/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2009 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#include "minc_1_rw.h"
#include <iostream>
#include <fstream>
#include <algorithm>

#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include <getopt.h>
#include <math.h>
#include <strtok.h>
#include <unistd.h>
#include "../gco/graph.cpp"
#include "../gco/maxflow.cpp"


using namespace minc;
 

void show_usage(const char *name)
{
  std::cerr 
      << "Usage: "<<name<<" <input> <init <output>  " << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--mask <mask.mnc> ROI mask"<<std::endl
      << "\t--debug print more info"<<std::endl
      << "\t--diagonal perform diagonal connection between voxels "<<std::endl
      << "\t--edge-smooth <f> use normalization for edge-map smoothing"<<std::endl
      << "\t--beta <f> apply factor to edge smoothing "<<std::endl;
}

class weigh_functor
{
public:
  virtual double operator()(double v1,double v2)=0;
};

class sym_inv_exp_distance:public weigh_functor
{
public:
  double sigma2;
  sym_inv_exp_distance(double sigma):sigma2(sigma*sigma)
  {}
  
  double operator()(double v1,double v2)
  {
    
    return ::exp(-(v1-v2)*(v1-v2)/sigma2);
  }
};

class bright_inv_exp_distance:public weigh_functor
{
public:
  double sigma2;
  bright_inv_exp_distance(double sigma):sigma2(sigma*sigma)
  {}
  
  double operator()(double v1,double v2)
  {
    
    return v1>v2?::exp(-(v1-v2)*(v1-v2)/sigma2):1.0;
  }
};

class dark_inv_exp_distance:public weigh_functor
{
public:
  double sigma2;
  dark_inv_exp_distance(double sigma):sigma2(sigma*sigma)
  {}
  
  double operator()(double v1,double v2)
  {
    
    return v2>v1?::exp(-(v1-v2)*(v1-v2)/sigma2):1.0;
  }
};

class linear_abs_distance:public weigh_functor
{
public:
  double sigma;
  linear_abs_distance(double sigma):sigma(sigma)
  {}
  
  double operator()(double v1,double v2)
  {
    double d=::fabs(v1-v2);
    return d>sigma?sigma/d:1.0;
  }
};

class squared_distance:public weigh_functor
{
public:
  double sigma2;
  squared_distance(double sigma):sigma2(sigma*sigma)
  {}
  
  double operator()(double v1,double v2)
  {
    double difference=(v1-v2)*(v1-v2);
    return difference>sigma2?sigma2/difference:1.0;
  }
};

int main(int argc,char **argv)
{
  int clobber=0;
  int verbose=0;
  int normalize=0;
  int debug=0;
  int expansion=0;
  int diagonal=0;
  int linear_distance=0;
  
  std::string mask_f;
  std::string init_f;
  std::string intensity_f;
  
  double weight=1.0;
  
  double quantization=1.0;
  double epsilon=1e-12;
  
  double edge_smooth=0.0;

  const double SOURCE_LINK=1e20,SINK_LINK=1e20;
  
  static struct option long_options[] =
  {
    {"verbose",     no_argument, &verbose, 1},
    {"debug",       no_argument, &debug, 1},
    {"quiet",       no_argument, &verbose, 0},
    {"clobber",     no_argument, &clobber, 1},
    {"diagonal",   no_argument, &diagonal, 1},
    {"mask",        required_argument, 0, 'm'},
    {"edge-smooth",     required_argument, 0       ,'E'},
    {0, 0, 0, 0}
  };
  
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "m:i:s:I:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
      case 0:
        break;
      case 'E':
        edge_smooth=atof(optarg);
        break;
      case 'm':
        mask_f=optarg;
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage(argv[0]);
        return 1;
    }
  }

  if ((argc - optind) < 3)
  {
    show_usage(argv[0]);
    return 1;
  }

  std::string output_f=argv[argc-1]; //last is output
  argc--;

  intensity_f=argv[optind];
  init_f=argv[optind+1];
    
  if (!clobber && !access (output_f.c_str (), F_OK))
  {
    std::cerr<<"File "<<output_f.c_str()<<" exists, use --clobber"<<std::endl;
    return 1;
  }
 
  try
  {
    minc_float_volume intensity;
    minc_byte_volume mask,init,cls;

    minc_1_reader rdr;
    if(verbose) std::cout<<"loading intensity:"<<intensity_f.c_str()<<std::endl;
    rdr.open(intensity_f.c_str());
    
    load_simple_volume(rdr,intensity);
    
    minc_1_reader rdr2;
    if(verbose) std::cout<<"loading init:"<<init_f.c_str()<<std::endl;
    rdr2.open(init_f.c_str());
    load_simple_volume(rdr2,init);
    if(init.size()!=intensity.size())
      REPORT_ERROR("Init size mismatch");

    if(!mask_f.empty())
    {
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading mask:"<<mask_f.c_str()<<std::endl;
      rdr.open(mask_f.c_str());
      load_simple_volume(rdr,mask);
      if(mask.size()!=intensity.size())
        REPORT_ERROR("Mask size mismatch");
    } else {
      mask.resize(intensity.size());
      mask=1; //all volume is selected
    }

    int num_voxels=mask.size().vol();
    int num_edges=num_voxels*(diagonal?16:6);

    typedef Graph<double,double,double> GraphType;

    GraphType gc( num_voxels, num_edges); 

    gc.add_node(intensity.dim(0)*intensity.dim(1)*intensity.dim(2));

    size_t num_sources=0;
    size_t num_sinks=0;

    for(int z=0;z<intensity.dim(2);z++)
      for(int y=0;y<intensity.dim(1);y++)
        for(int x=0;x<intensity.dim(0);x++)
        {
          int node0=x+y*mask.dim(0)+z*mask.dim(0)*mask.dim(1);
          if(init.get(x,y,z)==1) //SOURCE
          {
            num_sources++;
            gc.add_tweights(node0,SOURCE_LINK,0);
          } else if(init.get(x,y,z)==2) //sink 
          {
            num_sinks++;
            gc.add_tweights(node0,0,SINK_LINK);
          }
        }

    if(num_sinks==0)
    {
      std::cerr<<"No sink nodes defines ! (init image label 2)"<<std::endl;
      return 1;
    }

    if(num_sources==0)
    {
      std::cerr<<"No source nodes defined ! (init image label 1)"<<std::endl;
      return 1;
    }
    std::cout<<"num_sources="<<num_sources<<" num_sinks="<<num_sinks<<std::endl;
    
    int weight1=quantization;
    int weight2=quantization/sqrt(2.0); //(100/sqrt(2))
    int weight3=quantization/sqrt(3.0); //(100/sqrt(3))

    weigh_functor* dist_weight;
    if(linear_distance)
      dist_weight=new linear_abs_distance(edge_smooth);
    else
      dist_weight=new sym_inv_exp_distance(edge_smooth);


    if(edge_smooth==0.0) //use MAD to estimate the smoothing term
    {
      std::vector<double> gradients;
      for(int z=1;z<mask.dim(2);z++)
        for(int y=1;y<mask.dim(1);y++)
          for(int x=1;x<mask.dim(0);x++)
          {
            double i0=intensity.get(x,y,z);
            
            double i1=intensity.get(x-1,y,z);
            double i2=intensity.get(x,y-1,z);
            double i3=intensity.get(x,y,z-1);
            
            gradients.push_back(::fabs(i0-i1));
            gradients.push_back(::fabs(i0-i2));
            gradients.push_back(::fabs(i0-i3));
            
          }
      std::sort(gradients.begin(),gradients.end());
      edge_smooth=1.4826*gradients[gradients.size()/2];
      std::cout<<"edge_smooth="<<edge_smooth<<std::endl;
    }
    // now set up a grid neighborhood system, 6-connected neighborhood for now
    for(int z=1;z<mask.dim(2);z++)
      for(int y=1;y<mask.dim(1);y++)
        for(int x=1;x<mask.dim(0);x++)
    {
      int node0=x+y*mask.dim(0)+z*mask.dim(0)*mask.dim(1);

      int node1=x-1+y*mask.dim(0)+z*mask.dim(0)*mask.dim(1);
      int node2=x+(y-1)*mask.dim(0)+z*mask.dim(0)*mask.dim(1);
      int node3=x+y*mask.dim(0)+(z-1)*mask.dim(0)*mask.dim(1);
      
      double i0=intensity.get(x,y,z);
      
      double i1=intensity.get(x-1,y,z);
      double i2=intensity.get(x,y-1,z);
      double i3=intensity.get(x,y,z-1);
    
      gc.add_edge(node0,node1,weight1*(*dist_weight)(i1,i0),weight1*(*dist_weight)(i0,i1));
      gc.add_edge(node0,node2,weight1*(*dist_weight)(i2,i0),weight1*(*dist_weight)(i0,i2));
      gc.add_edge(node0,node3,weight1*(*dist_weight)(i3,i0),weight1*(*dist_weight)(i0,i3));

    
      if(diagonal)
      {
      
        int node1=(x-1)+(y-1)*mask.dim(0)+z*mask.dim(0)*mask.dim(1);
        int node2=(x-1)+y*mask.dim(0)+(z-1)*mask.dim(0)*mask.dim(1);
        int node3=x+(y-1)*mask.dim(0)+(z-1)*mask.dim(0)*mask.dim(1);
        int node4=(x-1)+(y-1)*mask.dim(0)+(z-1)*mask.dim(0)*mask.dim(1);
        
        double i1=intensity.get(x-1,y-1,z);
        double i2=intensity.get(x-1,y,z-1);
        double i3=intensity.get(x,y-1,z-1);
        double i4=intensity.get(x-1,y-1,z-1);

        gc.add_edge(node0,node1,weight2*(*dist_weight)(i1,i0),weight2*(*dist_weight)(i0,i1));
        gc.add_edge(node0,node2,weight2*(*dist_weight)(i2,i0),weight2*(*dist_weight)(i0,i2));
        gc.add_edge(node0,node3,weight2*(*dist_weight)(i3,i0),weight2*(*dist_weight)(i0,i3));

        gc.add_edge(node0,node4,weight3*(*dist_weight)(i4,i0),weight3*(*dist_weight)(i0,i4));
        
        if(y<(mask.dim(1)-1))
        {
          double i1=intensity.get(x-1,y+1,z);
          int node1=(x-1)+(y+1)*mask.dim(0)+z*mask.dim(0)*mask.dim(1);
          gc.add_edge(node0,node1,weight2*(*dist_weight)(i1,i0),weight2*(*dist_weight)(i0,i1));
        }
        
        if(z<(mask.dim(2)-1))
        {
          double i1=intensity.get(x-1,y,z+1);
          int node1=(x-1)+y*mask.dim(0)+(z+1)*mask.dim(0)*mask.dim(1);
          gc.add_edge(node0,node1,weight2*(*dist_weight)(i1,i0),weight2*(*dist_weight)(i0,i1));
        }
      }
    }

    std::cout<<"Flow="<<gc.maxflow()<<std::endl;

    cls.resize(intensity.size());

    for(int z=0;z<intensity.dim(2);z++)
      for(int y=0;y<intensity.dim(1);y++)
        for(int x=0;x<intensity.dim(0);x++)
        {
          int node0=x+y*mask.dim(0)+z*mask.dim(0)*mask.dim(1);
          cls.set(x,y,z,gc.what_segment(node0,GraphType::SINK)==GraphType::SOURCE?1:0);
        }

    save_volume(cls,output_f.c_str(),intensity_f.c_str());
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  return 0;
}
