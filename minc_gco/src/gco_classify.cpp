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
#include <time_stamp.h>    // for creating minc style history entry

#include <iostream>
#include <fstream>
#include <algorithm>

#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include <getopt.h>
#include <math.h>
#include <unistd.h>

#include "utils.h"

#include "strtok.h"
#include "gco/GCoptimization.h"

using namespace minc;


class weigh_functor
{
public:
  virtual double operator()(double difference)=0;
};


class inv_exp_distance:public weigh_functor
{
public:
  double sigma2;
  inv_exp_distance(double sigma):sigma2(sigma*sigma)
  {}
  double operator()(double difference)
  {
    return ::exp(-difference*difference/sigma2);
  }
};

class linear_abs_distance:public weigh_functor
{
public:
  double sigma;
  linear_abs_distance(double sigma):sigma(sigma)
  {}
  
  double operator()(double difference)
  {
    return fabs(difference)/sigma;
  }
};

class squared_distance:public weigh_functor
{
public:
  double sigma2;
  squared_distance(double sigma):sigma2(sigma*sigma)
  {}
  
  double operator()(double difference)
  {
    return difference*difference/sigma2;
  }
};


template<class T>void index_to_coord(T& array, size_t index, fixed_vec<3,int>& c)
{
  c[2]=index/(array.dim(1)*array.dim(2));
  index %= ( array.dim(1)*array.dim(2) );
  c[1]=index/(array.dim(1));
  index %= array.dim(1);
  c[0]=index;
}

void show_usage(const char *name)
{
  std::cerr 
      << "Usage: "<<name<<" <input_1> .. <input_n> <output>  " << std::endl
      << "\t--csv <input.csv> use files from csv file" << std::endl
      << "\t--ignore-missing ignore missing files, set probability to zero" << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--mask <mask.mnc> ROI mask"<<std::endl
      << "\t--iter <n> maximum number of iterations"<<std::endl
      << "\t--debug print more info"<<std::endl
      << "\t--smooth <f> smoothing term, alternative to --cooc"<<std::endl
      << "\t--cooc <label.txt> inter-label coocurance matrix -log , comma separated table"<<std::endl
      << "\t--epsilon <f> smallest thesholded probability"<<std::endl
      << "\t--expansion use alpha-expansion instead of swap, probably will not work if cooc matrix is used"<<std::endl
      << "\t--quant <n> quantization, used to convert floating point to integer, default 1000"<<std::endl
      << "\t--wdata <f> use this weight for data term, default 1.0"<<std::endl
      << "\t--wlabel <f> use this weight for label-interaction, default 1.0"<<std::endl
      << "\t--wintensity <f> use this weight for intensity gradient term, default 1.0"<<std::endl
      << "\t--diagonal perform diagonal connection between voxels "<<std::endl
      << "\t--intensity <input> use intensity difference for connection costs"<<std::endl
      << "\t--edge-smooth <f> sigma for edge weight map"<<std::endl
      << "\t--linear-distance use linear distance energy, instead of inv exp"<<std::endl
      << "\t--ball <n> use binary ball of the given radius for neighborhood"<<std::endl;
}

int main(int argc,char **argv)
{
  char *_history = time_stamp(argc, argv); //maybe we should free it afterwards
  std::string history=_history;
  free(_history);
  
  int clobber=0;
  int verbose=0;
  int normalize=0;
  int debug=0;
  int expansion=0;
  int diagonal=0;
  int linear_distance=0;
  int ignore_missing=0;
  
  int max_dist=1;
  int max_sq_dist=1;

  std::string csv_f;
  std::string mask_f;
  std::string init_f;
  std::string intensity_f;
  std::string label_energy_f;
  int number_of_classes=2;
  int maxiter=10;

  double wdata=1.0;
  double wlabel=1.0;
  double wintensity=1.0;
  
  double quantization=1000;
  double smoothing=0.301; //corresponds to -log(0.5)
  double epsilon=1e-12;
  double label_norm=1.0/6; // normalization factor for 6-connected neighborhood
  
  double edge_smooth=1.0;
  
  static struct option long_options[] =
  {
    {"verbose",     no_argument, &verbose, 1},
    {"debug",       no_argument, &debug, 1},
    {"quiet",       no_argument, &verbose, 0},
    {"clobber",     no_argument, &clobber, 1},
    {"expansion",   no_argument, &expansion, 1},
    {"diagonal",    no_argument, &diagonal, 1},
    {"linear-distance", no_argument, &linear_distance, 1},
    {"ignore-missing", no_argument, &ignore_missing, 1},
    {"mask",        required_argument, 0,       'm'},
    {"iter",        required_argument, 0,       'i'},
    {"smooth",      required_argument, 0,       's'},
    {"cooc",        required_argument, 0,       'C'},
    {"epsilon",     required_argument, 0       ,'e'},
    {"init" ,       required_argument, 0       ,'I'},
    {"wdata"  ,     required_argument, 0       ,'D'},
    {"wlabel" ,     required_argument, 0       ,'w'},
    {"quant" ,      required_argument, 0       ,'q'},
    {"wintensity" , required_argument, 0       ,'W'},
    {"intensity",   required_argument, 0       ,'N'},
    {"edge-smooth", required_argument, 0       ,'E'},
    {"csv",         required_argument, 0       ,'c'},
    {"ball",        required_argument, 0       ,'B'},
    {0, 0, 0, 0}
  };
  
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "m:i:s:I:E:w:W:C:q:D:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
      case 0:
        break;
      case 'N':
        intensity_f=optarg;
        break;
      case 'e':
        epsilon=atof(optarg);
        break;
      case 'E':
        edge_smooth=atof(optarg);
        break;
      case 'm':
        mask_f=optarg;
        break;
      case 'i':
        maxiter=atoi(optarg);
        break;
      case 's':
        smoothing=atof(optarg);
        break;
      case 'C':
        label_energy_f=optarg;
        break;
      case 'I':
        init_f=optarg;
        break;
      case 'c':
        csv_f=optarg;
        break;
      case 'q':
        quantization=atof(optarg);
        break;
      case 'D':
        wdata=atof(optarg);
        break;
      case 'w':
        wlabel=atof(optarg);
        break;
      case 'B':
        max_dist=atoi(optarg);
        max_sq_dist=max_dist*max_dist;
        break;
      case 'W':
        wintensity=atof(optarg);
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage(argv[0]);
        return 1;
    }
  }

  if ((argc - optind) < 1)
  {
    std::cout<<"argc="<<argc<<" optind="<<optind<<std::endl;
    show_usage(argv[0]);
    return 1;
  }

  std::string output_f=argv[argc-1]; //last is output
  argc--;

  if (!clobber && !access (output_f.c_str (), F_OK))
  {
    std::cerr<<"File "<<output_f.c_str()<<" exists, use --clobber"<<std::endl;
    return 1;
  }
  
  try
  {
    strings _input_f;
    
    if(!csv_f.empty())
    {
      string_table _table;
      read_table_n(csv_f.c_str(),_table,2);
      table_verify(_table);

      if(_table[0].size()!=2)
      {
        std::cerr<<"Unexpected input file format, expecting two columns in csv file"<<std::endl;
        return 1;
      }
      number_of_classes=0;
      for(size_t i=0;i<_table.size();i++)
      {
        int c=atoi(_table[i][0].c_str());
        if(number_of_classes <= c)
          number_of_classes=c+1;
      }
      _input_f.resize(number_of_classes);

      for(size_t i=0;i<_table.size();i++)
      {
        int c=atoi(_table[i][0].c_str());
        _input_f[c]=_table[i][1];
      }
    } else {
      for(int i=optind;i<argc;i++)
        _input_f.push_back(std::string(argv[i]));

      number_of_classes=_input_f.size();
    }
    if(verbose)
      std::cout<<"Number of classes:"<<number_of_classes<<std::endl;

    double label_norm=1.0/6; // normalization factor for 6-connected neighborhood

    if(diagonal)
    {
      label_norm=1.0/26; // 26-connected neighborhood
      max_sq_dist=3;
      max_dist=1;
    } else if(max_dist>1) {
      label_norm=1.0/(max_dist*max_dist*max_dist*8-1);
    }

    doubles_table classes_interaction_table;
    
    if(!label_energy_f.empty())
    {
      read_doubles_table(label_energy_f.c_str(), classes_interaction_table);
      if(!table_verify(classes_interaction_table))
        return 1;
    }
    
    //1  load all the volumes
    volumes input;
    if(verbose) std::cout<<"loading input files"<<std::endl;

    load_volumes(_input_f,input,0,verbose,ignore_missing);
    check_volumes(input,ignore_missing);
    
    minc_byte_volume mask,cls;
    
    if(!mask_f.empty())
    {
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading mask:"<<mask_f.c_str()<<std::endl;
      rdr.open(mask_f.c_str());
      load_simple_volume(rdr,mask);
      if(mask.size()!=input[0].size())
        REPORT_ERROR("Mask size mismatch");
    } else {
      mask.resize(input[0].size());
      mask=1; //all volume is selected
    }
    
    int num_voxels=mask.size().vol();


    if(verbose)
      std::cout<<"Using descrete quantization multiplier: "<<quantization<<std::endl
               <<"Label-energy weight: "<<wlabel<<std::endl
               <<"Data-energy weight: "<<wdata<<std::endl
               <<"Intensity-difference weight: "<<wintensity<<std::endl;
    
    // first set up the array for data costs
    int *data = new int[num_voxels*number_of_classes];
    for(int i=0;i<mask.c_buf_size();i++)
    {
      if(mask.c_buf()[i]==0)
      {
        data[i*number_of_classes]=0;
        for(int l=1;l<number_of_classes;l++)
            data[i*number_of_classes+l] = GCO_MAX_ENERGYTERM;
      } else {
        for(int l=0;l<number_of_classes;l++)
        {
          double v=input[l].c_buf()[i];
          data[i*number_of_classes+l] = v>epsilon?floor(-log(v)*wdata*quantization):GCO_MAX_ENERGYTERM;
        }
      }
    }

    
    int *smooth = new int[number_of_classes*number_of_classes];
    
    if(classes_interaction_table.empty())
    {
      if(verbose)
        std::cout<<"Using uniform label-interaction cost of "<<smoothing<<std::endl;
      
      // next set up the array for smooth costs
      for ( int l1 = 0; l1 < number_of_classes; l1++ )
        for (int l2 = 0; l2 < number_of_classes; l2++ )
          smooth[l1+l2*number_of_classes] = (l1==l2 ? 0: (1*smoothing*wlabel) );
    } else {
      if(verbose)
        std::cout<<"Using label-interaction from "<<label_energy_f.c_str()<<std::endl;
      
      if(classes_interaction_table.size()!=classes_interaction_table[0].size())
      {
        std::cerr<<"Label-interaction matrix must be square!"<<std::endl;
        return 1;
      }
      if(number_of_classes!=classes_interaction_table.size())
      {
        std::cerr<<"Label-interaction matrix expected to be "<<number_of_classes<<"x"<<number_of_classes<<std::endl;
        return 1;
      }
      
      for ( int l1 = 0; l1 < number_of_classes; l1++ )
      {
        for (int l2 = 0; l2 < number_of_classes; l2++ )
        {
          smooth[l1+l2*number_of_classes] = ::floor(quantization*classes_interaction_table[l1][l2]*wlabel);
        }
      }
    }
    
    // setup GC graph
    GCoptimizationGeneralGraph gc(num_voxels,number_of_classes);
    
    gc.setDataCost(data);
    gc.setSmoothCost(smooth);
    
    if(!init_f.empty())
    {
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading init:"<<init_f.c_str()<<std::endl;
      rdr.open(init_f.c_str());
      load_simple_volume(rdr,cls);
      if(cls.size()!=input[0].size())
        REPORT_ERROR("Init size mismatch");

      for(size_t i=0;i<mask.c_buf_size();i++)
          {
            gc.setLabel(i,cls.c_buf()[i]);
          }
      
    }
    
    minc_float_volume intensity;
    bool use_intensity=false;
    
    weigh_functor* dist_weight=NULL;
    
    if(!intensity_f.empty())
    {
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading intensity:"<<intensity_f.c_str()<<std::endl;
      rdr.open(intensity_f.c_str());
      
      load_simple_volume(rdr,intensity);
      
      
      if(linear_distance)
        dist_weight=new linear_abs_distance(edge_smooth);
      else
        dist_weight=new inv_exp_distance(edge_smooth);
      
      if(intensity.size()!=input[0].size())
        REPORT_ERROR("Init size mismatch");
      use_intensity=true;
    }
    
    std::vector<double> weight(max_sq_dist);
    for(size_t i=0;i<max_sq_dist;i++)
      weight[i]=quantization/::sqrt((double)(i+1));
    
    size_t pairs=0;
    size_t i=0;
    // now set up a grid neighborhood system
    for(int z=0;z<mask.dim(2);z++)
      for(int y=0;y<mask.dim(1);y++)
        for(int x=0;x<mask.dim(0);x++)
    {
      int i=x+y*mask.dim(0)+z*mask.dim(1)*mask.dim(0);

      //TODO: add mask ?
      for(int zz=-max_dist;zz<=max_dist;zz++)
        for(int yy=-max_dist;yy<=max_dist;yy++)
          for(int xx=-max_dist;xx<=max_dist;xx++)
          {
            int _x=x+xx;
            int _y=y+yy;
            int _z=z+zz;
            
            if(_x<0 || _x >= mask.dim(0) ||
               _y<0 || _y >= mask.dim(1) ||
               _z<0 || _z >= mask.dim(2) ) continue;
            
            int j= _x + _y*mask.dim(0) + _z*mask.dim(0)*mask.dim(1);

            if( j<=i ) continue;

            int dist=xx*xx+yy*yy+zz*zz;

            if(dist>max_sq_dist) continue;
      
            //pairs++;
            double weight_intensity=1.0;

            if(use_intensity)
            {
              weight_intensity=(*dist_weight)(intensity.c_buf()[j]-intensity.c_buf()[j]);
            }
            gc.setNeighbors(i, j, ::floor(wintensity*weight[dist-1]*weight_intensity));
          }
    }

   
    if(verbose)
      std::cout<<"Before optimization energy is "<<gc.compute_energy()<<std::endl;
    
    if(expansion)
      gc.expansion(maxiter);
    else
      gc.swap(maxiter);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);

    if(verbose)
      std::cout<<"After optimization energy is "<<gc.compute_energy()<<std::endl;

    cls.resize(input[0].size());
    
    for ( int  i = 0; i < num_voxels; i++ )
      cls.c_buf()[i] = gc.whatLabel(i);
    
    save_volume(cls,output_f.c_str(),_input_f[0].c_str(),history);
    //de-allocate data arrays
    delete[] data;
    delete[] smooth;
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  catch (GCException e){
           e.Report();
  }
}
