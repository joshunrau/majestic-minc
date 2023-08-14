/* ----------------------------- MNI Header -----------------------------------
 * @NAME       :  volume_similarity
 * @DESCRIPTION:  an example of calculating volume similarity metrics 
 * @COPYRIGHT  :
 *              Copyright 2011 Vladimir Fonov, McConnell Brain Imaging Centre, 
 *              Montreal Neurological Institute, McGill University.
 *              Permission to use, copy, modify, and distribute this
 *              software and its documentation for any purpose and without
 *              fee is hereby granted, provided that the above copyright
 *              notice appear in all copies.  The author and McGill University
 *              make no representations about the suitability of this
 *              software for any purpose.  It is provided "as is" without
 *              express or implied warranty.
 * ---------------------------------------------------------------------------- */
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <algorithm>
#include <unistd.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>

using namespace minc;

void show_usage (const char * prog)
{
  std::cout<<"Estimate negative log- inter-label co-occurance matrix"<<std::endl
           <<"Usage: "<<prog<<" <input1.mnc> [... <inputN.mnc>] <output.csv> "<<std::endl
           <<"\t--mask <msk.mnc> use mask"<<std::endl
           <<"\t--classes <n> set number of classes"<<std::endl
           <<"\t--epsilon <f> set minimal probability"<<std::endl;
}

typedef unsigned short label_type;

int main(int argc,char **argv)
{
  int verbose=0;
  int clobber=0;
  std::string mask_f;
  double max_energy=20.0;
  double epsilon=exp(-max_energy);
  
  int number_of_classes=4;
  
  static struct option long_options[] = {
    {"verbose",       no_argument,       &verbose, 1},
    {"quiet",         no_argument,       &verbose, 0},
    {"clobber",       no_argument,       &clobber, 1},
    {"mask",          required_argument, 0       ,'m'},
    {"classes",       required_argument, 0       ,'c'},
    {"epsilon",       required_argument, 0       ,'e'},
    {0, 0, 0, 0}
  };
  
  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    int c = getopt_long (argc, argv, "vqm:c:", long_options, &option_index);
    
    /* Detect the end of the options. */
    if (c == -1) break;
    
    switch (c)
    {
      case 0:
        break;
      case 'm':
        mask_f=optarg;
        break;
      case 'c':
        number_of_classes=atoi(optarg);
        break;
      case 'e':
        epsilon=atof(optarg);
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage (argv[0]);
        return 1;
    }
  }
  
  if((argc - optind) < 2) {
    show_usage (argv[0]);
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
    minc_byte_volume cls;
    minc_byte_volume mask;
    
    std::vector<double> label_interaction(number_of_classes*number_of_classes,0.0);
      
    double total_count=0.0;
    double avg_probability=0.0;
    
    for(int i=0;i<(argc-optind);i++)
    {
      
      minc_1_reader rdr;
      
      if(verbose) 
        std::cout<<"loading :"<<argv[i+optind]<<std::endl;
      
      rdr.open(argv[i+optind]);
      
      load_simple_volume(rdr,cls);
      
      if(i==0)
      {
        if(!mask_f.empty())
        {
          minc_1_reader rdr2;
          if(verbose) 
            std::cout<<"loading mask:"<<mask_f.c_str()<<std::endl;
          rdr2.open(mask_f.c_str());
          load_simple_volume(rdr2,mask);
          if(mask.size()!=cls.size())
            REPORT_ERROR("Mask size mismatch");
        } else {
          mask.resize(cls.size());
          mask=1; //all volume is selected
        }
      }
      
      
      for(int z=1;z<cls.dim(2);z++)
        for(int y=1;y<cls.dim(1);y++)
          for(int x=1;x<cls.dim(0);x++)
          {
            if(mask.get(x,y,z))
            {
              int l1=cls.get(x,y,z);
              int l2=0;
              
              if(l1>=number_of_classes)
              {
                std::cerr<<"found label "<<l1<<" skipping"<<std::endl;
                continue;
              }
              
              if(mask.get(x-1,y,z)) {
                l2=cls.get(x-1,y,z);
                if(l2>=number_of_classes)
                {
                  std::cerr<<"found label "<<l1<<" skipping"<<std::endl;
                  continue;
                }
                label_interaction[l1+l2*number_of_classes]+=1;
                label_interaction[l2+l1*number_of_classes]+=1;
                total_count+=2;

                if(l1!=l2)
                  avg_probability+=2;
              }

              if(mask.get(x,y-1,z)) {
                l2=cls.get(x,y-1,z);
                if(l2>=number_of_classes)
                {
                  std::cerr<<"found label "<<l1<<" skipping"<<std::endl;
                  continue;
                }
                label_interaction[l1+l2*number_of_classes]+=1;
                label_interaction[l2+l1*number_of_classes]+=1;
                total_count+=2;
                if(l1!=l2)
                  avg_probability+=2;
              }
              
              if(mask.get(x,y,z-1)) {
                l2=cls.get(x,y,z-1);
                if(l2>=number_of_classes)
                {
                  std::cerr<<"found label "<<l1<<" skipping"<<std::endl;
                  continue;
                }
                label_interaction[l1+l2*number_of_classes]+=1;
                label_interaction[l2+l1*number_of_classes]+=1;
                total_count+=2;
                if(l1!=l2)
                  avg_probability+=2;
              }
              
            }
          }
    }
    
    avg_probability/=total_count;
    
    for(int i=0;i<number_of_classes*number_of_classes;i++)
    {
      label_interaction[i]/=total_count;
    }
    
    std::ofstream out(output_f.c_str());
    
    for(int l1=0;l1<number_of_classes;l1++)
    {
      for(int l2=0;l2<number_of_classes;l2++)
      {
        label_interaction[l2+l1*number_of_classes]=
            label_interaction[l2+l1*number_of_classes]>epsilon?
          -log(label_interaction[l2+l1*number_of_classes]):max_energy;
        
        if(l1==l2)
          label_interaction[l2+l1*number_of_classes]=0.0;
        
        out<<label_interaction[l2+l1*number_of_classes];
        
        if(l2<(number_of_classes-1))
          out<<",";
      }
      out<<std::endl;
    }
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  return 0;
}
