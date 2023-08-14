#ifdef HAVE_CONFIG_H
#include <config.h>
#endif //HAVE_CONFIG_H


#include <stdio.h>
#include <stdlib.h>

#ifdef MT_USE_OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 1
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif


#include <spams.h>



int main(int argc,char **argv)
{
    float * PatchImg;
    float * PatchTemp;
    
    /*init with eye*/
    int patch_volume=5;
    int count=10;
    int i,j;
    
    float noise=0.1;
#ifdef MT_USE_OPENMP
    fprintf(stderr,"Using OpenMP, max number of threads=%d\n",omp_get_max_threads());
#endif
    
    PatchImg=(float*)malloc(patch_volume*count*sizeof(float));
    PatchTemp=(float*)malloc(patch_volume*1*sizeof(float));
    
    
    for(i=0;i<count;i++)
    {
        for(j=0;j<patch_volume;j++)
        {
            if(i==j)
                PatchImg[j+i*patch_volume]=1.0;
            else if(i==(j+patch_volume))
                PatchImg[j+i*patch_volume]=-1.0;
            else
                PatchImg[j+i*patch_volume]=0.0;
            
            /*corrupt with noise*/
            PatchImg[j+i*patch_volume]+=drand48()*noise;
        }
    }
    
    for(j=0;j<patch_volume;j++)
    {   
        PatchTemp[j]=0;
    } 
    
    PatchTemp[2]=1.0;
    
    for(j=0;j<patch_volume;j++)
    {   
        PatchTemp[j]+=drand48()*noise;
    } 
    
    
    
    SpMatrix<float> alpha1,alpha2,alpha3,alpha4;
    Matrix<float> M(PatchImg,   patch_volume, count);
    Matrix<float> X(PatchTemp,  patch_volume, 1 );
    
    M.print("M");
    X.print("X");
    
    float constraint=1.0;
    float lambda2=0.15;
    
    
    /*TEST L1COEFFS*/
    lasso<float>(M, X, alpha1, 
                    count, constraint, lambda2, L1COEFFS , true, false);
    alpha1.print("alpha L1COEFFS");
    
    /**/
    for(i=0;i<count;i++)
        printf("%f ",alpha1[i]);
    printf("\n");
    
    constraint=0.15;
    lambda2=0.0;
    /*TEST PENALTY*/
    lasso<float>(M, X, alpha2, 
                    count, constraint, lambda2, PENALTY , true, false);
    alpha2.print("alpha PENALTY");
    
    /*TEST SPARSITY*/
    lasso<float>(M, X, alpha3, 
                    count, constraint, lambda2, SPARSITY , true, false);
    alpha3.print("alpha SPARSITY");
    
    /*TEST PENALTY2*/
    lasso<float>(M, X, alpha4, 
                    count, constraint, lambda2, PENALTY2 , true, false);
    alpha4.print("alpha PENALTY2");
    
    
    free(PatchImg);
    free(PatchTemp);
    
    return 0;
}
