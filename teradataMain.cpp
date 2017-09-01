#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include <boost/random/random_device.hpp>
#include "normal.hpp"
#include <boost/random/beta_distribution.hpp>
#include <boost/random.hpp>
#include <boost/random/variate_generator.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include "utilities.hpp"

using namespace std;

int main (){
  
  size_t K  =    6; // Number of populations
  size_t Ns =    5; // Number of points/regions drawn from K-Simplex using dirichlet
  size_t I  =   10; // Number of individuals -- N  
  size_t L  =  100; // Number of SNP locations -- M
  double MN = 1e-6; // Minimum frequency
  
  /*==============================================================================================*/
  /* Read files fst and freq                                                                      */
  /*==============================================================================================*/
  ifstream fst("fst.txt",ios::in);
  ifstream freq("freq.txt",ios::in);

  if (!fst){
    cout << "Cannot open the F_st file.\n";
    return -1;
  }
  if (!freq){
    cout << "Cannot open the Freq file.\n";
    return -1;
  }

  vector<double> valuesfst;
  vector<double> valuesfreq;

  size_t fstsize = ReadFiles(valuesfst,fst);
  size_t freqsize = ReadFiles(valuesfreq,freq);

  if (fstsize != freqsize){
    cout << "The files are not of the same size! Please check and return!"<< endl;
    fst.close();
    freq.close();
    return -1;
  }
  
  size_t filesize = fstsize;
  
  fst.close();
  freq.close();
  

  /*==============================================================================================*/
  /* Betas computation                                                                            */
  /*==============================================================================================*/
 
  size_t BetaMatrixSize = L*K;
  double *B_kl = new double[BetaMatrixSize];
  size_t randindex;

  for (size_t i=0; i < L; i++){
    randindex = rand() % filesize; // Generate a random number in the range 0 to filesize-1
    ComputeBetas(B_kl, i, (double) valuesfst.at(randindex),(double) valuesfreq.at(randindex),K);
  }
  
  
  /*==============================================================================================*/              
  /* Thetas computation Schenario A                                                               */              
  /*==============================================================================================*/ 

  size_t ThetaMatrixSize = K*I;
  double *T_kn = new double[ThetaMatrixSize];
 
  // Dirichlet distribution options 
  double alpha = 0.2; // Dirichlet Parameter --> controls sparsity
  double gamma =  50; // Second level Dirichlet scale 
  
  ComputeThetasA(T_kn, alpha, gamma, MN, Ns, K, I);
  

  /*==============================================================================================*/              
  /* Thetas computation Schenario B                                                               */              
  /*==============================================================================================*/ 

  size_t s = 2;
  double a = 0.0, b;
  double *theta = (double *)malloc(I*K*sizeof(double));
  
  b = K + 1.00;
  ComputeThetasB(I,K,s,a,b,theta);

  /*==============================================================================================*/
  /* Cretate the files                                                                            */
  /*==============================================================================================*/
  
  // Create the file to write
  FILE  *f;
  // this one for writting output                                                                  
  f = fopen("output.txt", "w");
  if (f == NULL)
    {
      printf("Error opening writting file!\n");
      exit(EXIT_FAILURE);
    }


  double *S = (double *)malloc(K*sizeof(double));
  double *B = (double *)malloc(K*sizeof(double));

  int theta_index = 0, beta_index = 0;
  double p; double bino_dis;
  clock_t t;
  t = clock();
  for(int i=0;i<I;i++){

    memcpy(S,&theta[theta_index], K*sizeof(double) );
    p = 0.0;

    for(int j=0;j<L;j++){

      memcpy(B,&B_kl[beta_index], K*sizeof(double) );
      p = cblas_ddot(K, S, 1, B, 1);

      //bino_dis.param(BinomialDist::param_type(2,p));

      bino_dis=gen_from_binomial(2, p);
      //print to output file                                                           
      fprintf(f,"%d ",((int)bino_dis+1));
      //cout << "Bino_dis+1 = " << (int)bino_dis+1 << endl;

      //printf("prob is %f \n",p);                                                     
      //printf("number is %d\n",temp[j]);                                              
      p = 0.0;
      beta_index = beta_index + 1;
      //if(j==1) break;                                                                
    
    }

    
    beta_index = 0;
    theta_index = theta_index + K;
    fprintf(f,"\r\n");
    //if(i==0)break;
    break;
  }
  fclose(f);
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds                                          
  printf("it took %f seconds to execute \n", time_taken);

 
  delete[] B_kl;
  delete[] T_kn;
  
  return 0;
}
