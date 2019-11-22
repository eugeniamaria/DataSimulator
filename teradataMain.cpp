#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include <omp.h>
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

int main(int argc, char** argv)
{


  size_t K;         // Number of populations
  size_t Ns;        // Number of points/regions drawn from K-Simplex using dirichlet
  size_t I;         // Number of individuals -- N
  size_t L;         // Number of SNP locations -- M
  double MN;        // Minimum frequency
  int option;       // output type 1 for bed 0 for txt 
  char fname[1024]; // output file name
  
  int flag_K;
  int flag_Ns;
  int flag_I;
  int flag_L;
  int flag_MN;
  int flag_option;
  int flag_fname;
	
  //============================================================================                        
  // This part is taken from EVSL -- it allows us to give the arguments in any                          
  // particular order -- it requires a pre-defined value for each input argument                        
  //============================================================================                        
  int flg = findarg("help", NA, NULL, argc, argv);
  if (flg) {
    printf("\nUsage: ./GeneticDataSimulator -npop [int] -nregions [int] -nindividuals [int] -nSNP [int] -minfreq [double] -txtoutput [int] -filename [char]\n");
    printf("npop: # of simulated populations\n");
    printf("nregions: # of simulated regions\n");
    printf("nindividuals: # of simulated individuals \n");
    printf("nSNP: # of SNP markers\n");
    printf("minfreq: minimum frequency\n");
    printf("txtoutput: 0 for bed output (default), 1 for txt\n\n");
    printf("filename: name of output file");
    return 0;
  }
 
  flag_K      = findarg("npop",INT, &K, argc, argv);          if (!flag_K){K = 6;}
  flag_Ns     = findarg("nregions", INT, &Ns, argc, argv);    if (!flag_Ns){Ns = 5;}
  flag_I      = findarg("nindividuals", INT, &I, argc, argv); if (!flag_I){I = 100;}
  flag_L      = findarg("nSNP", INT, &L, argc, argv);         if (!flag_L){L = 10;}
  flag_MN     = findarg("minfreq", DOUBLE, &MN, argc, argv);  if (!flag_MN){MN = 1e-6;}
  flag_option = findarg("txtoutput",INT, &option, argc, argv);if (!flag_option){option = 0;}
  flag_fname  = findarg("filename",STR, fname, argc, argv);   
  //============================================================================
  
  clock_t bt;
  bt = omp_get_wtime();

  int padding = (int) log10 ((double) I)+1;
  int padding2 = (int) log10 ((double) L)+1;
 
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
  
  InitializeThread();

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
  /* Create the files                                                                            */
  /*==============================================================================================*/
  std::string fname1,fname2; 
  
  // Create the file to write
  
  if (flag_fname){
    if (option == 0){
      fname1 = string(fname) + ".ped";
      fname2 = string(fname) + ".map";
    }
    else{
      fname1 = string(fname) + ".txt";
    }
  }
  else{
    if (option == 0){
      fname1 = string("output_file") + ".ped";
      fname2 = string("output_file") + ".map";
    }
    else{
      fname1 = string("output_file") + ".txt";  
    }
  }


  FILE  *f,*f1; // f: ped or txt f1: map                                                                                        
  f = fopen(fname1.c_str(), "w+");
  
  if (f == NULL){
    printf("Error opening writing file!\n");
    exit(EXIT_FAILURE);
  }
  cout << "Files created!"<< endl;

  // Creating a row of the matrix
  double *S = (double *)malloc(K*sizeof(double)); // a row of the thetas (size K) 
  double *B = (double *)malloc(K*sizeof(double)); // a row of the betas  (size K)

  int theta_index = 0, beta_index = 0;
  double p; double bino_dis; int int_bino_dis;
    
  size_t chrnum, div, newdiv, pos;
  int flag = 1; 

  if (option == 0){
    f1  = fopen(fname2.c_str(), "w+");
    div = L/22; 
    newdiv = div;
    chrnum = 1;
  }
  
  pos = 750000;

  for(size_t i=0; i<I; i++){

    // the first 6 columns will include dump values
    if(option == 0){
      fprintf(f,"FAM%0*d ID%0*d %d %d %d %d ", padding, (int)i+1, padding, (int)i+1, 0, 0, 0, -9);
    }
   
    memcpy(S, &theta[theta_index], K*sizeof(double)); // get the row of thetas

    size_t count = 0; 
    for(size_t j=0; j<L; j++){ 
    
      memcpy(B,&B_kl[beta_index], K*sizeof(double)); // get the row of betas
      p = cblas_ddot(K, S, 1, B, 1);

      /* // Use Boost - No safe with threads
      bino_dis.param(BinomialDist::param_type(2,p));
      bino_dis = gen_from_binomial_BOOST(2, p);*/

      bino_dis = gen_from_binomial_GSL(2, p); //picks 0, 1 or 2
      int_bino_dis = (int) bino_dis;
     
      if(option == 0){
	if (flag == 1){
	  count++;
	  if (j < newdiv){
	    fprintf(f1, "%zu rs%0*d %d %zu", chrnum, padding2, (int)j+1, 0, pos);
	    pos++;
	  } 
	  else{
	    newdiv = newdiv + div; 
	    chrnum++;
	    fprintf(f1, "%zu rs%0*d %d %zu", chrnum, padding2, (int)j+1, 0, pos);
	  }
	}			
	if (int_bino_dis == 2)
	  fprintf(f,"%d %d",2,2);
	else if(int_bino_dis == 1) 
	  fprintf(f,"%d %d",1,2);
	else if(int_bino_dis == 0)
	  fprintf(f,"%d %d",1,1);
	else{
	  cout<<"unexpected sampled value "<<int_bino_dis<<" at i = "<<i<<" and j = "<<j<<endl;
	}			
      }
      else{
	fprintf(f,"%d",((int)bino_dis));
      }
      
      beta_index = beta_index + 1;
      
      if (j < L-1)
	fprintf(f," "); 
      if (option == 0 && flag == 1)
	fprintf(f1,"\r\n");   
    }
    
    beta_index = 0; 
    flag = 0; 
    theta_index = theta_index + K;
    fprintf(f,"\r\n");
  }
  
  std::cout << std::endl;
  fclose(f);
  
  clock_t et = omp_get_wtime();
  double time_taken = ((double)(et-bt)); // in seconds                                          
  printf("it took %8.5f seconds to execute \n", time_taken);

  delete[] B_kl;
  delete[] T_kn;
  
  return 0;
}
