#include "utilities.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
gsl_rng ** threadvec = new gsl_rng*[omp_get_max_threads()];

// parse command-line input parameters                                                                  
// this part is copied from http://www-users.cs.umn.edu/~saad/software/EVSL/                            
int findarg(const char *argname, ARG_TYPE type, void *val, int argc, char **argv) {
  int *outint;
  double *outdouble;
  char *outchar;
  int i;
  
  for (i=0; i<argc; i++) {
    if (argv[i][0] != '-') {
      continue;
    }
    if (!strcmp(argname, argv[i]+1)) {
      if (type == NA) {
        return 1;
      } else {
        if (i+1 >= argc /*|| argv[i+1][0] == '-'*/) {
          return 0;
        }

        switch (type) {
        case INT:
          outint = (int *) val;
          *outint = atoi(argv[i+1]);
          return 1;
          break;
        case DOUBLE:
          outdouble = (double *) val;
          *outdouble = atof(argv[i+1]);
          return 1;
          break;
        case STR:
          outchar = (char *) val;
          sprintf(outchar, "%s", argv[i+1]);
          //printf("%s",argv[i+1]);
          return 1;
          break;
        default:
          printf("unknown arg type\n");
        }
      }
    }
  }
  return 0;
}

/*===================================================================================*/
/* This computes the betas                                                           */
/*===================================================================================*/
void ComputeBetas(double *B_kl, size_t row, double fst, double freq, size_t K){
  //===============================================================//
  // fst : Wright's F_st for each SNP location                     //
  // freq : marginal allele frequency for each SNP location        //
  // K : number of populations                                     //
  // ==============================================================//

  size_t j;
  double parameter1, parameter2;
  
  // If fst is less than 1e-6 we replicate the frequence
  if(fst < 1e-6){
    for (j=0; j<K; j++){ B_kl[j+K*row] = freq;} 
  }
  else{
    parameter1 = freq*(1-fst)/fst;
    parameter2 = (1-freq)*(1-fst)/fst;

    for (j=0; j<K; j++){
      do{
	B_kl[j+K*row] = gen_from_beta(parameter1,parameter2);
      }
      while(B_kl[j+K*row]<=0.05 || B_kl[j+K*row]>=0.95); // redraw if B_kl is less than 0.05 or more than 0.95
    }
  }
}

 
/*===================================================================================*/
/* This computes the thetas (Scenario A)                                             */
/*===================================================================================*/
void ComputeThetasA(double *TB_kn, double alpha, double gamma, double MN, size_t NS, size_t K, size_t I ){

  size_t pos, i, j, k, l;
  size_t index        = 0;
  size_t probNSsize   = NS*K;
  size_t g            = I/NS; // uniform distribution of the individuals into the regions
  size_t probIsize    = K*g;
  size_t matsize      = K*I;
  
  double f            = MN/(1-(MN*K));
  double sumt;

  size_t *indexvector = new size_t[probIsize];
  double *A_ns        = new double[probNSsize];  
  double *res         = new double[probNSsize]; 
  double *res2        = new double[probIsize];
  double *TB_knt      = new double[matsize];
  double *A_ns2       = new double[K];
  double *buffer      = new double[g];
  double *temp        = new double[K];

  // Create the NS regions -- each region has proportions of the K ancestral populations
  for (i=0; i<probNSsize; i++){ A_ns[i] = alpha; }
  gen_from_dirichlet(NS,K,A_ns,res,temp);

  for (k=0; k<K; k++){
    for(l=0; l<g; l++){
      indexvector[index] = l*K+k;
      index++;}
  }


  // Generate thetas -- probability of individuals to be from populations
  for (i=0; i<NS; i++){                             // A_ns2 = K x 1 ; res = NS x K ; res2 = I/NS x K
    for (j=0; j<K; j++){ A_ns2[j] = res[i*K +j];}   // bring the i-th row of res (length of row=K)
    cblas_dscal(K,gamma,A_ns2,1);                   // create the new probability vector
    gen_from_dirichlet(g,K,A_ns2,res2,temp);        // draw I/NS items using A_ns2 as probability vector
    for (k=0; k<probIsize; k++){res2[k]=res2[k]+f;} // normalize by adding f
    
    for (l=0; l<g; l++){buffer[l] = 0;}
    colsum(g, K, res2, buffer);                     // sum of the columns
    normprob(g, K, res2 ,buffer);                   // create the probability
    
    
    for ( l=0; l<g; l++){
      sumt=0;
      for (k=0; k<K; k++){sumt = sumt + res2[l*K+k];}
    }
    

    for (l=0; l<probIsize; l++){
      pos = indexvector[l]; 
      TB_kn[i*probIsize+l] = res2[pos];             // store the result in TB_kn
      TB_knt[i*probIsize+l] = res2[l]; 
    }
  }

  delete[] res;
  delete[] temp;
  delete[] A_ns;
  delete[] res2;
  delete[] A_ns2;
  delete[] buffer;
  delete[] TB_knt;
  delete[] indexvector;
}

/*===================================================================================*/
/* This computes the thetas (Scenario B)                                             */
/*===================================================================================*/
void ComputeThetasB(size_t n, size_t k, size_t s, double a, double b, double *Q){
  
  double *xs = (double *)malloc(n*sizeof(double));
  double *q = (double *)malloc(k*sizeof(double));
  size_t i,m,qindex = 0;
  double sumq = 0.0;

  for(i=0;i<n;i++){
    /*                                                                                       
     specify where each individual falls based on [a,b]                           
     xs <- a + (0:(n-1))/(n-1)*(b-a)                                                 
    */
    xs[i] = a +  ((double)i/(((double)n-1))*(b-a));

    for (m=1;m<=k;m++){
      boost::math::normal_distribution<double> d(m,s);
      q[m-1] = pdf(d, xs[i]);
      sumq = sumq + q[m-1];
    }
    for (m=0;m<k;m++){
      Q[qindex] = q[m] / sumq;
      qindex++;
    }
    sumq = 0.0;
  }
  
  delete[] xs;
  delete[] q;
}


/*===================================================================================*/
/* Read the text files   -- reads multiple digit numbers                             */
/*===================================================================================*/
size_t ReadFiles(vector<double>& values, ifstream & filename)
{
  double num;
  size_t numoflines=0;

  while(filename >> num){
    values.push_back(num);
    numoflines++;
  } 

  return numoflines;
}


/*===================================================================================*/
/* This prints the read values                                                       */
/*===================================================================================*/
void printBlock(double *dbuffer,size_t bufferlength)
{
  size_t i;

  cout << ": ";
  for ( i=0;i<bufferlength;i++){
    cout << i << ':' << dbuffer[i] << '\n'; }
  cout << endl;
  cout << "The length " << bufferlength << endl; 
}
/*===================================================================================*/
/* Initialize Threading		                                                         */
/*===================================================================================*/
void InitializeThread(){
	
  gsl_rng_env_setup();	
  for (int b = 0; b < omp_get_max_threads(); b++){
    threadvec[b] = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(threadvec[b],b*clock());
  }
	
  std::cout << "number of threads : " << omp_get_max_threads() << std::endl;
	
}
/*===================================================================================*/
/*Random Number Generator from Beta Distribution--Using BOOST C++ Library            */
/*===================================================================================*/
double gen_from_beta(double param1, double param2)
{
  using namespace boost::random;

  // construct the distribution
  using boost::random::beta_distribution;
  beta_distribution<>mybeta(param1,param2);

  // construct the generator
  boost::random_device dev;
  boost::mt19937 rng(dev);  
  boost::variate_generator<boost::mt19937, beta_distribution<> > generator(rng,mybeta);

  // generate random values
  double r = generator();
  return r;
}



/*===================================================================================*/
/*Random Number Generator from Binomial Distribution--Using BOOST C++ Library        */
/*===================================================================================*/
double gen_from_binomial_BOOST(double param1, double param2)
{
  using namespace boost::random;

  // construct the distribution                                                                         
  using boost::random::binomial_distribution;
  binomial_distribution<>bino_dis(param1,param2);

  // construct the generator                                                                            
  boost::random_device dev;
  boost::mt19937 rng(dev);
  boost::variate_generator<boost::mt19937, binomial_distribution<> > generator(rng,bino_dis);

  // generate random values                                                                             
  double r = generator();
  return r;

}

/*===================================================================================*/                 
/*Random Number Generator from Binomial Distribution--Using GNU Scientific Library   */     
/*===================================================================================*/                 
double gen_from_binomial_GSL(double param1, double param2)                                            
{
  return gsl_ran_binomial(threadvec[omp_get_thread_num()], param2, (unsigned int)param1);
}

/*===================================================================================*/
/*Random Number Generator from Dirichlet Distribution--Using GNU Scientific Library  */
/*===================================================================================*/
void gen_from_dirichlet(size_t rows, size_t cols ,double *A_ns, double *res, double *temp)
{  
  size_t i,j; 
  double seed = time(NULL);
  
  // construct the generator
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
  
  // initialize the engine with the specified seed
  gsl_rng_set (rng, seed);  
  
  // generate random values
  for (i=0; i< rows; i++){
    gsl_ran_dirichlet(rng, cols, A_ns, temp);
    for (j = 0; j< cols; j++){
      if (temp[j] < 1e-16){
	res[i*cols + j]=0;
      } 
      else{
	res[i*cols + j] = temp[j];
      }
    }
  }

  // destroy the generator
  gsl_rng_free(rng);
 
}



/*===================================================================================*/
/*Sum of the columns of a matrix in row major                                        */
/*===================================================================================*/
void colsum(size_t rows , size_t cols, double *A, double *result)
{
  size_t i, j;
  
  for (i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      result[i] = result[i] + A[i*cols+j];
    }
  }
}


/*===================================================================================*/
/* Make the columns of a row-major matrix stochastic                                 */
/*===================================================================================*/
void normprob(size_t rows , size_t cols, double *A, double *sums)
{
  
  size_t i,j;
  double temps;
  
  for(i=0; i<rows; i++){
    for(j=0; j<cols; j++){
      temps =  A[i*cols+j];
      A[i*cols+j] = temps/sums[i];
    }
  }
}
