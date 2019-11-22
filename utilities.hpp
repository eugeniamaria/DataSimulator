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
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

using namespace std;
extern gsl_rng** threadvec;

/* types of user command-line input */
typedef enum {
  INT,
  DOUBLE,
  STR,
  NA
} ARG_TYPE;


int findarg(const char *argname, ARG_TYPE type, void *val, int argc, char **argv); 
void InitializeThread();
void printBlock(double *dbuffer,size_t bufferlength);
size_t ReadFiles(vector<double> & values, ifstream & filename);
void ComputeBetas(double *B_kl, size_t row, double fst, double freq, size_t K);
double gen_from_beta(double param1, double param2);
void ComputeThetasA(double *TB_kn, double alpha, double gamma, double MN, size_t NS, size_t K, size_t I);
void gen_from_dirichlet(size_t rows, size_t cols ,double *A_ns, double *res, double *temp);
void colsum(size_t rows, size_t cols, double *A, double *result);
void normprob(size_t rows, size_t cols, double *A, double *sums);
double gen_from_binomial_GSL(double param1, double param2);
void ComputeThetasB(size_t n, size_t k, size_t s, double a, double b, double *Q);
double gen_from_binomial_BOOST(double param1, double param2);
