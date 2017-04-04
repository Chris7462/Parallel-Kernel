#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "func.h"

#define REPS 100

using namespace std;

int main(int argc, char* argv[]){
  if (argc != 4) {
    fprintf(stderr, "Usage: kde num_threads sample_size kernel\n");
    exit(1);
  }
  const int m = atoi(argv[1]);	// number of threads
  const int n = atoi(argv[2]);  // number of sample size
  string ker = argv[3];		// specific which kernel

  if ( !(ker.compare("rect") == 0 || ker.compare("epan") == 0 || ker.compare("gauss") == 0 || ker.compare("quar") == 0) ){
    fprintf(stderr, "Wrong kernel input!\nUse: \"rect\", \"epan\", \"gauss\", or \"quar\". \n");
    exit(1);
  }

  vector<double> X(n);    // Array for saving generated normal random variables
  vector<double> eps(n);  // Array for saving generated measurment error
  vector<double> x(n);    // Grid for density estimation, equally spacing x

  /* mean and stdev for normal distribution */
  double mu = 2, sigma = 0.5;

  /* add random noise to the data */
  X = rnorm(n, mu, sigma);  // generate normal random variable
  eps = runif(n,-0.5,0.5);  // small random noise
  X = X+eps;

	/* get linsapce of X */
	x = linspace(X, n);
	//for_each(x.begin(), x.end(), [&](double &n){ cout << n << endl; });

  int i, j, rep;
  double* pdf = new double[n];
  double sum_ker;
  double t_start, t_end;
  double h = pow(n,-0.2);

  double acumtime = 0.0;
  for (rep = 0; rep < REPS; rep++){
    t_start = omp_get_wtime();
#pragma omp parallel num_threads(m) default(none) shared(x,X,i,pdf,ker,h) private(sum_ker,j)
    {
#pragma omp for
      for (i = 0; i < n; ++i ){
        sum_ker = 0;
        for (j = 0; j < n; ++j ){
          sum_ker = sum_ker + kernel(x[i],X[j], h, ker);
        }
        pdf[i] = sum_ker/n;
      }
    }
    t_end = omp_get_wtime();
    acumtime = t_end-t_start;
  }

  printf("avg %d run: m = %d\tn = %d\tt = %lf\n", REPS, m, n, acumtime/REPS);

  delete [] pdf;

  return 0;
}
