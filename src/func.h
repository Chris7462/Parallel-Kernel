#ifndef _FUNC_H_
#define _FUNC_H

#include <algorithm>
#include <random>
#include <string>
using namespace std;

vector<double> rnorm(const int, const double, const double);
vector<double> runif(const int, const double, const double);
vector<double> linspace(const vector<double>&, const int);
double kernel(double, double, double, string);
vector<double> operator+(const vector<double>&, const vector<double>&);

/* generate nromal random variable */
vector<double> rnorm(const int n, const double mu, const double sigma){
	vector<double> X(n);          /* Array for saving generated random variables */
  vector<double>::iterator it;

  default_random_engine generator;
  normal_distribution<double> distribution(mu,sigma);

  /* generate random variables */
  for ( it = X.begin(); it != X.end(); ++it ){
    *it = distribution(generator);
  }
  return(X);
}

/* generate uniform random variable */
vector<double> runif(const int n, const double a, const double b){
  vector<double> X(n);          /* Array for saving generated random variables */
  vector<double>::iterator it;

  default_random_engine generator;
  uniform_real_distribution<double> distribution(a,b);

  /* generate random variables */
  for ( it = X.begin(); it != X.end(); ++it ){
    *it = distribution(generator);
  }
  return(X);
}

/* find the grid points of X */
vector<double> linspace(const vector<double>& X, const int m){
	/* calculate the mean of X */
  double X_mu = accumulate(X.begin(), X.end(), 0.0)/X.size();

	/* calculate the stdev of X */
  double accum = 0.0;
  for_each(X.begin(), X.end(), [&](const double d){ accum += (d-X_mu)*(d-X_mu); });
  double X_stdev = sqrt(accum/(X.size()-1));

  /* calculate the min and max of X */
  double x_min, x_max, X_min, X_max;
  X_min = *min_element(X.begin(),X.end());
  X_max = *max_element(X.begin(),X.end());
  x_min = X_min-1.1*X_stdev;
  x_max = X_max+1.1*X_stdev;

  /* get grid points of X */
  vector<double> array;
  double step = (x_max-x_min)/(m-1);
  for ( int i = 0; i < m; ++i	){
    array.push_back(x_min+i*step);
  }
  return array;
}

double kernel(double x, double X, double h, string str){
  double u, val;
  u = (x-X)/h;

  if ( str.compare("rect") == 0 ){
    val = 0.5*static_cast<double>(abs(u)<=1);
  } else if ( str.compare("epan") == 0 ){
    val = 0.75*static_cast<double>(abs(u)<=1)*(1-u*u);
  } else if ( str.compare("gauss") == 0 ){
    val = 1/sqrt(2*M_PI)*exp(-0.5*u*u);
  } else if ( str.compare("quar") == 0 ){
    val = 0.9375*static_cast<double>(abs(u)<=1)*pow((1-u*u),2);
  }
  return(val);
}

vector<double> operator+(const vector<double>& a, const vector<double>& b){
  int n = a.size();
  vector<double> y;
  vector<double>::iterator it;

  for ( int i = 0 ; i < n; ++i ){
    y.push_back(a[i]+b[i]);
  }
  return y;
}

#endif
