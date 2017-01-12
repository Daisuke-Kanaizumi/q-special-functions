#include<iostream>
#include<cmath>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
using namespace std;
int main(void)
{
  // verification program for the q-Ramanujan function
  // reference
  // ENCYCLOPEDIA OF MATHEMATICS AND ITS APPLICATOINS 71
  // SPECIAL FUNCTIONS
  // author: G.E.ANDREWS, RICHARD ASKEY, RANJAN ROY
  // CAMBRIDGE, 1999
  // Page 551, Exercise 39

  cout.precision(16);
  int K=2000;
  // calculate Sum[q^(n^2)*(-x)^n/(q^2;q^2)_n*(xq^2;q^2)_n,{n,0,Infinity}]
  int j=1;
  kv::interval<double> q;
  kv::interval<double> x;  
  q=0.1;
  x=0.1;
  kv::interval<double> a;
  kv::interval<double> b;
  kv::interval<double> c;
  a=1.;
  b=1.;
  c=1.;
 for (int k=1; k<=K; k++){
    j = -1*j;
    b=b*(1-pow(q,2*k));
    c=c*(1-x*pow(q,2*k));
    a = a+j*pow(q,k*k)*pow(x,k)/(b*c);
  }
  b=b*(1-pow(q,2*K+2));
  c=c*(1-x*pow(q,2*K+2));
  double rad=abs(pow(q,(K+1)*(K+1))*pow(x,K+1)/(b*c)).upper();
  // .upper() allows to output supremum
  kv::interval<double> series=a+rad; // result must be an interval
  // series OK
  // calculate (xq^2;q^2)_Infinity
  kv::interval<double> z,q1;
  q1=q*q;
  z=x*q1;
  kv::interval<double> f;
  kv::interval<double> g;
  f=1/(1-z);
  g=1.;
  int m=1;
 // derived from the definition of q-Pochhammer symbol
  for (int l=1; l<=K; l++){
    m = -1*m;
    g=g*(1-pow(q1,l));
    f=f+m*pow(q1,l*(l+1)/2)/(g*(1-z*pow(q1,l)));
  }
  g=g*(1-pow(q1,K+1));
  double rad2=abs(pow(q1,(K+1)*(K+2)/2)/(g*(1-z*pow(q1,K+1)))).upper();
  kv::interval<double> series2=f+rad2; // result must be an interval
  // calculate the value of Euler function
  kv::interval<double> a1=1-q1;
  int j1=1;
  for (int k1=1; k1<=K; k1++){
    j1 = -1*j1;
  
    a1 = a1+j1*(1-pow(q1,2*k1+1))*pow(q1,k1*(3*k1+1)/2);
  }
  double rad3=abs((1-pow(q1,2*K+3))*pow(q1,(K+1)*(3*K+4)/2)).upper();
  // .upper() allows to output supremum
  kv::interval<double> euler=a1+rad3; //final result must be an interval
  kv::interval<double> qexp=euler/series2;
  kv::interval<double> res=qexp*series;
  // this identity can be proved by Heine transformation and q-binomial theorem
  cout << res << endl;
}
