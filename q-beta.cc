#include<iostream>
#include<cmath>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/beta.hpp>
using namespace std;
int main(void)
{
  // verification program for the q-beta function B_q(x,y) 
  // theorem by Fridrikh Israilevich Karpelevich was used
  // reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
  // Olshanetsky, Rogov, 1995
  cout.precision(16);
  int K=1500;
  // calculate gamma_q(x)
  int j=1;
  kv::interval<double> q; // q is any real number
  kv::interval<double> x;
  kv::interval<double> y;
  q=1.5;
  x=20.;
  y=10.;
  
  if (abs(q)<1){
    kv::interval<double> z;
    z=x;
    kv::interval<double> a;
    kv::interval<double> b;
    a=1/(1-pow(q,z));
    b=1.;
    // derived from the definition of q-Pochhammer symbol
    for (int k=1; k<=K; k++){
      j = -1*j;
      b=b*(1-pow(q,k));
      a=a+j*pow(q,k*(k+1)/2)/(b*(1-pow(q,z+k)));
    }
    b=b*(1-pow(q,K+1));
    double rad=abs(pow(-1,K+1)*pow(q,(K+1)*(K+2)/2)/(b*(1-pow(q,K+1+z)))).upper();
    kv::interval<double> series=a+rad; // result must be an interval
    kv::interval<double> c;
    c=pow(1-q,1-z);
    kv::interval<double> res1;
    res1=c*series;
    // calculate gamma_q(y)
 
    kv::interval<double> z2;

    z2=y;
    kv::interval<double> a2;
    kv::interval<double> b2;
    a2=1/(1-pow(q,z2));
    b2=1.;
    int j2;
    j2=1;
    // derived from the definition of q-Pochhammer symbol
    for (int k2=1; k2<=K; k2++){
      j2 = -1*j2;
      b2=b2*(1-pow(q,k2));
      a2=a2+j2*pow(q,k2*(k2+1)/2)/(b2*(1-pow(q,z2+k2)));
    }
    b2=b2*(1-pow(q,K+1));
    double rad2=abs(pow(-1,K+1)*pow(q,(K+1)*(K+2)/2)/(b2*(1-pow(q,K+1+z2)))).upper();
    kv::interval<double> series2=a2+rad2; // result must be an interval
    kv::interval<double> c2;
    c2=pow(1-q,1-z2);
    kv::interval<double> res2;
    res2=c2*series2;
    // calculate gamma_(x+y)
    kv::interval<double> z3;
    z3=x+y;
    kv::interval<double> a3;
    kv::interval<double> b3;
    a3=1/(1-pow(q,z3));
    b3=1.;
    int j3;
    j3=1;
    // derived from the definition of q-Pochhammer symbol
    for (int k3=1; k3<=K; k3++){
      j3 = -1*j3;
      b3=b3*(1-pow(q,k3));
      a3=a3+j3*pow(q,k3*(k3+1)/2)/(b3*(1-pow(q,z3+k3)));
    }
    b3=b3*(1-pow(q,K+1));
    double rad3=abs(pow(-1,K+1)*pow(q,(K+1)*(K+2)/2)/(b3*(1-pow(q,K+1+z3)))).upper();
    kv::interval<double> series3=a3+rad3; // result must be an interval
    kv::interval<double> c3;
    c3=pow(1-q,1-z3);
    kv::interval<double> res3;
    res3=c3*series3;
    // calculate the value of q-beta function
    kv::interval<double> res;
    res=res1*res2/res3;
    cout << res << endl;
  }
  else if(q==1){
    cout<<kv::beta(kv::interval<double>(x),kv::interval<double>(y))<<endl;
  }
  else if(abs(q)>1){
    kv::interval<double> z4;
    z4=x;
    kv::interval<double> q1;
    q1=1/q;
    kv::interval<double> a4;
    kv::interval<double> b4;
    a4=1/(1-pow(q1,z4));
    b4=1.;
    int j4=1;
    // derived from the definition of q-Pochhammer symbol
    for (int k4=1; k4<=K; k4++){
      j4 = -1*j4;
      b4=b4*(1-pow(q1,k4));
      a4=a4+j4*pow(q1,k4*(k4+1)/2)/(b4*(1-pow(q1,z4+k4)));
    }
    b4=b4*(1-pow(q1,K+1));
    double rad4=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b4*(1-pow(q1,K+1+z4)))).upper();
    kv::interval<double> series4=a4+rad4; // result must be an interval
    kv::interval<double> c4;
    c4=pow(q-1,1-z4)*pow(q,z4*(z4-1)/2);
    kv::interval<double> res4;
    res4=c4*series4;
    // calculate gamma_q(y)
 
    kv::interval<double> z5;

    z5=y;
    kv::interval<double> a5;
    kv::interval<double> b5;
    a5=1/(1-pow(q1,z5));
    b5=1.;
    int j5;
    j5=1;
    // derived from the definition of q-Pochhammer symbol
    for (int k5=1; k5<=K; k5++){
      j5 = -1*j5;
      b5=b5*(1-pow(q1,k5));
      a5=a5+j5*pow(q1,k5*(k5+1)/2)/(b5*(1-pow(q1,z5+k5)));
    }
    b5=b5*(1-pow(q1,K+1));
    double rad5=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b5*(1-pow(q1,K+1+z5)))).upper();
    kv::interval<double> series5=a5+rad5; // result must be an interval
    kv::interval<double> c5;
    c5=pow(q-1,1-z5)*pow(q,z5*(z5-1)/2);
    kv::interval<double> res5;
    res5=c5*series5;
    // calculate gamma_(x+y)
    kv::interval<double> z6;
    z6=x+y;
    kv::interval<double> a6;
    kv::interval<double> b6;
    a6=1/(1-pow(q1,z6));
    b6=1.;
    int j6;
    j6=1;
    // derived from the definition of q-Pochhammer symbol
    for (int k6=1; k6<=K; k6++){
      j6 = -1*j6;
      b6=b6*(1-pow(q1,k6));
      a6=a6+j6*pow(q1,k6*(k6+1)/2)/(b6*(1-pow(q1,z6+k6)));
    }
    b6=b6*(1-pow(q1,K+1));
    double rad6=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b6*(1-pow(q1,K+1+z6)))).upper();
    kv::interval<double> series6=a6+rad6; // result must be an interval
    kv::interval<double> c6;
    c6=pow(q-1,1-z6)*pow(q,z6*(z6-1)/2);
    kv::interval<double> res6;
    res6=c6*series6;
    // calculate the value of q-beta function
    kv::interval<double> ans;
    ans=res4*res5/res6;
    cout << ans << endl;
  }
}
 // definition of the q-gamma function is given in the reference below
 // Basic Hypergeometric Series Cambridge University Press (Gasper, Rahman, 1990)
