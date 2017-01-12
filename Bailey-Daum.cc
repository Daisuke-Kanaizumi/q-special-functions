#include<iostream>
#include<cmath>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

using namespace std;
int main(void)
{
  // verification program for the Bailey Daum sum
  // _2\phi_1(a,b;aq/b;q,-q/b)
  // reference
  // Encyclopedia of Mathematics and its Applications 98
  // CLASSICAL AND QUANTUM ORTHOGONAL POLYNOMIALS IN ONE VARIABLE
  // Mourad E.H.Ismail, 2005, Cambridge University Press
  // Theorem 12.5.3, page 315
cout.precision(16);
 int K=2000;

 kv::interval<double> A;
 kv::interval<double> B;
 kv::interval<double> q;
 A=0.5;
 B=0.2;
 q=0.1;

 // restriction:|q|<|b|
    
 // first, (-q;q)_Infinity
 int j=1;
 kv::interval<double> q1;
 kv::interval<double> z1;
 z1=-q;
 q1=q;
 kv::interval<double> a;
 kv::interval<double> b;
 a=1/(1-z1);
 b=1.;
 // derived from the definition of q-Pochhammer symbol
 for (int k=1; k<=K; k++){
   j = -1*j;
   b=b*(1-pow(q1,k));
   a=a+j*pow(q1,k*(k+1)/2)/(b*(1-z1*pow(q1,k)));
 }
 b=b*(1-pow(q1,K+1));
 double rad1=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b*(1-z1*pow(q1,K+1)))).upper();
 kv::interval<double> exp1=a+rad1; // result must be an interval
 // the value of Euler function (q,q)_Infinity will be cancelled once
 kv::interval<double> pro1=1/exp1; // pro1 OK
 // calculate  (aq,(aq^2)/b^2;q^2)_Infinity
 // first, calculate (aq,q^2)_Infinity 
 kv::interval<double> z2;
 kv::interval<double> q2;
 z2=A*q;
 q2=q*q;
 kv::interval<double> a2;
 kv::interval<double> b2;
 a2=1/(1-z2);
 b2=1.;
 // derived from the definition of q-Pochhammer symbol
 int j2=1;
 for (int k2=1; k2<=K; k2++){
   j2 = -1*j2;
   b2=b2*(1-pow(q2,k2));
   a2=a2+j2*pow(q2,k2*(k2+1)/2)/(b2*(1-z2*pow(q2,k2)));
 }
 b2=b2*(1-pow(q2,K+1));
 double rad2=abs(pow(-1,K+1)*pow(q2,(K+1)*(K+2)/2)/(b2*(1-z2*pow(q2,K+1)))).upper();
 kv::interval<double> exp2=a2+rad2;
 // calculate the value of Euler function
 // (q^2,q^2)_Infinity
 kv::interval<double> c1=1-q2;
 int m=1;
 for (int l=1; l<=K; l++){
   m = -1*m;
   c1 = c1+m*(1-pow(q2,2*l+1))*pow(q2,l*(3*l+1)/2);
 }
 double eulerrad1=abs(pow(-1,K+1)*(1-pow(q2,2*K+3))*pow(q2,(K+1)*(3*K+4)/2)).upper();
 kv::interval<double> q2euler=c1+eulerrad1;
 
 kv::interval<double> pro2=q2euler/exp2; // pro2 OK

 // next, ((a*q^2)/b^2;q^2)_Infinity
 kv::interval<double> z3;
 z3=A*q*q/(B*B);
 // reuse q2
 kv::interval<double> a3;
 kv::interval<double> b3;
 a3=1/(1-z3);
 b3=1.;
 // derived from the definition of q-Pochhammer symbol
 int j3=1;
 for (int k3=1; k3<=K; k3++){
   j3 = -1*j3;
   b3=b3*(1-pow(q2,k3));
   a3=a3+j3*pow(q2,k3*(k3+1)/2)/(b3*(1-z3*pow(q2,k3)));
 }
 b3=b3*(1-pow(q2,K+1));
 double rad3=abs(pow(-1,K+1)*pow(q2,(K+1)*(K+2)/2)/(b3*(1-z3*pow(q2,K+1)))).upper();
 kv::interval<double> exp3=a3+rad3; // result must be an interval   
 kv::interval<double> pro3=q2euler/exp3; // pro3 OK
 // value q2euler is reused
 // calculate (-q/b,aq/b;q)_Infinity
 // first, (-q/b;q)_Infinity
 kv::interval<double> z4;
 z4=-q/B;
 // reuse q1
 kv::interval<double> a4;
 kv::interval<double> b4;
 a4=1/(1-z4);
 b4=1.;
 // derived from the definition of q-Pochhammer symbol
 int j4=1;
 for (int k4=1; k4<=K; k4++){
   j4 = -1*j4;
   b4=b4*(1-pow(q1,k4));
   a4=a4+j4*pow(q1,k4*(k4+1)/2)/(b4*(1-z4*pow(q1,k4)));
 }
 b4=b4*(1-pow(q1,K+1));
 double rad4=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b4*(1-z4*pow(q1,K+1)))).upper();
 kv::interval<double> exp4=a4+rad4; // result must be an interval   
 // exp4 OK
 // next, (aq/b;q)_Infinity
 kv::interval<double> z5;
 z5=A*q/B;
 // reuse q1
 kv::interval<double> a5;
 kv::interval<double> b5;
 a5=1/(1-z5);
 b5=1.;
 // derived from the definition of q-Pochhammer symbol
 int j5=1;
 for (int k5=1; k5<=K; k5++){
   j5 = -1*j5;
   b5=b5*(1-pow(q1,k5));
   a5=a5+j5*pow(q1,k5*(k5+1)/2)/(b5*(1-z5*pow(q1,k5)));
 }
 b5=b5*(1-pow(q1,K+1));
 double rad5=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b5*(1-z5*pow(q1,K+1)))).upper();
 kv::interval<double> series=a5+rad5; // result must be an interval
 // exp5 OK
// calculate the value of Euler function
 kv::interval<double> c2=1-q1;
 int m1=1;
 for (int l1=1; l1<=K; l1++){
   m1 = -1*m1;
   c2 = c2+m1*(1-pow(q1,2*l1+1))*pow(q1,l1*(3*l1+1)/2);
  }
 double eulerrad2=abs(pow(-1,K+1)*(1-pow(q1,2*K+3))*pow(q1,(K+1)*(3*K+4)/2)).upper();
 kv::interval<double> euler=c2+eulerrad2;
 kv::interval<double> exp5=series/euler;
 kv::interval<double> res=(pro1*pro2*pro3)*(exp4*exp5);
 cout << res << endl;// final result  
}
