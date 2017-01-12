#include<iostream>
#include<cmath>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

using namespace std;
int main(void)
{
  // verification program for the Ramanujan _1\psi_1 sum
  // one type of bilateral basic hypergeometric series
  // the program is available when |b/a|<|z| and |z|<1
  // reference
  // Encyclopedia of Mathematics and its Applications 98
  // CLASSICAL AND QUANTUM ORTHOGONAL POLYNOMIALS IN ONE VARIABLE
  // Mourad E.H.Ismail, 2005, Cambridge University Press
  // Theorem 12.3.1, page 309
cout.precision(16);
 int K=2000;

 kv::interval<double> A;
 kv::interval<double> B;
 kv::interval<double> q;
 kv::interval<double> z;
 A=0.5; // non zero value
 B=0.2;
 q=0.1; //|q|<1
 z=0.6; // non zero value
 if(abs(B/A)>=1||abs(z)>=1)
   {
     cout<<"this program is unavailable" <<endl;
   }
 else
   {
     //calculate (b/a,q,q/az,az;q)_Infinity
     // first, (b/a;q)_Infinity
     int j=1;
     kv::interval<double> q1;
     kv::interval<double> z1;
     z1=B/A;
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
     // the value of Euler function will be cancelled
     kv::interval<double> pro1=1/exp1; // pro1 OK
     // next, calculate  (q/az;q)_Infinity
     kv::interval<double> z2;
     z2=q/(A*z);
     // reuse q1
     kv::interval<double> a2;
     kv::interval<double> b2;
     a2=1/(1-z2);
     b2=1.;
     // derived from the definition of q-Pochhammer symbol
     int j2=1;
     for (int k2=1; k2<=K; k2++){
       j2 = -1*j2;
       b2=b2*(1-pow(q1,k2));
       a2=a2+j2*pow(q1,k2*(k2+1)/2)/(b2*(1-z2*pow(q1,k2)));
     }
     b2=b2*(1-pow(q1,K+1));
     double rad2=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b2*(1-z2*pow(q1,K+1)))).upper();
     kv::interval<double> exp2=a2+rad2;
     kv::interval<double> pro2=1/exp2; // pro2 OK

     // next, (az;q)_Infinity
     kv::interval<double> z3;
     z3=A*z;
     // reuse q1
     kv::interval<double> a3;
     kv::interval<double> b3;
     a3=1/(1-z3);
     b3=1.;
     // derived from the definition of q-Pochhammer symbol
     int j3=1;
     for (int k3=1; k3<=K; k3++){
       j3 = -1*j3;
       b3=b3*(1-pow(q1,k3));
       a3=a3+j3*pow(q1,k3*(k3+1)/2)/(b3*(1-z3*pow(q1,k3)));
     }
     b3=b3*(1-pow(q1,K+1));
     double rad3=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b3*(1-z3*pow(q1,K+1)))).upper();
     kv::interval<double> exp3=a3+rad3; // result must be an interval   
     kv::interval<double> pro3=1/exp3; // pro3 OK
     // calculate (b,b/az,q/a,z;q)_Infinity
     // first, (b;q)_Infinity
     kv::interval<double> z4;
     z4=B;
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
     // next, (b/az;q)_Infinity
     kv::interval<double> z5;
     z5=B/(A*z);
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
     kv::interval<double> exp5=a5+rad5; // result must be an interval
     // exp5 OK
     // next, (q/a;q)_Infinity
     kv::interval<double> z6;
     z6=q/A;
     // reuse q1
     kv::interval<double> a6;
     kv::interval<double> b6;
     a6=1/(1-z6);
     b6=1.;
     // derived from the definition of q-Pochhammer symbol
     int j6=1;
     for (int k6=1; k6<=K; k6++){
       j6 = -1*j6;
       b6=b6*(1-pow(q1,k6));
       a6=a6+j6*pow(q1,k6*(k6+1)/2)/(b6*(1-z6*pow(q1,k6)));
     }
     b6=b6*(1-pow(q1,K+1));
     double rad6=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b6*(1-z6*pow(q1,K+1)))).upper();
     kv::interval<double> exp6=a6+rad6; // result must be an interval
     // next, (z;q)_Infinity    
     // reuse q1,z
     kv::interval<double> a7;
     kv::interval<double> b7;
     a7=1/(1-z);
     b7=1.;
     // derived from the definition of q-Pochhammer symbol
     int j7=1;
     for (int k7=1; k7<=K; k7++){
       j7 = -1*j7;
       b7=b7*(1-pow(q1,k7));
       a7=a7+j7*pow(q1,k7*(k7+1)/2)/(b7*(1-z*pow(q1,k7)));
     }
     b7=b7*(1-pow(q1,K+1));
     double rad7=abs(pow(-1,K+1)*pow(q1,(K+1)*(K+2)/2)/(b7*(1-z*pow(q1,K+1)))).upper();
     kv::interval<double> exp7=a7+rad7; // result must be an interval
     kv::interval<double> res=(pro1*pro2*pro3)*(exp4*exp5*exp6*exp7);
     cout << res << endl;// final result

   }
}
