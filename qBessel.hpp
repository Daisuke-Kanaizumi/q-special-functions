// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp

// verification program for q-Bessel functions
// March 9th, 2017

#ifndef QBESSEL_HPP
#define QBESSEL_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

namespace kv{
template <class T> interval<T>Jackson1(const interval<T>& z,const interval<T>& nu,const interval<T>& q){
// verification program for Jackson`s 1st q-Bessel function
interval<T>res,a,b,c,series;
int j,K;
T rad;
if(abs(z)>=2){
throw std::domain_error("Jackson`s 1st q-Bessel function is not defined");
}
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 if(z<0){
   throw std::domain_error("z must be positive");
 }
 if (pow(q,nu)<1){
j=1;
K=500;
a=1.;
b=1.;
c=1.;
// Leibniz criterion
 for (int k=1; k<=K; k++){
    j = -1*j;
    b=b*(1-pow(q,k));
    c=c*(1-pow(q,k+nu));
    a = a+j*pow(z/2,2*k)/(b*c);
  }
b=b*(1-pow(q,K+1));
c=c*(1-pow(q,K+nu+1));

rad=abs(pow(z/2,2*K+2)/(b*c)).upper();
series=a+rad*interval<T>(-1.,1.);
res=pow(z/2,nu)*series/Karpelevich(interval<T>(pow(q,nu+1)),interval<T>(q));
}
 else{
res=pow(z/2,nu)*infinite_qPochhammer(interval<T>(pow(q,nu+1)),interval<T>(q))*
  Heine(interval<T>(0),interval<T>(0),interval<T>(pow(q,nu+1)),interval<T>(q),interval<T>(-z*z/4))/Euler(interval<T>(q));
}
 return res;
}
template <class T> complex<interval<T> >Jackson1(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
complex<interval<T> >res;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 res=pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pow(q,nu+1)),interval<T>(q))*
  Heine(complex<interval<T> >(0),complex<interval<T> >(0),complex<interval<T> >(pow(q,nu+1)),interval<T>(q),complex<interval<T> >(-z*z/4))/Euler(interval<T>(q));
 return res;       
}

template <class T> complex<interval<T> >modified_qBessel1(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  // verification program for 1st modified q-Bessel
complex<interval<T> >res;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 res=pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pow(q,nu+1)),interval<T>(q))*
  Heine(complex<interval<T> >(0),complex<interval<T> >(0),complex<interval<T> >(pow(q,nu+1)),interval<T>(q),complex<interval<T> >(z*z/4))/Euler(interval<T>(q));
 return res;       
}


template <class T> interval<T>Jackson2(const interval<T>& z,const interval<T>& nu,const interval<T>& q){
// verification program for Jackson`s 2nd q-Bessel function
interval<T>res,a,b,c,series;
int j,K;
T rad;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 if (pow(q,nu)<1 && abs(z)>pow(q,nu)&&z<500){
   // this method will not converge when z>500
j=1;
K=500;
a=1.;
b=1.;
c=1.;
// Leibniz criterion
 for (int k=1; k<=K; k++){
    j = -1*j;
    b=b*(1-pow(q,k));
    c=c*(1-pow(q,k+nu));
    a = a+j*pow(q,k*(k-1))*pow(pow(q,nu+1)*z*z/4,k)/(b*c);
  }
  b=b*(1-pow(q,K+1));
  c=c*(1-pow(q,K+nu+1));
 
 rad=abs(pow(q,K*(K+1))*pow(pow(q,nu+1)*z*z/4,K+1)/(b*c)).upper();
 series=a+rad*interval<T>(-1.,1.);
 res=pow(z/2,nu)*series/Karpelevich(interval<T>(pow(q,nu+1)),interval<T>(q));
}
 else{
   res=pow(z/2,nu)*infinite_qPochhammer(interval<T>(pow(q,nu+1)),interval<T>(q))*
     _0phi_1(interval<T>(pow(q,nu+1)),interval<T>(q),interval<T>(-z*z*pow(q,nu+1)/4))/Euler(interval<T>(q));
 }
  return res;
}
template <class T> complex<interval<T> >Jackson2(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  complex<interval<T> >res,pq;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 pq=pow(q,nu+1);
 res=pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pq),interval<T>(q))*
   _0phi_1(complex<interval<T> >(pq),interval<T>(q),complex<interval<T> >(-z*z*pq/4))/Euler(interval<T>(q));
 return res;       
}
template <class T> complex<interval<T> >modified_qBessel2(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  //verification program for 2nd modified q-Bessel function
  complex<interval<T> >res,pq;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 pq=pow(q,nu+1);
 res=pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pq),interval<T>(q))*
   _0phi_1(complex<interval<T> >(pq),interval<T>(q),complex<interval<T> >(z*z*pq/4))/Euler(interval<T>(q));
 return res;       
}

template <class T> complex<interval<T> >Hahn_Exton(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  // verification program for Hahn-Exton q-Bessel function
  complex<interval<T> >res,pq;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 pq=pow(q,nu+1);
 res=pow(z,nu)*infinite_qPochhammer(complex<interval<T> >(pq),interval<T>(q))*
   _1phi_1(complex<interval<T> >(0),complex<interval<T> >(pq),interval<T>(q),complex<interval<T> >(z*z*q))/Euler(interval<T>(q));
 return res;       
}
template <class T> complex<interval<T> >little(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  // verification program for little q-Bessel function
  // reference
  // Koornwinder and Swarttouw, On q-analogues of the Fourier and Hankel transforms, 1992
  // Bouzeffour, New Addition Formula for the Little q-Bessel Functions, arXiv, 2013
  complex<interval<T> >res,pq;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 pq=pow(q,nu+1);
 res=pow(z,nu)*infinite_qPochhammer(complex<interval<T> >(pq),interval<T>(q))*
   _1phi_1(complex<interval<T> >(0),complex<interval<T> >(pq),interval<T>(q),complex<interval<T> >(z))/Euler(interval<T>(q));
 return res;       
}

}

#endif
