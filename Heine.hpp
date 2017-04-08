// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp
// Date: March 3rd, 2017

// This program was made in order to calculate the Heine hypergeometric function
// 2\phi1(a,b,c,q,z), this is a q extension of Gaussian hypergeometric function

// reference
// Fredrik Johansson, Computing hypergeometric functions rigorously, arXiv, 2016

#ifndef HEINE_HPP
#define HEINE_HPP
 
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/geoseries.hpp>
#include <algorithm>
#include <limits>
#include <cmath>
#include <kv/convert.hpp> // this was included to use complex numbers
#include <kv/complex.hpp>
namespace kv {
  
  //0<q<1, |z|<1
  template <class T> interval<T> Heine(const interval<T>& a, const interval<T>& b, const interval<T>& c,const interval<T>& q, const interval<T> & z) {
    int N;
    N=1000;
    interval<T> mid,rad,res,ratio,first;
    mid=1.;
    while(abs(c)>pow(1/q,N)){
      N=N+500;
      // throw std::domain_error("value of N not large enough");
    }
   if (q>=1){
     throw std::domain_error("value of q must be under 1");
   }
   if (q<=0){
     throw std::domain_error("q must be positive");
   }
    if (abs(z)>=1){
      throw std::domain_error("absolute value of z must be under 1");
    }
    for(int n=1;n<=N-1;n++){
      mid=mid+qPochhammer(interval<T>(a),interval<T>(q),int (n))*qPochhammer(interval<T>(b),interval<T>(q),int (n))*pow(z,n)
	/qPochhammer(interval<T>(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=abs(qPochhammer(interval<T>(a),interval<T>(q),int (N))*qPochhammer(interval<T>(b),interval<T>(q),int (N))*pow(z,N)
	      /qPochhammer(interval<T>(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
    ratio=abs(z)*(1+abs(pow(q,N)*(c-a)/(1-c*pow(q,N))))*(1+abs(pow(q,N)*(q-b)/(1-pow(q,N+1))));
    if(ratio<1){
      rad=(first/(1-ratio)).upper();
      res=mid+rad*interval<T>(-1.,1.);
      return res;
    }
    else{
      throw std::domain_error("ratio is more than 1");
    } 
 }
  template <class T> interval<T> q_log(const interval<T>& (q),const interval<T> &(z)){
    // verification program for q-logarithm l_q(1-z)
    // reference
    // Adrienne W. Kemp, C. David Kemp, The q-cluster distribution
    // Journal of Statistical Plannning and Inference, 2009

    // 0<q<1, 0<z<1
   if (q>=1){
     throw std::domain_error("value of q must be under 1");
   }
   if (q<=0){
     throw std::domain_error("q must be positive");
   }
   if (z>=1){
     throw std::domain_error("value of z must be under 1");
   }
   if (z<=0){
     throw std::domain_error("z must be positive");
   }
    interval<T>res;
    res=-z*Heine(interval<T>(q),interval<T>(q),interval<T>(q*q),interval<T>(q),interval<T>(z));
    return res;
  }
 template <class T> complex<interval<T> >Heine(const complex<interval<T> >& a, const complex<interval<T> >& b, const complex<interval<T> >& c,const interval<T>& q, const complex<interval<T> >& z) {
    int N;
    N=1000;
    complex<interval<T> > mid,res;
    interval<T>ratio,first;
    T rad;
    mid=1.;

    while(abs(c)>pow(1/q,N)){
      N=N+500;
      // throw std::domain_error("value of N not large enough");
    }
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    if (abs(z)>=1){
      throw std::domain_error("absolute value of z must be under 1");
    }
    for(int n=1;n<=N-1;n++){
      mid=mid+qPochhammer(complex<interval<T> >(a),interval<T>(q) ,int (n))*qPochhammer(complex<interval<T> >(b),interval<T>(q),int (n))*pow(z,n)
	/qPochhammer(complex<interval<T> >(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=abs(qPochhammer(complex<interval<T> >(a),interval<T>(q),int (N))*qPochhammer(complex<interval<T> >(b),interval<T>(q),int (N))*pow(z,N)
	      /qPochhammer(complex<interval<T> >(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
    ratio=abs(z)*(1+abs(pow(q,N)*(c-a)/(1-c*pow(q,N))))*(1+abs(pow(q,N)*(q-b)/(1-pow(q,N+1))));
    if(abs(ratio)<1){
      rad=(first/(1-ratio)).upper();
      res=complex_nbd(mid,rad);
      return res;
    }
    else{
      throw std::domain_error("ratio is more than 1");
    } 
}
  template <class T> complex<interval<T> >complex_nbd(const complex<interval<T> >& a,const T & rad){
    interval<T> rc,ic,rcc,icc;
    complex<interval<T> >res;
    rc=a.real();
    ic=a.imag();
    rcc=rc+rad*interval<T>(-1.,1.);
    icc=ic+rad*interval<T>(-1.,1.);
    res=complex<interval<T> >(rcc,icc);
    return res;
  }
}
#endif
