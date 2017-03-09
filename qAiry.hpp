// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp
// Date: March 4th, 2017

// verification program for q-Airy functions

#ifndef QAIRY_HPP
#define QAIRY_HPP
 
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/geoseries.hpp>
#include <algorithm>
#include <limits>
#include <cmath>
#include <kv/convert.hpp> // this was included to use complex numbers
#include <kv/complex.hpp>
namespace kv {

  // calculate basic hypergeometric series 1phi1
  
  // reference
  // Fredrik Johansson, Computing hypergeometric functions rigorously, arXiv, 2016

template <class T> interval<T> _1phi_1(const interval<T>& a, const interval<T>& c,const interval<T>& q, const interval<T> & z) {
  int N;
    N=1000;
    interval<T> mid,rad,res,ratio,first;
    mid=1.;
    for(int n=1;n<=N-1;n++){
      mid=mid+qPochhammer(interval<T>(a),interval<T>(q),int (n))*pow(-z,n)*pow(q,n*(n-1)/2)
	/qPochhammer(interval<T>(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=qPochhammer(interval<T>(a),interval<T>(q),int (N))*pow(-z,N)*pow(q,N*(N-1)/2)
      /qPochhammer(interval<T>(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N));
    ratio=abs(z)*abs(pow(q,N))*(1+abs(pow(q,N)*(c-a)/(1-c*pow(q,N))))*abs(1/(1-pow(q,N+1)));
    if(ratio<1){
      rad=first/(1-ratio);
      res=mid+rad;
      return res;
    }
    else{
      std::cout<<"ratio is more than 1"<<std::endl;
    }
 }
  template <class T> complex<interval<T> >_1phi_1(const complex<interval<T> >& a, const complex<interval<T> >& c,const interval<T>& q, const complex<interval<T> >& z) {
     int N;
    N=1000;
    complex<interval<T> > mid,rad,res,ratio,first;
    mid=1.;
    for(int n=1;n<=N-1;n++){
      mid=mid+qPochhammer(complex<interval<T> >(a),interval<T>(q) ,int (n))*pow(-z,n)*pow(q,n*(n-1)/2)
	/qPochhammer(complex<interval<T> >(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=qPochhammer(complex<interval<T> >(a),interval<T>(q),int (N))*pow(-z,N)*pow(q,N*(N-1)/2)
      /qPochhammer(complex<interval<T> >(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N));
    ratio=abs(z)*pow(q,N)*(1+abs(pow(q,N)*(c-a)/(1-c*pow(q,N))))*abs(1/(1-pow(q,N+1)));
    if(abs(ratio)<1){
      rad=first/(1-ratio);
      res=mid+rad;
      return res;
    }
    else{
      std::cout<<"ratio is more than 1"<<std::endl;
      }}

  // calculate basic hypergeometric series 0phi1
  
  // reference
  // Fredrik Johansson, Computing hypergeometric functions rigorously, arXiv, 2016


template <class T> interval<T> _0phi_1(const interval<T>& c,const interval<T>& q, const interval<T> & z) {
    int N;
    N=1000;
    interval<T> mid,rad,res,ratio,first;
    mid=1.;
    for(int n=1;n<=N-1;n++){
      mid=mid+pow(z,n)*pow(q,n*(n-1))
	/qPochhammer(interval<T>(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=pow(z,N)*pow(q,N*(N-1))
      /qPochhammer(interval<T>(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N));
    ratio=abs(z)*pow(q,2*N)*abs(1/(1-c*pow(q,N)))*abs(1/(1-pow(q,N+1)));
    if(ratio<1){
      rad=first/(1-ratio);
      res=mid+rad;
      return res;
    }
    else{
      std::cout<<"ratio is more than 1"<<std::endl;
    }
 }
  template <class T> complex<interval<T> >_0phi_1(const complex<interval<T> >& c,const interval<T>& q, const complex<interval<T> >& z) {
     int N;
    N=1000;
    complex<interval<T> > mid,rad,res,ratio,first;
    mid=1.;
    for(int n=1;n<=N-1;n++){
      mid=mid+pow(z,n)*pow(q,n*(n-1))
	/qPochhammer(complex<interval<T> >(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=pow(z,N)*pow(q,N*(N-1))
      /qPochhammer(complex<interval<T> >(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N));
    ratio=abs(z)*pow(q,2*N)*abs(1/(1-c*pow(q,N)))*abs(1/(1-pow(q,N+1)));
    if(abs(ratio)<1){
      rad=first/(1-ratio);
      res=mid+rad;
      return res;
    }
    else{
      std::cout<<"ratio is more than 1"<<std::endl;
      } 
 }
 template <class T> complex<interval<T> >HKW_qAiry(const interval<T>& q, const complex<interval<T> >& z) {
   // calculate Hamamoto-Kajiwara-Witte q-Airy function
   // reference
   // Hamamoto, Kajiwara, Witte, Hypergeometric solutions to the q-Painlev\'e equation of type $(A_1+A^{\prime}_1)~{(1)}$
   // International Mathematics Research Notices, 2006
   complex<interval<T> >res,qq;
   qq=q;
   res=_1phi_1(complex<interval<T> >(0),complex<interval<T> >(-qq),interval<T>(q),complex<interval<T> >(-z));
   return res;
 }
  template <class T> complex<interval<T> >KMNOY_qAiry(const interval<T>& q, const complex<interval<T> >& z) {
    // calculate Kajiwara-Masuda-Noumi-Ohta-Yamada q-Airy function
    // reference
    // Kajiwara, Masuda, Noumi, Ohta, Yamada, Hypergeometric Solutions to the q-Painlev\'e Equations
    // International Mathematics Research Notices, 2004
    complex<interval<T> >res,qq;
    qq=q;
    res=_1phi_1(complex<interval<T> >(0),complex<interval<T> >(-qq),interval<T>(q),complex<interval<T> >(-z*pow(q,0.5)));
    return res;
 }
template <class T> complex<interval<T> >Ramanujan_qAiry(const interval<T>& q, const complex<interval<T> >& z) {
  // calculate Ramanujan`s q-Airy function
  complex<interval<T> >res,qq;
   qq=q;
   res=_0phi_1(complex<interval<T> >(0),interval<T>(q),complex<interval<T> >(-qq*z));
   return res;
 }
}
#endif
