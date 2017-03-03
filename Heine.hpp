// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp
// Date: March 3rd, 2017

// This program was made in order to calculate the Heine hypergeometric function
// 2\phi1(a,b,c,q,z), this is a q extension of Gaussian hypergeometric function

// reference
// Fredik Johansson, Computing hypergeometric functions rigorously, arXiv, 2016

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
    for(int n=1;n<=N-1;n++){
      mid=mid+qPochhammer(interval<T>(a),interval<T>(q),int (n))*qPochhammer(interval<T>(b),interval<T>(q),int (n))*pow(z,n)
	/qPochhammer(interval<T>(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=qPochhammer(interval<T>(a),interval<T>(q),int (N))*qPochhammer(interval<T>(b),interval<T>(q),int (N))*pow(z,N)
      /qPochhammer(interval<T>(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N));
    ratio=abs(z)*(1+abs(pow(q,N)*(c-a)/(1-pow(q,N))))*(1+abs(pow(q,N)*(q-b)/(1-pow(q,N))));
    if(ratio<1){
      rad=first/(1-ratio);
      res=mid+rad;
      return res;
    }
    else{
      std::cout<<"ratio is more than 1"<<std::endl;
    } 
 }
}
#endif
