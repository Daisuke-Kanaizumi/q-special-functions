// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University

// Verified computation of infinite q-Pochhammer symbol by the Sokal's algorithm


// Reference
// Sokal, Alan. D. (2002). arXiv preprint math/0212035.

#ifndef SOKAL_HPP
#define SOKAL_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/complex.hpp>
#include <kv/constants.hpp>
#include <kv/Heine.hpp>
#include <kv/Pochhammer.hpp>
#include <cmath>

namespace kv{
  template <class T> complex<interval<T> >Sokal_qp(const complex<interval<T> >& z,const complex<interval<T> >& q, int N=50){
    if (abs(q)>=1){
      throw std::domain_error("value of q must be under 1");
    }
    complex<interval<T> >res,sum,zz;
    sum=0.;
    zz=1.;
    interval<T> gamma,pi;
    T rad;
    pi=kv::constants<interval<T> >::pi();
    gamma=1.;
    while(abs(q)<=exp(-gamma)){
      gamma=gamma+1;
    }
    while(abs(z)>=exp((N+1)*gamma)){
      N=N+10;
    }
    while(abs(z*pow(q,N+1))>=1){
      N=N+10;
    }
   
    for(int n=0;n<=N-1;n++){
      sum=sum+zz*pow(q,n*(n-1)*0.5)/qPochhammer(q,q,n);
   
      zz*=-z;
    }
    rad=(pow(abs(z),N)*exp(pi*pi/6./gamma-N*(N-1)*gamma*0.5)/(1-abs(z)*exp(-(N+1))*gamma)).upper();
    res=complex_nbd(sum,rad);
    return res;
  }
  
  template <class T> complex<interval<T> >Sokal_qp(const complex<interval<T> >& z,const interval<T>& q, int N=100){
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("value of q must be positive");
    }
    complex<interval<T> >res,sum,zz;
    sum=0.;
    zz=1.;
    interval<T> gamma, pi;
    pi=kv::constants<interval<T> >::pi();
    T rad;
    gamma=1.;
    while(abs(q)<=exp(-gamma)){
      gamma=gamma+1;
    }
    while(abs(z)>=exp((N+1)*gamma)){
      N=N+10;
    }
    while(abs(z)*pow(q,N+1)>=1){
      N=N+10;
    }
 
    for(int n=0;n<=N-1;n++){
      sum=sum+zz*pow(q,n*(n-1)*0.5)/qPochhammer(q,q,n);
      
      zz=-z*zz;
    }
    rad=(pow(abs(z),N)*exp(pi*pi/6./gamma-N*(N-1)*gamma*0.5)/(1-abs(z)*exp(-(N+1))*gamma)).upper();
    res=complex_nbd(sum,rad);
    return res;
  }
  
}
#endif

