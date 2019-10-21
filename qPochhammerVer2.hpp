// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University

// Verified computation of infinite q-Pochhammer symbol
// This implementation is the second version
// First Version:Pochhammer.hpp

// Reference
// Zhang, R. (2008). On asymptotics of q-Gamma functions. Journal of Mathematical Analysis and Applications, 339(2), 1313-1321.
// Lemma 1.1
#ifndef QPOCHHAMMERVER2_HPP
#define QPOCHHAMMERVER2_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/complex.hpp>
#include <kv/constants.hpp>
#include <kv/Heine.hpp>
#include <kv/Pochhammer.hpp>
#include <kv/Gatteschi.hpp>
#include <cmath>

namespace kv{
  template <class T> complex<interval<T> >qPVer2(const complex<interval<T> >& z,const interval<T>& q, int n=100, int K=100){
    complex<interval<T> >res,r,mid;
    mid=0.;
    T rad;
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("value of q must be positive");
    }
    if(q<0.9){
      while(abs(z)*pow(q,n)/(1-q)>=0.5){
	n=n+50;
      }
      rad=(2*pow(q,K*(K-1)*0.5)*pow(abs(z)*pow(q,n),K)
	   /qPochhammer(interval<T>(q),interval<T>(q),int (K))).upper();
      for(int k=0;k<=K-1;k++){
	mid=mid+pow(q,k*(k-1)*0.5)*pow(-z*pow(q,n),k)
	  /qPochhammer(interval<T>(q),interval<T>(q),int (k));
      }
      r=complex_nbd(mid,rad);
      res=qPochhammer(complex<interval<T> >(z),interval<T>(q),int(n))*r;
    }
    else{
      res=Gatteschi_qp(complex<interval<T> >(z),interval<T>(q));
    }
    return res;
  }
    template <class T> interval<T> qPVer2(const interval<T> & z,const interval<T>& q, int n=100, int K=100){
    interval<T> res,r,mid;
    mid=0.;
    T rad;
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("value of q must be positive");
    }
    if(q<0.9){
      while(abs(z)*pow(q,n)/(1-q)>=0.5){
	n=n+50;
      }
      rad=(2*pow(q,K*(K-1)*0.5)*pow(abs(z)*pow(q,n),K)
	   /qPochhammer(interval<T>(q),interval<T>(q),int (K))).upper();
      for(int k=0;k<=K-1;k++){
	mid=mid+pow(q,k*(k-1)*0.5)*pow(-z*pow(q,n),k)
	  /qPochhammer(interval<T>(q),interval<T>(q),int (k));
      }
      r=mid+rad*interval<T>(-1.,1.);
      res=qPochhammer(interval<T> (z),interval<T>(q),int(n))*r;
    }
    else{
      res=Gatteschi_qp(interval<T> (z),interval<T>(q));
    }
    return res;
}
  template <class T> interval<T> qPlarge(const interval<T> & z,const interval<T>& q, int n=5000){
    interval<T> res,qpn;
    // Ismail, M. E., & Zhang, R. (2006). Chaotic and periodic asymptotics for q-orthogonal polynomials. International Mathematics Research Notices, 2006.
    // z must be positive
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("value of q must be positive");
    }
    if (z<=0){
      throw std::domain_error("value of z must be positive");
    }
    T rad;
    qpn=qPochhammer(z,q,n);
    rad=(infinite_qPochhammer(interval<T>(-z*q*q),interval<T>(q))*z*pow(q,n)/(1-q)).upper();
    res=qpn*(1+rad*interval<T>(-1.,1.));
    return res;
  }
}
#endif
