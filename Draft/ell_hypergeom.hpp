// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Date: May 11th, 2017

#ifndef ELL_HYPERGEOM_HPP
#define ELL_HYPERGEOM_HPP

// verification program for the elliptic hypergeometric function

// Elliptic hypergeometric series were introduced by Frenkel & Turaev (1997) in their study of elliptic 6-j symbols (elliptic solutions of the Yang-Baxter eq.)

// For details of elliptic hypergeometric series, see Gasper & Rahman (2004) or Spiridonov (2008)
// Gasper, George; Rahman, Mizan (2004), Basic hypergeometric series, Encyclopedia of Mathematics and its Applications, 96 (2nd ed.), Cambridge University Press, 
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/complex.hpp>
#include <kv/Pochhammer.hpp>
#include <cmath>

namespace kv{
  // Before 2E1, the modified Jacobi theta function and the elliptic Pochhammer symbol will be implemented

  template <class T> interval<T> modified_Jacobi_theta(const interval<T> &a,const interval<T> &q){
    interval<T> res;
    res=infinite_qPochhammer(interval<T> (a),interval<T>(q))
      *infinite_qPochhammer(interval<T> (q/a),interval<T>(q));
    return res;
  } 
  template <class T> complex<interval<T> > modified_Jacobi_theta(const complex<interval<T> >&a,const interval<T> &q){
    complex<interval<T> >res;
    res=infinite_qPochhammer(complex<interval<T> >(a),interval<T>(q))
      *infinite_qPochhammer(complex<interval<T> >(q/a),interval<T>(q));
    return res;
  } 
  template <class T> interval<T>  elliptic_Pochhammer(const interval<T> &a,const interval<T> &q,const interval<T> &p,const int &n){
    interval<T> res,pro;
    int i;
    pro=1.;
    for(i=1;i<=n;i++){
      pro=pro*modified_Jacobi_theta(interval<T> (a*pow(q,i-1)),interval<T>(p));
    } 
    res=pro;
    return res;
  }
  template <class T> complex<interval<T> > elliptic_Pochhammer(const complex<interval<T> >&a,const interval<T> &q,const interval<T> &p,const int &n){
    complex<interval<T> >res,pro;
    int i;
    pro=1.;
    for(i=1;i<=n;i++){
      pro=pro*modified_Jacobi_theta(complex<interval<T> >(a*pow(q,i-1)),interval<T>(p));
    } 
    res=pro;
    return res;
  }
  template <class T> interval<T>  _2E_1(const interval<T> &a1,const interval<T> &a2,const interval<T> &b1,const interval<T> &q,const interval<T> &p,interval<T> (z)){
    // References
    // Zhang (2008), Plancherel-Rotach asymptotics for certain basic hypergeometric series, Advances in Mathematics, 217
    // Gasper, Rahman (2004), Basic hypergeometric series, Encyclopedia of Mathematics and its Applications, 96 (2nd ed.), Cambridge University Press, 

    int l,m,n,N;
    l=10;
    m=10;
    N=4;
    interval<T> res,mid,first;
    mid=1.;
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    if (p>=1){
      throw std::domain_error("value of p must be under 1");
   }
    if (p<=0){
      throw std::domain_error("p must be positive");
    }
    while(pow(q,N)>=abs(1/a1)){
      N=N+10;
    }
    while(pow(q,N)>=abs(1/a2)){
     N=N+10;
    }
    for(int n=1;n<=N-1;n++){
      mid=mid+elliptic_Pochhammer(interval<T>(a1),interval<T>(q),interval<T>(p),int (n))*elliptic_Pochhammer(interval<T>(a2),interval<T>(q),interval<T>(p),int (n))*pow(z,n)
	/elliptic_Pochhammer(interval<T>(b1),interval<T>(q),interval<T>(p),int (n))/elliptic_Pochhammer(interval<T>(q),interval<T>(q),interval<T>(p),int (n));
    }
    first=abs(elliptic_Pochhammer(interval<T>(a1),interval<T>(q),interval<T>(p),int (N))*elliptic_Pochhammer(interval<T>(a2),interval<T>(q),interval<T>(p),int (N))*pow(z,N)
	      /elliptic_Pochhammer(interval<T>(b1),interval<T>(q),interval<T>(p),int (N))/elliptic_Pochhammer(interval<T>(q),interval<T>(q),interval<T>(p),int (N)));
    
    while(pow(p,m)/(1-p)>=0.5){
      m=m+10;
    }
    while(abs(b1)*pow(p,m)/(1-p)>=0.5){
      m=m+10;
    }
    while(abs(a1)*pow(p,m)/(1-p)>=0.5){
      m=m+10;
    }
    while(abs(a2)*pow(p,m)/(1-p)>=0.5){
      m=m+10;
    }
    while(q*pow(p,l+1)>=0.5){
      l=l+10;
    }
    while(abs(b1)*pow(p,l+1)>=0.5){
      l=l+10;
    }
    while(abs(a1)*pow(p,l+1)>=0.5){
      l=l+10;
    }
    while(abs(a2)*pow(p,l+1)>=0.5){
      l=l+10;
    }
    interval<T>ratio1,ratio2,ratio3,ratio4,ratio;
    T rad;
    ratio1=(1+2*abs(b1)*pow(q,N)*pow(p,m)/(1-p))*(1+2*pow(q,N+1)*pow(p,m)/(1-p))
      *abs(qPochhammer(interval<T>(b1*pow(q,N)),interval<T>(p),int(m)))
      *abs(qPochhammer(interval<T>(pow(q,N+1)),interval<T>(p),int(m)));
    ratio2=(1+2*abs(a1)*pow(q,N)*pow(p,m)/(1-p))*(1+2*abs(a2)*pow(q,N)*pow(p,m)/(1-p))
      /abs(qPochhammer(interval<T>(a1*pow(q,N)),interval<T>(p),int(m)))
      /abs(qPochhammer(interval<T>(a2*pow(q,N)),interval<T>(p),int(m)));
    ratio3=pow(1-p,2)*q*abs(b1)/pow(Euler(interval<T>(p)),2)
      *(1+2*abs(b1)*pow(q,N)*pow(p,l+1))*(1+2*pow(q,N+1)*pow(p,l+1))
      *abs(qPochhammer(interval<T>(p*(1-p)*b1*pow(q,N)),interval<T>(p),int(l)))
      *abs(qPochhammer(interval<T>(p*(1-p)*pow(q,N+1)),interval<T>(p),int(l)));
    ratio4=abs(a1*a2)/(pow(Euler(interval<T>(p)),2)).lower()
      *(1+2*abs(a1)*pow(q,N)*pow(p,l+1))*(1+2*abs(a2)*pow(q,N)*pow(p,l+1))
      /abs(qPochhammer(interval<T>(p*(1-p)*a1*pow(q,N)),interval<T>(p),int(l)))
      /abs(qPochhammer(interval<T>(p*(1-p)*a2*pow(q,N)),interval<T>(p),int(l)));
    ratio=ratio1*ratio2*ratio3*ratio4*abs(z);
    if(ratio<1){
      rad=(first/(1-ratio)).upper();
      res=mid+rad*interval<T>(-1.,1.);
      return res;
    }
    else{
      throw std::domain_error("ratio is more than 1");
    } 
  }
  
}
#endif
