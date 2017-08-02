// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp

// Verification program for q-special functions 
// by using double exponential formula with verified error bounds.


// References

// Ismail, M. E., & Zhang, R. (2016).
// Integral and Series Representations of $ q $-Polynomials and Functions: Part I.
// arXiv preprint arXiv:1604.08441.

// Okayama, T. (2013). Error Estimates with Explicit Constants for Sinc Quadrature and Sinc Indefinite Integration over Infinite Intervals. arXiv preprint arXiv:1302.1314.

// Tanaka, K., Sugihara, M., Murota, K., Mori, M. (2007), Function Classes for Double Exponential Integration Formulas, Numerische Mathematik, 111(4), 631-655. 

// Okayama, T., Matsuo, T., Sugihara, M. (2013), Error Estimates with Explicit Constants for Sinc Quadrature and Sinc Indefinite Integration. Numerische Mathematik, 124(2), 361-394.

// Zhang, R. (2008). Plancherel–Rotach Asymptotics for Certain Basic Hypergeometric Series. Advances in Mathematics, 217(4), 1588-1613.

#ifndef QFUNC_DE_HPP
#define QFUNC_DE_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <kv/constants.hpp>
#include <kv/complex.hpp>
#include <kv/Pochhammer.hpp>
#include <kv/Heine.hpp>
#include <cmath>
#include <algorithm>
namespace kv{

 
  template <class T> complex<interval<T> >Ramanujan_qAiry_DE4(const interval<T>&q, const complex<interval<T> >&z){
    // Compute Ramanujan q-Airy function with DE formula
    complex<interval<T> >res,mid1,mid2,i,int1,int2;
    interval<T> d,beta,pi,K,C,h,e,Amax,xhat;
    T rad;
    int m,N;
    m=100;
    while(abs(z)*pow(q,m)/(1-q)>=0.5){
      m=m+50;
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    i=kv::complex<interval<T> >::i();
    pi=kv::constants<interval<T> >::pi();
    e=kv::constants<interval<T> >::e();
      if(q<=exp(-0.5)){

      beta=-1/(2*log(q));
      K=qPochhammer(interval<T>(-abs(z)*sqrt(q)*exp(pi*0.5)),interval<T>(q),int(m))*(1+2*abs(z)*exp(pi*0.5)*pow(q,m+0.5)/(1-q));
      N=100;
      d=1.;
      while(2*pi*d*N/beta<e){
	N=N+10;
      }
      xhat=1.;
      while(pi/2-d-2*exp(-xhat)*sin(d)<0){
	xhat=xhat+1;
      }
      interval<T> alpha;
      alpha=xhat-exp(-xhat);
      while(alpha<0){
	xhat=xhat+1;
      }
      Amax=std::max(std::max((exp(-exp(-xhat))+1)*(1+exp(xhat))/(exp(xhat)-exp(exp(-xhat)))/exp(beta*exp(-1/e))*exp(beta),4*exp(beta+1)),4/(1-exp(-alpha)));
      C=2*Amax/beta*(1+2*exp(-beta)/(1-exp(-beta*e)));
      h=log(2*pi*N*d/beta)/N;
      rad=(C*exp(-2*pi*d/h)).upper();
      mid1=0.;
      mid2=0.;
      for(int k1=-N;k1<=N;k1++){
	mid1=mid1+infinite_qPochhammer(complex<interval<T> >(z*sqrt(q)*exp(i*exp(k1*h-exp(-k1*h)))),interval<T>(q))
	  *exp(-beta*exp(2*k1*h-2*exp(-k1*h)))			       
	  *(1+exp(k1*h))*exp(-exp(-k1*h));
      }
      mid1=mid1*h;
      for(int k2=-N;k2<=N;k2++){
	mid2=mid2+infinite_qPochhammer(complex<interval<T> >(z*sqrt(q)*exp(-i*exp(k2*h-exp(-k2*h)))),interval<T>(q))
	  *exp(-beta*exp(2*k2*h-2*exp(-k2*h)))
	  *(1+exp(k2*h))*exp(-exp(-k2*h));
      }
      mid2=mid2*h;    
      
      int1=complex_nbd(mid1,rad);
      int2=complex_nbd(mid2,rad);
      
      res=(int1+int2)/sqrt(2*pi*log(1/q));
      return res;
      }
  
      else{
	throw std::domain_error("value of q must be smaller");
      }
  }
  template <class T> complex<interval<T> >Jackson2_DE4(const complex<interval<T> >&z, const int &nu,const interval<T>&q){
    // Compute Jackson's 2nd q-Bessel function with DE formula
    complex<interval<T> >res,mid1,mid2,i,int1,int2;
    interval<T> d,beta,pi,K,C,h,e,Amax,xhat;
    T rad;
    int m,N;
    m=100;
    while(pow(q,m+nu+0.5)/(1-q)>=0.5){
      m=m+50;
    }
    while(abs(z*z*pow(q,m+nu+0.5)/(1-q)*0.25)>=0.5){
      m=m+50;
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    if(q<=exp(-0.5)){
      beta=-1/(2*log(q));
      i=kv::complex<interval<T> >::i();
      pi=kv::constants<interval<T> >::pi();
      e=kv::constants<interval<T> >::e();
      K=qPochhammer(interval<T>(-abs(z*z)*pow(q,nu+0.5)*exp(pi*0.5)*0.25),interval<T>(q),int(m))*(1+abs(z*z)*pow(q,m+nu+0.5)*exp(pi*0.5)/(1-q)*0.5)
	/infinite_qPochhammer(pow(q,nu+0.5)*exp(-pi*0.5),interval<T>(q))*(1+2*exp(pi*0.5)*pow(q,m+nu+0.5)/(1-q));
      N=1000;

      d=1.;
      while(d>=abs((nu+0.5)*log(q))){
	d=d*0.5;
      }
      while(2*pi*d*N/beta<e){
	N=N+10;
      }
      xhat=1.;
      while(pi/2-d-2*exp(-xhat)*sin(d)<0){
	xhat=xhat+1;
      }
      interval<T> alpha;
      alpha=xhat-exp(-xhat);
      while(alpha<0){
	xhat=xhat+1;
      }
      Amax=std::max(std::max((exp(-exp(-xhat))+1)*(1+exp(xhat))/(exp(xhat)-exp(exp(-xhat)))/exp(beta*exp(-1/e))*exp(beta),4*exp(beta+1)),4/(1-exp(-alpha)));

      C=2*Amax/beta*(1+2*exp(-beta)/(1-exp(-beta*e)));
      h=log(2*pi*N*d/beta)/N;
      rad=(C*exp(-2*pi*d/h)).upper();
      mid1=0.;
      mid2=0.;
      for(int k1=-N;k1<=N;k1++){
	mid1=mid1+infinite_qPochhammer(complex<interval<T> >(z*z*pow(q,nu+0.5)*0.25*exp(i*exp(k1*h-exp(-k1*h)))),interval<T>(q))
	  /infinite_qPochhammer(complex<interval<T> >(-pow(q,nu+0.5)*exp(i*exp(k1*h-exp(-k1*h)))),interval<T>(q))
	  *exp(-beta*exp(2*k1*h-2*exp(-k1*h)))			       
	  *(1+exp(k1*h))*exp(-exp(-k1*h));
      }
      mid1=mid1*h;
      for(int k2=-N;k2<=N;k2++){
	mid2=mid2+infinite_qPochhammer(complex<interval<T> >(z*z*pow(q,nu+0.5)*0.25*exp(-i*exp(k2*h-exp(-k2*h)))),interval<T>(q))
	  /infinite_qPochhammer(complex<interval<T> >(-pow(q,nu+0.5)*exp(-i*exp(k2*h-exp(-k2*h)))),interval<T>(q))
	  *exp(-beta*exp(2*k2*h-2*exp(-k2*h)))
	  *(1+exp(k2*h))*exp(-exp(-k2*h));
      }
      mid2=mid2*h;    
      
      int1=complex_nbd(mid1,rad);
      int2=complex_nbd(mid2,rad);
      
      res=(int1+int2)/sqrt(2*pi*log(1/q))/Euler(interval<T>(q))*pow(0.5*z,nu);
      return res;
      }

      else{
	throw std::domain_error("value of q must be smaller");
      }
    }

template <class T> complex<interval<T> >Jackson2_DE4(const complex<interval<T> >&z, const interval<T> &nu,const interval<T>&q){
    // Compute Jackson's 2nd q-Bessel function with DE formula
    complex<interval<T> >res,mid1,mid2,i,int1,int2;
    interval<T> d,beta,pi,K,C,h,e,Amax,xhat;
    T rad;
    int m,N;
    m=100;
    while(pow(q,m+nu+0.5)/(1-q)>=0.5){
      m=m+50;
    }
    while(abs(z*z*pow(q,m+nu+0.5)/(1-q)*0.25)>=0.5){
      m=m+50;
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    int floornu;
    floornu=std::floor(mid(-2*nu));
    if(nu<0 && floornu%2==1 && mid(-2*nu)-floornu<=0){
      throw std::domain_error("singularity on real axis");
      // reject negative half integers
    }
    if(q<=exp(-0.5)){
      beta=-1/(2*log(q));
      i=kv::complex<interval<T> >::i();
      pi=kv::constants<interval<T> >::pi();
      e=kv::constants<interval<T> >::e();
      K=qPochhammer(interval<T>(-abs(z*z)*pow(q,nu+0.5)*exp(pi*0.5)*0.25),interval<T>(q),int(m))*(1+abs(z*z)*pow(q,m+nu+0.5)*exp(pi*0.5)/(1-q)*0.5)
	/infinite_qPochhammer(pow(q,nu+0.5)*exp(-pi*0.5),interval<T>(q))*(1+2*exp(pi*0.5)*pow(q,m+nu+0.5)/(1-q));
      N=1000;

      d=1.;
      while(d>=abs((nu+0.5)*log(q))){
	d=d*0.5;
      }
      while(2*pi*d*N/beta<e){
	N=N+10;
      }
      xhat=1.;
      while(pi/2-d-2*exp(-xhat)*sin(d)<0){
	xhat=xhat+1;
      }
      interval<T> alpha;
      alpha=xhat-exp(-xhat);
      while(alpha<0){
	xhat=xhat+1;
      }
      Amax=std::max(std::max((exp(-exp(-xhat))+1)*(1+exp(xhat))/(exp(xhat)-exp(exp(-xhat)))/exp(beta*exp(-1/e))*exp(beta),4*exp(beta+1)),4/(1-exp(-alpha)));

      C=2*Amax/beta*(1+2*exp(-beta)/(1-exp(-beta*e)));
      h=log(2*pi*N*d/beta)/N;
      rad=(C*exp(-2*pi*d/h)).upper();
      mid1=0.;
      mid2=0.;
      for(int k1=-N;k1<=N;k1++){
	mid1=mid1+infinite_qPochhammer(complex<interval<T> >(z*z*pow(q,nu+0.5)*0.25*exp(i*exp(k1*h-exp(-k1*h)))),interval<T>(q))
	  /infinite_qPochhammer(complex<interval<T> >(-pow(q,nu+0.5)*exp(i*exp(k1*h-exp(-k1*h)))),interval<T>(q))
	  *exp(-beta*exp(2*k1*h-2*exp(-k1*h)))			       
	  *(1+exp(k1*h))*exp(-exp(-k1*h));
      }
      mid1=mid1*h;
      for(int k2=-N;k2<=N;k2++){
	mid2=mid2+infinite_qPochhammer(complex<interval<T> >(z*z*pow(q,nu+0.5)*0.25*exp(-i*exp(k2*h-exp(-k2*h)))),interval<T>(q))
	  /infinite_qPochhammer(complex<interval<T> >(-pow(q,nu+0.5)*exp(-i*exp(k2*h-exp(-k2*h)))),interval<T>(q))
	  *exp(-beta*exp(2*k2*h-2*exp(-k2*h)))
	  *(1+exp(k2*h))*exp(-exp(-k2*h));
      }
      mid2=mid2*h;    
      
      int1=complex_nbd(mid1,rad);
      int2=complex_nbd(mid2,rad);
      
      res=(int1+int2)/sqrt(2*pi*log(1/q))/Euler(interval<T>(q))*pow(0.5*z,nu);
      return res;
      }
      else{
	throw std::domain_error("value of q must be smaller");
      }
    }
  template <class T> complex<interval<T> >Jackson2_DE4(const complex<interval<T> >&z, const complex<interval<T> >&nu,const interval<T>&q){
    // Compute Jackson's 2nd q-Bessel function with DE formula
    complex<interval<T> >res,mid1,mid2,i,int1,int2;
    interval<T> d,beta,pi,K,C,h,e,Amax,xhat;
    T rad;
    int m,N;
    m=100;
    while(abs(pow(q,m+nu+0.5))/(1-q)>=0.5){
      m=m+50;
    }
    while(abs(z*z*pow(q,m+nu+0.5)/(1-q)*0.25)>=0.5){
      m=m+50;
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    int floornu;
    floornu=std::floor(mid(-2*nu.real()));
    if(nu.real()<0 && floornu%2==1 && mid(-2*nu.real())-floornu<=0 && nu.imag()==0){
      throw std::domain_error("singularity on real axis");
      // reject negative half integers
    }
    if(q<=exp(-0.5)){
      beta=-1/(2*log(q));
      i=kv::complex<interval<T> >::i();
      pi=kv::constants<interval<T> >::pi();
      e=kv::constants<interval<T> >::e();
      K=qPochhammer(interval<T>(-abs(z*z*pow(q,nu+0.5))*exp(pi*0.5)*0.25),interval<T>(q),int(m))*(1+abs(z*z*exp(pi*0.5)*pow(q,m+nu+0.5))/(1-q)*0.5)
	/infinite_qPochhammer(abs(pow(q,nu+0.5)*exp(-pi*0.5)),interval<T>(q))*(1+2*exp(pi*0.5)*abs(pow(q,m+nu+0.5))/(1-q));
      N=1000;

      d=1.;
      while(d>=abs((nu.real()+0.5)*log(q))){
	d=d*0.5;
      }
      while(2*pi*d*N/beta<e){
	N=N+10;
      }
      xhat=1.;
      while(pi/2-d-2*exp(-xhat)*sin(d)<0){
	xhat=xhat+1;
      }
      interval<T> alpha;
      alpha=xhat-exp(-xhat);
      while(alpha<0){
	xhat=xhat+1;
      }
      Amax=std::max(std::max((exp(-exp(-xhat))+1)*(1+exp(xhat))/(exp(xhat)-exp(exp(-xhat)))/exp(beta*exp(-1/e))*exp(beta),4*exp(beta+1)),4/(1-exp(-alpha)));

      C=2*Amax/beta*(1+2*exp(-beta)/(1-exp(-beta*e)));
      h=log(2*pi*N*d/beta)/N;
      rad=(C*exp(-2*pi*d/h)).upper();
      mid1=0.;
      mid2=0.;
      for(int k1=-N;k1<=N;k1++){
	mid1=mid1+infinite_qPochhammer(complex<interval<T> >(z*z*pow(q,nu+0.5)*0.25*exp(i*exp(k1*h-exp(-k1*h)))),interval<T>(q))
	  /infinite_qPochhammer(complex<interval<T> >(-pow(q,nu+0.5)*exp(i*exp(k1*h-exp(-k1*h)))),interval<T>(q))
	  *exp(-beta*exp(2*k1*h-2*exp(-k1*h)))			       
	  *(1+exp(k1*h))*exp(-exp(-k1*h));
      }
      mid1=mid1*h;
      for(int k2=-N;k2<=N;k2++){
	mid2=mid2+infinite_qPochhammer(complex<interval<T> >(z*z*pow(q,nu+0.5)*0.25*exp(-i*exp(k2*h-exp(-k2*h)))),interval<T>(q))
	  /infinite_qPochhammer(complex<interval<T> >(-pow(q,nu+0.5)*exp(-i*exp(k2*h-exp(-k2*h)))),interval<T>(q))
	  *exp(-beta*exp(2*k2*h-2*exp(-k2*h)))
	  *(1+exp(k2*h))*exp(-exp(-k2*h));
      }
      mid2=mid2*h;    
      
      int1=complex_nbd(mid1,rad);
      int2=complex_nbd(mid2,rad);
      
      res=(int1+int2)/sqrt(2*pi*log(1/q))/Euler(interval<T>(q))*pow(0.5*z,nu);
      return res;
      }
      else{
	throw std::domain_error("value of q must be smaller");
      }
    }
 
}
#endif
    
