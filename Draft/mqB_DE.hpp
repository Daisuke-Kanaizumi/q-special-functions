// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp

// Verification program for modified q-Bessel functions 
// by using double exponential formula with verified error bounds.
// References

// Zhang, R. (2008). Plancherel–Rotach Asymptotics for Certain Basic Hypergeometric Series. Advances in Mathematics, 217(4), 1588-1613.
// Okayama, T. (2013). Error estimates with explicit constants for Sinc quadrature and Sinc indefinite integration over infinite intervals. arXiv preprint arXiv:1302.1314.
// Ismail, M. E. (1981). The basic Bessel functions and polynomials. SIAM Journal on Mathematical Analysis, 12(3), 454-468.
#ifndef MQB_DE_HPP
#define MQB_DE_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <kv/constants.hpp>
#include <kv/complex.hpp>
#include <kv/defint.hpp>
#include <kv/Pochhammer.hpp>
#include <kv/Heine.hpp>
#include <cmath>
#include <algorithm>
namespace kv{

  template <class TT> struct mqBreal {
    TT z, q, nu;
    int n; // Setting parameters
    mqBreal(TT z, TT q, TT nu,int n) : z(z),q(q),nu(nu),n(n) {}
    
    template <class T> T operator() (const T& t) {
      complex<T>  pro;
      T proreal;
      complex<T> i;
      pro=1.;
      i=complex<T>::i();
      for(int k=0;k<=n-1;k++){
	pro=pro/(1-z*0.5*pow(T(q),k)*exp(i*t))/(1-z*0.5*pow(T(q),k)*exp(-i*t));
      }
      pro=pro*cos(T(nu)*t);
      proreal=pro.real();
      return proreal;
    }
  };
  
  template <class TT> struct mqBimag {
    TT z, q, nu;
    int n; // Setting parameters
    mqBimag(TT z, TT q, TT nu,int n) : z(z),q(q),nu(nu),n(n) {}
    
    template <class T> T operator() (const T& t) {
      complex<T>  pro;
      T proimag;
      complex<T> i;
      pro=1.;
      i=complex<T>::i();
      for(int k=0;k<=n-1;k++){
	pro=pro/(1-z*0.5*pow(T(q),k)*exp(i*t))/(1-z*0.5*pow(T(q),k)*exp(-i*t));
      }
      pro=pro*cos(T(nu)*t);
      proimag=pro.imag();
      return proimag;
      
    }
  };

  template <class T> complex<interval<T> >modified_qBesselI2_DE(const interval<T> & z,const interval<T> & nu,const interval<T>& q){
    //verification program for 2nd modified q-Bessel function I2
    if(abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    if(nu<=0){
      throw std::domain_error("value of nu must be more than 0");
    }
    interval<T> pi,realint,imagint,second,e,d,h,K,C3,sum,intrad;
    complex<interval<T> >res,first;
    pi=constants<interval<T> >::pi();
    e=constants<interval<T> >::e();
    int m,n,M,N;
    m=100;
    while(abs(z*pow(q,m))*0.5/(1-q)>=0.5){
      m=m+10;
    }
    d=1.;
  
    n=150;
    M=n;
    h=log(4*d*n)/T(n);
    N=n-std::floor(mid(log(nu)/h));
    
    while(n<=nu*e*0.25*d){
      n=n+10;
    }
    intrad=pow(1+2*abs(z)*pow(q,m)/(1-q)*interval<T>(-1.,1.),2);
    realint=defint(mqBreal<interval<T> >(z,q,nu,m),interval<T>(0.),interval<T>(pi),10,10);
    imagint=defint(mqBimag<interval<T> >(z,q,nu,m),interval<T>(0.),interval<T>(pi),10,10);
    first=intrad*complex<interval<T> >(realint,imagint)/pi;
    
    K=(1+2*z*pow(q,m)/(1-q))/infinite_qPochhammer(interval<T>(abs(z)*0.5),interval<T>(q))
      /qPochhammer(interval<T>(z*0.5),interval<T>(q),int(m));

    C3=2*K*(2/(1-exp(-pi*e*0.5))/cos(d)/pow(cos(pi*0.5*sin(d)),nu+1)+exp(pi*nu*0.5));

    sum=0.;
    for(int k=-M;k<=N;k++){
      sum=sum+exp(-(nu)*log(1+exp(pi*sinh(k*h))))
	/infinite_qPochhammer(interval<T> (-z*0.5*(1.+exp(pi*sinh(k*h)))),interval<T>(q))
	/infinite_qPochhammer(interval<T> (-0.5*z/(1.+exp(pi*sinh(k*h)))),interval<T>(q))
	*pi*exp(pi*sinh(k*h))*cosh(k*h)/(1.+exp(pi*sinh(k*h)));
    }

    second=(h*sum+C3*exp(-2*pi*d/h)*interval<T>(-1.,1.))*sin(pi*nu)/pi;
    res=(first-second)*infinite_qPochhammer(interval<T> (z*z*0.25),interval<T>(q));
    return res;
    // to be repaired
  }
}
#endif
