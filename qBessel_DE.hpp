// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University

// Verification program for q-Bessel functions 
// by using double exponential formula with verified error bounds.


// References

// Zhang, R. (2008). Plancherel–Rotach Asymptotics for Certain Basic Hypergeometric Series. Advances in Mathematics, 217(4), 1588-1613.

// Okayama, T., Matsuo, T., & Sugihara, M. (2010). Error estimates with explicit constants for the tanh rule and the DE formula for indefinite integrals. JSIAM Letters, 2, 13-16.

// Okayama, T., Matsuo, T., & Sugihara, M. (2009). Error estimates with explicit constants for Sinc approximation, Sinc quadrature and Sinc indefinite integration.Mathematical Engineering Technical Reports. 2009-01, The University of Tokyo.
#ifndef QBESSEL_DE_HPP
#define QBESSEL_DE_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <kv/constants.hpp>
#include <kv/complex.hpp>
#include <kv/Pochhammer.hpp>
#include <kv/Heine.hpp>
#include <cmath>
#include <algorithm>
namespace kv{
  template <class T> complex<interval<T> >Jackson2_DE(const complex<interval<T> >&z, const complex<interval<T> >&nu,const interval<T>&q){
    // Compute Jackson's 2nd q-Bessel function with DE formula and finite integral representation
    // Rahman(1987), An Integral Representation and Some Transformation Properties of q-Bessel Functions, Journal of Mathematical Analysis and Applications 125
    complex<interval<T> >i,res,integral,sum;
    interval<T>e,pi,K,d,h,C3,ImMax;
    sum=0.;
    int N=100,n=100;
    T rad;
    if(nu.real()<=0){
      throw std::domain_error("real part of nu must be positive");
    }
    if(abs(z*pow(q,nu/2+n+0.5)*exp(ImMax))/2/(1-q)>=0.5){
      n=n+10;
    }
    if(abs(pow(q,nu+n)*exp(ImMax))/2/(1-q)>=0.5){
      n=n+10;
    }
    i=kv::complex<interval<T> >::i();
    pi=kv::constants<interval<T> >::pi();
    e=kv::constants<interval<T> >::e();
    d=1.;
    h=log(4*d*N)/N;
    C3=1/pow(cos(pi*0.5*sin(d)),2)/cos(d);
    ImMax=sin(pi*sin(d))/(cos(pi*sin(d))+1);
    K=pow((1+2*pow(q,n)*exp(2*ImMax)/(1-q))*qPochhammer(interval<T>(-exp(2*ImMax)),interval<T>(q),int(n)),2)
      *pow((1+abs(z*pow(q,n+(nu+1)*0.5))*exp(ImMax)/(1-q))*qPochhammer(interval<T>(-0.5*abs(z*pow(q,(nu+1)*0.5))*exp(ImMax)),interval<T>(q),int(n)),2)
      *pow((1+abs(pow(q,nu+n))*exp(2*ImMax)/(1-q))/qPochhammer(interval<T> (abs(pow(q,nu))*exp(-2*ImMax)),interval<T>(q),int(n)),2);
    rad=(2*K*pi*exp(-2*pi*d/h)*(exp(pi*0.5)+2*C3/(1-exp(-pi*e*0.5)))).upper();
    for(int k=-N;k<=N;k++){
      interval<T>psi;
      psi=pi*0.5*tanh(pi*0.5*sinh(k*h))+pi*0.5;
      sum=sum+infinite_qPochhammer(complex<interval<T> >(exp(2*i*psi)),interval<T>(q))
	*infinite_qPochhammer(complex<interval<T> >(exp(-2*i*psi)),interval<T>(q))
	/infinite_qPochhammer(complex<interval<T> >(exp(2*i*psi)*pow(q,nu)),interval<T>(q))
	/infinite_qPochhammer(complex<interval<T> >(exp(-2*i*psi)*pow(q,nu)),interval<T>(q))
	*infinite_qPochhammer(complex<interval<T> >(-i*z*0.5*pow(q,(nu+1)*0.5)*exp(i*psi)),interval<T>(q))
	*infinite_qPochhammer(complex<interval<T> >(-i*z*0.5*pow(q,(nu+1)*0.5)*exp(-i*psi)),interval<T>(q))
	*pi*pi*0.25*cosh(k*h)/pow(cosh(pi*0.5*sinh(k*h)),2);
    }
    sum=sum*h;
    integral=complex_nbd(sum,rad);
    res=integral*pow(z/2,nu)
      *infinite_qPochhammer(complex<interval<T> >(pow(q,2*nu)),interval<T>(q))
      /infinite_qPochhammer(complex<interval<T> >(pow(q,nu)),interval<T>(q))/(2.*pi);
       
    return res;
  }
}
#endif
