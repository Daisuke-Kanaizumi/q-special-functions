// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp
#ifndef JACKSOND_HPP
#define JACKSOND_HPP
#include <cmath>
namespace kv {
  template <class T, class F>
  interval<T>JacksonD(F f,const interval<T> & x,const interval<T>& q){
    interval<T>y;
    y=q*x;
    return (f(x)-f(y))/(1-q)/x;
  }
  template <class T, class F>
  complex<interval<T> >JacksonD(F f,const complex<interval<T> >& x,const interval<T>& q){
    complex<interval<T> >y;
    y=q*x;
    return (f(x)-f(y))/(1-q)/x;
  }

  template <class T, class F>
  interval<T>mJacksonD(F f,const interval<T> & x,const interval<T>& q){
    interval<T>y;
    y=q*x;
    return (f(x)-f(y))/(1-q)/x/(1+(1-q)*f(y));
    // modified with q-subtraction (Borges)
  }
  template <class T, class F>
  interval<T>modJacksonD(F f,const interval<T> & x,const interval<T>& q){
    interval<T>y;
    y=q*x;
    return (f(x)-q*f(y))/(1-q)/x;
    // modified with q-subtraction (Hahn)
  }
  template <class T, class F>
  interval<T>mmJacksonD(F f,const interval<T> & x,const interval<T>& q){
    interval<T>y;
    int n ;
    n=10000;
    y=q*x;
    return (f(x)*sqrt(1-(1-pow(q,1./n))*f(y)*f(y))-f(y)*sqrt(1-(1-pow(q,1./n))*f(x)*f(x)))/(1-q)/x;
    // modified with q-subtraction
  }
  template <class T, class F>
  interval<T>mmmJacksonD(F f,const interval<T> & x,const interval<T>& q){
    interval<T>y;
    int n ;
    n=10000;
    y=q*x;
    return (f(x)*sqrt(1-(1-pow(q,1./n))*f(y)*f(y)*f(y)*f(y))
	    -f(y)*sqrt(1-(1-pow(q,1./n))*f(x)*f(x)*f(x)*f(x)))/(1-q)/x/(1+(1-q)*f(x)*f(x)*f(y)*f(y));
    // modified with q-subtraction made by Euler's formal group
  }

  template <class T, class F>
  interval<T>JacksonD2(F f,const interval<T> &  x,const interval<T>& q){
    interval<T>y;
    y=q*x;
    return (JacksonD(f,x,q)-JacksonD(f,y,q))/(1-q)/x;
  }
  template <class T, class F>
  interval<T>JacksonD3(F f,const interval<T> &  x,const interval<T>& q){
    interval<T>y;
    y=q*x;
    return (JacksonD2(f,x,q)-JacksonD2(f,y,q))/(1-q)/x;
  }

  template <class T, class F>
  interval<T>Dh(F f,const interval<T> & x,const interval<T>& h){
    return (f(x+h)-f(x))/h;
  }
  template <class T, class F>
  interval<T>qmvf(F f,const interval<T> & x,const interval<T>& q){
    interval<T>m;
    m=mid(x);
    return f(m)+JacksonD(f,x,q)*(x-m);
  }
  template <class T, class F>
  interval<T>mqmvf(F f,const interval<T> & x,const interval<T>& q){
    interval<T>m;
    m=mid(x);
    return f(m)+mJacksonD(f,x,q)*(x-m);
    // test with x-sin(x), x-cos(x)
  }
 template <class T, class F>
  interval<T>modqmvf(F f,const interval<T> & x,const interval<T>& q){
    interval<T>m;
    m=mid(x);
    return f(m)+modJacksonD(f,x,q)*(x-m);
 // test with x*x*x-exp(x)
  }
  template <class T, class F>
  interval<T>mmqmvf(F f,const interval<T> & x,const interval<T>& q){
    interval<T>m;
    m=mid(x);
    return f(m)+mmJacksonD(f,x,q)*(x-m);
  }
  template <class T, class F>
  interval<T>mmmqmvf(F f,const interval<T> & x,const interval<T>& q){
    interval<T>m;
    m=mid(x);
    return f(m)+mmmJacksonD(f,x,q)*(x-m);
  }

  template <class T, class F>
  complex<interval<T> >qmvf(F f,const complex<interval<T> >& x,const interval<T> & q){
    complex<interval<T> >m;
    m=complex<interval<T> >(mid(x.real()),mid(x.imag()));
    return f(m)+JacksonD(f,x,q)*(x-m);
   
  }

  template <class T, class F>
  interval<T>qSchwarzD(F f,const interval<T> & x,const interval<T>& q){
    return JacksonD3(f,x,q)/JacksonD(f,x,q)
      -1.5*pow(JacksonD2(f,x,q)/JacksonD(f,x,q),2);
  }
  template <class T, class F>
  interval<T>qSmvf(F f,const interval<T> & x,const interval<T>& q){
    interval<T>m;
    m=mid(x);
    return f(m)+qSchwarzD(f,x,q)*(x-m);
  }
}
#endif 

