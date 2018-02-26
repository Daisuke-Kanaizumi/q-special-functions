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

  // next make higher order
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

