// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
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
  interval<T>HahnD(F f,const interval<T> & x,const interval<T>& q){
    interval<T>y,omega;
    y=q*x;
    omega=-2.5;
    return (f(x)-f(y+omega))/(x-y-omega);
  }
  // Hahn, W.: Über Orthogonalpolynome, die q-Differenzenlgleichungen genügen. Math. Nachr. 2, 4–34 (1949) 
  template <class T, class F>
  interval<T>hqmvf(F f,const interval<T> & x,const interval<T>& q){
    interval<T>m;
    m=mid(x);
    return f(m)+HahnD(f,x,q)*(x-m);
  }

  template <class T, class F>
  interval<T>RubinD(F f,const interval<T> & x,const interval<T>& q){
    interval<T>y;
    y=q*x;
    return (f(x/q)+f(-x/q)-2*f(-x)-f(y)+f(-y))/(1-q)/x*0.5;
  }
  // Rubin RL:A q2-analogue operator for q2-analogue Fourier analysis. J. Math. Anal. Appl. 1997, 212(2):571–582. 10.1006/jmaa.1997.5547
  // Rubin RL:Duhamel solutions of non-homogeneous q2-analogue wave equations. Proc. Am. Math. Soc. 2007, 135(3):777–785. 10.1090/S0002-9939-06-08525-X
  template <class T, class F>
  interval<T>rqmvf(F f,const interval<T> & x,const interval<T>& q){
    interval<T>m;
    m=mid(x);
    return f(m)+RubinD(f,x,q)*(x-m);
  }

  template <class T, class F>
  interval<T>mJacksonD(F f,const interval<T> & x,const interval<T>& q){
    interval<T>y;
    y=q*x;
    return (f(x)-f(y))/(1-q)/x/(1+(1-q)*f(y));
    // modified with q-subtraction (Borges, 2004)
  }
  template <class T, class F>
  interval<T>modJacksonD(F f,const interval<T> & x,const interval<T>& q){
    interval<T>y;
    y=q*x;
    return (f(x)-q*f(y))/(1-q)/x;
    // modified with q-subtraction (Jackson-Hahn-Cigler)
    // x*x*x*x*x-cos(x),x*x*x*x*x-sin(x),x*x*x*x*x-sqrt(x)
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
  interval<T>MJacksonD(F f,const interval<T> & x,const interval<T>& q){
    interval<T>y;
    y=q*x;
    return (f(x)-f(y))/(1-q)/x/(1-(1-q)*f(x)*f(y));
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
  interval<T>Mqmvf(F f,const interval<T> & x,const interval<T>& q){
    interval<T>m;
    m=mid(x);
    return f(m)+MJacksonD(f,x,q)*(x-m);
    // test with x-sin(x), x-cos(x)
  }

  template <class T, class F>
  interval<T>mmmqmvf(F f,const interval<T> & x,const interval<T>& q){
    interval<T>m;
    m=mid(x);
    return f(m)+mmmJacksonD(f,x,q)*(x-m);
  }

}
#endif 

