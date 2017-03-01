// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp

// verification program for Pochhammer symbol

#ifndef POCHHAMMER_HPP
#define POCHHAMMER_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

namespace kv{
  template <class T> interval<T> Pochhammer(const interval<T>& a, const int& n){
    interval<T> res;
    int i;
    i=1;
    if(n==0){
      return interval<T>(1.);
    }
    if(n==1){
      return a;
    }
    if(n<0){
      std::cout<<"Pochhammer symbol is not defined"<<std::endl;
    }
    else{
      res=a;
      for(i=1;i<=n-1;i++){
	res*=a+i;
      }
      return res;
    }
  }
 template <class T> interval<T> Euler(const interval<T>& q){
   // verify Euler function by using the pentagonal number theorem
   int K=500;
   int j=1;
   int k=1;
   interval<T>a,res;
   T rad;
   if (abs(q)>=1){
     std::cout<<"this program is unavailable"<<std::endl;
   }
   else{
     a=1-q;
     for (k=1; k<=K; k++){
       j = -1*j;
       
       a = a+j*(1-pow(q,2*k+1))*pow(q,k*(3*k+1)/2);
     }
     rad=abs((1-pow(q,2*K+3))*pow(q,(K+1)*(3*K+4)/2)).upper();
     // .upper() allows to output supremum
     res=a+rad;
     return res;
   }
 }
  template <class T> interval<T> Karpelevich(const interval<T>& z,const interval<T>& q){
  // reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
  // Olshanetsky, Rogov 1995
    interval<T>a,b,res;
    T rad;
    int K=500;
    int j=1;
    int k=1;
    if (abs(q)>=1){
      std::cout<<"this program is unavailable"<<std::endl;
    }
    if (abs(z)>=1){
      std::cout<<"this program is unavailable"<<std::endl;
    }
    else{
      a=1/(1-z);
      b=1.;
      // derived from the definition of q-Pochhammer symbol
      for (k=1; k<=K; k++){
	j = -1*j;
	b=b*(1-pow(q,k));
	a=a+j*pow(q,k*(k+1)/2)/(b*(1-z*pow(q,k)));
      }
      b=b*(1-pow(q,K+1));
      rad=abs(pow(q,(K+1)*(K+2)/2)/(b*(1-z*pow(q,K+1)))).upper();
      res=a+rad;
      return res;
    }
  }
  template <class T> complex<interval<T> > Karpelevich(const complex< interval<T> >& z,const interval<T>& q){
  // reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
  // Olshanetsky, Rogov 1995
    complex<interval<T> >res;
    T rad1,rad2;
    interval<T>x,y,a1,a2,b1,b2,real,imag;
    x=z.real();
    y=z.imag();
    int K=500;
    int j=1;
    int k=1;
    int kk=1;
    int jj=1;
    if (abs(q)>=1){
      std::cout<<"this program is unavailable"<<std::endl;
    }
    if (abs(z)>=1){
      std::cout<<"this program is unavailable"<<std::endl;
    }
    if (z.real()<0){
      std::cout<<"this program is unavailable"<<std::endl;
    }
    else{
      a1=(1-x)/(pow(1-x,2)+y*y);
      b1=1.;
      // derived from the definition of q-Pochhammer symbol
      for (k=1; k<=K; k++){
	j = -1*j;
	b1=b1*(1-pow(q,k));
	a1=a1+j*pow(q,k*(k+1)/2)*(1-x*pow(q,k))/(b1*(pow(1-x*pow(q,k),2)+pow(q,2*k)*y*y));
      }
      b1=b1*(1-pow(q,K+1));
      rad1=abs(pow(q,(K+2)*(K+1)/2)*(1-x*pow(q,K+1))/(b1*(pow(1-x*pow(q,K+1),2)+pow(q,2*(K+1))*y*y))).upper();
      real=a1+rad1;
      
      a2=y/(pow(1-x,2)+y*y);
      b2=1.;
      // derived from the definition of q-Pochhammer symbol
      for (kk=1; kk<=K; kk++){
	jj = -1*jj;
	b2=b2*(1-pow(q,kk));
		 a2=a2+jj*pow(q,kk*(kk+1)/2)*y*pow(q,kk)/(b2*(pow(1-x*pow(q,kk),2)+pow(q,2*kk)*y*y));
      }
      b2=b2*(1-pow(q,K+1));
      rad2=abs(pow(q,(K+2)*(K+1)/2)*y*pow(q,K+1)/(b2*(pow(1-x*pow(q,K+1),2)+pow(q,2*(K+1))*y*y))).upper();
      imag=a2+rad2;
      res=complex<interval<T> >(real,imag);
      return res;
    }
  }
  
  template <class T> interval<T> infinite_qPochhammer(const interval<T>& z,const interval<T>& q){
  // reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
  // Olshanetsky, Rogov 1995
    interval<T>res;
    res=Euler(interval<T>(q))/Karpelevich(interval<T>(z),interval<T>(q));
    return res;
  }

  template <class T> complex<interval<T> >infinite_qPochhammer(const complex< interval<T> >& z,const interval<T>& q){
  // reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
  // Olshanetsky, Rogov 1995
  complex<interval<T> >res;
  res=Euler(interval<T>(q))/Karpelevich(complex<interval<T> >(z),interval<T>(q));
  return res;
  }
  
  template <class T> interval<T> q_exp(const interval<T>& z,const interval<T>& q){
    // verification program for q-exponential function e_q(z)
    interval<T>res;
    res=1/infinite_qPochhammer(interval<T> (z),interval<T>(q));
    return res;
  }
  
  template <class T> complex<interval<T> >q_exp(const complex<interval<T> >& z,const interval<T>& q){
    // verification program for q-exponential function e_q(z)
    complex<interval<T> >res;
    res=1/infinite_qPochhammer(complex<interval<T> >(z),interval<T>(q));
    return res;
  }

  template <class T> interval<T> quantum_dilogarithm(const interval<T>& z,const interval<T>& q){
    // verification program for quantum dilogarithm
    interval<T>res;
    res=-log(infinite_qPochhammer(interval<T> (z),interval<T>(q)));
    return res;
  }

      template <class T> complex<interval<T> > quantum_dilogarithm (const complex<interval<T> >& z,const interval<T>& q){
    // verification program for quantum dilogarithm
      complex<interval<T> >res;
      res=-log(infinite_qPochhammer(complex<interval<T> >(z),interval<T>(q)));
    return res;
  }

  template <class T> interval<T> qAiry(const interval<T>& z,const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T> (-z*q), interval<T>(q*q));
    return res;
  }

  template <class T> interval<T> qPochhammer(const interval<T>& z,const interval<T>& q,const int& n){
    interval<T>res,qp;
    int k;
    qp=1.;
    if(n==0){
      res=1.;
      return res;
    }
    else {
      for(k=0;k<=n-1;k++){
	qp*=1-z*pow(q,k);
      }
      res=qp;
      return res;
    }
    }
  template <class T> complex<interval<T> >qPochhammer(const complex<interval<T> >& z,const interval<T>& q,const int& n){
    complex<interval<T> >res,qp;
    int k;
    qp=1.;
    if(n==0){
      res=1.;
      return res;
    }
    else {
      for(k=0;k<=n-1;k++){
	qp*=1-z*pow(q,k);
      }
      res=qp;
      return res;
    }
    }
}
#endif
