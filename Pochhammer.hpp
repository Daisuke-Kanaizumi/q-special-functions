// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp

// This header file includes following functions
// Pochhammer symbol
// Euler function (q;q)_{\infty}
// finite and infinite q-Pochhammer symbols
// q-exponential function e_q(z)
// q-sin, q-cos
// Jackson and Moak q-gamma, q-beta, symmetric q-gamma, symmetric q-beta (Brahim-Sidomou, 2013)
// quantum dilogarithm (Kirillov, 1994)

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
   // Verify Euler function by using the pentagonal number theorem
   // Reference:The truncated pentagonal number theorem
   // Andrews, Merca, 2012
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
     res=a+rad*interval<T>(-1.,1.);
     return res;
   }
 }
  template <class T> interval<T> Karpelevich(const interval<T>& z,const interval<T>& q){
  // Reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
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
      res=a+rad*interval<T>(-1.,1.);
      return res;
    }
  }
  template <class T> complex<interval<T> > Karpelevich(const complex< interval<T> >& z,const interval<T>& q){
  // Reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
  // Olshanetsky, Rogov 1995
    // available when re(z)>0,0<q<1,|z|<1
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
      a1=(1-x)/(pow(1-x,2)+y*y);
      b1=1.;
      // derived from the definition of q-Pochhammer symbol
      // Leibniz criterion
      for (k=1; k<=K; k++){
	j = -1*j;
	b1=b1*(1-pow(q,k));
	a1=a1+j*pow(q,k*(k+1)/2)*(1-x*pow(q,k))/(b1*(pow(1-x*pow(q,k),2)+pow(q,2*k)*y*y));
      }
      b1=b1*(1-pow(q,K+1));
      rad1=abs(pow(q,(K+2)*(K+1)/2)*(1-x*pow(q,K+1))/(b1*(pow(1-x*pow(q,K+1),2)+pow(q,2*(K+1))*y*y))).upper();
      real=a1+rad1*interval<T>(-1.,1.);
      
      a2=y/(pow(1-x,2)+y*y);
      b2=1.;
      // derived from the definition of q-Pochhammer symbol
      // Leibniz criterion
      for (kk=1; kk<=K; kk++){
	jj = -1*jj;
	b2=b2*(1-pow(q,kk));
		 a2=a2+jj*pow(q,kk*(kk+1)/2)*y*pow(q,kk)/(b2*(pow(1-x*pow(q,kk),2)+pow(q,2*kk)*y*y));
      }
      b2=b2*(1-pow(q,K+1));
      rad2=abs(pow(q,(K+2)*(K+1)/2)*y*pow(q,K+1)/(b2*(pow(1-x*pow(q,K+1),2)+pow(q,2*(K+1))*y*y))).upper();
      imag=a2+rad2*interval<T>(-1.,1.);
      res=complex<interval<T> >(real,imag);
      return res;
    
  }
template <class T> interval<T> infinite_qPochhammer(const interval<T>& z,const interval<T>& q){
  // Revised implementation of infinite q-Pochhammer symbols
  // March 9th, 2017

  // 0<q<1
  interval<T> res,r;
  T rad;
  int n;
  n=100;
  if(abs(z)<1){
    // reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
    // Olshanetsky, Rogov 1995
    // When absolute value of q is close to 1 (for example, q=0.999), zero division occurs
    res=Euler(interval<T>(q))/Karpelevich(interval<T>(z),interval<T>(q));
    return res;
  }
  else{
  // Reference: Plancherel-Rotach asymptotics for certain basic hypergeometric series
  // Zhang, 2008

  while(abs(z)*pow(q,n)/(1-q)>=0.5){
    n=n+500;
  }

  rad=(2*abs(z)*pow(q,n)/(1-q)).upper();
  r=interval<T>(1-rad,1+rad);
  res=qPochhammer(interval<T>(z),interval<T>(q),int(n))*r;
  return res;
  }
}

  template <class T> complex<interval<T> >infinite_qPochhammer(const complex<interval<T> >& z,const interval<T>& q){
    complex<interval<T> >res;
    interval<T> r;
    T rad;
    int n;
    n=100;
 if(abs(z)<1 && z.real()>0){
    // reference: The Modified q-Bessel Functions and the q-Bessel-Macdonald Functions
    // Olshanetsky, Rogov 1995
    // When absolute value of q is close to 1 (for example, q=0.999), zero division occurs
   res=Euler(interval<T>(q))/Karpelevich(complex<interval<T> >(z),interval<T>(q));
   return res;
  }
 else{
 // Reference: Plancherel-Rotach asymptotics for certain basic hypergeometric series
  // Zhang, 2008

  while(abs(z)*pow(q,n)/(1-q)>=0.5){
    n=n+500;
  }
  rad=(2*abs(z)*pow(q,n)/(1-q)).upper();
  r=interval<T>(1-rad,1+rad);
  res=qPochhammer(complex<interval<T> >(z),interval<T>(q),int(n))*r;
  return res;
 }
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

  template <class T> complex<interval<T> >q_sin(const complex<interval<T> >& z,const interval<T>& q){
    // verification program for q-sin function sin_q(z)
    complex<interval<T> >res,i;
    i=complex<interval<T> >::i();
    res=(q_exp(complex<interval<T> >(i*z),interval<T>(q))-q_exp(complex<interval<T> >(-i*z),interval<T>(q)))/(2*i);
    return res;
  }
  template <class T> complex<interval<T> >q_cos(const complex<interval<T> >& z,const interval<T>& q){
    // verification program for q-cos function cos_q(z)
    complex<interval<T> >res,i;
    i=complex<interval<T> >::i();
    res=(q_exp(complex<interval<T> >(i*z),interval<T>(q))+q_exp(complex<interval<T> >(-i*z),interval<T>(q)))/2;
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
 template <class T> interval<T> q_gamma(const interval<T>& z,const interval<T>& q){
   // q must be positive
   // verification program for q-gamma function
   interval<T>res;
   if(q<1 && q>0){
     if(pow(q,z)<1){
       res=pow(1-q,1-z)*Karpelevich(interval<T>(pow(q,z)),interval<T>(q));
     }
     else{
       res=Euler(interval<T>(q))*pow(1-q,1-z)/infinite_qPochhammer(interval<T>(pow(q,z)),interval<T>(q));
     }
   }
   if(q>1){ // Moak q-gamma function
     if(pow(q,-z)<1){
       res=pow(q-1,1-z)*pow(q,z*(z-1)/2)*Karpelevich(interval<T>(pow(q,-z)),interval<T>(1/q));
     }
     else{
       res=Euler(interval<T>(1/q))*pow(q-1,1-z)*pow(q,z*(z-1)/2)/infinite_qPochhammer(interval<T>(pow(q,-z)),interval<T>(1/q));
     }
   }
   return res;
 }
  template <class T> complex<interval<T> >q_gamma(const complex<interval<T> >& z,const interval<T>& q){
    complex<interval<T> >res;
    if(q<1 && q>0){
      res=Euler(interval<T>(q))*pow(1-q,1-z)/infinite_qPochhammer(complex<interval<T> >(pow(q,z)),interval<T>(q));
    }
    if(q>1){
      res=Euler(interval<T>(1/q))*pow(q-1,1-z)*pow(q,z*(z-1)/2)/infinite_qPochhammer(complex<interval<T> >(pow(q,-z)),interval<T>(1/q));
    }
    return res;
 }
 template <class T> interval<T> q_beta(const interval<T>& a,const interval<T>& b,const interval<T>& q){
   // q must be positive
   // verification program for q-beta function
   interval<T>res;
   res=q_gamma(interval<T>(a),interval<T>(q))*q_gamma(interval<T>(b),interval<T>(q))/q_gamma(interval<T>(a+b),interval<T>(q));
   return res;
 }
  template <class T> complex<interval<T> >q_beta(const complex<interval<T> >& a,const complex<interval<T> >& b,const interval<T>& q){
   // q must be positive
   // verification program for q-beta function
    complex<interval<T> >res;
    res=q_gamma(complex<interval<T> >(a),interval<T>(q))*q_gamma(complex<interval<T> >(b),interval<T>(q))/q_gamma(complex<interval<T> >(a+b),interval<T>(q));
    return res;
  }
 template <class T> interval<T> symmetric_q_gamma(const interval<T>& z,const interval<T>& q){
   // verification program for symmetric q-gamma function
   // reference
   // Brahim and Sidomou, On Some Symmetric q-Special Functions, 2013
   interval<T>res;
   res=pow(q,-(z-1)*(z-2)/2)*q_gamma(interval<T>(z),interval<T>(q*q));
   return res;
 }
 template <class T> interval<T> symmetric_q_beta(const interval<T>& a,const interval<T>& b,const interval<T>& q){
   // q,a,b must be positive
   // verification program for symmetric q-beta function
   // reference
   // Brahim and Sidomou, On Some Symmetric q-Special Functions, 2013
   interval<T>res;
   res=symmetric_q_gamma(interval<T>(a),interval<T>(q))*symmetric_q_gamma(interval<T>(b),interval<T>(q))/symmetric_q_gamma(interval<T>(a+b),interval<T>(q));
   return res;
 }
  template <class T> complex<interval<T> >symmetric_q_gamma(const complex<interval<T> >& z,const interval<T>& q){
   // verification program for symmetric q-gamma function
   // reference
   // Brahim and Sidomou, On Some Symmetric q-Special Functions, 2013
    complex<interval<T> >res;
    res=pow(q,-(z-1)*(z-2)/2)*q_gamma(complex<interval<T> >(z),interval<T>(q*q));
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
