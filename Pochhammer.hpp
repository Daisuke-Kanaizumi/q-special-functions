// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp

// This header file includes following functions
// Pochhammer symbol
// Euler function (q;q)_{\infty}
// finite and infinite q-Pochhammer symbols
// q-exponential function e_q(z)
// q-sin, q-cos
// Jackson and Moak q-gamma, q-beta, q-digamma, symmetric q-gamma, symmetric q-beta (Brahim-Sidomou, 2013)
// quantum dilogarithm, quantum polylogarithm (Kirillov, 1994)
// Ramanujan theta, Ramanujan psi sum

#ifndef POCHHAMMER_HPP
#define POCHHAMMER_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/complex.hpp>
#include <kv/constants.hpp>
#include <kv/dilogarithm.hpp>
#include <kv/Heine.hpp>
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
     throw std::domain_error("absolute value of q must be under 1");
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
      throw std::domain_error("absolute value of q must be under 1");
    }
    if (abs(z)>=1){
      throw std::domain_error("absolute value of z must be under 1");
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
    complex<interval<T> >res,logqp,logmid,logres;
    interval<T> r,pi;
    pi=constants<interval<T> >::pi();
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
 /*
 // abolished implementation (March 28th, 2017)
    if(abs(z)<1 && abs(arg(1-z))<pi && z.imag()!=0){
      // Reference: Thomas Prellberg (1995), Uniform q-series asymptotics for staircase polygons
      // Journal of Physics A: Mathematical and General
      logmid=dilogarithm(complex<interval<T> >(z))/log(q)+log(1-z)/2;
      // dilogarithm OK
      rad=(abs(log(q))*(log(abs(1-z))+atan(z.imag()/(1-z.real()))*z.real()/z.imag())/6).upper();
      logres=complex_nbd(logmid,rad);
      res=exp(logres);
      return res;   
    }
 */
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
  /*
template <class T> interval<T> infinite_qPochhammer(const interval<T>& z,const interval<T>& q){
  // abolished implementation (March 22nd, 2017)
  // reference: Ismail,Zhang (2016), Integral and Series Representations of q-Polynomials and Functions: Part I, arXiv
  // Chapter 6.1, Lemma 6.1
  interval<T>res;
  T rad;
  rad=(exp(abs(z)/(1-q))*abs(z)/(1-q)).upper();
  res=1+rad*interval<T>(-1,1);
  return res;
}
  */
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
  template <class T> interval<T> quantum_polylogarithm(const int& k, const interval<T>& z,const interval<T>& q){
    if (abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    if (abs(z)>=1){
      throw std::domain_error("absolute value of z must be under 1");
    }
    // verification program for quantum polylogarithm
    interval<T>res,mid,rad,ratio,first,power,qq,denom;
    if(k==2){
      res=quantum_dilogarithm(interval<T>(z),interval<T>(q));
    }
    if(k>=3){
      int N;
      N=1000;
      mid=z/pow(1-q,k-1);
      power=z;
      denom=1.;
      qq=q;
      for(int n=2;n<=N-1;n++){
	power=power*z;
	qq=qq*q;
	denom=denom+1;
	mid=mid+power/denom/pow(1-qq,k-1);
      }
      power=power*z;
      qq=qq*q;
      denom=denom+1;
      first=abs(power/denom/pow(1-qq,k-1));
      ratio=abs(z)*N*pow((1-qq)/(1-qq*q),k-1)/(N+1);
      if(ratio<1){
	rad=(first/(1-ratio)).upper();
	res=mid+rad*interval<T>(-1.,1.);
      }
      if (ratio>=1){
	throw std::domain_error("ratio must be under 1");
      }
    }
    return res;
  }
  template <class T> complex<interval<T> >quantum_polylogarithm(const int& k, const complex<interval<T> >& z,const interval<T>& q){
    if (abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    if (abs(z)>=1){
      throw std::domain_error("absolute value of z must be under 1");
    }
    // verification program for quantum polylogarithm
    complex<interval<T> >res,mid,power;
    interval<T>qq, denom,ratio,first;
    T rad;
    if(k==2){
      res=quantum_dilogarithm(complex<interval<T> >(z) ,interval<T>(q));
    }
    if(k>=3){
      int N;
      N=1000;
      mid=z/pow(1-q,k-1);
      power=z;
      denom=1.;
      qq=q;
      for(int n=2;n<=N-1;n++){
	power=power*z;
	qq=qq*q;
	denom=denom+1;
	mid=mid+power/denom/pow(1-qq,k-1);
      }
      power=power*z;
      qq=qq*q;      
      first=abs(power/N/pow(1-qq,k-1));
      ratio=abs(z)*N*pow((1-qq)/(1-qq*q),k-1)/(N+1);
      if(abs(ratio)<1){
	rad=(first/(1-ratio)).upper();
	res=complex_nbd(mid,rad);
      }
      if (abs(ratio)>=1){
	throw std::domain_error("ratio must be under 1");
      }
    }
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
 } template <class T> interval<T> q_digamma(const interval<T>& x,const interval<T>& q){
   // q,x must be positive
   // verification program for q-digamma function
   // Reference: Kamel Brahim (2009), Turan-Type Inequalities for some q-Special Functions
   // Journal of inequalities in pure and applied mathematics, Volume 10
   interval<T>res,sum,qq,first,ratio;
   T rad;
   int N=100;
   sum=0.;
   qq=1.;
   if (abs(q)>=1){
     throw std::domain_error("absolute value of q must be under 1");
   }
   if (x<=0){
     throw std::domain_error("implemented for positive x");
   }
   for(int n=1;n<=N-1;n++){
     qq=qq*q;
     sum=sum+pow(q,n*x)/(1-qq);
   }
   qq=qq*q;
   first=pow(q,N*x)/(1-qq);
   ratio=(1-qq)*pow(q,x)/(1-qq*q);
 if(abs(ratio)<1){
      rad=(first/(1-ratio)).upper();
      res=-log(1-q)+log(q)*(sum+rad*interval<T>(-1.,1.));
      return res;
    }
    else{
      std::cout<<"ratio is more than 1"<<std::endl;
    } 
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
  template <class T> interval<T>Ramanujan_theta(const interval<T>& a,const interval<T>& b){
    // verification program for Ramanujan theta function
    interval<T>res;
    if (abs(a*b)>=1){
      throw std::domain_error("absolute value of a*b must be under 1");
    }
    res=infinite_qPochhammer(interval<T>(-a),interval<T>(a*b))*infinite_qPochhammer(interval<T>(-b),interval<T>(a*b))*Euler(interval<T>(a*b));
    return res;
  }
  template <class T> complex<interval<T> >Ramanujan_psi_sum(const complex<interval<T> >& a,const complex<interval<T> >& b,const interval<T>(q),const complex<interval<T> >& z){
    // verification program for Ramanujan psi sum
    if(abs(b/a)<abs(z)&&abs(z)<1&&q<1){
    complex<interval<T> >res;
    res=infinite_qPochhammer(complex<interval<T> >(a*z),interval<T>(q))*infinite_qPochhammer(complex<interval<T> >(q/(a*z)),interval<T>(q))
      *Euler(interval<T>(q))*infinite_qPochhammer(complex<interval<T> >(b/a),interval<T>(q))
      /infinite_qPochhammer(complex<interval<T> >(z),interval<T>(q))/infinite_qPochhammer(complex<interval<T> >(b/(a*z)),interval<T>(q))
      /infinite_qPochhammer(complex<interval<T> >(b),interval<T>(q))/infinite_qPochhammer(complex<interval<T> >(q/a),interval<T>(q));
    return res;
    }
    else{
      throw std::domain_error("Ramanujan psi sum cannot be calculated");
    }
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
