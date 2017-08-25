// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp
 
#ifndef QGAMMA_HPP
#define QGAMMA_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/complex.hpp>
#include <kv/constants.hpp>
#include <limits>
#include <algorithm>
#include <kv/Heine.hpp>
#include <kv/Pochhammer.hpp>

namespace kv{
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
     /*if(abs(res).upper()==std::numeric_limits<T>::infinity()){
       // Use asymptotic expansion 
       // M Mansour (2006) An asymptotic expansion of the q-gamma function Γ q (x), Journal of Nonlinear Mathematical Physics, 13:4, 479-483, DOI: 10.2991/jnmp.2006.13.4.2
       res=sqrt(1+q)*pow(1-q,0.5-z)*Euler(interval<T>(q*q))*pow(1-q*q,0.5)/infinite_qPochhammer(interval<T>(pow(q*q,0.5)),interval<T>(q*q))*interval<T>(1.,(exp(pow(q,z)/(1-q-pow(q,z)))).upper());
       }*/

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
  template <class T> complex<interval<T> >qgamma_Gauss_multi(const complex<interval<T> >& z,const interval<T>& q, const int& p){
  if(q<1 && q>0){
    // M Mansour (2006) An asymptotic expansion of the q-gamma function Γ q (x), Journal of Nonlinear Mathematical Physics, 13:4, 479-483, DOI: 10.2991/jnmp.2006.13.4.2
    // G Gasper , M Rahman, Basic Hypergeometric Series 2nd Edition, Cambridge University Press, 2004.
    interval<T> pq,pro2;
    pq=(1-pow(q,p))/(1-q);//pq OK
    complex<interval<T> >res,pro1;    
    pro1=1.;
    pro2=1.;   
    for(int i=0;i<=p-1;i++){
      pro1=pro1*q_gamma(complex<interval<T> >((z+i)/p),interval<T>(pow(q,p)));
      // pro1 OK
    }
    for(int j=1;j<=p-1;j++){
      interval<T> jj;
      jj=j;
      pro2=pro2*q_gamma(interval<T> (jj/p),interval<T>(pow(q,p)));      
    }

    res=pro1*pow(pq,z-1)/pro2;
    return res;
  }  
  else{
    throw std::domain_error("implemented for 0<q<1");
  }
}
template <class T> complex<interval<T> >qgamma_Legendre(const complex<interval<T> >& z,const interval<T>& q){
  if(q<1 && q>0){
    interval<T>qg;
    qg=q_gamma(interval<T>(0.5),interval<T>(q*q));
    complex<interval<T> >res;    
    
    res=q_gamma(complex<interval<T> >(z*0.5),interval<T>(q*q))
      *q_gamma(complex<interval<T> >((z+1)*0.5),interval<T>(q*q))
      *pow(1+q,z-1)/qg;         
  
    return res;
  }  
  else{
    throw std::domain_error("implemented for 0<q<1");
  }
}

  template <class T> interval<T> q_digamma(const interval<T>& x,const interval<T>& q){
    // q,x must be positive
   // verification program for q-digamma function
   // Reference: Kamel Brahim (2009), Turan-Type Inequalities for some q-Special Functions
   // Journal of inequalities in pure and applied mathematics, Volume 10
   interval<T>res,sum,qq,first,ratio;
   T rad;
   int N=100;
   sum=0.;
   qq=1.;
    if (q>=1){
     throw std::domain_error("value of q must be under 1");
   }
   if (q<=0){
     throw std::domain_error("q must be positive");
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
  template <class T> complex<interval<T> >incomplete_q_gamma(const complex<interval<T> >& z,const complex<interval<T> >& a,const interval<T>& q){
    // verification program for incomplete q-gamma function
    // expansion formula is used
    // reference
    // Ahmed Salem, A q-analogue of the exponential integral, 2013
    // warning: "a" should neither be negative integer nor zero
    complex<interval<T> >res,qq;
    qq=q;
    res=pow(z*(1-q),a)*q_gamma(complex<interval<T> >(a),interval<T>(q))
      *Heine(complex<interval<T> >(z*(1-q)),complex<interval<T> >(pow(q,a)),complex<interval<T> >(0.),interval<T>(q),complex<interval<T> >(qq));
    return res;
  }
}
#endif
