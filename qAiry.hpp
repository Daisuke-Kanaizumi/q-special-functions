// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp
// Date: March 4th, 2017

// verification program for q-Airy functions

#ifndef QAIRY_HPP
#define QAIRY_HPP
 
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/constants.hpp>
#include <cmath>
#include <limits>
#include <kv/convert.hpp> // this was included to use complex numbers
#include <kv/complex.hpp>
namespace kv {

  // calculate basic hypergeometric series 1phi1
  
  // reference
  // Fredrik Johansson, Computing hypergeometric functions rigorously, arXiv, 2016

template <class T> interval<T> _1phi_1(const interval<T>& a, const interval<T>& c,const interval<T>& q, const interval<T> & z) {
  int N;
    N=1000;
    interval<T> mid,rad,res,ratio,first;
    mid=1.;
   while(abs(c)>pow(1/q,N)){
     N=N+500;
     //throw std::domain_error("value of N not large enough");
    }   if (q>=1){
     throw std::domain_error("value of q must be under 1");
   }
   if (q<=0){
     throw std::domain_error("q must be positive");
   }
    for(int n=1;n<=N-1;n++){
      mid=mid+qPochhammer(interval<T>(a),interval<T>(q),int (n))*pow(-z,n)*pow(q,n*(n-1)/2)
	/qPochhammer(interval<T>(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=abs(qPochhammer(interval<T>(a),interval<T>(q),int (N))*pow(-z,N)*pow(q,N*(N-1)/2)
	      /qPochhammer(interval<T>(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
    ratio=abs(z)*abs(pow(q,N))*(1+abs(pow(q,N)*(c-a)/(1-c*pow(q,N))))*abs(1/(1-pow(q,N+1)));
    if(ratio<1){
      rad=(first/(1-ratio)).upper();
      res=mid+rad*interval<T>(-1.,1.);
      return res;
    }
    else{
      throw std::domain_error("ratio is more than 1");
    } 
 }
  template <class T> complex<interval<T> >_1phi_1(const complex<interval<T> >& a, const complex<interval<T> >& c,const interval<T>& q, const complex<interval<T> >& z) {
     int N;
    N=1000;
  while(abs(c)>pow(1/q,N)){
     N=N+500;
     //throw std::domain_error("value of N not large enough");
    }
   if (q>=1){
     throw std::domain_error("value of q must be under 1");
   }
   if (q<=0){
     throw std::domain_error("q must be positive");
   }
    complex<interval<T> > mid,res;
    interval<T>ratio,first;
    T rad;
    mid=1.;
    for(int n=1;n<=N-1;n++){
      mid=mid+qPochhammer(complex<interval<T> >(a),interval<T>(q) ,int (n))*pow(-z,n)*pow(q,n*(n-1)/2)
	/qPochhammer(complex<interval<T> >(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=abs(qPochhammer(complex<interval<T> >(a),interval<T>(q),int (N))*pow(-z,N)*pow(q,N*(N-1)/2)
	      /qPochhammer(complex<interval<T> >(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
    ratio=abs(z)*pow(q,N)*(1+abs(pow(q,N)*(c-a)/(1-c*pow(q,N))))*abs(1/(1-pow(q,N+1)));
    if(abs(ratio)<1){
      rad=(first/(1-ratio)).upper();
      res=complex_nbd(mid,rad);
      return res;
    }
    else{
      throw std::domain_error("ratio is more than 1");
    } }

  // calculate basic hypergeometric series 0phi1
  
  // reference
  // Fredrik Johansson, Computing hypergeometric functions rigorously, arXiv, 2016


template <class T> interval<T> _0phi_1(const interval<T>& c,const interval<T>& q, const interval<T> & z) {
    int N;
    N=100;
    interval<T> mid,res,ratio,first;
    T rad;
  while(abs(c)>pow(1/q,N)){
     N=N+500;
     //throw std::domain_error("value of N not large enough");
    }
   if (q>=1){
     throw std::domain_error("value of q must be under 1");
   }
   if (q<=0){
     throw std::domain_error("q must be positive");
   }
   mid=1.;
    for(int n=1;n<=N-1;n++){
      mid=mid+pow(z,n)*pow(q,n*(n-1))
	/qPochhammer(interval<T>(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=abs(pow(z,N)*pow(q,N*(N-1))
	      /qPochhammer(interval<T>(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
    ratio=abs(z)*pow(q,2*N)*abs(1/(1-c*pow(q,N)))*abs(1/(1-pow(q,N+1)));
    if(ratio<1){
      rad=(first/(1-ratio)).upper();
      res=mid+rad*interval<T>(-1.,1.);
      //std::cout<<first<<std::endl;
      //std::cout<<ratio<<std::endl;
      //std::cout<<rad<<std::endl;
      return res;
    }
    else{
      throw std::domain_error("ratio is more than 1");
    } 
 }
  template <class T> complex<interval<T> >_0phi_1(const complex<interval<T> >& c,const interval<T>& q, const complex<interval<T> >& z) {
     int N;
    N=4000;
    complex<interval<T> > mid,res;
    interval<T>ratio,first;
    T rad;
    mid=1.;
  while(abs(c)>pow(1/q,N)){
     N=N+500;
     //throw std::domain_error("value of N not large enough");
    }
   if (q>=1){
     throw std::domain_error("value of q must be under 1");
   }
   if (q<=0){
     throw std::domain_error("q must be positive");
   }
    for(int n=1;n<=N-1;n++){
      mid=mid+pow(z,n)*pow(q,n*(n-1))
	/qPochhammer(complex<interval<T> >(c),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
    }
    first=abs(pow(z,N)*pow(q,N*(N-1))
	      /qPochhammer(complex<interval<T> >(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
    ratio=abs(z)*pow(q,2*N)*abs(1/(1-c*pow(q,N)))*abs(1/(1-pow(q,N+1)));
    if(abs(ratio)<1){
      rad=(first/(1-ratio)).upper();
      res=complex_nbd(mid,rad);
      return res;
    }
    else{
      throw std::domain_error("ratio is more than 1");
    } 
 }
 template <class T> complex<interval<T> >HKW_qAiry(const interval<T>& q, const complex<interval<T> >& z) {
   // calculate Hamamoto-Kajiwara-Witte q-Airy function
   // reference
   // Hamamoto, Kajiwara, Witte, Hypergeometric solutions to the q-Painlev\'e equation of type $(A_1+A^{\prime}_1)~{(1)}$
   // International Mathematics Research Notices, 2006
   complex<interval<T> >res,qq;
   qq=q;
   res=_1phi_1(complex<interval<T> >(0),complex<interval<T> >(-qq),interval<T>(q),complex<interval<T> >(-z));
   if((abs(res)).upper()==std::numeric_limits<T>::infinity()){
     // implementation to avoid numerical divergence
     // reference
     // T. H. Koornwinder and R. F. Swarttouw, On $q$-analogues of the Fourier and Hankel transforms
     // Transactions of the American Mathematical Society, 333 (1992), 445-461
     // symmetry relation is used
     res=infinite_qPochhammer(complex<interval<T> >(-z),interval<T> (q))*_1phi_1(complex<interval<T> >(0),complex<interval<T> >(-z),interval<T>(q),complex<interval<T> >(-qq))/infinite_qPochhammer(interval<T>(-q),interval<T>(q));
   }
   return res;
 }
  template <class T> complex<interval<T> >KMNOY_qAiry(const interval<T>& q, const complex<interval<T> >& z) {
    // calculate Kajiwara-Masuda-Noumi-Ohta-Yamada q-Airy function
    // reference
    // Kajiwara, Masuda, Noumi, Ohta, Yamada, Hypergeometric Solutions to the q-Painlev\'e Equations
    // International Mathematics Research Notices, 2004
    complex<interval<T> >res,qq;
    qq=q;
    res=_1phi_1(complex<interval<T> >(0),complex<interval<T> >(-qq),interval<T>(q),complex<interval<T> >(-z*pow(q,0.5)));
   if((abs(res)).upper()==std::numeric_limits<T>::infinity()){
     // implementation to avoid numerical divergence
     // reference
     // T. H. Koornwinder and R. F. Swarttouw, On $q$-analogues of the Fourier and Hankel transforms
     // Transactions of the American Mathematical Society, 333 (1992), 445-461
     // symmetry relation is used
     res=infinite_qPochhammer(complex<interval<T> >(-z*pow(qq,0.5)),interval<T> (q))*_1phi_1(complex<interval<T> >(0),complex<interval<T> >(-z*pow(qq,0.5)),interval<T>(q),complex<interval<T> >(-qq))/infinite_qPochhammer(interval<T>(-q),interval<T>(q));
   }
    return res;
 }
template <class T> interval<T> Ramanujan_qAiry(const interval<T>& q, const interval<T> & z) {
  // calculate Ramanujan`s q-Airy function
  interval<T> res;
  if(z>0){
    // reference: ENCYCLOPEDIA OF MATHEMATICS AND ITS APPLICATOINS 71, SPECIAL FUNCTIONS
    // G.E.ANDREWS, RICHARD ASKEY, RANJAN ROY
    // CAMBRIDGE, 1999, Page 551, Exercise 39
    int K=1000;
    int j=1;
    interval<T>a,b,c,series;
    a=1.;
    b=1.;
    c=1.;
    // Leibniz criterion
    for (int k=1; k<=K; k++){
    j = -1*j;
    b=b*(1-pow(q,2*k));
    c=c*(1-z*pow(q,2*k));
    a = a+j*pow(q,k*k)*pow(z,k)/(b*c);
    }
    b=b*(1-pow(q,2*K+2));
    c=c*(1-z*pow(q,2*K+2));
    T rad=abs(pow(q,(K+1)*(K+1))*pow(z,K+1)/(b*c)).upper();
    series=a+rad*interval<T>(-1.,1.);
    res=series*infinite_qPochhammer(interval<T>(z*q*q),interval<T>(q*q));
  }
  else{
   res=_0phi_1(interval<T> (0),interval<T>(q),interval<T> (-q*z));
  }
  if((abs(res)).upper()==std::numeric_limits<T>::infinity()){
    // implementation to avoid numerical divergence
    // reference
    // Mourad E. H. Ismail, Changgui Zhang, Zeros of entire functions and a problem of Ramanujan, Theorem 3.1
    // Advances in Mathematics 209 (2007) 363–380   
    
    // In this paper, the proof of theorem 3.1 is partially abbreviated.
    // By using the following identities, the proof will be completed.
    // Gasper and Rahman, Basic Hypergeometric Series, Cambridge University Press, 1990, Appendix I, (I.2), (I.5), (I.17), (I.25)
    // Weisstein, Eric W. "q-Pochhammer Symbol." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/q-PochhammerSymbol.html, formula (10)     
    
    res=(infinite_qPochhammer(interval<T>(q*z),interval<T>(q*q))*infinite_qPochhammer(interval<T>(q/z),interval<T>(q*q))
	 *_1phi_1(interval<T>(0.),interval<T>(q),interval<T>(q*q),interval<T>(pow(q,2)/z))-
	 q*infinite_qPochhammer(interval<T>(q*q*z),interval<T>(q*q))*infinite_qPochhammer(interval<T>(1/z),interval<T>(q*q))/(1-q)
	 *_1phi_1(interval<T>(0.),interval<T>(q*q*q),interval<T>(q*q),interval<T>(pow(q,3)/z)))
      /infinite_qPochhammer(interval<T>(q),interval<T>(q*q));
  }
   return res;
 }
template <class T> complex<interval<T> >Ramanujan_qAiry(const interval<T>& q, const complex<interval<T> >& z) {
  // calculate Ramanujan`s q-Airy function
  complex<interval<T> >res,qq;
   qq=q;
   res=_0phi_1(complex<interval<T> >(0),interval<T>(q),complex<interval<T> >(-qq*z));
 if((abs(res)).upper()==std::numeric_limits<T>::infinity()){
   // implementation to avoid numerical divergence
   // reference
   // Mourad E. H. Ismail, Changgui Zhang, Zeros of entire functions and a problem of Ramanujan, Theorem 3.1
   // Advances in Mathematics 209 (2007) 363–380
   
   // In this paper, the proof of theorem 3.1 is partially abbreviated.
   // By using the following identities, the proof will be completed.
   // Gasper and Rahman, Basic Hypergeometric Series, Cambridge University Press, 1990, Appendix I, (I.2), (I.5), (I.17), (I.25)
   // Weisstein, Eric W. "q-Pochhammer Symbol." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/q-PochhammerSymbol.html, formula (10)     
   res=(infinite_qPochhammer(complex<interval<T> >(qq*z),interval<T>(q*q))*infinite_qPochhammer(complex<interval<T> >(qq/z),interval<T>(q*q))
	*_1phi_1(complex<interval<T> >(0.),complex<interval<T> >(qq),interval<T>(q*q),complex<interval<T> >(pow(qq,2)/z))-
	q*infinite_qPochhammer(complex<interval<T> >(qq*qq*z),interval<T>(q*q))*infinite_qPochhammer(complex<interval<T> >(1/z),interval<T>(q*q))/(1-q)
	*_1phi_1(complex<interval<T> >(0.),complex<interval<T> >(qq*qq*qq),interval<T>(q*q),complex<interval<T> >(pow(qq,3)/z)))
      /infinite_qPochhammer(interval<T>(q),interval<T>(q*q));
  }
   return res;
 }
  template <class T> interval<T> IKMMS_qAiry(const interval<T>& q, const interval<T>& tau){
    // calculate Isojima-Konno-Mimura-Murata-Satsuma q-Airy function
    // reference
    // Isojima, Konno, Mimura, Murata, Satsuma, Ultradiscrete Painlev\`e II equation and a special function solution
    // JOURNAL OF PHYSICS A: MATHEMATICAL AND THEORETICAL, 2011
   if (q>=1){
     throw std::domain_error("value of q must be under 1");
   }
   if (q<=0){
     throw std::domain_error("q must be positive");
   }
   if (tau<=0){
     throw std::domain_error("tau must be positive");
   }
   int N;
   N=1000;
   interval<T> mid,res,ratio,first,pi,qq,pt;
   T rad;
   qq=1.;
   pt=1.;
   pi=constants<interval<T> >::pi();
   mid=sqrt(2)*sin(pi*(2*log(tau)/log(q)+1)*0.25);
  
  for(int n=1;n<=N-1;n++){
    qq*=q*q;
    pt*=tau;
    mid=mid+sqrt(2)*sin(pi*(2*n+2*log(tau)/log(q)+1)*0.25)*pow(q,n*(n+1)/2)*tau/(1-qq);
  }
  qq*=q*q;
  pt*=tau;
  first=abs(sqrt(2)*sin(pi*(2*N+2*log(tau)/log(q)+1)*0.25)*pow(q,N*(N+1)/2)*tau/(1-qq));
  ratio=abs(tau)*pow(q,N+1)/(1-pow(q,2*N+2));
  // ratio of majorant
  if(ratio<1){
    rad=(first/(1-ratio)).upper();
    res=mid+rad*interval<T>(-1.,1.);
    return res;
  }
  else{
    throw std::domain_error("ratio is more than 1");
  } 
  }
template <class T> interval<T> q_Bi(const interval<T>& q, const interval<T>& tau){
    // calculate q-Bi function
    // reference
    // Isojima, Konno, Mimura, Murata, Satsuma, Ultradiscrete Painlev\`e II equation and a special function solution
    // JOURNAL OF PHYSICS A: MATHEMATICAL AND THEORETICAL, 2011
   if (q>=1){
     throw std::domain_error("value of q must be under 1");
   }
   if (q<=0){
     throw std::domain_error("q must be positive");
   }
   if (tau<=0){
     throw std::domain_error("tau must be positive");
   }

   int N;
   N=1000;
   interval<T> mid,res,ratio,first,pi,qq,pt;
   T rad;
   qq=1.;
   pt=1.;
   pi=constants<interval<T> >::pi();
   mid=sqrt(2)*sin(pi*(2*log(tau)/log(q)+3)*0.25);
   for(int n=1;n<=N-1;n++){
    qq*=q*q;
    pt*=tau;
    mid=mid+sqrt(2)*sin(pi*(2*n+2*log(tau)/log(q)+3)*0.25)*pow(q,n*(n+1)/2)*tau/(1-qq);
   }
   qq*=q*q;
   pt*=tau;
   first=abs(sqrt(2)*sin(pi*(2*N+2*log(tau)/log(q)+3)*0.25)*pow(q,N*(N+1)/2)*tau/(1-qq));
   ratio=abs(tau)*pow(q,N+1)/(1-pow(q,2*N+2));
  // ratio of majorant
   if(ratio<1){
     rad=(first/(1-ratio)).upper();
     res=mid+rad*interval<T>(-1.,1.);
     return res;
   }
   else{
     throw std::domain_error("ratio is more than 1");
   } 
}
}
#endif
