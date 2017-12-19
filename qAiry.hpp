// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp
// Date: Dec 19th, 2017

// verification program for q-Airy functions

#ifndef QAIRY_HPP
#define QAIRY_HPP
 
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/constants.hpp>
#include <kv/Pochhammer.hpp>
#include <kv/qBessel.hpp>
#include <cmath>
#include <limits>
#include <kv/convert.hpp> // this was included to use complex numbers
#include <kv/complex.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ub = boost::numeric::ublas;
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
      // Alternative implementation
      // Reference
      // Koekoek-Swarttouw,The Askey-scheme of hypergeometric orthogonal polynomials and its q-analogue, arXiv (1996)
      // formula 0.6.12 and 0.6.13
        if((abs(res)).upper()==std::numeric_limits<T>::infinity()&&abs(a)<1&&abs(a)>0){
      res=infinite_qPochhammer(interval<T>(a),interval<T>(q))*infinite_qPochhammer(interval<T>(z),interval<T>(q))
	*Heine(interval<T>(c/a),interval<T>(0.),interval<T>(z),interval<T>(q),interval<T>(a))/infinite_qPochhammer(interval<T>(c),interval<T>(q));
	}
	if((abs(res)).upper()==std::numeric_limits<T>::infinity()&&abs(a*z/c)<1&&abs(a)>0){
	  res=infinite_qPochhammer(interval<T>(a*z/c),interval<T>(q))*Heine(interval<T>(c/a),interval<T>(0.),interval<T>(c),interval<T>(q),interval<T>(a*z/c));
	}

	return res;
    }
    else{
      throw std::domain_error("ratio is more than 1");
    } 
    
}template <class T> interval<T> _1phi_1_bs(const interval<T>& a, const interval<T>& c,const interval<T>& q, const interval<T> & z) {
  int N,nn;
  nn=10;
  N=std::pow(2,nn);
  interval<T> mid,rad,res,ratio,first;
  mid=0.;

  
    while(abs(c)>pow(1/q,N)){
     nn=nn+5;
     N=std::pow(2,nn);
     //throw std::domain_error("value of N not large enough");
    } 
    ub::vector<interval<T> >A(N),B(N),C(N);
    for(int i=0;i<=N/2-1;i++){
      A(i)=pow(-z,i)*pow(q,i*(i-1)*0.5)/qPochhammer(interval<T>(q),interval<T>(q),int (i));
    }
    if (q>=1){
     throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    for(int n=1;n<=N/2-1;n++){
      B(0)=qPochhammer(interval<T>(a),interval<T>(q),int (n));
      C(0)=qPochhammer(interval<T>(c),interval<T>(q),int (n));
      for(int j=1;j<=N/2-1;j++){
	B(j)=1.;
      }
      for(int k=1;k<=N/2-1;k++){
	C(k)=1.;
      }
      mid=mid+B(0)/C(0)*(A(2*n)*C(2*n+1)+B(2*n)*A(2*n+1));
	
    }
    first=abs(qPochhammer(interval<T>(a),interval<T>(q),int (N))*pow(-z,N)*pow(q,N*(N-1)/2)
	      /qPochhammer(interval<T>(c),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
    ratio=abs(z)*abs(pow(q,N))*(1+abs(pow(q,N)*(c-a)/(1-c*pow(q,N))))*abs(1/(1-pow(q,N+1)));
    if(ratio<1){
      rad=(first/(1-ratio)).upper();
      res=mid+rad*interval<T>(-1.,1.);
      // Alternative implementation
      // Reference
      // Koekoek-Swarttouw,The Askey-scheme of hypergeometric orthogonal polynomials and its q-analogue, arXiv (1996)
      // formula 0.6.12 and 0.6.13
        if((abs(res)).upper()==std::numeric_limits<T>::infinity()&&abs(a)<1&&abs(a)>0){
      res=infinite_qPochhammer(interval<T>(a),interval<T>(q))*infinite_qPochhammer(interval<T>(z),interval<T>(q))
	*Heine(interval<T>(c/a),interval<T>(0.),interval<T>(z),interval<T>(q),interval<T>(a))/infinite_qPochhammer(interval<T>(c),interval<T>(q));
	}
	if((abs(res)).upper()==std::numeric_limits<T>::infinity()&&abs(a*z/c)<1&&abs(a)>0){
	  res=infinite_qPochhammer(interval<T>(a*z/c),interval<T>(q))*Heine(interval<T>(c/a),interval<T>(0.),interval<T>(c),interval<T>(q),interval<T>(a*z/c));
	}
	if((abs(res)).upper()==std::numeric_limits<T>::infinity()&&a<1&&a>=0&&c>=0&&c<1){
  // Alternative implementation
  // Reference: Plancherel-Rotach asymptotics for certain basic hypergeometric series
  // Zhang, 2008, Advances in Mathematics
	  int nn,mu;
	  nn=5; 
	  mu=std::floor(nn/2);
	  interval<T>r,amid;
	  
	  T arad;
	  arad=(32*Euler(interval<T>(q))*infinite_qPochhammer(interval<T>(-q*pow(q,-nn)/abs(z)),interval<T>(q))
		*infinite_qPochhammer(interval<T>(-q*pow(q,nn)*abs(z)),interval<T>(q))
		/infinite_qPochhammer(interval<T>(a),interval<T>(q))
		*(pow(q,mu+1)/(1-q)+pow(q,mu*mu/2+mu/2)/pow(abs(z*pow(q,nn)),mu))).upper();
	  amid=Euler(interval<T>(q))*infinite_qPochhammer(interval<T>(q*pow(q,-nn)/abs(z)),interval<T>(q))
	    *infinite_qPochhammer(interval<T>(q*pow(q,nn)*abs(z)),interval<T>(q));
	  r=interval<T>(amid.lower()-arad,amid.upper()+arad);
	  res=Euler(interval<T>(q))*infinite_qPochhammer(interval<T>(c),interval<T>(q))*pow(-z*pow(q,nn),nn)*r
	    /infinite_qPochhammer(interval<T>(a),interval<T>(q))/pow(q,nn*(nn+1)/2);
	}  
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
       // Alternative implementation
       // Reference
       // Koekoek-Swarttouw,The Askey-scheme of hypergeometric orthogonal polynomials and its q-analogue, arXiv (1996)
       // formula 0.6.12 and 0.6.13
       if((abs(res)).upper()==std::numeric_limits<T>::infinity()&&abs(a)<1&&abs(a)>0){
	 res=infinite_qPochhammer(complex<interval<T> >(a),interval<T> (q))*infinite_qPochhammer(complex<interval<T> >(z),interval<T>(q))
	   *Heine(complex<interval<T> >(c/a),complex<interval<T> >(0.),complex<interval<T> >(z),interval<T>(q),complex<interval<T> >(a))/infinite_qPochhammer(complex<interval<T> >(c),interval<T>(q));
       }
       if((abs(res)).upper()==std::numeric_limits<T>::infinity()&&abs(a*z/c)<1&&abs(a)>0){
	 res=infinite_qPochhammer(complex<interval<T> >(a*z/c),interval<T>(q))*Heine(complex<interval<T> >(c/a),complex<interval<T> >(0.),complex<interval<T> >(c),interval<T>(q),complex<interval<T> >(a*z/c));
       }
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
    N=1000;
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
      if((abs(res)).upper()==std::numeric_limits<T>::infinity()&&abs(c)>0){
	// Alternative implementation
	// Reference
	// Koekoek-Swarttouw,The Askey-scheme of hypergeometric orthogonal polynomials and its q-analogue, arXiv (1996)
	// formula 0.6.8 and 0.6.9
	// H. T Koelink, Hansen-Lommel orthogonality Relations for Jackson`s q-Bessel functions, formula 3.2
	// Journal of Mathematical Analysis and Applications 175, 425-437 (1993)
	res=_1phi_1(interval<T>(z/c),interval<T>(0.),interval<T>(q),interval<T>(c))
	  /infinite_qPochhammer(interval<T>(c),interval<T>(q));
      }
    
      return res;
    }
    else{
      throw std::domain_error("ratio is more than 1");
    } 
 }
  template <class T> complex<interval<T> >_0phi_1(const complex<interval<T> >& c,const interval<T>& q, const complex<interval<T> >& z) {
    int N;
    N=1000;
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
      if((abs(res)).upper()==std::numeric_limits<T>::infinity()&&abs(c)>0){
	// Alternative implementation
	// Reference
	// Koekoek-Swarttouw,The Askey-scheme of hypergeometric orthogonal polynomials and its q-analogue, arXiv (1996)
	// formula 0.6.8 and 0.6.9
	// H. T Koelink, Hansen-Lommel orthogonality Relations for Jackson`s q-Bessel functions, formula 3.2
	// Journal of Mathematical Analysis and Applications 175, 425-437 (1993)
	res=_1phi_1(complex<interval<T> >(z/c),complex<interval<T> >(0.),interval<T>(q),complex<interval<T> >(c))
	  /infinite_qPochhammer(complex<interval<T> >(c),interval<T>(q));
 }
      
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
     res=infinite_qPochhammer(complex<interval<T> >(-z*pow(qq,0.5)),interval<T> (q))
       *_1phi_1(complex<interval<T> >(0),complex<interval<T> >(-z*pow(qq,0.5)),interval<T>(q),complex<interval<T> >(-qq))
       /infinite_qPochhammer(interval<T>(-q),interval<T>(q));
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
    complex<interval<T> >res,qq,i;
    qq=q;
    res=_0phi_1(complex<interval<T> >(0),interval<T>(q),complex<interval<T> >(-qq*z));
    if((abs(res)).upper()==std::numeric_limits<T>::infinity()){
      // implementation to avoid numerical divergence
      // reference
      // Mourad E. H. Ismail, Ruiming Zhang - Proceedings of the American Mathematical Society, 2017
      // Integral and series representations of $q$-polynomials and functions: Part II Schur polynomials and the Rogers-Ramanujan identities
     i=complex<interval<T> >::i(); 
     res=infinite_qPochhammer(complex<interval<T> >(-i*sqrt(z)*pow(q,0.25)),interval<T>(sqrt(q)))
       *_1phi_1(complex<interval<T> >(0.),complex<interval<T> >(-i*sqrt(z)*pow(q,0.25)),interval<T>(sqrt(q)),complex<interval<T> >(i*sqrt(z)*pow(q,0.25)));
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
      res=(infinite_qPochhammer(complex<interval<T> >(qq*z),interval<T>(q*q))*infinite_qPochhammer(complex<interval<T> >(qq/z),interval<T>(q*q))
	   *_1phi_1(complex<interval<T> >(0.),complex<interval<T> >(qq),interval<T>(q*q),complex<interval<T> >(pow(qq,2)/z))-
	   q*infinite_qPochhammer(complex<interval<T> >(qq*qq*z),interval<T>(q*q))*infinite_qPochhammer(complex<interval<T> >(1/z),interval<T>(q*q))/(1-q)
	   *_1phi_1(complex<interval<T> >(0.),complex<interval<T> >(qq*qq*qq),interval<T>(q*q),complex<interval<T> >(pow(qq,3)/z)))
	/infinite_qPochhammer(interval<T>(q),interval<T>(q*q));
       }
    return res;
  }

  template <class T> complex<interval<T> >Ramanujan_qAiry_Morita(const interval<T>& q, const complex<interval<T> >& z) {
    // Morita, T. (2011). A connection formula between the Ramanujan function and the $ q $-Airy function. arXiv preprint arXiv:1104.0755.
    complex<interval<T> >res,i,x;
    i=complex<interval<T> >::i(); 
    x=i*pow(q,0.75)/sqrt(z);
    res=(infinite_qPochhammer(complex<interval<T> >(-x/sqrt(q)),interval<T>(sqrt(q)))*infinite_qPochhammer(complex<interval<T> >(-q/x),interval<T>(sqrt(q)))*HKW_qAiry(interval<T>(sqrt(q)),complex<interval<T> >(-x))
	  +infinite_qPochhammer(complex<interval<T> >(x/sqrt(q)),interval<T>(sqrt(q)))*infinite_qPochhammer(complex<interval<T> >(q/x),interval<T>(sqrt(q)))*HKW_qAiry(interval<T>(sqrt(q)),complex<interval<T> >(x)))
      /infinite_qPochhammer(interval<T>(-1.),interval<T>(sqrt(q)));
    return res;
  }
  template <class T> complex<interval<T> >Ramanujan_qAiry_mqB(const interval<T>& q, const complex<interval<T> >& z) {
    complex<interval<T> >res;
    res=sqrt(q*z)/(1+q+q*q)
      *(modified_qBesselI2(complex<interval<T> >(2/(1+q+q*q)*pow(q*z,1.5)*(1-pow(q,1/3.))),complex<interval<T> >(-1/3.,0.),interval<T>(pow(q,1/3.)))-
	modified_qBesselI2(complex<interval<T> >(2/(1+q+q*q)*pow(q*z,1.5)*(1-pow(q,1/3.))),complex<interval<T> >(1/3.,0),interval<T>(pow(q,1/3.))));
    return res; 
    // to be updated
  }
  /*
  template <class T> complex<interval<T> >Ramanujan_qAiry_ae(const interval<T>& q, const complex<interval<T> >& Z) {
    // calculate Ramanujan`s q-Airy function with asymptotic formula
    // reference
    // Ismail, M. E., & Zhang, R. (2016).
    // Integral and Series Representations of $ q $-Polynomials and Functions: Part I.
    // arXiv preprint arXiv:1604.08441.
    complex<interval<T> >res,z,nbd;
    int n;
    n=10;
    
    z=Z*pow(q,2*n);
    
    res=pow(q,-n*n)*pow(-z,n)/Euler(interval<T>(q))*Euler(interval<T>(q*q))
      *infinite_qPochhammer(complex<interval<T> >(q*z),interval<T>(q*q))
      *infinite_qPochhammer(complex<interval<T> >(q/z),interval<T>(q*q));
    return res;
  }
  */
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
