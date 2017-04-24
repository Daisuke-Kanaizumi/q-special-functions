// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp

// verification program for q-Bessel functions
// March 9th, 2017

#ifndef QBESSEL_HPP
#define QBESSEL_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

namespace kv{
template <class T> interval<T>Jackson1(const interval<T>& z,const interval<T>& nu,const interval<T>& q){
// verification program for Jackson`s 1st q-Bessel function
interval<T>res,a,b,c,series;
int j,K;
T rad;
if(abs(z)>=2){
throw std::domain_error("Jackson`s 1st q-Bessel function is not defined");
}
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 if(z<0){
   throw std::domain_error("z must be positive");
 }
 if (pow(q,nu)<1){
j=1;
K=500;
a=1.;
b=1.;
c=1.;
// Leibniz criterion
 for (int k=1; k<=K; k++){
    j = -1*j;
    b=b*(1-pow(q,k));
    c=c*(1-pow(q,k+nu));
    a = a+j*pow(z/2,2*k)/(b*c);
  }
b=b*(1-pow(q,K+1));
c=c*(1-pow(q,K+nu+1));

rad=abs(pow(z/2,2*K+2)/(b*c)).upper();
series=a+rad*interval<T>(-1.,1.);
res=pow(z/2,nu)*series/Karpelevich(interval<T>(pow(q,nu+1)),interval<T>(q));
}
 else{
res=pow(z/2,nu)*infinite_qPochhammer(interval<T>(pow(q,nu+1)),interval<T>(q))*
  Heine(interval<T>(0),interval<T>(0),interval<T>(pow(q,nu+1)),interval<T>(q),interval<T>(-z*z/4))/Euler(interval<T>(q));
}
 return res;
}
template <class T> complex<interval<T> >Jackson1(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
complex<interval<T> >res;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 res=pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pow(q,nu+1)),interval<T>(q))*
  Heine(complex<interval<T> >(0),complex<interval<T> >(0),complex<interval<T> >(pow(q,nu+1)),interval<T>(q),complex<interval<T> >(-z*z/4))/Euler(interval<T>(q));
 return res;       
}

template <class T> complex<interval<T> >modified_qBessel1(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  // verification program for 1st modified q-Bessel
complex<interval<T> >res;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 res=pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pow(q,nu+1)),interval<T>(q))*
  Heine(complex<interval<T> >(0),complex<interval<T> >(0),complex<interval<T> >(pow(q,nu+1)),interval<T>(q),complex<interval<T> >(z*z/4))/Euler(interval<T>(q));
 return res;       
}


template <class T> interval<T>Jackson2(const interval<T>& z,const interval<T>& nu,const interval<T>& q){
// verification program for Jackson`s 2nd q-Bessel function
  interval<T>res,a,b,c,series,pq;
int j,K;
T rad;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 if (pow(q,nu)<1 && z*z<4*q){
j=1;
K=500;
a=1.;
b=1.;
c=1.;
// Leibniz criterion
 for (int k=1; k<=K; k++){
    j = -1*j;
    b=b*(1-pow(q,k));
    c=c*(1-pow(q,k+nu));
    a = a+j*pow(q,k*(k-1))*pow(pow(q,nu+1)*z*z/4,k)/(b*c);
  }
  b=b*(1-pow(q,K+1));
  c=c*(1-pow(q,K+nu+1));
 
 rad=abs(pow(q,K*(K+1))*pow(pow(q,nu+1)*z*z/4,K+1)/(b*c)).upper();
 series=a+rad*interval<T>(-1.,1.);
 res=pow(z/2,nu)*series/Karpelevich(interval<T>(pow(q,nu+1)),interval<T>(q));
}
 else{
   /* implementation by original definition
   res=pow(z*0.5,nu)*infinite_qPochhammer(interval<T>(pow(q,nu+1)),interval<T>(q))*
     _0phi_1(interval<T>(pow(q,nu+1)),interval<T>(q),interval<T>(-z*z*pow(q,nu+1)/4))/Euler(interval<T>(q));
   */
 // alternative implementaion
 // reference
 // H. T Koelink, Hansen-Lommel Orothogonality Relations for Jackson`s q-Bessel functions, formula 3.2
 // Journal of Mathematical Analysis and Applications 175, 425-437 (1993)
 pq=pow(q,nu+1);
 res=pow(z/2,nu)*_1phi_1(interval<T> (-z*z/4),interval<T> (0),interval<T>(q), interval<T> (pq))/Euler(interval<T>(q));

 }
  return res;
}
template <class T> complex<interval<T> >Jackson2(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  complex<interval<T> >res,pq;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 pq=pow(q,nu+1);
 /* implementation by original definition

 res=pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pq),interval<T>(q))*
   _0phi_1(complex<interval<T> >(pq),interval<T>(q),complex<interval<T> >(-z*z*pq/4))/Euler(interval<T>(q));
 */
 // alternative implementaion
 // reference
 // H. T Koelink, Hansen-Lommel Orothogonality Relations for Jackson`s q-Bessel functions, formula 3.2
 // Journal of Mathematical Analysis and Applications 175, 425-437 (1993)
 res=pow(z/2,nu)*_1phi_1(complex<interval<T> >(-z*z/4),complex<interval<T> >(0),interval<T>(q), complex<interval<T> >(pq))/Euler(interval<T>(q));
 return res;       
}

template <class T> interval<T>Jackson2_integral(const interval<T>& z,const int& nu,const interval<T>& q){
// verification program for Jackson`s 2nd q-Bessel function
// Integral representation is used

// References
// Rahman(1987), An Integral Representation and Some Transformation Properties of q-Bessel Functions, Journal of Mathematical Analysis and Applications 125
// Zhang(2008), Plancherel-Rotach asymptotics for certain basic hypergeometric series, Advances in Mathematics 217, Lemma 1.1
// Ismail, Zhang(2016), Integral and Series Representations of q-Polynomials and Functions: PartI, arXiv, Lemma 6.1

  if(nu<=0){
    throw std::domain_error("nu must be positive");
  }
  /* abolished implementation
     if(abs(q)>=1){
     throw std::domain_error("absolute value of q must be under 1");
     }
  
  interval<T> res,denom,num1,num2,integral;
  T denomrad,numrad1,numrad2;
 
  numrad1=(exp(1/(1-q))/(1-q)).upper();
  
  num1=(1.+numrad1*interval<T>(-1.,1.))*(1.+numrad1*interval<T>(-1.,1.));
  if(abs(z*sqrt(q)/2)*pow(q,nu/2)/(1-q)<0.5){
    numrad2=(abs(z)*pow(q,(nu+1)/2)/(1-q)).upper();
  }
  else{
    numrad2=(exp(abs(z*0.5*pow(q,(nu+1)/2))/(1-q))*abs(z*0.5*pow(q,(nu+1)/2))/(1-q)).upper();
  }
  num2=(1.+numrad2*interval<T>(-1.,1.))*(1.+numrad2*interval<T>(-1.,1.));
  //std::cout<<num2<<std::endl;
  if(pow(q,nu)/(1-q)<0.5){
    denomrad=(2*pow(q,nu)/(1-q)).upper();
  }
  else{
    throw std::domain_error("denominator cannot be evaluated");
  }
  denom=(1.+denomrad*interval<T>(-1.,1.))*(1.+denomrad*interval<T>(-1.,1.));
  integral=denom*num1*num2;
  res=integral*pow(z/2,nu)*infinite_qPochhammer(interval<T>(pow(q,2*nu)),interval<T>(q))/infinite_qPochhammer(interval<T>(pow(q,nu)),interval<T>(q))*0.5;
  return res;
*/
  if(abs(q)>=0.5){
    throw std::domain_error("implemented only for q under 0.5");
  }
  interval<T> res,denom,num1,num2,integral;
  T denomrad,numrad1,numrad2;
  if(q<1/3){
    numrad1=(2*q/(1-q)).upper();
  }
  else{
    numrad1=(2*q*q/(1-q)).upper();
  }
  num1=(1.+numrad1*interval<T>(-1.,1.))*(1.+numrad1*interval<T>(-1.,1.));

   if(abs(z*sqrt(q)/2)*pow(q,nu/2)/(1-q)<0.5){
    numrad2=(abs(z)*pow(q,(nu+1)/2)/(1-q)).upper();
  }
  else{
    numrad2=(exp(abs(z*0.5*pow(q,(nu+1)/2))/(1-q))*abs(z*0.5*pow(q,(nu+1)/2))/(1-q)).upper();
  }
  num2=(1.+numrad2*interval<T>(-1.,1.))*(1.+numrad2*interval<T>(-1.,1.));
  if(pow(q,nu)/(1-q)<0.5){
    denomrad=(2*pow(q,nu)/(1-q)).upper();
  }
  else{
    throw std::domain_error("denominator cannot be evaluated");
  }
  denom=(1.+denomrad*interval<T>(-1.,1.))*(1.+denomrad*interval<T>(-1.,1.));
  integral=denom*num1*num2;
  if(q<1/3){
    res=integral*pow(z/2,nu)*infinite_qPochhammer(interval<T>(pow(q,2*nu)),interval<T>(q))/infinite_qPochhammer(interval<T>(pow(q,nu)),interval<T>(q));
}
  else{
    res=integral*pow(z/2,nu)*infinite_qPochhammer(interval<T>(pow(q,2*nu)),interval<T>(q))/infinite_qPochhammer(interval<T>(pow(q,nu)),interval<T>(q))*(q*q+q+1);

  }
  return res;
}

template <class T> complex<interval<T> >modified_qBessel2(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  //verification program for 2nd modified q-Bessel function
  complex<interval<T> >res,pq;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 pq=pow(q,nu+1);
 res=pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pq),interval<T>(q))*
   _0phi_1(complex<interval<T> >(pq),interval<T>(q),complex<interval<T> >(z*z*pq/4))/Euler(interval<T>(q));
 return res;       
}
template <class T> interval<T> Hahn_Exton(const interval<T> & z,const interval<T> & nu,const interval<T>& q){
  // verification program for Hahn-Exton q-Bessel function
  interval<T> res,pq;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 pq=pow(q,nu+1);
 /* implementation by original definition
 res=pow(z,nu)*infinite_qPochhammer(interval<T> (pq),interval<T>(q))*
   _1phi_1(interval<T> (0),interval<T> (pq),interval<T>(q),interval<T> (z*z*q))/Euler(interval<T>(q));
 */
 res=pow(z,nu)*infinite_qPochhammer(interval<T> (z*z*q),interval<T>(q))*_1phi_1(interval<T> (0),interval<T> (z*z*q),interval<T>(q),interval<T> (pq))/Euler(interval<T> (q));
 return res;       
}
template <class T> complex<interval<T> >Hahn_Exton(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  // verification program for Hahn-Exton q-Bessel function
  complex<interval<T> >res,pq;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 pq=pow(q,nu+1);
 /* implementation by original definition
 res=pow(z,nu)*infinite_qPochhammer(complex<interval<T> >(pq),interval<T>(q))*
   _1phi_1(complex<interval<T> >(0),complex<interval<T> >(pq),interval<T>(q),complex<interval<T> >(z*z*q))/Euler(interval<T>(q));
 */
 // alternative implementation
 // A. B. Olde Daalhuis, Asymptotic Expansions for q-Gamma, q-Exponential and q-Bessel Functions, formula 4.6
 // Journal of Mathematical Analysis and Applications 186, 896-913 (1994)
 res=pow(z,nu)*infinite_qPochhammer(complex<interval<T> >(z*z*q),interval<T>(q))*_1phi_1(complex<interval<T> >(0),complex<interval<T> >(z*z*q),interval<T>(q), complex<interval<T> >(pq))/Euler(interval<T> (q));
 return res;       
}
template <class T> complex<interval<T> >little(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  // verification program for little q-Bessel function
  // reference
  // Koornwinder and Swarttouw, On q-analogues of the Fourier and Hankel transforms, 1992
  // Bouzeffour, New Addition Formula for the Little q-Bessel Functions, arXiv, 2013
  complex<interval<T> >res,pq;
if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 pq=pow(q,nu+1);
 res=pow(z,nu)*infinite_qPochhammer(complex<interval<T> >(pq),interval<T>(q))*
   _1phi_1(complex<interval<T> >(0),complex<interval<T> >(pq),interval<T>(q),complex<interval<T> >(z))/Euler(interval<T>(q));
 return res;       
}
template <class T> complex<interval<T> >qNeumann(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  // verification program for the q-Neumann function defined by Swarttouw
  // Neumann function = 2nd kind Bessel function
  // reference
  // A. B. Olde Daalhuis, Asymptotic Expansions for q-Gamma, q-Exponential and q-Bessel Functions, formula 4.6
  // Journal of Mathematical Analysis and Applications 186, 896-913 (1994)
  // R. F. Swarttouw, "The Hahn-Exton $q$-Bessel function", Ph.D. thesis, Delft University of Technology, 1992
  complex<interval<T> > res,qg;
  qg=q_gamma(complex<interval<T> >(0.5),interval<T> (q));
  res=q_gamma(complex<interval<T> >(nu),interval<T>(q))*q_gamma(complex<interval<T> >(1-nu),interval<T>(q))
    *(pow(q,nu*(nu+1)*0.5)*Hahn_Exton(complex<interval<T> >(z),complex<interval<T> >(nu),interval<T> (q))
      /q_gamma(complex<interval<T> >(0.5-nu),interval<T> (q))/q_gamma(complex<interval<T> >(0.5+nu),interval<T> (q))
      -Hahn_Exton(complex<interval<T> >(z*pow(q,-nu*0.5)),complex<interval<T> >(-nu),interval<T> (q))/pow(qg,2));
  return res;
}

}

#endif
