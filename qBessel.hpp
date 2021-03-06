// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University

// verification program for q-Bessel functions
// April 12th, 2018

#ifndef QBESSEL_HPP
#define QBESSEL_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/constants.hpp>
#include <kv/complex.hpp>
#include <kv/convert.hpp>
#include <kv/defint.hpp>
#include <kv/Heine.hpp>
#include <kv/Pochhammer.hpp>
#include <kv/QHypergeometric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <limits>
namespace ub = boost::numeric::ublas;
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
  complex<interval<T> >res,i;


if(abs(q)>=1){
throw std::domain_error("absolute value of q must be under 1");
}
 if(abs(z).upper()>=2){
throw std::domain_error("Jackson`s 1st q-Bessel function is not defined");
}

res=pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pow(q,nu+1)),interval<T>(q))*
  Heine(complex<interval<T> >(0),complex<interval<T> >(0),complex<interval<T> >(pow(q,nu+1)),interval<T>(q),complex<interval<T> >(-z*z/4))/Euler(interval<T>(q));
 
 return res;       
}

template <class T> complex<interval<T> >modified_qBesselI1(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  // verification program for 1st modified q-Bessel I1
  complex<interval<T> >res,i;
  if(abs(q)>=1){
    throw std::domain_error("absolute value of q must be under 1");
  }
  if(abs(z).upper()>=2){
    throw std::domain_error("1st modified q-Bessel function is not defined");
  }
  interval<T>pi;
  i=complex<interval<T> >::i();
  pi=constants<interval<T> >::pi();
  res=exp(-i*nu*pi/2)*Jackson1(complex<interval<T> >(i*z),complex<interval<T> >(nu),interval<T> (q));
  return res;       
}
  
  template <class T> complex<interval<T> >modified_qBesselK1(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
    // verification program for 1st modified q-Bessel K1, nu should not be an integer
    complex<interval<T> >res;
    if(abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    if(abs(z).upper()>=2){
      throw std::domain_error("1st modified q-Bessel function is not defined");
    }
    interval<T>pi;
    pi=constants<interval<T> >::pi();
    res=pi*0.5/sin(pi*nu)*
      (modified_qBesselI1(complex<interval<T> > (z), complex<interval<T> > (-nu), interval<T> (q))-modified_qBesselI1(complex<interval<T> > (z), complex<interval<T> > (nu), interval<T> (q)));
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

// if (pow(q,nu)<1 && z*z<4*q){
/*j=1;
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
*/ //}
// else{
   /* implementation by original definition*/
// res=pow(z*0.5,nu)*infinite_qPochhammer(interval<T>(pow(q,nu+1)),interval<T>(q))*
//  _0phi_1(interval<T>(pow(q,nu+1)),interval<T>(q),interval<T>(-z*z*pow(q,nu+1)/4))/Euler(interval<T>(q));
   
 // alternative implementaion
 // reference
 // H. T Koelink, Hansen-Lommel orthogonality Relations for Jackson`s q-Bessel functions, formula 3.2
 // Journal of Mathematical Analysis and Applications 175, 425-437 (1993)
 pq=pow(q,nu+1);
 res=pow(z/2,nu)*_1phi_1(interval<T> (-z*z/4),interval<T> (0),interval<T>(q), interval<T> (pq))/Euler(interval<T>(q));
 // }

  return res;
}
template <class T> complex<interval<T> >Jackson2(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
  complex<interval<T> >res,pq,i;
  i=complex<interval<T> >::i();
  if(abs(q)>=1){
    throw std::domain_error("absolute value of q must be under 1");
  }
  pq=pow(q,nu+1);
  /* implementation by original definition*/
     
  //  res=pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pq),interval<T>(q))*
  //  _0phi_1(complex<interval<T> >(pq),interval<T>(q),complex<interval<T> >(-z*z*pq/4))/Euler(interval<T>(q));
  
  // alternative implementaion
  // reference
  // H. T Koelink, Hansen-Lommel Orthogonality Relations for Jackson`s q-Bessel functions, formula 3.2
  // Journal of Mathematical Analysis and Applications 175, 425-437 (1993)
   try{
  res=pow(z/2,nu)*_1phi_1(complex<interval<T> >(-z*z/4),complex<interval<T> >(0),interval<T>(q), complex<interval<T> >(pq))/Euler(interval<T>(q));
  }
  // alternative implementaion
  // reference
  // Y. Chen , M. E. Ismail, K. A. Muttalib. Asymptotics of basic Bessel functions and q-Laguerre polynomials, Lemma 2
  // Journal of Computational and Applied Mathematics, 54(3), 263-272 (1994).
  catch(std::domain_error){
    ub::vector< complex<interval<T> > > a(3);
    ub::vector< complex<interval<T> > > b(2);
    ub::vector< complex<interval<T> > > c(2);
    a(0)=pow(q,(nu+0.5)*0.5);
    a(1)=-pow(q,(nu+0.5)*0.5);
    a(2)=0.;
    b(0)=-sqrt(q);
    b(1)=i*a(0)*z*0.5;
    c(0)=-sqrt(q);
    c(1)=-b(1);

    res=pow(z*0.5,nu)*infinite_qPochhammer(interval<T>(sqrt(q)),interval<T>(q))/Euler(interval<T>(q))*0.5
      *(infinite_qPochhammer(complex<interval<T> >(b(1)),interval<T>(sqrt(q)))
	*QHypergeom(ub::vector<complex<interval<T> > >(a),ub::vector<complex<interval<T> > >(b),interval<T>(sqrt(q)),complex<interval<T> >(sqrt(q),0))+
	infinite_qPochhammer(complex<interval<T> >(c(1)),interval<T>(sqrt(q)))
	*QHypergeom(ub::vector<complex<interval<T> > >(a),ub::vector<complex<interval<T> > >(c),interval<T>(sqrt(q)),complex<interval<T> >(sqrt(q),0)));
    	} 
  
    if((abs(res)).upper()==std::numeric_limits<T>::infinity()){  
    ub::vector< complex<interval<T> > > a(3);
    ub::vector< complex<interval<T> > > b(2);
    ub::vector< complex<interval<T> > > c(2);
    a(0)=pow(q,(nu+0.5)*0.5);
    a(1)=-pow(q,(nu+0.5)*0.5);
    a(2)=0.;
    b(0)=-sqrt(q);
    b(1)=i*a(0)*z*0.5;
    c(0)=-sqrt(q);
    c(1)=-b(1);

    res=pow(z*0.5,nu)*infinite_qPochhammer(interval<T>(sqrt(q)),interval<T>(q))/Euler(interval<T>(q))*0.5
      *(infinite_qPochhammer(complex<interval<T> >(b(1)),interval<T>(sqrt(q)))
	*QHypergeom(ub::vector<complex<interval<T> > >(a),ub::vector<complex<interval<T> > >(b),interval<T>(sqrt(q)),complex<interval<T> >(sqrt(q),0))+
	infinite_qPochhammer(complex<interval<T> >(c(1)),interval<T>(sqrt(q)))
	*QHypergeom(ub::vector<complex<interval<T> > >(a),ub::vector<complex<interval<T> > >(c),interval<T>(sqrt(q)),complex<interval<T> >(sqrt(q),0)));
     }
    //std::cout<<qPochhammer(complex<interval<T> >(pow(q,(nu+0.5)*0.5)),interval<T>(q),int(1000))
    //*qPochhammer(complex<interval<T> >(-pow(q,(nu+0.5)*0.5)),interval<T>(q),int(1000))<<std::endl;

    return res;       
}
  template <class TT> struct qBesselintegral_nu_int_real {
    TT x, q;
    int nu,n; // Setting parameters
    qBesselintegral_nu_int_real(TT x, TT q, int nu,int n) : x(x),q(q),nu(nu),n(n) {}
    
    template <class T> T operator() (const T& t) {
      complex<T>  pro;
      T proreal;
      complex<T> i;
      pro=1.;
      i=complex<T>::i();
      for(int k=0;k<=nu-1;k++){
	pro=pro*(1-pow(T(q),k)*exp(2*i*t))*(1-pow(T(q),k)*exp(-2*i*t));
      }
      for(int j=0;j<=n-1;j++){
	pro=pro*(1+i*T(x)*pow(T(q),T(nu)/2+0.5+j)*exp(i*t)/2)*(1+i*T(x)*pow(T(q),T(nu)/2+0.5+j)*exp(-i*t)/2);
      }
      
      proreal=pro.real();
      return proreal;	  
    }
  };
  template <class TT> struct qBesselintegral_nu_int_imag {
    TT x, q;
    int nu,n; // Setting parameters
    qBesselintegral_nu_int_imag(TT x, TT q, int nu,int n) : x(x),q(q),nu(nu),n(n) {}
    
    template <class T> T operator() (const T& t) {
      complex<T>  pro;
      T proimag;
      complex<T> i;
      pro=1.;
      i=complex<T>::i();
      for(int k=0;k<=nu-1;k++){
	pro=pro*(1-pow(T(q),k)*exp(2*i*t))*(1-pow(T(q),k)*exp(-2*i*t));
      }
      for(int j=0;j<=n-1;j++){
	pro=pro*(1+i*T(x)*pow(T(q),T(nu)/2+0.5+j)*exp(i*t)/2)*(1+i*T(x)*pow(T(q),T(nu)/2+0.5+j)*exp(-i*t)/2);
      }
      
      proimag=pro.imag();
      return proimag;	  
    }
  };
  template <class T> complex<interval<T> >Jackson2_integral(const interval<T>& z,const int & nu,const interval<T>& q){
      // verification program for Jackson`s 2nd q-Bessel function
      // Integral representation is used
    
      // References
      // Rahman(1987), An Integral Representation and Some Transformation Properties of q-Bessel Functions, Journal of Mathematical Analysis and Applications 125
      // Zhang(2008), Plancherel-Rotach asymptotics for certain basic hypergeometric series, Advances in Mathematics 217, Lemma 1.1
    int n;    
    n=100;
    if(nu<=0){
      throw std::domain_error("nu must be positive");
    }
    if(abs(z)*pow(q,nu/2+n+0.5)/2/(1-q)>=0.5){
      n=n+10;
    }
    complex<interval<T> >res,integral;
    interval<T>num,realint,imagint,pi;
    T numrad;
    numrad=(abs(z)*pow(q,nu/2+n+0.5)*2/(1-q)).upper();
    num=pow((1+numrad*interval<T>(-1.,1.)),2);
    pi=constants<interval<T> >::pi();
    realint=defint(qBesselintegral_nu_int_real<interval<T> >(z,q,nu,n),interval<T>(0.),interval<T>(pi),10,10);
    imagint=defint(qBesselintegral_nu_int_imag<interval<T> >(z,q,nu,n),interval<T>(0.),interval<T>(pi),10,10);
    integral=complex<interval<T> >(realint,imagint);
    
    res=num*integral*pow(z/2,nu)*infinite_qPochhammer(interval<T>(pow(q,2*nu)),interval<T>(q))/infinite_qPochhammer(interval<T>(pow(q,nu)),interval<T>(q))/(2*pi);
    
    
    return res;
  }
    template <class TT> struct qBesselintegral_nu_double_real {
      TT x,q,nu;
      int n; // Setting parameters
      qBesselintegral_nu_double_real(TT x, TT q, TT nu,int n) : x(x),q(q),nu(nu),n(n) {}
      
      template <class T> T operator() (const T& t) {
      complex<T>  pro;
      T proreal;
      complex<T> i;
      pro=1.;
      i=complex<T>::i();
      for(int k=0;k<=n-1;k++){
	pro=pro*(1-pow(T(q),k)*exp(2*i*t))*(1-pow(T(q),k)*exp(-2*i*t));
      }
      for(int j=0;j<=n-1;j++){
	pro=pro*(1+i*T(x)*pow(T(q),T(nu)/2+0.5+j)*exp(i*t)/2)*(1+i*T(x)*pow(T(q),T(nu)/2+0.5+j)*exp(-i*t)/2);
      }
      for(int l=0;l<=n-1;l++){
	pro=pro/(1-pow(T(q),l+nu)*exp(2*i*t))/(1-pow(T(q),l+nu)*exp(-2*i*t));
      }      
      proreal=pro.real();
      return proreal;	  
      }
  };
  template <class TT> struct qBesselintegral_nu_double_imag {
    TT x,q,nu;
    int n; // Setting parameters
    qBesselintegral_nu_double_imag(TT x, TT q, TT nu,int n) : x(x),q(q),nu(nu),n(n) {}
    
    template <class T> T operator() (const T& t) {
      complex<T>  pro;
      T proimag;
      complex<T> i;
      pro=1.;
      i=complex<T>::i();
      for(int k=0;k<=n-1;k++){
	pro=pro*(1-pow(T(q),k)*exp(2*i*t))*(1-pow(T(q),k)*exp(-2*i*t));
      }
      for(int j=0;j<=n-1;j++){
	pro=pro*(1+i*T(x)*pow(T(q),T(nu)/2+0.5+j)*exp(i*t)/2)*(1+i*T(x)*pow(T(q),T(nu)/2+0.5+j)*exp(-i*t)/2);
      }
      for(int l=0;l<=n-1;l++){
	pro=pro/(1-pow(T(q),l+nu)*exp(2*i*t))/(1-pow(T(q),l+nu)*exp(-2*i*t));
      }      
      
      proimag=pro.imag();
      return proimag;	  
    }
  };
  template <class T> complex<interval<T> >Jackson2_integral(const interval<T>& z,const interval<T> & nu,const interval<T>& q){
      // verification program for Jackson`s 2nd q-Bessel function
      // Integral representation is used
    
      // References
      // Rahman(1987), An Integral Representation and Some Transformation Properties of q-Bessel Functions, Journal of Mathematical Analysis and Applications 125
      // Zhang(2008), Plancherel-Rotach asymptotics for certain basic hypergeometric series, Advances in Mathematics 217, Lemma 1.1
    int n;    
    n=100;
    if(nu<=0){
      throw std::domain_error("nu must be positive");
    }
    if(abs(z)*pow(q,nu/2+n+0.5)/2/(1-q)>=0.5){
      n=n+10;
    }
    if(pow(q,nu+n)/2/(1-q)>=0.5){
      n=n+10;
    }

    complex<interval<T> >res,integral;
    interval<T>num1,num2,denom,realint,imagint,pi;
    T numrad1,numrad2,denomrad;
    numrad1=(abs(z)*pow(q,nu/2+n+0.5)*2/(1-q)).upper();
    num1=pow((1+numrad1*interval<T>(-1.,1.)),2);
    numrad2=(pow(q,n)*2/(1-q)).upper();
    num2=pow((1+numrad2*interval<T>(-1.,1.)),2);
    denomrad=(pow(q,nu+n)*2/(1-q)).upper();
    denom=pow((1+denomrad*interval<T>(-1.,1.)),2);
    pi=constants<interval<T> >::pi();
    realint=defint(qBesselintegral_nu_double_real<interval<T> >(z,q,nu,n),interval<T>(0.),interval<T>(pi),10,10);
    imagint=defint(qBesselintegral_nu_double_imag<interval<T> >(z,q,nu,n),interval<T>(0.),interval<T>(pi),10,10);
    integral=complex<interval<T> >(realint,imagint);
    
    res=num1*num2*denom*integral*pow(z/2,nu)*infinite_qPochhammer(interval<T>(pow(q,2*nu)),interval<T>(q))/infinite_qPochhammer(interval<T>(pow(q,nu)),interval<T>(q))/(2*pi);
    
    
    return res;
  }
    template <class TT> struct qBesselintegral_z_complex_real {
      TT q,nu;
      complex<TT> x;
      int n; // Setting parameters
      qBesselintegral_z_complex_real(complex<TT> x, TT q, TT nu,int n) : x(x),q(q),nu(nu),n(n) {}
      
      template <class T> T operator() (const T& t) {
      complex<T>  pro;
      T proreal;
      complex<T> i;
      pro=1.;
      i=complex<T>::i();
      for(int k=0;k<=n-1;k++){
	pro=pro*(1-pow(T(q),k)*exp(2*i*t))*(1-pow(T(q),k)*exp(-2*i*t));
      }
      for(int j=0;j<=n-1;j++){
	pro=pro*(1+i*complex<T>(x)*pow(T(q),T(nu)/2+0.5+j)*exp(i*t)/2)*(1+i*complex<T>(x)*pow(T(q),T(nu)/2+0.5+j)*exp(-i*t)/2);
      }
      for(int l=0;l<=n-1;l++){
	pro=pro/(1-pow(T(q),l+T(nu))*exp(2*i*t))/(1-pow(T(q),l+T(nu))*exp(-2*i*t));
      }      
      proreal=pro.real();
      return proreal;	  
      }
  };
  template <class TT> struct qBesselintegral_z_complex_imag {
    TT q,nu;
    complex<TT> x;
    int n; // Setting parameters
    qBesselintegral_z_complex_imag(complex<TT> x, TT q, TT nu,int n) : x(x),q(q),nu(nu),n(n) {}
    
    template <class T> T operator() (const T& t) {
      complex<T>  pro;
      T proimag;
      complex<T> i;
      pro=1.;
      i=complex<T>::i();
      for(int k=0;k<=n-1;k++){
	pro=pro*(1-pow(T(q),k)*exp(2*i*t))*(1-pow(T(q),k)*exp(-2*i*t));
      }
      for(int j=0;j<=n-1;j++){
	pro=pro*(1+i*complex<T>(x)*pow(T(q),T(nu)/2+0.5+j)*exp(i*t)/2)*(1+i*complex<T>(x)*pow(T(q),T(nu)/2+0.5+j)*exp(-i*t)/2);
      }
      for(int l=0;l<=n-1;l++){
	pro=pro/(1-pow(T(q),l+T(nu))*exp(2*i*t))/(1-pow(T(q),l+T(nu))*exp(-2*i*t));
      }      
      
      proimag=pro.imag();
      return proimag;	  
    }
  };
  template <class T> complex<interval<T> >Jackson2_integral(const complex<interval<T> >& z,const interval<T> & nu,const interval<T>& q){
      // verification program for Jackson`s 2nd q-Bessel function
      // Integral representation is used
    
      // References
      // Rahman(1987), An Integral Representation and Some Transformation Properties of q-Bessel Functions, Journal of Mathematical Analysis and Applications 125
      // Zhang(2008), Plancherel-Rotach asymptotics for certain basic hypergeometric series, Advances in Mathematics 217, Lemma 1.1
    int n;    
    n=100;
    if(nu<=0){
      throw std::domain_error("nu must be positive");
    }
    if(abs(z*pow(q,nu/2+n+0.5))/2/(1-q)>=0.5){
      n=n+10;
    }
    if(abs(pow(q,nu+n))/2/(1-q)>=0.5){
      n=n+10;
    }

    complex<interval<T> >res,integral,num1,denom;

    interval<T> num2;
    interval<T> pi,realint,imagint;
    T numrad1,numrad2,denomrad;
    numrad1=(abs(z*pow(q,nu/2+n+0.5))*2/(1-q)).upper();
    num1=pow(complex_nbd(complex<interval<T> >(1,0),numrad1),2);
    numrad2=(pow(q,n)*2/(1-q)).upper();
    num2=pow((1+numrad2*interval<T>(-1.,1.)),2);
    denomrad=(abs(pow(q,nu+n))*2/(1-q)).upper();
    denom=pow(complex_nbd(complex<interval<T> >(1,0),denomrad),2);
    pi=constants<interval<T> >::pi();
    realint=defint(qBesselintegral_z_complex_real<interval<T> >(z,q,nu,n),interval<T>(0.),interval<T>(pi),10,10);
    imagint=defint(qBesselintegral_z_complex_imag<interval<T> >(z,q,nu,n),interval<T>(0.),interval<T>(pi),10,10);
    integral=complex<interval<T> >(realint,imagint);
    
    res=num1*num2*denom*integral*pow(z/2,nu)*infinite_qPochhammer(complex<interval<T> >(pow(q,2*nu)),interval<T>(q))/infinite_qPochhammer(complex<interval<T> >(pow(q,nu)),interval<T>(q))/(2*pi);
    
    return res;
  }
  template <class T> interval<T> J2ratio(const interval<T> & z,const interval<T> & nu,const interval<T>& q){
    interval<T>res;
    res=Jackson2(interval<T>(z),interval<T>(nu),interval<T>(q))/Jackson2(interval<T>(z),interval<T>(nu-1),interval<T>(q));
    return res;
  }
  template <class T> complex<interval<T> >modified_qBesselI2(const complex<interval<T> >& z,const complex<interval<T> >& nu,const interval<T>& q){
    //verification program for 2nd modified q-Bessel function I2
    complex<interval<T> >res,i;
    if(abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    /*implementation by original definition
      interval<T>pi;
      i=complex<interval<T> >::i();
      pi=constants<interval<T> >::pi();
      res=exp(-i*nu*pi/2)*Jackson2(complex<interval<T> >(i*z),complex<interval<T> >(nu),interval<T> (q));
    */
    // alternative implementation
    // reference
    // Ismail, M. E., & Zhang, R. (2015). $ q $-Bessel Functions and Rogers-Ramanujan Type Identities.
    // arXiv preprint arXiv:1508.06861.
    res=pow(z/2,nu)/Euler(interval<T>(q))
      *_1phi_1(complex<interval<T> >(z*z/4),complex<interval<T> >(0.),interval<T>(q),complex<interval<T> >(pow(q,nu+1)));
    return res;       
  }
 template <class T> interval<T> modified_qBesselI2(const interval<T> & z,const interval<T> & nu,const interval<T>& q){
    //verification program for 2nd modified q-Bessel function I2
    interval<T> res;
    if(abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    /*implementation by original definition
      interval<T>pi;
      i=complex<interval<T> >::i();
      pi=constants<interval<T> >::pi();
      res=exp(-i*nu*pi/2)*Jackson2(complex<interval<T> >(i*z),complex<interval<T> >(nu),interval<T> (q));
    */
    // alternative implementation
    // reference
    // Ismail, M. E., & Zhang, R. (2015). $ q $-Bessel Functions and Rogers-Ramanujan Type Identities.
    // arXiv preprint arXiv:1508.06861.
    res=pow(z/2,nu)/Euler(interval<T>(q))
      *_1phi_1(interval<T> (z*z/4),interval<T> (0.),interval<T>(q),interval<T> (pow(q,nu+1)));
    return res;       
  }
  /*template <class T> interval<T> modified_qBesselI2_ae(const interval<T> & z,const interval<T> & nu,const interval<T>& q){
    //verification program for 2nd modified q-Bessel function I2, z>0
    // reference
    // Ismail, M. E., & Zhang, R. (2015). $ q $-Bessel Functions and Rogers-Ramanujan Type Identities.
    // arXiv preprint arXiv:1508.06861.
    interval<T> res;
    if(abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    res=pow(z/2.,nu)*infinite_qPochhammer(interval<T>(sqrt(q)),interval<T>(q))*0.5/Euler(interval<T>(q))
      *(infinite_qPochhammer(interval<T>(z*0.5*pow(q,(nu+0.5)*0.5)),interval<T>(sqrt(q)))+infinite_qPochhammer(interval<T>(-z*0.5*pow(q,(nu+0.5)*0.5)),interval<T>(sqrt(q))));
    return res;
    }*/
  template <class T> interval<T> Hahn_Exton(const interval<T> & z,const interval<T> & nu,const interval<T>& q){
    // verification program for Hahn-Exton q-Bessel function
    interval<T> res,pq;
    if(abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    pq=pow(q,nu+1);
    /* implementation by original definition*/
    /*    res=pow(z,nu)*infinite_qPochhammer(interval<T> (pq),interval<T>(q))*
	  _1phi_1(interval<T> (0),interval<T> (pq),interval<T>(q),interval<T> (z*z*q))/Euler(interval<T>(q));*/
       // alternative implementation
       // A. B. Olde Daalhuis, Asymptotic Expansions for q-Gamma, q-Exponential and q-Bessel Functions, formula 4.6
       // Journal of Mathematical Analysis and Applications 186, 896-913 (1994)
       
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
    /* implementation by original definition*/
    /*  res=pow(z,nu)*infinite_qPochhammer(complex<interval<T> >(pq),interval<T>(q))*
       _1phi_1(complex<interval<T> >(0),complex<interval<T> >(pq),interval<T>(q),complex<interval<T> >(z*z*q))/Euler(interval<T>(q));
    */
    // alternative implementation
    // A. B. Olde Daalhuis, Asymptotic Expansions for q-Gamma, q-Exponential and q-Bessel Functions, formula 4.6
    // Journal of Mathematical Analysis and Applications 186, 896-913 (1994)
    res=pow(z,nu)*infinite_qPochhammer(complex<interval<T> >(z*z*q),interval<T>(q))*_1phi_1(complex<interval<T> >(0),complex<interval<T> >(z*z*q),interval<T>(q), complex<interval<T> >(pq))/Euler(interval<T> (q));
     return res;       
  }
   template <class T> interval<T> HEratio(const interval<T> & z,const interval<T> & nu,const interval<T>& q){
    interval<T>res;
    res=Hahn_Exton(interval<T>(z),interval<T>(nu),interval<T>(q))/Hahn_Exton(interval<T>(z),interval<T>(nu-1),interval<T>(q));
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


}
  
#endif
