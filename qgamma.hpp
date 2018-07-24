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
#include <kv/qPochhammerVer2.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ub = boost::numeric::ublas;
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
template <class T> ub::matrix<interval<T> >MExp(const ub::matrix<interval<T> >& A){
  int n,M;
  M=100;
  n=A.size1();//A:square matrix
  ub::matrix< interval<T> > B(n, n),res(n, n),sum(n, n),pro(n,n);
  interval<T> error,norm;
  T b;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      sum(i,j)=0.;
      if(i==j)pro(i,j)=1.;
      else pro(i,j)=0.;
    }
  }
  for(int N=0;N<=M;N++){
    sum+=(1./Pochhammer(interval<T>(1.,1.),N))*pro;
    pro=prod(pro,A);
  }
  norm=abs(A(0,0));
  for(int i1=0;i1<n;i1++){
    for(int j1=0;j1<n;j1++){
      if (A(i1,j1)>abs(norm)) norm=abs(A(i1,j1));
    }
  }

  error=exp(norm)*pow(norm,M+1)/Pochhammer(interval<T>(1.,1.),M+1);
  b=(abs(error)).upper();

  for(int k1=0;k1<n;k1++){
    for(int l1=0;l1<n;l1++){
      B(k1,l1).assign(-1.,1.);
      B(k1,l1)=b*B(k1,l1);
    }
  }
  res=sum+B;
  return res;
}
template <class T> ub::matrix<interval<T> >q_gamma(const ub::matrix<interval<T> >& A,const interval<T>& q){
  int n;
  n=A.size1();//A:square matrix
  interval<T>buf;
  ub::matrix< interval<T> > I(n, n),B(n, n),res(n, n),exp(n,n),inv(n,n),qp(n,n),qa(n,n),AA(n,n);
  
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(i==j){
	I(i,j)=1.;
	inv(i,j)=1.;
      }
      else{ 
	I(i,j)=0.;
	inv(i,j)=0.;
      }
    }      
  }
  AA=log(q)*A;
  qa=MExp(AA);
  qp=infinite_qPochhammer(ub::matrix<interval<T> > (qa), interval<T> (q));
//std::cout<<qp<<std::endl;
  for(int i1=0;i1<n;i1++){
    buf=1./qp(i1,i1);
    for(int j1=0;j1<n;j1++){
      qp(i1,j1)*=buf;
      inv(i1,j1)*=buf;
    }
 
    for(int j2=0;j2<n;j2++){
      if(i1!=j2){
	buf=qp(j2,i1);
	for(int k=0;k<n;k++){
	  qp(j2,k)-=qp(i1,k)*buf;
	  inv(j2,k)-=inv(i1,k)*buf;
	}
      }
    }
  }
  B=I-A;
  B=log(1-q)*B;
  std::cout<<B<<std::endl;
  exp=MExp(B);
  res=infinite_qPochhammer(interval<T>(q),interval<T>(q))
    *prod(inv,exp);
  return res;
}
  template <class T> complex<interval<T> >qgamma_Gauss_multi(const complex<interval<T> >& z,const interval<T>& q, int  p=3){
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
  template <class T> complex<interval<T> >qgamma_shift(const complex<interval<T> >& z,const interval<T>& q, int p=3){
    // computing the q-gamma function with functional equation
    // G Gasper , M Rahman, Basic Hypergeometric Series 2nd Edition, Cambridge University Press, 2004.
    complex<interval<T> >pro;
    pro=1.;    
    for(int i=1;i<=p;i++){
      pro=pro*(1-pow(q,z-i))/(1-q);
    }
    pro=pro*q_gamma(complex<interval<T> >(z-p),interval<T>(q));
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
 template <class T> complex<interval<T> >elliptic_gamma(const complex<interval<T> >& z,const interval<T> & p ,const interval<T> & q){
    // verification program for elliptic gamma function
    // reference: M. A. Bershtein, A. I. Shechechkin (arXiv, 2016)
    // q-deformed Painlev\`e \tau function and q-deformed conformal blocks, Appendix A
    complex<interval<T> >res;
    /* if (abs(z)>=1){
      throw std::domain_error("implemented only for |z|<1");
      }*/
    if (abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    if (abs(p)>=1){
      throw std::domain_error("absolute value of p must be under 1");
    }
    res=inf_elliptic_Pochhammer(complex<interval<T> >(p*q/z),complex<interval<T> >(p),complex<interval<T> >(q))
      /inf_elliptic_Pochhammer(complex<interval<T> >(z),complex<interval<T> >(p),complex<interval<T> >(q));
    return res;
  }
  template <class T> complex<interval<T> > modified_Jacobi_theta(const complex<interval<T> >&a,const interval<T> &q){
    complex<interval<T> >res;
    res=infinite_qPochhammer(complex<interval<T> >(a),interval<T>(q))
      *infinite_qPochhammer(complex<interval<T> >(q/a),interval<T>(q));
    return res;
  } 
  template <class T> interval<T>  modified_Jacobi_theta(const interval<T> &a,const interval<T> &q){
    interval<T> res;
    res=qPVer2(interval<T> (a),interval<T>(q))
      *qPVer2(interval<T> (q/a),interval<T>(q));
    return res;
  }
  template <class T> complex<interval<T> >elliptic_gamma_tilde(const complex<interval<T> >& z,const interval<T> & p ,const interval<T> & q){
    complex<interval<T> >res;
    if (abs(z)>=1){
      throw std::domain_error("implemented only for |z|<1");
    }
    if (abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    if (abs(p)>=1){
      throw std::domain_error("absolute value of p must be under 1");
    }
    res=infinite_qPochhammer(interval<T> (q),interval<T> (q))/infinite_qPochhammer(interval<T> (p),interval<T> (p))
      *pow(modified_Jacobi_theta(q,p),1-log(z)/log(q))*elliptic_gamma(complex<interval<T> >(z),interval<T>(p),interval<T>(q));
    return res;
  }
  template <class T> complex<interval<T> >elliptic_gamma_shift(const complex<interval<T> >& z,const interval<T> & p ,const interval<T> & q,int n){
    complex<interval<T> >res,pro1,pro2;
    interval<T> r;
    r=pow(q,n);
    if (abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    if (abs(p)>=1){
      throw std::domain_error("absolute value of p must be under 1");
    }
    for(int i=1;i<=n-1;i++){
      pro1=pro1*elliptic_gamma_tilde(complex<interval<T> >(i/T(n)),interval<T>(p),interval<T>(r));
    }
    for(int j=0;j<=n-1;j++){
      pro2=pro2*elliptic_gamma_tilde(complex<interval<T> >((z+j)/T(n)),interval<T>(p),interval<T>(r));
    }
    res=pow(modified_Jacobi_theta(r,p)/modified_Jacobi_theta(q,p),z-1)*pro2/pro1;
    return res;
  }
 template <class T> complex<interval<T> >elliptic_gamma_shift2(const complex<interval<T> >& z,const interval<T> & p ,const interval<T> & q,int n){
    complex<interval<T> >res,pro1,pro2;
    interval<T> r;
    r=pow(q,n);
    if (abs(q)>=1){
      throw std::domain_error("absolute value of q must be under 1");
    }
    if (abs(p)>=1){
      throw std::domain_error("absolute value of p must be under 1");
    }
    for(int i=1;i<=n-1;i++){
      pro1=pro1*elliptic_gamma_tilde(complex<interval<T> >(T(i)/T(n)),interval<T>(p),interval<T>(r));
    }
    for(int j=0;j<=n-1;j++){
      pro2=pro2*elliptic_gamma_shift(complex<interval<T> >((z+T(j))/T(n)),interval<T>(p),interval<T>(r),int(n));
    }
    res=pow(modified_Jacobi_theta(r,p)/modified_Jacobi_theta(q,p),z-1)*pro2/pro1;
    return res;
  }
}

#endif
