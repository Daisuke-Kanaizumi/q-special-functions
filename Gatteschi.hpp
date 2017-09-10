// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp

// Verified computation of infinite q-Pochhammer symbol by the Gatteschi algorithm
// Reference: Gabutti, B. and Allasia, G., 2008. Evaluation of q-gamma function and q-analogues by iterative algorithms. Numerical Algorithms, 49(1), pp.159-168.
#ifndef GATTESCHI_HPP
#define GATTESCHI_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/complex.hpp>
#include <kv/Heine.hpp>
#include <cmath>

namespace kv{
  template <class T> complex<interval<T> >Gatteschi_qp(const complex<interval<T> >& z,const interval<T>& q, int N=100){
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("value of q must be positive");
    }
    
    int N0,m;
    N0=std::floor((0.5-log(abs(z))/log(q)).upper());
    while (abs(z*pow(q,N))>=1){
      N=N+50;
    }
    // This condition is not mentioned in the reference above, but necessary for the convergence of the error term
    if(N<N0){
      N=N0;
    }
    // Optimal value of N mentioned in the reference above
    m=std::floor((sqrt(2*log(std::pow(10.,-18.))/log(q))).upper());
    complex<interval<T> > mid,res,qn,sum,d,x;
    interval<T> qpower;
    qn=1.;
    sum=1.;
    d=1.;
    x=1.;
    for(int i=0;i<=N-1;i++){
      qn=qn*(1-z*pow(q,i));
    }
    // The reference above has mistaken Q_N(x) and G(x) (formula 16) 
    for(int k=1;k<=m-1;k++){
      d=d*pow(q,k-1)/(1-pow(q,k))*(1-q);
      x=x*(-z);
      qpower=qpower*pow(q,N);
      sum=sum+d*x*qpower;
    }
    mid=qn*sum;
    T rad;
    rad=(abs(qn)*pow(abs(z*pow(q,N)),m)*pow(q,m*(m-1)*0.5)/(1-abs(z*pow(q,N)))).upper();
    res=complex_nbd(mid,rad);
    return res;
  }
  template <class T> interval<T> Gatteschi_qp(const interval<T> & z,const interval<T>& q,int N=100){
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("value of q must be positive");
    }
    
    int N0,m;
    N0=std::floor((0.5-log(abs(z))/log(q)).upper());
    while (abs(z*pow(q,N))>=1){
      N=N+50;
    }
    // This condition is not mentioned in the reference above, but necessary for the convergence of the error term
    if(N<N0){
      N=N0;
    }
    // Optimal value of N mentioned in the reference above
    m=std::floor((sqrt(2*log(std::pow(10.,-18.))/log(q))).upper());
    interval<T>  mid,res,qn,sum,d,x,qpower;
    qn=1.;
    sum=1.;
    d=1.;
    x=1.;
    for(int i=0;i<=N-1;i++){
      qn=qn*(1-z*pow(q,i));
    }
    // The reference above has mistaken Q_N(x) and G(x) (formula 16) 
    for(int k=1;k<=m-1;k++){
      d=d*pow(q,k-1)/(1-pow(q,k))*(1-q);
      x=x*(-z);
      qpower=qpower*pow(q,N);
      sum=sum+d*x*qpower;
    }
    mid=qn*sum;
    T rad;
    rad=(abs(qn)*pow(abs(z*pow(q,N)),m)*pow(q,m*(m-1)*0.5)/(1-abs(z*pow(q,N)))).upper();
    res=mid+rad*interval<T>(-1.,1.);
    return res;
  }
}
#endif
