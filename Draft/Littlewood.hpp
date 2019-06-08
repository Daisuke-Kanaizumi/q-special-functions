// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp

// Verified computation of infinite q-Pochhammer symbol by using Littlewood's theorem
// Reference
// J. E. Littlewood, On the Asymptotic Approximation to Integral Functions of Zero Order,
// Proceedings of the London Mathematical Society, (2)5  (1907) 361-410.
// This is reprinted in Collected Papers, Vol2 (Oxford Univ. Press, 1970) 1059-1108.
// Corrected in
// Y. Chen , M. E. Ismail, K. A. Muttalib. Asymptotics of basic Bessel functions and q-Laguerre polynomials, Theorem 5
// Journal of Computational and Applied Mathematics, 54(3), 263-272 (1994).

#ifndef LITTLEWOOD_HPP
#define LITTLEWOOD_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/complex.hpp>
#include <kv/constants.hpp>
#include <kv/Heine.hpp>
#include <kv/Pochhammer.hpp>
#include <kv/Gatteschi.hpp>
#include <cmath>

namespace kv{
  template <class T> complex<interval<T> >qPL(const complex<interval<T> >& z,const interval<T>& q){
    int n,s;
    T rad;
    interval<T>r,theta,omega,beta;
    complex<interval<T> >a,res,i,series,sum;
    i=complex<interval<T> >::i();
    s=120;
    theta=arg(-z);
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("value of q must be positive");
    }
    if(abs(z)<=1){
      throw std::domain_error("value of z must be larger");
    }
    r=abs(z);
    omega=-log(q);
    //n=std::floor(mid(log(r)/omega)-1);
    n=1;
    while((n+1)*omega<log(r)){
      n=n+1;
    }
    //std::cout<<n<<std::endl; //OK 
    beta=n+1-log(r)/omega;
    a=beta-i*theta/omega;
    // std::cout<<a<<std::endl; //OK
    series=0.;
    for(int t=1;t<=s-1;t++){
      series=series+1./t/pow(z,t)/(1-pow(q,t));
    }
    rad=(abs(1./s/pow(z,s)/(1-pow(q,s)))/(1-abs(1/z/q))).upper();
    sum=complex_nbd(series,rad);
    res=(1-z)*exp(log(-z)*log(-z)*0.5/omega-log(-z)*0.5+0.5*omega*a-0.5*omega*a*a
		  +log(infinite_qPochhammer(complex<interval<T> >(-pow(q,1-a)),interval<T>(q)))
		  +log(infinite_qPochhammer(complex<interval<T> >(-pow(q,a)),interval<T>(q)))-sum);
    //std::cout<<beta<<std::endl; //OK
    //std::cout<<sum<<std::endl; //OK
    //std::cout<<infinite_qPochhammer(complex<interval<T> >(-pow(q,a)),interval<T>(q))<<std::endl; 
    return res;
   
  }
}
#endif
