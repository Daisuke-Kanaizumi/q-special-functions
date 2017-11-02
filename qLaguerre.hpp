// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp


// verification program for q-Laguerre polynomial
// References

// Ismail, M. E., & Zhang, R. (2016).
// Integral and Series Representations of $ q $-Polynomials and Functions: Part I.
// arXiv preprint arXiv:1604.08441.
// Y. Chen , M. E. Ismail, K. A. Muttalib. Asymptotics of basic Bessel functions and q-Laguerre polynomials, Lemma 2
// Journal of Computational and Applied Mathematics, 54(3), 263-272 (1994).
#ifndef QLAGUERRE_HPP
#define QLAGUERRE_HPP
 
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/constants.hpp>
#include <kv/Pochhammer.hpp>
#include <kv/QHypergeometric.hpp>
#include <cmath>
#include <limits>
#include <kv/convert.hpp> // this was included to use complex numbers
#include <kv/complex.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ub = boost::numeric::ublas;
namespace kv {
template <class T> complex<interval<T> > qLpoly(const complex<interval<T> >& x,const interval<T>& q,const int n,const complex<interval<T> >& alpha){
complex<interval<T> >res,sum;
if(n<=1000){
sum=0.;
for(int k=0;k<=n;k++){
sum=sum+pow(q,alpha*k+k*k)*pow(-x,k)
  /qPochhammer(q,q,k)/qPochhammer(q,q,n-k)/qPochhammer(pow(q,alpha+1),q,k);
}
sum=sum*qPochhammer(pow(q,alpha+1),q,n);
res=sum;
}
 else{
ub::vector<complex<interval<T> > >a(1),b(2);
a(0)=pow(q,alpha+n+1);
b(0)=pow(q,alpha+1);
b(1)=-x*pow(q,n+alpha+1);

res=qPochhammer(pow(q,alpha+1),q,n)/qPochhammer(q,q,n)
  *infinite_qPochhammer(complex<interval<T> >(-x*pow(q,alpha+n+1)),interval<T>(q))
  *QHypergeom(ub::vector<complex<interval<T> > >(a),ub::vector<complex<interval<T> > >(b),interval<T>(q),complex<interval<T> >(-x*pow(q,alpha+1)));

//std::cout<<QHypergeom(ub::vector<complex<interval<T> > >(a),ub::vector<complex<interval<T> > >(b),interval<T>(q),complex<interval<T> >(-x*pow(q,alpha+1)))<<std::endl;
}
if((abs(res)).upper()==std::numeric_limits<T>::infinity()){
res=infinite_qPochhammer(complex<interval<T> >(-x*pow(q,n+alpha+1)),interval<T> (q))
  /infinite_qPochhammer(complex<interval<T> >(pow(q,n+alpha+1)),interval<T> (q))/qPochhammer(q,q,n)
  *_1phi_1(complex<interval<T> >(-x),complex<interval<T> >(-x*pow(q,n+alpha+1)),interval<T>(q),complex<interval<T> >(pow(q,alpha+1)));
  }
return res;
}
}
#endif
