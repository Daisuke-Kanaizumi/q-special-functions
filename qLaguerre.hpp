// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University
// Email: daisuke15@asagi.waseda.jp


// verification program for q-Laguerre polynomial
// References

// Ismail, M. E., & Zhang, R. (2016).
// Integral and Series Representations of $ q $-Polynomials and Functions: Part I.
// arXiv preprint arXiv:1604.08441.
// Y. Chen , M. E. Ismail, & K. A. Muttalib. Asymptotics of basic Bessel functions and q-Laguerre polynomials, Lemma 2
// Journal of Computational and Applied Mathematics, 54(3), 263-272 (1994).

//Koekoek, R., & Swarttouw, R. F. (1996). The Askey-scheme of hypergeometric orthogonal polynomials and its q-analogue. arXiv preprint math/9602214.

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
#include <kv/psa.hpp>
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
template <class T> interval<T>  qLpoly(const interval<T> & x,const interval<T>& q,const int n,const interval<T> & alpha){
interval<T> res,sum;
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
ub::vector<interval<T>  >a(1),b(2);
a(0)=pow(q,alpha+n+1);
b(0)=pow(q,alpha+1);
b(1)=-x*pow(q,n+alpha+1);

res=qPochhammer(pow(q,alpha+1),q,n)/qPochhammer(q,q,n)
  *infinite_qPochhammer(interval<T> (-x*pow(q,alpha+n+1)),interval<T>(q))
  *QHypergeom(ub::vector<interval<T>  >(a),ub::vector<interval<T>  >(b),interval<T>(q),interval<T> (-x*pow(q,alpha+1)));

//std::cout<<QHypergeom(ub::vector<complex<interval<T> > >(a),ub::vector<complex<interval<T> > >(b),interval<T>(q),complex<interval<T> >(-x*pow(q,alpha+1)))<<std::endl;
}
if((abs(res)).upper()==std::numeric_limits<T>::infinity()){
res=infinite_qPochhammer(interval<T> (-x*pow(q,n+alpha+1)),interval<T> (q))
  /infinite_qPochhammer(interval<T> (pow(q,n+alpha+1)),interval<T> (q))/qPochhammer(q,q,n)
  *_1phi_1(interval<T> (-x),interval<T> (-x*pow(q,n+alpha+1)),interval<T>(q),interval<T> (pow(q,alpha+1)));
}
return res;
}
template <class T> interval<T>  qLpoly_psa(const interval<T> & x,const interval<T>& q,const int n,const interval<T> & alpha){
interval<T> res;

kv::psa<interval<T> >a,b,c;
a.v.resize(n+1);
b.v.resize(n+1);
c.v.resize(n+1);
for(int i=0;i<=n;i++){
a.v(i)=qPochhammer(-x,q,i)*std::pow(-1,i)*pow(q,i*(i+1)*0.5)*pow(q,alpha*i)/qPochhammer(q,q,i);
}
for(int j=0;j<=n;j++){
b.v(j)=std::pow(-1,j)*pow(q,j*(j-1)*0.5)/qPochhammer(q,q,j);
}
c=a/b;
res=c.v(n);
return res;
}
}
#endif
