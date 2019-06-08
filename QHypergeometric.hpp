// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University

// verification program for the q-Hypergeometric function

#ifndef QHYPERGEOMETRIC_HPP
#define QHYPERGEOMETRIC_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/complex.hpp>
#include <kv/qAiry.hpp>
#include <kv/Heine.hpp>
#include <kv/Pochhammer.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ub = boost::numeric::ublas;
namespace kv{
template <class T> interval<T> qPochhammer(const ub::vector<interval<T> >& a,const interval<T>& q,int n){
interval<T> res,pro;
int r;
r=a.size();
pro=1.;
for(int i=0;i<=r-1;i++){
pro=pro*qPochhammer(interval<T>(a(i)),interval<T>(q),int (n));
}
res=pro;
return res;
}
template <class T> complex<interval<T> >qPochhammer(const ub::vector<complex<interval<T> > >& a,const interval<T>& q,int n){
complex<interval<T> >res,pro;
int r;
r=a.size();
pro=1.;
for(int i=0;i<=r-1;i++){
  pro=pro*qPochhammer(complex<interval<T> >(a(i)),interval<T>(q),int (n));
}
res=pro;
return res;
}

template <class T> interval<T> QHypergeom(const ub::vector<interval<T> >& a,const ub::vector<interval<T> >& b,const interval<T>& q,const interval<T>& z){
  interval<T>res,mid,first,ratio,pro1,pro2;
  T rad;
  int r,s;
  r=a.size();
  s=b.size();
  mid=1.;
  pro1=1.;
  pro2=1.;
  int N;
  N=1000;
  for(int i=0;i<=s-1;i++){
    while(abs(b(i))>pow(1/q,N)){
      N=N+500;
    }
  }
  if (q>=1){
    throw std::domain_error("value of q must be under 1");
  }
  if (q<=0){
    throw std::domain_error("q must be positive");
  }
  if (r>s+1){
    throw std::domain_error("r>s+1 is not implemented");
  }
  if(r==s+1){
    if (abs(z)>=1){
      throw std::domain_error("absolute value of z must be under 1");
    }
    if(r==2&&s==1){
      // Heine hypergeometric function
      res=Heine(interval<T>(a(0)),interval<T>(a(1)),interval<T>(b(0)),interval<T>(q),interval<T>(z));
    }
    else{
      for(int n=1;n<=N-1;n++){
	mid=mid+qPochhammer(ub::vector<interval<T> >(a),interval<T>(q),int (n))*pow(z,n)
	  /qPochhammer(ub::vector<interval<T> >(b),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
      }
      first=abs(qPochhammer(ub::vector<interval<T> >(a),interval<T>(q),int (N))*pow(z,N)
		/qPochhammer(ub::vector<interval<T> >(b),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
      for(int j=0;j<=s-1;j++){
	pro1=pro1*(1+abs(b(j)-a(j))*pow(q,N)/abs(1-b(j)*pow(q,N)));
      }
      ratio=abs(z)*pro1*(1+pow(q,N)*abs(q-a(r-1))/abs(1-pow(q,N+1)));
      if(ratio<1){
	rad=(first/(1-ratio)).upper();
	res=mid+rad*interval<T>(-1.,1.);
      }
      else{
	throw std::domain_error("ratio is more than 1");
      }
      
    }
  }
  if(r<=s){
    if(r==0&&s==1){
      res=_0phi_1(interval<T>(b(0)),interval<T>(q),interval<T>(z));
    }
    if(r==1&&s==1){
      res=_1phi_1(interval<T>(a(0)),interval<T>(b(0)),interval<T>(q),interval<T>(z));
    }
    
    else{
      for(int n=1;n<=N-1;n++){
	mid=mid+qPochhammer(ub::vector<interval<T> >(a),interval<T>(q),int (n))*std::pow(-1,1+s-r)*pow(q,n*(n-1)*0.5*(1+s-r))*pow(z,n)
	  /qPochhammer(ub::vector<interval<T> >(b),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
      }
      first=abs(qPochhammer(ub::vector<interval<T> >(a),interval<T>(q),int (N))*std::pow(-1,1+s-r)*pow(q,N*(N-1)*0.5*(1+s-r))*pow(z,N)
		/qPochhammer(ub::vector<interval<T> >(b),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
      for(int k=0;k<=r-1;k++){
	pro1=pro1*(1+abs(b(k)-a(k))*pow(q,N)/abs(1-b(k)*pow(q,N)));
      }
      for(int h=r;h<=s-1;h++){
	pro2=pro2*pow(q,N*(1+s-r))/abs(1-b(h)*pow(q,N));
      }
      ratio=abs(z)*pro1*pro2/abs(1-pow(q,N));
      if(ratio<1){
	rad=(first/(1-ratio)).upper();
	res=mid+rad*interval<T>(-1.,1.);
      }
      else{
	throw std::domain_error("ratio is more than 1");
      }
    }
  }
return res;
}
  
  template <class T> complex<interval<T> >QHypergeom(const ub::vector<complex<interval<T> > >& a,const ub::vector<complex<interval<T> > >& b,const interval<T>& q,const complex<interval<T> >& z){
    complex<interval<T> >res,mid;
    interval<T>first,ratio,pro1,pro2;
    T rad;
    int r,s;
    r=a.size();
    s=b.size();
    mid=1.;
    pro1=1.;
    pro2=1.;
    int N;
    N=1000;
    for(int i=0;i<=s-1;i++){
      while(abs(b(i))>pow(1/q,N)){
	N=N+500;
      }
    }
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    if (r>s+1){
      throw std::domain_error("r>s+1 is not implemented");
    }
    if(r==s+1){
      if (abs(z)>=1){
	throw std::domain_error("absolute value of z must be under 1");
      }
      if(r==2&&s==1){
	// Heine hypergeometric function
	res=Heine(complex<interval<T> >(a(0)),complex<interval<T> >(a(1)),complex<interval<T> >(b(0)),interval<T>(q),complex<interval<T> >(z));
}
      else{
	for(int n=1;n<=N-1;n++){
	  mid=mid+qPochhammer(ub::vector<complex<interval<T> > >(a),interval<T>(q),int (n))*pow(z,n)
	    /qPochhammer(ub::vector<complex<interval<T> > >(b),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
	}
	first=abs(qPochhammer(ub::vector<complex<interval<T> > >(a),interval<T>(q),int (N))*pow(z,N)
		  /qPochhammer(ub::vector<complex<interval<T> > >(b),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
	for(int j=0;j<=s-1;j++){
	  pro1=pro1*(1+abs(b(j)-a(j))*pow(q,N)/abs(1-b(j)*pow(q,N)));
	}
	ratio=abs(z)*pro1*(1+pow(q,N)*abs(q-a(r-1))/abs(1-pow(q,N+1)));
	if(ratio<1){
	  rad=(first/(1-ratio)).upper();
	  res=complex_nbd(mid,rad);
	}
	else{
	  throw std::domain_error("ratio is more than 1");
	}
	
      }
    }
    if(r<=s){
      if(r==0&&s==1){
	res=_0phi_1(complex<interval<T> >(b(0)),interval<T>(q),complex<interval<T> >(z));
      }
      if(r==1&&s==1){
	res=_1phi_1(complex<interval<T> >(a(0)),complex<interval<T> >(b(0)),interval<T>(q),complex<interval<T> >(z));
      }
      
      else{
	for(int n=1;n<=N-1;n++){
	  mid=mid+qPochhammer(ub::vector<complex<interval<T> > >(a),interval<T>(q),int (n))*std::pow(-1,1+s-r)*pow(q,n*(n-1)*0.5*(1+s-r))*pow(z,n)
	    /qPochhammer(ub::vector<complex<interval<T> > >(b),interval<T>(q),int (n))/qPochhammer(interval<T>(q),interval<T>(q),int (n));
	}
	first=abs(qPochhammer(ub::vector<complex<interval<T> > >(a),interval<T>(q),int (N))*std::pow(-1,1+s-r)*pow(q,N*(N-1)*0.5*(1+s-r))*pow(z,N)
		  /qPochhammer(ub::vector<complex<interval<T> > >(b),interval<T>(q),int (N))/qPochhammer(interval<T>(q),interval<T>(q),int (N)));
	for(int k=0;k<=r-1;k++){
	  pro1=pro1*(1+abs(b(k)-a(k))*pow(q,N)/abs(1-b(k)*pow(q,N)));
	}
	for(int h=r;h<=s-1;h++){
	  pro2=pro2*pow(q,N*(1+s-r))/abs(1-b(h)*pow(q,N));
	}
	ratio=abs(z)*pro1*pro2/abs(1-pow(q,N));
	if(ratio<1){
	  rad=(first/(1-ratio)).upper();
	  res=complex_nbd(mid,rad);
	}
	else{
	  throw std::domain_error("ratio is more than 1");
	}
      }
    }
    return res;
  }
  template <class T> complex<interval<T> >qAppell1(const complex<interval<T> >& a,const complex<interval<T> >& b,const complex<interval<T> >& bp,const complex<interval<T> >& c ,
						     const interval<T>& q,const complex<interval<T> >& x,const complex<interval<T> >& y){
    // verification program for the first q-Appell function
    // reference:DLMF http://dlmf.nist.gov/17.11 formula 17.11.1
    complex<interval<T> >res;
    ub::vector<complex<interval<T> > >v1(3),v2(2);
    v1(0)=c/a;v1(1)=x;v1(2)=y;
    v2(0)=b*x;v2(1)=bp*y;
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    if (abs(x).upper()>=1){
      throw std::domain_error("absolute value of x must be under 1");
    }
    if (abs(y).upper()>=1){
      throw std::domain_error("absolute value of y must be under 1");
    }

    res=infinite_qPochhammer(complex<interval<T> >(a),interval<T>(q))
      *infinite_qPochhammer(complex<interval<T> >(b*x),interval<T>(q))
      *infinite_qPochhammer(complex<interval<T> >(bp*y),interval<T>(q))
      /infinite_qPochhammer(complex<interval<T> >(c),interval<T>(q))
      /infinite_qPochhammer(complex<interval<T> >(x),interval<T>(q))
      /infinite_qPochhammer(complex<interval<T> >(y),interval<T>(q))
      *QHypergeom(ub::vector<complex<interval<T> > >(v1),ub::vector<complex<interval<T> > >(v2),interval<T>(q),complex<interval<T> >(a));
    return res;
  }

  template <class T> interval<T> infinite_qPochhammer(const ub::vector<interval<T> >& a,const interval<T>& q){
    interval<T> res,pro;
    int r;
    r=a.size();
    pro=1.;
    for(int i=0;i<=r-1;i++){
      pro=pro*infinite_qPochhammer(interval<T>(a(i)),interval<T>(q));
    }
    res=pro;
    return res;
  }
  template <class T> complex<interval<T> >infinite_qPochhammer(const ub::vector<complex<interval<T> > >& a,const interval<T>& q){
    complex<interval<T> >res,pro;
    int r;
    r=a.size();
    pro=1.;
    for(int i=0;i<=r-1;i++){
      pro=pro*infinite_qPochhammer(complex<interval<T> >(a(i)),interval<T>(q));
    }
    res=pro;
    return res;
  }
  template <class T> complex<interval<T> >qLauricellaD(const complex<interval<T> >& a,const ub::vector<complex<interval<T> > >& b,const complex<interval<T> >& c ,
						     const interval<T>& q,const ub::vector<complex<interval<T> > >& x){
    // verification program for the q-Lauricella function type D
    // reference: Andrews (1972), Gasper-Rahman (2004) (Page 300, Exercise 10.17)
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    if(b.size()!=x.size()){
      throw std::domain_error("size of b and x must be same");
    }
    int r;
    r=b.size();
    for(int k=0;k<=r-1;k++){
      if(abs(x(k)).upper()>=1){
      throw std::domain_error("value of x must be under 1");
      }
    }
    complex<interval<T> >res;
    ub::vector<complex<interval<T> > >v1(r),v2(r);
    for(int i=0;i<=r-1;i++){
      v1(i)=b(i)*x(i);
    }
    v2(0)=c/a;
    for(int j=1;j<=r-1;j++){
      v2(j)=x(j);
    }
    res=infinite_qPochhammer(complex<interval<T> >(a),interval<T>(q))/infinite_qPochhammer(complex<interval<T> >(c),interval<T>(q))
      *infinite_qPochhammer(ub::vector<complex<interval<T> > >(v1),interval<T>(q))/infinite_qPochhammer(ub::vector<complex<interval<T> > >(x),interval<T>(q))
      *QHypergeom(ub::vector<complex<interval<T> > >(v2),ub::vector<complex<interval<T> > >(v1),interval<T>(q),complex<interval<T> >(a));
    return res;
  }
  template <class T> complex<interval<T> >VWP(const complex<interval<T> >& a,const ub::vector<complex<interval<T> > >& b,const interval<T>& q,const complex<interval<T> >& x){
    // verification program for the very-well-poised q-hypergeometric function
    if (q>=1){
      throw std::domain_error("value of q must be under 1");
    }
    if (q<=0){
      throw std::domain_error("q must be positive");
    }
    if (abs(x)>=1){
      throw std::domain_error("absolute value of x must be under 1");
    }
    int r;
    r=b.size()+3;
    complex<interval<T> >res;
    ub::vector<complex<interval<T> > >v1(r+1),v2(r);
    v1(0)=a;v1(1)=q*sqrt(a);v1(2)=-q*sqrt(a);
    for(int i=3;i<=r;i++){
      v1(i)=b(i-3);
    }
    v2(0)=sqrt(a);v2(1)=-sqrt(a);
    for(int j=2;j<=r-1;j++){
      v2(j)=a*q/b(j-2);
    }

    res=QHypergeom(ub::vector<complex<interval<T> > >(v1),ub::vector<complex<interval<T> > >(v2),interval<T>(q),complex<interval<T> >(x));
    return res;
  }
  template <class T> interval<T> Euler_qlog(const interval<T>& q,const interval<T> & x){
    // verification program for Euler's q-logarithm
    // reference: Koelink-Van Assche
    interval<T> res;
    ub::vector<interval<T> >v1(3),v2(2);
    v1(0)=q;v1(1)=q;v1(2)=q*x;
    v2(0)=q*q;v2(1)=0;
    res=QHypergeom(ub::vector<interval<T>  >(v1),ub::vector<interval<T>  >(v2),interval<T>(q),interval<T> (q));
    return res;
    
  }
}
#endif
