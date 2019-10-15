//#include <vcp/imats.hpp>
//#include <vcp/matrix.hpp>
//#include <vcp/matrix_assist.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include "kv/interval.hpp"
// define rounding operations for double
// if "rdouble.hpp" is not included, result of computation is not "verified" 
#include "kv/rdouble.hpp"
typedef kv::interval<double> itv;
namespace ub = boost::numeric::ublas;
itv factorial(int k){
    itv sum;
    sum=1.;
    for (int i = 1; i <= k; ++i)
    {
        sum *= i;
    }
    return sum;
}
int main(void){
std::cout.precision(17);
  int n=2,M=42;
  /*  vcp::matrix< itv, vcp::imats< double > > A,B,res,sum; 
  sum.zeros(n);
  A(0,0)=1.;
  A(0,1)=-3.;
  A(1,0)=-2.;
  A(1,1)=5.;
  for(int N=0;N<=M;N++){
    sum=sum*A+1/factorial(N);
  }
  norm=max(abs(A));
  error=exp(norm)*pow(norm,M+1)/factorial(M+1);
  b=max(abs(error));
  B.ones(n);
  B=B*itv(-b,b);
  */
  // use ublas
  ub::matrix< itv > A(n, n),B(n, n),res(n, n),sum(n, n),pro(n,n);
  itv error,norm;
  double b;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      sum(i,j)=0.;
      if(i==j)pro(i,j)=1.;
      else pro(i,j)=0.;
    }
  }
	A(0, 0)=1.; A(0, 1)=-3.;
	A(1, 0)=-2; A(1, 1)=5.;
  for(int N=0;N<=M;N++){
   sum+=(1/factorial(N))*pro;
   pro=prod(pro,A);
  }
  //fprintf( stderr, "Check\n" );
  //norm=max(abs(A));
  norm=abs(A(0,0));
  for(int i1=0;i1<n;i1++){
    for(int j1=0;j1<n;j1++){
      if (A(i1,j1)>abs(norm)) norm=abs(A(i1,j1));
    }
  }

  error=exp(norm)*pow(norm,M+1)/factorial(M+1);
  b=(abs(error)).upper();

  for(int k1=0;k1<n;k1++){
    for(int l1=0;l1<n;l1++){
      B(k1,l1).assign(-1.,1.);
      B(k1,l1)=b*B(k1,l1);
    }
  }
  res=sum+B;

  std::cout<<res<<std::endl;
  return 0;
}
