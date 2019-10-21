// Author: Daisuke Kanaizumi
// Affiliation: Department of Applied Mathematics, Waseda University

// verification program for q-series
// reference: Wolfram Mathworld

#ifndef QSERIES_HPP
#define QSERIES_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/Pochhammer.hpp>
#include <cmath>

namespace kv{
  template <class T> interval<T>Rogers_Ramanujan1(const interval<T>& q){
    interval<T>res;
    res=1/infinite_qPochhammer(interval<T>(q),interval<T>(pow(q,5)))
      /infinite_qPochhammer(interval<T>(pow(q,4)),interval<T>(pow(q,5)));
    return res;
  }

  template <class T> interval<T>Rogers_Ramanujan2(const interval<T>& q){
    interval<T>res;
    res=1/infinite_qPochhammer(interval<T>(q*q),interval<T>(pow(q,5)))
      /infinite_qPochhammer(interval<T>(pow(q,3)),interval<T>(pow(q,5)));
    return res;
  }
  template <class T> interval<T>Rogers_Ramanujan_continued_fraction(const interval<T>& q){
    interval<T>res;
    res=pow(q,0.2)*Rogers_Ramanujan2(interval<T>(q))/Rogers_Ramanujan1(interval<T>(q));
    return res;
}
  template <class T> interval<T>Rogers_Selberg1(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(pow(q,3)),interval<T>(pow(q,7)))
      *infinite_qPochhammer(interval<T>(pow(q,4)),interval<T>(pow(q,7)))
      *Euler(interval<T>(pow(q,7)))/Euler(interval<T>(q*q));
    return res;
  }
  template <class T> interval<T>Rogers_Selberg2(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(pow(q,2)),interval<T>(pow(q,7)))
      *infinite_qPochhammer(interval<T>(pow(q,5)),interval<T>(pow(q,7)))
      *Euler(interval<T>(pow(q,7)))/Euler(interval<T>(q*q));
    return res;
  }
  template <class T> interval<T>Rogers_Selberg3(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(q),interval<T>(pow(q,7)))
      *infinite_qPochhammer(interval<T>(pow(q,6)),interval<T>(pow(q,7)))
      *Euler(interval<T>(pow(q,7)))/Euler(interval<T>(q*q));
    return res;
  }
  template <class T> interval<T>Jackson_Slater(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(q),interval<T>(pow(q,8)))
      *infinite_qPochhammer(interval<T>(pow(q,7)),interval<T>(pow(q,8)))
      *Euler(interval<T>(pow(q,8)))
      *infinite_qPochhammer(interval<T>(pow(q,6)),interval<T>(pow(q,16)))
      *infinite_qPochhammer(interval<T>(pow(q,10)),interval<T>(pow(q,16)))
      /Euler(interval<T>(q));
    return res;
  }
  template <class T> interval<T>Gollnitz_Gordon1(const interval<T>& q){
    interval<T>res;
    res=1/infinite_qPochhammer(interval<T>(q),interval<T>(pow(q,8)))
      /infinite_qPochhammer(interval<T>(pow(q,4)),interval<T>(pow(q,8)))
      /infinite_qPochhammer(interval<T>(pow(q,7)),interval<T>(pow(q,8)));
    return res;
  }
  template <class T> interval<T>Gollnitz_Gordon2(const interval<T>& q){
    interval<T>res;
    res=1/infinite_qPochhammer(interval<T>(pow(q,3)),interval<T>(pow(q,8)))
      /infinite_qPochhammer(interval<T>(pow(q,4)),interval<T>(pow(q,8)))
      /infinite_qPochhammer(interval<T>(pow(q,5)),interval<T>(pow(q,8)));
    return res;
  }
  template <class T> interval<T>Bailey_Mod9_1(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(pow(q,4)),interval<T>(pow(q,9)))
      *infinite_qPochhammer(interval<T>(pow(q,5)),interval<T>(pow(q,9)))
      *Euler(interval<T>(pow(q,9)))/Euler(interval<T>(pow(q,3)));
    return res;
  }
  template <class T> interval<T>Bailey_Mod9_2(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(q*q),interval<T>(pow(q,9)))
      *infinite_qPochhammer(interval<T>(pow(q,7)),interval<T>(pow(q,9)))
      *Euler(interval<T>(pow(q,9)))/Euler(interval<T>(pow(q,3)));
    return res;
  }
  template <class T> interval<T>Bailey_Mod9_3(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(q),interval<T>(pow(q,9)))
      *infinite_qPochhammer(interval<T>(pow(q,8)),interval<T>(pow(q,9)))
      *Euler(interval<T>(pow(q,9)))/Euler(interval<T>(pow(q,3)));
    return res;
  }
  template <class T> interval<T>Rogers_Mod14_1(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(pow(q,6)),interval<T>(pow(q,14)))
      *infinite_qPochhammer(interval<T>(pow(q,8)),interval<T>(pow(q,14)))
      *Euler(interval<T>(pow(q,14)))/Euler(interval<T>(q));
    return res;
  }
  template <class T> interval<T>Rogers_Mod14_2(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(pow(q,4)),interval<T>(pow(q,14)))
      *infinite_qPochhammer(interval<T>(pow(q,10)),interval<T>(pow(q,14)))
      *Euler(interval<T>(pow(q,14)))/Euler(interval<T>(q));
    return res;
  }
  template <class T> interval<T>Rogers_Mod14_3(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(pow(q,2)),interval<T>(pow(q,14)))
      *infinite_qPochhammer(interval<T>(pow(q,12)),interval<T>(pow(q,14)))
      *Euler(interval<T>(pow(q,14)))/Euler(interval<T>(q));
    return res;
  }
  template <class T> interval<T>Dyson_Mod27_1(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(pow(q,12)),interval<T>(pow(q,27)))
      *infinite_qPochhammer(interval<T>(pow(q,15)),interval<T>(pow(q,27)))
      *Euler(interval<T>(pow(q,27)))/Euler(interval<T>(q));
    return res;
  }
  template <class T> interval<T>Dyson_Mod27_2(const interval<T>& q){
    interval<T>res;
    res=Euler(interval<T>(pow(q,9)))/Euler(interval<T>(q));
    return res;
  }
  template <class T> interval<T>Dyson_Mod27_3(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(pow(q,6)),interval<T>(pow(q,27)))
      *infinite_qPochhammer(interval<T>(pow(q,21)),interval<T>(pow(q,27)))
      *Euler(interval<T>(pow(q,27)))/Euler(interval<T>(q));
    return res;
  }
  template <class T> interval<T>Dyson_Mod27_4(const interval<T>& q){
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(pow(q,3)),interval<T>(pow(q,27)))
      *infinite_qPochhammer(interval<T>(pow(q,24)),interval<T>(pow(q,27)))
      *Euler(interval<T>(pow(q,27)))/Euler(interval<T>(q));
    return res;
  }
  template <class T> interval<T>Gessel_Stanton(const interval<T>& q){
    // reference: Savage, Sills, 2009
    // On an identity of Gessel and Stanton and new little Gollnitz identities 
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(-pow(q,3)),interval<T>(pow(q,8)))
      *infinite_qPochhammer(interval<T>(-pow(q,5)),interval<T>(pow(q,8)))
      *infinite_qPochhammer(interval<T>(-q*q),interval<T>(q*q));
    return res;
  }
  template <class T> interval<T>Lebesgue(const interval<T>& a,const interval<T>& q){
    // reference: Savage, Sills, 2009
    // On an identity of Gessel and Stanton and new little Gollnitz identities 
    interval<T>res;
    res=infinite_qPochhammer(interval<T>(a*q),interval<T>(q*q))
      /infinite_qPochhammer(interval<T>(q),interval<T>(q*q));
    return res;
  }
}
#endif
