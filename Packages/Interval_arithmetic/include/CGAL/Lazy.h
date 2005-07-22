// Copyright (c) 2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_LAZY_H
#define CGAL_LAZY_H

#include <CGAL/basic.h>
#include <CGAL/Handle.h>
#include <CGAL/Object.h>
#include <CGAL/Lazy_exact_nt.h>
#include <boost/static_assert.hpp>
CGAL_BEGIN_NAMESPACE

template <typename AT, typename ET, typename EFT, typename E2A> class Lazy;
template <typename ET> class Lazy_exact_nt;


template <typename AT, typename ET, typename EFT, typename E2A>
inline
AT
approx(const Lazy<AT,ET, EFT, E2A>& l)
{
  return l.approx();
}

template <typename ET>
inline
const Interval_nt<true>&
approx(const Lazy_exact_nt<ET>& l)
{
  return l.approx();
}


inline
double
approx(double d)
{
  return d;
}

inline
float
approx(float f)
{
  return f;
}

inline
int
approx(int i)
{
  return i;
}

inline
unsigned int
approx(unsigned int i)
{
  return i;
}

inline
Null_vector
approx(const Null_vector nv)
{
  return nv;
}

inline
Null_vector
exact(const Null_vector nv)
{
  return nv;
}

inline
Origin
approx(const Origin nv)
{
  return nv;
}

inline
Origin
exact(const Origin nv)
{
  return nv;
}

inline
Orientation
approx(const Orientation nv)
{
  return nv;
}

inline
Orientation
exact(const Orientation nv)
{
  return nv;
}

template <typename AT, typename ET, typename EFT, typename E2A>
inline
ET
exact(const Lazy<AT,ET,EFT,E2A>& l)
{
  return l.exact();
}

template <typename ET>
inline
ET
exact(const Lazy_exact_nt<ET>& l)
{
  return l.exact();
}

inline
double
exact(double d)
{
  return d;
}

inline
float
exact(float f)
{
  return f;
}

inline
int
exact(int i)
{
  return i;
}

inline
unsigned int
exact(unsigned int i)
{
  return i;
}


template <typename AT, typename ET, typename EFT, typename E2A>
inline
void
print(const Lazy<AT,ET,EFT,E2A>& l, std::ostream& os, int level)
{
  l.print(os, level);
}

template <typename ET>
inline
void
print(const Lazy_exact_nt<ET>& l, std::ostream& os, int level)
{
  //os << "TBD" << std::endl;
  l.print(os, level);
}

inline
void
print(double d, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++){
    os << "    ";
  }
  os << d << std::endl;
}


void
msg(std::ostream& os, int level, char* s)
  {
    int i;
    for(i = 0; i < level; i++){
      os << "    ";
    }
    os << s << std::endl;
  }

inline
void
print(const Null_vector& nv, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++){
    os << "    ";
  }
  os << "Null_vector" << std::endl;
}

inline
void
print(const Origin& nv, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++){
    os << "    ";
  }
  os << "Origin" << std::endl;
}




//____________________________________________________________
// The rep for the leaf node

template <typename AT, typename ET, typename E2A> 
class Lazy_construct_rep_0 : public Lazy_construct_rep<AT, ET, E2A>
{

  typedef Lazy_construct_rep<AT, ET, E2A> Base;
public:

  void
  update_exact()
  {
    this->et = new ET();
  }

  Lazy_construct_rep_0()
    : Lazy_construct_rep<AT,ET, E2A>()
  {}

  Lazy_construct_rep_0(const AT& a, const ET& e)
    : Lazy_construct_rep<AT,ET,E2A>(a, e)
  {}
  
  Lazy_construct_rep_0(const AT& a, void*)
    : Lazy_construct_rep<AT,ET,E2A>(a)
  {}

  Lazy_construct_rep_0(const ET& e)
    : Lazy_construct_rep<AT,ET,E2A>(AT(), e)
  {
    E2A e2a;
    this->at = e2a(e); 
  }

  void
  print(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
  }
};


//____________________________________________________________

template <typename AC, typename EC, typename E2A, typename L1>
class Lazy_construct_rep_1 : public Lazy_construct_rep<typename AC::result_type , typename EC::result_type, E2A>
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_construct_rep<AT, ET, E2A> Base;

  EC ec_;
  L1 l1_;

public:

  void
  update_exact()
  {
    this->et = new ET(ec_(CGAL::exact(l1_)));
    // Prune lazy tree
    l1_ = L1(); 
  }


  Lazy_construct_rep_1(const AC& ac, const EC& ec, const L1& l1)
    : Lazy_construct_rep<AT,ET, E2A>(ac(CGAL::approx(l1))), ec_(ec), l1_(l1)
  {}

  void
  print(std::ostream& os, int level) const 
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      CGAL::msg(os, level, "One child node:");
      CGAL::print(l1_, os, level+1);
    }
  }
};



//____________________________________________________________

template <typename AC, typename EC, typename E2A, typename L1, typename L2>
class Lazy_construct_rep_2 : public Lazy_construct_rep<typename AC::result_type , typename EC::result_type, E2A>
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_construct_rep<AT, ET, E2A> Base;

  EC ec_;
  L1 l1_;
  L2 l2_;

public:

  void
  update_exact()
  {
    this->et = new ET(ec_(CGAL::exact(l1_), CGAL::exact(l2_)));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
  }


  Lazy_construct_rep_2(const AC& ac, const EC& ec, const L1& l1, const L2& l2)
    : Lazy_construct_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2))), l1_(l1), l2_(l2)
  {}

  void
  print(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      CGAL::msg(os, level, "Two child nodes:");
      CGAL::print(l1_, os, level+1);
      CGAL::print(l2_, os, level+1);
    }
  }
};



//____________________________________________________________

template <typename AC, typename EC, typename E2A, typename L1, typename L2, typename L3>
class Lazy_construct_rep_3 : public Lazy_construct_rep<typename AC::result_type , typename EC::result_type, E2A>
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_construct_rep<AT, ET, E2A> Base;

  EC ec_;
  L1 l1_;
  L2 l2_;
  L3 l3_;

public:

  void
  update_exact()
  {
    this->et = new ET(ec_(CGAL::exact(l1_), CGAL::exact(l2_), CGAL::exact(l3_)));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
    l3_ = L3();
  }


  Lazy_construct_rep_3(const AC& ac, const EC& ec, const L1& l1, const L2& l2, const L3& l3)
    : Lazy_construct_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2), CGAL::approx(l3))), l1_(l1), l2_(l2), l3_(l3)
  {}

  void
  print(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      CGAL::msg(os, level, "Three child nodes:");
      CGAL::print(l1_, os, level+1);
      CGAL::print(l2_, os, level+1);
      CGAL::print(l3_, os, level+1);
    }
  }
};


//____________________________________________________________

template <typename AC, typename EC, typename E2A, typename L1, typename L2, typename L3, typename L4>
class Lazy_construct_rep_4 : public Lazy_construct_rep<typename AC::result_type , typename EC::result_type, E2A>
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_construct_rep<AT, ET, E2A> Base;

  EC ec_;
  L1 l1_;
  L2 l2_;
  L3 l3_;
  L4 l4_;

public:

  void
  update_exact()
  {
    this->et = new ET(ec_(CGAL::exact(l1_), CGAL::exact(l2_), CGAL::exact(l3_), CGAL::exact(l4_)));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
    l3_ = L3();
    l4_ = L4();
  }


  Lazy_construct_rep_4(const AC& ac, const EC& ec, const L1& l1, const L2& l2, const L3& l3, const L4& l4)
    : Lazy_construct_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2), CGAL::approx(l3), CGAL::approx(l4))), l1_(l1), l2_(l2), l3_(l3), l4_(l4)
  {}

  void
  print(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    
    if(this->is_lazy()){
      CGAL::msg(os, level, "Four child nodes:");
      CGAL::print(l1_, os, level+1);
      CGAL::print(l2_, os, level+1);
      CGAL::print(l3_, os, level+1);
      CGAL::print(l4_, os, level+1);
    }
  }
};

//____________________________________________________________


template <typename AC, typename EC, typename E2A, typename L1, typename L2, typename L3, typename L4, typename L5>
 class Lazy_construct_rep_5 : public Lazy_construct_rep<typename AC::result_type , typename EC::result_type, E2A>
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_construct_rep<AT, ET, E2A> Base;

  EC ec_;
  L1 l1_;
  L2 l2_;
  L3 l3_;
  L4 l4_;
  L5 l5_;

public:

  void
  update_exact()
  {
    this->et = new ET(ec_(CGAL::exact(l1_), CGAL::exact(l2_), CGAL::exact(l3_), CGAL::exact(l4_), CGAL::exact(l5_)));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
    l3_ = L3();
    l4_ = L4();
    l5_ = L5();
  }

 Lazy_construct_rep_5(const AC& ac, const EC& ec, const L1& l1, const L2& l2, const L3& l3, const L4& l4, const L5& l5)
    : Lazy_construct_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2), CGAL::approx(l3), CGAL::approx(l4), CGAL::approx(l5))), l1_(l1), l2_(l2), l3_(l3), l4_(l4), l5_(l5)
  {}

  void
  print(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    
    if(this->is_lazy()){
      CGAL::msg(os, level, "Five child nodes:");
      CGAL::print(l1_, os, level+1);
      CGAL::print(l2_, os, level+1);
      CGAL::print(l3_, os, level+1);
      CGAL::print(l4_, os, level+1);
      CGAL::print(l5_, os, level+1);
    }
  }
};

struct Approx_converter {
  template < typename T >
  const typename T::AT
  operator()(const T&t) const
  { return t.approx(); }
};

struct Exact_converter {
  template < typename T >
  const typename T::ET
  operator()(const T&t) const
  { return t.exact(); }
};


//____________________________________________________________
// The handle class
template <typename AT_, typename ET_, typename EFT, typename E2A>
class Lazy : public Handle
{
public :
  typedef AT_ AT;
  typedef ET_ ET;
  typedef Lazy<AT, ET, EFT, E2A> Self;
  typedef Lazy_construct_rep<AT, ET, E2A> Self_rep;

  typedef Self Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  Lazy()
  {
    PTR = new Lazy_construct_rep_0<AT,ET, E2A>();
  }
  
  Lazy (Self_rep *r)
  { 
    PTR = r; 
  }
  
  AT approx() const 
  { return ptr()->approx(); }


  ET  exact() const
  { return ptr()->exact(); }

  void
  print(std::ostream& os, int level=1) const
  {
    ptr()->print(os, level);
  }

private:
  Self_rep * ptr() const { return (Self_rep*) PTR; }

};


template <typename AT, typename ET, typename EFT, typename E2A>
std::ostream&
operator<<(std::ostream& os, const Lazy<AT,ET,EFT, E2A>& lazy)
{
  if(is_pretty(os)){
    lazy.print(os);
  } else {
    os << lazy.approx();
  }
  return os;
} 

template <typename AT, typename ET, typename EFT, typename E2A>
bool
operator==(const Lazy<AT,ET,EFT,E2A>& a, const Lazy<AT,ET,EFT,E2A>& b)
{
  try
  {
    return a.approx() == b.approx();
  }
  catch (Interval_nt<false>::unsafe_comparison)
  {
    // std::cerr << "Interval filter failure (==)" << std::endl;
    return a.exact() == b.exact();
  }
}

template <typename AT, typename ET, typename EFT, typename E2A>
bool
operator!=(const Lazy<AT,ET,EFT,E2A>& a, const Lazy<AT,ET,EFT,E2A>& b)
{
  return ! (a == b);
}





//____________________________________________________________
// A helper class to select the return type.

// AF: Regarder le tutorial metaprogramming

template <typename AK, typename EK, typename AT, typename ET, typename EFT, typename E2A>
struct Lazy_construction_return_type {
  typedef Lazy<AT, ET, EFT, E2A> result_type;
};




// As people will write bool b = lazy_assign(..) 
// Or are there cases where we need a lazy bool in the kernel???
template <typename AK, typename EK, typename ET, typename EFT, typename E2A>
struct Lazy_construction_return_type<AK,EK,bool, ET,EFT,E2A>  {
  typedef bool  result_type;
};


template <typename AK, typename EK, typename AC, typename EC, typename EFT, typename E2A >
struct Lazy_construction_bbox {
  typedef typename AC::result_type result_type;

  AC ac;
  EC ec;
  template <typename L1>
  result_type operator()(const L1& l1) const
  {
    try {
      return ac(CGAL::approx(l1));
    } catch (Interval_nt_advanced::unsafe_comparison) {
      return ec(CGAL::exact(l1));
    }
  }
};


template <typename AK, typename EK, typename AC, typename EC, typename EFT, typename E2A >
struct Lazy_construction_nt {

  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_exact_nt<ET> result_type;

  AC ac;
  EC ec;
  template <typename L1>
  result_type operator()(const L1& l1) const
  {
    try {
      return new Lazy_construct_rep_1<AC, EC, To_interval<ET>, L1>(ac, ec, l1);
    } catch (Interval_nt_advanced::unsafe_comparison) {
      return new Lazy_construct_rep_0<AT,ET,To_interval<ET> >(ec(CGAL::exact(l1)));
    }
  }
};




//____________________________________________________________
// The functor that creates all kinds of Lazy<T>

template <typename AK, typename EK, typename AC, typename EC, typename EFT, typename E2A>
struct Lazy_construction {


  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef typename Lazy_construction_return_type<AK, EK, AT, ET, EFT, E2A>::result_type result_type;

  AC ac;
  EC ec;

public:

  result_type
  operator()() const
  {
    return new Lazy_construct_rep_0<AT,ET,E2A>();
  }


  template <typename L1>
  result_type
  operator()(const L1& l1) const
  {
    try {
      return  new Lazy_construct_rep_1<AC, EC, E2A, L1>(ac, ec, l1);
    } catch (Interval_nt_advanced::unsafe_comparison) {
      return new Lazy_construct_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1)));
    }
  }

  template <typename L1, typename L2>
  result_type
  operator()(const L1& l1, const L2& l2) const
  {
    try {
      return new Lazy_construct_rep_2<AC, EC, E2A, L1, L2>(ac, ec, l1, l2);
    } catch (Interval_nt_advanced::unsafe_comparison) {
      return new Lazy_construct_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2)));
    }
  }
  

  template <typename L1, typename L2, typename L3>
  result_type
  operator()(const L1& l1, const L2& l2, const L3& l3) const
  {
    try {
      return new Lazy_construct_rep_3<AC, EC, E2A, L1, L2, L3>(ac, ec, l1, l2, l3);
    } catch (Interval_nt_advanced::unsafe_comparison) {
      return new Lazy_construct_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3)));
    }
  }

  template <typename L1, typename L2, typename L3, typename L4>
  result_type
  operator()(const L1& l1, const L2& l2, const L3& l3, const L4& l4) const
  {
    try {
    return new Lazy_construct_rep_4<AC, EC, E2A, L1, L2, L3, L4>(ac, ec, l1, l2, l3, l4);
    } catch (Interval_nt_advanced::unsafe_comparison) {
      return new Lazy_construct_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3), CGAL::exact(l4)));
    }
  }

  template <typename L1, typename L2, typename L3, typename L4, typename L5>
  result_type
  operator()(const L1& l1, const L2& l2, const L3& l3, const L4& l4, const L5& l5) const
  {
    try {
    return new Lazy_construct_rep_5<AC, EC, E2A, L1, L2, L3, L4, L5>(ac, ec, l1, l2, l3, l4, l5);
    } catch (Interval_nt_advanced::unsafe_comparison) {
      return new Lazy_construct_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3), CGAL::exact(l4), CGAL::exact(l5)));
    }
  }
};


CGAL_END_NAMESPACE


#endif // CGAL_LAZY_H




