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
const AT&
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
const double&
approx(const double& d)
{
  return d;
}

inline
const float&
approx(const float& f)
{
  return f;
}

inline
const int&
approx(const int& i)
{
  return i;
}

inline
const unsigned int&
approx(const unsigned int& i)
{
  return i;
}

inline
const Null_vector&
approx(const Null_vector& nv)
{
  return nv;
}

inline
const Null_vector&
exact(const Null_vector& nv)
{
  return nv;
}

inline
const Origin&
approx(const Origin& nv)
{
  return nv;
}

inline
const Origin&
exact(const Origin& nv)
{
  return nv;
}

inline
const Orientation&
approx(const Orientation& nv)
{
  return nv;
}

inline
const Orientation&
exact(const Orientation& nv)
{
  return nv;
}

template <typename AT, typename ET, typename EFT, typename E2A>
inline
const ET&
exact(const Lazy<AT,ET,EFT,E2A>& l)
{
  return l.exact();
}

template <typename ET>
inline
const ET&
exact(const Lazy_exact_nt<ET>& l)
{
  return l.exact();
}

inline
const double&
exact(const double& d)
{
  return d;
}

inline
const float&
exact(const float& f)
{
  return f;
}

inline
const int&
exact(const int& i)
{
  return i;
}

inline
const unsigned int&
exact(const unsigned int& i)
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

inline
void
print(double d, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++){
    os << "    ";
  }
  os << d << std::endl;
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
  const typename T::AT&
  operator()(const T&t) const
  { return t.approx(); }
};

struct Exact_converter {
  template < typename T >
  const typename T::ET&
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
  
  const AT& approx() const 
  { return ptr()->approx(); }


  const ET&  exact() const
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

template <typename AK, typename EK, typename AT, typename ET, typename EFT, typename E2A>
struct Lazy_construction_return_type {
  typedef Lazy<AT, ET, EFT, E2A> result_type;
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

  template <typename L1, typename L2>
  result_type operator()(const L1& l1, const L2& l2) const
  {
    try {
      return new Lazy_construct_rep_2<AC, EC, To_interval<ET>, L1,L2>(ac, ec, l1,l2);
    } catch (Interval_nt_advanced::unsafe_comparison) {
      return new Lazy_construct_rep_0<AT,ET,To_interval<ET> >(ec(CGAL::exact(l1), exact(l2)));
    }
  }

  template <typename L1, typename L2, typename L3>
  result_type operator()(const L1& l1, const L2& l2, const L3& l3) const
  {
    try {
      return new Lazy_construct_rep_3<AC, EC, To_interval<ET>, L1,L2,L3>(ac, ec, l1,l2,l3);
    } catch (Interval_nt_advanced::unsafe_comparison) {
      return new Lazy_construct_rep_0<AT,ET,To_interval<ET> >(ec(CGAL::exact(l1), exact(l2), exact(l3)));
    }
  }
};

template <typename T>
struct Object_cast {
  typedef T result_type;

  const T&
  operator()(const Object& o) const
  {
    return object_cast<T>(o);
  }
};



// The following functor returns an Object with a Lazy<Something> inside
// As the nested kernels return Objects of AK::Something and EK::Something
// we have to unwrap them from the Object, and wrap them in a Lazy<Something>
//
// TODO: write operators for other than two arguments. For the current kernel we only need two for Interscet_2

template <typename LK, typename AK, typename EK, typename AC, typename EC, typename EFT, typename E2A>
struct Lazy_construction_object {


  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Object result_type;

  typedef Lazy<Object, Object, EFT, E2A> Lazy_object;
  AC ac;
  EC ec;

public:

  template <typename L1, typename L2>
  result_type
  operator()(const L1& l1, const L2& l2) const
  {
    try {
      Lazy_object lo(new Lazy_construct_rep_2<AC, EC, E2A, L1, L2>(ac, ec, l1, l2));
 
      if(object_cast<typename AK::Point_2>(& (lo.approx()))){
	typedef Lazy_construct_rep_1<Object_cast<typename AK::Point_2>, Object_cast<typename EK::Point_2>, E2A, Lazy_object> Lcr;
	Lcr * lcr = new Lcr(Object_cast<typename AK::Point_2>(), Object_cast<typename EK::Point_2>(), lo); 
	return make_object(typename LK::Point_2(lcr));
      } else if(object_cast<typename AK::Segment_2>(& (lo.approx()))){
	typedef Lazy_construct_rep_1<Object_cast<typename AK::Segment_2>, Object_cast<typename EK::Segment_2>, E2A, Lazy_object> Lcr;
	Lcr * lcr = new Lcr(Object_cast<typename AK::Segment_2>(), Object_cast<typename EK::Segment_2>(), lo); 
	return make_object(typename LK::Segment_2(lcr));
      } else if(object_cast<typename AK::Ray_2>(& (lo.approx()))){
	typedef Lazy_construct_rep_1<Object_cast<typename AK::Ray_2>, Object_cast<typename EK::Ray_2>, E2A, Lazy_object> Lcr;
	Lcr * lcr = new Lcr(Object_cast<typename AK::Ray_2>(), Object_cast<typename EK::Ray_2>(), lo); 
	return make_object(typename LK::Ray_2(lcr));
      } else if(object_cast<typename AK::Line_2>(& (lo.approx()))){
	typedef Lazy_construct_rep_1<Object_cast<typename AK::Line_2>, Object_cast<typename EK::Line_2>, E2A, Lazy_object> Lcr;
	Lcr * lcr = new Lcr(Object_cast<typename AK::Line_2>(), Object_cast<typename EK::Line_2>(), lo); 
	return make_object(typename LK::Line_2(lcr));
      }  else if(object_cast<typename AK::Circle_2>(& (lo.approx()))){
	typedef Lazy_construct_rep_1<Object_cast<typename AK::Circle_2>, Object_cast<typename EK::Circle_2>, E2A, Lazy_object> Lcr;
	Lcr * lcr = new Lcr(Object_cast<typename AK::Circle_2>(), Object_cast<typename EK::Circle_2>(), lo); 
	return make_object(typename LK::Circle_2(lcr));
      } else if(object_cast<typename AK::Iso_rectangle_2>(& (lo.approx()))){
	typedef Lazy_construct_rep_1<Object_cast<typename AK::Iso_rectangle_2>, Object_cast<typename EK::Iso_rectangle_2>, E2A, Lazy_object> Lcr;
	Lcr * lcr = new Lcr(Object_cast<typename AK::Iso_rectangle_2>(), Object_cast<typename EK::Iso_rectangle_2>(), lo); 
	return make_object(typename LK::Iso_rectangle_2(lcr));
      }  else if(object_cast<typename AK::Triangle_2>(& (lo.approx()))){
	typedef Lazy_construct_rep_1<Object_cast<typename AK::Triangle_2>, Object_cast<typename EK::Triangle_2>, E2A, Lazy_object> Lcr;
	Lcr * lcr = new Lcr(Object_cast<typename AK::Triangle_2>(), Object_cast<typename EK::Triangle_2>(), lo); 
	return make_object(typename LK::Triangle_2(lcr));
      }   else {
	std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's" << std::endl;
      }
    } catch (Interval_nt_advanced::unsafe_comparison) {
      ET eto = ec(CGAL::exact(l1), CGAL::exact(l2));

      if(const typename EK::Point_2* ptr = object_cast<typename EK::Point_2>(&eto)){
	make_object(typename LK::Point_2(new Lazy_construct_rep_0<typename AK::Point_2, typename EK::Point_2, E2A>(*ptr)));
      } else if(const typename EK::Segment_2* ptr = object_cast<typename EK::Segment_2>(&eto)){
	make_object(typename LK::Segment_2(new Lazy_construct_rep_0<typename AK::Segment_2, typename EK::Segment_2, E2A>(*ptr)));
      } else if(const typename EK::Ray_2* ptr = object_cast<typename EK::Ray_2>(&eto)){
	make_object(typename LK::Ray_2(new Lazy_construct_rep_0<typename AK::Ray_2, typename EK::Ray_2, E2A>(*ptr)));
      } else if(const typename EK::Line_2* ptr = object_cast<typename EK::Line_2>(&eto)){
	make_object(typename LK::Line_2(new Lazy_construct_rep_0<typename AK::Line_2, typename EK::Line_2, E2A>(*ptr)));
      } else if(const typename EK::Triangle_2* ptr = object_cast<typename EK::Triangle_2>(&eto)){
	make_object(typename LK::Triangle_2(new Lazy_construct_rep_0<typename AK::Triangle_2, typename EK::Triangle_2, E2A>(*ptr)));
      } else if(const typename EK::Iso_rectangle_2* ptr = object_cast<typename EK::Iso_rectangle_2>(&eto)){
	make_object(typename LK::Iso_rectangle_2(new Lazy_construct_rep_0<typename AK::Iso_rectangle_2, typename EK::Iso_rectangle_2, E2A>(*ptr)));
      } else if(const typename EK::Circle_2* ptr = object_cast<typename EK::Circle_2>(&eto)){
	make_object(typename LK::Circle_2(new Lazy_construct_rep_0<typename AK::Circle_2, typename EK::Circle_2, E2A>(*ptr)));
      } else{
	std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's" << std::endl;
      }      
    }
    return Object();
  }
  

};

template <typename T>
struct First {
  typedef typename T::first_type result_type;

  const typename T::first_type&
  operator()(const T& p) const
  {
    return p.first;
  }
};

template <typename T>
struct Second {
  typedef typename T::second_type result_type;

  const typename T::second_type&
  operator()(const T& p) const
  {
    return p.second;
  }
};

template <typename AK, typename EK, typename AC, typename EC, typename EFT, typename E2A>
struct Pair_of_lazy_construction {


  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy<AT, ET, EFT, E2A> Lazy_pair;

  typedef Lazy<typename AT::first_type, typename ET::first_type, EFT, E2A> Lazy_pair_first_type;
  typedef Lazy<typename AT::second_type, typename ET::second_type, EFT, E2A> Lazy_pair_second_type;

  typedef std::pair<Lazy_pair_first_type, Lazy_pair_second_type > Pair_of_lazy; 
  typedef Pair_of_lazy result_type;
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
      Lazy_pair lazy_pair(new Lazy_construct_rep_2<AC, EC, E2A, L1, L2>(ac, ec, l1, l2));
      Lazy_pair_first_type first(new Lazy_construct_rep_1<First<AT>, First<ET>, E2A, Lazy_pair>(First<AT>(), First<ET>(), lazy_pair)); 
      Lazy_pair_second_type second(new Lazy_construct_rep_1<Second<AT>, Second<ET>, E2A, Lazy_pair>(Second<AT>(), Second<ET>(), lazy_pair));
      return std::make_pair(first, second);
    } catch (Interval_nt_advanced::unsafe_comparison) {
      ET et = ec(l1.exact(), l2.exact());
      return std::make_pair(Lazy_pair_first_type(new Lazy_construct_rep_0<typename AT::first_type, typename ET::first_type, E2A>(et.first)),
			    Lazy_pair_second_type(new Lazy_construct_rep_0<typename AT::second_type, typename ET::second_type, E2A>(et.second))) ;
    }
  }

};


//____________________________________________________________
// The functor that has Lazy<Something> as result type

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




