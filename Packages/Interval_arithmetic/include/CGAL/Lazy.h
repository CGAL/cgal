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

CGAL_BEGIN_NAMESPACE

template <typename AT, typename ET, typename EFT, typename E2A> class Lazy;
template <typename AK, typename EK, typename E2A> class Lazy_exact_nt;


template <typename AT, typename ET, typename EFT, typename E2A>
inline
AT
approx(const Lazy<AT,ET, EFT, E2A>& l)
{
  return l.approx();
}

template <typename AK, typename EK, typename E2A>
inline
typename AK::FT
approx(const Lazy_exact_nt<AK,EK, E2A>& l)
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

template <typename AK, typename EK, typename E2A>
inline
typename EK::FT
exact(const Lazy_exact_nt<AK,EK,E2A>& l)
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

template <typename AK, typename EK, typename E2A>
inline
void
print(const Lazy_exact_nt<AK,EK,E2A>& l, std::ostream& os, int level)
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
// Abstract base class
template <typename AT, typename ET, typename E2A>
struct Lazy_construct_rep : public Rep
{
  AT at;
  ET *et;

  Lazy_construct_rep ()
      : at(), et(NULL) {}

  Lazy_construct_rep (const AT& a)
      : at(a), et(NULL) 
  {}

  Lazy_construct_rep (const AT& a, const ET& e)
      : at(a), et(new ET(e)) 
  {}


private:
  Lazy_construct_rep (const Lazy_construct_rep&) { abort(); } // cannot be copied.
public:

  AT& approx() //  This is not const AT&, because otherwise K::Assign_2::operator(T& t, const Object& o) does not compile 
  {
      return at;
  }

  ET & exact()
  {
      if (et==NULL) {
          update_exact();
	  E2A e2a;
	  e2a(*et); // improve approximation
      }
      return *et;
  }

  bool is_lazy() const
  {
    return et == NULL;
  }


  virtual
  void
  print(std::ostream& os, int level) const = 0;

  void
  print_at_et(std::ostream& os, int level) const
  {
    int i;
    for(i = 0; i < level; i++){
      os << "    ";
    }
    os << "Approximation: " << at << std::endl;
    if(et!= NULL){
      for(i = 0; i < level; i++){
	os << "    ";
      }
      os << "Exact        : " << *et << std::endl;
    }
  }

  virtual void update_exact() = 0;
  virtual ~Lazy_construct_rep () { delete et; };
};


//_____________________________________________________________
// the base class for lazy numbers
template <typename ET, typename E2A>
struct Lazy_exact_rep : public Lazy_construct_rep<Interval_nt<true>, ET, E2A >
{
  typedef Lazy_construct_rep<Interval_nt<true>, ET, E2A > Base;

  Lazy_exact_rep (const Interval_nt<true> & i)
      : Base(i) {}

private:
  Lazy_exact_rep (const Lazy_exact_rep&) { abort(); } // cannot be copied.
public:

  virtual void update_exact(){}

  void
  print(std::ostream& os, int level) const
  {
    print_at_et(os, level);
  }
};

// int constant
template <typename ET, typename E2A>
struct Lazy_exact_Int_Cst : public Lazy_exact_rep<ET,E2A>
{
  typedef Lazy_exact_rep<ET,E2A> Base;

  Lazy_exact_Int_Cst (int i)
      : Lazy_exact_rep<ET,E2A>(double(i)) {}

  void update_exact()  { this->et = new ET((int)this->at.inf()); }

  void
  print(std::ostream& os, int level) const
  {
    print_at_et(os, level);
  }
};

// double constant
template <typename ET, typename E2A>
struct Lazy_exact_Cst : public Lazy_exact_rep<ET,E2A>
{
  typedef Lazy_exact_rep<ET,E2A> Base;

  Lazy_exact_Cst (double d)
      : Lazy_exact_rep<ET,E2A>(d) {}

  void update_exact()  { this->et = new ET(this->at.inf()); }

  void
  print(std::ostream& os, int level) const
  {
    print_at_et(os, level);
  }
};

// Exact constant
template <typename ET, typename E2A>
struct Lazy_exact_Ex_Cst : public Lazy_exact_rep<ET,E2A>
{
  typedef Lazy_exact_rep<ET,E2A> Base;
  Lazy_exact_Ex_Cst (const ET & e)
      : Lazy_exact_rep<ET,E2A>(to_interval(e))
  {
    this->et = new ET(e);
  }

  void update_exact()  { CGAL_assertion(false); }

  void
  print(std::ostream& os, int level) const
  {
    print_at_et(os, level);
  }
};

/*
// Construction from a Lazy_exact_nt<ET1> (which keeps the lazyness).
template <typename ET, typename ET1>
struct Lazy_lazy_exact_Cst : public Lazy_exact_rep<ET>
{
  Lazy_lazy_exact_Cst (const Lazy_exact_nt<ET1> &x)
      : Lazy_exact_rep<ET>(x.approx()), l(x) {}

  void update_exact()  { this->et = new ET(l.exact()); }

  Lazy_exact_nt<ET1> l;
};
*/
// Unary  operations: abs, sqrt, square.
// Binary operations: +, -, *, /, min, max.

// Base unary operation
template <typename AK, typename EK, typename E2A>
struct Lazy_exact_unary : public Lazy_exact_rep<typename EK::FT,E2A>
{
  typedef Lazy_exact_rep<typename EK::FT,E2A> Base;
  typedef Lazy_exact_nt<AK,EK,E2A> L;
  const  L op1;


  Lazy_exact_unary (const Interval_nt<true> &i, const L  &a)
      : Lazy_exact_rep<typename EK::FT,E2A>(i), op1(a) 
  {}

  void
  print(std::ostream& os, int level) const 
  {
    print_at_et(os, level);
    if(is_lazy()){
      CGAL::msg(os, level, "One child node:");
      op1.print(os, level+1);
    }
  }
};

// Base binary operation
template <typename AK, typename EK, typename E2A>
struct Lazy_exact_binary : public Lazy_exact_unary<AK,EK,E2A>
{
  typedef Lazy_exact_rep<typename EK::FT,E2A> BaseBase;
  typedef Lazy_exact_nt<AK,EK,E2A> L;
  const L op2;

  Lazy_exact_binary (const Interval_nt<true> &i,
		     const L &a, const L &b)
      : Lazy_exact_unary<AK,EK,E2A>(i, a), op2(b) 
  {}

  void
  print(std::ostream& os, int level) const 
  {
    print_at_et(os, level);
    if(is_lazy()){
      CGAL::msg(os, level, "Two child nodes:");
      op1.print(os, level+1);
      op2.print(os, level+1);
    }
  }
};

// Here we could use a template class for all operations (STL provides
// function objects plus, minus, multiplies, divides...).  But it would require
// a template template parameter, and GCC 2.95 seems to crash easily with them.

// Macro for unary operations
#define CGAL_LAZY_UNARY_OP(OP, NAME)                                  \
template <typename AK, typename EK, typename E2A>                     \
struct NAME : public Lazy_exact_unary<AK,EK,E2A>                      \
{                                                                     \
  typedef Lazy_exact_rep<typename EK::FT,E2A> BaseBase;               \
  typedef Lazy_exact_nt<AK,EK,E2A> L;                                 \
  NAME (const L &a)                                                  \
      : Lazy_exact_unary<AK,EK,E2A>(OP(a.approx()), a) {}              \
                                                                     \
  void                                                                \
  print(std::ostream& os, int level) const                           \
  {                                                                  \
    print_at_et(os, level);                                      \
    if(is_lazy()){                                                   \
      CGAL::msg(os, level, "One child node:");                 \
      op1.print(os, level+1);                                 \
    }                                                                \
  }                                                                   \
  void update_exact()  { this->et = new typename EK::FT(OP(this->op1.exact())); } \
};



CGAL_LAZY_UNARY_OP(CGAL::opposite,  Lazy_exact_Opp)
CGAL_LAZY_UNARY_OP(CGAL_NTS abs,    Lazy_exact_Abs)
CGAL_LAZY_UNARY_OP(CGAL_NTS square, Lazy_exact_Square)
CGAL_LAZY_UNARY_OP(CGAL::sqrt,      Lazy_exact_Sqrt)

// A macro for +, -, * and /
#define CGAL_LAZY_BINARY_OP(OP, NAME)                                 \
template <typename AK, typename EK, typename E2A>                                                \
struct NAME : public Lazy_exact_binary<AK,EK,E2A>                            \
{                                                                     \
  typedef Lazy_exact_rep<typename EK::FT,E2A> BaseBase;               \
  typedef Lazy_exact_nt<AK,EK,E2A> L;          \
  NAME (const L &a, const L &b)                                       \
    : Lazy_exact_binary<AK,EK,E2A>(a.approx() OP b.approx(), a, b) {}        \
                                                                      \
                                                                     \
  void                                                                \
  print(std::ostream& os, int level) const                           \
  {                                                                  \
    print_at_et(os, level);                                      \
    if(is_lazy()){                                                   \
     CGAL::msg(os, level, "Two child node:");                 \
      op1.print(os, level+1);                                 \
      op2.print(os, level+1);                                 \
    }                                                                \
  }                                                                   \
  void update_exact()                                                 \
  {                                                                   \
    this->et = new typename EK::FT(this->op1.exact() OP this->op2.exact());        \
  }                                                                   \
};

CGAL_LAZY_BINARY_OP(+, Lazy_exact_Add)
CGAL_LAZY_BINARY_OP(-, Lazy_exact_Sub)
CGAL_LAZY_BINARY_OP(*, Lazy_exact_Mul)
CGAL_LAZY_BINARY_OP(/, Lazy_exact_Div)

// Minimum
template <typename AK, typename EK, typename E2A>
struct Lazy_exact_Min : public Lazy_exact_binary<AK,EK,E2A>
{
  typedef Lazy_exact_nt<AK,EK,E2A> L;
  Lazy_exact_Min (const L &a, const L &b)
    : Lazy_exact_binary<AK,EK,E2A>(min(a.approx(), b.approx()), a, b) {}

  void update_exact()
  {
    this->et = new typename EK::FT(min(this->op1.exact(), this->op2.exact()));
  }
};

// Maximum
template <typename AK, typename EK, typename E2A>
struct Lazy_exact_Max : public Lazy_exact_binary<AK,EK,E2A>
{
  typedef Lazy_exact_nt<AK,EK,E2A> L;
  Lazy_exact_Max (const L &a, const L &b)
    : Lazy_exact_binary<AK,EK,E2A>(max(a.approx(), b.approx()), a, b) {}

  void update_exact()
  {
    this->et = new typename EK::FT(max(this->op1.exact(), this->op2.exact()));

  }
};



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
    at = e2a(e); 
  }

  void
  print(std::ostream& os, int level) const
  {
    print_at_et(os, level);
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
    print_at_et(os, level);
    if(is_lazy()){
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
    print_at_et(os, level);
    if(is_lazy()){
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
    print_at_et(os, level);
    if(is_lazy()){
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
    print_at_et(os, level);
    
    if(is_lazy()){
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
    print_at_et(os, level);
    
    if(is_lazy()){
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

///_______________________________________________________________

template <typename AK, typename EK, typename E2A>
class Lazy_exact_nt : public Handle
{
public :
  typedef typename AK::FT AT;
  typedef typename EK::FT ET;
  typedef typename Number_type_traits<ET>::Has_gcd      Has_gcd;
  typedef typename Number_type_traits<ET>::Has_division Has_division;
  typedef typename Number_type_traits<ET>::Has_sqrt     Has_sqrt;

  typedef Lazy_exact_nt<AK,EK,E2A> Self;
  typedef Lazy_construct_rep<AT, ET, E2A> Self_rep;

  // Lazy_exact_nt () {} // Handle is not such a nice stuff...  at the moment.

  Lazy_exact_nt (Self_rep *r)
  { PTR = r; }

  // Operations
  Lazy_exact_nt (double d)
  { PTR = new Lazy_exact_Cst<ET,E2A>(d); }

  Lazy_exact_nt (int i = 0)
  { PTR = new Lazy_exact_Int_Cst<ET,E2A>(i); }

  Lazy_exact_nt (const ET & e)
  { PTR = new Lazy_exact_Ex_Cst<ET,E2A>(e); }

  /*
  template <class ET1>
  Lazy_exact_nt (const Lazy_exact_nt<ET1> &x)
  { PTR = new Lazy_lazy_exact_Cst<ET, ET1>(x); }
  */
  Self operator- () const
  { return new Lazy_exact_Opp<AK,EK,E2A>(*this); }

  const Interval_nt<true>& approx() const
  { return ptr()->approx(); }

  Interval_nt<false> interval() const
  { 
    const Interval_nt<true>& i = ptr()->approx();
    return Interval_nt<false>(i.inf(), i.sup());
  }

  Interval_nt_advanced approx_adv() const
  { return ptr()->approx(); }

  const ET & exact() const
  { return ptr()->exact(); }

  static const double & get_relative_precision_of_to_double()
  {
      return relative_precision_of_to_double;
  }

  static void set_relative_precision_of_to_double(const double & d)
  {
      CGAL_assertion(d > 0 && d < 1);
      relative_precision_of_to_double = d;
  }

  void
  print(std::ostream& os, int level) const
  {
    ptr()->print(os, level);
  }

private:
  Self_rep * ptr() const { return (Self_rep*) PTR; }

  static double relative_precision_of_to_double;
};

template <typename AK, typename EK, typename E2A >
double Lazy_exact_nt<AK,EK,E2A>::relative_precision_of_to_double = 0.00001;


template <typename AK, typename EK, typename E2A>
bool
operator<(const Lazy_exact_nt<AK,EK,E2A>& a, const Lazy_exact_nt<AK,EK,E2A>& b)
{
  try
  {
    return a.approx() < b.approx();
  }
  catch (Interval_nt<false>::unsafe_comparison)
  {
    // std::cerr << "Interval filter failure (<)" << std::endl;
    return a.exact() < b.exact();
  }
}

template <typename AK, typename EK, typename E2A>
bool
operator==(const Lazy_exact_nt<AK,EK,E2A>& a, const Lazy_exact_nt<AK,EK,E2A>& b)
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

template <typename AK, typename EK, typename E2A>
inline
bool
operator>(const Lazy_exact_nt<AK,EK,E2A>& a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return b < a; }

template <typename AK, typename EK, typename E2A>
inline
bool
operator<=(const Lazy_exact_nt<AK,EK,E2A>& a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return ! (b < a); }

template <typename AK, typename EK, typename E2A>
inline
bool
operator>=(const Lazy_exact_nt<AK,EK,E2A>& a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return ! (a < b); }

template <typename AK, typename EK, typename E2A>
inline
bool
operator!=(const Lazy_exact_nt<AK,EK,E2A>& a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return ! (a == b); }


// Mixed operators with int.
template <typename AK, typename EK, typename E2A>
bool
operator<(int a, const Lazy_exact_nt<AK,EK,E2A>& b)
{
  try {
    return a < b.approx();
  }
  catch (Interval_nt<false>::unsafe_comparison) {
    return a < b.exact();
  }
}

template <typename AK, typename EK, typename E2A>
bool
operator<(const Lazy_exact_nt<AK,EK,E2A>& a, int b)
{
  try {
    return a.approx() < b;
  }
  catch (Interval_nt<false>::unsafe_comparison) {
    return a.exact() < b;
  }
}

template <typename AK, typename EK, typename E2A>
bool
operator==(int a, const Lazy_exact_nt<AK,EK,E2A>& b)
{
  try {
    return a == b.approx();
  }
  catch (Interval_nt<false>::unsafe_comparison) {
    return a == b.exact();
  }
}

template <typename AK, typename EK, typename E2A>
inline
bool
operator==(const Lazy_exact_nt<AK,EK,E2A>& a, int b)
{ return b == a; }

template <typename AK, typename EK, typename E2A>
inline
bool
operator>(int a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return b < a; }

template <typename AK, typename EK, typename E2A>
inline
bool
operator>(const Lazy_exact_nt<AK,EK,E2A>& a, int b)
{ return b < a; }

template <typename AK, typename EK, typename E2A>
inline
bool
operator<=(int a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return ! (b < a); }

template <typename AK, typename EK, typename E2A>
inline
bool
operator<=(const Lazy_exact_nt<AK,EK,E2A>& a, int b)
{ return ! (b < a); }

template <typename AK, typename EK, typename E2A>
inline
bool
operator>=(int a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return ! (a < b); }

template <typename AK, typename EK, typename E2A>
inline
bool
operator>=(const Lazy_exact_nt<AK,EK,E2A>& a, int b)
{ return ! (a < b); }

template <typename AK, typename EK, typename E2A>
inline
bool
operator!=(int a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return ! (a == b); }

template <typename AK, typename EK, typename E2A>
inline
bool
operator!=(const Lazy_exact_nt<AK,EK,E2A>& a, int b)
{ return ! (b == a); }


template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator+(const Lazy_exact_nt<AK,EK,E2A>& a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return new Lazy_exact_Add<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator-(const Lazy_exact_nt<AK,EK,E2A>& a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return new Lazy_exact_Sub<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator*(const Lazy_exact_nt<AK,EK,E2A>& a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return new Lazy_exact_Mul<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator/(const Lazy_exact_nt<AK,EK,E2A>& a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return new Lazy_exact_Div<AK,EK,E2A>(a, b); }

// mixed operators
template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator+(const Lazy_exact_nt<AK,EK,E2A>& a, int b)
{ return new Lazy_exact_Add<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator-(const Lazy_exact_nt<AK,EK,E2A>& a, int b)
{ return new Lazy_exact_Sub<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator*(const Lazy_exact_nt<AK,EK,E2A>& a, int b)
{ return new Lazy_exact_Mul<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator/(const Lazy_exact_nt<AK,EK,E2A>& a, int b)
{ return new Lazy_exact_Div<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator+(int a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return new Lazy_exact_Add<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator-(int a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return new Lazy_exact_Sub<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator*(int a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return new Lazy_exact_Mul<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
Lazy_exact_nt<AK,EK,E2A>
operator/(int a, const Lazy_exact_nt<AK,EK,E2A>& b)
{ return new Lazy_exact_Div<AK,EK,E2A>(a, b); }





template <typename AK, typename EK, typename E2A>
double
to_double(const Lazy_exact_nt<AK,EK,E2A> & a)
{
    const Interval_nt<true>& app = a.approx();
    if (app.sup() == app.inf())
	return app.sup();

    // If it's precise enough, then OK.
    if ((app.sup() - app.inf())
	    < Lazy_exact_nt<AK,EK,E2A>::get_relative_precision_of_to_double()
	      * std::max(std::fabs(app.inf()), std::fabs(app.sup())) )
        return CGAL::to_double(app);

    // Otherwise we trigger exact computation first,
    // which will refine the approximation.
    a.exact();
    return CGAL::to_double(a.approx());
}


template <typename AK, typename EK, typename E2A>
inline
std::pair<double,double>
to_interval(const Lazy_exact_nt<AK,EK,E2A> & a)
{
    return a.approx().pair();
}

template <typename AK, typename EK, typename E2A>
inline
Sign
sign(const Lazy_exact_nt<AK,EK,E2A> & a)
{
  try
  {
    return CGAL_NTS sign(a.approx());
  }
  catch (Interval_nt<false>::unsafe_comparison)
  {
    // std::cerr << "Interval filter failure (sign)" << std::endl;
    return CGAL_NTS sign(a.exact());
  }
}

template <typename AK, typename EK, typename E2A>
inline
Comparison_result
compare(const Lazy_exact_nt<AK,EK,E2A> & a, const Lazy_exact_nt<AK,EK,E2A> & b)
{
  try
  {
    return CGAL_NTS compare(a.approx(), b.approx());
  }
  catch (Interval_nt<false>::unsafe_comparison)
  {
    // std::cerr << "Interval filter failure (compare)" << std::endl;
    return CGAL_NTS compare(a.exact(), b.exact());
  }
}

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A>
abs(const Lazy_exact_nt<AK,EK,E2A> & a)
{ return new Lazy_exact_Abs<AK,EK,E2A>(a); }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A>
square(const Lazy_exact_nt<AK,EK,E2A> & a)
{ return new Lazy_exact_Square<AK,EK,E2A>(a); }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A>
sqrt(const Lazy_exact_nt<AK,EK,E2A> & a)
{ return new Lazy_exact_Sqrt<AK,EK,E2A>(a); }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A>
min(const Lazy_exact_nt<AK,EK,E2A> & a, const Lazy_exact_nt<AK,EK,E2A> & b)
{ return new Lazy_exact_Min<AK,EK,E2A>(a, b); }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A>
max(const Lazy_exact_nt<AK,EK,E2A> & a, const Lazy_exact_nt<AK,EK,E2A> & b)
{ return new Lazy_exact_Max<AK,EK,E2A>(a, b); }


template <typename AK, typename EK, typename E2A>
std::ostream &
operator<< (std::ostream & os, const Lazy_exact_nt<AK,EK,E2A> & a)
{ 
  if(is_pretty(os)){
    a.print(os, 0);
  }   else {
    os << CGAL::to_double(a); 
  }
  return os;
}





template <typename AK, typename EK, typename E2A>
std::istream &
operator>> (std::istream & is, Lazy_exact_nt<AK,EK,E2A> & a)
{
  typename EK::FT e;
  is >> e;
  a = e;
  return is;
}

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A> &
operator+=(Lazy_exact_nt<AK,EK,E2A> & a, const Lazy_exact_nt<AK,EK,E2A> & b)
{ return a = a + b; }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A> &
operator-=(Lazy_exact_nt<AK,EK,E2A> & a, const Lazy_exact_nt<AK,EK,E2A> & b)
{ return a = a - b; }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A> &
operator*=(Lazy_exact_nt<AK,EK,E2A> & a, const Lazy_exact_nt<AK,EK,E2A> & b)
{ return a = a * b; }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A> &
operator/=(Lazy_exact_nt<AK,EK,E2A> & a, const Lazy_exact_nt<AK,EK,E2A> & b)
{ return a = a / b; }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A> &
operator+=(Lazy_exact_nt<AK,EK,E2A> & a, int b)
{ return a = a + b; }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A> &
operator-=(Lazy_exact_nt<AK,EK,E2A> & a, int b)
{ return a = a - b; }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A> &
operator*=(Lazy_exact_nt<AK,EK,E2A> & a, int b)
{ return a = a * b; }

template <typename AK, typename EK, typename E2A>
inline
Lazy_exact_nt<AK,EK,E2A> &
operator/=(Lazy_exact_nt<AK,EK,E2A> & a, int b)
{ return a = a / b; }

template <typename AK, typename EK, typename E2A>
inline
bool
is_finite(const Lazy_exact_nt<AK,EK,E2A> & a)
{
  return is_finite(a.approx()) || is_finite(a.exact());
}

template <typename AK, typename EK, typename E2A>
inline
bool
is_valid(const Lazy_exact_nt<AK,EK,E2A> & a)
{
  return is_valid(a.approx()) || is_valid(a.exact());
}

template <typename AK, typename EK, typename E2A>
inline
io_Operator
io_tag (const Lazy_exact_nt<AK,EK,E2A>&)
{ return io_Operator(); }

/* TODO  
template <typename ET, typename AK, typename EK, typename E2A>
struct converter<ET, Lazy_exact_nt<AK,EK,E2A> >
{
    static inline typename ET do_it (const Lazy_exact_nt<AK,EK,E2A> & z)
    {
        return z.exact();
    }
};
*/






//____________________________________________________________
// A helper class to select the return type.

// AF: Regarder le tutorial metaprogramming

template <typename AK, typename EK, typename AT, typename ET, typename EFT, typename E2A>
struct Lazy_construction_return_type {
  typedef Lazy<AT, ET, EFT, E2A> result_type;
};

template <typename AK, typename EK, typename ET, typename EFT, typename E2A>
struct Lazy_construction_return_type<AK,EK,Interval_nt<>, ET,EFT,E2A>  {
  typedef Lazy_exact_nt<AK,EK,E2A>  result_type;
};


// Bboxes remain Bboxes and do not become lazy.
template <typename AK, typename EK, typename ET, typename EFT, typename E2A>
struct Lazy_construction_return_type<AK,EK,Bbox_2, ET,EFT,E2A>  {
  typedef Bbox_2  result_type;
};


// As people will write bool b = lazy_assign(..) 
// Or are there cases where we need a lazy bool in the kernel???
template <typename AK, typename EK, typename ET, typename EFT, typename E2A>
struct Lazy_construction_return_type<AK,EK,bool, ET,EFT,E2A>  {
  typedef bool  result_type;
};


template <typename AK, typename EK, typename AC, typename EC, typename EFT, typename E2A >
struct Lazy_construction_bbox_2 {
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





//____________________________________________________________
// The magic functor

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
      return new Lazy_construct_rep_1<AC, EC, E2A, L1>(ac, ec, l1);
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




