// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL$
// $Id$
//
//
// Author(s)     : Andreas Fabri, Sylvain Pion

#ifndef CGAL_LAZY_H
#define CGAL_LAZY_H

#include <CGAL/basic.h>
#include <CGAL/Handle.h>
#include <CGAL/Object.h>
#include <CGAL/Kernel/Type_mapper.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/min_max_n.h>
#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <vector>
#include <CGAL/Default.h>

#ifdef CGAL_HAS_THREADS
#  include <boost/thread/tss.hpp>
#endif

namespace CGAL {

template <typename AT, typename ET, typename EFT, typename E2A> class Lazy;

template <typename ET_>
class Lazy_exact_nt;

template <typename AT, typename ET, typename EFT, typename E2A>
inline
const AT&
approx(const Lazy<AT,ET, EFT, E2A>& l)
{
  return l.approx();
}

// Where is this one (non-const) needed ?  Is it ?
template <typename AT, typename ET, typename EFT, typename E2A>
inline
AT&
approx(Lazy<AT,ET, EFT, E2A>& l)
{
  return l.approx();
}


template <typename AT, typename ET, typename EFT, typename E2A>
inline
const ET&
exact(const Lazy<AT,ET,EFT,E2A>& l)
{
  return l.exact();
}


template <typename AT, typename ET, typename EFT, typename E2A>
inline
unsigned
depth(const Lazy<AT,ET,EFT,E2A>& l)
{
  return l.depth();
}


#define CGAL_LAZY_FORWARD(T) \
  inline const T & approx(const T& d) { return d; } \
  inline const T & exact (const T& d) { return d; } \
  inline unsigned  depth (const T&  ) { return 0; }


CGAL_LAZY_FORWARD(double)
CGAL_LAZY_FORWARD(float)
CGAL_LAZY_FORWARD(int)
CGAL_LAZY_FORWARD(unsigned)
CGAL_LAZY_FORWARD(Return_base_tag)
CGAL_LAZY_FORWARD(Null_vector)
CGAL_LAZY_FORWARD(Origin)
CGAL_LAZY_FORWARD(Orientation)
CGAL_LAZY_FORWARD(Bbox_2)
CGAL_LAZY_FORWARD(Bbox_3)



#ifdef CGAL_LAZY_KERNEL_DEBUG
template <class T>
void
print_at(std::ostream& os, const T& at)
{
  os << at;
}

template <class T>
void
print_at(std::ostream& os, const std::vector<T>& at)
{
  os << "std::vector";
}

template <>
void
print_at(std::ostream& os, const Object& o)
{
  os << "Object";
}

template <class T1, class T2>
void
print_at(std::ostream& os, const std::pair<T1,T2> & at)
{
  os << "[ " << at.first << " | " << at.second << " ]" << std::endl ;
}


template <typename AT, typename ET, typename EFT, typename E2A>
inline
void
print_dag(const Lazy<AT,ET,EFT,E2A>& l, std::ostream& os, int level = 0)
{
  l.print_dag(os, level);
}

inline
void
print_dag(double d, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++)
    os << "    ";
  os << d << std::endl;
}

inline
void
msg(std::ostream& os, int level, char* s)
{
    for(int i = 0; i < level; i++)
      os << "    ";
    os << s << std::endl;
}

inline
void
print_dag(const Null_vector& nv, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++)
    os << "    ";
  os << "Null_vector" << std::endl;
}

inline
void
print_dag(const Origin& nv, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++)
    os << "    ";
  os << "Origin" << std::endl;
}
#endif


struct Depth_base {
#ifdef CGAL_PROFILE
  unsigned depth_;
  Depth_base() { set_depth(0); }
  unsigned depth() const { return depth_; }
  void set_depth(unsigned i)
  {
    depth_ = i;
    CGAL_HISTOGRAM_PROFILER(std::string("[Lazy_kernel DAG depths]"), i);
                            //(unsigned) ::log2(double(i)));
  }
#else
  unsigned depth() const { return 0; }
  void set_depth(unsigned) {}
#endif
};


// Abstract base class for lazy numbers and lazy objects
template <typename AT_, typename ET, typename E2A>
class Lazy_rep : public Rep, public Depth_base
{
  Lazy_rep (const Lazy_rep&); // cannot be copied.

public:

  typedef AT_ AT;

  mutable AT at;
  mutable ET *et;

  Lazy_rep ()
    : at(), et(NULL) {}

  Lazy_rep (const AT& a)
      : at(a), et(NULL) {}

  Lazy_rep (const AT& a, const ET& e)
      : at(a), et(new ET(e)) {}

  const AT& approx() const
  {
      return at;
  }

  AT& approx()
  {
      return at;
  }

  const ET & exact() const
  {
    if (et==NULL)
      update_exact();
    return *et;
  }

  ET & exact()
  {
    if (et==NULL)
      update_exact();
    return *et;
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void print_at_et(std::ostream& os, int level) const
  {
    for(int i = 0; i < level; i++){
      os << "    ";
    }
    os << "Approximation: ";
    print_at(os, at);
    os << std::endl;
    if(! is_lazy()){
      for(int i = 0; i < level; i++){
	os << "    ";
      }
      os << "Exact: ";
      print_at(os, *et);
      os << std::endl;
    }
  }

  virtual void print_dag(std::ostream& os, int level) const {}
#endif

  bool is_lazy() const { return et == NULL; }
  virtual void update_exact() const = 0;
  virtual ~Lazy_rep() { delete et; }
};


//____________________________________________________________
// The rep for the leaf node
// FIXME TODO : Factorize all the Lazy_rep_[0-8] !!!

template <typename AT, typename ET, typename E2A>
class Lazy_rep_0 : public Lazy_rep<AT, ET, E2A>
{

  typedef Lazy_rep<AT, ET, E2A> Base;
public:

  void
  update_exact() const
  {
    this->et = new ET();
  }

  Lazy_rep_0()
    : Lazy_rep<AT,ET, E2A>() {}

  Lazy_rep_0(const AT& a, const ET& e)
    : Lazy_rep<AT,ET,E2A>(a, e) {}

  Lazy_rep_0(const AT& a, void*)
    : Lazy_rep<AT,ET,E2A>(a) {}

  Lazy_rep_0(const ET& e)
    : Lazy_rep<AT,ET,E2A>(E2A()(e), e) {}

  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
  }
};


//____________________________________________________________

template <typename AC, typename EC, typename E2A, typename L1>
class Lazy_rep_1
  : public Lazy_rep<typename AC::result_type, typename EC::result_type, E2A>
  , private EC
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET(ec()(CGAL::exact(l1_)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
  }

  Lazy_rep_1(const AC& ac, const EC& ec, const L1& l1)
    : Lazy_rep<AT,ET, E2A>(ac(CGAL::approx(l1))), EC(ec), l1_(l1)
  {
    this->set_depth(CGAL::depth(l1_) + 1);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with one child node:");
      CGAL::print_dag(l1_, os, level+1);
    }
  }
#endif
};



//____________________________________________________________

template <typename AC, typename EC, typename E2A, typename L1, typename L2>
class Lazy_rep_2
  : public Lazy_rep<typename AC::result_type, typename EC::result_type, E2A>
  , private EC
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET(ec()(CGAL::exact(l1_), CGAL::exact(l2_)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
  }

  Lazy_rep_2(const AC& ac, const EC& /*ec*/, const L1& l1, const L2& l2)
    : Lazy_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2))),
      l1_(l1), l2_(l2)
  {
    this->set_depth(max_n(CGAL::depth(l1_), CGAL::depth(l2_)) + 1);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with two child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
    }
  }
#endif
};


//____________________________________________________________

template <typename AC, typename EC, typename E2A,
          typename L1, typename L2, typename L3>
class Lazy_rep_3
  : public Lazy_rep<typename AC::result_type, typename EC::result_type, E2A>
  , private EC
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;
  mutable L3 l3_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET(ec()(CGAL::exact(l1_), CGAL::exact(l2_),
                           CGAL::exact(l3_)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
    l3_ = L3();
  }

  Lazy_rep_3(const AC& ac, const EC& /*ec*/,
             const L1& l1, const L2& l2, const L3& l3)
    : Lazy_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2),
                             CGAL::approx(l3))),
      l1_(l1), l2_(l2), l3_(l3)
  {
    this->set_depth(max_n(CGAL::depth(l1_),
                          CGAL::depth(l2_),
                          CGAL::depth(l3_)) + 1);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with three child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
      CGAL::print_dag(l3_, os, level+1);
    }
  }
#endif
};


//____________________________________________________________

template <typename AC, typename EC, typename E2A,
          typename L1, typename L2, typename L3, typename L4>
class Lazy_rep_4
  : public Lazy_rep<typename AC::result_type, typename EC::result_type, E2A>
  , private EC
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;
  mutable L3 l3_;
  mutable L4 l4_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET(ec()(CGAL::exact(l1_), CGAL::exact(l2_),
                           CGAL::exact(l3_), CGAL::exact(l4_)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
    l3_ = L3();
    l4_ = L4();
  }

  Lazy_rep_4(const AC& ac, const EC& /*ec*/,
             const L1& l1, const L2& l2, const L3& l3, const L4& l4)
   : Lazy_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2),
                            CGAL::approx(l3), CGAL::approx(l4))),
     l1_(l1), l2_(l2), l3_(l3), l4_(l4)
  {
    this->set_depth(max_n(CGAL::depth(l1_),
                          CGAL::depth(l2_),
                          CGAL::depth(l3_),
                          CGAL::depth(l4_)) + 1);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);

    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with four child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
      CGAL::print_dag(l3_, os, level+1);
      CGAL::print_dag(l4_, os, level+1);
    }
  }
#endif
};

//____________________________________________________________


template <typename AC, typename EC, typename E2A,
          typename L1, typename L2, typename L3, typename L4, typename L5>
class Lazy_rep_5
  : public Lazy_rep<typename AC::result_type, typename EC::result_type, E2A>
  , private EC
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;
  mutable L3 l3_;
  mutable L4 l4_;
  mutable L5 l5_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET(ec()(CGAL::exact(l1_), CGAL::exact(l2_),
                           CGAL::exact(l3_), CGAL::exact(l4_),
                           CGAL::exact(l5_)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
    l3_ = L3();
    l4_ = L4();
    l5_ = L5();
  }

 Lazy_rep_5(const AC& ac, const EC& /*ec*/,
            const L1& l1, const L2& l2, const L3& l3, const L4& l4,
            const L5& l5)
    : Lazy_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2),
                             CGAL::approx(l3), CGAL::approx(l4),
                             CGAL::approx(l5))),
      l1_(l1), l2_(l2), l3_(l3), l4_(l4), l5_(l5)
  {
    this->set_depth(max_n(CGAL::depth(l1_),
                          CGAL::depth(l2_),
                          CGAL::depth(l3_),
                          CGAL::depth(l4_),
                          CGAL::depth(l5_)) + 1);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);

    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with five child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
      CGAL::print_dag(l3_, os, level+1);
      CGAL::print_dag(l4_, os, level+1);
      CGAL::print_dag(l5_, os, level+1);
    }
  }
#endif
};

template <typename AC, typename EC, typename E2A,
          typename L1, typename L2, typename L3, typename L4,
          typename L5, typename L6>
class Lazy_rep_6
  : public Lazy_rep<typename AC::result_type, typename EC::result_type, E2A>
  , private EC
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;
  mutable L3 l3_;
  mutable L4 l4_;
  mutable L5 l5_;
  mutable L6 l6_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET(ec()(CGAL::exact(l1_), CGAL::exact(l2_),
                           CGAL::exact(l3_), CGAL::exact(l4_),
                           CGAL::exact(l5_), CGAL::exact(l6_)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
    l3_ = L3();
    l4_ = L4();
    l5_ = L5();
    l6_ = L6();
  }

 Lazy_rep_6(const AC& ac, const EC& /*ec*/,
            const L1& l1, const L2& l2, const L3& l3, const L4& l4,
            const L5& l5, const L6& l6)
    : Lazy_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2),
                             CGAL::approx(l3), CGAL::approx(l4),
                             CGAL::approx(l5), CGAL::approx(l6))),
      l1_(l1), l2_(l2), l3_(l3), l4_(l4), l5_(l5), l6_(l6)
  {
    this->set_depth(max_n(CGAL::depth(l1_),
                          CGAL::depth(l2_),
                          CGAL::depth(l3_),
                          CGAL::depth(l4_),
                          CGAL::depth(l5_),
                          CGAL::depth(l6_)) + 1);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);

    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with 6 child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
      CGAL::print_dag(l3_, os, level+1);
      CGAL::print_dag(l4_, os, level+1);
      CGAL::print_dag(l5_, os, level+1);
      CGAL::print_dag(l6_, os, level+1);
    }
  }
#endif
};

template <typename AC, typename EC, typename E2A,
          typename L1, typename L2, typename L3, typename L4,
          typename L5, typename L6, typename L7>
class Lazy_rep_7
  : public Lazy_rep<typename AC::result_type, typename EC::result_type, E2A>
  , private EC
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;
  mutable L3 l3_;
  mutable L4 l4_;
  mutable L5 l5_;
  mutable L6 l6_;
  mutable L7 l7_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET(ec()(CGAL::exact(l1_), CGAL::exact(l2_),
                           CGAL::exact(l3_), CGAL::exact(l4_),
                           CGAL::exact(l5_), CGAL::exact(l6_),
                           CGAL::exact(l7_)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
    l3_ = L3();
    l4_ = L4();
    l5_ = L5();
    l6_ = L6();
    l7_ = L7();
  }

 Lazy_rep_7(const AC& ac, const EC& /*ec*/,
            const L1& l1, const L2& l2, const L3& l3, const L4& l4,
            const L5& l5, const L6& l6, const L7& l7)
    : Lazy_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2),
                             CGAL::approx(l3), CGAL::approx(l4),
                             CGAL::approx(l5), CGAL::approx(l6),
                             CGAL::approx(l7))),
      l1_(l1), l2_(l2), l3_(l3), l4_(l4), l5_(l5), l6_(l6), l7_(l7)
  {
    this->set_depth(max_n(CGAL::depth(l1_),
                          CGAL::depth(l2_),
                          CGAL::depth(l3_),
                          CGAL::depth(l4_),
                          CGAL::depth(l5_),
                          CGAL::depth(l6_),
                          CGAL::depth(l7_)) + 1);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);

    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with 7 child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
      CGAL::print_dag(l3_, os, level+1);
      CGAL::print_dag(l4_, os, level+1);
      CGAL::print_dag(l5_, os, level+1);
      CGAL::print_dag(l6_, os, level+1);
      CGAL::print_dag(l7_, os, level+1);
    }
  }
#endif
};

template <typename AC, typename EC, typename E2A,
          typename L1, typename L2, typename L3, typename L4,
          typename L5, typename L6, typename L7, typename L8>
class Lazy_rep_8
  : public Lazy_rep<typename AC::result_type, typename EC::result_type, E2A>
  , private EC
{
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;
  mutable L3 l3_;
  mutable L4 l4_;
  mutable L5 l5_;
  mutable L6 l6_;
  mutable L7 l7_;
  mutable L8 l8_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET(ec()(CGAL::exact(l1_), CGAL::exact(l2_),
                           CGAL::exact(l3_), CGAL::exact(l4_),
                           CGAL::exact(l5_), CGAL::exact(l6_),
                           CGAL::exact(l7_), CGAL::exact(l8_)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
    l3_ = L3();
    l4_ = L4();
    l5_ = L5();
    l6_ = L6();
    l7_ = L7();
    l8_ = L8();
  }

 Lazy_rep_8(const AC& ac, const EC& /*ec*/,
            const L1& l1, const L2& l2, const L3& l3, const L4& l4,
            const L5& l5, const L6& l6, const L7& l7, const L8& l8)
    : Lazy_rep<AT,ET,E2A>(ac(CGAL::approx(l1), CGAL::approx(l2),
                             CGAL::approx(l3), CGAL::approx(l4),
                             CGAL::approx(l5), CGAL::approx(l6),
                             CGAL::approx(l7), CGAL::approx(l8))),
       l1_(l1), l2_(l2), l3_(l3), l4_(l4), l5_(l5), l6_(l6), l7_(l7), l8_(l8)
  {
    this->set_depth(max_n(CGAL::depth(l1_),
                          CGAL::depth(l2_),
                          CGAL::depth(l3_),
                          CGAL::depth(l4_),
                          CGAL::depth(l5_),
                          CGAL::depth(l6_),
                          CGAL::depth(l7_),
                          CGAL::depth(l8_)) + 1);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);

    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with 8 child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
      CGAL::print_dag(l3_, os, level+1);
      CGAL::print_dag(l4_, os, level+1);
      CGAL::print_dag(l5_, os, level+1);
      CGAL::print_dag(l6_, os, level+1);
      CGAL::print_dag(l7_, os, level+1);
      CGAL::print_dag(l8_, os, level+1);
    }
  }
#endif
};


template < typename K1, typename K2 >
struct Approx_converter
{
  typedef K1         Source_kernel;
  typedef K2         Target_kernel;
  //typedef Converter  Number_type_converter;

  template < typename T >
  const typename T::AT&
  operator()(const T&t) const
  { return t.approx(); }

  const Null_vector&
  operator()(const Null_vector& n) const
  { return n; }
  
  const Bbox_2&
  operator()(const Bbox_2& b) const
  { return b; }
  
  const Bbox_3&
  operator()(const Bbox_3& b) const
  { return b; }
};

template < typename K1, typename K2 >
struct Exact_converter
{
  typedef K1         Source_kernel;
  typedef K2         Target_kernel;
  //typedef Converter  Number_type_converter;

  template < typename T >
  const typename T::ET&
  operator()(const T&t) const
  { return t.exact(); }

  const Null_vector&
  operator()(const Null_vector& n) const
  { return n; }
  
  const Bbox_2&
  operator()(const Bbox_2& b) const
  { return b; }
  
  const Bbox_3&
  operator()(const Bbox_3& b) const
  { return b; }
};

//____________________________________________________________



template <typename AC, typename EC, typename E2A, typename L1>
class Lazy_rep_with_vector_1
  : public Lazy_rep<std::vector<Object>, std::vector<Object>, E2A>
  , private EC
{
  typedef std::vector<Object> AT;
  typedef std::vector<Object> ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const 
  {
// TODO : This looks really unfinished...
    std::vector<Object> vec;
    this->et = new ET();
    //this->et->reserve(this->at.size());
    ec()(CGAL::exact(l1_), std::back_inserter(*(this->et)));
    if(this->et==NULL)
    E2A()(*(this->et));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
  }

  Lazy_rep_with_vector_1(const AC& ac, const EC& /*ec*/, const L1& l1)
    : l1_(l1)
  {
    ac(CGAL::approx(l1), std::back_inserter(this->at));
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    os << "A Lazy_rep_with_vector_1 of size " <<  this->at.size() << std::endl;
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with one child node:");
      CGAL::print_dag(l1_, os, level+1);

    }
  }
#endif
};


template <typename AC, typename EC, typename E2A, typename L1, typename L2>
class Lazy_rep_with_vector_2
  : public Lazy_rep<std::vector<Object>, std::vector<Object>, E2A>
  , private EC
{
  typedef std::vector<Object> AT;
  typedef std::vector<Object> ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET();
    this->et->reserve(this->at.size());
    ec()(CGAL::exact(l1_), CGAL::exact(l2_), std::back_inserter(*(this->et)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
  }

  Lazy_rep_with_vector_2(const AC& ac, const EC& /*ec*/, const L1& l1, const L2& l2)
    : l1_(l1), l2_(l2)
  {
    ac(CGAL::approx(l1), CGAL::approx(l2), std::back_inserter(this->at));
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    os << "A Lazy_rep_with_vector_2 of size " <<  this->at.size() << std::endl;
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with two child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
    }
  }
#endif
};


template <typename AC, typename EC, typename E2A, typename L1, typename L2, typename R1>
class Lazy_rep_2_1
  : public Lazy_rep<typename R1::AT, typename R1::ET, E2A>
  , private EC
{
  typedef typename R1::AT AT;
  typedef typename R1::ET ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET();
    ec()(CGAL::exact(l1_), CGAL::exact(l2_), *(this->et));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
  }

  Lazy_rep_2_1(const AC& ac, const EC& /*ec*/, const L1& l1, const L2& l2)
    : Lazy_rep<AT,ET,E2A>(), l1_(l1), l2_(l2)
  {
    ac(CGAL::approx(l1), CGAL::approx(l2), this->at);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    os << "A Lazy_rep_2_1" << std::endl;
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with two child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
    }
  }
#endif
};


//____________________________________________________________________________________
// The following rep class stores two non-const reference parameters of type R1 and R2

template <typename AC, typename EC, typename E2A, typename L1, typename L2, typename R1, typename R2>
class Lazy_rep_2_2
  : public Lazy_rep<std::pair<typename R1::AT,typename R2::AT>, std::pair<typename R1::ET, typename R2::ET>, E2A>
  , private EC
{
  typedef std::pair<typename R1::AT, typename R2::AT> AT;
  typedef std::pair<typename R1::ET, typename R2::ET> ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET();
    ec()(CGAL::exact(l1_), CGAL::exact(l2_), this->et->first, this->et->second );
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
  }

  Lazy_rep_2_2(const AC& ac, const EC& /*ec*/, const L1& l1, const L2& l2)
    : Lazy_rep<AT,ET,E2A>(), l1_(l1), l2_(l2)
  {
    ac(CGAL::approx(l1), CGAL::approx(l2), this->at.first, this->at.second);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    os << "A Lazy_rep_2_2"  << std::endl;
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with two child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
    }
  }
#endif
};


//____________________________________________________________
// The handle class
template <typename AT_, typename ET_, typename EFT, typename E2A>
class Lazy : public Handle
{
public :

  typedef Lazy<AT_, ET_, EFT, E2A>  Self;
  typedef Lazy_rep<AT_, ET_, E2A>   Self_rep;

  typedef AT_ AT; // undocumented
  typedef ET_ ET; // undocumented

  typedef AT  Approximate_type;
  typedef ET  Exact_type;

/*
  typedef Self Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }
*/

  Lazy()
    : Handle(zero()) {}

  // Before Lazy::zero() used Boost.Thread, the definition of Lazy() was:
  //   Lazy()
  //   #ifndef CGAL_HAS_THREAD
  //     : Handle(zero()) {}
  //   #else
  //   {
  //     PTR = new Lazy_rep_0<AT, ET, E2A>();
  //   }
  //   #endif

  Lazy(Self_rep *r)
  {
    PTR = r;
  }

  Lazy(const ET& e)
  {
    PTR = new Lazy_rep_0<AT,ET,E2A>(e);
  }

  const AT& approx() const
  { return ptr()->approx(); }

  const ET& exact() const
  { return ptr()->exact(); }

  AT& approx()
  { return ptr()->approx(); }

  ET& exact()
  { return ptr()->exact(); }

  unsigned depth() const
  {
    return ptr()->depth();
  }

  void print_dag(std::ostream& os, int level) const
  {
    ptr()->print_dag(os, level);
  }

private:

  // We have a static variable for optimizing the default constructor,
  // which is in particular heavily used for pruning DAGs.
  static const Self & zero()
  {
#ifdef CGAL_HAS_THREADS
    static boost::thread_specific_ptr<Self> z;
    if (z.get() == NULL) {
        z.reset(new Self(new Lazy_rep_0<AT, ET, E2A>()));
    }
    return * z.get();
#else
    static const Self z = new Lazy_rep_0<AT, ET, E2A>();
    return z;
#endif
  }

  Self_rep * ptr() const { return (Self_rep*) PTR; }
};





// The magic functor for Construct_bbox_[2,3], as there is no Lazy<Bbox>

template <typename LK, typename AC, typename EC>
struct Lazy_construction_bbox
{
  static const bool Protection = true;
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename AC::result_type result_type;

  AC ac;
  EC ec;

  template <typename L1>
  result_type operator()(const L1& l1) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
    Protect_FPU_rounding<Protection> P;
    try {
      return ac(CGAL::approx(l1));
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return ec(CGAL::exact(l1));
    }
  }
};

template <typename LK, typename AC, typename EC>
struct Lazy_construction_nt {

  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::E2A E2A;
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy_exact_nt<ET> result_type;

  AC ac;
  EC ec;

  template <typename L1>
  result_type operator()(const L1& l1) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return new Lazy_rep_1<AC, EC, To_interval<ET>, L1>(ac, ec, l1);
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return new Lazy_rep_0<AT,ET,To_interval<ET> >(ec(CGAL::exact(l1)));
    }
  }

  template <typename L1, typename L2>
  result_type operator()(const L1& l1, const L2& l2) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return new Lazy_rep_2<AC, EC, To_interval<ET>, L1,L2>(ac, ec, l1,l2);
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return new Lazy_rep_0<AT,ET,To_interval<ET> >(ec(CGAL::exact(l1), CGAL::exact(l2)));
    }
  }

  template <typename L1, typename L2, typename L3>
  result_type operator()(const L1& l1, const L2& l2, const L3& l3) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return new Lazy_rep_3<AC, EC, To_interval<ET>, L1,L2,L3>(ac, ec, l1,l2,l3);
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return new Lazy_rep_0<AT,ET,To_interval<ET> >(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3)));
    }
  }

  template <typename L1, typename L2, typename L3, typename L4>
  result_type operator()(const L1& l1, const L2& l2, const L3& l3, const L4& l4) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return new Lazy_rep_4<AC, EC, To_interval<ET>, L1,L2,L3,L4>(ac, ec, l1,l2,l3,l4);
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return new Lazy_rep_0<AT,ET,To_interval<ET> >(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3), CGAL::exact(l4)));
    }
  }

  template <typename L1, typename L2, typename L3, typename L4, typename L5>
  result_type operator()(const L1& l1, const L2& l2, const L3& l3, const L4& l4, const L5& l5) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return new Lazy_rep_5<AC, EC, To_interval<ET>, L1,L2,L3,L4,L5>(ac, ec, l1,l2,l3,l4,l5);
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return new Lazy_rep_0<AT,ET,To_interval<ET> >(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3), CGAL::exact(l4), CGAL::exact(l5)));
    }
  }
};


template <typename LK>
Object
make_lazy(const Object& eto)
{
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::E2A E2A;

  if (eto.is_empty())
    return Object();

#define CGAL_Kernel_obj(X) \
  if (const typename EK::X* ptr = object_cast<typename EK::X>(&eto)) \
    return make_object(typename LK::X(new Lazy_rep_0<typename AK::X, typename EK::X, E2A>(*ptr)));

#include <CGAL/Kernel/interface_macros.h>

  std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's (#2)" << std::endl;
  std::cerr << "dynamic type of the Object : " << eto.type().name() << std::endl;

  return Object();
}


// This functor selects the i'th element in a vector of Object's
// and casts it to what is in the Object

template <typename T2>
struct Ith {
  typedef T2 result_type;

  // We keep a Sign member object
  // for future utilisation, in case
  // we have pairs of 2 T2 objects e.g.
  // for a numeric_point vector returned
  // from a construction of a possible
  // lazy algebraic kernel

  int i;
  Sign sgn;

  Ith(int i_)
    : i(i_)
  {sgn=NEGATIVE;}

  Ith(int i_, bool b_)
    : i(i_)
  { sgn= (b_) ? POSITIVE : ZERO;}

  const T2&
  operator()(const std::vector<Object>& v) const
  {
    if(sgn==NEGATIVE)
    return *object_cast<T2>(&v[i]);

    typedef std::pair<T2,unsigned int >         Pair_type_1;
    typedef std::pair<T2,std::pair<bool,bool> > Pair_type_2;

    if(const Pair_type_1 *p1 = object_cast<Pair_type_1>(&v[i]))
    	return p1->first;
    else if(const Pair_type_2 *p2 = object_cast<Pair_type_2>(&v[i]))
        return p2->first;

    CGAL_error_msg( " Unexpected encapsulated type ");
  }
};


template <typename LK, typename AC, typename EC>
struct Lazy_cartesian_const_iterator_2
{
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::Cartesian_const_iterator_2 result_type;

  AC ac;
  EC ec;

public:

  template < typename L1>
  result_type
  operator()(const L1& l1) const
  {
    return result_type(&l1);
  }

  template < typename L1>
  result_type
  operator()(const L1& l1, int) const
  {
    return result_type(&l1,2);
  }

};


template <typename LK, typename AC, typename EC>
struct Lazy_cartesian_const_iterator_3
{
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::Cartesian_const_iterator_3 result_type;

  AC ac;
  EC ec;

public:

  template < typename L1>
  result_type
  operator()(const L1& l1) const
  {
    return result_type(&l1);
  }

  template < typename L1>
  result_type
  operator()(const L1& l1, int) const
  {
    return result_type(&l1,3);
  }

};


// This is the magic functor for functors that write their result in a  reference argument
// In a first version we assume that the references are of type Lazy<Something>,
// and that the result type is void

template <typename LK, typename AK, typename EK, typename AC, typename EC, typename EFT, typename E2A>
struct Lazy_functor_2_1
{
  static const bool Protection = true;
  typedef void result_type;

  AC ac;
  EC ec;

public:

  template <typename L1, typename L2, typename R1>
  void
  operator()(const L1& l1, const L2& l2, R1& r1) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      // we suppose that R1 is a Lazy<Something>
      r1 = R1(new Lazy_rep_2_1<AC, EC, E2A, L1, L2, R1>(ac, ec, l1, l2));
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      typename R1::ET et;
      ec(CGAL::exact(l1), CGAL::exact(l2), et);
      r1 = R1(new Lazy_rep_0<typename R1::AT,typename R1::ET,E2A>(et));
    }
  }
};


template <typename T>
struct First
{
   typedef typename T::first_type result_type;

   const typename T::first_type&
   operator()(const T& p) const
   {
     return p.first;
   }
 };

template <typename T>
struct Second
{
  typedef typename T::second_type result_type;

  const typename T::second_type&
  operator()(const T& p) const
  {
    return p.second;
  }
};

// This is the magic functor for functors that write their result in a reference argument
// In a first version we assume that the references are of type Lazy<Something>,
// and that the result type is void

//template <typename LK, typename AK, typename EK, typename AC, typename EC, typename EFT, typename E2A>
template <typename LK, typename AC, typename EC>
struct Lazy_functor_2_2
{
  static const bool Protection = true;

  typedef void result_type;
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;

  AC ac;
  EC ec;

public:

  template <typename L1, typename L2, typename R1, typename R2>
  void
  operator()(const L1& l1, const L2& l2, R1& r1, R2& r2) const
  {
    typedef Lazy<typename R1::AT, typename R1::ET, EFT, E2A> Handle_1;
    typedef Lazy<typename R2::AT, typename R2::ET, EFT, E2A> Handle_2;
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      typedef Lazy<std::pair<typename R1::AT, typename R2::AT>, std::pair<typename R1::ET, typename R2::ET>, EFT, E2A> Lazy_pair;
      Lazy_pair lv(new Lazy_rep_2_2<AC, EC, E2A, L1, L2, R1, R2>(ac, ec, l1, l2));
      // lv->approx() is a std::pair<R1::AT, R2::AT>;
      r1 = R1(Handle_1(new Lazy_rep_1<First<std::pair<typename R1::AT, typename R2::AT> >, First<std::pair<typename R1::ET, typename R2::ET> >, E2A, Lazy_pair>(First<std::pair<typename R1::AT, typename R2::AT> >(), First<std::pair<typename R1::ET, typename R2::ET> >(), lv)));
      r2 = R2(Handle_2(new Lazy_rep_1<Second<std::pair<typename R1::AT, typename R2::AT> >, Second<std::pair<typename R1::ET, typename R2::ET> >, E2A, Lazy_pair>(Second<std::pair<typename R1::AT, typename R2::AT> >(), Second<std::pair<typename R1::ET, typename R2::ET> >(), lv)));
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      typename R1::ET et1, et2;
      ec(CGAL::exact(l1), CGAL::exact(l2), et1, et2);
      r1 = R1(Handle_1(new Lazy_rep_0<typename R1::AT,typename R1::ET,E2A>(et1)));
      r2 = R2(Handle_2(new Lazy_rep_0<typename R2::AT,typename R2::ET,E2A>(et2)));
    }
  }
};


// This is the magic functor for functors that write their result as Objects into an output iterator

template <typename LK, typename AC, typename EC>
struct Lazy_intersect_with_iterators
{
  static const bool Protection = true;
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;
  typedef void result_type;
  typedef Lazy<Object, Object, EFT, E2A> Lazy_object;
  typedef Lazy<std::vector<Object>, std::vector<Object>, EFT, E2A> Lazy_vector;

  AC ac;
  EC ec;

public:

  // In the example we intersect two Lazy<Segment>s
  // and write into a back_inserter(list<Object([Lazy<Point>,Lazy<Segment>]) >)
  template <typename L1, typename L2, typename OutputIterator>
  OutputIterator
  operator()(const L1& l1, const L2& l2, OutputIterator it) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      Lazy_vector lv(new Lazy_rep_with_vector_2<AC, EC, E2A, L1, L2>(ac, ec, l1, l2));
      // lv.approx() is a std::vector<Object([AK::Point_2,AK::Segment_2])>
      // that is, when we get here we have constructed all approximate results
      for (unsigned int i = 0; i < lv.approx().size(); i++) {
// FIXME : I'm not sure how this work...
#define CGAL_Kernel_obj(X) if (object_cast<typename AK::X>(& (lv.approx()[i]))) { \
	  *it++ = make_object(typename LK::X(new Lazy_rep_1<Ith<typename AK::X>, \
                                                                      Ith<typename EK::X>, E2A, Lazy_vector> \
                                                 (Ith<typename AK::X>(i), Ith<typename EK::X>(i), lv))); \
          continue; \
	}

#include <CGAL/Kernel/interface_macros.h>

        std::cerr << "we need  more casts" << std::endl;
      }

    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      // TODO: Instead of using a vector, write an iterator adapter
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      std::vector<Object> exact_objects;
      ec(CGAL::exact(l1), CGAL::exact(l2), std::back_inserter(exact_objects));
      for (std::vector<Object>::const_iterator oit = exact_objects.begin();
	   oit != exact_objects.end();
	   ++oit){
	*it++ = make_lazy<LK>(*oit);
      }
    }
    return it;
  }
};


template <typename T>
struct Object_cast
{
  typedef T result_type;

  const T&
  operator()(const Object& o) const
  {
    return *object_cast<T>(&o);
  }
};

// The following functor returns an Object with a Lazy<Something> inside
// As the nested kernels return Objects of AK::Something and EK::Something
// we have to unwrap them from the Object, and wrap them in a Lazy<Something>
//
// TODO: write operators for other than two arguments. For the current kernel we only need two for Intersect_2

template <typename LK, typename AC, typename EC>
struct Lazy_construction_object
{
  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Object result_type;

  typedef Lazy<Object, Object, EFT, E2A> Lazy_object;

  AC ac;
  EC ec;

public:

  template <typename L1>
  result_type
  operator()(const L1& l1) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      Lazy_object lo(new Lazy_rep_1<AC, EC, E2A, L1>(ac, ec, l1));

      if(lo.approx().is_empty())
        return Object();

#define CGAL_Kernel_obj(X) \
      if (object_cast<typename AK::X>(& (lo.approx()))) { \
	typedef Lazy_rep_1<Object_cast<typename AK::X>, Object_cast<typename EK::X>, E2A, Lazy_object> Lcr; \
	Lcr * lcr = new Lcr(Object_cast<typename AK::X>(), Object_cast<typename EK::X>(), lo); \
	return make_object(typename LK::X(lcr)); \
      }

#include <CGAL/Kernel/interface_macros.h>

      std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's (#1)" << std::endl;
      std::cerr << "dynamic type of the Object : " << lo.approx().type().name() << std::endl;

    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      ET eto = ec(CGAL::exact(l1));
      return make_lazy<LK>(eto);
    }
    return Object();
  }

  template <typename L1, typename L2>
  result_type
  operator()(const L1& l1, const L2& l2) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      Lazy_object lo(new Lazy_rep_2<AC, EC, E2A, L1, L2>(ac, ec, l1, l2));

      if(lo.approx().is_empty())
        return Object();

#define CGAL_Kernel_obj(X) \
      if (object_cast<typename AK::X>(& (lo.approx()))) { \
	typedef Lazy_rep_1<Object_cast<typename AK::X>, Object_cast<typename EK::X>, E2A, Lazy_object> Lcr; \
	Lcr * lcr = new Lcr(Object_cast<typename AK::X>(), Object_cast<typename EK::X>(), lo); \
	return make_object(typename LK::X(lcr)); \
      }

#include <CGAL/Kernel/interface_macros.h>

      std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's (#1)" << std::endl;
      std::cerr << "dynamic type of the Object : " << lo.approx().type().name() << std::endl;

    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      ET eto = ec(CGAL::exact(l1), CGAL::exact(l2));
      return make_lazy<LK>(eto);
    }
    return Object();
  }

};



//____________________________________________________________
// The magic functor that has Lazy<Something> as result type

template <typename LK, typename AC, typename EC, typename E2A_ = Default>
struct Lazy_construction
{
  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename EK::FT EFT;
  typedef typename Default::Get<E2A_, typename LK::E2A>::type E2A;
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Lazy<AT, ET, EFT, E2A> Handle;
  typedef typename Type_mapper<AT,AK,LK>::type result_type;

  AC ac;
  EC ec;

public:

  result_type
  operator()() const
  {
    return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>()) );
  }

  template <typename L1>
  result_type
  operator()(const L1& l1) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return  result_type( Handle(new Lazy_rep_1<AC, EC, E2A, L1>(ac, ec, l1)) );
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1)))) );
    }
  }

  template <typename L1, typename L2>
  result_type
  operator()(const L1& l1, const L2& l2) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return result_type( Handle(new Lazy_rep_2<AC, EC, E2A, L1, L2>(ac, ec, l1, l2)) );
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2)))) );
    }
  }

  template <typename L1, typename L2, typename L3>
  result_type
  operator()(const L1& l1, const L2& l2, const L3& l3) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
      Protect_FPU_rounding<Protection> P;
    try {
      return result_type( Handle(new Lazy_rep_3<AC, EC, E2A, L1, L2, L3>(ac, ec, l1, l2, l3)) );
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3)))) );
    }
  }

  template <typename L1, typename L2, typename L3, typename L4>
  result_type
  operator()(const L1& l1, const L2& l2, const L3& l3, const L4& l4) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return result_type( Handle(new Lazy_rep_4<AC, EC, E2A, L1, L2, L3, L4>(ac, ec, l1, l2, l3, l4)) );
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3), CGAL::exact(l4)))) );
    }
  }

  template <typename L1, typename L2, typename L3, typename L4, typename L5>
  result_type
  operator()(const L1& l1, const L2& l2, const L3& l3, const L4& l4, const L5& l5) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return result_type( Handle(new Lazy_rep_5<AC, EC, E2A, L1, L2, L3, L4, L5>(ac, ec, l1, l2, l3, l4, l5)) );
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3), CGAL::exact(l4), CGAL::exact(l5)))) );
    }
  }

  template <typename L1, typename L2, typename L3, typename L4, typename L5, typename L6>
  result_type
  operator()(const L1& l1, const L2& l2, const L3& l3, const L4& l4, const L5& l5, const L6& l6) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return result_type( Handle(new Lazy_rep_6<AC, EC, E2A, L1, L2, L3, L4, L5, L6>(ac, ec, l1, l2, l3, l4, l5, l6)) );
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3), CGAL::exact(l4), CGAL::exact(l5), CGAL::exact(l6)))) );
    }
  }

  template <typename L1, typename L2, typename L3, typename L4, typename L5, typename L6, typename L7>
  result_type
  operator()(const L1& l1, const L2& l2, const L3& l3, const L4& l4, const L5& l5, const L6& l6, const L7& l7) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return result_type( Handle(new Lazy_rep_7<AC, EC, E2A, L1, L2, L3, L4, L5, L6, L7>(ac, ec, l1, l2, l3, l4, l5, l6, l7)) );
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3), CGAL::exact(l4), CGAL::exact(l5), CGAL::exact(l6), CGAL::exact(l7)))) );
    }
  }

  template <typename L1, typename L2, typename L3, typename L4, typename L5, typename L6, typename L7, typename L8>
  result_type
  operator()(const L1& l1, const L2& l2, const L3& l3, const L4& l4, const L5& l5, const L6& l6, const L7& l7, const L8& l8) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      return result_type( Handle(new Lazy_rep_8<AC, EC, E2A, L1, L2, L3, L4, L5, L6, L7, L8>(ac, ec, l1, l2, l3, l4, l5, l6, l7, l8)) );
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>(ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3), CGAL::exact(l4), CGAL::exact(l5), CGAL::exact(l6), CGAL::exact(l7), CGAL::exact(l8)))) );
    }
  }

};

} //namespace CGAL

#endif // CGAL_LAZY_H
