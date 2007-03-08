// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sebastian Limbach  <slimbach@mpi-sb.mpg.de>

#ifndef CGAL_UTILS_CLASSES_H
#define CGAL_UTILS_CLASSES_H
#include <CGAL/number_type_basic.h>

CGAL_BEGIN_NAMESPACE

template < class NT1, class NT2 > struct Equal_to;
template < class NT1, class NT2 > struct Not_equal_to;
template < class NT1, class NT2 > struct Greater;
template < class NT1, class NT2 > struct Less;
template < class NT1, class NT2 > struct Greater_equal;
template < class NT1, class NT2 > struct Less_equal;


template < class NT, class Compare = std::less< NT > >
struct Min :public Binary_function< NT, NT, NT > {
 Min() {}
 Min(const Compare& c_) : c(c_) {}
 NT operator()( const NT& x, const NT& y) const
    { return (std::min) BOOST_PREVENT_MACRO_SUBSTITUTION ( x, y, c); }
protected:
 Compare c;
};

template < class NT, class Compare = std::less< NT > >
struct Max :public Binary_function< NT, NT, NT > {
 Max() {}
 Max(const Compare& c_) : c(c_) {}
 NT operator()( const NT& x, const NT& y) const
    { return (std::max) BOOST_PREVENT_MACRO_SUBSTITUTION ( x, y, c); }
protected:
 Compare c;
};

template< class Number_type >
class Is_valid 
  : public Unary_function< Number_type, bool > {
  public:
    bool operator()( const Number_type& ) const {
      return true;
    };
};




template < class NT1, class NT2 = NT1 >
struct Equal_to : public Binary_function< NT1, NT2, bool > {
  bool operator()( const NT1& x, const NT2& y) const
  { return x == y; }
};

template < class NT1, class NT2 = NT1 >
struct Not_equal_to : public Binary_function< NT1, NT2, bool > {
  bool operator()( const NT1& x, const NT2& y) const
  { return x != y; }
};

template < class NT1, class NT2 = NT1 >
struct Greater : public Binary_function< NT1, NT2, bool > {
  bool operator()( const NT1& x, const NT2& y) const
  { return x > y; }
};

template < class NT1, class NT2 = NT1 >
struct Less : public Binary_function< NT1, NT2, bool > {
  bool operator()( const NT1& x, const NT2& y) const
  { return x < y; }
};

template < class NT1, class NT2 = NT1 >
struct Greater_equal : public Binary_function< NT1, NT2, bool > {
  bool operator()( const NT1& x, const NT2& y) const
  { return x >= y; }
};

template < class NT1, class NT2 = NT1 >
struct Less_equal : public Binary_function< NT1, NT2, bool > {
  bool operator()( const NT1& x, const NT2& y) const
  { return x <= y; }
};

CGAL_END_NAMESPACE

#endif // CGAL_UTILS_CLASSES_H
