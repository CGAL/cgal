// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/certified_numeric_predicates.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_CERTIFIED_NUMERIC_PREDICATES_H
#define CGAL_CERTIFIED_NUMERIC_PREDICATES_H

#include <CGAL/number_utils.h>
#include <CGAL/Interval_arithmetic.h>

#include <boost/optional.hpp>
#include <boost/none.hpp>

CGAL_BEGIN_NAMESPACE

//
// NOTE: There is a "Certified" namespace in "Interval_arithmetic.h" which contains
// certified versions of the interval_nt<> operators and some predicates.
// We use a slightly different name here "certified" to avoid coallisions.
//
namespace certified {

// This should really be in optional itself
template<class T> optional<T> make_optional( T v ) { return optional<T>(v); }

template<class T> optional<T> make_optional( std::pair<T,bool> v ) 
{
  if ( v.second )
       return optional<T>(v); 
  else return boost::none ;  
}


//
// NOTE: CGAL::Certified defines certified comparison operators for interval_nt<>.
// In order to allow the certified predicate versions in this header to simply forward 
// to the non-certified counterparts, a using directive is used to bring about the certified
// comparison operators for the interval type.
// This allows the code to not need a specialization for interval_nt<>, except for the
// cases of sign() and compare().
// 
template <class NT>
inline
optional<bool>
is_zero(const NT& x)
{ 
  using namespace CGAL::Certified ;
  return make_optional( CGAL::is_zero(x) ; 
}

template <class NT>
inline
optional<bool>
is_one(const NT& x)
{
  using namespace CGAL::Certified ;
  return make_optional( CGAL::is_one(x) ; 
}

template <class NT>
inline
optional<bool>
is_negative(const NT& x)
{ 
  using namespace CGAL::Certified ;
  return make_optional( CGAL::is_negative(x) ; 
}

template <class NT>
inline
optional<bool>
is_positive(const NT& x)
{ 
  using namespace CGAL::Certified ;
  return make_optional( CGAL::is_positive(x) ; 
}

template <class NT>
inline
optional<Sign>
sign(const NT& x)
{ return make_optional(CGAL::sign(x)); }

template <class NT1, class NT2>
inline
optional<Comparison_result>
compare(const NT1& n1, const NT2& n2)
{ return make_optional(CGAL::compare(n1,n2)); }


//
// Interval_nt<> specializations.
//
template<bool Protected>
inline
optional<Sign>
sign(const Interval_nt<Protected>& x)
{ return make_optional(CGAL::Certified::sign(x)); }

template<bool Protected>
inline
optional<Comparison_result>
compare(const Interval_nt<Protected>& n1, const Interval_nt<Protected>& n2 )
{ return make_optional(CGAL::Certified::compare(n1,n2)); }

#define CGAL_CERTIFIED_NTS CGAL_NTS :: certified ::

CGAL_END_NAMESPACE

#endif // CGAL_CERTIFIED_NUMERIC_PREDICATES_H
 
