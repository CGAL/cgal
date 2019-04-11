// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_UTILS_H
#define CGAL_KSR_UTILS_H

#include <sstream>

#define CGAL_KSR_ASSERT_POINTS_ALMOST_EQUAL(a,b) \
  CGAL_assertion_msg (CGAL::approximate_sqrt(CGAL::squared_distance((a), (b))) < 1e-15, \
                      std::string("Points " + CGAL::KSR::to_string(a) + " and " \
                                  + CGAL::KSR::to_string(b) + " should be almost equal").c_str())

namespace CGAL
{
namespace KSR
{

// Size type
#ifdef CGAL_KSR_USE_STD_SIZE_T_AS_SIZE_TYPE
typedef std::size_t size_t;
#else
typedef boost::uint32_t size_t;
#endif

// Use -1 as no element identifier
inline size_t no_element() { return size_t(-1); }

// Use -2 as special invalid identifier
inline size_t invalid() { return size_t(-2); }

// Vector normalization
template <typename Vector>
inline Vector normalize (const Vector& v)
{
  return v / CGAL::approximate_sqrt(v*v);
}

template <typename Type1, typename Type2, typename ResultType>
inline bool intersection_2 (const Type1& t1, const Type2& t2, ResultType& result)
{
  typedef typename Kernel_traits<Type1>::Kernel::Intersect_2 Intersect_2;
    
  typename cpp11::result_of<Intersect_2(Type1, Type2)>::type
    inter = intersection (t1, t2);
  if (!inter)
    return false;
        
  if (const ResultType* typed_inter = boost::get<ResultType>(&*inter))
  {
    result = *typed_inter;
    return true;
  }
  return false;
}

template <typename ResultType, typename Type1, typename Type2>
inline ResultType intersection_2 (const Type1& t1, const Type2& t2)
{
  ResultType out;
  bool okay = intersection_2 (t1, t2, out);
  CGAL_assertion_msg (okay, "Intersection not found");
  return out;
}
  
template <typename Type1, typename Type2, typename ResultType>
inline bool intersection_3 (const Type1& t1, const Type2& t2, ResultType& result)
{
  typedef typename Kernel_traits<Type1>::Kernel::Intersect_3 Intersect_3;
    
  typename cpp11::result_of<Intersect_3(Type1, Type2)>::type
    inter = intersection (t1, t2);
  if (!inter)
    return false;
        
  if (const ResultType* typed_inter = boost::get<ResultType>(&*inter))
  {
    result = *typed_inter;
    return true;
  }
  return false;
}

template <typename ResultType, typename Type1, typename Type2>
inline ResultType intersection_3 (const Type1& t1, const Type2& t2)
{
  ResultType out;
  bool okay = intersection_3 (t1, t2, out);
  CGAL_assertion_msg (okay, "Intersection not found");
  return out;
}

template <typename Point>
std::string to_string (const Point& p)
{
  std::ostringstream oss;
  oss.precision(18);
  oss << p;
  return oss.str();
}

}
} 


#endif // CGAL_KSR_UTILS_H
