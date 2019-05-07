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
#include <cmath>

// Line discretization
#define CGAL_KSR_SAME_VECTOR_TOLERANCE 0.99999
#define CGAL_KSR_SAME_POINT_TOLERANCE 1e-10

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
using std::vector;
#else

typedef boost::uint32_t size_t;

template <typename ValueType>
class vector
{
public:
  typedef ValueType value_type;
  typedef std::vector<ValueType> Base;
  typedef typename Base::const_iterator const_iterator;
  typedef typename Base::iterator iterator;
private:
  std::vector<ValueType> m_data;
public:

  vector (KSR::size_t size = 0)
    : m_data (size) { }
  vector (KSR::size_t size, const ValueType& def)
    : m_data (size, def) { }

  const_iterator begin() const { return m_data.begin(); }
  const_iterator end() const { return m_data.end(); }
  iterator begin() { return m_data.begin(); }
  iterator end() { return m_data.end(); }

  KSR::size_t size() const { return static_cast<KSR::size_t>(m_data.size()); }
  bool empty() const { return m_data.empty(); }

  void reserve (const KSR::size_t& size) { m_data.reserve(std::size_t(size)); }
  void resize (const KSR::size_t& size) { m_data.resize(std::size_t(size)); }

  const ValueType& operator[] (const KSR::size_t& idx) const { return m_data[std::size_t(idx)]; }
  ValueType& operator[] (const KSR::size_t& idx) { return m_data[std::size_t(idx)]; }

  void erase (iterator it) { m_data.erase (it); }
  void insert (iterator it, const ValueType& v) { m_data.insert (it, v); }

  const ValueType& front() const { return m_data.front(); }
  ValueType& front() { return m_data.front(); }
  const ValueType& back() const { return m_data.back(); }
  ValueType& back() { return m_data.back(); }

  void push_back (const ValueType& v) { m_data.push_back (v); }

  void swap (vector& other) { m_data.swap (other.m_data); }
};

#endif

typedef vector<KSR::size_t> Idx_vector;
typedef typename Idx_vector::iterator Idx_iterator;

typedef std::set<KSR::size_t> Idx_set;
typedef typename Idx_set::iterator Idx_set_iterator;

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
  bool intersection_found = intersection_2 (t1, t2, out);
  CGAL_assertion_msg (intersection_found, "Intersection not found");
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
  bool intersection_found = intersection_3 (t1, t2, out);
  CGAL_assertion_msg (intersection_found, "Intersection not found");
  return out;
}

template <typename Point_3>
bool do_intersect (const vector<Point_3>& a, const vector<Point_3>& b)
{
  typedef typename Kernel_traits<Point_3>::Kernel::Triangle_3 Triangle_3;
  
  for (KSR::size_t i = 1; i < a.size() - 1; ++ i)
  {
    Triangle_3 ta (a[0], a[i], a[i+1]);
    for (KSR::size_t j = 1; j < b.size() - 1; ++ j)
    {
      Triangle_3 tb (b[0], b[j], b[j+1]);
      if (CGAL::do_intersect (ta, tb))
        return true;
    }
  }
  
  return false;
}

template <typename Line_3, typename Point_3>
inline bool intersection_3 (const Line_3& seg, const vector<Point_3>& polygon, Point_3& result)
{
  typedef typename Kernel_traits<Point_3>::Kernel::Triangle_3 Triangle_3;

  for (KSR::size_t i = 1; i < polygon.size() - 1; ++ i)
  {
    Triangle_3 triangle (polygon[0], polygon[i], polygon[i+1]);
    if (intersection_3 (seg, triangle, result))
      return true;
  }
  
  return false;
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
