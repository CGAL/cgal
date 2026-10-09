// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_SMDS_3_UTILITIES_H
#define CGAL_SMDS_3_UTILITIES_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/functional.h>
#include <CGAL/Has_timestamp.h>
#include <CGAL/tags.h>

#include <iterator>
#include <string>
#include <sstream>
#include <type_traits>
#include <utility>

namespace CGAL {
namespace SMDS_3 {
namespace internal {

struct Debug_messages_tools {
  template <typename Vertex_handle>
  static std::string disp_vert(Vertex_handle v, Tag_true) {
    std::stringstream ss;
    ss.precision(17);
    ss << (void*)(&*v) << "[ts=" << v->time_stamp() << "]"
       << "(" << v->point() <<")";
    return ss.str();
  }

  template <typename Vertex_handle>
  static std::string disp_vert(Vertex_handle v, Tag_false) {
    std::stringstream ss;
    ss.precision(17);
    ss << (void*)(&*v) << "(" << v->point() <<")";
    return ss.str();
  }

  template <typename Vertex_handle>
  static std::string disp_vert(Vertex_handle v)
  {
    typedef typename std::iterator_traits<Vertex_handle>::value_type Vertex;
    return disp_vert(v, CGAL::internal::Has_timestamp<Vertex>());
  }
};

/**
 * @class First_of
 * Function object which returns the first element of a pair
 */
template <typename Pair>
struct First_of :
  public CGAL::cpp98::unary_function<Pair, const typename Pair::first_type&>
{
  typedef CGAL::cpp98::unary_function<Pair, const typename Pair::first_type&> Base;
  typedef typename Base::result_type                                  result_type;
  typedef typename Base::argument_type                                argument_type;

  result_type operator()(const argument_type& p) const { return p.first; }
}; // end class First_of


/**
 * @class Ordered_pair
 * Stores two elements in an ordered manner, i.e. first() < second()
 */
template <typename T>
class Ordered_pair
{
public:
  Ordered_pair(const T& t1, const T& t2)
  : data_(t1,t2)
  {
    if ( ! (t1 < t2) )
    {
      data_.second = t1;
      data_.first = t2;
    }
  }

  const T& first() const { return data_.first; }
  const T& second() const { return data_.second; }

  bool operator<(const Ordered_pair& rhs) const { return data_ < rhs.data_; }

private:
  std::pair<T,T> data_;
};


/**
 * @class Iterator_not_in_complex
 * @brief A class to filter elements which do not belong to the complex
 */
template < typename C3T3 >
class Iterator_not_in_complex
{
  const C3T3& c3t3_;
public:
  Iterator_not_in_complex(const C3T3& c3t3) : c3t3_(c3t3) { }

  template <typename Iterator>
  bool operator()(Iterator it) const { return ! c3t3_.is_in_complex(*it); }
}; // end class Iterator_not_in_complex


} // end namespace internal
} // end namespace SMDS_3

namespace SMDS_3_internal {
template <typename T, typename = void>
struct Has_in_dimension : std::false_type
{};

template <typename T>
struct Has_in_dimension<T, std::void_t<decltype(std::declval<T>().in_dimension())>>
  : std::true_type
{};

template <typename T, typename = void>
struct Has_is_corner : std::false_type
{};

template <typename T>
struct Has_is_corner<T, std::void_t<decltype(std::declval<T>().is_corner())>>
  : std::true_type
{};

template <typename Tr>
bool is_corner(const typename Tr::Vertex_handle v, const Tr&)
{
  using V = typename Tr::Triangulation_data_structure::Vertex;

  if constexpr(Has_in_dimension<V>::value)
    return v->in_dimension() == 0;
  else if constexpr(Has_is_corner<V>::value)
    return v->ccdt_3_data().is_corner();
  else
    return false;
}

} // namespace SMDS_3_internal
} //namespace CGAL

#endif // CGAL_SMDS_3_UTILITIES_H
