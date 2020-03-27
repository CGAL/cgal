// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno

#ifndef GEOMETRY_CONTAINER_H
#define GEOMETRY_CONTAINER_H
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <boost/geometry/io/wkt/write.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/shared_ptr.hpp>

struct Dummy_deleter{
  template<class T>
  void operator()(T*){
  }
};

namespace CGAL{
namespace internal{

/* \brief The Geometry_container struct is a wrapper
 that provides an easy way to use the read and write functions
 of `Polygon_IO`.

 \tparam Range is a model of `RandomAccessRange`
 \tparam TAG is the tag corresponding to the wkt structure you want
 to be read or written with this `CGAL::Geometry_container`
*/
template <typename Range, typename TAG>
struct Geometry_container{
  typedef std::pair<Range, TAG> type;
  typedef typename Range::difference_type difference_type;
  typedef typename Range::iterator iterator;
  typedef typename Range::const_iterator const_iterator;
  typedef typename Range::reverse_iterator reverse_iterator;
  typedef typename Range::const_reverse_iterator const_reverse_iterator;
  typedef typename Range::size_type size_type;
  typedef typename Range::value_type value_type;
  boost::shared_ptr<Range>  range;
  bool must_delete;
  //
  // Default constructor.
  // Creates a new internal Range.
  // De-allocate memory after usage.
  Geometry_container():range(new Range()), must_delete(true)
  {
  }
  /*
   Copy constructor.
   Memory NOT de-allocated after usage.
  */
  Geometry_container(Range& range)
    :range(&range, Dummy_deleter()), must_delete(false){}

  iterator begin()
  { return range->begin(); }

  iterator end()
  { return range->end(); }

  reverse_iterator rbegin()
  { return range->rbegin(); }

  reverse_iterator rend()
  { return range->rend(); }

  const_iterator begin()const
  { return range->begin(); }

  const_iterator end()const
  { return range->end(); }

  const_reverse_iterator rbegin()const
  { return range->rbegin(); }

  const_reverse_iterator rend()const
  { return range->rend(); }

  void clear(){ range->clear(); }

  template<typename size_t>
  void resize(size_t n){ range->resize(n); }

  template<typename T>
  void push_back(const T& t){ range->push_back(t); }

  size_type size() { return range->size(); }

  value_type operator[](size_type i)
  { return range[i]; }
};//end Geometry_container
} //end internal
}//end CGAL

namespace boost{

template< class T, typename TAG >
struct range_iterator<CGAL::internal::Geometry_container<T, TAG> >
{ typedef typename T::iterator type; };

template< class T, typename TAG >
struct range_iterator<const CGAL::internal::Geometry_container<T, TAG> >
{ typedef typename T::const_iterator type; };

template< class T, typename TAG >
struct range_mutable_iterator<CGAL::internal::Geometry_container<T, TAG> >
{ typedef typename range_mutable_iterator<T>::type type; };

}//end boost
#endif // GEOMETRY_CONTAINER_H
#endif
