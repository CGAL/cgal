// Copyright (c) 2009 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_PROPERTY_MAP_H
#define CGAL_PROPERTY_MAP_H

#include <boost/property_map.hpp>
#include <boost/tuple/tuple.hpp>

namespace CGAL {

 template <typename T>
 struct identity_property_map
    : public boost::put_get_helper<T, 
        identity_property_map<T> >
  {
    typedef T key_type;
    typedef T value_type;
    typedef T reference;
    typedef boost::readable_property_map_tag category;

    inline value_type operator[](const key_type& v) const { return v; }
  };


 template <typename Pair>
 struct first_of_pair_property_map
    : public boost::put_get_helper<typename Pair::first_type, 
        first_of_pair_property_map<Pair> >
  {
    typedef Pair key_type;
    typedef typename Pair::first_type value_type;
    typedef typename Pair::first_type& reference;
    typedef boost::read_write_property_map_tag category;

    inline value_type operator[](const key_type& v) const { return v.first; }
    inline reference operator[](const key_type& v) { return v.first; }
  };


  template <typename Tuple, int N>
 struct nth_of_tuple_property_map
    : public boost::put_get_helper<typename boost::tuples::element<N,Tuple>::type, 
        nth_of_tuple_property_map<Tuple,N> >
  {
    typedef Tuple key_type;
    typedef typename boost::tuples::element<N,Tuple>::type value_type;
    typedef value_type& reference;
    typedef boost::read_write_property_map_tag category;
    
  };

};

namespace boost
{
  
  template <typename Tuple, int N>
  typename boost::tuples::element<N,Tuple>::type
  get(CGAL::nth_of_tuple_property_map<Tuple,N> pm, const Tuple& k)
  {
    return get<N>(k);
  }
  
  template <typename Tuple, int N>
  typename boost::tuples::element<N,Tuple>::type&
  get(CGAL::nth_of_tuple_property_map<Tuple,N> pm, Tuple& k)
  {
    return get<N>(k);
  }

  template <typename Tuple, int N>
  void
  put(CGAL::nth_of_tuple_property_map<Tuple,N> pm, Tuple& k, const typename boost::tuples::element<N,Tuple>::type& v)
  {
    get<N>(k) = v;
  }

}


#endif // CGAL_PROPERTY_MAP_H
