// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_COMMON_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_COMMON_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/algorithm.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/optional/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <functional>
#include <utility>
#include <vector>
#include <vector>
#include <set>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class Handle>
inline bool handle_assigned(Handle h) { Handle null; return h != null; }

template<class Iterator, class Handle>
bool handle_exists(Iterator begin, Iterator end, Handle h)
{
  if(handle_assigned(h))
  {
    while(begin != end)
      if(begin++ == h)
        return true;
  }

  return false;
}

template <class TM>
struct No_constrained_edge_map
{
  typedef typename boost::graph_traits<TM>::edge_descriptor         key_type;
  typedef bool                                                      value_type;
  typedef value_type                                                reference;
  typedef boost::readable_property_map_tag                          category;
  friend bool get(No_constrained_edge_map, key_type) { return false; }
};

} // namespace Surface_mesh_simplification

template<class N>
inline std::string n_to_string(const N& n) {
  return boost::str(boost::format("%|5.19g|") % n);
}

template<class XYZ>
inline std::string xyz_to_string(const XYZ& xyz) {
  return boost::str(boost::format("(%|5.19g|,%|5.19g|,%|5.19g|)") % xyz.x() % xyz.y() % xyz.z());
}

template<class Matrix>
inline std::string matrix_to_string(const Matrix& m) {
  return boost::str(boost::format("[%1%|%2%|%3%]") % xyz_to_string(m.r0()) % xyz_to_string(m.r1()) % xyz_to_string(m.r2()));
}

template<class T>
inline std::string optional_to_string(const boost::optional<T>& o) {
  if(o)
    return boost::str(boost::format("%1%") % *o);
  else return std::string("NONE");
}

} // namespace CGAL

#if defined(CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE) \
  || defined(CGAL_SURFACE_SIMPLIFICATION_ENABLE_LT_TRACE)
#define CGAL_SMS_ENABLE_TRACE
#endif

#ifdef CGAL_SMS_ENABLE_TRACE
  #include<string>
  #include<iostream>
  #include<sstream>

  namespace internal { namespace { bool cgal_enable_sms_trace = true; } }
  #define CGAL_SMS_TRACE_IMPL(m) \
    if(::internal::cgal_enable_sms_trace) { \
    std::ostringstream ss; ss << m; std::string s = ss.str(); \
    /*Surface_simplification_external_trace(s)*/ std::cerr << s << std::endl; \
    }

  #define CGAL_SMS_DEBUG_CODE(code) code
#else
  #define CGAL_SMS_DEBUG_CODE(code)
#endif

#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_LT_TRACE
  #define CGAL_SMS_LT_TRACE(l,m) if((l) <= CGAL_SURFACE_SIMPLIFICATION_ENABLE_LT_TRACE) CGAL_SMS_TRACE_IMPL(m)
#else
  #define CGAL_SMS_LT_TRACE(l,m)
#endif

#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
  #define CGAL_SMS_TRACE_IF(c,l,m) if((c) && ((l) <= CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE)) CGAL_SMS_TRACE_IMPL(m)
  #define CGAL_SMS_TRACE(l,m)      if((l) <= CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE) CGAL_SMS_TRACE_IMPL(m)
#else
  #define CGAL_SMS_TRACE_IF(c,l,m)
  #define CGAL_SMS_TRACE(l,m)
#endif
#undef CGAL_SMS_ENABLE_TRACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_COMMON_H
