// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_FUNCTORS_FOR_FACE_GRAPH_WRAPPER_H
#define CGAL_FUNCTORS_FOR_FACE_GRAPH_WRAPPER_H 1

#include <CGAL/license/Surface_mesh_topology.h>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/helpers.h>

////////////////////////////////////////////////////////////////////////////////
/** This file contains the following functors for Face_graph_wrapper:
 * Get_beta<typename HEG, unsigned int i>::
         operator() (const HEG& heg, Dart_const_handle dh)
 * Index_from_halfedge_descriptor<Mesh>::run(m, h)
 * Halfedge_descriptor_from_index<Mesh>::run(m, i)
 * Is_index_used<Mesh>::run(m, i)
*/
////////////////////////////////////////////////////////////////////////////////
namespace CGAL {
////////////////////////////////////////////////////////////////////////////////
template<typename P>
class Surface_mesh;
////////////////////////////////////////////////////////////////////////////////
namespace internal {
////////////////////////////////////////////////////////////////////////////////
/// Get_beta
template<typename HEG, unsigned int i>
struct Get_beta
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

  static Dart_const_handle value(const HEG& /*heg*/, Dart_const_handle /*dh*/)
  {
    std::cout<<"ERROR Get_beta<HEG, "<<i<<">"<<std::endl;
    CGAL_assertion(false);
    return nullptr;
  }
};
template<typename HEG>
struct Get_beta<HEG, 0>
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  static Dart_const_handle value(const HEG& heg, Dart_const_handle dh)
  { return prev(dh, heg); }
};
template<typename HEG>
struct Get_beta<HEG, 1>
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  static Dart_const_handle value(const HEG& heg, Dart_const_handle dh)
  { return next(dh, heg); }
};
template<typename HEG>
struct Get_beta<HEG, 2>
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  static Dart_const_handle value(const HEG& heg, Dart_const_handle dh)
  { return opposite(dh, heg); }
};
////////////////////////////////////////////////////////////////////////////////
template<typename Mesh>
struct Index_from_halfedge_descriptor
{
  typedef boost::uint32_t size_type;
  typedef typename boost::template graph_traits<Mesh>::halfedge_descriptor
  halfedge_descriptor;

  static size_type run(const Mesh& m, halfedge_descriptor h)
  {
    size_type res=0;
    for (typename boost::template graph_traits<Mesh>::halfedge_iterator
           it=halfedges(m).begin(), itend=halfedges(m).end(); it!=itend; ++it, ++res)
    { if ((*it)==h) { return res; } }
    return (std::numeric_limits<size_type>::max)();
  }
};
template<typename P>
struct Index_from_halfedge_descriptor<CGAL::Surface_mesh<P> >
{
  using Mesh=CGAL::Surface_mesh<P>;
  typedef boost::uint32_t size_type;
  typedef typename boost::template graph_traits<Mesh>::halfedge_descriptor
  halfedge_descriptor;

  static size_type run(const Mesh& /*m*/, halfedge_descriptor h)
  { return (size_type)(h); }
};
////////////////////////////////////////////////////////////////////////////////
template<typename Mesh>
struct Halfedge_descriptor_from_index
{
  typedef boost::uint32_t size_type;
  typedef typename boost::template graph_traits<Mesh>::halfedge_descriptor
  halfedge_descriptor;

  static halfedge_descriptor run(const Mesh& m, size_type i)
  {
    for (typename boost::template graph_traits<Mesh>::halfedge_iterator
           it=halfedges(m).begin(), itend=halfedges(m).end(); it!=itend; ++it, --i)
    { if (i==0) { return *it; } }
    return halfedge_descriptor(nullptr);
  }
};
template<typename P>
struct Halfedge_descriptor_from_index<CGAL::Surface_mesh<P> >
{
  using Mesh=CGAL::Surface_mesh<P>;
  typedef boost::uint32_t size_type;
  typedef typename boost::template graph_traits<Mesh>::halfedge_descriptor
  halfedge_descriptor;

  static halfedge_descriptor run(const Mesh& /*m*/, size_type i)
  { return halfedge_descriptor(i); }
};
////////////////////////////////////////////////////////////////////////////////
template<typename Mesh>
struct Is_index_used
{
  typedef boost::uint32_t size_type;

  static bool run(const Mesh& m, size_type i)
  { return i<m.size_of_halfedges(); }
};
template<typename P>
struct Is_index_used<CGAL::Surface_mesh<P> >
{
  using Mesh=CGAL::Surface_mesh<P>;
  typedef boost::uint32_t size_type;

  static bool run(const Mesh& m, size_type i)
  { return i<(m.number_of_halfedges()+m.number_of_removed_halfedges()) &&
             !m.is_removed(typename Mesh::Halfedge_index(i)); }
};
////////////////////////////////////////////////////////////////////////////////
} // namespace internal
} // namespace CGAL

#endif // CGAL_FUNCTORS_FOR_FACE_GRAPH_WRAPPER_H //
// EOF //
