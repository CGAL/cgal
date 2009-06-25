// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Surface_mesh_simplification/edge_collapse.h $
// $Id: edge_collapse.h 49869 2009-06-10 15:51:32Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_VISITOR_BASE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_VISITOR_BASE_H

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

CGAL_BEGIN_NAMESPACE

namespace Surface_mesh_simplification
{

template<class ECM_>
struct Edge_collapse_visitor_base
{
  typedef ECM_ ECM ;
  
  typedef Edge_profile<ECM> Profile ;
  
  typedef boost::graph_traits  <ECM> GraphTraits ; 
  typedef halfedge_graph_traits<ECM> HalfedgeGraphTraits ; 
  
  typedef typename GraphTraits::edges_size_type   size_type ;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor ;
  typedef typename HalfedgeGraphTraits::Point     Point ;
  typedef typename Kernel_traits<Point>::Kernel   Kernel ;
  typedef typename Kernel::FT                     FT ;
  
  virtual void OnStarted( ECM& ) const {} 
  
  virtual void OnFinished ( ECM& ) const {} 
  
  virtual void OnStopConditionReached( Profile const& ) const {} 
  
  virtual void OnCollected( Profile const&, boost::optional<FT> const& ) const {}                
  
  virtual void OnSelected( Profile const&, boost::optional<FT> const&, size_type, size_type ) const {}                
  
  virtual void OnCollapsing(Profile const&, boost::optional<Point> const& ) const {}                
  
  virtual void OnCollapsed( Profile const&, vertex_descriptor const& ) const {}

  virtual void OnNonCollapsable(Profile const& ) const {}                
} ;

} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_VISITOR_BASE_H //
// EOF //
 
