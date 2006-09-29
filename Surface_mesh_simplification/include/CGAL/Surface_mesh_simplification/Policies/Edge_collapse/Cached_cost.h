// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Surface_mesh_simplification/Policies/Edge_length_cost.h $
// $Id: Edge_length_cost.h 32188 2006-07-03 17:19:46Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_CACHED_COST_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_CACHED_COST_H

#include <CGAL/Surface_mesh_simplification/Detail/ECMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Surface_mesh_simplification
{

template<class ECM_>
class Cached_cost 
{
  
public:
  
  typedef ECM_ ECM ;
  
  typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
  
  typedef typename halfedge_graph_traits<ECM>::Point Point ;
  typedef typename Kernel_traits<Point>::Kernel      Kernel ;
  typedef typename Kernel::FT                        FT ;
  
  typedef optional<FT> result_type ;

  typedef void Params ;
  
public :

  template<class Cache>    
  result_type operator()( edge_descriptor const& aEdge
                        , ECM&                   aSurface
                        , Cache const&           aCache
                        , Params const*          aParams
                        ) const
  {
    return aCache.cost();
  }
};

} // namespace Surface_mesh_simplification


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_CACHED_COST_H
// EOF //
 
