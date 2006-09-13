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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_CACHED_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_CACHED_PLACEMENT_H

#include <CGAL/Surface_mesh_simplification/Detail/TSMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse
{

template<class TSM_>
class Cached_placement
{
  
public:
  
  typedef TSM_ TSM ;
  
  typedef typename boost::graph_traits<TSM>::edge_descriptor edge_descriptor ;
  
  typedef typename Geometric_graph_traits<TSM>::Point Point_3 ;
  
  typedef optional<Point_3> result_type ;

  typedef void Params ;
  
public :
    
  template<class Collapse_data>
  result_type operator()( edge_descriptor const& aEdge
                        , TSM&                   aSurface
                        , Collapse_data const&   aData
                        , Params const*          aParams
                        ) const
  {
    return aData.placement();
  }
};

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_CACHED_PLACEMENT_H
// EOF //
 
