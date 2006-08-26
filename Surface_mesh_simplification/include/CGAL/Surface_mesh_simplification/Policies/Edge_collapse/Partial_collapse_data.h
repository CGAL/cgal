// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL: $
// $Id: $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_PARTIAL_COLLAPSE_DATA_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_PARTIAL_COLLAPSE_DATA_H

#include <CGAL/Surface_mesh_simplification/Detail/TSMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse
{

template<class TSM_>    
class Partial_collapse_data
{
public:

  typedef TSM_ TSM ;

  typedef typename Surface_geometric_traits<TSM>::FT FT ;

  typedef optional<FT> Optional_cost_type ;
  
public :

  Partial_collapse_data() {}

  Partial_collapse_data( Optional_cost_type aCost )
   :
    mCost     (aCost) 
  {}
  
  Optional_cost_type cost() const { return mCost ; }
  
private :

  Optional_cost_type mCost ;
};    

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_PARTIAL_COLLAPSE_DATA_H
// EOF //
 
