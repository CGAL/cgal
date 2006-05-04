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
// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_COUNT_RATIO_STOP_PRED_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_COUNT_RATIO_STOP_PRED_H 1

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{

//*******************************************************************************************************************
//                                -= stopping condition predicate =-
//
// Determines whether the simplification has finished.
// The arguments are (current_cost,vertex,vertex,is_edge,initial_pair_count,current_pair_count,surface) and the result is bool
//
//*******************************************************************************************************************

// 
// Stops when the ratio of initial to current vertex pairs is below some value.
//
template<class TSM_>    
class Count_ratio_stop_condition
{
public:

    typedef TSM_ TSM ;
    
    typedef typename boost::graph_traits<TSM>::vertex_descriptor vertex_descriptor ;
    typedef typename boost::graph_traits<TSM>::edges_size_type   size_type ;
    
    typedef typename Surface_geometric_traits<TSM>::FT FT ;

public :
  
  Count_ratio_stop_condition( double aRatio ) : mRatio(aRatio) {}
  
  bool operator()( FT                // aCurrentCost
                 , vertex_descriptor // p
                 , vertex_descriptor // q
                 , bool              // aIsEdge
                 , size_type         aInitialCount
                 , size_type         aCurrentCount
                 , TSM&              // aSurface
                 ) const 
  {
    return ( static_cast<double>(aCurrentCount) / static_cast<double>(aInitialCount) ) < mRatio ;
  }
  
private:
  
  double mRatio ;
};    

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_COUNT_RATIO_STOP_PRED_H //
// EOF //
 
