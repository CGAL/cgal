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
// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_MINIMAL_COLLAPSE_DATA_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_MINIMAL_COLLAPSE_DATA_H

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{

template<class TSM_>    
class Minimal_collapse_data
{
public:

  typedef TSM_ TSM ;
  
  typedef typename boost::graph_traits<TSM>::vertex_descriptor vertex_descriptor ;
  typedef typename boost::graph_traits<TSM>::edge_descriptor   edge_descriptor ;
  typedef typename boost::graph_traits<TSM>::edges_size_type   size_type ;
  
  typedef typename Surface_geometric_traits<TSM>::Point_3 Point_3 ;
  typedef typename Surface_geometric_traits<TSM>::FT      FT ;

public :

  Minimal_collapse_data()
    :
     mP       ()
    ,mQ       ()
    ,mIsPFixed(false)
    ,mIsQFixed(false)
    ,mEdge    ()
    ,mSurface (0)                             
  {}
  
  Minimal_collapse_data ( vertex_descriptor const& aP 
                        , vertex_descriptor const& aQ
                        , bool                     aIsPFixed
                        , bool                     aIsQFixed
                        , edge_descriptor   const& aEdge 
                        , TSM&                     aSurface 
                        )
    :
     mP       (aP)
    ,mQ       (aQ)
    ,mIsPFixed(aIsPFixed)
    ,mIsQFixed(aIsQFixed)
    ,mEdge    (aEdge)
    ,mSurface (addressof(aSurface))                             
  {} 

  void set( vertex_descriptor const& aP 
          , vertex_descriptor const& aQ
          , bool                     aIsPFixed
          , bool                     aIsQFixed
          , edge_descriptor   const& aEdge 
          , TSM&                     aSurface 
          )
  {
    mP        = aP ;
    mQ        = aQ ;
    mIsPFixed = aIsPFixed ;
    mIsQFixed = aIsQFixed ;
    mEdge     = aEdge ;
    mSurface  = addressof(aSurface);
  } 
  
  vertex_descriptor const& p()       const { return mP ; }
  vertex_descriptor const& q()       const { return mQ ; }
  edge_descriptor const&   edge()    const { return mEdge ; }
  TSM&                     surface() const { return *mSurface ; }
  
  bool is_p_fixed   () const { return mIsPFixed ; }
  bool is_q_fixed   () const { return mIsQFixed ; }
  bool is_edge_fixed() const { return mIsPFixed && mIsQFixed ; }
  
protected:

  vertex_descriptor mP ;
  vertex_descriptor mQ ;
  bool              mIsPFixed ;
  bool              mIsQFixed ;
  edge_descriptor   mEdge ;
  TSM*              mSurface ;
  
};    

} } // namespace Triangulated_surface_mesh::Simplification


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_MINIMAL_COLLAPSE_DATA_H
// EOF //
 
