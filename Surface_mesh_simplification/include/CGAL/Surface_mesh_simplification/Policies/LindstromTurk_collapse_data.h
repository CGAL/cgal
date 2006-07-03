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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROMTURK_COLLAPSE_DATA_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROMTURK_COLLAPSE_DATA_H 1

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Minimal_collapse_data.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{


template<class TSM_>    
class LindstromTurk_collapse_data : public Minimal_collapse_data<TSM_>
{
public:

  typedef Minimal_collapse_data<TSM_> Base ;
  
  typedef typename Base::TSM               TSM ;
  typedef typename Base::vertex_descriptor vertex_descriptor;
  typedef typename Base::edge_descriptor   edge_descriptor;
  typedef typename Base::size_type         size_type ;
  typedef typename Base::Point_3           Point_3 ;
  typedef typename Base::FT                FT ;

  typedef optional<Point_3> Optional_point_3 ;  
  typedef optional<FT>      Optional_FT ;
  
  struct Params
  {
    Params()
      :
      VolumeWeight  ( FT(1) / FT(2) )
     ,BoundaryWeight( FT(1) / FT(2) )
     ,ShapeWeight   ( 0 )
    {}
    
    Params( FT const& aVolumeWeight, FT const& aBoundaryWeight, FT const& aShapeWeight )
      :
      VolumeWeight  (aVolumeWeight)
     ,BoundaryWeight(aBoundaryWeight)
     ,ShapeWeight   (aShapeWeight)
    {}
      
    FT VolumeWeight ;
    FT BoundaryWeight ;
    FT ShapeWeight ;
  };
  
public :

  LindstromTurk_collapse_data() {}
  
  LindstromTurk_collapse_data ( vertex_descriptor const&  aP 
                              , vertex_descriptor const&  aQ
                              , bool                      aIsPFixed
                              , bool                      aIsQFixed
                              , edge_descriptor const&    aEdge 
                              , TSM&                      aSurface 
                              )
    :
     Base(aP,aQ,aIsPFixed,aIsQFixed,aEdge,aSurface)
  {} 
  
  LindstromTurk_collapse_data ( vertex_descriptor const&  aP 
                              , vertex_descriptor const&  aQ
                              , bool                      aIsPFixed
                              , bool                      aIsQFixed
                              , edge_descriptor const&    aEdge 
                              , TSM&                      aSurface 
                              , Optional_FT const&        aCost
                              , Optional_point_3 const&   aVertexPoint
                              )
    :
     Base(aP,aQ,aIsPFixed,aIsQFixed,aEdge,aSurface)
    ,mCost        (aCost)
    ,mVertexPoint (aVertexPoint)
  {} 

  Optional_FT      cost        () const { return mCost ; }
  Optional_point_3 vertex_point() const { return mVertexPoint ; }
  
private:

  Optional_FT      mCost ;
  Optional_point_3 mVertexPoint ;
  
};    

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROMTURK_COLLAPSE_DATA_H //
// EOF //
 
