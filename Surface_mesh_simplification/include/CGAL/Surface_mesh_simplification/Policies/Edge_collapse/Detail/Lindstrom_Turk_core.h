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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H 1

#include <vector>

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>
#include <CGAL/Surface_mesh_simplification/Policies/LindstromTurk_collapse_data.h>

CGAL_BEGIN_NAMESPACE

//
// This should be in 
//
// Implementation of the collapsing cost and placement strategy from:
//
//  "Fast and Memory Efficient Polygonal Symplification"
//  Peter Lindstrom, Greg Turk
//

namespace Triangulated_surface_mesh { namespace Simplification  { namespace Edge_collapse
{

template<class TSM_>
class LindstromTurkCore
{
public:
    
  typedef typename TSM_TSM ;
  
  typedef Collapse_data_ Collapse_data ;
  
  typedef boost::graph_traits<TSM> GraphTraits ;
  
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor ;
  typedef typename GraphTraits::edge_descriptor   edge_descriptor ;
  typedef typename GraphTraits::in_edge_iterator  in_edge_iterator ;
  
  typedef LindstromTurk_params Params ;
  
  typedef typename Surface_geometric_traits<TSM>::Kernel Kernel ;
  
  typedef typename Kernel::Point_3  Point ;
  typedef typename Kernel::Vector_3 Vector ;
  typedef typename Kernel::FT       FT ;
 
  typedef optional<FT>     Optional_FT ;
  typedef optional<Point>  Optional_point ;
  typedef optional<Vector> Optional_vector ;
  
  typedef MatrixC33<Kernel> Matrix ;
  
  typedef tuple<Optional_FT,Optional_point> result_type ;
    
public:
  
  LindstromTurkCore( Params const&          aParams
                   , edge_descriptor const& aP_Q
                   , edge_descriptor const& aQ_P
                   , TSM&                   aSurface 
                   ) ;
    
  result_type compute() ;
  
private :

  struct Triangle
  {
    Triangle() {}
    
    Triangle( Vector const& aNormalV, FT const& aNormalL ) : NormalV(aNormalV), NormalL(aNormalL) {}
    
    Vector NormalV ;
    FT     NormalL ;
  } ;
  
  typedef std::vector<Triangle>          Triangles ;
  typedef std::vector<vertex_descriptor> Link ;
  typedef std::vector<edge_descriptor>   edge_descriptor_vector ;
  
  struct BoundaryEdge
  {
    BoundaryEdge ( Point s_, Point t_, Vector const& v_, Vector const& n_ ) : s(s_), t(t_), v(v_), n(n_) {}

    Point  s, t ;      
    Vector v, n ;
  } ;
  typedef std::vector<BoundaryEdge> BoundaryEdges ;
  
  class Constrians
  {
  public:
  
    Constrians() : n(0), A(NULL_MATRIX), b(NULL_VECTOR) {}

    void Add_if_alpha_compatible( Vector const& Ai, FT const& bi ) ;
  
    void Add_from_gradient ( Matrix const& H, Vector const& c ) ;
    
    int    n ;
    Matrix A ;
    Vector b ;
    
  private:
  
    // alpha = 1 degree  
    static double squared_cos_alpha() { return 0.999695413509  ; }
    static double squared_sin_alpha() { return 3.04586490453e-4; }
  } ;
  
private :
    
  void Add_boundary_preservation_constrians( BoundaryEdges const& aBdry ) ;
  void Add_volume_preservation_constrians( Triangles const& aTriangles );
  void Add_boundary_and_volume_optimization_constrians( BoundaryEdges const& aBdry, Triangles const& aTriangles ) ;
  void Add_shape_optimization_constrians( Link const& aLink ) ;

  FT Compute_boundary_cost( Vector const& v, BoundaryEdges const& aBdry ) ;
  FT Compute_volume_cost  ( Vector const& v, Triangles const& aTriangles ) ;
  FT Compute_shape_cost   ( Point const& p, Link const& aLink ) ;

  bool is_border ( edge_descriptor const& edge ) const
  {
    edge_is_border_t is_border_property ;
    return get(is_border_property,mSurface,edge) ;
  }    
  bool is_undirected_edge_a_border ( edge_descriptor const& edge ) const
  {
    return is_border(edge) || is_border(opposite_edge(edge,mSurface)) ;
  }    

  Point get_point( vertex_descriptor const& v ) const
  {
    vertex_point_t vertex_point_property ;
    return get(vertex_point_property,mSurface,v) ;
  }
  
  
  static Vector Point_cross_product ( Point const& a, Point const& b ) 
  {
    return cross_product(a-ORIGIN,b-ORIGIN); 
  }

  // This is the (uX)(Xu) product described in the Lindstrom-Turk paper
  static Matrix LT_product( Vector const& u ) 
  {
    FT a00 = ( u.y()*u.y() ) + ( u.z()*u.z() ) ;
    FT a01 = -u.x()*u.y();
    FT a02 = -u.x()*u.z();
  
    FT a10 = a01 ;
    FT a11 = ( u.x()*u.x() ) + ( u.z()*u.z() ) ;
    FT a12 = - u.y() * u.z();
  
    FT a20 = a02 ;
    FT a21 = a12 ;
    FT a22 =  ( u.x()*u.x() ) + ( u.y()*u.y() ) ;
  
    return Matrix(a00,a01,a02
                 ,a10,a11,a12
                 ,a20,a21,a22
                 );
  }
    
  Triangle Get_triangle ( vertex_descriptor const& v0
                        , vertex_descriptor const& v1
                        , vertex_descriptor const& v2 
                        ) ;
                       
  void Extract_triangle( vertex_descriptor const& v0
                       , vertex_descriptor const& v1
                       , vertex_descriptor const& v2 
                       , edge_descriptor   const& e02
                       , Triangles&               rTriangles
                       ) ;
                       
  void Extract_triangles_and_link( Triangles& rTriangles, Link& rLink );
  
  void Extract_boundary_edge ( edge_descriptor edge, BoundaryEdges& rBdry ) ;
  void Extract_boundary_edges( vertex_descriptor const& v
                             , edge_descriptor_vector&  rCollected
                             , BoundaryEdges&           rBdry
                             ) ;
  void Extract_boundary_edges( BoundaryEdges& rBdry ) ;

  
private:    

  Params const&            mParams ; 
  vertex_descriptor const& mP ;
  vertex_descriptor const& mQ ;
  bool                     mIsPFixed ;
  bool                     mIsQFixed ;
  edge_descriptor const&   mP_Q ;
  edge_descriptor const&   mQ_P ;
  TSM&                     mSurface ;

private:    

  Constrians mConstrians ;
  
  Point mV ;
  
};

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse

CGAL_END_NAMESPACE

#include <CGAL/Surface_mesh_simplification/Policies/Detail/Lindstrom_Turk_core_impl.h>

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H //
// EOF //
 
