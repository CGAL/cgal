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

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_params.h>

CGAL_BEGIN_NAMESPACE

//
// This should be in 
//
// Implementation of the collapsing cost and placement strategy from:
//
//  "Fast and Memory Efficient Polygonal Symplification"
//  Peter Lindstrom, Greg Turk
//

namespace Surface_mesh_simplification
{

template<class ECM_>
class LindstromTurkCore
{
public:
    
  typedef ECM_ ECM ;
  
  typedef Edge_profile<ECM> Profile ;
  
  typedef boost::graph_traits<ECM> GraphTraits ;
  
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor ;
  typedef typename GraphTraits::edge_descriptor   edge_descriptor ;
  typedef typename GraphTraits::in_edge_iterator  in_edge_iterator ;
  
  typedef LindstromTurk_params Params ;
  
  typedef typename halfedge_graph_traits<ECM>::Point Point ;
  
  typedef typename Kernel_traits<Point>::Kernel Kernel ;
  
  typedef typename Kernel::Vector_3 Vector ;
  typedef typename Kernel::FT       FT ;
 
  typedef optional<FT>     Optional_FT ;
  typedef optional<Point>  Optional_point ;
  typedef optional<Vector> Optional_vector ;
  
  typedef MatrixC33<Kernel> Matrix ;
  
  typedef typename Profile::Triangle                 Triangle ;
  typedef typename Profile::vertex_descriptor_vector vertex_descriptor_vector ;
  
  typedef typename Profile::Triangle_vector       ::const_iterator const_triangle_iterator ;
  typedef typename Profile::edge_descriptor_vector::const_iterator const_border_edge_iterator ;
  
public:
  
  LindstromTurkCore( Params const& aParams, Profile const& aProfile ) ;
    
  Optional_point compute_placement() ;
  Optional_FT    compute_cost( Optional_point const& p ) ;
  
private :

  struct Triangle_data
  {
    Triangle_data( Vector const& aNormalV, FT const& aNormalL ) : NormalV(aNormalV), NormalL(aNormalL) {}
    
    Vector NormalV ;
    FT     NormalL ;
  } ;
  struct Boundary_data
  {
    Boundary_data ( Point s_, Point t_, Vector const& v_, Vector const& n_ ) : s(s_), t(t_), v(v_), n(n_) {}

    Point  s, t ;      
    Vector v, n ;
  } ;
  typedef std::vector<Triangle_data> Triangle_data_vector ;
  typedef std::vector<Boundary_data> Boundary_data_vector ;
  
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
    static FT squared_cos_alpha() { return FT(0.999695413509)  ; }
    static FT squared_sin_alpha() { return FT(3.04586490453e-4); }
  } ;
  
private :
    
  void Extract_triangle_data();
  void Extract_boundary_data();
  
  void Add_boundary_preservation_constrians( Boundary_data_vector const& aBdry ) ;
  void Add_volume_preservation_constrians( Triangle_data_vector const& aTriangles );
  void Add_boundary_and_volume_optimization_constrians( Boundary_data_vector const& aBdry, Triangle_data_vector const& aTriangles ) ;
  void Add_shape_optimization_constrians( vertex_descriptor_vector const& aLink ) ;

  FT Compute_boundary_cost( Vector const& v, Boundary_data_vector const& aBdry ) ;
  FT Compute_volume_cost  ( Vector const& v, Triangle_data_vector const& aTriangles ) ;
  FT Compute_shape_cost   ( Point  const& p, vertex_descriptor_vector const& aLink ) ;

  Point const& get_point ( vertex_descriptor const& v ) const 
  {
    return get(vertex_point,surface(),v);
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

  static FT big_value() { return static_cast<FT>((std::numeric_limits<double>::max)()) ; }
  
  // Returns 'n' if it is finite, CGALi::infinite otherwise.
  static FT made_finite ( FT n ) { return ( CGAL_NTS is_finite(n) ) ? n : big_value() ; }
  
  static Vector made_finite ( Vector const& v ) { return Vector( made_finite(v.x()), made_finite(v.y()), made_finite(v.z()) ) ; }
  
  ECM& surface() const { return mProfile.surface() ; }
  
private:    

  Params const&  mParams ; 
  Profile const& mProfile ;

private:    

  Triangle_data_vector mTriangle_data ;
  Boundary_data_vector mBdry_data ;
  
  Constrians mConstrians ;
};

} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core_impl.h>

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H //
// EOF //
 
