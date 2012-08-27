// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

#include <CGAL/Cartesian_converter.h>

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_params.h>

namespace CGAL {

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

template<class ECM_, class Kernel_ = typename Kernel_traits< typename halfedge_graph_traits<ECM_>::Point>::Kernel >
class LindstromTurkCore
{
public:
    
  typedef ECM_    ECM ;
  typedef Kernel_ Kernel ;
  
  typedef Edge_profile<ECM> Profile ;
  
  typedef boost::graph_traits<ECM> GraphTraits ;
  
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor ;
  typedef typename GraphTraits::edge_descriptor   edge_descriptor ;
  typedef typename GraphTraits::in_edge_iterator  in_edge_iterator ;
  
  typedef LindstromTurk_params Params ;
  
  typedef typename Kernel::Point_3 Point ;

  typedef typename halfedge_graph_traits<ECM>::Point ECM_Point ;
  
  typedef typename Kernel_traits<ECM_Point>::Kernel ECM_Kernel ;
  
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

  Point get_point ( vertex_descriptor const& v ) const 
  {
    return convert(get(vertex_point,surface(),v));
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

  static bool is_finite ( FT     const& n ) { return CGAL_NTS is_finite(n) ; }
  static bool is_finite ( Point  const& p ) { return is_finite(p.x())  && is_finite(p.y())  && is_finite(p.z()) ;  }
  static bool is_finite ( Vector const& v ) { return is_finite(v.x())  && is_finite(v.y())  && is_finite(v.z()) ;  }
  static bool is_finite ( Matrix const& m ) { return is_finite(m.r0()) && is_finite(m.r1()) && is_finite(m.r2()) ; }

  template<class T>
  static optional<T> filter_infinity ( T const& n ) { return is_finite(n) ? optional<T>(n) : optional<T>() ; }

  ECM& surface() const { return mProfile.surface() ; }
  
private:    

  Params const&  mParams ; 
  Profile const& mProfile ;

    void Add_constraint_if_alpha_compatible( Vector const& Ai, FT const& bi ) ;
  
    void Add_constraint_from_gradient ( Matrix const& H, Vector const& c ) ;

private:    

  Triangle_data_vector mTriangle_data ;
  Boundary_data_vector mBdry_data ;
  
  int    mConstraints_n ;
  Matrix mConstraints_A ;
  Vector mConstraints_b ;

  Cartesian_converter<ECM_Kernel,Kernel> convert ;

  FT mSquared_cos_alpha;
  FT mSquared_sin_alpha;
};

} // namespace Surface_mesh_simplification

} //namespace CGAL

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core_impl.h>

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H //
// EOF //
 
