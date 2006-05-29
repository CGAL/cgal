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
// $URL: $
// $Id: $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROM_TURK_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROM_TURK_H 1

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh_simplification/TSMS_common.h>
#include <CGAL/Surface_mesh_simplification/Policies/LindstromTurk_collapse_data.h>
#include <CGAL/Cartesian/MatrixC33.h>

CGAL_BEGIN_NAMESPACE

//
// This should be in 
//
// Implementation of the collapsing cost and placement strategy from:
//
//  "Fast and Memory Efficient Polygonal Symplification"
//  Peter Lindstrom, Greg Turk
//

namespace Triangulated_surface_mesh { namespace Simplification 
{

template<class Collapse_data_>
class LindstromTurkImpl
{
public:
    
  typedef Collapse_data_ Collapse_data ;
  
  typedef typename Collapse_data::TSM               TSM ;
  typedef typename Collapse_data::vertex_descriptor vertex_descriptor ;
  typedef typename Collapse_data::edge_descriptor   edge_descriptor ;
  typedef typename Collapse_data::Params            Params ;
     
  typedef typename Surface_geometric_traits<TSM>::Kernel Kernel ;
  
  typedef typename Kernel::Point_3  Point ;
  typedef typename Kernel::Vector_3 Vector ;
  typedef typename Kernel::FT       FT ;
 
  typedef optional<FT>    Optional_FT ;
  typedef optional<Point> Optional_point ;
  
  typedef Simple_cartesian<FT> CKernel ;
  typedef MatrixC33<CKernel>   CMatrix ;
  typedef Vector_3<CKernel>    CVector ; 
  
    
public:
  
  typedef shared_ptr<Collapse_data> result_type ;
    
  LindstromTurkImpl( Params const&            aParams
                   , vertex_descriptor const& aP
                   , vertex_descriptor const& aQ
                   , edge_descriptor const&   aP_Q
                   , edge_descriptor const&   aQ_P
                   , TSM&                     aSurface 
                   ) ;
    
  result_type result() const { return mResult ; }   
  
private :

  struct Triangle
  {
    Triangle() {}
    
    Triangle( CVector const& aNormalV, FT const& aNormalL ) : NormalV(aNormalV), NormalL(aNormalL) {}
    
    CVector NormalV ;
    FT      NormalL ;
  } ;
  
  typedef std::vector<Triangle>          Triangles ;
  typedef std::vector<vertex_descriptor> Link ;

  struct Boundary
  {
    Boundary ( CVector const& op_
             , CVector const& opN_
             , CVector const& pq_
             , CVector const& pqN_
             , CVector const& qr_
             , CVector const& qrN_
             )
      :
      op(op_), opN(opN_), pq(pq_), pqN(pqN_), qr(qr_), qrN(qrN_)
    {}
      
    CVector op, opN, pq, pqN, qr, qrN  ;
  } ;
  typedef optional<Boundary> OptionalBoundary ;
  
  class Constrians
  {
  public:
  
    Constrians() : n(0), A(NULL_MATRIX), b(NULL_VECTOR) {}

    void Add_if_alpha_compatible( CVector const& Ai, FT const& bi ) ;
  
    void Add_from_gradient ( CMatrix const& H, CVector const& c ) ;
    
    int     n ;
    CMatrix A ;
    CVector b ;
    
  private:
  
    // alpha = 1 degree  
    static double squared_cos_alpha() { return 0.999695413509  ; }
    static double squared_sin_alpha() { return 3.04586490453e-4; }
  } ;
  
private :
    
  void Add_boundary_preservation_constrians( OptionalBoundary const& aBdry ) ;
  void Add_volume_preservation_constrians( Triangles const& aTriangles );
  void Add_boundary_and_volume_optimization_constrians( OptionalBoundary const& aBdry, Triangles const& aTriangles ) ;
  void Add_shape_optimization_constrians( Link const& aLink ) ;

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
  
  static CVector toCVector ( Vector const& v )
  {
    return CVector(v.x(),v.y(),v.z());
  }
  
  static Vector toVector ( CVector const& cv )
  {
    Rational_traits<FT> rat ;
    
    if (   (rat.denominator(cv.x()) != rat.denominator(cv.y()) )
        || (rat.denominator(cv.x()) != rat.denominator(cv.z()) ) 
       )
         return Vector(rat.numerator  (cv.x()) * rat.denominator(cv.y()) * rat.denominator(cv.z())
                      ,rat.numerator  (cv.y()) * rat.denominator(cv.x()) * rat.denominator(cv.z())
                      ,rat.numerator  (cv.z()) * rat.denominator(cv.x()) * rat.denominator(cv.y())
                      ,rat.denominator(cv.x()) * rat.denominator(cv.y()) * rat.denominator(cv.z())
                      );
    else
         return Vector(rat.numerator(cv.x())
                      ,rat.numerator(cv.y())
                      ,rat.numerator(cv.z())
                      ,rat.denominator(cv.x())
                      );
  }
  
  static CVector Point_cross_product ( Point const& a, Point const& b ) 
  {
    return toCVector(cross_product(a-ORIGIN,b-ORIGIN)); 
  }

  // This is the (uX)(Xu) product described in the Lindstrom-Turk paper
  static CMatrix LT_product( CVector const& u ) 
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
  
    return CMatrix(a00,a01,a02
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
  
  OptionalBoundary Extract_boundary() ;

private:    

  Params const&            mParams ; 
  vertex_descriptor const& mP ;
  vertex_descriptor const& mQ ;
  edge_descriptor const&   mP_Q ;
  edge_descriptor const&   mQ_P ;
  TSM&                     mSurface ;

private:    

  Constrians mConstrians ;
  
  result_type mResult;   

};

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#include <CGAL/Surface_mesh_simplification/Policies/Detail/Lindstrom_Turk_impl.h>

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROM_TURK_H //
// EOF //
 
