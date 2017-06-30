// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
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
//
// Author(s)     : Yin Xu, Andreas Fabri and Ilker O. Yaz

#ifndef CGAL_PMP_WEIGHTS_H
#define CGAL_PMP_WEIGHTS_H

#include <CGAL/license/Polygon_mesh_processing/core.h>

/// @cond CGAL_DOCUMENT_INTERNAL

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace internal {

// Returns the cotangent value of half angle v0 v1 v2
// using formula in -[Meyer02] Discrete Differential-Geometry Operators for- page 19
// The potential problem with previous one (Cotangent_value) is that it does not produce symmetric results
// (i.e. for v0, v1, v2 and v2, v1, v0 returned cot weights can be slightly different)
// This one provides stable results.
template<class PolygonMesh>
struct Cotangent_value_Meyer_impl
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  template <class VertexPointMap>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, const VertexPointMap& ppmap)
  {
    typedef typename Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type >::Kernel::Vector_3 Vector;
    
    Vector a = get(ppmap, v0) - get(ppmap, v1);
    Vector b = get(ppmap, v2) - get(ppmap, v1);
    
    double dot_ab = to_double(a*b);
    // rewritten for safer fp arithmetic
    //double dot_aa = a.squared_length();
    //double dot_bb = b.squared_length();
    //double divider = CGAL::sqrt( dot_aa * dot_bb - dot_ab * dot_ab );

    Vector cross_ab = CGAL::cross_product(a, b);
    double divider = to_double(CGAL::approximate_sqrt(cross_ab*cross_ab));

    if(divider == 0 /*|| divider != divider*/) 
    {
      CGAL::collinear(get(ppmap, v0), get(ppmap, v1), get(ppmap, v2)) ? 
        CGAL_warning(!"Infinite Cotangent value with degenerate triangle!") :
        CGAL_warning(!"Infinite Cotangent value due to floating point arithmetic!");
      

      return dot_ab > 0 ? (std::numeric_limits<double>::max)() :
                         -(std::numeric_limits<double>::max)();
    }
    
    return dot_ab / divider;
  }
};

// Same as above but with a different API
template<class PolygonMesh
  , class VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
class Cotangent_value_Meyer
{
protected:
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef VertexPointMap Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3  Vector;

  PolygonMesh& pmesh_;
  Point_property_map ppmap_;

public:

  Cotangent_value_Meyer(PolygonMesh& pmesh_, VertexPointMap vpmap_)
    : pmesh_(pmesh_)
    , ppmap_(vpmap_)
  {}

  PolygonMesh& pmesh()
  {
    return pmesh_;
  }

  Point_property_map& ppmap()
  {
    return ppmap_;
  }

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    return Cotangent_value_Meyer_impl<PolygonMesh>()(v0, v1, v2, ppmap());
  }
};


// imported from skeletonization
template<class PolygonMesh,
         class VertexPointMap =
          typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
class Cotangent_value_Meyer_secure
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef VertexPointMap Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3  Vector;

  PolygonMesh& pmesh_;
  Point_property_map ppmap_;

public:

  Cotangent_value_Meyer_secure(PolygonMesh& pmesh_, VertexPointMap vpmap_)
    : pmesh_(pmesh_)
    , ppmap_(vpmap_)
  {}

  PolygonMesh& pmesh()
  {
    return pmesh_;
  }

  Point_property_map& ppmap()
  {
    return ppmap_;
  }

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    Vector a = get(ppmap(), v0) - get(ppmap(), v1);
    Vector b = get(ppmap(), v2) - get(ppmap(), v1);

    double dot_ab = a*b;
    double dot_aa = a.squared_length();
    double dot_bb = b.squared_length();
    double lb = -0.999, ub = 0.999;
    double cosine = dot_ab / CGAL::sqrt(dot_aa) / CGAL::sqrt(dot_bb);
    cosine = (cosine < lb) ? lb : cosine;
    cosine = (cosine > ub) ? ub : cosine;
    double sine = std::sqrt(1.0 - cosine * cosine);
    return cosine / sine;
  }
};

// Returns the cotangent value of half angle v0 v1 v2 by clamping between [1, 89] degrees
// as suggested by -[Friedel] Unconstrained Spherical Parameterization-
template<class PolygonMesh
  , class VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type
  , class CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Cotangent_value_clamped : CotangentValue
{
  Cotangent_value_clamped()
  {}
public:

  Cotangent_value_clamped(PolygonMesh& pmesh_, VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  {}

  PolygonMesh& pmesh()
  {
    return CotangentValue::pmesh();
  }

  VertexPointMap& ppmap()
  {
    return CotangentValue::ppmap();
  }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    const double cot_1 = 57.289962;
    const double cot_89 = 0.017455;
    double value = CotangentValue::operator()(v0, v1, v2);
    return (std::max)(cot_89, (std::min)(value, cot_1));
  }
};

template<class PolygonMesh
  , class VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type
  , class CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Cotangent_value_clamped_2 : CotangentValue
{
  Cotangent_value_clamped_2()
  {}

public:

  Cotangent_value_clamped_2(PolygonMesh& pmesh_, VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  {}

  PolygonMesh& pmesh()
  {
    return CotangentValue::pmesh();
  }

  VertexPointMap& ppmap()
  {
    return CotangentValue::ppmap();
  }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    const double cot_5 = 5.671282;
    const double cot_175 = -cot_5;
    double value = CotangentValue::operator()(v0, v1, v2);
    return (std::max)(cot_175, (std::min)(value, cot_5));
  }
};

template<class PolygonMesh
  , class VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type
  , class CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Cotangent_value_minimum_zero : CotangentValue
{
  Cotangent_value_minimum_zero()
  {}
public:
  Cotangent_value_minimum_zero(PolygonMesh& pmesh_, VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  {}

  PolygonMesh& pmesh()
  {
    return CotangentValue::pmesh();
  }

  VertexPointMap& ppmap()
  {
    return CotangentValue::ppmap();
  }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    double value = CotangentValue::operator()(v0, v1, v2);
    return (std::max)(0.0, value);
  }
};

template<class PolygonMesh,
         class CotangentValue = Cotangent_value_Meyer_impl<PolygonMesh> >
struct Cotangent_value_minimum_zero_impl : CotangentValue
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  template <class VertexPointMap>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, const VertexPointMap& ppmap)
  {
    double value = CotangentValue::operator()(v0, v1, v2,ppmap);
    return (std::max)(0.0, value);
  }
};

template<class PolygonMesh
    , class VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type
    , class CotangentValue
           = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Voronoi_area : CotangentValue
{
  //Voronoi_area()
  //{}
  
public:
  Voronoi_area(PolygonMesh& pmesh_, VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  {}

  //Voronoi_area(PolygonMesh& pmesh_)
  //  : CotangentValue(pmesh_, get(CGAL::vertex_point, pmesh_))
  //{}

  PolygonMesh& pmesh()
  {
    return CotangentValue::pmesh();
  }

  VertexPointMap& ppmap()
  {
    return CotangentValue::ppmap();
  }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::in_edge_iterator in_edge_iterator;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  typedef VertexPointMap Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3  Vector;

  double operator()(vertex_descriptor v0) {

    //return 1.0;
    double voronoi_area = 0.0;
    BOOST_FOREACH(halfedge_descriptor he,
                  halfedges_around_target( halfedge(v0,pmesh()), pmesh()) )
    {
      if( is_border(he,pmesh()) ) { continue; }

      CGAL_assertion(CGAL::is_triangle_mesh(pmesh()));
      CGAL_assertion( v0 == target(he, pmesh()) );
      vertex_descriptor v1 = source(he, pmesh());
      vertex_descriptor v_op = target(next(he, pmesh()), pmesh());

      const Point& v0_p = get(ppmap(), v0);
      const Point& v1_p = get(ppmap(), v1);
      const Point& v_op_p = get(ppmap(), v_op);

      // (?) check if there is a better way to predicate triangle is obtuse or not
      CGAL::Angle angle0 = CGAL::angle(v1_p, v0_p, v_op_p);
      CGAL::Angle angle1 = CGAL::angle(v_op_p, v1_p, v0_p);
      CGAL::Angle angle_op = CGAL::angle(v0_p, v_op_p, v1_p);

      bool obtuse = (angle0 == CGAL::OBTUSE) || (angle1 == CGAL::OBTUSE) || (angle_op == CGAL::OBTUSE);

      if(!obtuse) {
        double cot_v1 = CotangentValue::operator()(v_op, v1, v0);
        double cot_v_op = CotangentValue::operator()(v0, v_op, v1);

        double term1 = cot_v1   * to_double((v_op_p - v0_p).squared_length());
        double term2 = cot_v_op * to_double((v1_p  - v0_p).squared_length());
        voronoi_area += (1.0 / 8.0) * (term1 + term2);
      }
      else {
        double area_t = to_double(CGAL::approximate_sqrt(squared_area(v0_p, v1_p, v_op_p)));
        if(angle0 == CGAL::OBTUSE) {
          voronoi_area += area_t / 2.0;
        }
        else {
          voronoi_area += area_t / 4.0;
        }
      }
    }
    CGAL_warning(voronoi_area != 0 && "Zero voronoi area!");
    return voronoi_area;
  }
};
// Returns the cotangent value of half angle v0 v1 v2 by dividing the triangle area
// as suggested by -[Mullen08] Spectral Conformal Parameterization-
template<class PolygonMesh
  , class VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type
  , class CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Cotangent_value_area_weighted : CotangentValue
{
  Cotangent_value_area_weighted()
  {}

public:

  Cotangent_value_area_weighted(PolygonMesh& pmesh_, VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  {}

  PolygonMesh& pmesh()
  {
    return CotangentValue::pmesh();
  }

  VertexPointMap& ppmap()
  {
    return CotangentValue::ppmap();
  }

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    return CotangentValue::operator()(v0, v1, v2)
      / CGAL::sqrt(squared_area(v0->point(), v1->point(), v2->point()));
  }
};
/////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////// Edge Weight Calculators ///////////////////////////////////
// Cotangent weight calculator 
// Cotangent_value:               as suggested by -[Sorkine07] ARAP Surface Modeling-
// Cotangent_value_area_weighted: as suggested by -[Mullen08] Spectral Conformal Parameterization-
template< class PolygonMesh,
          class CotangentValue = Cotangent_value_minimum_zero_impl<PolygonMesh> >
struct Cotangent_weight_impl : CotangentValue
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  // Returns the cotangent weight of specified halfedge_descriptor
  // Edge orientation is trivial
  template<class VertexPointMap>
  double operator()(halfedge_descriptor he, PolygonMesh& pmesh, const VertexPointMap& ppmap)
  {
    vertex_descriptor v0 = target(he, pmesh);
    vertex_descriptor v1 = source(he, pmesh);
     // Only one triangle for border edges
    if (is_border_edge(he, pmesh))
     {
       
       halfedge_descriptor he_cw = opposite( next(he, pmesh) , pmesh );
       vertex_descriptor v2 = source(he_cw, pmesh);
       if (is_border_edge(he_cw, pmesh) )
       {
         halfedge_descriptor he_ccw = prev( opposite(he, pmesh) , pmesh );
         v2 = source(he_ccw, pmesh);
       }
       return ( CotangentValue::operator()(v0, v2, v1, ppmap)/2.0 );
     }
     else
     {
       halfedge_descriptor he_cw = opposite( next(he, pmesh) , pmesh );
       vertex_descriptor v2 = source(he_cw, pmesh);     
       halfedge_descriptor he_ccw = prev( opposite(he, pmesh) , pmesh );
       vertex_descriptor v3 = source(he_ccw, pmesh);

        return ( CotangentValue::operator()(v0, v2, v1, ppmap)/2.0 +
                 CotangentValue::operator()(v0, v3, v1, ppmap)/2.0 );
     }
  }
};

template<class PolygonMesh
  , class VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type
  , class CotangentValue
           = Cotangent_value_minimum_zero<PolygonMesh, VertexPointMap> >
class Cotangent_weight : CotangentValue
{
  Cotangent_weight()
  {}

public:
  Cotangent_weight(PolygonMesh& pmesh_, VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  {}

  Cotangent_weight(PolygonMesh& pmesh_)
    : CotangentValue(pmesh_, get(CGAL::vertex_point, pmesh_))
  {}

  PolygonMesh& pmesh()
  {
    return CotangentValue::pmesh();
  }

  VertexPointMap& ppmap()
  {
    return CotangentValue::ppmap();
  }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<PolygonMesh,vertex_point_t>::type Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3  Vector;

  // Returns the cotangent weight of specified halfedge_descriptor
  // Edge orientation is trivial
  double operator()(halfedge_descriptor he)
  {
    vertex_descriptor v0 = target(he, pmesh());
    vertex_descriptor v1 = source(he, pmesh());
     // Only one triangle for border edges
    if (is_border_edge(he, pmesh()))
     {
       
       halfedge_descriptor he_cw = opposite( next(he, pmesh()) , pmesh() );
       vertex_descriptor v2 = source(he_cw, pmesh());
       if (is_border_edge(he_cw, pmesh()) )
       {
         halfedge_descriptor he_ccw = prev( opposite(he, pmesh()) , pmesh() );
         v2 = source(he_ccw, pmesh());
       }
       return ( CotangentValue::operator()(v0, v2, v1)/2.0 );
     }
     else
     {
       halfedge_descriptor he_cw = opposite( next(he, pmesh()) , pmesh() );
       vertex_descriptor v2 = source(he_cw, pmesh());     
       halfedge_descriptor he_ccw = prev( opposite(he, pmesh()) , pmesh() );
       vertex_descriptor v3 = source(he_ccw, pmesh());

        return ( CotangentValue::operator()(v0, v2, v1)/2.0 + CotangentValue::operator()(v0, v3, v1)/2.0 );
     }
  }
};

// Single cotangent from -[Chao10] Simple Geometric Model for Elastic Deformation
template<class PolygonMesh,
         class CotangentValue = Cotangent_value_Meyer_impl<PolygonMesh> >
struct Single_cotangent_weight_impl : CotangentValue
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  // Returns the cotangent of the opposite angle of the edge
  // 0 for border edges (which does not have an opposite angle)
  template <class VertexPointMap>
  double operator()(halfedge_descriptor he, PolygonMesh& pmesh, const VertexPointMap& ppmap)
  {
    if(is_border(he, pmesh)) { return 0.0;}
     
    vertex_descriptor v0 = target(he, pmesh);
    vertex_descriptor v1 = source(he, pmesh);

    vertex_descriptor v_op = target(next(he, pmesh), pmesh);
     return CotangentValue::operator()(v0, v_op, v1, ppmap);
  }
};

template<class PolygonMesh
  , class VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type
  , class CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Single_cotangent_weight : CotangentValue
{
  Single_cotangent_weight()
  {}
public:
  Single_cotangent_weight(PolygonMesh& pmesh_, VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  {}

  PolygonMesh& pmesh()
  {
    return CotangentValue::pmesh();
  }

  VertexPointMap& ppmap()
  {
    return CotangentValue::ppmap();
  }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<PolygonMesh,vertex_point_t>::type Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3  Vector;

  // Returns the cotangent of the opposite angle of the edge
  // 0 for border edges (which does not have an opposite angle)
  double operator()(halfedge_descriptor he)
  {
    if(is_border(he, pmesh())) { return 0.0;}
     
    vertex_descriptor v0 = target(he, pmesh());
    vertex_descriptor v1 = source(he, pmesh());

    vertex_descriptor v_op = target(next(he, pmesh()), pmesh());
     return CotangentValue::operator()(v0, v_op, v1);
  }
};


template<class PolygonMesh
  , class VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type
  , class CotangentValue = Cotangent_value_Meyer<PolygonMesh, VertexPointMap>
>
class Cotangent_weight_with_triangle_area : CotangentValue
{
  typedef PolygonMesh                                         PM;
  typedef VertexPointMap                                      VPMap;
  typedef typename boost::property_traits<VPMap>::value_type  Point;

  typedef typename boost::graph_traits<PM>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::vertex_descriptor     vertex_descriptor;

  Cotangent_weight_with_triangle_area()
  {}
public:
  Cotangent_weight_with_triangle_area(PolygonMesh& pmesh_, VertexPointMap vpmap_)
    : CotangentValue(pmesh_, vpmap_)
  {}

  PolygonMesh& pmesh()
  {
    return CotangentValue::pmesh();
  }

  VertexPointMap& ppmap()
  {
    return CotangentValue::ppmap();
  }

  double operator()(halfedge_descriptor he)
  {
    vertex_descriptor v0 = target(he, pmesh());
    vertex_descriptor v1 = source(he, pmesh());

     // Only one triangle for border edges
    if (is_border_edge(he, pmesh()))
    {
      halfedge_descriptor he_cw = opposite( next(he, pmesh()) , pmesh() );
      vertex_descriptor v2 = source(he_cw, pmesh());
      if (is_border_edge(he_cw, pmesh()) )
      {
        halfedge_descriptor he_ccw = prev( opposite(he, pmesh()) , pmesh() );
        v2 = source(he_ccw, pmesh());
      }

      const Point& v0_p = get(ppmap(), v0);
      const Point& v1_p = get(ppmap(), v1);
      const Point& v2_p = get(ppmap(), v2);
      double area_t = to_double(CGAL::sqrt(squared_area(v0_p, v1_p, v2_p)));

      return ( CotangentValue::operator()(v0, v2, v1) / area_t );
    }
    else
    {
      halfedge_descriptor he_cw = opposite( next(he, pmesh()) , pmesh() );
      vertex_descriptor v2 = source(he_cw, pmesh());
      halfedge_descriptor he_ccw = prev( opposite(he, pmesh()) , pmesh() );
      vertex_descriptor v3 = source(he_ccw, pmesh());

      const Point& v0_p = get(ppmap(), v0);
      const Point& v1_p = get(ppmap(), v1);
      const Point& v2_p = get(ppmap(), v2);
      const Point& v3_p = get(ppmap(), v3);
      double area_t1 = to_double(CGAL::sqrt(squared_area(v0_p, v1_p, v2_p)));
      double area_t2 = to_double(CGAL::sqrt(squared_area(v0_p, v1_p, v3_p)));

      return ( CotangentValue::operator()(v0, v2, v1) / area_t1
             + CotangentValue::operator()(v0, v3, v1) / area_t2);
     }

    return 0.;
  }
};

// Mean value calculator described in -[Floater04] Mean Value Coordinates-
template<class PolygonMesh
       , class VertexPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
class Mean_value_weight
{
  //Mean_value_weight()
  //{}

  PolygonMesh& pmesh_;
  VertexPointMap vpmap_;

public:
  Mean_value_weight(PolygonMesh& pmesh_, VertexPointMap vpmap)
    : pmesh_(pmesh_), vpmap_(vpmap)
  {}

  PolygonMesh& pmesh()
  {
    return pmesh_;
  }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;

  typedef VertexPointMap Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3  Vector;

  // Returns the mean-value coordinate of specified halfedge_descriptor
  // Returns different value for different edge orientation (which is a normal behaivour according to formula)
  double operator()(halfedge_descriptor he)
  {
    vertex_descriptor v0 = target(he, pmesh());
    vertex_descriptor v1 = source(he, pmesh());
    Vector vec = v0->point() - v1->point();
    double norm = CGAL::sqrt( vec.squared_length() );

    // Only one triangle for border edges
    if ( is_border_edge(he, pmesh()) )
    {
      halfedge_descriptor he_cw = opposite( next(he, pmesh()) , pmesh() );
      vertex_descriptor v2 = source(he_cw, pmesh());
      if ( is_border_edge(he_cw, pmesh()) )
      {
        halfedge_descriptor he_ccw = prev( opposite(he, pmesh()) , pmesh() );
        v2 = source(he_ccw, pmesh());
      }

      return ( half_tan_value_2(v1, v0, v2)/norm);
    }
    else
    {
      halfedge_descriptor he_cw = opposite( next(he, pmesh()) , pmesh() );
      vertex_descriptor v2 = source(he_cw, pmesh());     
      halfedge_descriptor he_ccw = prev( opposite(he, pmesh()) , pmesh() );
      vertex_descriptor v3 = source(he_ccw, pmesh());

      return ( half_tan_value_2(v1, v0, v2)/norm + half_tan_value_2(v1, v0, v3)/norm);
    }
  }

private:
  // Returns the tangent value of half angle v0_v1_v2/2
  double half_tan_value(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    Vector vec0 = v1->point() - v2->point();
    Vector vec1 = v2->point() - v0->point();
    Vector vec2 = v0->point() - v1->point();
    double e0_square = vec0.squared_length();
    double e1_square = vec1.squared_length();
    double e2_square = vec2.squared_length();
    double e0 = CGAL::sqrt(e0_square);
    double e2 = CGAL::sqrt(e2_square);
    double cos_angle = ( e0_square + e2_square - e1_square ) / 2.0 / e0 / e2;
    cos_angle = (std::max)(-1.0, (std::min)(1.0, cos_angle)); // clamp into [-1, 1]
    double angle = acos(cos_angle);
    
    return ( tan(angle/2.0) );
  }

  // My deviation built on Meyer_02
  double half_tan_value_2(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    Vector a = v0->point() - v1->point();
    Vector b = v2->point() - v1->point();
    double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    double dot_aa = a.squared_length();
    double dot_bb = b.squared_length();
    double dot_aa_bb = dot_aa * dot_bb;

    double cos_rep = dot_ab;
    double sin_rep = CGAL::sqrt(dot_aa_bb  - dot_ab * dot_ab);
    double normalizer = CGAL::sqrt(dot_aa_bb); // |a|*|b|

    return (normalizer - cos_rep) / sin_rep; // formula from [Floater04] page 4
                                             // tan(Q/2) = (1 - cos(Q)) / sin(Q)
  }
};

template< class PolygonMesh, 
          class PrimaryWeight = Cotangent_weight<PolygonMesh>,
          class SecondaryWeight = Mean_value_weight<PolygonMesh> >
class Hybrid_weight : public PrimaryWeight, SecondaryWeight
{
  PrimaryWeight primary;
  SecondaryWeight secondary;

  Hybrid_weight()
  {}

public:
  Hybrid_weight(PolygonMesh& pmesh_)
    : primary(pmesh_), secondary(pmesh_)
  {}

  PolygonMesh& pmesh()
  {
    return primary.pmesh();
  }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;

  double operator()(halfedge_descriptor he)
  {
    double weight = primary(he);
    //if(weight < 0) { std::cout << "Negative weight" << std::endl; }
    return (weight >= 0) ? weight : secondary(he);
  }
};

// Trivial uniform weights (created for test purposes)
template<class PolygonMesh>
class Uniform_weight
{
public:
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;

  double operator()(halfedge_descriptor /*e*/)
  { return 1.0; }
};

////////////////////////////////////////////////////////////////////////////
//                              FAIRING                                   //
template<class PolygonMesh>
class Scale_dependent_weight_fairing
{
  PolygonMesh& pmesh_;
public:
  Scale_dependent_weight_fairing(PolygonMesh& pmesh_)
    : pmesh_(pmesh_)
  {}

  PolygonMesh& pmesh()
  {
    return pmesh_;
  }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<PolygonMesh,vertex_point_t>::type Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3  Vector;

  double w_i(vertex_descriptor /*v_i*/) { return 1.0; }

  double w_ij(halfedge_descriptor he)
  {
    Vector v = target(he, pmesh())->point() - source(he, pmesh())->point();
    double divider = CGAL::sqrt(v.squared_length());
    if(divider == 0.0) {
      CGAL_warning(!"Scale dependent weight - zero length edge.");
      return (std::numeric_limits<double>::max)();
    }
    return 1.0 / divider;
  }
};

template<class PolygonMesh
  , class VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type
>
class Cotangent_weight_with_voronoi_area_fairing
{
  typedef PolygonMesh PM;
  typedef VertexPointMap VPMap;
  Voronoi_area<PM, VPMap> voronoi_functor;
  Cotangent_weight<PM, VPMap, Cotangent_value_Meyer<PM, VPMap> > cotangent_functor;

public:
  Cotangent_weight_with_voronoi_area_fairing(PM& pmesh_)
    : voronoi_functor(pmesh_, get(CGAL::vertex_point, pmesh_))
    , cotangent_functor(pmesh_, get(CGAL::vertex_point, pmesh_))
  {}

  Cotangent_weight_with_voronoi_area_fairing(PM& pmesh_, VPMap vpmap_)
    : voronoi_functor(pmesh_, vpmap_)
    , cotangent_functor(pmesh_, vpmap_)
  {}

  PM& pmesh()
  {
    return voronoi_functor.pmesh();
  }

  typedef typename boost::graph_traits<PM>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;

  double w_i(vertex_descriptor v_i)
  {
    return 0.5 / voronoi_functor(v_i);
  }

  double w_ij(halfedge_descriptor he) {

    return cotangent_functor(he) * 2.0;
  }
};

template<class PolygonMesh>
class Uniform_weight_fairing
{
public:
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  Uniform_weight_fairing(PolygonMesh&)
  {}

  double w_ij(halfedge_descriptor /*e*/) { return 1.0; }

  double w_i(vertex_descriptor /*v_i*/) { return 1.0; }
};
////////////////////////////////////////////////////////////////////////////

}//namespace internal


}//namespace CGAL
/// @endcond
#endif //CGAL_PMP_WEIGHTS_H
