// Copyright (c) 2011-2013 GeometryFactory
// All rights reserved.
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
// $URL:$
// $Id:$
//
// Author(s)     : Yin Xu, Andreas Fabri and Ilker O. Yaz

#ifndef CGAL_SURFACE_MODELING_WEIGHTS_H
#define CGAL_SURFACE_MODELING_WEIGHTS_H
/// @cond CGAL_DOCUMENT_INTERNAL

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel/global_functions_3.h>

namespace CGAL {
namespace internal {
/////////////////////////////////////////////////////////////////////////////////////////
// Returns the cotangent value of half angle v0 v1 v2
template<class Polyhedron>
class Cotangent_value
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
  typedef typename Polyhedron::Traits::Vector_3  Vector;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    Vector vec0 = v1->point() - v2->point();
    Vector vec1 = v2->point() - v0->point();
    Vector vec2 = v0->point() - v1->point();
    double e0_square = vec0.squared_length();
    double e1_square = vec1.squared_length();
    double e2_square = vec2.squared_length();
    double e0 = std::sqrt(e0_square); 
    double e2 = std::sqrt(e2_square);
    double cos_angle = ( e0_square + e2_square - e1_square ) / 2.0 / e0 / e2;
    double sin_angle = std::sqrt(1-cos_angle*cos_angle);

    return (cos_angle/sin_angle);
  }
};

// Returns the cotangent value of half angle v0 v1 v2
// using formula in -[Meyer02] Discrete Differential-Geometry Operators for- page 19
// The potential problem with previous one (Cotangent_value) is that it does not produce symmetric results
// (i.e. for v0, v1, v2 and v2, v1, v0 returned cot weights can be slightly different)
// This one provides stable results.
template<class Polyhedron>
class Cotangent_value_Meyer
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
  typedef typename Polyhedron::Traits::Vector_3  Vector;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    Vector a = v0->point() - v1->point();
    Vector b = v2->point() - v1->point();

    
    double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    // rewritten for safer fp arithmetic
    //double dot_aa = a.squared_length();
    //double dot_bb = b.squared_length();
    //double divider = std::sqrt( dot_aa * dot_bb - dot_ab * dot_ab );

    Vector cross_ab = CGAL::cross_product(a, b);
    double divider = std::sqrt(cross_ab*cross_ab);

    if(divider == 0 /*|| divider != divider*/) 
    {
      CGAL::collinear(v0->point(), v1->point(), v2->point()) ? 
        CGAL_warning(!"Infinite Cotangent value with degenerate triangle!") :
        CGAL_warning(!"Infinite Cotangent value due to floating point arithmetic!");
      

      return dot_ab > 0 ? (std::numeric_limits<double>::max)() :
                         -(std::numeric_limits<double>::max)();
    }
    
    return dot_ab / divider;
  }
};

// Returns the cotangent value of half angle v0 v1 v2 by clamping between [1, 89] degrees
// as suggested by -[Friedel] Unconstrained Spherical Parameterization-
template<class Polyhedron, class CotangentValue = Cotangent_value_Meyer<Polyhedron> >
class Cotangent_value_clamped : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    const double cot_1 = 57.289962;
    const double cot_89 = 0.017455;
    double value = CotangentValue::operator()(v0, v1, v2);
    return (std::max)(cot_89, (std::min)(value, cot_1));
  }
};

template<class Polyhedron, class CotangentValue = Cotangent_value_Meyer<Polyhedron> >
class Cotangent_value_clamped_2 : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    const double cot_5 = 5.671282;
    const double cot_175 = -cot_5;
    double value = CotangentValue::operator()(v0, v1, v2);
    return (std::max)(cot_175, (std::min)(value, cot_5));
  }
};

template<class Polyhedron, class CotangentValue = Cotangent_value_Meyer<Polyhedron> >
class Cotangent_value_minimum_zero : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    double value = CotangentValue::operator()(v0, v1, v2);
    return (std::max)(0.0, value);
  }
};

template<class Polyhedron, 
         class CotangentValue = Cotangent_value_Meyer<Polyhedron> >
class Voronoi_area : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator in_edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
  typedef typename Polyhedron::Traits::Point_3   Point;

  double operator()(vertex_descriptor v0, Polyhedron& polyhedron) {
    //return 1.0;
    double voronoi_area = 0.0;
    in_edge_iterator e, e_end;
    for (boost::tie(e,e_end) = boost::in_edges(v0, polyhedron); e != e_end; e++)
    {
      halfedge_descriptor he = halfedge(*e,polyhedron);
      if( is_border(he,polyhedron) ) { continue; }

      vertex_descriptor v1 = source(he, polyhedron);
      vertex_descriptor v_op = target(next(he, polyhedron), polyhedron);

      const Point& v0_p = v0->point();
      const Point& v1_p = v1->point();
      const Point& v_op_p = v_op->point();

      // (?) check if there is a better way to predicate triangle is obtuse or not
      CGAL::Angle angle0 = CGAL::angle(v1_p, v0_p, v_op_p);
      CGAL::Angle angle1 = CGAL::angle(v_op_p, v1_p, v0_p);
      CGAL::Angle angle_op = CGAL::angle(v0_p, v_op_p, v1_p);

      bool obtuse = (angle0 == CGAL::OBTUSE) || (angle1 == CGAL::OBTUSE) || (angle_op == CGAL::OBTUSE);

      if(!obtuse) {
        double cot_v1 = CotangentValue::operator()(v_op, v1, v0);
        double cot_v_op = CotangentValue::operator()(v0, v_op, v1);

        double term1 = cot_v1   * (v_op_p - v0_p).squared_length();
        double term2 = cot_v_op * (v1_p  - v0_p).squared_length();
        voronoi_area += (1.0 / 8.0) * (term1 + term2);
      }
      else {
        double area_t = std::sqrt(squared_area(v0_p, v1_p, v_op_p));
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
template<class Polyhedron, 
         class CotangentValue = Cotangent_value_Meyer<Polyhedron> >
class Cotangent_value_area_weighted : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    return CotangentValue::operator()(v0, v1, v2)
      / std::sqrt(squared_area(v0->point(), v1->point(), v2->point()));
  }
};
/////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////// Edge Weight Calculators ///////////////////////////////////
// Cotangent weight calculator 
// Cotangent_value:               as suggested by -[Sorkine07] ARAP Surface Modeling-
// Cotangent_value_area_weighted: as suggested by -[Mullen08] Spectral Conformal Parameterization-
template<class Polyhedron, 
         class CotangentValue = Cotangent_value_minimum_zero<Polyhedron> >
class Cotangent_weight : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  typedef typename Polyhedron::Traits::Vector_3  Vector;
  typedef typename Polyhedron::Traits::Point_3   Point;

  // Returns the cotangent weight of specified halfedge_descriptor
  // Edge orientation is trivial
  double operator()(halfedge_descriptor he, Polyhedron& polyhedron)
  {
     vertex_descriptor v0 = target(he, polyhedron);
     vertex_descriptor v1 = source(he, polyhedron);
     // Only one triangle for border edges
     if (is_border_edge(he, polyhedron))
     {
       
       halfedge_descriptor he_cw = opposite( next(he, polyhedron) , polyhedron );
       vertex_descriptor v2 = source(he_cw, polyhedron);
       if (is_border_edge(he_cw, polyhedron) )
       {
          halfedge_descriptor he_ccw = prev( opposite(he, polyhedron) , polyhedron );
          v2 = source(he_ccw, polyhedron);
       }
       return ( CotangentValue::operator()(v0, v2, v1)/2.0 );
     }
     else
     {
        halfedge_descriptor he_cw = opposite( next(he, polyhedron) , polyhedron );
        vertex_descriptor v2 = source(he_cw, polyhedron);     
        halfedge_descriptor he_ccw = prev( opposite(he, polyhedron) , polyhedron );
        vertex_descriptor v3 = source(he_ccw, polyhedron);

        return ( CotangentValue::operator()(v0, v2, v1)/2.0 + CotangentValue::operator()(v0, v3, v1)/2.0 );
     }
  }
};

// Single cotangent from -[Chao10] Simple Geometric Model for Elastic Deformation
template<class Polyhedron, 
         class CotangentValue = Cotangent_value_Meyer<Polyhedron> >
class Single_cotangent_weight : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  typedef typename Polyhedron::Traits::Vector_3  Vector;
  typedef typename Polyhedron::Traits::Point_3   Point;

  // Returns the cotangent of the opposite angle of the edge
  // 0 for border edges (which does not have an opposite angle)
  double operator()(halfedge_descriptor he, Polyhedron& polyhedron)
  {
     if(is_border(he, polyhedron)) { return 0.0;}
     
     vertex_descriptor v0 = target(he, polyhedron);
     vertex_descriptor v1 = source(he, polyhedron);

     vertex_descriptor v_op = target(CGAL::next_edge(he, polyhedron), polyhedron);
     return CotangentValue::operator()(v0, v_op, v1);
  }
};

// Mean value calculator described in -[Floater04] Mean Value Coordinates-
template<class Polyhedron>
class Mean_value_weight
{
public:
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  typedef typename Polyhedron::Traits::Vector_3  Vector;
  typedef typename Polyhedron::Traits::Point_3   Point;

  // Returns the mean-value coordinate of specified halfedge_descriptor
  // Returns different value for different edge orientation (which is a normal behaivour according to formula)
  double operator()(halfedge_descriptor he, Polyhedron& polyhedron)
  {
    vertex_descriptor v0 = target(he, polyhedron);
    vertex_descriptor v1 = source(he, polyhedron);
    Vector vec = v0->point() - v1->point();
    double norm = std::sqrt( vec.squared_length() );

    // Only one triangle for border edges
    if ( is_border_edge(he, polyhedron) )
    {
      halfedge_descriptor he_cw = opposite( next(he, polyhedron) , polyhedron );
      vertex_descriptor v2 = source(he_cw, polyhedron);
      if ( is_border_edge(he_cw, polyhedron) )
      {
        halfedge_descriptor he_ccw = prev( opposite(he, polyhedron) , polyhedron );
        v2 = source(he_ccw, polyhedron);
      }

      return ( half_tan_value_2(v1, v0, v2)/norm);
    }
    else
    {
      halfedge_descriptor he_cw = opposite( next(he, polyhedron) , polyhedron );
      vertex_descriptor v2 = source(he_cw, polyhedron);     
      halfedge_descriptor he_ccw = prev( opposite(he, polyhedron) , polyhedron );
      vertex_descriptor v3 = source(he_ccw, polyhedron);

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
    double e0 = std::sqrt(e0_square); 
    double e2 = std::sqrt(e2_square);
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
    double sin_rep = std::sqrt(dot_aa_bb  - dot_ab * dot_ab);
    double normalizer = std::sqrt(dot_aa_bb); // |a|*|b|

    return (normalizer - cos_rep) / sin_rep; // formula from [Floater04] page 4
                                             // tan(Q/2) = (1 - cos(Q)) / sin(Q)
  }
};

template< class Polyhedron, 
          class PrimaryWeight = Cotangent_weight<Polyhedron>,
          class SecondaryWeight = Mean_value_weight<Polyhedron> >
class Hybrid_weight : public PrimaryWeight, SecondaryWeight
{
public:
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor   halfedge_descriptor;

  double operator()(halfedge_descriptor he, Polyhedron& polyhedron)
  {
    double weight = PrimaryWeight::operator()(he, polyhedron);
    //if(weight < 0) { std::cout << "Negative weight" << std::endl; }
    return (weight >= 0) ? weight : SecondaryWeight::operator()(he, polyhedron);
  }
};

// Trivial uniform weights (created for test purposes)
template<class Polyhedron>
class Uniform_weight
{
public:
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor   halfedge_descriptor;

  double operator()(halfedge_descriptor /*e*/, Polyhedron& /*polyhedron*/)
  { return 1.0; }
};

////////////////////////////////////////////////////////////////////////////
//                              FAIRING                                   //
template<class Polyhedron>
class Scale_dependent_weight_fairing
{
public:
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
  typedef typename Polyhedron::Traits::Vector_3  Vector;

  double w_i(vertex_descriptor /*v_i*/, Polyhedron& /*polyhedron*/) { return 1.0; }

  double w_ij(halfedge_descriptor he, Polyhedron& polyhedron)
  {
    Vector v = target(he, polyhedron)->point() - source(he, polyhedron)->point();
    double divider = std::sqrt(v.squared_length());
    if(divider == 0.0) {
      CGAL_warning(!"Scale dependent weight - zero length edge.");
      return (std::numeric_limits<double>::max)();
    }
    return 1.0 / divider;
  }
};

template<class Polyhedron>
class Cotangent_weight_with_voronoi_area_fairing {
public:
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  double w_i(vertex_descriptor v_i, Polyhedron& polyhedron) {
    Voronoi_area<Polyhedron> voronoi_functor;
    return 0.5 / voronoi_functor(v_i, polyhedron);
  }

  double w_ij(halfedge_descriptor he, Polyhedron& polyhedron) {
    Cotangent_weight<Polyhedron, Cotangent_value_Meyer<Polyhedron> > cotangent_functor;
    return cotangent_functor(he, polyhedron) * 2.0;
  }
};

template<class Polyhedron>
class Uniform_weight_fairing
{
public:
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  double w_ij(halfedge_descriptor /*e*/, Polyhedron& /*polyhedron*/) { return 1.0; }

  double w_i(vertex_descriptor /*v_i*/, Polyhedron& /*polyhedron*/) { return 1.0; }
};
////////////////////////////////////////////////////////////////////////////

}//namespace internal
/// @endcond



}//namespace CGAL
#endif //CGAL_SURFACE_MODELING_WEIGHTS_H
