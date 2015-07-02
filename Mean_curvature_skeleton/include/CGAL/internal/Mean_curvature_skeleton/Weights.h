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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_MCS_WEIGHTS_H
#define CGAL_MCS_WEIGHTS_H
/// @cond CGAL_DOCUMENT_INTERNAL

#include <boost/graph/graph_traits.hpp>
#include <cmath>

namespace CGAL {
namespace internal {

struct Vector{
  double coords[3];
  template<class Traits>
  Vector(const typename Traits::Point_3& b, const typename Traits::Point_3& a, const Traits& traits) {
    coords[0] = traits.compute_x_3_object()(a) - traits.compute_x_3_object()(b);
    coords[1] = traits.compute_y_3_object()(a) - traits.compute_y_3_object()(b);
    coords[2] = traits.compute_z_3_object()(a) - traits.compute_z_3_object()(b);
  }
  double& operator[](int i)       { return coords[i]; }
  double  operator[](int i) const { return coords[i]; }
  double squared_length() const {
    return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
  }
  double length() const {
    return std::sqrt(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]);
  }
  bool normalize() {
    double len = length();
    if (len < 1e-10)
    {
      return false;
    }
    coords[0] /= len;
    coords[1] /= len;
    coords[2] /= len;
    return true;
  }
  double dot(const Vector& b) {
    return coords[0] * b.coords[0] + coords[1] * b.coords[1] + coords[2] * b.coords[2];
  }
};

/////////////////////////////////////////////////////////////////////////////////////////
// Returns the cotangent value of half angle v0 v1 v2
template<class TriangleMesh>
class Cotangent_value
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  template<class Traits>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, const Traits& traits)
  {
    Vector vec0(v2->point(), v1->point(), traits);
    Vector vec1(v0->point(), v2->point(), traits);
    Vector vec2(v1->point(), v0->point(), traits);
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
template<class TriangleMesh>
class Cotangent_value_Meyer
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  template<class Traits>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, const Traits& traits)
  {
    Vector a(v1->point(), v0->point(), traits);
    Vector b(v1->point(), v2->point(), traits);
    double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    double dot_aa = a.squared_length();
    double dot_bb = b.squared_length();
    return dot_ab / std::sqrt( dot_aa * dot_bb - dot_ab * dot_ab );
  }
};

template<class TriangleMesh>
class Cotangent_value_Meyer_secure
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  template<class Traits>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, const Traits& traits)
  {
    Vector a(v1->point(), v0->point(), traits);
    Vector b(v1->point(), v2->point(), traits);
    double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    double dot_aa = a.squared_length();
    double dot_bb = b.squared_length();
    double lb = -0.999, ub = 0.999;
    double cosine = dot_ab / std::sqrt(dot_aa) / std::sqrt(dot_bb);
    cosine = (cosine < lb) ? lb : cosine;
    cosine = (cosine > ub) ? ub : cosine;
    double sine = std::sqrt(1.0 - cosine * cosine);
    return cosine / sine;
  }
};

// Returns the cotangent value of half angle v0 v1 v2 by clamping between [1, 89] degrees
// as suggested by -[Friedel] Unconstrained Spherical Parameterization-
template<class TriangleMesh, class CotangentValue = Cotangent_value_Meyer<TriangleMesh> >
class Cotangent_value_clamped : CotangentValue
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  template<class Traits>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, const Traits& traits)
  {
    const double cot_1 = 57.289962;
    const double cot_89 = 0.017455;
    double value = CotangentValue::operator()(v0, v1, v2, traits);
    return (std::max)(cot_89, (std::min)(value, cot_1));
  }
};

template<class TriangleMesh, class CotangentValue = Cotangent_value_Meyer<TriangleMesh> >
class Cotangent_value_minimum_zero : CotangentValue
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  template <class Traits>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, const Traits& traits)
  {
    double value = CotangentValue::operator()(v0, v1, v2, traits);
    return (std::max)(0.0, value);
  }
};

// Returns the cotangent value of half angle v0 v1 v2 by dividing the triangle area
// as suggested by -[Mullen08] Spectral Conformal Parameterization-
template<class TriangleMesh,
         class CotangentValue = Cotangent_value_Meyer<TriangleMesh> >
class Cotangent_value_area_weighted : CotangentValue
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  template <class Traits>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, const Traits& traits)
  {
    return CotangentValue::operator()(v0, v1, v2, traits)
      / traits.compute_area_3_object()(v0->point(), v1->point(), v2->point());
  }
};
/////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////// Edge Weight Calculators ///////////////////////////////////
// Cotangent weight calculator
// Cotangent_value:               as suggested by -[Sorkine07] ARAP Surface Modeling-
// Cotangent_value_area_weighted: as suggested by -[Mullen08] Spectral Conformal Parameterization-
template<class TriangleMesh,
         class CotangentValue = Cotangent_value_minimum_zero<TriangleMesh> >
class Cotangent_weight : CotangentValue
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  //  typedef typename TriangleMesh::Traits::Point_3   Point;

  // Returns the cotangent weight of specified halfedge_descriptor
  // Edge orientation is trivial
  template <class Traits>
  double operator()(halfedge_descriptor e, TriangleMesh& hg, const Traits& traits)
  {
     vertex_descriptor v0 = target(e, hg);
     vertex_descriptor v1 = source(e, hg);
     // Only one triangle for border edges
     if ( is_border(e, hg) ||
          is_border( opposite(e, hg), hg) )
     {
       halfedge_descriptor e_cw = opposite(next(e, hg), hg);
       vertex_descriptor v2 = source(e_cw, hg);
       if ( is_border(e_cw, hg) ||
            is_border( opposite(e_cw, hg), hg) )
       {
          halfedge_descriptor e_ccw = prev(opposite(e, hg), hg);
          v2 = source(e_ccw, hg);
       }
       return ( CotangentValue::operator()(v0, v2, v1, traits)/2.0 );
     }
     else
     {
        halfedge_descriptor e_cw = opposite(next(e, hg), hg);
        vertex_descriptor v2 = source(e_cw, hg);
        halfedge_descriptor e_ccw = prev(opposite(e, hg), hg);
        vertex_descriptor v3 = source(e_ccw, hg);

        return ( CotangentValue::operator()(v0, v2, v1, traits)/2.0 + CotangentValue::operator()(v0, v3, v1, traits)/2.0 );
     }
  }
};

// Single cotangent from -[Chao10] Simple Geometric Model for Elastic Deformation
template<class TriangleMesh,
         class CotangentValue = Cotangent_value_Meyer<TriangleMesh> >
class Single_cotangent_weight : CotangentValue
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  //  typedef typename TriangleMesh::Traits::Point_3   Point;

  // Returns the cotangent of the opposite angle of the edge
  // 0 for border edges (which does not have an opposite angle)
  template <class Traits>
  double operator()(halfedge_descriptor e, TriangleMesh& hg, const Traits& traits)
  {
     if( is_border(e, hg) ) { return 0.0;}

     vertex_descriptor v0 = target(e, hg);
     vertex_descriptor v1 = source(e, hg);

     vertex_descriptor v_op = target(next(e, hg), hg);
     return CotangentValue::operator()(v0, v_op, v1, traits);
  }
};

// Mean value calculator described in -[Floater04] Mean Value Coordinates-
template<class TriangleMesh>
class Mean_value_weight
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  typedef typename TriangleMesh::Traits::Point_3   Point;

  // Returns the mean-value coordinate of specified halfedge_descriptor
  // Returns different value for different edge orientation (which is a normal behaivour according to formula)
  template <class Traits>
  double operator()(halfedge_descriptor e, TriangleMesh& hg, const Traits& traits)
  {
    vertex_descriptor v0 = target(e, hg);
    vertex_descriptor v1 = source(e, hg);
    Vector vec(v1->point(), v0->point(), traits);
    double norm = std::sqrt( vec.squared_length() );

    // Only one triangle for border edges
    if ( is_border(e, hg) ||
         is_border( opposite(e, hg), hg) )
    {
      halfedge_descriptor e_cw = opposite(next(e, hg), hg);
      vertex_descriptor v2 = source(e_cw, hg);
      if ( is_border(e_cw, hg) ||
           is_border( opposite(e_cw, hg), hg) )
      {
        halfedge_descriptor e_ccw = prev(opposite(e, hg), hg);
        v2 = source(e_ccw, hg);
      }

      return ( half_tan_value_2(v1, v0, v2, traits)/norm);
    }
    else
    {
      halfedge_descriptor e_cw = opposite(next(e, hg), hg);
      vertex_descriptor v2 = source(e_cw, hg);
      halfedge_descriptor e_ccw = prev(opposite(e, hg), hg);
      vertex_descriptor v3 = source(e_ccw, hg);

      return ( half_tan_value_2(v1, v0, v2, traits)/norm + half_tan_value_2(v1, v0, v3, traits)/norm);
    }
  }

private:
  // Returns the tangent value of half angle v0_v1_v2/2
  template <class Traits>
  double half_tan_value(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, const Traits& traits)
  {
    Vector vec0(v2->point(), v1->point(), traits);
    Vector vec1(v0->point(), v2->point(), traits);
    Vector vec2(v1->point(), v0->point(), traits);
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
  template <class Traits>
  double half_tan_value_2(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, const Traits& traits)
  {
    Vector a(v1->point(), v0->point(), traits);
    Vector b(v1->point(), v2->point(), traits);
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

template< class TriangleMesh,
          class PrimaryWeight = Cotangent_weight<TriangleMesh>,
          class SecondaryWeight = Mean_value_weight<TriangleMesh> >
class Hybrid_weight : public PrimaryWeight, SecondaryWeight
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;

  template<class Traits>
  double operator()(halfedge_descriptor e, TriangleMesh& hg, const Traits& traits)
  {
    double weight = PrimaryWeight::operator()(e, hg);
    //if(weight < 0) { std::cout << "Negative weight" << std::endl; }
    return (weight >= 0) ? weight : SecondaryWeight::operator()(e, hg, traits);
  }
};

// Trivial uniform weights (created for test purposes)
template<class TriangleMesh>
class Uniform_weight
{
public:
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;

  template<class Traits>
  double operator()(halfedge_descriptor /*e*/, TriangleMesh& /*hg*/,const Traits&)
  { return 1.0; }
};



}//namespace internal
/// @endcond
}//namespace CGAL
#endif //CGAL_MCS_WEIGHTS_H
