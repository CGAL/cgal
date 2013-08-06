#ifndef CGAL_SURFACE_MODELING_WEIGHTS_H
#define CGAL_SURFACE_MODELING_WEIGHTS_H
/// @cond CGAL_DOCUMENT_INTERNAL

#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<double>::Vector_3 Vector;

namespace CGAL {
namespace internal {

template<class Point>
Vector to_vector(const Point& b, const Point& a) {
  return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

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

  template<class VertexPointMap>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, VertexPointMap vpm)
  {
    const Vector& a = to_vector(get(vpm, v1), get(vpm, v0));
    const Vector& b = to_vector(get(vpm, v1), get(vpm, v2));

    double dot_ab = a*b;
    Vector cross_ab = CGAL::cross_product(a, b);
    double divider = std::sqrt(cross_ab*cross_ab);

    if(divider == 0 /*|| divider != divider*/) 
    {
      this->collinear(get(vpm, v0), get(vpm, v1), get(vpm, v2)) ? 
        CGAL_warning(!"Infinite Cotangent value with degenerate triangle!") :
        CGAL_warning(!"Infinite Cotangent value due to floating point arithmetic!");

      return dot_ab > 0 ? (std::numeric_limits<double>::max)() :
        -(std::numeric_limits<double>::max)();
    }

    return dot_ab / divider;
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  // WARNING: this two functions are just used when cotangent weight turns out to be +-inf,
  //          just for raising a proper warning message (i.e nothing functional)
  template<class Point>
  bool collinear(const Point&, const Point&, const Point&) {
    return true; 
  }
  template<class Kernel>
  bool collinear(const CGAL::Point_3<Kernel>& a, const CGAL::Point_3<Kernel>& b, const CGAL::Point_3<Kernel>& c) {
    return CGAL::collinear(a, b, c);
  }
  ///////////////////////////////////////////////////////////////////////////////////////
  
};

// Returns the cotangent value of half angle v0 v1 v2 by clamping between [1, 89] degrees
// as suggested by -[Friedel] Unconstrained Spherical Parameterization-
template<class Polyhedron, class CotangentValue = Cotangent_value_Meyer<Polyhedron> >
class Cotangent_value_clamped : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  template<class VertexPointMap>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, VertexPointMap vpm)
  {
    const double cot_1 = 57.289962;
    const double cot_89 = 0.017455;
    double value = CotangentValue::operator()(v0, v1, v2, vpm);
    return (std::max)(cot_89, (std::min)(value, cot_1));
  }
};

template<class Polyhedron, class CotangentValue = Cotangent_value_Meyer<Polyhedron> >
class Cotangent_value_minimum_zero : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  template<class VertexPointMap>
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2, VertexPointMap vpm)
  {
    double value = CotangentValue::operator()(v0, v1, v2, vpm);
    return (std::max)(0.0, value);
  }
};


///////////////////////////// Edge Weight Calculators ///////////////////////////////////
// Cotangent weight calculator 
// Cotangent_value:               as suggested by -[Sorkine07] ARAP Surface Modeling-
// Cotangent_value_area_weighted: as suggested by -[Mullen08] Spectral Conformal Parameterization-
template<class Polyhedron, 
         class CotangentValue = Cotangent_value_minimum_zero<Polyhedron> >
class Cotangent_weight : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  // Returns the cotangent weight of specified edge_descriptor
  // Edge orientation is trivial
  template<class VertexPointMap>
  double operator()(edge_descriptor e, Polyhedron& polyhedron, VertexPointMap vpm)
  {
     vertex_descriptor v0 = boost::target(e, polyhedron);
     vertex_descriptor v1 = boost::source(e, polyhedron);
     // Only one triangle for border edges
     if (boost::get(CGAL::edge_is_border, polyhedron, e) ||
         boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(e, polyhedron)))
     {       
       edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
       vertex_descriptor v2 = boost::source(e_cw, polyhedron);
       if (boost::get(CGAL::edge_is_border, polyhedron, e_cw) ||
           boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(e_cw, polyhedron)) )
       {
          edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
          v2 = boost::source(e_ccw, polyhedron);
       }
       return ( CotangentValue::operator()(v0, v2, v1, vpm)/2.0 );
     }
     else
     {
        edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
        vertex_descriptor v2 = boost::source(e_cw, polyhedron);     
        edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
        vertex_descriptor v3 = boost::source(e_ccw, polyhedron);

        return ( CotangentValue::operator()(v0, v2, v1, vpm)/2.0 + CotangentValue::operator()(v0, v3, v1, vpm)/2.0 );
     }
  }
};

// Single cotangent from -[Chao10] Simple Geometric Model for Elastic Deformation
template<class Polyhedron, 
         class CotangentValue = Cotangent_value_Meyer<Polyhedron> >
class Single_cotangent_weight : CotangentValue
{
public:
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  // Returns the cotangent of the opposite angle of the edge
  // 0 for border edges (which does not have an opposite angle)
  template<class VertexPointMap>
  double operator()(edge_descriptor e, Polyhedron& polyhedron, VertexPointMap vpm)
  {
     if(boost::get(CGAL::edge_is_border, polyhedron, e)) { return 0.0;}
     
     vertex_descriptor v0 = boost::target(e, polyhedron);
     vertex_descriptor v1 = boost::source(e, polyhedron);

     vertex_descriptor v_op = boost::target(CGAL::next_edge(e, polyhedron), polyhedron);
     return CotangentValue::operator()(v0, v_op, v1, vpm);
  }
};

// Mean value calculator described in -[Floater04] Mean Value Coordinates-
// WARNING: Need to be updated to use point pmap
template<class Polyhedron>
class Mean_value_weight
{
public:
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  // Returns the mean-value coordinate of specified edge_descriptor
  // Returns different value for different edge orientation (which is a normal behaivour according to formula)
  double operator()(edge_descriptor e, Polyhedron& polyhedron)
  {
    vertex_descriptor v0 = boost::target(e, polyhedron);
    vertex_descriptor v1 = boost::source(e, polyhedron);
    Vector vec(v1->point(), v0->point());
    double norm = std::sqrt( vec.squared_length() );

    // Only one triangle for border edges
    if (boost::get(CGAL::edge_is_border, polyhedron, e) ||
        boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(e, polyhedron)))
    {
      edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
      vertex_descriptor v2 = boost::source(e_cw, polyhedron);
      if (boost::get(CGAL::edge_is_border, polyhedron, e_cw) || 
          boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(e_cw, polyhedron)) )
      {
        edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
        v2 = boost::source(e_ccw, polyhedron);
      }

      return ( half_tan_value_2(v1, v0, v2)/norm);
    }
    else
    {
      edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
      vertex_descriptor v2 = boost::source(e_cw, polyhedron);     
      edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
      vertex_descriptor v3 = boost::source(e_ccw, polyhedron);

      return ( half_tan_value_2(v1, v0, v2)/norm + half_tan_value_2(v1, v0, v3)/norm);
    }
  }

private:
  // Returns the tangent value of half angle v0_v1_v2/2
  double half_tan_value(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    Vector vec0(v2->point(), v1->point());
    Vector vec1(v0->point(), v2->point());
    Vector vec2(v1->point(), v0->point());
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
    Vector a(v1->point(), v0->point());
    Vector b(v1->point(), v2->point());
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
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;

  double operator()(edge_descriptor e, Polyhedron& polyhedron)
  {
    double weight = PrimaryWeight::operator()(e, polyhedron);
    //if(weight < 0) { std::cout << "Negative weight" << std::endl; }
    return (weight >= 0) ? weight : SecondaryWeight::operator()(e, polyhedron);
  }
};

// Trivial uniform weights (created for test purposes)
template<class Polyhedron>
class Uniform_weight
{
public:
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;

  double operator()(edge_descriptor /*e*/, Polyhedron& /*polyhedron*/)
  { return 1.0; }
};



}//namespace internal
/// @endcond
}//namespace CGAL
#endif //CGAL_SURFACE_MODELING_WEIGHTS_H