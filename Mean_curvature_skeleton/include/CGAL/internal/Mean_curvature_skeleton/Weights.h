#ifndef CGAL_SURFACE_MODELING_WEIGHTS_H
#define CGAL_SURFACE_MODELING_WEIGHTS_H
/// @cond CGAL_DOCUMENT_INTERNAL

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
namespace CGAL {
namespace internal {

struct Vector{
  double coords[3];
  template<class Point>
  Vector(const Point& b, const Point& a) {
    coords[0] = a[0] - b[0];
    coords[1] = a[1] - b[1];
    coords[2] = a[2] - b[2];
  }
  double& operator[](int i)       { return coords[i]; }
  double  operator[](int i) const { return coords[i]; }
  double squared_length() const {
    return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
  }
  double length() const {
    return sqrtf(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]);
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
template<class Polyhedron>
class Cotangent_value
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
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

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    Vector a(v1->point(), v0->point());
    Vector b(v1->point(), v2->point());
    double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    double dot_aa = a.squared_length();
    double dot_bb = b.squared_length();
    return dot_ab / std::sqrt( dot_aa * dot_bb - dot_ab * dot_ab );
  }
};

template<class Polyhedron>
class Cotangent_value_Meyer_secure
{
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    Vector a(v1->point(), v0->point());
    Vector b(v1->point(), v2->point());
    double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    double dot_aa = a.squared_length();
    double dot_bb = b.squared_length();
    double lb = -0.999, ub = 0.999;
    double cosine = dot_ab / sqrtf(dot_aa) / sqrtf(dot_bb);
    cosine = (cosine < lb) ? lb : cosine;
    cosine = (cosine > ub) ? ub : cosine;
    double sine = sqrtf(1.0 - cosine * cosine);
    return cosine / sine;
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
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  typedef typename Polyhedron::Traits::Point_3   Point;

  // Returns the cotangent weight of specified edge_descriptor
  // Edge orientation is trivial
  double operator()(edge_descriptor e, Polyhedron& polyhedron)
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
       return ( CotangentValue::operator()(v0, v2, v1)/2.0 );
     }
     else
     {
        edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
        vertex_descriptor v2 = boost::source(e_cw, polyhedron);     
        edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
        vertex_descriptor v3 = boost::source(e_ccw, polyhedron);

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
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  typedef typename Polyhedron::Traits::Point_3   Point;

  // Returns the cotangent of the opposite angle of the edge
  // 0 for border edges (which does not have an opposite angle)
  double operator()(edge_descriptor e, Polyhedron& polyhedron)
  {
     if(boost::get(CGAL::edge_is_border, polyhedron, e)) { return 0.0;}
     
     vertex_descriptor v0 = boost::target(e, polyhedron);
     vertex_descriptor v1 = boost::source(e, polyhedron);

     vertex_descriptor v_op = boost::target(CGAL::next_edge(e, polyhedron), polyhedron);
     return CotangentValue::operator()(v0, v_op, v1);
  }
};

// Mean value calculator described in -[Floater04] Mean Value Coordinates-
template<class Polyhedron>
class Mean_value_weight
{
public:
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  typedef typename Polyhedron::Traits::Point_3   Point;

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
