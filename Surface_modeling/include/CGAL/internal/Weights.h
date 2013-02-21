#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

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
// Returns the cotangent value of half angle v0 v1 v2 by dividing the triangle area
template<class Polyhedron>
class Cotangent_value_area_weighted
  : Cotangent_value<Polyhedron>
{
public:
  double operator()(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    return Cotangent_value<Polyhedron>::operator()(v0, v1, v2)
      / std::sqrt(squared_area(v0->point(), v1->point(), v2->point()));
  }
};

// Cotangent weight calculator 
// Cotangent_value:               as suggested by -[Sorkine07] ARAP Surface Modeling-
// Cotangent_value_area_weighted: as suggested by -[Mullen08] Spectral Conformal Parameterization-
template<class Polyhedron, class CotangentValue>
class Cotangent_weight
{
public:
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  typedef typename Polyhedron::Traits::Vector_3  Vector;
  typedef typename Polyhedron::Traits::Point_3   Point;

  Cotangent_weight(Polyhedron& polyhedron) : polyhedron(polyhedron)
  { }

  // Returns the cotangent weight of specified edge_descriptor
  double operator()(edge_descriptor e)
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
       return ( cotangent_value(v0, v2, v1)/2.0 );
     }
     else
     {
        edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
        vertex_descriptor v2 = boost::source(e_cw, polyhedron);     
        edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
        vertex_descriptor v3 = boost::source(e_ccw, polyhedron);

        return ( cotangent_value(v0, v2, v1)/2.0 + cotangent_value(v0, v3, v1)/2.0 );
     }
  }

private:
  Polyhedron& polyhedron;	
  CotangentValue cotangent_value;
};


// Mean value calculator described in -[Floater04] Mean Value Coordinates-
template<class Polyhedron>
class Mean_value_weight
{
public:
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  typedef typename Polyhedron::Traits::Vector_3  Vector;
  typedef typename Polyhedron::Traits::Point_3   Point;

  Mean_value_weight(Polyhedron& polyhedron) : polyhedron(polyhedron)
  { }
  // Returns the mean-value coordinate of specified edge_descriptor
  double operator()(edge_descriptor e)
  {
    vertex_descriptor v0 = boost::target(e, polyhedron);
    vertex_descriptor v1 = boost::source(e, polyhedron);
    Vector vec = v0->point() - v1->point();
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

      return ( half_tan_value(v1, v0, v2)/norm );
    }
    else
    {
      edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
      vertex_descriptor v2 = boost::source(e_cw, polyhedron);     
      edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
      vertex_descriptor v3 = boost::source(e_ccw, polyhedron);

      return ( half_tan_value(v1, v0, v2)/norm + half_tan_value(v1, v0, v3)/norm );
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
    double angle = acos(cos_angle);

    return ( tan(angle/2.0) );
  }

  Polyhedron& polyhedron;	
};