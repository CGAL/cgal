// Copyright (c) 2018 GeometryFactory (France).
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
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_MESH_SMOOTHING_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_MESH_SMOOTHING_IMPL_H

#include <math.h>
#include <utility>
#include <iterator>

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/property_map.h>
#include <CGAL/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>


namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename PolygonMesh, typename VertexPointMap, typename VertexConstraintMap, typename GeomTraits>
class Compatible_remesher
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_iterator vertex_iterator;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Segment_3 Segment;
  typedef typename GeomTraits::Triangle_3 Triangle;

  typedef std::vector<Triangle> Triangle_list;
  typedef std::pair<halfedge_descriptor, halfedge_descriptor> he_pair;
  typedef std::map<halfedge_descriptor, he_pair> Edges_around_map;

  typedef CGAL::AABB_triangle_primitive<GeomTraits, typename Triangle_list::iterator> AABB_Primitive;
  typedef CGAL::AABB_traits<GeomTraits, AABB_Primitive> AABB_Traits;
  typedef CGAL::AABB_tree<AABB_Traits> Tree;

public:
  Compatible_remesher(PolygonMesh& pmesh, VertexPointMap& vpmap, VertexConstraintMap& vcmap) :
    mesh_(pmesh), vpmap_(vpmap), vcmap_(vcmap)
  {}

  ~Compatible_remesher()
  {
    delete tree_ptr_;
  }

  template<typename FaceRange>
  void init_smoothing(const FaceRange& face_range)
  {
    check_vertex_range(face_range);

    BOOST_FOREACH(face_descriptor f, face_range)
      input_triangles_.push_back(triangle(f));

    tree_ptr_ = new Tree(input_triangles_.begin(), input_triangles_.end());
    tree_ptr_->accelerate_distance_queries();
  }

  std::size_t remove_degenerate_faces()
  {
    std::size_t nb_removed_faces = 0;
    nb_removed_faces = CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh_);
    return nb_removed_faces;
  }

  void angle_relaxation()
  {
    std::map<vertex_descriptor, Point> barycenters;
    std::map<vertex_descriptor, Vector> n_map;


    BOOST_FOREACH(vertex_descriptor v, vrange_)
    {
      if(!is_border(v, mesh_) && !is_constrained(v))
      {
        // compute normal to v
        Vector vn = compute_vertex_normal(v, mesh_,
                         Polygon_mesh_processing::parameters::vertex_point_map(vpmap_)
                         .geom_traits(traits_));
        n_map[v] = vn;

        Edges_around_map he_map;
        typename Edges_around_map::iterator it;
        BOOST_FOREACH(halfedge_descriptor hi, halfedges_around_source(v, mesh_))
          he_map[hi] = he_pair( next(hi, mesh_), prev(opposite(hi, mesh_), mesh_) );

        // measure initial angles
        measure_angles(he_map);

        // calculate movement
        Vector move = calc_move(he_map);

        barycenters[v] = (get(vpmap_, v) + move) ;

      } // not on border
    } // for each v

    // compute locations on tangent plane
    typedef typename std::map<vertex_descriptor, Point>::value_type VP;
    std::map<vertex_descriptor, Point> new_locations;
    BOOST_FOREACH(const VP& vp, barycenters)
    {
      Point p = get(vpmap_, vp.first);
      Point q = vp.second;
      Vector n = n_map[vp.first];

      new_locations[vp.first] = q + ( n * Vector(q, p) ) * n ;
    }

    // update location
    std::size_t moved_points = 0;
    BOOST_FOREACH(const VP& vp, new_locations)
    {
      // iff movement impoves all angles
      if(does_it_impove(vp.first, vp.second))
      {
        moved_points++;
        put(vpmap_, vp.first, vp.second);
      }
    }

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout<<"moved: "<< moved_points <<" points based on angle."<<std::endl;
    std::cout<<"not imporved min angle: "<< vrange_.size() - moved_points <<" times."<<std::endl;
#endif
  }

  void area_relaxation(const double& precision)
  {
    std::size_t moved_points = 0;
    BOOST_FOREACH(vertex_descriptor v, vrange_)
    {
       if(!is_border(v, mesh_) && !is_constrained(v))
       {
         if (gradient_descent(v, precision))
           moved_points++;
       }
    }

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout<<"moved : "<<moved_points<<" points based on area."<<std::endl;
    std::cout<<"non convex energy found: "<<vrange_.size() - moved_points<<" times."<<std::endl;
#endif
  }

  void project_to_surface()
  {
    BOOST_FOREACH( vertex_descriptor v, vertices(mesh_))
    {
      if(!is_border(v, mesh_) && !is_constrained(v))
      {
        Point p_query = get(vpmap_, v);
        Point projected = tree_ptr_->closest_point(p_query);
        put(vpmap_, v, projected);
      }
    }
  }

private:
  // helper functions
  // ----------------
  Triangle triangle(face_descriptor f) const
  {
    halfedge_descriptor h = halfedge(f, mesh_);
    vertex_descriptor v1 = target(h, mesh_);
    vertex_descriptor v2 = target(next(h, mesh_), mesh_);
    vertex_descriptor v3 = target(next(next(h, mesh_), mesh_), mesh_);
    return Triangle(get(vpmap_, v1), get(vpmap_, v2), get(vpmap_, v3));
  }

  double sqlength(const vertex_descriptor& v1, const vertex_descriptor& v2) const
  {
    return to_double(CGAL::squared_distance(get(vpmap_, v1), get(vpmap_, v2)));
  }

  double sqlength(const halfedge_descriptor& h) const
  {
   vertex_descriptor v1 = target(h, mesh_);
   vertex_descriptor v2 = source(h, mesh_);
   return sqlength(v1, v2);
  }

  double sqlength(const edge_descriptor& e) const
  {
   return sqlength(halfedge(e, mesh_));
  }

  // angle bisecting functions
  // -------------------------
  Vector calc_move(const Edges_around_map& he_map)
  {
    Vector move = CGAL::NULL_VECTOR;
    double weights_sum = 0;
    typename Edges_around_map::const_iterator it;
    for(it = he_map.begin(); it != he_map.end(); ++it)
    {
      halfedge_descriptor main_he = it->first;
      he_pair incident_pair = it->second;

      // avoid zero angles
      Point pt = get(vpmap_, source(incident_pair.first, mesh_));
      Point p1 = get(vpmap_, target(incident_pair.first, mesh_));
      Point p2 = get(vpmap_, source(incident_pair.second, mesh_));
      CGAL_assertion(target(incident_pair.second, mesh_) == source(incident_pair.first, mesh_));
      Vector e1(pt, p1);
      Vector e2(pt, p2);
      double angle = get_angle(e1, e2);
      if(angle < 1e-5)
        continue;

      // rotate
      Vector rotated_edge = rotate_edge(main_he, incident_pair);

      // small angles carry more weight
      double weight = 1.0 / (angle*angle);
      weights_sum += weight;

       move += weight * rotated_edge;
    }

    if(weights_sum != 0)
     move /= weights_sum;
    return move;
  }

  Vector rotate_edge(const halfedge_descriptor& main_he, const he_pair& incd_edges)
  {
    // get common vertex around which the edge is rotated
    Point pt = get(vpmap_, target(main_he, mesh_));

    // ps is the vertex that is being moved
    Point ps = get(vpmap_, source(main_he, mesh_));

    // get "equidistant" points - in fact they are at equal angles
    Point equidistant_p1 = get(vpmap_, target(incd_edges.first, mesh_));
    Point equidistant_p2 = get(vpmap_, source(incd_edges.second, mesh_));
    CGAL_assertion(target(incd_edges.second, mesh_) == source(incd_edges.first, mesh_));

    Vector edge1(pt, equidistant_p1);
    Vector edge2(pt, equidistant_p2);
    Vector vec_main_he(pt, ps);

    // check degenerate cases
    double precision = 1e-3;

    if ( edge1.squared_length()               < precision ||
       edge2.squared_length()                 < precision ||
       sqlength(main_he)                      < precision ||
       (edge1 - vec_main_he).squared_length() < precision ||
       (edge2 - vec_main_he).squared_length() < precision )
    {
      return CGAL::NULL_VECTOR;
    }

    CGAL_assertion(vec_main_he.squared_length() > precision);

    // find bisector
    Vector bisector = CGAL::NULL_VECTOR;
    internal::normalize(edge1, traits_);
    internal::normalize(edge2, traits_);
    bisector = edge1 + edge2;

    // special case in angle bisecting: If no edge is actually degenerate
    // but edge1 and edge2 are almost parallel, then by adding them the bisector is
    // almost zero.
    if( bisector.squared_length() < precision )
    {
      // angle is (almost) 180 degrees, take the perpendicular
      // which is normal to edge and tangent to the surface
      Vector normal_vec = find_perpendicular(edge1, pt, ps);

      CGAL_assertion(normal_vec != CGAL::NULL_VECTOR);
      CGAL_assertion(CGAL::scalar_product(edge1, normal_vec) < precision);

      Segment b_segment(pt, pt + normal_vec);
      Point b_segment_end = b_segment.target();

      if(CGAL::angle(b_segment_end, pt, ps) == CGAL::OBTUSE)
      {
        b_segment = b_segment.opposite();
      }
      bisector = Vector(b_segment);
    }

    correct_bisector(bisector, main_he);

    double target_length = CGAL::sqrt(sqlength(main_he));
    double bisector_length = CGAL::sqrt(bisector.squared_length());

    CGAL_assertion(   ( target_length - precision    <   bisector_length     ) &&
             ( bisector_length        <   target_length + precision )    );
    return bisector;
  }

  void correct_bisector(Vector& bisector_vec, const halfedge_descriptor& main_he)
  {
    // get common vertex around which the edge is rotated
    Point pt = get(vpmap_, target(main_he, mesh_));

    // create a segment to be able to translate
    Segment bisector(pt, pt + bisector_vec);

    // scale
    double scale_factor = CGAL::sqrt(  sqlength(main_he) / bisector.squared_length() );
    typename GeomTraits::Aff_transformation_3 t_scale(CGAL::SCALING, scale_factor);
    bisector = bisector.transform(t_scale);

    // translate
    Vector vec(bisector.source(), pt);
    typename GeomTraits::Aff_transformation_3 t_translate(CGAL::TRANSLATION, vec);
    bisector = bisector.transform(t_translate);

    // take the opposite so that their sum is the overall displacement
    bisector_vec = -Vector(bisector);
  }

  Vector find_perpendicular(const Vector& input_vec, const Point& s, const Point& pv)
  {
    Vector s_pv(s, pv);
    Vector aux_normal = CGAL::cross_product(input_vec, s_pv);
    return CGAL::cross_product(aux_normal, input_vec);
  }

  // angle measurement & evaluation
  // ------------------------------
  void measure_angles(const Edges_around_map& he_map)
  {
    min_angle_ = 2 * CGAL_PI;
    typename Edges_around_map::const_iterator it;
    for(it = he_map.begin(); it != he_map.end(); ++it)
    {
      halfedge_descriptor main_he = it->first;
      he_pair incident_pair = it->second;
      calc_angles(main_he, incident_pair);
    }
  }

  void calc_angles(const halfedge_descriptor& main_he, const he_pair& incd_edges)
  {
    // get common vertex around which the edge is rotated
    Point pt = get(vpmap_, target(main_he, mesh_));
    // ps is the vertex that is being moved
    Point ps = get(vpmap_, source(main_he, mesh_));
    // get "equidistant" points - in fact they are at equal angles
    Point equidistant_p1 = get(vpmap_, target(incd_edges.first, mesh_));
    Point equidistant_p2 = get(vpmap_, source(incd_edges.second, mesh_));
    CGAL_assertion(target(incd_edges.second, mesh_) == source(incd_edges.first, mesh_));

    Vector edge1(pt, equidistant_p1);
    Vector edge2(pt, equidistant_p2);
    Vector vec_main_he(pt, ps);

    // angles in [0, 2pi]
    double a1 = get_angle(edge1, vec_main_he);
    double a2 = get_angle(edge2, vec_main_he);

    min_angle_ = std::min(min_angle_, std::min(a1, a2));
  }

  bool does_it_impove(const vertex_descriptor& v, const Point& new_location)
  {
    Edges_around_map he_map;
    typename Edges_around_map::iterator it;
    BOOST_FOREACH(halfedge_descriptor hi, halfedges_around_source(v, mesh_))
      he_map[hi] = he_pair( next(hi, mesh_), prev(opposite(hi, mesh_), mesh_) );
    return evaluate_angles(he_map, new_location);
  }

  bool evaluate_angles(const Edges_around_map& he_map, const Point& new_location)
  {
    typename Edges_around_map::const_iterator it;
    for(it = he_map.begin(); it != he_map.end(); ++it)
    {
      halfedge_descriptor main_he = it->first;
      he_pair incd_edges = it->second;

      // get common vertex around which the edge is rotated
      Point pt = get(vpmap_, target(main_he, mesh_));
      // ps is the vertex that is being moved
      Point new_point = new_location;
      // get "equidistant" points
      Point equidistant_p1 = get(vpmap_, target(incd_edges.first, mesh_));
      Point equidistant_p2 = get(vpmap_, source(incd_edges.second, mesh_));
      CGAL_assertion(target(incd_edges.second, mesh_) == source(incd_edges.first, mesh_));

      Vector edge1(pt, equidistant_p1);
      Vector edge2(pt, equidistant_p2);
      Vector vec_main_he(pt, new_point);

      double a1 = get_angle(edge1, vec_main_he);
      double a2 = get_angle(edge2, vec_main_he);

      if(a1 < min_angle_ || a2 < min_angle_)
        return false;
    }
    return true;
  }

  double get_angle(const Vector& e1, const Vector& e2)
  {
    //double rad_to_deg = 180. / CGAL_PI;
    double cos_angle = (e1 * e2)
     / std::sqrt(e1.squared_length() * e2.squared_length());

    return std::acos(cos_angle); //* rad_to_deg;
  }

  // gradient descent
  // ----------------
  bool gradient_descent(const vertex_descriptor& v, const double& precision)
  {
    bool move_flag;
    double x, y, z, x_new, y_new, z_new, drdx, drdy, drdz;
    x = get(vpmap_, v).x();
    y = get(vpmap_, v).y();
    z = get(vpmap_, v).z();

    double S_av = compute_average_area_around(v);
    double energy = measure_energy(v, S_av);

    // if the adjacent areas are absolutely equal
    if(energy == 0)
      return false;

    double energy_new = 0;
    double relative_energy = 1;
    unsigned int t = 1;
    double eta0 = 0.01;
    //double power_t = 0.25;
    double t0 = 0.001;
    double eta = eta0 / (1 + t0*t);

    while(relative_energy > precision)
    {
      drdx=0, drdy=0, drdz=0;
      compute_derivatives(drdx, drdy, drdz, v, S_av);

      x_new = x - eta * drdx;
      y_new = y - eta * drdy;
      z_new = z - eta * drdz;

      Point moved(x_new, y_new, z_new);
      energy_new = measure_energy(v, S_av, moved);

      if(energy_new < energy)
      {
        put(vpmap_, v, moved);
        move_flag = true;
      }
      else
        return false;

      relative_energy = CGAL::to_double( (energy - energy_new) / energy );

      // update
      x = x_new;
      y = y_new;
      z = z_new;
      energy = energy_new;
      t++;

      //eta = eta0 / pow(t, power_t);
      eta = eta0 / (1 + t0 * t);
    }

    return move_flag;
  }

  void compute_derivatives(double& drdx, double& drdy, double& drdz, const vertex_descriptor& v, const double& S_av)
  {
    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(v, mesh_))
    {
      vertex_descriptor pi = source(next(h, mesh_), mesh_);
      vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
      double S = element_area(v, pi, pi1);

      Vector vec(get(vpmap_, pi), get(vpmap_, pi1));

      // minimize r:
      // r = Σ(S-S_av)^2
      // dr/dx = 2 Σ(S - S_av) dS/dx
      // area of triangle with respect to (x_a, y_a, z_a) =
      // (1/2) [(v_z - v_y)x_a + (v_x - v_z)y_a + (v_y - v_x)z_a + constants]
      // vector v is (x_c - x_b, y_c - y_b, z_c - z_b)
      drdx += (S - S_av) * 0.5 * (vec.z() - vec.y());
      drdy += (S - S_av) * 0.5 * (vec.x() - vec.z());
      drdz += (S - S_av) * 0.5 * (vec.y() - vec.x());
    }

    drdx *= 2;
    drdy *= 2;
    drdz *= 2;
  }

  double element_area(const vertex_descriptor& p1,
            const vertex_descriptor& p2,
            const vertex_descriptor& p3) const
  {
    return to_double(CGAL::approximate_sqrt(
               traits_.compute_squared_area_3_object()(
                  get(vpmap_, p1),
                  get(vpmap_, p2),
                  get(vpmap_, p3))));
  }

  double element_area(const Point& P,
            const vertex_descriptor& p2,
            const vertex_descriptor& p3) const
  {
    return to_double(CGAL::approximate_sqrt(
               traits_.compute_squared_area_3_object()(
                  P,
                  get(vpmap_, p2),
                  get(vpmap_, p3))));
  }

  double compute_average_area_around(const vertex_descriptor& v)
  {
    double sum_areas = 0;
    unsigned int number_of_edges = 0;

    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(v, mesh_))
    {
      // opposite vertices
      vertex_descriptor pi = source(next(h, mesh_), mesh_);
      vertex_descriptor pi1 = target(next(h, mesh_), mesh_);

      double S = element_area(v, pi, pi1);
      sum_areas += S;
      number_of_edges++;
    }
    
    return sum_areas / number_of_edges;
  }

  double measure_energy(const vertex_descriptor& v, const double& S_av)
  {
    double energy = 0;
    unsigned int number_of_edges = 0;

    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(v, mesh_))
    {
      vertex_descriptor pi = source(next(h, mesh_), mesh_);
      vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
      double S = element_area(v, pi, pi1);

      energy += (S - S_av)*(S - S_av);
      number_of_edges++;
    }

    return to_double( energy / number_of_edges );
  }

  double measure_energy(const vertex_descriptor& v, const double& S_av, const Point& new_P)
  {
    double energy = 0;
    unsigned int number_of_edges = 0;
    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(v, mesh_))
    {
      vertex_descriptor pi = source(next(h, mesh_), mesh_);
      vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
      double S = element_area(new_P, pi, pi1);

      energy += (S - S_av)*(S - S_av);
      number_of_edges++;
    }

    return to_double( energy / (2 * number_of_edges) );
  }

  bool is_constrained(const vertex_descriptor& v)
  {
    return get(vcmap_, v);
  }

  template<typename FaceRange>
  void check_vertex_range(const FaceRange& face_range)
  {
   BOOST_FOREACH(face_descriptor f, face_range)
   {
    BOOST_FOREACH(vertex_descriptor v, vertices_around_face(halfedge(f, mesh_), mesh_))
     vrange_.insert(v);
   }
  }

private:

  // data members
  // ------------
  PolygonMesh& mesh_;
  VertexPointMap& vpmap_;
  VertexConstraintMap vcmap_;
  Triangle_list input_triangles_;
  Tree* tree_ptr_;
  GeomTraits traits_;
  std::set<vertex_descriptor> vrange_;
  double min_angle_;

};

} // namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL




#endif // CGAL_POLYGON_MESH_PROCESSING_MESH_SMOOTHING_IMPL_H
