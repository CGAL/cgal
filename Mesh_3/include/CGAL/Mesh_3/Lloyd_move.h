// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : Lloyd move function
//******************************************************************************

#ifndef CGAL_MESH_3_LLOYD_MOVE_H
#define CGAL_MESH_3_LLOYD_MOVE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>
#include <CGAL/Mesh_3/config.h>
#include <CGAL/Mesh_3/Uniform_sizing_field.h>

#include <CGAL/Time_stamper.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/ch_graham_andrew.h>

#include <boost/unordered_set.hpp>

#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

namespace CGAL {

namespace Mesh_3 {

template <typename C3T3,
          typename SizingField = Uniform_sizing_field<typename C3T3::Triangulation> >
class Lloyd_move
{
  typedef typename C3T3::Triangulation                        Tr;
  typedef typename Tr::Geom_traits                            Gt;

  typedef typename Tr::Vertex_handle                          Vertex_handle;
  typedef typename Tr::Edge                                   Edge;
  typedef typename Tr::Facet                                  Facet;
  typedef typename Tr::Cell_handle                            Cell_handle;

  typedef typename Tr::Bare_point                             Bare_point;
  typedef typename Tr::Weighted_point                         Weighted_point;

  typedef typename std::vector<Facet>                         Facet_vector;
  typedef typename std::vector<Cell_handle>                   Cell_vector;

  typedef typename Gt::FT                                     FT;
  typedef typename Gt::Point_2                                Point_2;
  typedef typename Gt::Vector_3                               Vector_3;
  typedef typename Gt::Tetrahedron_3                          Tetrahedron_3;
  typedef typename Gt::Plane_3                                Plane_3;
  typedef typename Gt::Aff_transformation_3                   Aff_transformation_3;

public:
  typedef SizingField                                         Sizing_field;

  /**
   * @brief Return the move to apply on \c v according to Lloyd optimization
   * function.
   */
  Vector_3 operator()(const Vertex_handle& v,
                      const Cell_vector& incident_cells,
                      const C3T3& c3t3,
                      const Sizing_field& sizing_field = Sizing_field() ) const
  {
#ifdef CGAL_MESH_3_OPTIMIZER_DEBUG_VERBOSE
    std::cout << "computing move of: " << &*v
              << " pos: " << c3t3.triangulation().point(v)
              << " dim: " << c3t3.in_dimension(v) << std::endl;
#endif

    switch ( c3t3.in_dimension(v) )
    {
      case 3:
        return lloyd_move_inside_domain(v,incident_cells,c3t3,sizing_field);
        break;
      case 2:
        return lloyd_move_on_boundary(v,c3t3,sizing_field);
        break;
      case 1:
      case 0:
      case -1:
        // Don't move edge or corner vertices
        // N.B.: dimension = -1 is possible if we added points on a far sphere
        //       during initialization
        return CGAL::NULL_VECTOR;
        break;
      default:
        // Should not happen
        CGAL_assertion(false);
        return CGAL::NULL_VECTOR;
        break;
    }

    return CGAL::NULL_VECTOR;
  }

#if defined(CGAL_MESH_3_OPTIMIZER_VERBOSE) \
 || defined (CGAL_MESH_3_EXPORT_PERFORMANCE_DATA)
  static std::string name() { return std::string("Lloyd"); }
#endif

private:
  /**
   * Project_on_plane defines `operator()` to project a point object on the plane `plane`.
   */
  struct Project_on_plane
  {
    Project_on_plane(const Plane_3& plane, const C3T3& c3t3)
      : plane_(plane), c3t3_(c3t3)
    { }

    Bare_point operator()(const Bare_point& p) const
    { return c3t3_.triangulation().geom_traits().
          construct_projected_point_3_object()(plane_, p); }

  private:
    const Plane_3& plane_;
    const C3T3& c3t3_;
  };

  /**
   * `To_2d` defines `operator()` to transform from `Bare_point` to `Point_2`.
   */
  struct To_2d
  {
    To_2d(const Aff_transformation_3& to_2d) : to_2d_(to_2d) {}

    Point_2 operator()(const Bare_point& p) const
    { return Point_2(to_2d_.transform(p).x(), to_2d_.transform(p).y()); }

  private:
    const Aff_transformation_3& to_2d_;
  };

  /**
   * `To_3d` defines `operator()` to transform from `Point_2` to `Bare_point`.
   */
  struct To_3d
  {
    To_3d(const Aff_transformation_3& to_3d) : to_3d_(to_3d) {}

    Bare_point operator()(const Point_2& p) const
    { return to_3d_.transform((Bare_point(p.x(),p.y(),0))); }

  private:
    const Aff_transformation_3& to_3d_;
  };


  /**
   * Return the move for the inside vertex \c v.
   */
  Vector_3 lloyd_move_inside_domain(const Vertex_handle& v,
                                    const Cell_vector& incident_cells,
                                    const C3T3& c3t3,
                                    const Sizing_field& sizing_field) const
  {
    // Move data
    Vector_3 move = CGAL::NULL_VECTOR;
    FT sum_masses(0);

    // -----------------------------------
    // Tessellate Voronoi cell
    // -----------------------------------

    // Stores vertex that have already been used
    typedef CGAL::Hash_handles_with_or_without_timestamps   Hash_fct;
    typedef boost::unordered_set<Vertex_handle, Hash_fct>   Vertex_container;
    typedef typename Vertex_container::iterator             VC_it;

    Vertex_container treated_vertices;

    for (typename Cell_vector::const_iterator cit = incident_cells.begin();
         cit != incident_cells.end();
         ++cit)
    {
      const int& k = (*cit)->index(v);

      // For each vertex of the cell
      for ( int i=1 ; i<4 ; ++i )
      {
        const Vertex_handle& v1 = (*cit)->vertex((k+i)&3);

        std::pair<VC_it, bool> is_insert_successful = treated_vertices.insert(v1);
        if ( ! is_insert_successful.second ) // vertex has already been treated
          continue;

        // Vertex has not been treated: turn around edge(v,v1)
        turn_around_edge(v, Edge(*cit,k,(k+i)&3), c3t3,
                         move, sum_masses, sizing_field);
      }
    }

    CGAL_assertion(sum_masses != 0.0);
    return move / sum_masses;
  }

  /**
   * Return the move for the on-boundary vertex \c v.
   */
  Vector_3 lloyd_move_on_boundary(const Vertex_handle& v,
                                  const C3T3& c3t3,
                                  const Sizing_field& sizing_field) const
  {
    CGAL_precondition(c3t3.in_dimension(v) == 2);

    // get all surface delaunay ball point
    std::vector<Bare_point> points = extract_lloyd_boundary_points(v,c3t3);

    switch(points.size())
    {
      case 0: // could happen if there is an isolated surface point into mesh
      case 1: // don't do anything, as the point is already on the surface
      {
        return CGAL::NULL_VECTOR;
        break;
      }
      case 2: // segment centroid
      {
        const Bare_point& a = points.front();
        const Bare_point& b = points.back();
        return centroid_segment_move(v, a, b, c3t3, sizing_field);
        break;
      }
      case 3: // triangle centroid
      {
        const Bare_point& a = points.at(0);
        const Bare_point& b = points.at(1);
        const Bare_point& c = points.at(2);
        return centroid_triangle_move(v, a, b, c, c3t3, sizing_field);
        break;
      }
      default: // >= 4 points, centroid + projection
        return centroid_general_move(v, points.begin(), points.end(), c3t3, sizing_field);
        break;
    }

    return CGAL::NULL_VECTOR;
  }

  /**
   * Returns a vector containing the surface delaunay ball centers of the surface
   * facets that are incident to vertex \c v.
   */
  std::vector<Bare_point> extract_lloyd_boundary_points(const Vertex_handle& v,
                                                        const C3T3& c3t3) const
  {
    const Tr& tr = c3t3.triangulation();

    typename Gt::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

    Facet_vector incident_facets;
    incident_facets.reserve(64);
    tr.finite_incident_facets(v, std::back_inserter(incident_facets));

    std::vector<Bare_point> points;
    points.reserve(64);

    const Weighted_point& position = tr.point(v);

    // Get c3t3's facets incident to v, and add their surface delaunay ball
    // center to output
    for ( typename Facet_vector::iterator fit = incident_facets.begin() ;
                                          fit != incident_facets.end() ;
                                          ++fit )
    {
      if ( ! c3t3.is_in_complex(*fit) )
        continue;

      const Cell_handle cell = fit->first;
      const int i = fit->second;

      // In the case of a periodic triangulation, the incident facets of a point
      // do not necessarily have the same offsets. Worse, the surface centers
      // might not have the same offset as their facet. Thus, no solution except
      // calling a function 'get_closest_point(p, q)' that simply returns q
      // for a non-periodic triangulation, and checks all possible offsets for
      // periodic triangulations
      points.push_back(tr.get_closest_point(cp(position),
                                            cell->get_facet_surface_center(i)));
    }

    return points;
  }

  /**
   * Return the move from \c v to the centroid of the segment [a,b].
   */
  Vector_3 centroid_segment_move(const Vertex_handle& v,
                                 const Bare_point& a,
                                 const Bare_point& b,
                                 const C3T3& c3t3,
                                 const Sizing_field& sizing_field) const
  {
    typename Gt::Construct_point_3 cp = c3t3.triangulation().geom_traits().construct_point_3_object();
    typename Gt::Construct_vector_3 vector = c3t3.triangulation().geom_traits().construct_vector_3_object();

    const Weighted_point position = c3t3.triangulation().point(v);
    const Bare_point& p = cp(position);

    FT da = density_1d(a, v, sizing_field);
    FT db = density_1d(b, v, sizing_field);

    CGAL_assertion( !is_zero(da + db) );
    return ( (vector(p,a)*da + vector(p,b)*db) / (da + db) );
  }

  /**
   * Return the move from \c v to the centroid of triangle [a,b,c].
   */
  Vector_3 centroid_triangle_move(const Vertex_handle& v,
                                  const Bare_point& a,
                                  const Bare_point& b,
                                  const Bare_point& c,
                                  const C3T3& c3t3,
                                  const Sizing_field& sizing_field) const
  {
    typename Gt::Construct_point_3 cp = c3t3.triangulation().geom_traits().construct_point_3_object();
    typename Gt::Construct_vector_3 vector = c3t3.triangulation().geom_traits().construct_vector_3_object();

    const Weighted_point& position = c3t3.triangulation().point(v);
    const Bare_point& p = cp(position);

    FT da = density_2d<true>(a,v,sizing_field);
    FT db = density_2d<false>(b,v,sizing_field);
    FT dc = density_2d<false>(c,v,sizing_field);

    CGAL_assertion( !is_zero(da+db+dc) );
    return ( (da*vector(p,a) + db*vector(p,b) + dc*vector(p,c)) / (da+db+dc) );
  }

  /**
   * Compute the approximate centroid of the intersection between the 3D voronoi
   * cell and the boundary. The input is the set of intersection points between
   * Voronoi edges and the boundary.
   */
  template <typename ForwardIterator>
  Vector_3 centroid_general_move(const Vertex_handle& v,
                                 ForwardIterator first,
                                 ForwardIterator last,
                                 const C3T3& c3t3,
                                 const Sizing_field& sizing_field) const
  {
    CGAL_assertion(std::distance(first,last) > 3);

    // Fit plane using point-based PCA: compute least square fitting plane
    Plane_3 plane;
    Bare_point point;
    CGAL::linear_least_squares_fitting_3(first, last, plane, point, Dimension_tag<0>(),
                                         c3t3.triangulation().geom_traits(),
                                         Default_diagonalize_traits<FT, 3>());

    // Project all points to the plane
    std::transform(first, last, first, Project_on_plane(plane, c3t3));
    CGAL_assertion(std::distance(first, last) > 3);

    // Get 2D-3D transformations
    Aff_transformation_3 to_3d = compute_to_3d_transform(plane, *first, c3t3);
    Aff_transformation_3 to_2d = to_3d.inverse();

    // Transform to 2D
    std::vector<Point_2> points_2d;
    points_2d.reserve(std::distance(first,last));
    std::transform(first, last, std::back_inserter(points_2d), To_2d(to_2d));

    // Compute 2D convex hull
    CGAL_assertion(points_2d.size() > 3);
    std::vector<Point_2> ch_2d;
    // AF: We have to debug CGAL::convex_hull_2 = ch_akl_toussaint
    //     as it triggers filter failures unnecessarily
    CGAL::ch_graham_andrew(points_2d.begin(),points_2d.end(),
                           std::back_inserter(ch_2d));

    // Lift back convex hull to 3D
    std::vector<Bare_point> polygon_3d;
    polygon_3d.reserve(ch_2d.size());
    std::transform(ch_2d.begin(), ch_2d.end(),
                   std::back_inserter(polygon_3d), To_3d(to_3d));

    // Compute centroid using quadrature sizing
    return centroid_3d_polygon_move(v, polygon_3d.begin(), polygon_3d.end(),
                                    c3t3, sizing_field);
  }

  /**
   * Return the move from \c v to the centroid of polygon[first,last].
   * The polygon has to be convex.
   */
  template <typename ForwardIterator>
  Vector_3 centroid_3d_polygon_move(const Vertex_handle& v,
                                    ForwardIterator first,
                                    ForwardIterator last,
                                    const C3T3& c3t3,
                                    const Sizing_field& sizing_field) const
  {
    CGAL_precondition(std::distance(first,last) >= 3);

    typename Gt::Compute_area_3 area = c3t3.triangulation().geom_traits().compute_area_3_object();
    typename Gt::Construct_centroid_3 centroid = c3t3.triangulation().geom_traits().construct_centroid_3_object();
    typename Gt::Construct_point_3 cp = c3t3.triangulation().geom_traits().construct_point_3_object();
    typename Gt::Construct_vector_3 vector = c3t3.triangulation().geom_traits().construct_vector_3_object();

    // Vertex current position
    const Weighted_point& vertex_weighted_position = c3t3.triangulation().point(v);
    const Bare_point& vertex_position = cp(vertex_weighted_position);

    // Use as reference point to triangulate
    const Bare_point& a = *first++;
    const Bare_point* b = &(*first++);

    // Treat first point (optimize density_2d call)
    const Bare_point& c = *first++;

    const Bare_point& triangle_centroid = centroid(a,*b,c);
    FT density = density_2d<true>(triangle_centroid, v, sizing_field);

    FT sum_masses = density * area(a,*b,c);
    Vector_3 move = sum_masses * vector(vertex_position, triangle_centroid);

    b = &c;

    // Next points
    while ( first != last )
    {
      const Bare_point& c = *first++;

      const Bare_point& triangle_centroid = centroid(a,*b,c);
      FT density = density_2d<false>(triangle_centroid, v, sizing_field);
      FT mass = density * area(a,*b,c);

      move = move + mass * vector(vertex_position, triangle_centroid);
      sum_masses += mass;

      // Go one step forward
      b = &c;
    }
    CGAL_assertion(sum_masses != 0);

    return move / sum_masses;
  }

  /**
   * Return the transformation from `reference_point` to `plane`.
   */
  Aff_transformation_3 compute_to_3d_transform(const Plane_3& plane,
                                               const Bare_point& reference_point,
                                               const C3T3& c3t3) const
  {
    typename Gt::Construct_base_vector_3 base = c3t3.triangulation().geom_traits().construct_base_vector_3_object();
    typename Gt::Construct_orthogonal_vector_3 orthogonal_vector = c3t3.triangulation().geom_traits().construct_orthogonal_vector_3_object();

    Vector_3 u = base(plane, 1);
    u = u / CGAL::sqrt(u*u);

    Vector_3 v = base(plane, 2);
    v = v / CGAL::sqrt(v*v);

    Vector_3 w = orthogonal_vector(plane);
    w = w / CGAL::sqrt(w*w);

    return Aff_transformation_3(u.x(),v.x(),w.x(),reference_point.x(),
                                u.y(),v.y(),w.y(),reference_point.y(),
                                u.z(),v.z(),w.z(),reference_point.z());
  }

  /**
   * Return density_1d
   */
  template <typename Sizing_field>
  FT density_1d(const Bare_point& p,
                const Vertex_handle& v,
                const Sizing_field& sizing_field) const
  {
    FT s = sizing_field(p,v);
    CGAL_assertion(!is_zero(s));

    // s^(d+2)
    return ( 1/(s*s*s) );
  }

  /**
   * Return density_2d
   */
  template <bool use_v, typename Sizing_field>
  FT density_2d(const Bare_point& p,
                const Vertex_handle& v,
                const Sizing_field& sizing_field) const
  {
    FT s = use_v ? sizing_field(p,v) : sizing_field(p);
    CGAL_assertion(!is_zero(s));

    // s^(d+2)
    return ( 1/(s*s*s*s) );
  }

  /**
   * Return density_3d
   */
  template <typename Sizing_field>
  FT density_3d(const Bare_point& p,
                const Cell_handle& cell,
                const Sizing_field& sizing_field) const
  {
    FT s = sizing_field(p,cell);
    CGAL_assertion(!is_zero(s));

    // s^(d+2)
    return ( 1/(s*s*s*s*s) );
  }

  /**
   * Turn around the edge \c edge and add the values computed from tets made by
   * `v` and the circumcenters of cells incident to \c edge.
   *
   * Note that this function abundantly uses dual() calls and using a cell base
   * which stores the circumcenter thus improves its efficiency.
   */
  void turn_around_edge(const Vertex_handle& v,
                        const Edge& edge,
                        const C3T3& c3t3,
                        Vector_3& move,
                        FT& sum_masses,
                        const Sizing_field& sizing_field) const
  {
    CGAL_precondition(c3t3.in_dimension(v) == 3);

    typedef typename Tr::Cell_circulator Cell_circulator;

    const Tr& tr = c3t3.triangulation();

    typename Gt::Construct_centroid_3 centroid = tr.geom_traits().construct_centroid_3_object();
    typename Gt::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();
    typename Gt::Construct_tetrahedron_3 tetrahedron = tr.geom_traits().construct_tetrahedron_3_object();
    typename Gt::Construct_translated_point_3 translate = tr.geom_traits().construct_translated_point_3_object();
    typename Gt::Construct_vector_3 vector = tr.geom_traits().construct_vector_3_object();
    typename Gt::Compute_volume_3 volume = tr.geom_traits().compute_volume_3_object();

    Cell_circulator current_cell = tr.incident_cells(edge);
    Cell_circulator done = current_cell;

    // a & b are fixed points
    const Weighted_point& wa = tr.point(v);
    const Bare_point& a = cp(wa);

    Bare_point b = tr.dual(current_cell);
    const Weighted_point& a_b = tr.point(current_cell, current_cell->index(v));
    Vector_3 ba = Vector_3(cp(a_b), b);
    ++current_cell;
    CGAL_assertion(current_cell != done);

    // c & d are moving points
    Bare_point c = tr.dual(current_cell);
    const Weighted_point& a_c = tr.point(current_cell, current_cell->index(v));
    Vector_3 ca = Vector_3(cp(a_c), c);
    ++current_cell;
    CGAL_assertion(current_cell != done);

    while ( current_cell != done )
    {
      Bare_point d = tr.dual(current_cell);
      const Weighted_point& a_d = tr.point(current_cell, current_cell->index(v));
      Vector_3 da = Vector_3(cp(a_d), d);

      Tetrahedron_3 tet = tetrahedron(a, translate(a, ba), translate(a, ca), translate(a, da));
      Bare_point tet_centroid = centroid(tet);

      // Compute mass
      FT density = density_3d(tet_centroid, current_cell, sizing_field);
      FT abs_volume = CGAL::abs(volume(tet));
      FT mass = abs_volume * density;

      move = move + mass * vector(a, tet_centroid);
      sum_masses += mass;

      ++current_cell;
      c = d;
      ca = da;
    }
  }
};

} // end namespace Mesh_3


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_LLOYD_MOVE_H
