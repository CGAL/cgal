// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : Lloyd move function
//******************************************************************************

#ifndef CGAL_MESH_3_LLOYD_MOVE_H
#define CGAL_MESH_3_LLOYD_MOVE_H

#include <CGAL/Mesh_3/config.h>

#include <CGAL/convex_hull_2.h>
#include <CGAL/ch_graham_andrew.h>

#include <CGAL/Mesh_3/Uniform_sizing_field.h>

#include <string>

namespace CGAL {

namespace Mesh_3 {

template <typename C3T3,
  typename SizingField = Uniform_sizing_field<typename C3T3::Triangulation> >
class Lloyd_move
{
  typedef typename C3T3::Triangulation  Tr;
  typedef typename Tr::Geom_traits      Gt;
  
  typedef typename Tr::Cell_handle    Cell_handle;
  typedef typename Tr::Vertex_handle  Vertex_handle;
  typedef typename Tr::Edge           Edge;
  typedef typename Tr::Facet          Facet;
  typedef typename Tr::Point          Point_3;  
  typedef typename Point_3::Point     Bare_point_3;
  
  typedef typename std::vector<Cell_handle>   Cell_vector;
  typedef typename std::vector<Vertex_handle> Vertex_vector;
  typedef typename std::vector<Facet>         Facet_vector;
  
  typedef typename Gt::FT             FT;
  typedef typename Gt::Vector_3       Vector_3;
  typedef typename Gt::Tetrahedron_3  Tetrahedron_3;
  typedef typename Gt::Plane_3        Plane_3;
  typedef typename Gt::Point_2        Point_2;
  
  typedef typename Gt::Aff_transformation_3   Aff_transformation_3;
  
public:
  typedef SizingField Sizing_field;
  
  /**
   * @brief Return move to apply on \c v according to Lloyd optimization 
   * function
   */
  Vector_3 operator()(const Vertex_handle& v,
                      const Cell_vector& incident_cells,
                      const C3T3& c3t3,
                      const Sizing_field& sizing_field = Sizing_field() ) const
  {
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
   * Project_on_plane defines operator() to project Point_3 object on plane.
   */
  struct Project_on_plane
  {
    Project_on_plane(const Plane_3& plane) : plane_(plane) {}
    
    Point_3 operator()(const Point_3& p) const
    { return Gt().construct_projected_point_3_object()(plane_,p); }
    
  private:
    const Plane_3& plane_;
  };
  
  /**
   * To_2d defines operator() to transform Point_3 into Point_2
   */
  struct To_2d
  {
    To_2d(const Aff_transformation_3& to_2d) : to_2d_(to_2d) {}
    
    Point_2 operator()(const Point_3& p) const
    { return Point_2(to_2d_.transform(p).x(), to_2d_.transform(p).y()); }
    
  private:
    const Aff_transformation_3& to_2d_;
  };
  
  /**
   * To_3d defines operator() to transform Point_2 into Point_3
   */
  struct To_3d
  {
    To_3d(const Aff_transformation_3& to_3d) : to_3d_(to_3d) {}
    
    Point_3 operator()(const Point_2& p) const
    { return to_3d_.transform(Point_3(Bare_point_3(p.x(),p.y(),0))); }
    
  private:
    const Aff_transformation_3& to_3d_;
  };
  
 
  /**
   * Return move for inside vertex \c v
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
    // get cells incident to v
    const Tr& tr = c3t3.triangulation();
    
    // Stores vertex that have already been used
    std::set<Vertex_handle> treated_vertex;
    
    for (typename Cell_vector::const_iterator cit = incident_cells.begin();
         cit != incident_cells.end();
         ++cit)
    {
      const int& k = (*cit)->index(v);
      
      // For each vertex of the cell
      for ( int i=1 ; i<4 ; ++i )
      {
        const Vertex_handle& v1 = (*cit)->vertex((k+i)&3);
        if ( treated_vertex.find(v1) != treated_vertex.end() )
          continue;
        
        // Vertex has not been treated: turn around edge(v,v1)
        treated_vertex.insert(v1);
        turn_around_edge(v, Edge(*cit,k,(k+i)&3), tr,
                         move, sum_masses, sizing_field);
      }
    }
    
    CGAL_assertion(sum_masses != 0.0);
    return move / sum_masses;
  }
  
  /**
   * Return move for on boundary vertex \c v
   */
  Vector_3 lloyd_move_on_boundary(const Vertex_handle& v,
                                  const C3T3& c3t3,
                                  const Sizing_field& sizing_field) const
  {
    CGAL_precondition(c3t3.in_dimension(v) == 2);
    
    // get all surface delaunay ball point
    std::vector<Point_3> points = extract_lloyd_boundary_points(v,c3t3);
    
    switch(points.size())
    {
      case 0: // could happen if there is an isolated surface point into mesh
      case 1: // don't do anything, as v->point() is already on the surface
      {
        return CGAL::NULL_VECTOR;
        break;
      }
      case 2: // centroid 
      {
        const Point_3& a = points.front();
        const Point_3& b = points.back();
        return centroid_segment_move(v,a,b,sizing_field);
        break;
      }
      case 3: // triangle centroid 
      {
        const Point_3& a = points.at(0);
        const Point_3& b = points.at(1);
        const Point_3& c = points.at(2);
        return centroid_triangle_move(v,a,b,c,sizing_field);
        break;
      }
      default: // >= 4 points, centroid + projection
        return centroid_general_move(v,points.begin(),points.end(),sizing_field);
        break;
    }
    
    return CGAL::NULL_VECTOR;
  }
  
  /**
   * Returns a vector containing surface delaunay ball center of surface
   * facets incident to vertex \c v
   */
  std::vector<Point_3> extract_lloyd_boundary_points(const Vertex_handle& v,
                                                     const C3T3& c3t3) const
  {
    Facet_vector incident_facets;
    incident_facets.reserve(64);
    c3t3.triangulation().finite_incident_facets(v,std::back_inserter(incident_facets));
    
    std::vector<Point_3> points;
    points.reserve(64);
    
    // Get c3t3's facets incident to v, and add their surface delaunay ball
    // center to output
    for ( typename Facet_vector::iterator fit = incident_facets.begin() ;
         fit != incident_facets.end() ;
         ++fit )
    {
      if ( ! c3t3.is_in_complex(*fit) )
        continue;
      
      points.push_back(fit->first->get_facet_surface_center(fit->second));
    }
    
    return points;
  }
  
  /**
   * Return move from \c v to centroid of segment [a,b]
   */
  Vector_3 centroid_segment_move(const Vertex_handle& v,
                                 const Point_3& a,
                                 const Point_3& b,
                                 const Sizing_field& sizing_field) const
  {
    typename Gt::Construct_vector_3 vector =
      Gt().construct_vector_3_object();
    
    const Point_3& p = v->point();
    
    FT da = density_1d(a,v,sizing_field);
    FT db = density_1d(b,v,sizing_field);

    CGAL_assertion( !is_zero(da+db) ); 
    return ( (vector(p,a)*da + vector(p,b)*db) / (da+db) );
  }
  
  /**
   * Return move from \c v to centroid of triangle [a,b,c]
   */
  Vector_3 centroid_triangle_move(const Vertex_handle& v,
                                  const Point_3& a,
                                  const Point_3& b,
                                  const Point_3& c,
                                  const Sizing_field& sizing_field) const
  {
    typename Gt::Construct_vector_3 vector =
      Gt().construct_vector_3_object();
    
    const Point_3& p = v->point();
    
    FT da = density_2d<true>(a,v,sizing_field);
    FT db = density_2d<false>(b,v,sizing_field);
    FT dc = density_2d<false>(c,v,sizing_field);

    CGAL_assertion( !is_zero(da+db+dc) );
    return ( (da*vector(p,a) + db*vector(p,b) + dc*vector(p,c)) / (da+db+dc) );
  }
  
  /**
   * compute approximate centroid of intersection between 3D voronoi cell and
   * boundary. The input is the set of intersection points between Voronoi 
   * edges and the boundary.
   */
  template <typename ForwardIterator>
  Vector_3 centroid_general_move(const Vertex_handle& v,
                                 ForwardIterator first,
                                 ForwardIterator last,
                                 const Sizing_field& sizing_field) const
  {
    CGAL_assertion(std::distance(first,last) > 3);
    
    // Fit plane using point-based PCA: compute least square fitting plane
    Plane_3 plane;
    CGAL::linear_least_squares_fitting_3(first,last,plane,Dimension_tag<0>());
    
    // Project all points to the plane
    std::transform(first, last, first, Project_on_plane(plane));
    CGAL_assertion(std::distance(first,last) > 3);
    
    // Get 2D-3D transformations
    Aff_transformation_3 to_3d = compute_to_3d_transform(plane, *first);
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
    std::vector<Point_3> polygon_3d;
    polygon_3d.reserve(ch_2d.size());
    std::transform(ch_2d.begin(), ch_2d.end(),
                   std::back_inserter(polygon_3d), To_3d(to_3d));
    
    // Compute centroid using quadrature sizing
    return centroid_3d_polygon_move(v,
                                    polygon_3d.begin(), polygon_3d.end(),
                                    sizing_field);
  }
  
  /**
   * Return move from \c v to centroid of polygon[first,last]
   * Polygon has to be convex
   */
  template <typename ForwardIterator>
  Vector_3 centroid_3d_polygon_move(const Vertex_handle& v,
                                    ForwardIterator first,
                                    ForwardIterator last,
                                    const Sizing_field& sizing_field) const
  {
    CGAL_precondition(std::distance(first,last) >= 3);
    
    typename Gt::Construct_vector_3 vector =
      Gt().construct_vector_3_object();
    
    typename Gt::Construct_centroid_3 centroid =
      Gt().construct_centroid_3_object();
    
    typename Gt::Compute_area_3 area = 
      Gt().compute_area_3_object();
    
    // Vertex current position
    const Point_3& vertex_position = v->point();
    
    // Use as reference point to triangulate
    const Point_3& a = *first++;
    const Point_3* b = &(*first++);
    
    // Treat first point (optimize density_2d call)
    const Point_3& c = *first++;
    
    Point_3 triangle_centroid = centroid(a,*b,c);
    FT density = density_2d<true>(triangle_centroid, v, sizing_field);
    
    FT sum_masses = density * area(a,*b,c);
    Vector_3 move = sum_masses * vector(vertex_position, triangle_centroid);
    
    b = &c;
    
    // Next points
    while ( first != last )
    {
      const Point_3& c = *first++;

      Point_3 triangle_centroid = centroid(a,*b,c);
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
   * Returns the transformation from reference_point to plane
   */
  Aff_transformation_3 compute_to_3d_transform(const Plane_3& plane,
                                               const Point_3& reference_point) const
  {
    typename Gt::Construct_base_vector_3 base =
      Gt().construct_base_vector_3_object();
    
    typename Gt::Construct_orthogonal_vector_3 orthogonal_vector =
      Gt().construct_orthogonal_vector_3_object();
    
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
   * returns density_1d
   */
  template <typename Sizing_field>
  FT density_1d(const Point_3& p,
                const Vertex_handle& v,
                const Sizing_field& sizing_field) const
  {
    FT s = sizing_field(p,v);
    CGAL_assertion(!is_zero(s));

    // s^(d+2)
    return ( 1/(s*s*s) );
  }
  
  /**
   * returns density_2d
   */
  template <bool use_v, typename Sizing_field>
  FT density_2d(const Point_3& p,
                const Vertex_handle& v,
                const Sizing_field& sizing_field) const
  {
    FT s = use_v ? sizing_field(p,v) : sizing_field(p);
    CGAL_assertion(!is_zero(s));

    // s^(d+2)
    return ( 1/(s*s*s*s) );
  }
  
  /**
   * returns density_3d
   */
  template <typename Sizing_field>
  FT density_3d(const Point_3& p,
                const Cell_handle& cell,
                const Sizing_field& sizing_field) const
  {
    FT s = sizing_field(p,cell);
    CGAL_assertion(!is_zero(s));

    // s^(d+2)
    return ( 1/(s*s*s*s*s) );
  }

  /**
   * Turns around edge \c edge and add values computed from tets made by
   * v->point() and circumcenters of \c edge incident cells
   *
   * Note that this function uses lots of dual() calls, so using a cell base
   * which stores the circumcenter improve efficiency.
   */
  void turn_around_edge(const Vertex_handle& v,
                        const Edge& edge,
                        const Tr& tr,
                        Vector_3& move,
                        FT& sum_masses,
                        const Sizing_field& sizing_field) const
  {
    typedef typename Tr::Cell_circulator Cell_circulator;
    
    typename Gt::Compute_volume_3 volume = 
      Gt().compute_volume_3_object();
    
    typename Gt::Construct_centroid_3 centroid =
      Gt().construct_centroid_3_object();
    
    typename Gt::Construct_vector_3 vector =
      Gt().construct_vector_3_object();
    
    Cell_circulator current_cell = tr.incident_cells(edge);
    Cell_circulator done = current_cell;
    
    // a & b are fixed points
    const Point_3& a = v->point();
    const Point_3 b = tr.dual(current_cell++);
    CGAL_assertion(current_cell != done);
    
    // c & d are moving points
    Point_3 c = tr.dual(current_cell++);
    CGAL_assertion(current_cell != done);
    
    while ( current_cell != done )
    {
      const Point_3 d = tr.dual(current_cell++);
      
      Point_3 tet_centroid = centroid(a,b,c,d);
      
      // Compute mass
      FT density = density_3d(tet_centroid, current_cell, sizing_field);
      FT abs_volume = CGAL::abs(volume(a,b,c,d));
      FT mass = abs_volume * density;
      
      move = move + mass * vector(a,tet_centroid);
      sum_masses += mass;
      
      c = d;
    }
  }
};

} // end namespace Mesh_3 


} //namespace CGAL

#endif // CGAL_MESH_3_LLOYD_MOVE_H
