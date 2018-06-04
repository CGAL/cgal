// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot
//

#ifndef CGAL_STRUCTURE_POINT_SET_3_H
#define CGAL_STRUCTURE_POINT_SET_3_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/assertions.h>
#include <CGAL/intersections.h>

#include <CGAL/centroid.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits_3.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/iterator/counting_iterator.hpp>

#include <iterator>
#include <list>
#include <limits>

namespace CGAL {

/*!
\ingroup PkgPointSetProcessingAlgorithms

\brief A 3D point set with structure information based on a set of
detected planes.

Given a point set in 3D space along with a set of fitted planes, this
class stores a simplified and structured version of the point
set. Each output point is assigned to one, two or more primitives
(depending wether it belongs to a planar section, an edge or a if it
is a vertex). The implementation follow \cgalCite{cgal:la-srpss-13}.

\tparam Kernel a model of `ShapeDetectionTraits` that must provide in
addition a function `Intersect_3 intersection_3_object() const` and a
functor `Intersect_3` with:
- `boost::optional< boost::variant< Traits::Plane_3, Traits::Line_3 > > operator()(typename Traits::Plane_3, typename Traits::Plane_3)`
- `boost::optional< boost::variant< Traits::Line_3, Traits::Point_3 > > operator()(typename Traits::Line_3, typename Traits::Plane_3)`

*/
template <typename Kernel>
class Point_set_with_structure
{
  typedef Point_set_with_structure<Kernel> Self;

  typedef typename Kernel::FT FT;
  typedef typename Kernel::Segment_3 Segment;
  typedef typename Kernel::Line_3 Line;
  typedef typename Kernel::Point_2 Point_2;

  enum Point_status { POINT, RESIDUS, PLANE, EDGE, CORNER, SKIPPED };

public:


  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Plane_3 Plane;

  /// Tag classifying the coherence of a triplet of points with
  /// respect to an inferred surface
  enum Coherence_type
    {
      INCOHERENT = -1, ///< Incoherent (facet violates the underlying structure)
      FREEFORM = 0,    ///< Free-form coherent (facet is between 3 free-form points)
      VERTEX = 1,      ///< Structure coherent, facet adjacent to a vertex
      CREASE = 2,      ///< Structure coherent, facet adjacent to an edge
      PLANAR = 3       ///< Structure coherent, facet inside a planar section
    };
  
private:

  class My_point_property_map{
    const std::vector<Point>& points;
  public:
    typedef Point value_type;
    typedef const value_type& reference;
    typedef std::size_t key_type;
    typedef boost::lvalue_property_map_tag category;  
    My_point_property_map (const std::vector<Point>& pts) : points (pts) {}
    reference operator[] (key_type k) const { return points[k]; }
    friend inline reference get (const My_point_property_map& ppmap, key_type i) 
    { return ppmap[i]; }
  };

  struct Edge
  {
    CGAL::cpp11::array<std::size_t, 2> planes;
    std::vector<std::size_t> indices; // Points belonging to intersection
    Line support;
    bool active;

    Edge (std::size_t a, std::size_t b)
      : support (Point (FT(0.), FT(0.), FT(0.)),
                 Vector (FT(0.), FT(0.), FT(0.)))
      , active(true)
    { planes[0] = a; planes[1] = b; }
  };
  struct Corner
  {
    std::vector<std::size_t> planes;
    std::vector<std::size_t> edges;
    std::vector<Vector> directions;
    Point support;
    bool active;

    Corner (std::size_t p1, std::size_t p2, std::size_t p3,
            std::size_t e1, std::size_t e2, std::size_t e3)
    {
      planes.resize (3); planes[0] = p1; planes[1] = p2; planes[2] = p3;
      edges.resize (3); edges[0] = e1; edges[1] = e2; edges[2] = e3;
      active = true;
    }
  };
      

  std::vector<Point> m_points;
  std::vector<Vector> m_normals;
  std::vector<std::size_t> m_indices;
  std::vector<Point_status> m_status;
    
  std::vector<Plane> m_planes;
  std::vector<std::vector<std::size_t> > m_indices_of_assigned_points;
  std::vector<Edge> m_edges;
  std::vector<Corner> m_corners;
    
public:


  /*!
    Constructs a structured point set based on the input points and the
    associated shape detection object.

    \tparam PointRange is a model of `ConstRange`. The value type of
    its iterator is the key type of the named parameter `point_map`.
    \tparam PlaneRange is a model of `ConstRange`. The value type of
    its iterator is the key type of the named parameter `plane_map`.

    \param points input point range.
    \param planes input plane range.
    \param epsilon size parameter.
    \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

    \cgalNamedParamsBegin
      \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `Kernel::Point_3`.
      If this parameter is omitted, `CGAL::Identity_property_map<Kernel::Point_3>` is used.\cgalParamEnd
      \cgalParamBegin{normal_map} a model of `ReadablePropertyMap` with value type
      `Kernel::Vector_3`.\cgalParamEnd
      \cgalParamBegin{plane_index_map} a model of `ReadablePropertyMap` with value type `int`.
      Associates the index of a point in the input range to the index of plane (-1 if point does is not assigned to
      a plane).\cgalParamEnd
      \cgalParamBegin{plane_map} a model of `ReadablePropertyMap` with value type
      `Kernel::Plane_3`. If this parameter is omitted, `CGAL::Identity_property_map<Kernel::Plane_3>`
      is used.\cgalParamEnd
      \cgalParamBegin{attraction_factor} multiple of `epsilon` used to connect simplices.\cgalParamEnd
    \cgalNamedParamsEnd

  */
  template <typename PointRange,
            typename PlaneRange,
            typename NamedParameters>
  Point_set_with_structure (const PointRange& points,
                            const PlaneRange& planes,
                            double epsilon,
                            const NamedParameters& np)
  {
    init (points, planes, epsilon, np);
  }

  /// \cond SKIP_IN_MANUAL
  // deprecated
  template <typename PointRange,
            typename PointMap,
            typename NormalMap,
            typename PlaneRange,
            typename PlaneMap,
            typename IndexMap>
  CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::Point_set_with_structure(), please update your code")
  Point_set_with_structure (const PointRange& points,
                            PointMap point_map,
                            NormalMap normal_map,
                            const PlaneRange& planes,
                            PlaneMap plane_map,
                            IndexMap index_map,
                            double epsilon,
                            double attraction_factor = 3.)
  {
    init (points, planes, epsilon,
          CGAL::parameters::point_map (point_map).
          normal_map (normal_map).
          plane_map (plane_map).
          plane_index_map (index_map).
          attraction_factor (attraction_factor));
  }

  template <typename PointRange,
            typename PlaneRange,
            typename NamedParameters>
  void init (const PointRange& points,
             const PlaneRange& planes,
             double epsilon,
             const NamedParameters& np)
  {
    using boost::choose_param;

    // basic geometric types
    typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
    typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;
    typedef typename Point_set_processing_3::GetPlaneMap<PlaneRange, NamedParameters>::type PlaneMap;
    typedef typename Point_set_processing_3::GetPlaneIndexMap<NamedParameters>::type PlaneIndexMap;

    CGAL_static_assertion_msg(!(boost::is_same<NormalMap,
                                typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value),
                              "Error: no normal map");
    CGAL_static_assertion_msg(!(boost::is_same<PlaneIndexMap,
                                typename Point_set_processing_3::GetPlaneIndexMap<NamedParameters>::NoMap>::value),
                              "Error: no plane index map");

    PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
    NormalMap normal_map = choose_param(get_param(np, internal_np::normal_map), NormalMap());
    PlaneMap plane_map = choose_param(get_param(np, internal_np::plane_map), PlaneMap());
    PlaneIndexMap index_map = choose_param(get_param(np, internal_np::plane_index_map), PlaneIndexMap());
    double attraction_factor = choose_param(get_param(np, internal_np::attraction_factor), 3.);
    
    m_points.reserve(points.size());
    m_normals.reserve(points.size());
    m_indices_of_assigned_points.resize (planes.size());

    m_indices.resize (points.size (), (std::numeric_limits<std::size_t>::max)());
    m_status.resize (points.size (), POINT);

    std::size_t idx = 0;
    for (typename PointRange::const_iterator it = points.begin();
         it != points.end(); ++ it)
    {
      m_points.push_back (get(point_map, *it));
      m_normals.push_back (get(normal_map, *it));
      int plane_index = get (index_map, idx);
      if (plane_index != -1)
      {
        m_indices_of_assigned_points[std::size_t(plane_index)].push_back(idx);
        m_indices[idx] = std::size_t(plane_index);
        m_status[idx] = PLANE;
      }
      ++ idx;
    }


    m_planes.reserve (planes.size());
    for (typename PlaneRange::const_iterator it = planes.begin();
         it != planes.end(); ++ it)
      m_planes.push_back (get (plane_map, *it));

    run (epsilon, attraction_factor);
    clean ();
  }
  /// \endcond

  std::size_t size () const { return m_points.size (); }
  std::pair<Point, Vector> operator[] (std::size_t i) const
  { return std::make_pair (m_points[i], m_normals[i]); }
  const Point& point (std::size_t i) const { return m_points[i]; }
  const Vector& normal (std::size_t i) const { return m_normals[i]; }

  /*!

    Returns all `Plane_shape` objects that are adjacent to the point
    with index `i`.

    \note Points not adjacent to any plane are free-form points,
    points adjacent to 1 plane are planar points, points adjacent to 2
    planes are edge points and points adjacent to 3 or more planes are
    vertices.

   */
  template <typename OutputIterator>
  void adjacency (std::size_t i, OutputIterator output) const
  {
    if (m_status[i] == PLANE || m_status[i] == RESIDUS)
      *(output ++) = m_planes[m_indices[i]];
    else if (m_status[i] == EDGE)
      {
        *(output ++) = m_planes[m_edges[m_indices[i]].planes[0]];
        *(output ++) = m_planes[m_edges[m_indices[i]].planes[1]];
      }
    else if (m_status[i] == CORNER)
      {
        for (std::size_t j = 0; j < m_corners[m_indices[i]].planes.size(); ++ j)
          *(output ++) = m_planes[m_corners[m_indices[i]].planes[j]];
      }
  }

  /*!

    Computes the coherence of a facet between the 3 points indexed by
    `f` with respect to the underlying structure.

   */
  Coherence_type facet_coherence (const CGAL::cpp11::array<std::size_t, 3>& f) const
  {
    // O- FREEFORM CASE
    if (m_status[f[0]] == POINT &&
        m_status[f[1]] == POINT &&
        m_status[f[2]] == POINT)
      return FREEFORM;
      
    // 1- PLANAR CASE
    if (m_status[f[0]] == PLANE &&
        m_status[f[1]] == PLANE &&
        m_status[f[2]] == PLANE)
      {
        if (m_indices[f[0]] == m_indices[f[1]] &&
            m_indices[f[0]] == m_indices[f[2]])
          return PLANAR;
        else
          return INCOHERENT;
      }

    for (std::size_t i = 0; i < 3; ++ i)
      {
        Point_status sa = m_status[f[(i+1)%3]];
        Point_status sb = m_status[f[(i+2)%3]];
        Point_status sc = m_status[f[(i+3)%3]];
        std::size_t a = m_indices[f[(i+1)%3]];
        std::size_t b = m_indices[f[(i+2)%3]];
        std::size_t c = m_indices[f[(i+3)%3]];

        // O- FREEFORM CASE
        if (sa == POINT && sb == POINT && sc == PLANE)
          return FREEFORM;
        if (sa == POINT && sb == PLANE && sc == PLANE)
          {
            if (b == c)
              return FREEFORM;
            else
              return INCOHERENT;
          }
          
        // 2- CREASE CASES
        if (sa == EDGE && sb == EDGE && sc == PLANE)
          {
            if ((c == m_edges[a].planes[0] ||
                 c == m_edges[a].planes[1]) &&
                (c == m_edges[b].planes[0] ||
                 c == m_edges[b].planes[1]))
              return CREASE;
            else
              return INCOHERENT;
          }

        if (sa == EDGE && sb == PLANE && sc == PLANE)
          {
            if (b == c &&
                (b == m_edges[a].planes[0] ||
                 b == m_edges[a].planes[1]))
              return CREASE;
            else
              return INCOHERENT;
          }


        // 3- CORNER CASES
        if (sc == CORNER)
          {
            if (sa == EDGE && sb == EDGE)
              {
                bool a0 = false, a1 = false, b0 = false, b1 = false;

                if ((m_edges[a].planes[0] != m_edges[b].planes[0] &&
                     m_edges[a].planes[0] != m_edges[b].planes[1] &&
                     m_edges[a].planes[1] != m_edges[b].planes[0] &&
                     m_edges[a].planes[1] != m_edges[b].planes[1]))
                  return INCOHERENT;
                  
                for (std::size_t j = 0; j < m_corners[c].planes.size (); ++ j)
                  {
                    if (m_corners[c].planes[j] == m_edges[a].planes[0])
                      a0 = true;
                    else if (m_corners[c].planes[j] == m_edges[a].planes[1])
                      a1 = true;
                    if (m_corners[c].planes[j] == m_edges[b].planes[0])
                      b0 = true;
                    else if (m_corners[c].planes[j] == m_edges[b].planes[1])
                      b1 = true;
                  }
                if (a0 && a1 && b0 && b1)
                  return VERTEX;
                else
                  return INCOHERENT;
              }
            else if (sa == PLANE && sb == PLANE)
              {
                if (a != b)
                  return INCOHERENT;

                for (std::size_t j = 0; j < m_corners[c].planes.size (); ++ j)
                  if (m_corners[c].planes[j] == a)
                    return VERTEX;
                  
                return INCOHERENT;
              }
            else if (sa == PLANE && sb == EDGE)
              {
                bool pa = false, b0 = false, b1 = false;
                if (a != m_edges[b].planes[0] && a != m_edges[b].planes[1])
                  return INCOHERENT;
                  
                for (std::size_t j = 0; j < m_corners[c].planes.size (); ++ j)
                  {
                    if (m_corners[c].planes[j] == a)
                      pa = true;
                    if (m_corners[c].planes[j] == m_edges[b].planes[0])
                      b0 = true;
                    else if (m_corners[c].planes[j] == m_edges[b].planes[1])
                      b1 = true;
                  }
                if (pa && b0 && b1)
                  return VERTEX;
                else
                  return INCOHERENT;
              }
            else if (sa == EDGE && sb == PLANE)
              {
                bool a0 = false, a1 = false, pb = false;
                if (b != m_edges[a].planes[0] && b != m_edges[a].planes[1])
                  return INCOHERENT;
                  
                for (std::size_t j = 0; j < m_corners[c].planes.size (); ++ j)
                  {
                    if (m_corners[c].planes[j] == b)
                      pb = true;
                    if (m_corners[c].planes[j] == m_edges[a].planes[0])
                      a0 = true;
                    else if (m_corners[c].planes[j] == m_edges[a].planes[1])
                      a1 = true;
                  }
                if (a0 && a1 && pb)
                  return VERTEX;
                else
                  return INCOHERENT;
              }
            else
              return INCOHERENT;
          }
      }


    return INCOHERENT;
  }


  /// \cond SKIP_IN_MANUAL  
private:


  void clean ()
  {
    std::vector<Point> points;
    std::vector<Vector> normals;
    std::vector<std::size_t> indices;
    std::vector<Point_status> status;
      
    for (std::size_t i = 0; i < m_points.size (); ++ i)
      if (m_status[i] != SKIPPED)
        {
          points.push_back (m_points[i]);
          normals.push_back (m_normals[i]);
          status.push_back (m_status[i]);
          if (m_status[i] == RESIDUS)
            status.back () = PLANE;
          indices.push_back (m_indices[i]);
        }
      
    m_points.swap (points);
    m_normals.swap (normals);
    m_indices.swap (indices);
    m_status.swap (status);
  }


  void run (double epsilon, double attraction_factor = 3.)
  {
    if (m_planes.empty ())
      return;
      
    double radius = epsilon * attraction_factor;

#ifdef CGAL_PSP3_VERBOSE
    std::cerr << "Computing planar points... " << std::endl;
#endif
      
    project_inliers ();
    resample_planes (epsilon);
      
#ifdef CGAL_PSP3_VERBOSE
    std::cerr << " -> Done" << std::endl;
    std::cerr << "Finding adjacent primitives... " << std::endl;
#endif
      
    find_pairs_of_adjacent_primitives (radius);

#ifdef CGAL_PSP3_VERBOSE
    std::cerr << " -> Found " << m_edges.size () << " pair(s) of adjacent primitives." << std::endl;
    std::cerr << "Computing edges... " << std::endl;
#endif
      
    compute_edges (epsilon);

#ifdef CGAL_PSP3_VERBOSE
    std::cerr << " -> Done" << std::endl;
    std::cerr << "Creating edge-anchor points... " << std::endl;
    {
      std::size_t size_before = m_points.size ();
#endif

      create_edge_anchor_points (radius, epsilon);

#ifdef CGAL_PSP3_VERBOSE
      std::cerr << " -> " << m_points.size () - size_before << " anchor point(s) created." << std::endl;
    }

    std::cerr << "Computating first set of corners... " << std::endl;
#endif
      
    compute_corners (radius);

#ifdef CGAL_PSP3_VERBOSE
    std::cerr << " -> Found " << m_corners.size () << " triple(s) of adjacent primitives/edges." << std::endl;
    std::cerr << "Merging corners... " << std::endl;
    {
      std::size_t size_before = m_points.size ();
#endif
        
      merge_corners (radius);

#ifdef CGAL_PSP3_VERBOSE
      std::cerr << " -> " << m_points.size () - size_before << " corner point(s) created." << std::endl;
    }

    std::cerr << "Computing corner directions... " << std::endl;
#endif
      
    compute_corner_directions (epsilon);

#ifdef CGAL_PSP3_VERBOSE
    std::cerr << " -> Done" << std::endl;
    std::cerr << "Refining sampling... " << std::endl;
#endif
      
    refine_sampling (epsilon);

#ifdef CGAL_PSP3_VERBOSE
    std::cerr << " -> Done" << std::endl;

    std::cerr << "Cleaning data set... " << std::endl;
#endif
      
    clean ();

#ifdef CGAL_PSP3_VERBOSE
    std::cerr << " -> Done" << std::endl;
#endif
  }

  void project_inliers ()
  {
    for(std::size_t i = 0; i < m_indices_of_assigned_points.size (); ++ i)
      for (std::size_t j = 0; j < m_indices_of_assigned_points[i].size(); ++ j)
        {
          std::size_t ind = m_indices_of_assigned_points[i][j];
          m_points[ind] = m_planes[i].projection (m_points[ind]);
        }
  }

  void resample_planes (double epsilon)
  {
    double grid_length = epsilon * (std::sqrt(2.) - 1e-3);

    for (std::size_t c = 0; c < m_planes.size (); ++ c)
      {
        //plane attributes and 2D projection vectors
        const Plane& plane = m_planes[c];
        Vector vortho = plane.orthogonal_vector();
        Vector b1 = plane.base1();
        Vector b2 = plane.base2();
			
        b1 = b1 / std::sqrt (b1 * b1);
        b2 = b2 / std::sqrt (b2 * b2);

        std::vector<Point_2> points_2d;

        //storage of the 2D points in "pt_2d"
        for (std::size_t j = 0; j < m_indices_of_assigned_points[c].size(); ++ j)
          {
            std::size_t ind = m_indices_of_assigned_points[c][j];
            const Point& pt = m_points[ind];
            points_2d.push_back (Point_2 (b1.x() * pt.x() + b1.y() * pt.y() + b1.z() * pt.z(),
                                          b2.x() * pt.x() + b2.y() * pt.y() + b2.z() * pt.z()));
          }


        //creation of a 2D-grid with cell width = grid_length, and image structures
        CGAL::Bbox_2 box_2d = CGAL::bbox_2 (points_2d.begin(), points_2d.end());
        std::size_t Nx = static_cast<std::size_t>((box_2d.xmax() - box_2d.xmin()) / grid_length) + 1;
        std::size_t Ny = static_cast<std::size_t>((box_2d.ymax() - box_2d.ymin()) / grid_length) + 1;
          
        std::vector<std::vector<bool> > Mask (Nx, std::vector<bool> (Ny, false));
        std::vector<std::vector<bool> > Mask_border (Nx, std::vector<bool> (Ny, false));
        std::vector<std::vector<std::vector<std::size_t> > >
          point_map (Nx, std::vector<std::vector<std::size_t> > (Ny, std::vector<std::size_t>()));

        //storage of the points in the 2D-grid "point_map"
        for (std::size_t i = 0; i < points_2d.size(); ++ i)
          {
            std::size_t ind_x = static_cast<std::size_t>((points_2d[i].x() - box_2d.xmin()) / grid_length);
            std::size_t ind_y = static_cast<std::size_t>((points_2d[i].y() - box_2d.ymin()) / grid_length);
            Mask[ind_x][ind_y] = true;
            point_map[ind_x][ind_y].push_back (m_indices_of_assigned_points[c][i]);
          }

        //hole filing in Mask in 4-connexity
        for (std::size_t j = 1; j < Ny - 1; ++ j)
          for (std::size_t i = 1; i < Nx - 1; ++ i)
            if( !Mask[i][j]
                && Mask[i-1][j] && Mask[i][j-1]
                && Mask[i][j+1] && Mask[i+1][j] )
              Mask[i][j]=true;
					
        //finding mask border in 8-connexity	
        for (std::size_t j = 1; j < Ny - 1; ++ j)
          for (std::size_t i = 1; i < Nx - 1; ++ i)
            if( Mask[i][j] &&
                ( !Mask[i-1][j-1] || !Mask[i-1][j] ||
                  !Mask[i-1][j+1] || !Mask[i][j-1] ||
                  !Mask[i][j+1] || !Mask[i+1][j-1] ||
                  !Mask[i+1][j]|| !Mask[i+1][j+1] ) )
              Mask_border[i][j]=true;
          
        for (std::size_t j = 0; j < Ny; ++ j)
          {
            if (Mask[0][j])
              Mask_border[0][j]=true;
            if (Mask[Nx-1][j])
              Mask_border[Nx-1][j]=true;
          }

        for (std::size_t i = 0; i < Nx; ++ i)
          {
            if(Mask[i][0])
              Mask_border[i][0]=true;
            if(Mask[i][Ny-1])
              Mask_border[i][Ny-1]=true;
          }

        //saving of points to keep
        for (std::size_t j = 0; j < Ny; ++ j)
          for (std::size_t i = 0; i < Nx; ++ i)
            if( point_map[i][j].size()>0)
              {
                //inside: recenter (cell center) the first point of the cell and desactivate the others points 
                if (!Mask_border[i][j] && Mask[i][j])
                  {
                    double x2pt = (i+0.5) * grid_length + box_2d.xmin();
                    double y2pt = (j+0.4) * grid_length + box_2d.ymin();
							
                    if (i%2 == 1)
                      {
                        x2pt = (i+0.5) * grid_length + box_2d.xmin();
                        y2pt = (j+0.6) * grid_length + box_2d.ymin();
                      }

                    FT X1 = x2pt * b1.x() + y2pt * b2.x() - plane.d() * vortho.x();
                    FT X2 = x2pt * b1.y() + y2pt * b2.y() - plane.d() * vortho.y();
                    FT X3 = x2pt * b1.z() + y2pt * b2.z() - plane.d() * vortho.z();

                    std::size_t index_pt = point_map[i][j][0];
                    m_points[index_pt] = Point (X1, X2, X3);
                    m_normals[index_pt] = m_planes[c].orthogonal_vector();
                    m_status[index_pt] = PLANE;

                    for (std::size_t np = 1; np < point_map[i][j].size(); ++ np)
                      m_status[point_map[i][j][np]] = SKIPPED;
                  }

                //border: recenter (barycenter) the first point of the cell and desactivate the others points
                else if (Mask_border[i][j] && Mask[i][j])
                  {
                    std::vector<Point> pts;
                    for (std::size_t np = 0; np < point_map[i][j].size(); ++ np)
                      pts.push_back (m_points[point_map[i][j][np]]);

                    m_points[point_map[i][j][0]] = CGAL::centroid (pts.begin (), pts.end ());
                    m_status[point_map[i][j][0]] = PLANE;
                    for (std::size_t np = 1; np < point_map[i][j].size(); ++ np)
                      m_status[point_map[i][j][np]] = SKIPPED;
                  }
              }
        // point use to filling 4-connexity holes are store in HPS_residus
            else if (point_map[i][j].size()==0 && !Mask_border[i][j] && Mask[i][j])
              {
                double x2pt = (i+0.5) * grid_length + box_2d.xmin();
                double y2pt = (j+0.49) * grid_length + box_2d.ymin();
                if(i%2==1)
                  {
                    x2pt = (i+0.5) * grid_length + box_2d.xmin();
                    y2pt = (j+0.51) * grid_length + box_2d.ymin();
                  }
                FT X1 = x2pt * b1.x() + y2pt * b2.x() - plane.d() * vortho.x();
                FT X2 = x2pt * b1.y() + y2pt * b2.y() - plane.d() * vortho.y();
                FT X3 = x2pt * b1.z() + y2pt * b2.z() - plane.d() * vortho.z();

                m_points.push_back (Point (X1, X2, X3));
                m_normals.push_back (m_planes[c].orthogonal_vector());
                m_indices.push_back (c);
                m_status.push_back (RESIDUS);
              }
      }

  }

  void find_pairs_of_adjacent_primitives (double radius)
  {
    typedef typename CGAL::Search_traits_3<Kernel> Search_traits_base;
    typedef Search_traits_adapter <std::size_t, typename Pointer_property_map<Point>::type, Search_traits_base> Search_traits;
    typedef CGAL::Kd_tree<Search_traits> Tree;
    typedef CGAL::Fuzzy_sphere<Search_traits> Fuzzy_sphere;

    typename Pointer_property_map<Point>::type pmap = make_property_map(m_points);

    Tree tree (boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t> (0),
               boost::counting_iterator<std::size_t, boost::use_default, std::ptrdiff_t> (m_points.size()),
               typename Tree::Splitter(),
               Search_traits (pmap));

    std::vector<std::vector<bool> > adjacency_table (m_planes.size (),
                                                     std::vector<bool> (m_planes.size (), false));

    //compute a basic adjacency relation (two primitives are neighbors
    //if at least one point of the primitive 1 is a k-nearest neighbor
    //of a point of the primitive 2 and vice versa)
    for (std::size_t i = 0; i < m_points.size (); ++ i)
      {
        std::size_t ind_i = m_indices[i];

        if (ind_i == (std::numeric_limits<std::size_t>::max)())
          continue;

        Fuzzy_sphere query (i, radius, 0., tree.traits());
          
        std::vector<std::size_t> neighbors;
        tree.search (std::back_inserter (neighbors), query);

          
        for (std::size_t k = 0; k < neighbors.size(); ++ k)
          {
            std::size_t ind_k = m_indices[neighbors[k]];
            if (ind_k != (std::numeric_limits<std::size_t>::max)() && ind_k != ind_i)
              adjacency_table[ind_i][ind_k] = true;
          }
      }

    //verify the symmetry and store the pairs of primitives in
    //m_edges
    for (std::size_t i = 0; i < adjacency_table.size() - 1; ++ i)
      for (std::size_t j = i + 1; j < adjacency_table[i].size(); ++ j)
        if ((adjacency_table[i][j]) && (adjacency_table[j][i]))
          m_edges.push_back (Edge (i, j));

  }

  void compute_edges (double epsilon)
  {
    for (std::size_t i = 0; i < m_edges.size(); ++ i)
      {
        const Plane& plane1 = m_planes[m_edges[i].planes[0]];
        const Plane& plane2 = m_planes[m_edges[i].planes[1]];       

        double angle_A = std::acos (CGAL::abs (plane1.orthogonal_vector() * plane2.orthogonal_vector()));
        double angle_B = CGAL_PI - angle_A;

        typename cpp11::result_of<typename Kernel::Intersect_3(Plane, Plane)>::type
          result = CGAL::intersection(plane1, plane2);

        if (!result)
          {
#ifdef CGAL_PSP3_VERBOSE
            std::cerr << "Warning: bad plane/plane intersection" << std::endl;
#endif
            continue;
          }

        if (const Line* l = boost::get<Line>(&*result))
          m_edges[i].support = *l;
        else
          {
#ifdef CGAL_PSP3_VERBOSE
            std::cerr << "Warning: bad plane/plane intersection" << std::endl;
#endif
            continue;
          }
        
        Vector direction_p1 (0., 0., 0.);
        for (std::size_t k = 0; k < m_indices_of_assigned_points[m_edges[i].planes[0]].size(); ++ k)
          {
            std::size_t index_point = m_indices_of_assigned_points[m_edges[i].planes[0]][k];
              
            const Point& point = m_points[index_point];
            Point projected = m_edges[i].support.projection (point);
            if (std::sqrt (CGAL::squared_distance (point, projected))
                < 2 * (std::min) (4., 1 / std::sin (angle_A)) * epsilon
                && m_status[index_point] != SKIPPED)
              direction_p1 = direction_p1 + Vector (projected, point);
          }
        if (direction_p1.squared_length() > 0)
          direction_p1 = direction_p1 / std::sqrt (direction_p1 * direction_p1);

        Vector direction_p2 (0., 0., 0.);
        for (std::size_t k = 0; k < m_indices_of_assigned_points[m_edges[i].planes[1]].size(); ++ k)
          {
            std::size_t index_point = m_indices_of_assigned_points[m_edges[i].planes[1]][k];
              
            const Point& point = m_points[index_point];
            Point projected = m_edges[i].support.projection (point);
            if (std::sqrt (CGAL::squared_distance (point, projected))
                < 2 * (std::min) (4., 1 / std::sin (angle_A)) * epsilon
                && m_status[index_point] != SKIPPED)
              direction_p2 = direction_p2 + Vector (projected, point);
          }
        if (direction_p2.squared_length() > 0)
          direction_p2 = direction_p2 / std::sqrt (direction_p2 * direction_p2);

        double angle = std::acos (direction_p1 * direction_p2);
      
        if (direction_p1.squared_length() == 0
            || direction_p2.squared_length() == 0
            || (CGAL::abs (angle - angle_A) > 1e-2
                && CGAL::abs (angle - angle_B) > 1e-2 ))
          {
            m_edges[i].active = false;
          }
      }
  }

  void create_edge_anchor_points (double radius, double epsilon)
  {
    double d_DeltaEdge = std::sqrt (2.) * epsilon;
    double r_edge = d_DeltaEdge / 2.;
      
    for (std::size_t i = 0; i < m_edges.size(); ++ i)
      {
        const Plane& plane1 = m_planes[m_edges[i].planes[0]];
        const Plane& plane2 = m_planes[m_edges[i].planes[1]];

        const Line& line = m_edges[i].support;

        if (!(m_edges[i].active))
          {
            continue;
          }

        Vector normal = 0.5 * plane1.orthogonal_vector () + 0.5 * plane2.orthogonal_vector();
							
        //find set of points close (<attraction_radius) to the edge and store in intersection_points
        std::vector<std::size_t> intersection_points;
        for (std::size_t k = 0; k < m_indices_of_assigned_points[m_edges[i].planes[0]].size(); ++ k)
          {
            std::size_t index_point = m_indices_of_assigned_points[m_edges[i].planes[0]][k];
            const Point& point = m_points[index_point];
            Point projected = line.projection (point);
            if (CGAL::squared_distance (point, projected) < radius * radius)
              intersection_points.push_back (index_point);
          }
        for (std::size_t k = 0; k < m_indices_of_assigned_points[m_edges[i].planes[1]].size(); ++ k)
          {
            std::size_t index_point = m_indices_of_assigned_points[m_edges[i].planes[1]][k];
            const Point& point = m_points[index_point];
            Point projected = line.projection (point);
            if (CGAL::squared_distance (point, projected) < radius * radius)
              intersection_points.push_back (index_point);
          }

        if (intersection_points.empty ())
          {
            continue;
          }

        const Point& t0 = m_points[intersection_points[0]];
        Point t0p = line.projection (t0);
        double dmin = 0.;
        double dmax = 0.;
        Point Pmin = t0p;
        Point Pmax = t0p;
        Vector dir = line.to_vector ();
          
        //compute the segment of the edge
        for (std::size_t k = 0; k < intersection_points.size(); ++ k)
          {
            std::size_t ind = intersection_points[k];
            const Point& point = m_points[ind];
            Point projected = line.projection (point);
            double d = Vector (t0p, projected) * dir;
                  
            if (d < dmin)
              {
                dmin = d;
                Pmin = projected;
              }
            else if (d > dmax)
              {
                dmax = d;
                Pmax = projected;
              }
          }

        // make a partition in a 1D image by voting if at the same
        // time at least one point of plane1 and one of point2 fall in
        // the same cell (same step as for planes)
        Segment seg (Pmin,Pmax);
        std::size_t number_of_division = static_cast<std::size_t>(std::sqrt (seg.squared_length ()) / d_DeltaEdge) + 1;
        std::vector<std::vector<std::size_t> > division_tab (number_of_division);

        for (std::size_t k = 0; k < intersection_points.size(); ++ k)
          {
            std::size_t ind = intersection_points[k];
            const Point& point = m_points[ind];
            Point projected = line.projection (point);

            std::size_t tab_index = static_cast<std::size_t>(std::sqrt (CGAL::squared_distance (seg[0], projected))
                                                             / d_DeltaEdge);

            division_tab[tab_index].push_back (ind);
          }

        //C1-CREATE the EDGE
        std::vector<int> index_of_edge_points;
        for (std::size_t j = 0; j < division_tab.size(); ++ j)
          {
            bool p1found = false, p2found = false;
            for (std::size_t k = 0; k < division_tab[j].size () && !(p1found && p2found); ++ k)
              {
                if (m_indices[division_tab[j][k]] == m_edges[i].planes[0])
                  p1found = true;
                if (m_indices[division_tab[j][k]] == m_edges[i].planes[1])
                  p2found = true;
              }

            if (!(p1found && p2found))
              {
                division_tab[j].clear();
                continue;
              }

            Point perfect (seg[0].x() + (seg[1].x() - seg[0].x()) * (j + 0.5) / double(number_of_division),
                           seg[0].y() + (seg[1].y() - seg[0].y()) * (j + 0.5) / double(number_of_division),
                           seg[0].z() + (seg[1].z() - seg[0].z()) * (j + 0.5) / double(number_of_division));

            // keep closest point, replace it by perfect one and skip the others
            double dist_min = (std::numeric_limits<double>::max)();
            std::size_t index_best = 0;

            for (std::size_t k = 0; k < division_tab[j].size(); ++ k)
              {
                std::size_t inde = division_tab[j][k];

                if (CGAL::squared_distance (line, m_points[inde]) < d_DeltaEdge * d_DeltaEdge)
                  m_status[inde] = SKIPPED; // Deactive points too close (except best, see below)
                  
                double distance = CGAL::squared_distance (perfect, m_points[inde]);
                if (distance < dist_min)
                  {
                    dist_min = distance;
                    index_best = inde;
                  }
              }

            m_points[index_best] = perfect;
            m_normals[index_best] = normal;
            m_status[index_best] = EDGE;
            m_indices[index_best] = i;
            m_edges[i].indices.push_back (index_best);

          }

        //C2-CREATE the ANCHOR
        Vector direction_p1(0,0,0);
        Vector direction_p2(0,0,0);

        for (std::size_t j = 0; j < division_tab.size() - 1; ++ j)
          {
            if (division_tab[j].empty () || division_tab[j+1].empty ())
              continue;
            Point anchor (seg[0].x() + (seg[1].x() - seg[0].x()) * (j + 1) / double(number_of_division),
                          seg[0].y() + (seg[1].y() - seg[0].y()) * (j + 1) / double(number_of_division),
                          seg[0].z() + (seg[1].z() - seg[0].z()) * (j + 1) / double(number_of_division));
              
            Plane ortho = seg.supporting_line().perpendicular_plane(anchor); 

            std::vector<Point> pts1, pts2;
            //Computation of the permanent angle and directions
            for (std::size_t k = 0; k < division_tab[j].size(); ++ k)
              { 
                std::size_t inde = division_tab[j][k];
                std::size_t plane = m_indices[inde];
                if (plane == m_edges[i].planes[0])
                  pts1.push_back (m_points[inde]);
                else if (plane == m_edges[i].planes[1])
                  pts2.push_back (m_points[inde]);
              }

            typename cpp11::result_of<typename Kernel::Intersect_3(Plane, Plane)>::type
              result = CGAL::intersection (plane1, ortho);
            if (result)
              {
                if (const Line* l = boost::get<Line>(&*result))
                  {
                    if (!(pts1.empty()))
                      {
                        Vector vecp1 = l->to_vector();
                        vecp1 = vecp1/ std::sqrt (vecp1 * vecp1);
                        Vector vtest1 (anchor, CGAL::centroid (pts1.begin (), pts1.end ()));
                        if (vtest1 * vecp1<0)
                          vecp1 = -vecp1;

                        direction_p1 = direction_p1+vecp1;

                        Point anchor1 = anchor + vecp1 * r_edge;
                        m_points.push_back (anchor1);
                        m_normals.push_back (m_planes[m_edges[i].planes[0]].orthogonal_vector());
                        m_indices.push_back (m_edges[i].planes[0]);
                        m_status.push_back (PLANE);
                      }
                  }
                else
                  {
#ifdef CGAL_PSP3_VERBOSE
                    std::cerr<<"Warning: bad plane/plane intersection"<<std::endl;
#endif
                  }
              }
            else
              {
#ifdef CGAL_PSP3_VERBOSE
                std::cerr<<"Warning: bad plane/plane intersection"<<std::endl;
#endif
              }

            
            result = CGAL::intersection (plane2,ortho);
            if (result)
              {
                if (const Line* l = boost::get<Line>(&*result))
                  {
                    if (!(pts2.empty()))
                      {
                        Vector vecp2 = l->to_vector();
                        vecp2 = vecp2 / std::sqrt (vecp2 * vecp2);
                        Vector vtest2 (anchor, CGAL::centroid (pts2.begin (), pts2.end ()));
                        if (vtest2 * vecp2 < 0)
                          vecp2 =- vecp2;

                        direction_p2 = direction_p2+vecp2;

                        Point anchor2 = anchor + vecp2 * r_edge;
                        m_points.push_back (anchor2);
                        m_normals.push_back (m_planes[m_edges[i].planes[1]].orthogonal_vector());
                        m_indices.push_back (m_edges[i].planes[1]);
                        m_status.push_back (PLANE);
                      }
                  }
                else
                  {
#ifdef CGAL_PSP3_VERBOSE
                    std::cerr<<"Warning: bad plane/plane intersection"<<std::endl;
#endif
                  }
              }
            else
              {
#ifdef CGAL_PSP3_VERBOSE
                std::cerr<<"Warning: bad plane/plane intersection"<<std::endl;
#endif
              }
          }
        
        //if not information enough (not enough edges to create
        //anchor) we unactivate the edge, else we update the angle
        //and directions
        if ( !(direction_p1.squared_length()>0 || direction_p2.squared_length()>0) )
          {
            m_edges[i].active = false;
            for (std::size_t j = 0; j < m_edges[i].indices.size (); ++ j)
              m_status[m_edges[i].indices[j]] = SKIPPED;
          }
      }
  }

  void compute_corners (double radius)
  {
    if (m_edges.size () < 3)
      return;

    std::vector<std::vector<std::size_t> > plane_edge_adj (m_planes.size());
    for (std::size_t i = 0; i < m_edges.size (); ++ i)
      if (m_edges[i].active)
        {
          plane_edge_adj[m_edges[i].planes[0]].push_back (i);
          plane_edge_adj[m_edges[i].planes[1]].push_back (i);
        }

    std::vector<std::set<std::size_t> > edge_adj (m_edges.size ());

    for (std::size_t i = 0; i < plane_edge_adj.size (); ++ i)
      {
        if (plane_edge_adj[i].size () < 2)
          continue;
          
        for (std::size_t j = 0; j < plane_edge_adj[i].size ()- 1; ++ j)
          for (std::size_t k = j + 1; k < plane_edge_adj[i].size (); ++ k)
            {
              edge_adj[plane_edge_adj[i][j]].insert (plane_edge_adj[i][k]);
              edge_adj[plane_edge_adj[i][k]].insert (plane_edge_adj[i][j]);
            }
      }

    for (std::size_t i = 0; i < edge_adj.size (); ++ i)
      {
        if (edge_adj[i].size () < 2)
          continue;

        std::set<std::size_t>::iterator end = edge_adj[i].end();
        end --;
          
        for (std::set<std::size_t>::iterator jit = edge_adj[i].begin ();
             jit != end; ++ jit)
          {
            std::size_t j = *jit;
            if (j < i)
              continue;

            std::set<std::size_t>::iterator begin = jit;
            begin ++;
            for (std::set<std::size_t>::iterator kit = begin;
                 kit != edge_adj[i].end (); ++ kit)
              {
                std::size_t k = *kit;
                if (k < j)
                  continue;
                  
                std::set<std::size_t> planes;
                planes.insert (m_edges[i].planes[0]);
                planes.insert (m_edges[i].planes[1]);
                planes.insert (m_edges[j].planes[0]);
                planes.insert (m_edges[j].planes[1]);
                planes.insert (m_edges[k].planes[0]);
                planes.insert (m_edges[k].planes[1]);

                if (planes.size () == 3) // Triple found
                  {
                    std::vector<std::size_t> vecplanes (planes.begin (), planes.end ());
                    m_corners.push_back (Corner (vecplanes[0], vecplanes[1], vecplanes[2],
                                                 i, j, k));
                  }
              }
          }
      }


    for (std::size_t i = 0; i < m_corners.size (); ++ i)
      {
        //calcul pt d'intersection des 3 plans
        const Plane& plane1 = m_planes[m_corners[i].planes[0]];
        const Plane& plane2 = m_planes[m_corners[i].planes[1]];
        const Plane& plane3 = m_planes[m_corners[i].planes[2]];

        typename cpp11::result_of<typename Kernel::Intersect_3(Plane, Plane)>::type
          result = CGAL::intersection(plane1, plane2);
        
        if (result)
          {
            if (const Line* l = boost::get<Line>(&*result))
              {
                typename cpp11::result_of<typename Kernel::Intersect_3(Line, Plane)>::type
                  result2 = CGAL::intersection(*l, plane3);

                if (result2)
                  {
                    if (const Point* p = boost::get<Point>(&*result2))
                      m_corners[i].support = *p;
                    else
                      {
#ifdef CGAL_PSP3_VERBOSE
                        std::cerr << "Warning: bad plane/line intersection" << std::endl;
#endif
                        m_corners[i].active = false;
                        continue;
                      }
                  }
                else
                  {
#ifdef CGAL_PSP3_VERBOSE
                    std::cerr << "Warning: bad plane/line intersection" << std::endl;
#endif
                    m_corners[i].active = false;
                    continue;

                  }
              }
            else
              {
#ifdef CGAL_PSP3_VERBOSE
                std::cerr << "Warning: bad plane/plane intersection" << std::endl;
#endif
                m_corners[i].active = false;
                continue;
              }
          }
        else
          {
#ifdef CGAL_PSP3_VERBOSE
            std::cerr << "Warning: bad plane/plane intersection" << std::endl;
#endif
            m_corners[i].active = false;
            continue;
          }

        // test if point is in bbox + delta
        CGAL::Bbox_3 bbox = CGAL::bbox_3 (m_points.begin (), m_points.end ());
          
        double margin_x = 0.1 * (bbox.xmax() - bbox.xmin());
        double X_min = bbox.xmin() - margin_x;
        double X_max = bbox.xmax() + margin_x; 
        double margin_y = 0.1 * (bbox.ymax() - bbox.ymin());
        double Y_min = bbox.ymin() - margin_y;
        double Y_max = bbox.ymax() + margin_y; 
        double margin_z = 0.1* (bbox.zmax() - bbox.zmin());
        double Z_min = bbox.zmin() - margin_z;
        double Z_max = bbox.zmax() + margin_z;
          
        if ((m_corners[i].support.x() < X_min) || (m_corners[i].support.x() > X_max)
            || (m_corners[i].support.y() < Y_min) || (m_corners[i].support.y() > Y_max)
            || (m_corners[i].support.z() < Z_min) || (m_corners[i].support.z() > Z_max))
          {
            m_corners[i].active = false;
            continue;
          }

        // test if corner is in neighborhood of at least one point each of the 3 planes
        std::vector<bool> neighborhood (3, false);

        for (std::size_t k = 0; k < 3; ++ k)
          {
            for (std::size_t j = 0; j < m_edges[m_corners[i].edges[k]].indices.size(); ++ j)
              {
                const Point& p = m_points[m_edges[m_corners[i].edges[k]].indices[j]];

                if (CGAL::squared_distance (m_corners[i].support, p) < radius * radius)
                  {
                    neighborhood[k] = true;
                    break;
                  }
              }
          }

        if ( !(neighborhood[0] && neighborhood[1] && neighborhood[2]) )
          m_corners[i].active = false;
      }
  }

  void merge_corners (double radius)
  {
    for (std::size_t k = 0; k < m_corners.size(); ++ k)
      {
        if (!(m_corners[k].active))
          continue;

        int count_plane_number=3;
        
        for (std::size_t kb = k + 1; kb < m_corners.size(); ++ kb)
          {
            if (!(m_corners[kb].active))
              continue;

            int count_new_plane = 0;

            if (CGAL::squared_distance (m_corners[kb].support, m_corners[k].support) >= radius * radius)
              continue;

            for (std::size_t i = 0; i < m_corners[kb].planes.size (); ++ i)
              {
                bool testtt = true; 
                for (std::size_t l = 0; l < m_corners[k].planes.size(); ++ l)
                  if (m_corners[kb].planes[i] == m_corners[k].planes[l])
                    {
                      testtt = false;
                      break;
                    }
                if (!testtt)
                  continue;

                m_corners[k].planes.push_back (m_corners[kb].planes[i]);
                ++ count_new_plane;
                m_corners[kb].active = false;

                std::vector<bool> is_edge_in (3, false);
                for (std::size_t l = 0; l < m_corners[k].edges.size(); ++ l)
                  {
                    for (std::size_t j = 0; j < 3; ++ j)
                      if (m_corners[k].edges[l] == m_corners[kb].edges[j])
                        is_edge_in[j] = true;
                  }
                for (std::size_t j = 0; j < 3; ++ j)
                  if (!(is_edge_in[j]))
                    m_corners[k].edges.push_back (m_corners[kb].edges[j]);

              }
              
            //update barycenter
            m_corners[k].support = CGAL::barycenter (m_corners[k].support, count_plane_number,
                                                     m_corners[kb].support, count_new_plane);
            count_plane_number += count_new_plane;
          }

        // Compute normal vector
        Vector normal (0., 0., 0.);
        for (std::size_t i = 0; i < m_corners[k].planes.size(); ++ i)
          normal = normal + (1. / (double)(m_corners[k].planes.size()))
            * m_planes[m_corners[k].planes[i]].orthogonal_vector();
          
        m_points.push_back (m_corners[k].support);
        m_normals.push_back (normal);
        m_indices.push_back (k);
        m_status.push_back (CORNER);
      }
  }

  void compute_corner_directions (double epsilon)
  {
    for (std::size_t k = 0; k < m_corners.size(); ++ k)
      {
        for (std::size_t ed = 0; ed < m_corners[k].edges.size(); ++ ed)
          {
            if (m_corners[k].edges[ed] < m_edges.size())
              {  
                const Edge& edge = m_edges[m_corners[k].edges[ed]];

                Vector direction (0., 0., 0.);
                for (std::size_t i = 0; i < edge.indices.size(); ++ i)
                  {
                    std::size_t index_pt = edge.indices[i];
                    if (std::sqrt (CGAL::squared_distance (m_corners[k].support,
                                                           m_points[index_pt])) < 5 * epsilon)
                      direction = direction + Vector (m_corners[k].support, m_points[index_pt]);
                  }

                if (direction.squared_length() > 1e-5)
                  m_corners[k].directions.push_back (direction / std::sqrt (direction * direction));
                else
                  m_corners[k].directions.push_back (Vector (0., 0., 0.));
              }
            else
              m_corners[k].directions.push_back (Vector (0., 0., 0.));
          }
      }
  }
    
  void refine_sampling (double epsilon)
  {
    double d_DeltaEdge = std::sqrt (2.) * epsilon;

    for (std::size_t k = 0; k < m_corners.size(); ++ k)
      {
        if (!(m_corners[k].active))
          continue;
          
        for (std::size_t ed = 0; ed < m_corners[k].edges.size(); ++ ed)
          {
            const Edge& edge = m_edges[m_corners[k].edges[ed]];

            for (std::size_t i = 0; i < edge.indices.size(); ++ i)
              {
                //if too close from a corner, ->remove
                if (CGAL::squared_distance (m_corners[k].support, m_points[edge.indices[i]])
                    < d_DeltaEdge * d_DeltaEdge)
                  m_status[edge.indices[i]] = SKIPPED;
				
                //if too close from a corner (non dominant side), ->remove
                if (m_corners[k].directions[ed].squared_length() > 0
                    && (m_corners[k].directions[ed]
                        * Vector (m_corners[k].support, m_points[edge.indices[i]]) < 0)
                    && (CGAL::squared_distance (m_corners[k].support, m_points[edge.indices[i]])
                        < 4 * d_DeltaEdge * d_DeltaEdge))
                  m_status[edge.indices[i]] = SKIPPED;
              }
              
          }
      }

    for (std::size_t k = 0; k < m_corners.size(); ++ k)
      {
        if (!(m_corners[k].active))
          continue;
		
        for (std::size_t ed = 0; ed < m_corners[k].edges.size(); ++ ed)
          {
            if (m_corners[k].directions[ed].squared_length() <= 0.)
              continue;
              
            Edge& edge = m_edges[m_corners[k].edges[ed]];

            //rajouter un edge a epsilon du cote dominant si pas de point entre SS_edge/2 et 3/2*SS_edge
            bool is_in_interval = false;
            for (std::size_t i = 0; i < edge.indices.size(); ++ i)
              {
                std::size_t index_pt = edge.indices[i];
                double dist = CGAL::squared_distance (m_corners[k].support,
                                                      m_points[index_pt]);
                if (m_status[index_pt] != SKIPPED
                    && dist < 1.5 * d_DeltaEdge && dist > d_DeltaEdge / 2)
                  {
                    Vector move (m_corners[k].support,
                                 m_points[index_pt]);
                    if (move * m_corners[k].directions[ed] > 0.)
                      {
                        is_in_interval = true;
                        break;
                      }
                  }
              }

            //rajouter un edge a 1 epsilon du cote dominant si pas de point entre SS_edge/2 et 3/2*SS_edge
            if (!is_in_interval)
              {
                Point new_edge = m_corners[k].support + m_corners[k].directions[ed] * d_DeltaEdge;
                m_points.push_back (new_edge);
                m_normals.push_back (0.5 * m_planes[m_edges[m_corners[k].edges[ed]].planes[0]].orthogonal_vector()
                                     + 0.5 * m_planes[m_edges[m_corners[k].edges[ed]].planes[1]].orthogonal_vector());
                m_status.push_back (EDGE);
                m_indices.push_back (m_corners[k].edges[ed]);
                edge.indices.push_back (m_points.size() - 1);
              }
						
            //rajouter un edge a 1/3 epsilon du cote dominant
            Point new_edge = m_corners[k].support + m_corners[k].directions[ed] * d_DeltaEdge / 3;
            m_points.push_back (new_edge);
            m_normals.push_back (0.5 * m_planes[m_edges[m_corners[k].edges[ed]].planes[0]].orthogonal_vector()
                                 + 0.5 * m_planes[m_edges[m_corners[k].edges[ed]].planes[1]].orthogonal_vector());
            m_status.push_back (EDGE);
            m_indices.push_back (m_corners[k].edges[ed]);
            edge.indices.push_back (m_points.size() - 1);
          }
      }

  }
  /// \endcond    
};


  


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/** 
   \ingroup PkgPointSetProcessingAlgorithms

   This is an implementation of the Point Set Structuring algorithm. This
   algorithm takes advantage of a set of detected planes: it detects adjacency
   relationships between planes and resamples the detected planes, edges and
   corners to produce a structured point set.

   The size parameter `epsilon` is used both for detecting adjacencies and for
   setting the sampling density of the structured point set.

   For more details, please refer to \cgalCite{cgal:la-srpss-13}.

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam PlaneRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `plane_map`.
   \tparam OutputIterator Type of the output iterator. The type of the
   objects put in it is `std::pair<Kernel::Point_3, Kernel::Vector_3>`.
   Note that the user may use a
   <A HREF="http://www.boost.org/libs/iterator/doc/function_output_iterator.html">function_output_iterator</A>
   to match specific needs.

   \param points input point range.
   \param planes input plane range.
   \param output output iterator where output points are written
   \param epsilon size parameter.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadablePropertyMap` with value type
     `geom_traits::Vector_3`.\cgalParamEnd
     \cgalParamBegin{plane_index_map} a model of `ReadablePropertyMap` with value type `int`.
     Associates the index of a point in the input range to the index of plane (-1 if point does is not assigned to
     a plane).\cgalParamEnd
     \cgalParamBegin{plane_map} a model of `ReadablePropertyMap` with value type
     `geom_traits::Plane_3`. If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Plane_3>`
     is used.\cgalParamEnd
     \cgalParamBegin{attraction_factor} multiple of `epsilon` used to connect simplices.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

*/
template <typename PointRange,
          typename PlaneRange,
          typename OutputIterator,
          typename NamedParameters
          >
OutputIterator
structure_point_set (const PointRange& points,
                     const PlaneRange& planes,
                     OutputIterator output,
                     double epsilon,
                     const NamedParameters& np)
{
  using boost::choose_param;

  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  Point_set_with_structure<Kernel> pss (points, planes, epsilon, np);

  for (std::size_t i = 0; i < pss.size(); ++ i)
    *(output ++) = pss[i];

  return output;
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange,
          typename PlaneRange,
          typename OutputIterator>
OutputIterator
structure_point_set (const PointRange& points, ///< range of points.
                     const PlaneRange& planes, ///< range of planes.
                     OutputIterator output, ///< output iterator where output points are written.
                     double epsilon) ///< size parameter.
{
  return structure_point_set
    (points, planes, output, epsilon,
     CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond

#ifndef CGAL_NO_DEPRECATED_CODE
/// \cond SKIP_IN_MANUAL

namespace Shape_detection_3{
//Forward declarations
template <class Traits>
class Efficient_RANSAC;
template <typename Traits>
class Plane_map;
template <typename Traits>
class Point_to_shape_index_map;
} // end of namespace Shape_detection_3

// deprecated API
template <typename Traits,
          typename OutputIterator
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::structure_point_set(), please update your code")
OutputIterator
structure_point_set (typename Traits::Input_range::iterator first,  ///< iterator over the first input point.
                     typename Traits::Input_range::iterator beyond, ///< past-the-end iterator over the input points.
                     typename Traits::Point_map point_map, ///< property map: value_type of InputIterator -> Point_3.
                     typename Traits::Normal_map normal_map, ///< property map: value_type of InputIterator -> Vector_3.
                     OutputIterator output, ///< output iterator where output points are written
                     Shape_detection_3::Efficient_RANSAC<Traits>&
                     shape_detection, ///< shape detection object
                     double epsilon, ///< size parameter
                     double attraction_factor = 3.) ///< attraction factor
{
  typename Shape_detection_3::Efficient_RANSAC<Traits>::Plane_range planes = shape_detection.planes();
  return structure_point_set (CGAL::make_range(first, beyond),
                              planes,
                              output,
                              epsilon, // epsilon for structuring points
                              CGAL::parameters::point_map (point_map).
                              normal_map (normal_map).
                              plane_map (CGAL::Shape_detection_3::Plane_map<Traits>()).
                              plane_index_map (Shape_detection_3::Point_to_shape_index_map<Traits>(CGAL::make_range(first, beyond), planes)).
                              attraction_factor(attraction_factor));
}

// deprecated API
template <typename Traits,
          typename OutputIterator
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::structure_point_set(), please update your code")
OutputIterator
structure_point_set (typename Traits::Input_range::iterator first,  ///< iterator over the first input point.
                     typename Traits::Input_range::iterator beyond, ///< past-the-end iterator over the input points.
                     OutputIterator output, ///< output iterator where output points are written
                     Shape_detection_3::Efficient_RANSAC<Traits>&
                     shape_detection, ///< shape detection object
                     double epsilon, ///< size parameter
                     double attraction_factor = 3.) ///< attraction factor
{
  return structure_point_set (first, beyond,
                              typename Traits::Point_map(),
                              typename Traits::Normal_map(),
                              output,
                              shape_detection,
                              epsilon,
                              attraction_factor);
}
/// \endcond
#endif // CGAL_NO_DEPRECATED_CODE


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_STRUCTURE_POINT_SET_3_H

