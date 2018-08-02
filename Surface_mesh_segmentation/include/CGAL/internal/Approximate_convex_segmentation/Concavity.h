// Copyright (c) 2018  GeometryFactory Sarl (France).
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
// Author(s)     : Liubomyr Piadyk

#ifndef CGAL_INTERNAL_CONCAVITY_H
#define CGAL_INTERNAL_CONCAVITY_H

#include <CGAL/license/Surface_mesh_segmentation.h>

#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for_each.h>
#include <tbb/mutex.h>
#endif

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

namespace CGAL
{
namespace internal
{

  template <class TriangleMesh, class Vpm, class GeomTraits, class ConcurrencyTag, class Mesh = TriangleMesh>
  class Concavity
  {
    // predefined structs
    struct Intersection_functor;

    // typedefs
    typedef typename GeomTraits::Point_3 Point_3;
    typedef typename GeomTraits::Vector_3 Vector_3;
    typedef typename GeomTraits::Ray_3 Ray_3;
    
    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

    typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator vertex_iterator;
    
    typedef std::map<vertex_descriptor, Vector_3> Normals_map;

    typedef CGAL::Face_filtered_graph<TriangleMesh> Filtered_graph;

    typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> AABB_primitive;
    typedef CGAL::AABB_tree<CGAL::AABB_traits<GeomTraits, AABB_primitive> > AABB_tree;

    typedef boost::optional<typename AABB_tree::template Intersection_and_primitive_id<Ray_3>::Type> Ray_intersection;
  
  public:
    Concavity(const TriangleMesh& mesh, const Vpm& vpm, const GeomTraits& traits, bool shortest_method = false)
    : m_mesh(mesh)
    , m_vpm(vpm)
    , m_traits(traits)
    , m_shortest_method(shortest_method)
    , m_normals_computed(false)
    {}
    
    Concavity(const TriangleMesh& mesh, const Vpm& vpm, const GeomTraits& traits, bool shortest_method, const Normals_map& normals_map)
    : m_mesh(mesh)
    , m_vpm(vpm)
    , m_traits(traits)
    , m_shortest_method(shortest_method)
    , m_normals_map(normals_map)
    , m_normals_computed(true)
    {}

    /**
    * Computes concavity value of a segment of the mesh which id is specified.
    */
    template <class FacetPropertyMap, class DistancesMap>
    double compute(FacetPropertyMap facet_ids, std::size_t segment_id, DistancesMap& distances)
    {
      compute_normals();
      
      Filtered_graph filtered_mesh(m_mesh, segment_id, facet_ids);

      Concavity<Filtered_graph, Vpm, GeomTraits, ConcurrencyTag, TriangleMesh> concavity(filtered_mesh, m_vpm, m_traits, m_shortest_method, m_normals_map);
      return concavity.compute(distances);
    }

    /**
    * Computes concavity value of the whole mesh along with the distance for each vertex.
    */
    template <class DistancesMap>
    double compute(DistancesMap& distances)
    {
      CGAL_assertion(!CGAL::is_empty(m_mesh));

      Mesh conv_hull;
      std::vector<Point_3> pts;

      if (num_vertices(m_mesh) <= 4) return 0;

      // extract the list points of the mesh
      BOOST_FOREACH(vertex_descriptor vert, vertices(m_mesh))
      {
        pts.push_back(get(m_vpm, vert));
      }

      // compute convex hull
      CGAL::convex_hull_3(pts.begin(), pts.end(), conv_hull); 
      
      if (m_shortest_method)
      {
        return compute_shortest(vertices(m_mesh), conv_hull, distances);
      }
      return compute_projected(vertices(m_mesh), conv_hull, distances);
    }
    
    double compute()
    {
      boost::unordered_map<vertex_descriptor, double> distances;
      return compute(boost::make_assoc_property_map(distances));
    }

    /**
    * Constructs list of vertices from the list of faces and computes concavity value with the convex hull provided.
    * Faces list is a subset of all faces in the mesh.
    */
    template <class DistancesMap>
    double compute(const std::vector<face_descriptor>& faces, const Mesh& conv_hull, DistancesMap distances)
    {
      boost::unordered_set<vertex_descriptor> pts;

      BOOST_FOREACH(face_descriptor face, faces)
      {
        BOOST_FOREACH(vertex_descriptor vert, vertices_around_face(halfedge(face, m_mesh), m_mesh))
        {
          pts.insert(vert);
        }
      }

      if (pts.size() <= 4) return 0;

      if (m_shortest_method)
      {
        return compute_shortest(std::make_pair(pts.begin(), pts.end()), conv_hull, distances);
      }
      return compute_projected(std::make_pair(pts.begin(), pts.end()), conv_hull, distances);
    }
    
    double compute(const std::vector<face_descriptor>& faces, const Mesh& conv_hull)
    {
      boost::unordered_map<vertex_descriptor, double> distances;
      return compute(faces, conv_hull, boost::make_assoc_property_map(distances));
    }

    void compute_normals()
    {
      if (m_normals_computed) return; // if the normals are already computed, then skip

      CGAL::Polygon_mesh_processing::compute_vertex_normals(m_mesh, boost::associative_property_map<Normals_map>(m_normals_map));
      m_normals_computed = true;
    }

  private:
    const TriangleMesh& m_mesh;
    Vpm m_vpm;
    const GeomTraits& m_traits;

    bool m_shortest_method;

    Normals_map m_normals_map;
    bool m_normals_computed;
    
    /**
    * Computes concavity value projecting vertices from a list onto a convex hull.
    */
    template <class iterator, class DistancesMap>
    double compute_projected(const std::pair<iterator, iterator>& verts, const Mesh& conv_hull, DistancesMap distances)
    {
      // compute normals if normals are not computed
      compute_normals();

      // construct AABB for fast computations of intersections between ray and convex hull
      AABB_tree tree(faces(conv_hull).begin(), faces(conv_hull).end(), conv_hull);

      // compute intersections and select the largest projection length
      double result = 0;

      // functor that computes intersection, its projection length from a vertex and maximizes the result variable
      struct Intersection_functor
      {
        Intersection_functor(const TriangleMesh& mesh, const Vpm& vpm, const Normals_map& normals_map, const AABB_tree& tree, DistancesMap& distances, double& result
#ifdef CGAL_LINKED_WITH_TBB
        , tbb::mutex& mutex
#endif
        )
        : m_mesh(mesh), m_vpm(vpm), m_normals_map(normals_map), m_tree(tree), m_distances(distances), m_result(result)
#ifdef CGAL_LINKED_WITH_TBB
       , m_mutex(mutex)
#endif
       {}

        void operator() (const vertex_descriptor& vert) const
        {
          Point_3 origin = get(m_vpm, vert);
          Ray_3 ray(origin, m_normals_map.at(vert));
          
          Ray_intersection intersection = m_tree.first_intersection(ray);

#ifdef CGAL_LINKED_WITH_TBB
          m_mutex.lock();
#endif
          put(m_distances, vert, 0.);
          if (intersection)
          {
            const Point_3* intersection_point =  boost::get<Point_3>(&(intersection->first));
            if (intersection_point)
            {
              double d = CGAL::squared_distance(origin, *intersection_point);
              put(m_distances, vert, d);
              m_result = std::max(m_result, d);
            }
          }
#ifdef CGAL_LINKED_WITH_TBB
          m_mutex.unlock();
#endif
        }

      private:
        const TriangleMesh& m_mesh;
        const Vpm& m_vpm;
        const Normals_map& m_normals_map;
        const AABB_tree& m_tree;
        DistancesMap& m_distances;
        double& m_result;
#ifdef CGAL_LINKED_WITH_TBB
        tbb::mutex& m_mutex;
#endif
      };

#ifdef CGAL_LINKED_WITH_TBB
      tbb::mutex mutex;
#endif
      Intersection_functor intersection_functor(m_mesh, m_vpm, m_normals_map, tree, distances, result
#ifdef CGAL_LINKED_WITH_TBB
      , mutex
#endif
      );

#ifdef CGAL_LINKED_WITH_TBB
      if (boost::is_convertible<ConcurrencyTag, Parallel_tag>::value)
      {
        tbb::parallel_for_each(verts.first, verts.second, intersection_functor);
      }
      else
#endif
      {
        BOOST_FOREACH(vertex_descriptor vert, verts)
        {
          intersection_functor(vert);
        }
      }

      return CGAL::sqrt(result);
    }
    
    /**
    * Computes concavity values as the shortest distances from the point of a vertex on the convex hull.
    */
    template <class iterator, class DistancesMap>
    double compute_shortest(const std::pair<iterator, iterator>& verts, const Mesh& conv_hull, DistancesMap distances)
    {
      // construct AABB for fast the shortest distance computations 
      AABB_tree tree(faces(conv_hull).begin(), faces(conv_hull).end(), conv_hull);
      tree.build();
      tree.accelerate_distance_queries();

      // compute intersections and select the largest projection length
      double result = 0;

      // functor that computes the shortest distance and maximizes the result variable
      struct Shortest_distance_functor
      {
        Shortest_distance_functor(const Vpm& vpm, const AABB_tree& tree, DistancesMap& distances, double& result
#ifdef CGAL_LINKED_WITH_TBB
        , tbb::mutex& mutex
#endif
        )
        : m_vpm(vpm), m_tree(tree), m_distances(distances), m_result(result)
#ifdef CGAL_LINKED_WITH_TBB
       , m_mutex(mutex)
#endif
       {}

        void operator() (const vertex_descriptor& vert) const
        {
          double distance = static_cast<double>(m_tree.squared_distance(get(m_vpm, vert)));

#ifdef CGAL_LINKED_WITH_TBB
          m_mutex.lock();
#endif
          put(m_distances, vert, distance);
          m_result = std::max(m_result, distance);
#ifdef CGAL_LINKED_WITH_TBB
          m_mutex.unlock();
#endif
        }

      private:
        const Vpm& m_vpm;
        const AABB_tree& m_tree;
        DistancesMap& m_distances;
        double& m_result;
#ifdef CGAL_LINKED_WITH_TBB
        tbb::mutex& m_mutex;
#endif
      };

#ifdef CGAL_LINKED_WITH_TBB
      tbb::mutex mutex;
#endif
      Shortest_distance_functor shortest_functor(m_vpm, tree, distances, result
#ifdef CGAL_LINKED_WITH_TBB
      , mutex
#endif
      );

#ifdef CGAL_LINKED_WITH_TBB
      if (boost::is_convertible<ConcurrencyTag, Parallel_tag>::value)
      {
        tbb::parallel_for_each(verts.first, verts.second, shortest_functor);
      }
      else
#endif
      {
        BOOST_FOREACH(vertex_descriptor vert, verts)
        {
          shortest_functor(vert);
        }
      }

      return CGAL::sqrt(result);
    }
  };
}
}

#endif // CGAL_INTERNAL_CONCAVITY_H
