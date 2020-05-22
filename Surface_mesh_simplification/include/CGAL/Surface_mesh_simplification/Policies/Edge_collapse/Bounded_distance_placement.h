// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno,
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/assertions.h>
#include <CGAL/Default.h>
#include <CGAL/intersections.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/optional.hpp>

#include <vector>
#include <type_traits>

namespace CGAL {
namespace Surface_mesh_simplification {

// An AABB tree can also be passed to the placement instead, see undocumented specialization below
template<typename BasePlacement, typename GeomTraits>
class Bounded_distance_placement
{
  typedef GeomTraits                                                          Geom_traits;
  typedef typename Geom_traits::FT                                            FT;

  typedef typename GeomTraits::Triangle_3                                     Triangle;
  typedef std::vector<Triangle>                                               Triangle_container;
  typedef typename Triangle_container::iterator                               TC_iterator;

  typedef CGAL::AABB_triangle_primitive<GeomTraits, TC_iterator>              Primitive;
  typedef CGAL::AABB_traits<GeomTraits, Primitive>                            Traits;
  typedef CGAL::AABB_tree<Traits>                                             AABB_tree;

private:
  template <typename Profile>
  void initialize_tree(const Profile& profile) const
  {
    CGAL_static_assertion((std::is_same<GeomTraits, typename Profile::Geom_traits>::value));

    typedef typename Profile::Triangle_mesh                                   Triangle_mesh;
    typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
    typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

    CGAL_precondition(m_tree_ptr == nullptr);

    const Triangle_mesh& tm = profile.surface_mesh();
    const Geom_traits& gt = profile.geom_traits();

    m_input_triangles.reserve(faces(tm).size());
    for(face_descriptor f : faces(profile.surface_mesh()))
    {
      halfedge_descriptor h = halfedge(f, tm);
      CGAL_assertion(!is_border(h, tm));

      m_input_triangles.emplace_back(gt.construct_triangle_3_object()(
                                       get(profile.vertex_point_map(), source(h, tm)),
                                       get(profile.vertex_point_map(), target(h, tm)),
                                       get(profile.vertex_point_map(), target(next(h, tm), tm))));
    }

    m_tree_ptr = new AABB_tree(m_input_triangles.begin(), m_input_triangles.end());
    const_cast<AABB_tree*>(m_tree_ptr)->build();
    const_cast<AABB_tree*>(m_tree_ptr)->accelerate_distance_queries();
  }

public:
  Bounded_distance_placement(const FT dist,
                             const BasePlacement& placement = BasePlacement())
    :
      m_sq_threshold_dist(CGAL::square(dist)),
      m_tree_ptr(nullptr),
      m_base_placement(placement)
  { }

  ~Bounded_distance_placement()
  {
    if(m_tree_ptr != nullptr)
      delete m_tree_ptr;
  }

  template <typename Profile>
  boost::optional<typename Profile::Point>
  operator()(const Profile& profile) const
  {
    typedef typename Profile::Point                                           Point;

    boost::optional<typename Profile::Point> op = m_base_placement(profile);
    if(op)
    {
      if(m_tree_ptr == nullptr)
        initialize_tree(profile);

      CGAL_assertion(m_tree_ptr != nullptr);
      CGAL_assertion(!m_tree_ptr->empty());

      const Point& p = *op;

      const Point& cp = m_tree_ptr->best_hint(p).first;

      // We could do better by having access to the internal kd-tree
      // and call search_any_point with a fuzzy_sphere.

      // if no input vertex is closer than the threshold, then
      // any face closer than the threshold is intersected by
      // the sphere (avoid the inclusion of the mesh into the threshold sphere)
      if(CGAL::compare_squared_distance(p, cp, m_sq_threshold_dist) != LARGER ||
         m_tree_ptr->do_intersect(CGAL::Sphere_3<Geom_traits>(p, m_sq_threshold_dist)))
        return op;

      return boost::optional<Point>();
    }

    return op;
  }

private:
  const FT m_sq_threshold_dist;
  mutable const AABB_tree* m_tree_ptr;
  mutable std::vector<Triangle> m_input_triangles;

  const BasePlacement m_base_placement;
};

// Undocumented specizalization where an _already built_ AABB tree is passed
template<typename BasePlacement, typename AABBTraits>
class Bounded_distance_placement<BasePlacement, CGAL::AABB_tree<AABBTraits> >
{
  typedef CGAL::AABB_tree<AABBTraits>                                         AABB_tree;
  typedef typename AABB_tree::AABB_traits::FT                                 FT;

public:
  Bounded_distance_placement(const FT dist,
                             const AABB_tree& tree,
                             const BasePlacement& placement = BasePlacement())
    :
      m_sq_threshold_dist(CGAL::square(dist)),
      m_tree_ptr(&tree),
      m_base_placement(placement)
  { }

  template <typename Profile>
  boost::optional<typename Profile::Point>
  operator()(const Profile& profile) const
  {
    typedef typename Profile::Geom_traits                                     Geom_traits;
    typedef typename Profile::Point                                           Point;

    boost::optional<typename Profile::Point> op = m_base_placement(profile);
    if(op)
    {
      CGAL_assertion(m_tree_ptr != nullptr);
      CGAL_assertion(!m_tree_ptr->empty());

      const Point& p = *op;

      const Point& cp = m_tree_ptr->best_hint(p).first;

      if(CGAL::compare_squared_distance(p, cp, m_sq_threshold_dist) != LARGER ||
         m_tree_ptr->do_intersect(CGAL::Sphere_3<Geom_traits>(p, m_sq_threshold_dist)))
        return op;

      return boost::optional<Point>();
    }

    return op;
  }

private:
  const FT m_sq_threshold_dist;
  mutable const AABB_tree* m_tree_ptr;

  const BasePlacement m_base_placement;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H
