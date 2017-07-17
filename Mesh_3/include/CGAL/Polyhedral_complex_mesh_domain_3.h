// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
// Copyright (c) 2014-2017 GeometryFactory Sarl (France)
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_POLYHEDRAL_COMPLEX_MESH_DOMAIN_3_H
#define CGAL_POLYHEDRAL_COMPLEX_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/Random.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Mesh_3/Polyline_with_context.h>
#include <CGAL/Polygon_mesh_processing/Detect_features_in_polyhedra.h>
#include <CGAL/Mesh_3/properties_Polyhedron_3.h>

#include <CGAL/enum.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/dynamic_bitset.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/foreach.hpp>

#include <string>
#include <vector>
#include <fstream>


namespace CGAL {
/// @cond DEVELOPERS
namespace internal {
namespace Mesh_3 {

template <typename Graph>
void dump_graph_edges(std::ostream& out, const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

  out.precision(17);
  BOOST_FOREACH(edge_descriptor e, edges(g))
  {
    vertex_descriptor s = source(e, g);
    vertex_descriptor t = target(e, g);
    out << "2 " << g[s] << " " << g[t] << "\n";
  }
}

template <typename Graph>
void dump_graph_edges(const char* filename, const Graph& g)
{
  std::ofstream out(filename);
  dump_graph_edges(out, g);
}

template <typename Kernel>
struct Angle_tester
{
  template <typename vertex_descriptor, typename Graph>
  bool operator()(vertex_descriptor& v, const Graph& g) const
  {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
    if (out_degree(v, g) != 2)
      return true;
    else
    {
      out_edge_iterator out_edge_it, out_edges_end;
      boost::tie(out_edge_it, out_edges_end) = out_edges(v, g);

      vertex_descriptor v1 = target(*out_edge_it++, g);
      vertex_descriptor v2 = target(*out_edge_it++, g);
      CGAL_assertion(out_edge_it == out_edges_end);

      const typename Kernel::Point_3& p = g[v];
      const typename Kernel::Point_3& p1 = g[v1];
      const typename Kernel::Point_3& p2 = g[v2];

      return (CGAL::angle(p1, p, p2) == CGAL::ACUTE);
    }
  }
};

template <typename Polyhedron>
struct Is_featured_edge {
  const Polyhedron* polyhedron;
  Is_featured_edge() : polyhedron(0) {} // required by boost::filtered_graph
  Is_featured_edge(const Polyhedron& polyhedron) : polyhedron(&polyhedron) {}

  bool operator()(typename boost::graph_traits<Polyhedron>::edge_descriptor e) const {
    return halfedge(e, *polyhedron)->is_feature_edge();
  }
}; // end Is_featured_edge<Polyhedron>

template <typename Polyhedron>
struct Is_border_edge {
  const Polyhedron* polyhedron;
  Is_border_edge() : polyhedron(0) {} // required by boost::filtered_graph
  Is_border_edge(const Polyhedron& polyhedron) : polyhedron(&polyhedron) {}

  bool operator()(typename boost::graph_traits<Polyhedron>::edge_descriptor e) const {
    return is_border(halfedge(e, *polyhedron), *polyhedron) ||
      is_border(opposite(halfedge(e, *polyhedron), *polyhedron), *polyhedron);
  }
}; // end Is_featured_edge<Polyhedron>

template<typename Polyhedral_mesh_domain,
         typename Polyline_with_context,
         typename Graph>
struct Extract_polyline_with_context_visitor
{
  typedef typename Polyhedral_mesh_domain::Polyhedron_type Polyhedron;
  std::vector<Polyline_with_context>& polylines;
  const Graph& graph;

  Extract_polyline_with_context_visitor
  (const Graph& graph,
   typename std::vector<Polyline_with_context>& polylines)
    : polylines(polylines), graph(graph)
  {}

  void start_new_polyline()
  {
    polylines.push_back(Polyline_with_context());
  }

  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor vd)
  {
    if(polylines.back().polyline_content.empty()) {
      polylines.back().polyline_content.push_back(graph[vd]);
    }
  }

  void add_edge(typename boost::graph_traits<Graph>::edge_descriptor ed)
  {
    typename boost::graph_traits<Graph>::vertex_descriptor
      s = source(ed, graph),
      t = target(ed, graph);
    Polyline_with_context& polyline = polylines.back();
    CGAL_assertion(!polyline.polyline_content.empty());
    if(polyline.polyline_content.back() != graph[s]) {
      polyline.polyline_content.push_back(graph[s]);
    } else if(polyline.polyline_content.back() != graph[t]) {
      // if the edge is zero-length, it is ignored
      polyline.polyline_content.push_back(graph[t]);
    }
    const typename boost::edge_bundle_type<Graph>::type &
      set_of_indices = graph[ed];
    polyline.context.adjacent_patches_ids.insert(set_of_indices.begin(),
                                                 set_of_indices.end());
  }

  void end_polyline()
  {
    // ignore degenerated polylines
    if(polylines.back().polyline_content.size() < 2)
      polylines.resize(polylines.size() - 1);
    // else {
    //   std::cerr << "Polyline with " << polylines.back().polyline_content.size()
    //             << " vertices, incident to "
    //             << polylines.back().context.adjacent_patches_ids.size()
    //             << " patches:\n ";
    //   for(auto p: polylines.back().polyline_content)
    //     std::cerr << " " << p;
    //   std::cerr << "\n";
    // }
  }
};


} // end CGAL::internal::Mesh_3
} // end CGAL::internal

/// @endcond

/*!
\ingroup PkgMesh_3Domains

The class `Polyhedral_complex_mesh_domain_3` implements a domain
defined by a collection of polyhedral surfaces, forming a complex.

The constraints on the complex are:
  - a polyhedral surface of the complex cannot self-intersect,
  - two polyhedral surfaces of the complex are either disjoint, or
share a subset of their border edges.

It is a model of the concept `MeshDomainWithFeatures_3`. It also provides a
member function to automatically detect sharp features and boundaries from
the input polyhedral surface(s).

The union of the polyhedral surfaces is a non-manifold surface, called the
2D surface of the domain. It is locally manifold, with or without borders,
but at the intersections of polyhedral surfaces.

The complement of the 2D surface is decomposed into:
 - zero, one, or many sub-domains of the mesh domain,
 - plus the exterior of the mesh domain.

If the domain has sub-domains, each one must be the union of one or many
connected components of the complement of the 2D surface. The sub-domains
have indices, of integral type `Subdomain_index`, and the exterior of the
mesh domain is associated with the subdomain index `0`, like for any mesh
domain in CGAL.

Each polyhedral surface is oriented, and has two sides. The positive side
is union of the positive side of all of its facets, usually named the
"exterior" of the surface. The negative side is the other side. The use of
`Polyhedral_complex_mesh_domain_3` assumes that for each polyhedral
surface, the sub-domain indices on both sides are known.

\tparam Polyhedron stands for the type of the input polyhedral surface(s).
The only requirements for this type is that the triangles of the surfaces
must be accessible through an object of the class
`TriangleAccessor`. @todo Document the requirements.

\tparam IGT_ stands for a geometric traits class
providing the types and functors required to implement
the intersection tests and intersection computations
for polyhedral boundary surfaces. This parameter has to be instantiated
with a model of the concept `IntersectionGeometricTraits_3`.

\tparam TriangleAccessor provides access to the triangles
of the input polyhedral
surface. It must be a model of the concept
`TriangleAccessor_3`. It defaults to
`Triangle_accessor_3<Polyhedron,IGT_>`. The type `IGT_::Triangle_3` must
be identical to the type `TriangleAccessor::Triangle_3`.

\cgalModels `MeshDomainWithFeatures_3`

\sa `TriangleAccessor_3`
\sa `IntersectionGeometricTraits_3`
\sa `CGAL::Triangle_accessor_3<Polyhedron_3<K>,K>`
\sa `CGAL::make_mesh_3()`.
\sa `CGAL::Mesh_domain_with_polyline_features_3<MeshDomain>`
\sa `CGAL::Polyhedral_mesh_domain_3<Polyhedron,IGT_,TriangleAccessor>`
\sa `CGAL::Mesh_polyhedron_3<IGT_>`
*/
template < class IGT_,
           class Polyhedron = typename Mesh_polyhedron_3<IGT_>::type,
           class TriangleAccessor=Triangle_accessor_3<Polyhedron,IGT_>
           >
class Polyhedral_complex_mesh_domain_3
  : public Mesh_domain_with_polyline_features_3<
      Polyhedral_mesh_domain_3< Polyhedron,
                                IGT_,
                                TriangleAccessor,
                                Tag_true,   //Use_patch_id_tag
                                Tag_true > >//Use_exact_intersection_tag
{
  /// @cond DEVELOPERS
protected:
  typedef boost::adjacency_list<
    boost::setS, // this avoids parallel edges
    boost::vecS,
    boost::undirectedS,
    typename Polyhedron::Point,
    typename Polyhedron::Vertex::Set_of_indices> Featured_edges_copy_graph;

private:
  typedef IGT_ IGT;
  typedef Polyhedral_mesh_domain_3<Polyhedron, IGT_, TriangleAccessor,
                                   Tag_true, Tag_true >       BaseBase;
  typedef Polyhedral_complex_mesh_domain_3<IGT_, Polyhedron,
                                           TriangleAccessor>  Self;
  /// @endcond

public:
  /// The base class
  typedef Mesh_domain_with_polyline_features_3<
    Polyhedral_mesh_domain_3<
      Polyhedron, IGT_, TriangleAccessor,
      Tag_true, Tag_true > > Base;
  /*!
  Numerical type.
  */
  typedef typename Base::FT        FT;

  /// The polyhedron type
  typedef Polyhedron Polyhedron_type;
  
  /// \name Index types
  /// @{
  /// The types are `int` or types compatible with `int`.
  typedef typename Base::Corner_index         Corner_index;
  typedef typename Base::Curve_segment_index  Curve_segment_index;
  typedef typename Base::Surface_patch_index  Surface_patch_index;
  typedef typename Base::Subdomain_index      Subdomain_index;
  /// @}

  /// @cond DEVELOPERS
  typedef typename Base::Ray_3                Ray_3;
  typedef typename Base::Index                Index;
  typedef Surface_patch_index Patch_id;

  typedef typename Base::Subdomain            Subdomain;
  typedef typename Base::Bounding_box         Bounding_box;
  typedef typename Base::AABB_tree            AABB_tree;
  typedef typename Base::AABB_primitive       AABB_primitive;
  typedef typename Base::AABB_primitive_id    AABB_primitive_id;
  // Backward compatibility
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index                 Surface_index;
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX

  typedef typename Base::R         R;
  typedef typename Base::Point_3   Point_3;

  typedef CGAL::Tag_true           Has_features;

  typedef std::vector<Point_3> Bare_polyline;
  typedef Mesh_3::Polyline_with_context<Surface_patch_index, Curve_segment_index,
                                        Bare_polyline > Polyline_with_context;
  /// @endcond

  /// Constructor
  /*! Constructs a domain defined by a set of polyhedral surfaces,
  describing a polyhedral complex.
  @param begin first iterator on the input polyhedral surfaces
  @param end past the end iterator on the input polyhedral surfaces
  @param indices_begin first iterator on the pairs of subdomain indices
             (two subdomain indices per input polyhedral surface),
             corresponding to the first input polyhedral surface
  @param indices_end past the end iterator on the pairs of subdomain indices

  @tparam InputPolyhedraIterator model of `InputIterator`, holding `Polyhedron`'s
  @tparam InputPairOfSubdomainIndicesIterator model of `InputIterator`, holding
              `std::pair<Subdomain_index, Subdomain_index>`

  @pre `std::distance(begin, end) == std::distance(indices_begin, indices_end)`
  */
  template <typename InputPolyhedraIterator,
            typename InputPairOfSubdomainIndicesIterator>
  Polyhedral_complex_mesh_domain_3
  ( InputPolyhedraIterator begin,
    InputPolyhedraIterator end,
    InputPairOfSubdomainIndicesIterator indices_begin,
    InputPairOfSubdomainIndicesIterator indices_end
#ifndef DOXYGEN_RUNNING
    , CGAL::Random* p_rng = NULL
#endif
    )
    : Base(p_rng)
    , patch_indices(indices_begin, indices_end)
    , borders_detected_(false)
  {
    stored_polyhedra.reserve(std::distance(begin, end));
    CGAL_assertion(stored_polyhedra.capacity() ==
                   std::size_t(std::distance(indices_begin, indices_end)));
    for (; begin != end; ++begin) {
      stored_polyhedra.push_back(*begin);
      this->add_primitives(stored_polyhedra.back());
    }
    this->build();
  }

  /// @cond DEVELOPERS
  Polyhedral_complex_mesh_domain_3
    (
    CGAL::Random* p_rng = NULL
    )
    : Base(p_rng)
    , borders_detected_(false)
  {}

  const std::vector<Polyhedron>& polyhedra() const {
    return stored_polyhedra;
  }
  /// @endcond

  /// @cond DEVELOPERS
  /*!
  * construct_initial_points_object() is one of the very first functions called
  * when make_mesh_3 starts
  * BEFORE make_mesh_3 starts, we have to make sure that (at least) the borders have
  * been detected, and the polyhedral complex internal data structures initialized
  * So, this function is overloaded to make sure they are (checking that
  * borders_detected_ is false is enough)
  * Then, call the base class function
  */
  typename BaseBase::Construct_initial_points construct_initial_points_object() const
  {
    Self* self = const_cast<Self*>(this);
    if (!borders_detected_)
      self->detect_borders(self->stored_polyhedra, true);
    return this->BaseBase::construct_initial_points_object();
  }

  void detect_features(FT angle_in_degree,
                       std::vector<Polyhedron_type>& p,
                       const bool dont_protect);//if true, features will not be protected
  /// @endcond
  /*!
  Detects sharp features and boundaries of the polyhedral components of the complex
  (including potential internal polyhedra),
  and inserts them as features of the domain. `angle_bound` gives the maximum
  angle (in degrees) between the two normal vectors of adjacent triangles.
  For an edge of the polyhedron, if the angle between the two normal vectors of its
  incident facets is bigger than the given bound, then the edge is considered as
  a feature edge, and inserted as a feature of the domain.
  */ 
  void detect_features(FT angle_bound = FT(60)) {
    detect_features(angle_bound, stored_polyhedra, false/*do protect*/);
  }

  /// @cond DEVELOPERS
  void detect_borders(std::vector<Polyhedron_type>& p, const bool dont_protect);
  /// @endcond
  /*!
  Detects border edges of the polyhedral components of the complex,
  and inserts them as features of the domain.
  This function should be called alone only, and not before or after `detect_features()`.
  */
  void detect_borders() {
    detect_borders(stored_polyhedra, false/*do protect*/);
  }

  /// @cond DEVELOPERS
  template <typename Surf_p_index>
  void reindex_patches(const std::vector<Surf_p_index>& map);

  template<typename PointSet>
  void merge_duplicated_points(const PointSet& duplicated_points);

  const std::pair<Subdomain_index, Subdomain_index>&
  incident_subdomains_indices(Surface_patch_index patch_id) const
  {
    CGAL_assertion(patch_id > 0 &&
                   std::size_t(patch_id) < patch_id_to_polyhedron_id.size());

    const std::size_t polyhedron_id = this->patch_id_to_polyhedron_id[patch_id];
    return this->patch_indices[polyhedron_id];
  }

  const std::vector<Surface_patch_index>& boundary_patches() const
  {
    return this->boundary_patches_ids;
  }

  const std::vector<std::size_t>& inside_polyhedra() const
  {
    return this->inside_polyhedra_ids;
  }

  const std::vector<std::size_t>& boundary_polyhedra() const
  {
    return this->boundary_polyhedra_ids;
  }

  void compute_boundary_patches()
  {
    if (!this->boundary_patches_ids.empty())
      return;
    //patches are numbered from 1 to n
    //patch_id_to_polyhedron_id has size n+1
    for (Surface_patch_index patch_id = 1;
      static_cast<std::size_t>(patch_id) < patch_id_to_polyhedron_id.size();
      ++patch_id)
    {
      const std::pair<Subdomain_index, Subdomain_index>& subdomains
        = incident_subdomains_indices(patch_id);
      if (subdomains.first == 0 || subdomains.second == 0)
        this->boundary_patches_ids.push_back(patch_id);
    }
    for(std::size_t poly_id = 0, end = stored_polyhedra.size();
        poly_id < end; ++poly_id)
    {
      const std::pair<Subdomain_index, Subdomain_index>& subdomains
        = this->patch_indices[poly_id];
      if (subdomains.first != 0 && subdomains.second != 0)
        this->inside_polyhedra_ids.push_back(poly_id);
      else
        this->boundary_polyhedra_ids.push_back(poly_id);
    }
  }
  /// @endcond

  /// @cond DEVELOPERS
  template <typename C3t3>
  void add_vertices_to_c3t3_on_patch_without_feature_edges(C3t3& c3t3) const {
#if CGAL_MESH_3_VERBOSE
    std::cout << "add_vertices_to_c3t3_on_patch_without_feature_edges...";
    std::cout.flush();
#endif
    CGAL::Random random(0);
    typedef typename C3t3::Triangulation Tr;
    Tr& tr = c3t3.triangulation();
    typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
    typename Tr::Geom_traits::Construct_weighted_point_3 cwp
      = tr.geom_traits().construct_weighted_point_3_object();

    const std::size_t nb_of_patch_plus_one = this->nb_of_patch_plus_one();
    const std::size_t nb_of_extra_vertices_per_patch = 20;
    std::vector<std::vector<Vertex_const_handle> >
      several_vertices_on_patch(nb_of_patch_plus_one);
    BOOST_FOREACH(const Polyhedron& p, this->stored_polyhedra)
    {
      for (typename Polyhedron::Vertex_const_iterator
        vit = p.vertices_begin(), end = p.vertices_end();
        vit != end; ++vit)
      {
        if (vit->is_feature_vertex()) { continue; }
        const Patch_id patch_id = vit->halfedge()->face()->patch_id();
        CGAL_assertion(static_cast<std::size_t>(patch_id) <= nb_of_patch_plus_one);
        typename Tr::Vertex_handle tr_v = tr.nearest_power_vertex(vit->point());
        if (tr_v != typename Tr::Vertex_handle()) {
          typedef typename IGT::Sphere_3 Sphere_3;
          const Sphere_3 sphere(tr_v->point().point(), tr_v->point().weight());
          if (!sphere.has_on_negative_side(vit->point())) continue;
        }
        if (several_vertices_on_patch[patch_id].size() <
          nb_of_extra_vertices_per_patch)
        {
          several_vertices_on_patch[patch_id].push_back(vit);
        }
        else {
          int i = random.uniform_smallint(0,
            static_cast<int>(nb_of_extra_vertices_per_patch - 1));
          several_vertices_on_patch[patch_id][i] = vit;
        }
      }
    }
    for (Patch_id patch_id = 1; std::size_t(patch_id) < nb_of_patch_plus_one;
      ++patch_id)
    {
      if (this->patch_has_featured_edges.test(patch_id)) {
        if (!several_vertices_on_patch[patch_id].empty()) {
          Vertex_const_handle v =
            several_vertices_on_patch[patch_id].front();
          typename Tr::Vertex_handle tr_v = tr.nearest_power_vertex(v->point());
          FT sq_dist_v = CGAL::squared_distance(v->point(),
            tr_v->point().point());
          for (std::size_t i = 1,
            size = several_vertices_on_patch[patch_id].size();
            i < size; ++i)
          {
            Vertex_const_handle other_v =
              several_vertices_on_patch[patch_id][i];
            tr_v = tr.nearest_power_vertex(other_v->point());
            const FT sq_dist_other_v =
              CGAL::squared_distance(other_v->point(),
              tr_v->point().point());
            if (sq_dist_other_v > sq_dist_v) {
              v = other_v;
              sq_dist_v = sq_dist_other_v;
            }
          }
          tr_v = tr.insert(cwp(v->point()));
          c3t3.set_dimension(tr_v, 2);
          c3t3.set_index(tr_v, patch_id);
        } // end if several_vertices_on_patch is empty for patch_id
      } // end if patch has featured edges
      else { // the patch is closed
        BOOST_FOREACH(Vertex_const_handle v,
          several_vertices_on_patch[patch_id])
        {
          typename Tr::Vertex_handle tv = tr.insert(cwp(v->point()));
          c3t3.set_dimension(tv, 2);
          c3t3.set_index(tv, patch_id);
        }
      }
    }
#if CGAL_MESH_3_VERBOSE
    std::cout << "\badd_vertices_to_c3t3_on_patch_without_feature_edges done.";
    std::cout << std::endl;
#endif
  }

  struct Is_in_domain
  {
    Is_in_domain(const Polyhedral_complex_mesh_domain_3& domain)
      : r_domain_(domain) {}

    boost::optional<AABB_primitive_id> shoot_a_ray_1(const Ray_3 ray) const {
      return r_domain_.bounding_aabb_tree_ptr()->
        first_intersected_primitive(ray);
    }

#if USE_ALL_INTERSECTIONS
    boost::optional<AABB_primitive_id> shoot_a_ray_2(const Ray_3 ray) const {
      const Point_3& p = ray.source();
      typedef typename AABB_tree::
        template Intersection_and_primitive_id<Ray_3>::Type Inter_and_prim;
      std::vector<Inter_and_prim> all_intersections;
      r_domain_.bounding_aabb_tree_ptr()->
        all_intersections(ray, std::back_inserter(all_intersections));
      if(all_intersections.empty())
        return boost::none;
      else {
        for(const Inter_and_prim& i_p: all_intersections) {
          if(boost::get<Point_3>( &i_p.first) == 0) return AABB_primitive_id();
        }
        auto it = std::min_element
          (all_intersections.begin(), all_intersections.end(),
           [p](const Inter_and_prim& a,
               const Inter_and_prim& b)
           {
             const Point_3& pa = boost::get<Point_3>(a.first);
             const Point_3& pb = boost::get<Point_3>(b.first);
             return compare_distance_to_point(p, pa, pb)
             == CGAL::SMALLER;
           });
        return it->second;
      }
    }
#endif // USE_ALL_INTERSECTIONS

    Subdomain operator()(const Point_3& p) const {
      if(r_domain_.bounding_aabb_tree_ptr() == 0) return Subdomain();
      const Bounding_box& bbox = r_domain_.bounding_aabb_tree_ptr()->bbox();

      if(   p.x() < bbox.xmin() || p.x() > bbox.xmax()
            || p.y() < bbox.ymin() || p.y() > bbox.ymax()
            || p.z() < bbox.zmin() || p.z() > bbox.zmax() )
      {
        return Subdomain();
      }
  
      // Shoot ray
      typename IGT::Construct_ray_3 ray = IGT().construct_ray_3_object();
      typename IGT::Construct_vector_3 vector = IGT().construct_vector_3_object();

      while(true) {
        Random_points_on_sphere_3<Point_3> random_point(1.);

        const Ray_3 ray_shot = ray(p, vector(CGAL::ORIGIN,*random_point));

#define USE_ALL_INTERSECTIONS 0
#if USE_ALL_INTERSECTIONS
        boost::optional<AABB_primitive_id> opt = shoot_a_ray_2(ray_shot);
#else // first_intersected_primitive
        boost::optional<AABB_primitive_id> opt = shoot_a_ray_1(ray_shot);
#endif // first_intersected_primitive
        // for(int i = 0; i < 20; ++i) {
        //   const Ray_3 ray_shot2 = ray(p, vector(CGAL::ORIGIN,*random_point));
        //   boost::optional<AABB_primitive_id> opt2 = shoot_a_ray_1(ray_shot2);
        //   if(opt != opt2) {
        //     if(!opt  && *opt2 == 0) continue;
        //     if(!opt2 && *opt  == 0) continue;
        //     std::cerr << "Not the same result for:\n  "
        //               << ray_shot
        //               << "\n  " << ray_shot2 << std::endl;
        //     abort();
        //   }
        // }
        if(!opt)
          return Subdomain();
        else {
          typename Polyhedron::Facet_const_handle fh = *opt;
          if(fh == AABB_primitive_id()) continue; // loop
          const std::pair<Subdomain_index, Subdomain_index>& pair =
            r_domain_.incident_subdomains_indices(fh->patch_id());
          typename Polyhedron::Halfedge_const_handle he = fh->halfedge();
          const Point_3& a = he->vertex()->point();
          const Point_3& b = he->next()->vertex()->point();
          const Point_3& c = he->next()->next()->vertex()->point();
          switch(orientation(a, b, c, p)) {
          case NEGATIVE: // inner region
            return pair.first == 0 ?
              Subdomain() :
              Subdomain(Subdomain_index(pair.first));
          case POSITIVE: // outer region
            return pair.second == 0 ?
              Subdomain() :
              Subdomain(Subdomain_index(pair.second));
          default: /* COPLANAR */ continue; // loop
          } // end switch on the orientation
        } // opt
      } // end while(true)
    } // end operator()
  private:
    const Polyhedral_complex_mesh_domain_3& r_domain_;
  };
  Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }
  /// @endcond

protected:
  /// @cond DEVELOPERS
  void initialize_ts(Polyhedron_type& p) const;

  void add_features_from_split_graph_into_polylines(Featured_edges_copy_graph& graph);

  template <typename Edge_predicate, typename P2Vmap>
  void add_featured_edges_to_graph(const Polyhedron& poly,
                                   const Edge_predicate& pred,
                                   Featured_edges_copy_graph& graph,
                                   P2Vmap& p2vmap);

  std::size_t nb_of_patch_plus_one() const
  {
    return patch_id_to_polyhedron_id.size();
  }
  /// @endcond

protected:
  /// @cond DEVELOPERS
  std::vector<Polyhedron> stored_polyhedra;
  std::vector<std::pair<Subdomain_index, Subdomain_index> > patch_indices;
  std::vector<std::size_t> patch_id_to_polyhedron_id;
  boost::dynamic_bitset<> patch_has_featured_edges;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  std::vector<std::vector<Vertex_handle> > several_vertices_on_patch;
  std::vector<Surface_patch_index> boundary_patches_ids;
  std::vector<std::size_t> inside_polyhedra_ids;
  std::vector<std::size_t> boundary_polyhedra_ids;
  mutable bool borders_detected_;
  /// @endcond

private:
  // Disabled copy constructor & assignment operator
  Polyhedral_complex_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Polyhedral_complex_mesh_domain_3


///@cond DEVELOPERS
template < typename GT_, typename P_, typename TA_>
void
Polyhedral_complex_mesh_domain_3<GT_,P_,TA_>::
initialize_ts(Polyhedron_type& p) const
{
  std::size_t ts = 0;
  for(typename Polyhedron_type::Vertex_iterator v = p.vertices_begin(),
      end = p.vertices_end() ; v != end ; ++v)
  {
    v->set_time_stamp(ts++);
  }
  for(typename Polyhedron_type::Facet_iterator fit = p.facets_begin(),
       end = p.facets_end() ; fit != end ; ++fit )
  {
    fit->set_time_stamp(ts++);
  }
  for(typename Polyhedron_type::Halfedge_iterator hit = p.halfedges_begin(),
       end = p.halfedges_end() ; hit != end ; ++hit )
  {
    hit->set_time_stamp(ts++);
  }
}

template < typename GT_, typename P_, typename TA_>
void
Polyhedral_complex_mesh_domain_3<GT_, P_, TA_>::
detect_borders(std::vector<Polyhedron_type>& poly, const bool dont_protect)
{
  if (borders_detected_)
    return;//border detection has already been done

  detect_features(180, poly, dont_protect);

  borders_detected_ = true;
}

template < typename GT_, typename P_, typename TA_>
void
Polyhedral_complex_mesh_domain_3<GT_,P_,TA_>::
detect_features(FT angle_in_degree,
                std::vector<Polyhedron_type>& poly,
                const bool dont_protect)
{
  CGAL_assertion(!borders_detected_);
  if (borders_detected_)
    return;//prevent from not-terminating

  typedef Featured_edges_copy_graph G_copy;
  G_copy g_copy;
  typedef typename boost::graph_traits<G_copy>::vertex_descriptor vertex_descriptor;
  typedef std::map<typename Polyhedron_type::Point,
                   vertex_descriptor> P2vmap;
  // TODO: replace this map by and unordered_map
  P2vmap p2vmap;

  namespace PMP = CGAL::Polygon_mesh_processing;
  PMP::Detect_features_in_polyhedra<Polyhedron_type, Surface_patch_index> detect_features;
  BOOST_FOREACH(Polyhedron_type& p, poly)
  {
    initialize_ts(p);

#if CGAL_MESH_3_VERBOSE
    std::size_t poly_id = &p-&poly[0];
    std::cerr << "Polyhedron #" << poly_id << " :\n";
    std::cerr << "  material #" << patch_indices[poly_id].first << "\n";
    std::cerr << "  material #" << patch_indices[poly_id].second << "\n";
#endif // CGAL_MESH_3_VERBOSE

    // Get sharp features
    detect_features.detect_sharp_edges(p, angle_in_degree);
    detect_features.detect_surface_patches(p);
    detect_features.detect_vertices_incident_patches(p);

    internal::Mesh_3::Is_featured_edge<Polyhedron_type> is_featured_edge(p);

    add_featured_edges_to_graph(p, is_featured_edge, g_copy, p2vmap);
  }
  const std::size_t nb_of_patch_plus_one =
    detect_features.maximal_surface_patch_index()+1;
  this->patch_id_to_polyhedron_id.resize(nb_of_patch_plus_one);
  this->patch_has_featured_edges.resize(nb_of_patch_plus_one);
  this->several_vertices_on_patch.resize(nb_of_patch_plus_one);
#if CGAL_MESH_3_VERBOSE
  std::cerr << "Number of patches: " << (nb_of_patch_plus_one - 1) << std::endl;
#endif
  BOOST_FOREACH(Polyhedron_type& p, poly)
  {
    const std::size_t polyhedron_id = &p - &poly[0];
    BOOST_FOREACH(typename Polyhedron_type::Facet_const_handle fh, faces(p))
    {
      patch_id_to_polyhedron_id[fh->patch_id()] = polyhedron_id;
    }
    for(typename Polyhedron_type::Halfedge_iterator
          heit = p.halfedges_begin(), end  = p.halfedges_end();
        heit != end; ++heit)
    {
      if(is_border(heit, p) || !heit->is_feature_edge()) continue;
      patch_has_featured_edges.set(heit->face()->patch_id());
    }
    for(typename Polyhedron_type::Vertex_iterator
          vit = p.vertices_begin(), end  = p.vertices_end();
        vit != end; ++vit)
    {
      if( vit->is_feature_vertex() ) { continue; }
      const Patch_id patch_id = vit->halfedge()->face()->patch_id();
      if(patch_has_featured_edges.test(patch_id)) continue;
      several_vertices_on_patch[patch_id].push_back(vit);
    }
  }
  for(Patch_id patch_id = 1; std::size_t(patch_id) < nb_of_patch_plus_one;
      ++patch_id)
  {
    CGAL_assertion(patch_has_featured_edges.test(patch_id) ==
                   several_vertices_on_patch[patch_id].empty() );
    if(several_vertices_on_patch[patch_id].empty()) continue;
    std::random_shuffle(several_vertices_on_patch[patch_id].begin(),
                        several_vertices_on_patch[patch_id].end());
    if(several_vertices_on_patch[patch_id].size()>20)
      several_vertices_on_patch[patch_id].resize(20);
#if __cplusplus > 201103L
    several_vertices_on_patch.shrink_to_fit();
#endif
  }
  if (!dont_protect)
    add_features_from_split_graph_into_polylines(g_copy);

  borders_detected_ = true;/*done by Mesh_3::detect_features*/

  this->compute_boundary_patches();
}

template < typename GT_, typename P_, typename TA_>
template < typename PointSet >
void
Polyhedral_complex_mesh_domain_3<GT_,P_,TA_>::
merge_duplicated_points(const PointSet& duplicated_points)
{
  typedef typename Polyhedron_type::Point Point;
  typedef typename Polyhedron_type::Halfedge_around_vertex_const_circulator HVcirc;
  typedef std::pair<Point, const Polyhedron_type*> Point_and_mesh;
  typedef std::multimap<Point_and_mesh, Patch_id> Patch_multimap;
  typedef typename Patch_multimap::iterator                   Patch_iterator;
  typedef typename Patch_multimap::value_type                 Pt_patch_pair;

  if (duplicated_points.empty())
    return;

  Patch_multimap patches;
  BOOST_FOREACH(const Polyhedron_type& p, stored_polyhedra)
  {
    for (typename Polyhedron_type::Vertex_const_iterator
      vit = p.vertices_begin(), end = p.vertices_end();
      vit != end; ++vit)
    {
      if (duplicated_points.find(vit->point()) == duplicated_points.end())
        continue;

      HVcirc it = vit->vertex_begin();
      HVcirc itend = it;
      do {
        if(!it->is_border()) {
          CGAL_assertion(static_cast<std::size_t>(it->face()->patch_id())
                         < this->nb_of_patch_plus_one());
          patches.insert(Pt_patch_pair(Point_and_mesh(vit->point(), &p),
                                       it->face()->patch_id()));
        }
        ++it;
      } while (it != itend);
    }
  }

  //collect ids to be modified 
  std::vector<Patch_id> new_ids(this->nb_of_patch_plus_one());
  for (std::size_t i = 0; i < this->nb_of_patch_plus_one(); ++i)
    new_ids[i] = Patch_id(i);

  for (Patch_iterator pit = patches.begin();
       pit != patches.end();
       /*pit = range_end inside loop*/)
  {
    const Patch_iterator range_begin = pit;

    // First loop: find the minimal patch index around `p` in a given
    // polyhedron, and the past-the-end iterator of the equal-range (the
    // loop will end with the first iterator when the point differs from
    // the point of the range).
    Patch_iterator it = range_begin;
    CGAL_assertion(static_cast<std::size_t>(it->second) < new_ids.size());
    Patch_id min_id_around_p = new_ids[it->second];
    ++it;
    for (; it != patches.end() && it->first == range_begin->first; ++it)
    {
      CGAL_assertion(static_cast<std::size_t>(it->second) < new_ids.size());
      min_id_around_p = (std::min)(min_id_around_p, new_ids[it->second]);
    }
    const Patch_iterator range_end = it;

    // In a second loop on the equal-range, update new_ids around p
    for (it = range_begin; it != range_end; ++it)
    {
      CGAL_assertion(static_cast<std::size_t>(it->second) < new_ids.size());
      new_ids[it->second] = min_id_around_p;
    }
    pit = range_end;
  }
  reindex_patches(new_ids);
}

template < typename GT_, typename P_, typename TA_>
void
Polyhedral_complex_mesh_domain_3<GT_,P_,TA_>::
add_features_from_split_graph_into_polylines(Featured_edges_copy_graph& g_copy)
{
  std::vector<Polyline_with_context> polylines;

  internal::Mesh_3::Extract_polyline_with_context_visitor<
    Polyhedral_complex_mesh_domain_3,
    Polyline_with_context,
    Featured_edges_copy_graph
    > visitor(g_copy, polylines);
  internal::Mesh_3::Angle_tester<GT_> angle_tester;
  split_graph_into_polylines(g_copy, visitor, angle_tester);

  this->add_features_with_context(polylines.begin(),
                                  polylines.end());

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  {//DEBUG
    std::ofstream og("polylines_graph.polylines.txt");
    og.precision(17);
    BOOST_FOREACH(const Polyline_with_context& poly, polylines)
    {
      og << poly.polyline_content.size() << " ";
      BOOST_FOREACH(const Point_3& p, poly.polyline_content)
        og << p << " ";
      og << std::endl;
    }
    og.close();
  }
#endif // CGAL_MESH_3_PROTECTION_DEBUG & 2

}

template < typename GT_, typename P_, typename TA_>
template <typename Edge_predicate, typename P2vmap>
void
Polyhedral_complex_mesh_domain_3<GT_,P_,TA_>::
add_featured_edges_to_graph(const Polyhedron_type& p,
                            const Edge_predicate& pred,
                            Featured_edges_copy_graph& g_copy,
                            P2vmap& p2vmap)
{
  typedef boost::filtered_graph<Polyhedron_type,
                                Edge_predicate > Featured_edges_graph;
  Featured_edges_graph orig_graph(p, pred);

  typedef Featured_edges_graph Graph;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Graph_vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor Graph_edge_descriptor;
  typedef Featured_edges_copy_graph G_copy;
  typedef typename boost::graph_traits<G_copy>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<G_copy>::edge_descriptor edge_descriptor;

  const Featured_edges_graph& graph = orig_graph;

  BOOST_FOREACH(Graph_vertex_descriptor v, vertices(graph)){
    vertex_descriptor vc;
    typename P2vmap::iterator it = p2vmap.find(v->point());
    if(it == p2vmap.end()) {
      vc = add_vertex(g_copy);
      g_copy[vc] = v->point();
      p2vmap[v->point()] = vc;
    }
  }

  BOOST_FOREACH(Graph_edge_descriptor e, edges(graph)){
    vertex_descriptor vs = p2vmap[source(e,graph)->point()];
    vertex_descriptor vt = p2vmap[target(e,graph)->point()];
    CGAL_warning_msg(vs != vt, "ignore self loop");
    if(vs != vt) {
      const std::pair<edge_descriptor, bool> pair = add_edge(vs,vt,g_copy);
      typename Polyhedron_type::Halfedge_handle he = halfedge(e, p);
      if(!is_border(he, p)) {
        g_copy[pair.first].insert(he->face()->patch_id());;
      }
      he = he->opposite();
      if(!is_border(he, p)) {
        g_copy[pair.first].insert(he->face()->patch_id());;
      }
    }
  }

#if CGAL_MESH_3_PROTECTION_DEBUG > 1
  {// DEBUG
    internal::Mesh_3::dump_graph_edges("edges-graph.polylines.txt", g_copy);
  }
#endif
}

template < typename GT_, typename P_, typename TA_>
template < typename Surf_p_index>
void
Polyhedral_complex_mesh_domain_3<GT_,P_,TA_>::
reindex_patches(const std::vector<Surf_p_index>& map)
{
  for(std::size_t i = 0, end = stored_polyhedra.size(); i < end; ++i) {
    Polyhedron_type& poly = stored_polyhedra[i];
    for(typename Polyhedron_type::Facet_iterator fit = poly.facets_begin(),
          end = poly.facets_end(); fit != end; ++fit)
    {
      fit->set_patch_id(map[fit->patch_id()]);
    }
  }
  BOOST_FOREACH(Surface_patch_index& id,
                boundary_patches_ids)
  {
    id = map[id];
  }
  Base::reindex_patches(map);
}
/// @endcond

} //namespace CGAL


#endif // CGAL_POLYHEDRAL_COMPLEX_MESH_DOMAIN_3_H
