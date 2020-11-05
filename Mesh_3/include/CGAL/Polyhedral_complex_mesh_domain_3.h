// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
// Copyright (c) 2014-2017 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_POLYHEDRAL_COMPLEX_MESH_DOMAIN_3_H
#define CGAL_POLYHEDRAL_COMPLEX_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/config.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Mesh_3/Polyline_with_context.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/internal/Mesh_3/helpers.h>

#include <CGAL/enum.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Random.h>
#include <CGAL/Union_find.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <string>
#include <vector>
#include <fstream>

namespace CGAL {

/*!
\ingroup PkgMesh3Domains

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

\tparam Polyhedron stands for the type of the input polyhedral surface(s),
model of `FaceListGraph`.

\tparam IGT_ stands for a geometric traits class
providing the types and functors required to implement
the intersection tests and intersection computations
for polyhedral boundary surfaces. This parameter has to be instantiated
with a model of the concept `IntersectionGeometricTraits_3`.

\cgalModels `MeshDomainWithFeatures_3`

\sa `IntersectionGeometricTraits_3`
\sa `CGAL::make_mesh_3()`.
\sa `CGAL::Mesh_domain_with_polyline_features_3<MeshDomain>`
\sa `CGAL::Polyhedral_mesh_domain_3<Polyhedron,IGT_,TriangleAccessor>`
\sa `CGAL::Mesh_polyhedron_3<IGT_>`
*/
template < class IGT_,
           class Polyhedron = typename Mesh_polyhedron_3<IGT_>::type,
           class TriangleAccessor=CGAL::Default
           >
class Polyhedral_complex_mesh_domain_3
  : public Mesh_domain_with_polyline_features_3<
      Polyhedral_mesh_domain_3< Polyhedron,
                                IGT_,
                                CGAL::Default,
                                int,   //Use_patch_id_tag
                                Tag_true > >//Use_exact_intersection_tag
{
public:
  /// The base class
  typedef Mesh_domain_with_polyline_features_3<
    Polyhedral_mesh_domain_3<
      Polyhedron, IGT_, CGAL::Default,
      int, Tag_true > > Base;
  /// @cond DEVELOPERS
private:
  typedef IGT_ IGT;
  typedef Polyhedral_mesh_domain_3<Polyhedron, IGT_, CGAL::Default,
                                   int, Tag_true >       BaseBase;
  typedef Polyhedral_complex_mesh_domain_3<IGT_, Polyhedron>  Self;

protected:
  typedef typename Base::Surface_patch_index Patch_id;
  typedef typename boost::property_map<Polyhedron,
                                       CGAL::vertex_incident_patches_t<Patch_id>
                                       >::type VIPMap;
  typedef typename boost::property_traits<VIPMap>::value_type Set_of_indices;

  typedef boost::adjacency_list<
    boost::setS, // this avoids parallel edges
    boost::vecS,
    boost::undirectedS,
    typename IGT::Point_3,
    Set_of_indices> Featured_edges_copy_graph;

  /// @endcond

public:
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
  typedef typename Base::Curve_index          Curve_index;
  typedef typename Base::Surface_patch_index  Surface_patch_index;
  typedef typename Base::Subdomain_index      Subdomain_index;
  /// @}

  /// @cond DEVELOPERS
  typedef typename Base::Ray_3                Ray_3;
  typedef typename Base::Index                Index;

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
  typedef Mesh_3::Polyline_with_context<Surface_patch_index, Curve_index,
                                        Bare_polyline > Polyline_with_context;
  /// @endcond

  /// Constructor
  /*! Constructs a domain defined by a set of polyhedral surfaces,
  describing a polyhedral complex.
  @param begin first iterator on the input polyhedral surfaces
  @param end past the end iterator on the input polyhedral surfaces
  @param indices_begin first iterator on the pairs of subdomain indices
             (two subdomain indices per input polyhedral surface),
             corresponding to the first input polyhedral surface.
             Each pair should be ordered as follows : the first (resp. second)
             element is the index of the subdomain that lies on the
             positively oriented (resp. negatively oriented) side of the
             corresponding input polyhedral surface.
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
    , CGAL::Random* p_rng = nullptr
#endif
    )
    : Base(p_rng)
    , patch_indices(indices_begin, indices_end)
    , borders_detected_(false)
  {
    patch_id_to_polyhedron_id.resize(std::distance(begin, end)+1);
    stored_polyhedra.reserve(patch_id_to_polyhedron_id.size()-1);
    CGAL_assertion(stored_polyhedra.capacity() ==
                   std::size_t(std::distance(indices_begin, indices_end)));

    Surface_patch_index sp_index = 1;

    for (; begin != end; ++begin) {
      typedef boost::graph_traits<Polyhedron_type> Graph_traits;
      typedef typename Graph_traits::face_descriptor face_descriptor;
      stored_polyhedra.push_back(*begin);
      patch_id_to_polyhedron_id[sp_index] = sp_index - 1;

      typedef typename boost::property_map<Polyhedron_type,
                                           CGAL::face_patch_id_t<Patch_id>
                                           >::type PIDMap;
      PIDMap pid_map = get(face_patch_id_t<Patch_id>(), stored_polyhedra.back());
      for(face_descriptor fd :
                    faces(stored_polyhedra.back())) {
        put(pid_map, fd, sp_index);
      }
      ++sp_index;
      this->add_primitives(stored_polyhedra.back());
    }
    this->build();
  }

  /// @cond DEVELOPERS
  Polyhedral_complex_mesh_domain_3
    (
    CGAL::Random* p_rng = nullptr
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
#ifdef CGAL_MESH_3_VERBOSE
    std::cout << "add_vertices_to_c3t3_on_patch_without_feature_edges...";
    std::cout.flush();
#endif
    CGAL::Random random(0);

    typedef typename C3t3::Triangulation                  Tr;
    typedef typename Tr::Weighted_point                   Weighted_point;
    typedef typename IGT::Sphere_3                        Sphere_3;
    typedef typename Polyhedron::Vertex_const_handle      Vertex_const_handle;

    Tr& tr = c3t3.triangulation();

    typename Tr::Geom_traits::Compute_weight_3 cw =
      tr.geom_traits().compute_weight_3_object();
    typename Tr::Geom_traits::Construct_point_3 cp =
      tr.geom_traits().construct_point_3_object();
    typename Tr::Geom_traits::Construct_weighted_point_3 cwp =
      tr.geom_traits().construct_weighted_point_3_object();

    const std::size_t nb_of_patch_plus_one = this->nb_of_patch_plus_one();
    const std::size_t nb_of_extra_vertices_per_patch = 20;
    std::vector<std::vector<Vertex_const_handle> >
      several_vertices_on_patch(nb_of_patch_plus_one);

    // For a patch, a "free" vertex is a vertex that is not on a feature,
    // and with a point that is not inside the set of protecting balls.
    std::vector<std::size_t> nb_of_free_vertices_on_patch(nb_of_patch_plus_one);

    // First path to count of the number of free vertices per patch...
    for(const Polyhedron& p : this->stored_polyhedra)
    {
      for (typename Polyhedron::Vertex_const_iterator
        vit = p.vertices_begin(), end = p.vertices_end();
        vit != end; ++vit)
      {
        if (vit->is_feature_vertex()) { continue; }
        const Patch_id patch_id = vit->halfedge()->face()->patch_id();
        CGAL_assertion(std::size_t(patch_id) <= nb_of_patch_plus_one);
        typename Tr::Vertex_handle tr_v = tr.nearest_power_vertex(vit->point());
        if (tr_v != typename Tr::Vertex_handle())
        {
          const Weighted_point& trv_wp = tr.point(tr_v);
          const Sphere_3 sphere(cp(trv_wp), cw(trv_wp));
          if (!sphere.has_on_unbounded_side(vit->point()))
            continue;
        }
        ++nb_of_free_vertices_on_patch[patch_id];
      }
    }
    std::vector<std::size_t>
      needed_vertices_on_patch = nb_of_free_vertices_on_patch;
    for(std::size_t i = 0, end = needed_vertices_on_patch.size(); i < end; ++i)
    {
      needed_vertices_on_patch[i] = (std::min)(nb_of_extra_vertices_per_patch,
                                               needed_vertices_on_patch[i]);
    }

    // Then a second path to fill `several_vertices_on_patch`...
    // The algorithm is adapted from SGI `random_sample_n`
    for(const Polyhedron& p : this->stored_polyhedra)
    {
      for (typename Polyhedron::Vertex_const_iterator
        vit = p.vertices_begin(), end = p.vertices_end();
        vit != end; ++vit)
      {
        if (vit->is_feature_vertex()) { continue; }
        const Patch_id patch_id = vit->halfedge()->face()->patch_id();
        CGAL_assertion(std::size_t(patch_id) <= nb_of_patch_plus_one);

        // If needed_vertices_on_patch is null, no need to proceed with the
        // rest of the loop.
        if(needed_vertices_on_patch[patch_id] == 0) continue;

        typename Tr::Vertex_handle tr_v = tr.nearest_power_vertex(vit->point());
        if (tr_v != typename Tr::Vertex_handle())
        {
          const Weighted_point& trv_wp = tr.point(tr_v);
          const Sphere_3 sphere(cp(trv_wp), cw(trv_wp));
          if (!sphere.has_on_unbounded_side(vit->point()))
            continue;
        }

        // here we have a new free vertex on patch #`patch_id`

        if(random.uniform_smallint(
               boost::uint32_t(0),
               boost::uint32_t(nb_of_free_vertices_on_patch[patch_id]))
           < boost::uint32_t(needed_vertices_on_patch[patch_id]))
        {
          several_vertices_on_patch[patch_id].push_back(vit);
          --needed_vertices_on_patch[patch_id];
        }
        --nb_of_free_vertices_on_patch[patch_id];
      }
    }
    for (Patch_id patch_id = 1; std::size_t(patch_id) < nb_of_patch_plus_one;
      ++patch_id)
    {
      CGAL_assertion(several_vertices_on_patch[patch_id].size()
                     <= nb_of_extra_vertices_per_patch);
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
        for(Vertex_const_handle v :
          several_vertices_on_patch[patch_id])
        {
          typename Tr::Vertex_handle tv = tr.insert(cwp(v->point()));
          c3t3.set_dimension(tv, 2);
          c3t3.set_index(tv, patch_id);
        }
      }
    }
#ifdef CGAL_MESH_3_VERBOSE
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
          if(AABB_primitive_id() == *opt) continue; // loop
          Surface_patch_index face_id = r_domain_.make_surface_index(*opt);
          const std::pair<Subdomain_index, Subdomain_index>& pair =
            r_domain_.incident_subdomains_indices(face_id);
          const typename IGT::Triangle_3 triangle =
            BaseBase::template Primitive_type<Polyhedron_type>::datum(*opt);
          const Point_3& a = triangle[0];
          const Point_3& b = triangle[1];
          const Point_3& c = triangle[2];
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
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
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
  typedef Polyhedron_type Polyhedron;
  typedef typename boost::property_map<Polyhedron,vertex_time_stamp_t>::type Vtmap;
  typedef typename boost::property_map<Polyhedron,halfedge_time_stamp_t>::type Htmap;
  typedef typename boost::property_map<Polyhedron,face_time_stamp_t>::type Ftmap;
  Vtmap vtm = get(vertex_time_stamp,p);
  Htmap htm = get(halfedge_time_stamp,p);
  Ftmap ftm = get(face_time_stamp,p);

  std::size_t ts = 0;
  typedef boost::graph_traits<Polyhedron> Graph_traits;
  for(typename Graph_traits::vertex_descriptor vd : vertices(p))
  {
    put(vtm,vd,ts++);
  }

  for(typename Graph_traits::face_descriptor fd : faces(p))
  {
    put(ftm,fd,ts++);
  }

  for(typename Graph_traits::halfedge_descriptor hd : halfedges(p))
  {
    put(htm,hd,ts++);
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
  typedef typename boost::graph_traits<G_copy>::vertex_descriptor
                                                      graph_vertex_descriptor;
  typedef std::map<typename IGT::Point_3,
                   graph_vertex_descriptor> P2vmap;
  // TODO: replace this map by and unordered_map
  P2vmap p2vmap;

  typedef typename boost::graph_traits<Polyhedron_type>::vertex_descriptor
                                                            vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron_type>::face_descriptor
                                                              face_descriptor;
  typedef typename boost::graph_traits<Polyhedron_type>::halfedge_descriptor
                                                          halfedge_descriptor;

  typedef typename boost::property_map<Polyhedron_type,
                                       CGAL::face_patch_id_t<Patch_id>
                                       >::type                        PIDMap;
  typedef typename boost::property_map<Polyhedron_type,
                                       CGAL::vertex_incident_patches_t<Patch_id>
                                       >::type                        VIPMap;
  typedef typename boost::property_map<Polyhedron_type,
                                       CGAL::edge_is_feature_t
                                       >::type                        EIFMap;
  typedef typename boost::property_map<Polyhedron_type,
                                       CGAL::vertex_feature_degree_t
                                       >::type                        VFDMap;
  namespace PMP = CGAL::Polygon_mesh_processing;
  std::size_t nb_of_patch_plus_one = 1;
  for(Polyhedron_type& p : poly)
  {
    initialize_ts(p);

#ifdef CGAL_MESH_3_VERBOSE
    std::size_t poly_id = &p-&poly[0];
    std::cerr << "Polyhedron #" << poly_id << " :\n";
    std::cerr << "  material #" << patch_indices[poly_id].first << "\n";
    std::cerr << "  material #" << patch_indices[poly_id].second << "\n";
#endif // CGAL_MESH_3_VERBOSE

    // Get sharp features
    PIDMap pid_map = get(face_patch_id_t<Patch_id>(), p);
    VIPMap vip_map = get(vertex_incident_patches_t<Patch_id>(), p);
    EIFMap eif = get(CGAL::edge_is_feature, p);
    VFDMap vertex_feature_degree_map = get(CGAL::vertex_feature_degree, p);
    nb_of_patch_plus_one +=PMP::sharp_edges_segmentation(p, angle_in_degree
      , eif
      , pid_map
      , PMP::parameters::first_index(nb_of_patch_plus_one)
                        .face_index_map(get_initialized_face_index_map(p))
                        .vertex_incident_patches_map(vip_map)
                        .vertex_feature_degree_map(vertex_feature_degree_map));

    Mesh_3::internal::Is_featured_edge<Polyhedron_type> is_featured_edge(p);

    add_featured_edges_to_graph(p, is_featured_edge, g_copy, p2vmap);
  }
  this->patch_id_to_polyhedron_id.resize(nb_of_patch_plus_one);
  this->patch_has_featured_edges.resize(nb_of_patch_plus_one);
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "Number of patches: " << (nb_of_patch_plus_one - 1) << std::endl;
#endif
  for(Polyhedron_type& p : poly)
  {
    PIDMap pid_map = get(face_patch_id_t<Patch_id>(), p);
    EIFMap eif = get(CGAL::edge_is_feature, p);
    const std::size_t polyhedron_id = &p - &poly[0];
    for(face_descriptor f : faces(p))
    {
      patch_id_to_polyhedron_id[get(pid_map, f)] = polyhedron_id;
    }
    for(halfedge_descriptor he : halfedges(p))
    {
      if(is_border(he, p) || !get(eif, edge(he, p))) continue;
      patch_has_featured_edges.set(get(pid_map, face(he, p)));
    }
    VFDMap vertex_feature_degree_map = get(CGAL::vertex_feature_degree, p);
    for(vertex_descriptor v : vertices(p))
    {
      if( get(vertex_feature_degree_map, v) != 0 ) { continue; }
      const Patch_id patch_id = get(pid_map, face(halfedge(v, p), p));
      if(patch_has_featured_edges.test(patch_id)) continue;
    }
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

  //collect ids to be modified
  typedef CGAL::Union_find<Patch_id> Union_find_t;
  Union_find_t union_find;
  std::vector<typename Union_find_t::handle> handles;
  handles.resize(this->nb_of_patch_plus_one());
  for (std::size_t i = 0; i < this->nb_of_patch_plus_one(); ++i) {
    handles[i] = union_find.make_set(Patch_id(i));
  }

  Patch_multimap patches;
  for(const Polyhedron_type& p : stored_polyhedra)
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
          CGAL_assertion(std::size_t(it->face()->patch_id()) <
                         this->nb_of_patch_plus_one());
          patches.insert(Pt_patch_pair(Point_and_mesh(vit->point(), &p),
                                       it->face()->patch_id()));
        }
        ++it;
      } while (it != itend);
    }
  }

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
    ++it;
    for (; it != patches.end() && it->first == range_begin->first; ++it)
    {
    }
    const Patch_iterator range_end = it;
#if CGAL_MESH_3_VERBOSE > 10
    std::cerr << "Point " << range_begin->first.first << " is duplicated, in "
              << "the following patches:\n";
#endif // CGAL_MESH_3_VERBOSE
    typename Union_find_t::handle first_handle =
      handles[range_begin->second];
    // In a second loop on the equal-range, update new_ids around p
    for (it = boost::next(range_begin); it != range_end; ++it)
    {
#if CGAL_MESH_3_VERBOSE > 10
      std::cerr << " - #" << it->second << "\n";
#endif // CGAL_MESH_3_VERBOSE
      union_find.unify_sets(first_handle,
                            handles[it->second]);
    }
    pit = range_end;
  }
  std::vector<Patch_id> new_ids(this->nb_of_patch_plus_one());
  for (std::size_t i = 0; i < this->nb_of_patch_plus_one(); ++i) {
    new_ids[i] = *union_find.find(handles[i]);
  }
  reindex_patches(new_ids);
}

template < typename GT_, typename P_, typename TA_>
void
Polyhedral_complex_mesh_domain_3<GT_,P_,TA_>::
add_features_from_split_graph_into_polylines(Featured_edges_copy_graph& g_copy)
{
  std::vector<Polyline_with_context> polylines;

  Mesh_3::internal::Extract_polyline_with_context_visitor<
    Polyline_with_context,
    Featured_edges_copy_graph
    > visitor(g_copy, polylines);
  Mesh_3::internal::Angle_tester<GT_> angle_tester;
  split_graph_into_polylines(g_copy, visitor, angle_tester);

  this->add_features_with_context(polylines.begin(),
                                  polylines.end());

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  {//DEBUG
    std::ofstream og("polylines_graph.polylines.txt");
    og.precision(17);
    for(const Polyline_with_context& poly : polylines)
    {
      og << poly.polyline_content.size() << " ";
      for(const Point_3& p : poly.polyline_content)
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
  typedef Polyhedron_type Polyhedron;
  typedef Patch_id P_id;
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

  typedef typename boost::property_map<Polyhedron,vertex_point_t>::const_type Vpm;
  Vpm vpm = get(vertex_point, p);
  for(Graph_vertex_descriptor v : make_range(vertices(graph))){
    vertex_descriptor vc;
    typename P2vmap::iterator it = p2vmap.find(get(vpm,v));
    if(it == p2vmap.end()) {
      vc = add_vertex(g_copy);
      g_copy[vc] = get(vpm, v);
      p2vmap[get(vpm,v)] = vc;
    }
  }

  typedef typename boost::property_map<Polyhedron,face_patch_id_t<P_id> >::type Face_patch_id_pmap;
  Face_patch_id_pmap fpm = get(face_patch_id_t<P_id>(),p);

  for(Graph_edge_descriptor e : make_range(edges(graph))){
    vertex_descriptor vs = p2vmap[get(vpm,source(e,graph))];
    vertex_descriptor vt = p2vmap[get(vpm,target(e,graph))];
    CGAL_warning_msg(vs != vt, "ignore self loop");
    if(vs != vt) {
      const std::pair<edge_descriptor, bool> pair = add_edge(vs,vt,g_copy);
      typename boost::graph_traits<Polyhedron>::halfedge_descriptor he = halfedge(e, p);
      if(!is_border(he, p)) {
        g_copy[pair.first].insert(get(fpm, face(he,p)));
      }
      he = opposite(he,p);
      if(!is_border(he, p)) {
        g_copy[pair.first].insert(get(fpm, face(he,p)));
      }
    }
  }

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  {// DEBUG
    Mesh_3::internal::dump_graph_edges("edges-graph.polylines.txt", g_copy);
  }
#endif
}

template < typename GT_, typename P_, typename TA_>
template < typename Surf_p_index>
void
Polyhedral_complex_mesh_domain_3<GT_,P_,TA_>::
reindex_patches(const std::vector<Surf_p_index>& map)
{
#if CGAL_MESH_3_VERBOSE > 10
  for(Surf_p_index i = 0, end = Surf_p_index(map.size()); i < end; ++i) {
    if(i != map[i])
      std::cerr << "patch #" << i << " is reindexed to " << map[i] << "\n";
  }
#endif // CGAL_MESH_3_VERBOSE
  typedef typename boost::graph_traits<Polyhedron_type>::face_descriptor
    Face_descriptor;
  for(std::size_t i = 0, end = stored_polyhedra.size(); i < end; ++i) {
    Polyhedron_type& poly = stored_polyhedra[i];
    typename boost::property_map<Polyhedron_type,
                                 face_patch_id_t<Patch_id> >::type
      face_pid_pmap = get(face_patch_id_t<Patch_id>(), poly);
    for(Face_descriptor fd : faces(poly))
    {
      const Surf_p_index id = get(face_pid_pmap, fd);
      const Surf_p_index new_id = map[id];
      put(face_pid_pmap, fd, new_id);
#if CGAL_MESH_3_VERBOSE > 11
      if(id != new_id) {
        std::cerr << "On polyhedron #" << i << ", patch #" << id
                  << " turned into patch #" << new_id << "\n";
      }
#endif // CGAL_MESH_3_VERBOSE
    }
  }
  for(Surface_patch_index& id :
                boundary_patches_ids)
  {
    id = map[id];
  }
  Base::reindex_patches(map);
}
/// @endcond

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYHEDRAL_COMPLEX_MESH_DOMAIN_3_H
