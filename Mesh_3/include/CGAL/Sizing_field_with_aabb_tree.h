// Copyright (c) 2010-2012-2016 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_MESH_3_SIZING_FIELD_WITH_AABB_TREE_H
#define CGAL_MESH_3_SIZING_FIELD_WITH_AABB_TREE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Profile_counter.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Mesh_3/Protect_edges_sizing_field.h> // for weight_modifier

#include <cstddef>
#include <memory>
#include <limits>

#include <CGAL/Mesh_3/experimental/Facet_patch_id_map.h>
#include <CGAL/Mesh_3/experimental/Get_curve_index.h>
#include <CGAL/Mesh_3/experimental/AABB_filtered_projection_traits.h>

#include <boost/container/flat_set.hpp>
#if defined(CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY) || defined(CGAL_NO_ASSERTIONS) == 0
#  include <sstream>
#  include <boost/format.hpp>
#endif  // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY || (! CGAL_NO_ASSERTIONS)

namespace CGAL {
/*!
* @ingroup PkgMesh3DomainFields
*
* The class `Sizing_field_with_aabb_tree` is a model of concept `MeshDomainField_3`.
* It provides a sizing field to be used in the meshing process of a polyhedral surface with features.
*
* At each query point `p`, the field value is the minimum of the distances to
*   - the surface patches `p` does not belong to, and
*   - for each curve `C` it does belong to, the intersection points of `C` with
*     the plane orthogonal to `C` and passing through `p`.
*
* This sizing field is designed to be used for the edges of the mesh complex,
* in `Mesh_edge_criteria_3` or as `edge_size` in `Mesh_criteria_3`.
* Using this sizing field for complex edges provides a high-quality hint
* to the protecting balls placement algorithm, since it ensures that no
* protecting ball will intersect a surface patch to which the
* corresponding vertex does not belong.
*
* @tparam GT is the geometric traits class. It must match the type `Triangulation::Geom_traits`,
* where `Triangulation` is the nested type of the model of `MeshComplex_3InTriangulation_3` used
* in the meshing process.
* @tparam MeshDomain is the type of the domain. It must be a model of `MeshDomainWithFeatures_3`,
* and derive from `CGAL::Polyhedral_mesh_domain_with_features_3`


* @cgalModels{MeshDomainField_3}
*/
template <typename GT,
          typename MeshDomain
#ifndef DOXYGEN_RUNNING
        , typename Input_facets_AABB_tree_ = CGAL::Default
        , typename Get_curve_index_ = CGAL::Default
        , typename Facet_patch_id_map_ = CGAL::Default
#endif
          >
struct Sizing_field_with_aabb_tree
{
  using FT      = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Index   = typename MeshDomain::Index;

private:
  typedef typename MeshDomain::Corner_index        Corner_index;
  typedef typename MeshDomain::Curve_index         Curve_index;
  typedef typename MeshDomain::Surface_patch_index Patch_index;

  typedef boost::container::flat_set<Curve_index> Curves_ids;
  typedef boost::container::flat_set<Patch_index> Patches_ids;

  typedef std::map<Point_3, Corner_index> Corners_indices;
  typedef std::vector<Point_3>     Corners;
  typedef std::vector<Curves_ids>  Corners_incident_curves;
  typedef std::vector<Patches_ids> Corners_incident_patches;
  typedef std::vector<Patches_ids> Curves_incident_patches;

  typedef GT Kernel_;
  typedef CGAL::Delaunay_triangulation_3<Kernel_> Dt;

  using Input_facets_AABB_tree = typename CGAL::Default::Get<
    Input_facets_AABB_tree_,
    typename MeshDomain::AABB_tree
    >::type;
  using AABB_traits                                 = typename Input_facets_AABB_tree::AABB_traits;
  using Point_and_primitive_id                      = typename AABB_traits::Point_and_primitive_id;
  typedef typename Input_facets_AABB_tree::Primitive  Input_facets_AABB_tree_primitive_;
  typedef typename MeshDomain::Curves_AABB_tree       Input_curves_AABB_tree_;
  typedef typename Input_curves_AABB_tree_::Primitive Input_curves_AABB_tree_primitive_;

  using Get_curve_index = typename CGAL::Default::Get<
    Get_curve_index_,
    CGAL::Mesh_3::Get_curve_index<typename MeshDomain::Curves_AABB_tree::Primitive>
    >::type;
  using Facet_patch_id_map = typename CGAL::Default::Get<
    Facet_patch_id_map_,
    CGAL::Mesh_3::Facet_patch_id_map<MeshDomain,
                                     typename Input_facets_AABB_tree::Primitive>
    >::type;

public:
  /// \name Creation
  /// @{

  /*!
  * Constructor for the sizing field.
  * @param d the maximal value returned by `operator()`, corresponding
  * to an upper bound on feature edge lengths in the output mesh.
  * @param domain the mesh domain to be meshed
  */
  Sizing_field_with_aabb_tree(const FT& d, const MeshDomain& domain)
    : Sizing_field_with_aabb_tree(d,
                                  domain,
                                  domain.aabb_tree(),
                                  Get_curve_index(),
                                  Facet_patch_id_map())
  {}

  /// @}

  struct Face_ids_traversal_traits {
    using Limits = std::numeric_limits<Patch_index>;
    Patch_index min{(Limits::max)()};
    Patch_index max{Limits::lowest()};
    Facet_patch_id_map index_map{};

    bool go_further() const { return true; }
    template <typename T> bool do_intersect(std::nullptr_t, const T&) { return true; }

    template <typename P> void intersection(std::nullptr_t, const P& primitive) {
      const Patch_index id = get(index_map, primitive.id());
      if (id < min) min = id;
      if (id > max) max = id;
    }
  };

#ifndef DOXYGEN_RUNNING
  Sizing_field_with_aabb_tree
  (const typename Kernel_::FT d,
   const MeshDomain& domain,
   const Input_facets_AABB_tree& aabb_tree,
   Get_curve_index get_curve_index = Get_curve_index(),
   Facet_patch_id_map facet_patch_id_map = Facet_patch_id_map()
   )
    : d_ptr{std::make_shared<Private_data>(d,
                                           domain,
                                           aabb_tree,
                                           get_curve_index,
                                           facet_patch_id_map)}
  {
    if(!d_ptr->aabb_tree.empty()) {
      // Initialize the per-patch kd-trees
      //  - First compute the number of patches and min/max patch ids
      Face_ids_traversal_traits traversal_traits;
      d_ptr->aabb_tree.traversal(nullptr, traversal_traits);
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      std::cerr << "min: " << traversal_traits.min << ", max: " << traversal_traits.max << '\n';
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      d_ptr->min_patch_id = traversal_traits.min;
      d_ptr->kd_trees_ptrs.resize(traversal_traits.max - d_ptr->min_patch_id + 1);

      using Node = CGAL::AABB_node<AABB_traits>;
      using Primitive = typename AABB_traits::Primitive;
      using Point_and_primitive_id = typename AABB_traits::Point_and_primitive_id;

      std::vector<std::set<Point_and_primitive_id>> vertices_per_patch;
      vertices_per_patch.resize(d_ptr->kd_trees_ptrs.size());

      //  - then fill sets of vertices per patch
      auto push_vertex = [&](const auto& primitive) {
        const Patch_index id = get(d_ptr->facet_patch_id_map, primitive.id());
        const Point_3& p = primitive.reference_point();
        vertices_per_patch[id - d_ptr->min_patch_id].emplace(p, primitive.id());
      };

      struct Visit_all_primitives {
        decltype(push_vertex) register_vertex;
        bool go_further() const { return true; }
        bool do_intersect(std::nullptr_t, const Node&) { return true; }
        void intersection(std::nullptr_t, const Primitive& primitive)
        {
          register_vertex(primitive);
        }
      } visit_all_primitives{push_vertex};
      d_ptr->aabb_tree.traversal(nullptr, visit_all_primitives);

      //  - then create the kd-trees
      for(std::size_t i = 0; i < vertices_per_patch.size(); ++i) {
        if(!vertices_per_patch[i].empty()) {
          d_ptr->kd_trees_ptrs[i] = std::make_unique<Kd_tree>(vertices_per_patch[i].begin(),
                                                              vertices_per_patch[i].end());
        }
      }
    }
    {
      Corner_index maximal_corner_index = 0;
      typedef std::pair<Corner_index, Point_3> Corner_index_and_point;
      std::vector<Corner_index_and_point > corners_tmp;
      d_ptr->domain.get_corners(std::back_inserter(corners_tmp));
      for(const Corner_index_and_point& pair : corners_tmp)
      {
        // Fill `corners_indices`
        if(pair.first > maximal_corner_index) maximal_corner_index = pair.first;
        d_ptr->corners_indices.
          insert(typename Corners_indices::value_type(pair.second, pair.first));
      }
      d_ptr->corners.resize(maximal_corner_index+1);
      for(const Corner_index_and_point& pair : corners_tmp)
      {
        // Fill `corners`
        d_ptr->corners[pair.first] = pair.second;
      }
    }

    //fill incidences of corners with curves
    d_ptr->corners_incident_curves.resize(d_ptr->corners.size());
    for(const typename Corners_indices::value_type& pair : d_ptr->corners_indices) {
      d_ptr->dt.insert(pair.first);

      // Fill `corners_incident_curves[corner_id]`
      Curves_ids& incident_curves = d_ptr->corners_incident_curves[pair.second];
      d_ptr->domain.get_corner_incident_curves(pair.second,
                                        std::inserter(incident_curves,
                                                      incident_curves.end()));
      // For each incident loops, insert a point on the loop, as far as
      // possible.
      for(Curve_index curve_index : incident_curves) {
        if(domain.is_loop(curve_index)) {
          FT curve_length = d_ptr->domain.curve_length(curve_index);
          Point_3 other_point =
            d_ptr->domain.construct_point_on_curve(pair.first,
                                                   curve_index,
                                                   curve_length / 2);
          d_ptr->dt.insert(other_point);
        }
      }
    }

    if (d_ptr->aabb_tree.empty())
      return;//there is no surface --> no surface patches

    //fill incidences with patches
    d_ptr->curves_incident_patches.resize(domain.maximal_curve_index()+1);
    d_ptr->corners_incident_patches.resize(d_ptr->corners.size());
    for (Curve_index curve_id = 1;
         std::size_t(curve_id) < d_ptr->curves_incident_patches.size();
         ++curve_id)
    {
      const typename MeshDomain::Surface_patch_index_set& incident_patches =
        d_ptr->domain.get_incidences(curve_id);

      // Fill `curves_incident_patches[curve_id]`
      d_ptr->curves_incident_patches[curve_id].
        insert(boost::container::ordered_unique_range,
               incident_patches.begin(), incident_patches.end());
    }
    for(const typename Corners_indices::value_type& pair : d_ptr->corners_indices) {
      // Get `corners_incident_curves[corner_id]`
      Curves_ids& incident_curves = d_ptr->corners_incident_curves[pair.second];

      // Fill `corners_incident_patches[corner_id]`
      Patches_ids& incident_patches = d_ptr->corners_incident_patches[pair.second];
      for(Curve_index curve_index : incident_curves) {
        const Patches_ids& curve_incident_patches =
          d_ptr->curves_incident_patches[curve_index];
        incident_patches.insert(boost::container::ordered_unique_range,
                                curve_incident_patches.begin(),
                                curve_incident_patches.end());
      }
    }
  }

  std::optional<Point_and_primitive_id>
  closest_point_on_surfaces(const Point_3& p,
                            const Patches_ids& patch_ids_to_ignore) const
  {
    std::optional<Point_and_primitive_id> result{};
    if(d_ptr->aabb_tree.empty()) return result;
    for(std::size_t i = 0; i < d_ptr->kd_trees_ptrs.size(); ++i) {
      const auto patch_id = static_cast<Patch_index>(i + d_ptr->min_patch_id);
      if(patch_ids_to_ignore.find(patch_id) != patch_ids_to_ignore.end()) continue;
      if(d_ptr->kd_trees_ptrs[i]) {
        const Kd_tree& kd_tree = *d_ptr->kd_trees_ptrs[i];
        const Point_and_primitive_id closest_point = kd_tree.closest_point(p);
        if(!result || CGAL::compare_distance(p, closest_point.first, result->first) == CGAL::SMALLER)
        {
          result = closest_point;
        }
      }
    }
    if(!result) return result;
    CGAL::Mesh_3::Filtered_projection_traits<
      AABB_traits,
      Facet_patch_id_map
      > projection_traits(*result,
                          patch_ids_to_ignore.begin(), patch_ids_to_ignore.end(),
                          d_ptr->aabb_tree.traits(),
                          d_ptr->facet_patch_id_map);
    d_ptr->aabb_tree.traversal(p, projection_traits);
    result = projection_traits.closest_point_and_primitive();
    return result;
  }
#endif // DOXYGEN_RUNNING

public:
  /// \name Operations
  /// @{

  /*!
  * returns the value of the sizing field at the point `p`,
  * assumed to be included in the input complex feature with dimension `dimension`
  * and mesh subcomplex index `id`.
  */
  double operator()(const Point_3& p,
                    const int dim,
                    const Index& id) const
  {
    CGAL_PROFILER("Sizing field");
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
    if(dim <= 1) {
      std::cerr << "Sizing("  << p << ", dim=" << dim
                << ", index=#" << CGAL::IO::oformat(id) << "): ";
    }
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
    FT result = d_ptr->d_;
    if(dim == 0) {
      if(d_ptr->dt.dimension() < 1) {
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        std::cerr << result << "(dt.dimension() < 1)\n";
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        return result;
      }
      typename Dt::Point nearest;

      typename Dt::Locate_type lt;
      int li, lj;
      const typename Dt::Cell_handle ch = d_ptr->dt.locate(p, lt, li, lj);
      if(lt == Dt::VERTEX) {
//         std::cerr << "lt == Dt::VERTEX\n";
        const typename Dt::Vertex_handle vh = ch->vertex(li);
        std::vector<typename Dt::Vertex_handle> vs;
        vs.reserve(32);
        d_ptr->dt.finite_adjacent_vertices(vh, std::back_inserter(vs));
        CGAL_assertion(!vs.empty());
        nearest = d_ptr->dt.point(vs[0]);
//         std::cerr << "sq_dist = " << CGAL::squared_distance(p, nearest)
//                   << std::endl;
        typename Kernel_::Compare_distance_3 compare_dist;
        for (typename std::vector<typename Dt::Vertex_handle>::const_iterator
               it = vs.begin(); it != vs.end(); ++it)
        {
//           std::cerr << "sq_dist = " << CGAL::squared_distance(p, dt.point(*it))
//                   << std::endl;
          if(compare_dist(p, d_ptr->dt.point(*it), nearest) == CGAL::SMALLER) {
//             std::cerr << "  nearest!\n";
            nearest = d_ptr->dt.point(*it);
          }
        }
      } else {
//         std::cerr << "lt=" << lt << std::endl;
        const typename Dt::Vertex_handle vh = d_ptr->dt.nearest_vertex(p, ch);
        nearest = d_ptr->dt.point(vh);
      }
      const FT dist = CGAL_NTS sqrt(CGAL::squared_distance( nearest, p));
      // std::cerr << (std::min)(dist / FT(1.5), d_) << "\n";
      result = (std::min)(dist * FT(0.5), result);

      // now search in the AABB tree
      typename Corners_indices::const_iterator ids_it =
          d_ptr->corners_indices.find(p);
      if(ids_it == d_ptr->corners_indices.end()) {
        std::cerr << "ERROR at " << __FILE__ << " line " << __LINE__ << "\n";
      }
      else if(!d_ptr->aabb_tree.empty()) {
        const Patches_ids& ids = d_ptr->corners_incident_patches[ids_it->second];
        CGAL_assertion(! ids.empty());

        const auto closest_point_and_primitive = closest_point_on_surfaces(p, ids);

        if(closest_point_and_primitive != std::nullopt) {
          result =
            (std::min)(FT(0.9 / CGAL::sqrt(CGAL::Mesh_3::internal::weight_modifier) *
                       CGAL_NTS
                       sqrt(CGAL::squared_distance(p,
                                                   closest_point_and_primitive->first))),
                       result);
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
          {
            std::stringstream s;

            using boost::io::group;
            using std::setprecision;
            s << boost::format("\nSizing field is %1% at point (%2%)"
                               " on corner!\n"
                               "Closest face id: %3%\n"
                               "Closest point: %4%\n"
                               "Ids are { ")
              % group(setprecision(17),result)
              % group(setprecision(17),p)
              % CGAL::IO::oformat(get(d_ptr->facet_patch_id_map,
                                  closest_point_and_primitive->second))
              % group(setprecision(17),
                      closest_point_and_primitive->first);
            for(Patch_index i : ids) {
              s << CGAL::IO::oformat(i) << " ";
            }
            s << "}\n";
            std::cerr << s.str();
            std::cerr << result << " (result)\n";
          }
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        } else {
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
          std::cerr << result << " (projection not found)\n";
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        }
      } // end if(corners.find(p) != corners.end()) and !aabb_tree.empty()
    } // end if(dim == 0)
    else if(dim != 1) {
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      std::cerr << result << "\n";
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      return result;
    }
    else { // dim == 1
      const typename MeshDomain::Curve_index& curve_id =
        d_ptr->domain.curve_index(id);
      const Patches_ids& ids = d_ptr->curves_incident_patches[curve_id];
      CGAL_assertion(! ids.empty());

      if(!d_ptr->aabb_tree.empty()) {
        //Compute distance to surface patches
        const auto closest_point_and_primitive = closest_point_on_surfaces(p, ids);
        if(closest_point_and_primitive == std::nullopt) {
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
          std::cerr << result << " (projection not found) ids:";
          for(Patch_index i : ids) {
            std::cerr << " " << CGAL::IO::oformat(i);
          }
          std::cerr << "\n";

#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
          return result;
        }

        CGAL_assertion(ids.count(get(d_ptr->facet_patch_id_map,
                                     closest_point_and_primitive->second)) == 0);

        result =
          (std::min)(FT(0.9 / CGAL::sqrt(CGAL::Mesh_3::internal::weight_modifier) *
                     CGAL_NTS
                     sqrt(CGAL::squared_distance(p,
                                                 closest_point_and_primitive->first))),
                     result);

#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        {
          std::stringstream s;

          s << boost::format("\nSizing field is %1% at point (%2%)"
                             " on curve #%3% !\n"
                             "Closest face id: %4%\n"
                             "Ids are { ")
            % result % p % curve_id
            % CGAL::IO::oformat(get(d_ptr->facet_patch_id_map,
                                closest_point_and_primitive->second));
          for(Patch_index i : ids) {
            s << CGAL::IO::oformat(i) << " ";
          }
          s << "}\n";
          std::cerr << s.str();
          std::cerr << result << " (result)\n";
        }
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
#ifndef CGAL_NO_ASSERTIONS
        if(result <= 0) {
          std::stringstream s;

          s << boost::format("Sizing field is %1% at point (%2%)"
                             " on curve #%3% !\n"
                             "Closest face id: %4%\n"
                             "Ids are { ")
            % result % p % curve_id
            % CGAL::IO::oformat(get(d_ptr->facet_patch_id_map,
                                closest_point_and_primitive->second));
          for(Patch_index i : ids) {
            s << CGAL::IO::oformat(i) << " ";
          }
          s << "}\n";
          CGAL_assertion_msg(result <=0, s.str().c_str());
        }
#ifdef PROTECTION_DEBUG
        else if (result <= (d_ / 1e7)) {
          std::stringstream s;
          s << boost::format("Sizing field is %1% at point (%2%)"
                             " on curve #%3% !\n"
                             "Closest face id: %4%\n"
                             "Ids are { ")
            % result % p % curve_id
            % closest_point_and_primitive->second->patch_id();
          for(Patch_index i : ids) {
            s << CGAL::IO::oformat(i) << " ";
          }
          s << "}\n";
          std::cerr << "ERROR at " << __FILE__ << " line " << __LINE__ << " :\n"
                    << s.str() << std::endl;
        }
#endif // PROTECTION_DEBUG
#endif // CGAL_NO_ASSERTIONS
      } // end if(!aabb_tree.empty())
      //Compute distance to the curves, and exclude the one on which p lies
      CGAL::Mesh_3::Filtered_projection_traits<typename Input_curves_AABB_tree_::AABB_traits,
                                               Get_curve_index >
        curves_projection_traits(curve_id,
                                 d_ptr->domain.curves_aabb_tree().traits(),
                                 d_ptr->get_curve_index);

      d_ptr->domain.curves_aabb_tree().traversal(p, curves_projection_traits);

      //Compute distance to the curve on which p lies
      typedef typename GT::Segment_3                        Segment_3;
      typedef typename GT::Plane_3                          Plane_3;

      const typename Input_curves_AABB_tree_::Point_and_primitive_id& ppid
        = d_ptr->domain.curves_aabb_tree().closest_point_and_primitive(p);

      Segment_3 curr_segment(*ppid.second.second, *(ppid.second.second + 1));
      Plane_3 curr_ortho_plane(p, curr_segment.to_vector()/*normal*/);
      Input_curves_AABB_tree_primitive_ curr_prim(ppid.second);

      std::vector<Input_curves_AABB_tree_primitive_> prims;
      d_ptr->domain.curves_aabb_tree().
          all_intersected_primitives(curr_ortho_plane, std::back_inserter(prims));

#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      std::cerr << std::endl;
      std::cerr << "p = " << p << std::endl;
      std::cerr << "curr_ortho_plane = " << curr_ortho_plane << std::endl;
      std::cerr << "PRIMITIVES FOUND : " << prims.size() << std::endl;
#endif

      Point_3 closest_intersection;
      Input_curves_AABB_tree_primitive_ closest_primitive = prims[0];
      FT sqd_intersection = -1;
      for(Input_curves_AABB_tree_primitive_ prim : prims)
      {
        if (prim.id() == curr_prim.id())
          continue;//curr_prim is the closest primitive

        if (curve_id != prim.id().first->first)
          continue;//don't deal with the same curves as what is done above

        const auto int_res = CGAL::intersection(prim.datum(), curr_ortho_plane);
        if (int_res)
        {
          if (const Point_3* pp = std::get_if<Point_3>(&*int_res))
          {
            FT new_sqd = CGAL::squared_distance(p, *pp);
            FT dist = CGAL::abs(d_ptr->domain.signed_geodesic_distance(p, *pp, curve_id));

#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
            std::cerr << "Intersection point : Point_3(" << *pp << ") ";
            std::cerr << "\n  new_sqd = " << new_sqd ;
            std::cerr << "\n  dist = " << dist << "\n";
#endif
            if (new_sqd * 1e10 < CGAL::squared_distance(curr_segment.source(),
                                                        curr_segment.target()))
            {
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
              std::cerr << "  too close, compared to possible rounding errors, "
                        << "SKIPPED\n";
#endif
              continue;
            }
            if (CGAL_NTS sqrt(new_sqd) > 0.9 * dist)
              continue;
            if (sqd_intersection == -1 || new_sqd < sqd_intersection)
            {
              sqd_intersection = new_sqd;
              closest_intersection = *pp;
              closest_primitive = prim;
            }
          }
          else
            continue;// intersection is a segment : collinear case
        }
        else
          CGAL_assertion(false);//prim was returned as an intersected primitive
      }

      //compare closest_projection and closest_intersection, and keep the closest
      if (curves_projection_traits.found())
      {
        FT tmp_sqd = CGAL::squared_distance(p, curves_projection_traits.closest_point());
        if (sqd_intersection == -1 || tmp_sqd < sqd_intersection)
        {
          sqd_intersection = tmp_sqd;
          closest_intersection = curves_projection_traits.closest_point();
          closest_primitive = Input_curves_AABB_tree_primitive_(
            curves_projection_traits.closest_point_and_primitive().second);
        }
      }
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      std::cout << " curve_id = " << curve_id
                << " proj_cid = " << closest_primitive.id().first->first
                << " (" << get(d_ptr->get_curve_index, closest_primitive.id()) << ")"
                << std::endl;
      std::cerr << " --- domain.curves_aabb_tree().traversal \n";
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      if (sqd_intersection > 0)
      {
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        std::cerr << "FOUND!\n";
        std::cerr << "  closest_point: " << closest_intersection << "\n"
                  << "  distance = " << CGAL_NTS sqrt(sqd_intersection) << std::endl;
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        FT new_result =
          (std::min)(FT(0.45 / CGAL::sqrt(CGAL::Mesh_3::internal::weight_modifier) *
                     CGAL_NTS sqrt(sqd_intersection)),
                     d_ptr->d_);

#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        std::cerr << "result     = " << result << "\n";
        std::cerr << "new_result = " << new_result << "\n";
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        if(result > new_result) {
          result = new_result;
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
          std::stringstream s;

          s << boost::format("\nSizing field is %1% at point (%2%)"
                             " on curve #%3% !\n"
                             "Closest CURVE id: %4%\n"
                             "Ids are { ")
            % result % p % curve_id
            % closest_primitive.id().first->first;
          for(int i : ids) {
            s << i << " ";
          }
          s << "}\n";
          std::cerr << s.str();
          std::cerr << result << " (result)\n";
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        }
      }
    } // end dim == 1
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
    std::cerr << result << std::endl;
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
    return result;
  }

  /// @}

private:
  using Kd_tree = CGAL::AABB_search_tree<AABB_traits>;
  struct Private_data {
    using FT = typename Kernel_::FT;
    Private_data(FT d,
                 const MeshDomain& domain,
                 const Input_facets_AABB_tree& aabb_tree,
                 Get_curve_index get_curve_index,
                 Facet_patch_id_map facet_patch_id_map)
      : d_(d)
      , domain(domain)
      , aabb_tree(aabb_tree)
      , get_curve_index(get_curve_index)
      , facet_patch_id_map(facet_patch_id_map)
    {}
    FT d_;
    const MeshDomain& domain;
    const Input_facets_AABB_tree& aabb_tree;
    Get_curve_index     get_curve_index;
    Facet_patch_id_map  facet_patch_id_map;
    Dt dt{};
    Corners          corners{};
    Corners_indices  corners_indices{};
    Curves_incident_patches   curves_incident_patches{};
    Corners_incident_patches  corners_incident_patches{};
    Corners_incident_curves   corners_incident_curves{};
    std::vector<std::unique_ptr<Kd_tree>> kd_trees_ptrs{};
    Patch_index min_patch_id{};
  };
  std::shared_ptr<Private_data> d_ptr;
};

}//end namespace CGAL

#endif // CGAL_MESH_3_SIZING_FIELD_WITH_AABB_TREE_H
