// Copyright (c) 2010-2012-2016 GeometryFactory Sarl (France).
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
// Author(s)     : Laurent Rineau
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************


#ifndef CGAL_MESH_3_SIZING_FIELD_WITH_AABB_TREE_H
#define CGAL_MESH_3_SIZING_FIELD_WITH_AABB_TREE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Profile_counter.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include "Get_facet_patch_id.h"
#include "Get_curve_index.h"
#include <CGAL/Mesh_3/Protect_edges_sizing_field.h> // for weight_modifier

#include <boost/foreach.hpp>
#if defined(CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY) || defined(CGAL_NO_ASSERTIONS) == 0
#  include <sstream>
#  include <boost/format.hpp>
#endif  // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY || (! CGAL_NO_ASSERTIONS)

#include "AABB_filtered_projection_traits.h"

template <typename GeomTraits, typename MeshDomain, typename Patches_ids,
          typename Input_facets_AABB_tree = typename MeshDomain::AABB_tree,
          typename Get_curve_index_ = CGAL::Default,
          typename Get_facet_patch_id_ = CGAL::Default
          >
struct Sizing_field_with_aabb_tree
{
  typedef GeomTraits Kernel_;
  typedef typename MeshDomain::Index Index;
  typedef std::vector<Patches_ids> Patches_ids_map;
  typedef CGAL::Delaunay_triangulation_3<Kernel_> Dt;
  typedef typename Kernel_::FT FT;
  typedef typename Kernel_::Point_3 Point_3;
  typedef std::map<Point_3, Patches_ids> Corners;

  typedef Input_facets_AABB_tree Input_facets_AABB_tree_;
  typedef typename Input_facets_AABB_tree_::Primitive Input_facets_AABB_tree_primitive_;
  typedef typename MeshDomain::Curves_AABB_tree Input_curves_AABB_tree_;
  typedef typename Input_curves_AABB_tree_::Primitive Input_curves_AABB_tree_primitive_;

  typedef typename CGAL::Default::Get<
    Get_curve_index_,
    CGAL::Mesh_3::Get_curve_index<typename MeshDomain::Curves_AABB_tree::Primitive>
    >::type Get_curve_index;
  typedef typename CGAL::Default::Get<
    Get_facet_patch_id_,
    CGAL::Mesh_3::Get_facet_patch_id<typename Input_facets_AABB_tree::Primitive>
    >::type Get_facet_patch_id;

  template <typename InputCornersIterator>
  Sizing_field_with_aabb_tree(typename Kernel_::FT d, 
                              InputCornersIterator corners_begin,
                              InputCornersIterator corners_end,
                              const Input_facets_AABB_tree_& aabb_tree,
                              const MeshDomain& domain,
                              const Patches_ids_map& patches_ids_map,
                              Get_curve_index get_curve_index
                              = Get_curve_index(),
                              Get_facet_patch_id get_facet_patch_id
                              = Get_facet_patch_id())
    : d_(d), aabb_tree(aabb_tree), 
      domain(domain), patches_ids_map(patches_ids_map),
      get_curve_index(get_curve_index),
      get_facet_patch_id(get_facet_patch_id),
      dt(),
      corners(corners_begin, corners_end)
  {
    typedef typename 
      std::iterator_traits<InputCornersIterator>::value_type Pair;

    BOOST_FOREACH(const Pair& p, std::make_pair(corners_begin, corners_end))
    {
      dt.insert(p.first);
    }
  }

  double operator()(const Point_3& p, 
                    const int dim, 
                    const Index& id) const
  { 
    CGAL_PROFILER("Sizing field");
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
    if(dim <= 1) {
      std::cerr << "Sizing("  << p << ", dim=" << dim 
                << ", index=#" << id << "): ";
    }
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
    double result;
    if(dim == 0) {
      if(dt.dimension() < 1) {
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        std::cerr << d_ << "(dt.dimension() < 1)\n";
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        return d_;
      }
      typename Dt::Point nearest;

      typename Dt::Locate_type lt;
      int li, lj;
      const typename Dt::Cell_handle ch = dt.locate(p, lt, li, lj);
      if(lt == Dt::VERTEX) {
// 	std::cerr << "lt == Dt::VERTEX\n";
        const typename Dt::Vertex_handle vh = ch->vertex(li);
        std::vector<typename Dt::Vertex_handle> vs;
        vs.reserve(32);
        dt.finite_adjacent_vertices(vh, std::back_inserter(vs));
        CGAL_assertion(!vs.empty());
        nearest = vs[0]->point();
// 	std::cerr << "sq_dist = " << CGAL::squared_distance(p, nearest)
// 		  << std::endl;
        typename Kernel_::Compare_distance_3 compare_dist;
        for (typename std::vector<typename Dt::Vertex_handle>::const_iterator
               it = vs.begin(); it != vs.end(); ++it) 
        {
// 	  std::cerr << "sq_dist = " << CGAL::squared_distance(p, (*it)->point())
// 		  << std::endl;
          if(compare_dist(p, (*it)->point(), nearest) == CGAL::SMALLER) {
// 	    std::cerr << "  nearest!\n";
            nearest =  (*it)->point();
          }
        }
      } else {
// 	std::cerr << "lt=" << lt << std::endl;
        const typename Dt::Vertex_handle vh = dt.nearest_vertex(p, ch);
        nearest = vh->point();
      }
      const FT dist = CGAL_NTS sqrt(CGAL::squared_distance( nearest, p));
      // std::cerr << (std::min)(dist / FT(1.5), d_) << "\n";
      result = (std::min)(dist / FT(2), d_);

      // now search in the AABB tree
      typename Corners::const_iterator ids_it = corners.find(p);
      if(ids_it == corners.end()) {
        std::cerr << "ERROR at " << __FILE__ << " line " << __LINE__ << "\n";
      }
      else {
        const Patches_ids& ids = ids_it->second;
        CGAL_assertion(! ids.empty());

        CGAL::Mesh_3::Filtered_projection_traits<
          typename Input_facets_AABB_tree_::AABB_traits,
          Get_facet_patch_id 
          > projection_traits(ids.begin(), ids.end(),
                              aabb_tree.traits(), 
                              get_facet_patch_id);

        aabb_tree.traversal(p, projection_traits);

        if(projection_traits.found()) {
          result = 
            (std::min)(0.9 / CGAL::sqrt(CGAL::Mesh_3::internal::weight_modifier) * 
                       CGAL_NTS 
                       sqrt(CGAL::squared_distance(p, 
                                                   projection_traits.closest_point())),
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
              % get(get_facet_patch_id,
                    projection_traits.closest_point_and_primitive().second)
              % group(setprecision(17),
                      projection_traits.closest_point_and_primitive().first);
            BOOST_FOREACH(int i, ids) {
              s << i << " ";
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
      }
    }
    else if(dim != 1) {
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      std::cerr << d_ << "\n";
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      return d_;
    } 
    else { // dim == 1
      const typename MeshDomain::Curve_segment_index& curve_id = 
        domain.curve_segment_index(id);
      const Patches_ids& ids = patches_ids_map[curve_id];

      CGAL_assertion(! ids.empty());

      //Compute distance to surface patches
      CGAL::Mesh_3::Filtered_projection_traits
      <
        typename Input_facets_AABB_tree_::AABB_traits
        , Get_facet_patch_id
        > projection_traits(ids.begin(), ids.end(),
                            aabb_tree.traits(),
                            get_facet_patch_id);

      aabb_tree.traversal(p, projection_traits);

      if(!projection_traits.found()) {
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        std::cerr << d_ << " (projection not found)\n";
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        return d_;
      }

      CGAL_assertion(ids.count(get(get_facet_patch_id,
                                   projection_traits.closest_point_and_primitive().second)) == 0);

      result = 
        (std::min)(0.9 / CGAL::sqrt(CGAL::Mesh_3::internal::weight_modifier) * 
                   CGAL_NTS 
                   sqrt(CGAL::squared_distance(p, 
                                               projection_traits.closest_point())),
                   d_);
      
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      {
        std::stringstream s;
       
        s << boost::format("\nSizing field is %1% at point (%2%)"
                           " on curve #%3% !\n"
                           "Closest face id: %4%\n"
                           "Ids are { ")
          % result % p % curve_id
          % get(get_facet_patch_id,
                projection_traits.closest_point_and_primitive().second);
        BOOST_FOREACH(int i, ids) {
          s << i << " ";
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
          % get(get_facet_patch_id,
                projection_traits.closest_point_and_primitive().second);
        BOOST_FOREACH(int i, ids) {
          s << i << " ";
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
          % get(get_facet_patch_id,
                projection_traits.closest_point_and_primitive().second);
        BOOST_FOREACH(int i, ids) {
          s << i << " ";
        }
        s << "}\n";
        std::cerr << "ERROR at " << __FILE__ << " line " << __LINE__ << " :\n"
                  << s.str() << std::endl;
      }
#endif // PROTECTION_DEBUG
#endif // CGAL_NO_ASSERTIONS

      //Compute distance to the curves, and exclude the one on which p lies
      CGAL::Mesh_3::Filtered_projection_traits<typename Input_curves_AABB_tree_::AABB_traits,
                                               Get_curve_index >
        curves_projection_traits(curve_id,
                                 domain.curves_aabb_tree().traits(),
                                 get_curve_index);

      domain.curves_aabb_tree().traversal(p, curves_projection_traits);

      //Compute distance to the curve on which p lies
      typedef typename GeomTraits::Segment_3                        Segment_3;
      typedef typename GeomTraits::Plane_3                          Plane_3;
      typedef typename CGAL::cpp11::result_of<
        typename GeomTraits::Intersect_3(Segment_3, Plane_3)>::type Intersection_result;

      const typename Input_curves_AABB_tree_::Point_and_primitive_id& ppid
        = domain.curves_aabb_tree().closest_point_and_primitive(p);

      Segment_3 curr_segment(*ppid.second.second, *(ppid.second.second + 1));
      Plane_3 curr_ortho_plane(p, curr_segment.to_vector()/*normal*/);
      Input_curves_AABB_tree_primitive_ curr_prim(ppid.second);

      std::vector<Input_curves_AABB_tree_primitive_> prims;
      domain.curves_aabb_tree().all_intersected_primitives(curr_ortho_plane,
                                                           std::back_inserter(prims));

#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
      std::cout << std::endl;
      std::cout << "P = " << p << std::endl;
      std::cout << "PRIMITIVES FOUND : " << prims.size() << std::endl;
#endif

      Point_3 closest_intersection;
      Input_curves_AABB_tree_primitive_ closest_primitive = prims[0];
      FT sqd_intersection = -1;
      BOOST_FOREACH(Input_curves_AABB_tree_primitive_ prim, prims)
      {
        if (prim.id() == curr_prim.id())
          continue;//curr_prim is the closest primitive

        if (curve_id != prim.id().first->first)
          continue;//don't deal with the same curves as what is done above

        Intersection_result int_res
          = CGAL::intersection(prim.datum(), curr_ortho_plane);
        if (int_res)
        {
          if (const Point_3* pp = boost::get<Point_3>(&*int_res))
          {
            FT new_sqd = CGAL::squared_distance(p, *pp);
            FT gdist = (std::min)(CGAL::abs(domain.geodesic_distance(p, *pp, curve_id)),
              CGAL::abs(domain.geodesic_distance(*pp, p, curve_id)));

#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
            std::cout << "Intersection point : Point_3(" << *pp << ") ";
            std::cout << " new_sqd = " << new_sqd ;
            std::cout << " gdist = " << gdist << std::endl;;
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
            if (CGAL_NTS sqrt(new_sqd) > 0.9 * gdist)
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
                << " (" << get(get_curve_index, closest_primitive.id()) << ")"
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
        double new_result = 
          (std::min)(0.45 / CGAL::sqrt(CGAL::Mesh_3::internal::weight_modifier) * 
                     CGAL_NTS sqrt(sqd_intersection),
                     d_);

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
          BOOST_FOREACH(int i, ids) {
            s << i << " ";
          }
          s << "}\n";
          std::cerr << s.str();
          std::cerr << result << " (result)\n";
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
        }
      }
    }
#ifdef CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
    std::cerr << result << std::endl;
#endif // CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY
    return result;
  }
private:
  typename Kernel_::FT d_;
  const Input_facets_AABB_tree_& aabb_tree;
  const MeshDomain& domain;
  const Patches_ids_map& patches_ids_map;
  Get_curve_index get_curve_index;
  Get_facet_patch_id get_facet_patch_id;
  Dt dt;
  Corners corners;
};

#endif // CGAL_MESH_3_SIZING_FIELD_WITH_AABB_TREE_H
