// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_AUTOREFINEMENT_H
#define CGAL_POLYGON_MESH_PROCESSING_AUTOREFINEMENT_H

#include <CGAL/license/Polygon_mesh_processing/geometric_repair.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// output
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#ifdef CGAL_DEBUG_PMP_AUTOREFINE
// debug
#include <CGAL/Surface_mesh.h>
#endif

#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

#ifndef DOXYGEN_RUNNING
namespace autorefine_impl {

template <class EK>
void generate_subtriangles(const typename EK::Triangle_3& t,
                           const std::vector<typename EK::Segment_3>& segments,
                           const std::vector<typename EK::Point_3>& points,
                           std::vector<typename EK::Triangle_3>& new_triangles)
{
  typedef CGAL::Projection_traits_3<EK> P_traits;
  typedef CGAL::Exact_intersections_tag Itag;
  //typedef CGAL::No_constraint_intersection_requiring_constructions_tag Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, Default, Itag> CDT_base;
  typedef CGAL::Constrained_triangulation_plus_2<CDT_base> CDT;

  typename EK::Vector_3 n = normal(t[0], t[1], t[2]);
  //~ bool orientation_flipped = false;
  //~ if (n.x() < 0)
  //~ {
    //~ orientation_flipped = true;
    //~ n = -n;
  //~ }
  //~ else
  //~ {
    //~ if (n.x()==0)
    //~ {
      //~ if (n.y() < 0)
      //~ {
        //~ orientation_flipped = true;
        //~ n = -n;
      //~ }
      //~ else
      //~ {
        //~ if (n.y()==0)
        //~ {
          //~ if (n.z() < 0)
          //~ {
            //~ orientation_flipped = true;
            //~ n = -n;
          //~ }
        //~ }
      //~ }
    //~ }
  //~ }

  P_traits cdt_traits(n);
  CDT cdt(cdt_traits);

  cdt.insert_outside_affine_hull(t[0]);
  cdt.insert_outside_affine_hull(t[1]);
  typename CDT::Vertex_handle v = cdt.tds().insert_dim_up(cdt.infinite_vertex(), false);
  v->set_point(t[2]);

  for (const typename EK::Segment_3& s : segments)
    cdt.insert_constraint(s[0], s[1]);

  for (const typename EK::Point_3& p : points)
    cdt.insert(p);

  //~ if (orientation_flipped)
    //~ for (typename CDT::Face_handle fh : cdt.finite_face_handles())
    //~ {
      //~ new_triangles.emplace_back(fh->vertex(0)->point(),
                                 //~ fh->vertex(cdt.cw(0))->point(),
                                 //~ fh->vertex(cdt.ccw(0))->point());
    //~ }
  //~ else
#ifdef CGAL_DEBUG_PMP_AUTOREFINE_DUMP_TRIANGULATIONS
    static int k = 0;
    std::stringstream buffer;
    buffer.precision(17);
    int nbt=0;
#endif
    for (typename CDT::Face_handle fh : cdt.finite_face_handles())
    {
      new_triangles.emplace_back(fh->vertex(0)->point(),
                                 fh->vertex(cdt.ccw(0))->point(),
                                 fh->vertex(cdt.cw(0))->point());
#ifdef CGAL_DEBUG_PMP_AUTOREFINE_DUMP_TRIANGULATIONS
      ++nbt;
      buffer << fh->vertex(0)->point() << "\n";
      buffer << fh->vertex(cdt.ccw(0))->point() << "\n";
      buffer << fh->vertex(cdt.cw(0))->point() << "\n";
#endif
    }

#ifdef CGAL_DEBUG_PMP_AUTOREFINE_DUMP_TRIANGULATIONS
    std::ofstream dump("triangulation_"+std::to_string(k)+".off");
    dump << "OFF\n" << 3*nbt << " " << nbt << " 0\n";
    dump << buffer.str();
    for (int i=0; i<nbt; ++i)
      dump << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
    ++k;
#endif
}

template <class EK>
struct Intersection_visitor
{
  std::vector<typename EK::Segment_3>& all_segments_1;
  std::vector<typename EK::Segment_3>& all_segments_2;
  std::vector<typename EK::Point_3>& all_points_1;
  std::vector<typename EK::Point_3>& all_points_2;

  Intersection_visitor(std::vector<typename EK::Segment_3>& all_segments_1,
                       std::vector<typename EK::Segment_3>& all_segments_2,
                       std::vector<typename EK::Point_3>& all_points_1,
                       std::vector<typename EK::Point_3>& all_points_2)
    : all_segments_1(all_segments_1)
    , all_segments_2(all_segments_2)
    , all_points_1(all_points_1)
    , all_points_2(all_points_2)
  {}

  typedef void result_type;
  void operator()(const typename EK::Point_3& p)
  {
    all_points_1.push_back(p);
    all_points_2.push_back(p);
  }

  void operator()(const typename EK::Segment_3& s)
  {
    all_segments_1.push_back(s);
    all_segments_2.push_back(s);
  }

  void operator()(const typename EK::Triangle_3& t)
  {
    for (std::size_t i=0; i<3; ++i)
    {
      typename EK::Segment_3 s(t[i], t[(i+1)%3]);
      all_segments_1.push_back(s);
      all_segments_2.push_back(s);
    }

  }

  void operator()(const std::vector<typename EK::Point_3>& poly)
  {
    std::size_t nbp = poly.size();
    for (std::size_t i=0; i<nbp; ++i)
    {
      typename EK::Segment_3 s(poly[i], poly[(i+1)%nbp]);
      all_segments_1.push_back(s);
      all_segments_2.push_back(s);
    }
  }
};

template <class EK>
bool is_output_valid(std::vector<Point_3<EK>> soup_points,
                     std::vector<std::array<std::size_t, 3> > soup_triangles)
{
  typedef Surface_mesh<typename EK::Point_3> Exact_mesh;
  Exact_mesh etm;
  orient_polygon_soup(soup_points, soup_triangles);
  polygon_soup_to_polygon_mesh(soup_points, soup_triangles, etm);

#ifdef CGAL_DEBUG_PMP_AUTOREFINE
  std::cerr << std::setprecision(17);
  auto verbose_fail_msg = [&](int i, const std::pair<typename Exact_mesh::Face_index, typename Exact_mesh::Face_index>& p)
  {
    typename Exact_mesh::Halfedge_index h1 = halfedge(p.first, etm), h2 = halfedge(p.second, etm);
    std::cerr << "DEBUG: failing at check #" << i << "\n";
    std::cerr << "DEBUG: " << etm.point(source(h1, etm)) << " " <<  etm.point(target(h1, etm)) << " " <<  etm.point(target(next(h1, etm), etm)) << " " << etm.point(source(h1, etm)) << "\n";
    std::cerr << "DEBUG: " << etm.point(source(h2, etm)) << " " <<  etm.point(target(h2, etm)) << " " <<  etm.point(target(next(h2, etm), etm)) << " " << etm.point(source(h2, etm)) << "\n";
  };
#endif

  //TODO: double check me
  auto skip_faces = [&](const std::pair<typename Exact_mesh::Face_index, typename Exact_mesh::Face_index>& p)
  {
    typename Exact_mesh::Halfedge_index h1 = etm.halfedge(p.first), h2=etm.halfedge(p.second);

    boost::container::small_vector<typename Exact_mesh::Halfedge_index, 3> v1;
    v1.push_back(prev(h1, etm));
    v1.push_back(h1);
    v1.push_back(next(h1, etm));

    boost::container::small_vector<typename Exact_mesh::Halfedge_index, 3> v2;
    v2.push_back(prev(h2, etm));
    v2.push_back(h2);
    v2.push_back(next(h2, etm));

    //collect identical vertices
    boost::container::small_vector<std::pair<typename Exact_mesh::Halfedge_index, typename Exact_mesh::Halfedge_index>, 3> common;
    for(typename Exact_mesh::Halfedge_index h1 : v1)
      for(typename Exact_mesh::Halfedge_index h2 : v2)
        if (etm.point(target(h1, etm))==etm.point(target(h2,etm)))
          common.push_back(std::make_pair(h1,h2));

    if (common.empty())
    {
#ifdef CGAL_DEBUG_PMP_AUTOREFINE
      verbose_fail_msg(1, p);
#endif
      return false;
    }

    switch (common.size())
    {
      case 1:
      {
        // geometric check if the opposite segments intersect the triangles
        const typename EK::Triangle_3 t1(etm.point(source(h1,etm)),
                                         etm.point(target(h1,etm)),
                                         etm.point(target(next(h1,etm),etm)));
        const typename EK::Triangle_3 t2(etm.point(source(h2,etm)),
                                         etm.point(target(h2,etm)),
                                         etm.point(target(next(h2,etm),etm)));

        const typename EK::Segment_3 s1(etm.point(source(common[0].first,etm)), etm.point(target(next(common[0].first,etm),etm)));
        const typename EK::Segment_3 s2(etm.point(source(common[0].second,etm)), etm.point(target(next(common[0].second,etm),etm)));

        if(do_intersect(t1, s2) || do_intersect(t2, s1))
        {
#ifdef CGAL_DEBUG_PMP_AUTOREFINE
          verbose_fail_msg(2, p);
#endif
          return false;
        }
        return true;
      }
      case 2:
      {
        // shared edge
        h1 = next(common[0].first, etm) == common[1].first ? common[1].first : common[0].first;
        h2 = next(common[0].second, etm) == common[1].second ? common[1].second : common[0].second;

        if( CGAL::coplanar(etm.point(source(h1,etm)),
                           etm.point(target(h1,etm)),
                           etm.point(target(etm.next(h1),etm)),
                           etm.point(target(etm.next(h2),etm))) &&
            CGAL::coplanar_orientation(etm.point(source(h1,etm)),
                           etm.point(target(h1,etm)),
                           etm.point(target(etm.next(h1),etm)),
                           etm.point(target(etm.next(h2),etm))) == CGAL::POSITIVE)
        {
#ifdef CGAL_DEBUG_PMP_AUTOREFINE
          verbose_fail_msg(3, p);
#endif
          return false;
        }
        return true;
      }
      default: // size == 3
        return true;
    }
  };

  std::vector< std::pair<typename Exact_mesh::Face_index, typename Exact_mesh::Face_index> > si_faces;
  self_intersections(etm,
                     CGAL::filter_output_iterator(std::back_inserter(si_faces),
                     skip_faces));

  return si_faces.empty();
}

} // end of autorefine_impl

template <class TriangleMesh, class Point_3, class NamedParameters = parameters::Default_named_parameters>
void autorefine_soup_output(const TriangleMesh& tm,
                            std::vector<Point_3>& soup_points,
                            std::vector<std::array<std::size_t, 3> >& soup_triangles,
                            const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, tm));

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef std::pair<face_descriptor, face_descriptor> Pair_of_faces;

  std::vector<Pair_of_faces> si_pairs;

  // collect intersecting pairs of triangles
  self_intersections<Concurrency_tag>(tm, std::back_inserter(si_pairs), np);

  if (si_pairs.empty()) return;

  // mark degenerate faces so that we can ignore them
  typedef CGAL::dynamic_face_property_t<bool> Degen_property_tag;
  typedef typename boost::property_map<TriangleMesh, Degen_property_tag>::const_type Is_degen_map;
  Is_degen_map is_degen = get(Degen_property_tag(), tm);

  for(face_descriptor f : faces(tm))
    put(is_degen, f, is_degenerate_triangle_face(f, tm, np));

  // assign an id per triangle involved in an intersection
  // + the faces involved in the intersection
  typedef CGAL::dynamic_face_property_t<int> TID_property_tag;
  typedef typename boost::property_map<TriangleMesh, TID_property_tag>::const_type Triangle_id_map;

  Triangle_id_map tid_map = get(TID_property_tag(), tm);
  for (face_descriptor f : faces(tm))
    put(tid_map, f, -1);

  std::vector<face_descriptor> intersected_faces;
  int tid=-1;
  for (const Pair_of_faces& p : si_pairs)
  {
    if (get(tid_map, p.first)==-1 && !get(is_degen, p.first))
    {
      put(tid_map, p.first, ++tid);
      intersected_faces.push_back(p.first);
    }
    if (get(tid_map, p.second)==-1 && !get(is_degen, p.second))
    {
      put(tid_map, p.second, ++tid);
      intersected_faces.push_back(p.second);
    }
  }

  // init the vector of triangles used for the autorefinement of triangles
  typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
  std::vector< EK::Triangle_3 > triangles(tid+1);
  Cartesian_converter<GT, EK> to_exact;

  for(face_descriptor f : intersected_faces)
  {
    halfedge_descriptor h = halfedge(f, tm);
    triangles[get(tid_map, f)]= EK::Triangle_3(
      to_exact( get(vpm, source(h, tm)) ),
      to_exact( get(vpm, target(h, tm)) ),
      to_exact( get(vpm, target(next(h, tm), tm)) ) );
  }

  std::vector< std::vector<EK::Segment_3> > all_segments(triangles.size());
  std::vector< std::vector<EK::Point_3> > all_points(triangles.size());


  typename EK::Intersect_3 intersection = EK().intersect_3_object();
  for (const Pair_of_faces& p : si_pairs)
  {
    int i1 = get(tid_map, p.first),
        i2 = get(tid_map, p.second);

    if (i1==-1 || i2==-1) continue; //skip degenerate faces

    const EK::Triangle_3& t1 = triangles[i1];
    const EK::Triangle_3& t2 = triangles[i2];

    auto inter = intersection(t1, t2);

    if (inter != boost::none)
    {
      autorefine_impl::Intersection_visitor<EK> intersection_visitor(all_segments[i1],  all_segments[i2],
                                                                     all_points[i1], all_points[i2]);

      boost::apply_visitor(intersection_visitor, *inter);
    }
  }

  // now refine triangles
  std::vector<EK::Triangle_3> new_triangles;
  for(std::size_t ti=0; ti<triangles.size(); ++ti)
  {
    if (all_segments[ti].empty() && all_points[ti].empty())
      new_triangles.push_back(triangles[ti]);
    else
      autorefine_impl::generate_subtriangles<EK>(triangles[ti], all_segments[ti], all_points[ti], new_triangles);
  }


  // brute force output: create a soup, orient and to-mesh
  Cartesian_converter<EK, GT> to_input;
  std::map<EK::Point_3, std::size_t> point_id_map;
#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
  std::vector<EK::Point_3> exact_soup_points;
#endif

  for (vertex_descriptor v : vertices(tm))
  {
    if (point_id_map.insert(std::make_pair(to_exact(get(vpm,v)), soup_points.size())).second)
    {
      soup_points.push_back(get(vpm,v));
#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
      exact_soup_points.push_back(to_exact(get(vpm,v)));
#endif
    }
  }

  auto get_point_id = [&](const typename EK::Point_3& pt)
  {
    auto insert_res = point_id_map.insert(std::make_pair(pt, soup_points.size()));
    if (insert_res.second)
    {
      soup_points.push_back(to_input(pt));
#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
      exact_soup_points.push_back(pt);
#endif
    }
    return insert_res.first->second;
  };

  for (face_descriptor f : faces(tm))
  {
    if (get(is_degen, f)) continue; //skip degenerate faces

    int tid = get(tid_map, f);
    if (tid == -1)
    {
      halfedge_descriptor h = halfedge(f, tm);
      soup_triangles.emplace_back(
        CGAL::make_array(get_point_id(to_exact(get(vpm,source(h, tm)))),
                         get_point_id(to_exact(get(vpm,target(h, tm)))),
                         get_point_id(to_exact(get(vpm,target(next(h, tm), tm)))))
      );
    }
  }
  for (const typename EK::Triangle_3& t : new_triangles)
  {
    soup_triangles.emplace_back(CGAL::make_array(get_point_id(t[0]), get_point_id(t[1]), get_point_id(t[2])));
  }


#ifndef CGAL_NDEBUG
  CGAL_assertion( autorefine_impl::is_output_valid(exact_soup_points, soup_triangles) );
#endif

#ifdef CGAL_DEBUG_PMP_AUTOREFINE
  autorefine_impl::is_output_valid(exact_soup_points, soup_triangles);
    throw std::runtime_error("invalid output");
#endif

}
#endif

/**
 * \ingroup PMP_corefinement_grp
 * \link coref_def_subsec autorefines \endlink `tm`. Refines a triangle mesh
 * so that no triangles intersects in their interior.
 * Self-intersection edges will be marked as constrained. If an edge that was marked as
 * constrained is split, its sub-edges will be marked as constrained as well.
 *
 * @tparam TriangleMesh a model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tm input triangulated surface mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalParamNBegin{geom_traits}
 *   \cgalParamDescription{an instance of a geometric traits class}
 *   \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *   \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *   \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 * \cgalParamNEnd
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *  \cgalParamNEnd
 *
 *   \cgalParamNBegin{edge_is_constrained_map}
 *     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
 *                    as key type and `bool` as value type}
 *     \cgalParamDefault{a constant property map returning `false` for any edge}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{face_index_map}
 *     \cgalParamDescription{a property map associating to each face of `tm` a unique index between `0` and `num_faces(tm) - 1`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{an automatically indexed internal map}
 *     \cgalParamExtra{If the property map is writable, the indices of the faces of `tm1` and `tm2`
 *                     will be set after the corefinement is done.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{visitor}
 *     \cgalParamDescription{a visitor used to track the creation of new faces}
 *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
 *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
void
autorefine(      TriangleMesh& tm,
           const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  std::vector<typename GT::Point_3> soup_points;
  std::vector<std::array<std::size_t, 3> > soup_triangles;

  autorefine_soup_output(tm, soup_points, soup_triangles, np);

  clear(tm);
  orient_polygon_soup(soup_points, soup_triangles);
  polygon_soup_to_polygon_mesh(soup_points, soup_triangles, tm);
}


} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_AUTOREFINEMENT_H
