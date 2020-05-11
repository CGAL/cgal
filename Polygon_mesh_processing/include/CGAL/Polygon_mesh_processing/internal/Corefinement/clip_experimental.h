// Copyright (c) 2019-2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé
//                 Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CLIP_EXPERIMENTAL_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CLIP_EXPERIMENTAL_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/assertions.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Exact_kernel_selector.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename PolygonMesh>
void output_mesh(const char* filename, const PolygonMesh& pmesh)
{
  std::ofstream output(filename);
  output.precision(17);
  output << pmesh;
  output.close();
}

template <typename TriangleMesh, typename EVPM,
          typename NamedParameters1, typename NamedParameters2>
bool clip_mesh_exactly(TriangleMesh& cc,
                       EVPM cc_evpm,
                       TriangleMesh clipper, // intentional copy
                       const NamedParameters1& np_cc,
                       const NamedParameters2& np_c)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;

  // These typedefs concern the VPM of the CLIPPER mesh
  typedef typename GetK<TriangleMesh, NamedParameters2>::Kernel                 Clipper_kernel;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters2>::type      Clipper_VPM;

#if 1
  typedef CGAL::Exact_predicates_exact_constructions_kernel                     Exact_kernel;
  typedef typename Exact_kernel::Point_3                                        Exact_point;

  typedef typename CGAL::Cartesian_converter<Clipper_kernel, Exact_kernel>      C2E;
#else
  typedef typename CGAL::Exact_kernel_selector<Clipper_kernel>::Exact_kernel    Exact_kernel;
  typedef typename Exact_kernel::Point_3                                        Exact_point;

  typedef typename CGAL::Exact_kernel_selector<Clipper_kernel>::C2E             C2E;
#endif

  typedef CGAL::dynamic_vertex_property_t<Exact_point>                          EP_property_tag;

  // must have equal exact point types because of corefinement.h
  CGAL_static_assertion((std::is_same<EVPM, typename boost::property_map<TriangleMesh, EP_property_tag>::type>::value));

  using parameters::get_parameter;
  using parameters::choose_parameter;

  // Initialize an exact vpm for the clipper
  Clipper_VPM clipper_vpm = choose_parameter(get_parameter(np_c, internal_np::vertex_point),
                                             get_property_map(CGAL::vertex_point, clipper));

  C2E to_exact;
  EVPM clipper_evpm = get(EP_property_tag(), clipper);
  for(vertex_descriptor vd : vertices(clipper))
    put(clipper_evpm, vd, to_exact(get(clipper_vpm, vd)));

#ifdef CGAL_DEBUG_CLIPPING
  const bool valid_input = is_valid(cc) && is_valid(clipper) &&
                           !does_self_intersect(cc, parameters::vertex_point_map(cc_evpm)) &&
                           !does_self_intersect(clipper, parameters::vertex_point_map(clipper_evpm)) &&
                           does_bound_a_volume(clipper, parameters::vertex_point_map(clipper_evpm));
  if(!valid_input)
  {
    std::cerr << "Invalid input for clip()" << std::endl;
    std::cerr << "is cc valid: " << cc.is_valid() << std::endl;
    std::cerr << "is clipper valid: " << clipper.is_valid() << std::endl;
    std::cerr << "does part self intersect? " << does_self_intersect(cc, parameters::vertex_point_map(cc_evpm)) << std::endl;
    std::cerr << "does clipper self intersect? " << does_self_intersect(clipper, parameters::vertex_point_map(clipper_evpm)) << std::endl;
    std::cerr << "clipper bounds a volume? " << does_bound_a_volume(clipper) << std::endl;
    return false;
  }
#endif

  bool clip_volume = choose_parameter(get_parameter(np_cc, internal_np::clip_volume), true);
  bool use_compact_clipper = choose_parameter(get_parameter(np_cc, internal_np::use_compact_clipper), false);

  // @todo is it possible to forward clip_volume/use_compact_clipper/visitor?
  clip_volume = false;
  use_compact_clipper = false;

  bool res = clip(cc, clipper,
                  parameters::vertex_point_map(cc_evpm)
                             .clip_volume(clip_volume)
                             .use_compact_clipper(use_compact_clipper)
                             .throw_on_self_intersection(true),
                  parameters::vertex_point_map(clipper_evpm));

  CGAL_postcondition(is_valid_polygon_mesh(cc));
  CGAL_postcondition(is_valid_polygon_mesh(clipper));

  return res;
}

template <typename TriangleMesh, typename VPM, typename EVPM>
void fill_triangle_mesh(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
                        const TriangleMesh& tm,
                        VPM tm_vpm,
                        EVPM tm_evpm,
                        TriangleMesh& si_face,
                        EVPM si_face_evpm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;

  vertex_descriptor vd0 = add_vertex(si_face);
  vertex_descriptor vd1 = add_vertex(si_face);
  vertex_descriptor vd2 = add_vertex(si_face);

  // Note that we don't even fill the internal VPM
  put(si_face_evpm, vd0, get(tm_evpm, source(hd, tm)));
  put(si_face_evpm, vd1, get(tm_evpm, target(hd, tm)));
  put(si_face_evpm, vd2, get(tm_evpm, target(next(hd, tm), tm)));

  CGAL::Euler::add_face(std::initializer_list<vertex_descriptor>{vd0, vd1, vd2}, si_face);
}

template <typename TriangleMesh, typename FacePairRange,
          typename NamedParameters1, typename NamedParameters2>
bool clip_single_self_intersecting_cc(const FacePairRange& self_intersecting_faces,
                                      TriangleMesh& cc,
                                      TriangleMesh& clipper,
                                      const NamedParameters1& np_cc,
                                      const NamedParameters2& np_c)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor           face_descriptor;

  typedef typename GetK<TriangleMesh, NamedParameters1>::Kernel                 CC_kernel;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters1>::type      CC_VPM;

#if 1
  typedef CGAL::Exact_predicates_exact_constructions_kernel                     Exact_kernel;
  typedef typename Exact_kernel::Point_3                                        Exact_point;

  typedef typename CGAL::Cartesian_converter<CC_kernel, Exact_kernel>           C2E;
  typedef typename CGAL::Cartesian_converter<Exact_kernel, CC_kernel>           E2C;
#else
  typedef typename CGAL::Exact_kernel_selector<CC_kernel>::Exact_kernel         Exact_kernel;
  typedef typename Exact_kernel::Point_3                                        Exact_point;

  typedef typename CGAL::Exact_kernel_selector<CC_kernel>::C2E                  C2E;
  typedef typename CGAL::Exact_kernel_selector<CC_kernel>::E2C                  E2C;
#endif

  typedef CGAL::dynamic_vertex_property_t<Exact_point>                          EP_property_tag;
  typedef typename boost::property_map<TriangleMesh, EP_property_tag>::type     EVPM;

  typedef std::pair<face_descriptor, face_descriptor>                           Pair_of_faces;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  CGAL_precondition(does_self_intersect(cc));

#ifdef CGAL_DEBUG_CLIPPING
  std::cout << "CC self-intersects" << std::endl;
#endif

  bool res = true;

  // 1. Gather the problematic faces
  // 2. Copy (one-by-one) the problematic faces into independent meshes
  // 3. Remove the problematic faces from 'cc'
  // 4. Clip (cc - problematic_faces)
  // 2bis. Clip the problematic faces
  // 5. Re-add the clipped problematic faces from step 2.

  // 1. ------------------------------
  std::set<face_descriptor> si_faces;
  for(const Pair_of_faces& fp : self_intersecting_faces)
  {
    si_faces.insert(fp.first);
    si_faces.insert(fp.second);
  }

  std::size_t independent_faces_n = si_faces.size();
#ifdef CGAL_DEBUG_CLIPPING
  std::cout << independent_faces_n << " independent faces" << std::endl;
#endif

  // 2. ------------------------------
  // Switch to an exact VPM
  CC_VPM cc_vpm = choose_parameter(get_parameter(np_cc, internal_np::vertex_point),
                                   get_property_map(CGAL::vertex_point, cc));

  C2E to_exact;
  EVPM cc_evpm = get(EP_property_tag(), cc);
  for(vertex_descriptor vd : vertices(cc))
    put(cc_evpm, vd, to_exact(get(cc_vpm, vd)));

  // forced to do that because clip() wants both meshes to be of the same type (rather than just requiring
  // that the two VPMs have matching types)
  std::vector<TriangleMesh> clipped_si_faces;
  clipped_si_faces.reserve(independent_faces_n);
  std::vector<EVPM> si_face_evpms;
  si_face_evpms.reserve(independent_faces_n);

  std::size_t counter = 0;
  for(face_descriptor fd : si_faces)
  {
#ifdef CGAL_DEBUG_CLIPPING
    std::cout << "  Build single mesh for face: " << fd << " (" << counter+1 << "/" << independent_faces_n << ")" << std::endl;
#endif

    const halfedge_descriptor hd = halfedge(fd, cc);
    if(!is_degenerate_triangle_face(face(hd, cc), cc))
    {
      ++counter;

      clipped_si_faces.resize(counter);
      TriangleMesh& si_face = clipped_si_faces.back();
      si_face_evpms.resize(counter);
      EVPM& si_face_evpm = si_face_evpms.back();
      si_face_evpm = get(EP_property_tag(), si_face);

      fill_triangle_mesh(hd, cc, cc_vpm, cc_evpm, si_face, si_face_evpm);
    }
  }
#ifdef CGAL_DEBUG_CLIPPING
  std::cout << clipped_si_faces.size() << " problematic faces to clip" << std::endl;
#endif

  // 3. ------------------------------
  for(const face_descriptor fd : si_faces)
    CGAL::Euler::remove_face(halfedge(fd, cc), cc);

  // 3bis. iteratively remove faces incident to non-manifold vertices
  std::vector<halfedge_descriptor> non_manifold_cones;
  CGAL::Polygon_mesh_processing::non_manifold_vertices(cc, std::back_inserter(non_manifold_cones));
  std::set<vertex_descriptor> nm_vertices;

  std::unordered_set<face_descriptor> nm_faces;
  for(halfedge_descriptor nm_hd : non_manifold_cones)
  {
    // keep the faces incident to the nm vertex for only one (the first one) of the cones
    if(nm_vertices.insert(target(nm_hd, cc)).second)
      continue;

    for(halfedge_descriptor hd : CGAL::halfedges_around_target(nm_hd, cc))
      if(cc.face(hd) != face_descriptor())
        nm_faces.insert(cc.face(hd));
  }

  clipped_si_faces.reserve(clipped_si_faces.size() + nm_faces.size());
  for(face_descriptor fd : nm_faces)
  {
    clipped_si_faces.resize(clipped_si_faces.size() + 1);
    si_face_evpms.resize(si_face_evpms.size() + 1);

    TriangleMesh& si_face = clipped_si_faces.back();
    EVPM si_face_evpm = si_face_evpms.back();

    const halfedge_descriptor hd = halfedge(fd, cc);
    fill_triangle_mesh(hd, cc, cc_vpm, cc_evpm, si_face, si_face_evpm);
    CGAL::Euler::remove_face(hd, cc);
  }

  CGAL_postcondition(is_valid(cc));
  CGAL_postcondition(!does_self_intersect(cc));
#ifdef CGAL_DEBUG_CLIPPING
  std::cout << "Pruned mesh: " << vertices(cc).size() << " nv "  << faces(cc).size() << " nf " << std::endl;
  output_mesh("pruned_mesh.off", cc);
#endif

  // 4. ------------------------------
#ifdef CGAL_DEBUG_CLIPPING
  std::cout << "Clipping CC's main part (w/o self-intersecting faces)" << std::endl;
#endif

  if(!clip_mesh_exactly(cc, cc_evpm, clipper, np_cc, np_c))
    res = false;

#ifdef CGAL_DEBUG_CLIPPING
  output_mesh("pruned_mesh_clipped.off", cc);
#endif

  // 5. ------------------------------
  counter = 0;
  for(std::size_t i=0, csfs=clipped_si_faces.size(); i<csfs; ++i)
  {
    TriangleMesh& clipped_si_face = clipped_si_faces[i];
    EVPM si_face_evpm = si_face_evpms[i];

#ifdef CGAL_DEBUG_CLIPPING
    std::cout << "  Clipping face: (" << counter++ << "/" << independent_faces_n << ")" << std::endl;
    std::cout << "    Clipped: " << num_vertices(clipped_si_face) << " nv " << num_faces(clipped_si_face) << " nf" << std::endl;
    std::cout << "    Clipper: " << num_vertices(clipper) << " nv " << num_faces(clipper) << " nf" << std::endl;
#endif

    if(!clip_mesh_exactly(clipped_si_face, si_face_evpm, clipper, np_cc, np_c))
      res = false;

    if(is_empty(clipped_si_face))
      continue;

#ifdef CGAL_DEBUG_CLIPPING
    std::cout << "  Re-attaching part!" << std::endl;
    std::cout << "    " << num_vertices(clipped_si_face) << " nv " << num_faces(clipped_si_face) << " nf" << std::endl;
    std::cout << "    copying..." << std::endl;
#endif

    CGAL::copy_face_graph(clipped_si_face, cc,
                          parameters::vertex_point_map(si_face_evpm),
                          parameters::vertex_point_map(cc_evpm));

#ifdef CGAL_DEBUG_CLIPPING
    std::cout << "\t" << vertices(cc).size() << " nv " << faces(cc).size() << " nf" << std::endl;
    std::cout << "    stitching..." << std::endl;
#endif

    stitch_borders(cc, parameters::vertex_point_map(cc_evpm)); // @todo do something local

#ifdef CGAL_DEBUG_CLIPPING
    std::cout << "\t" << vertices(cc).size() << " nv " << faces(cc).size() << " nf" << std::endl;
#endif
  }

  // 5bis. (Bring everyone back to the real world)
  E2C to_input;
  for(vertex_descriptor vd : vertices(cc))
    put(cc_vpm, vd, to_input(get(cc_evpm, vd)));

#ifdef CGAL_DEBUG_CLIPPING
  std::cout << "Done with this CC!" << std::endl;
#endif

  return res;
}

template <typename TriangleMesh,
          typename NamedParameters1, typename NamedParameters2>
bool clip_single_cc(TriangleMesh& cc,
                    TriangleMesh& clipper,
                    const NamedParameters1& np_cc,
                    const NamedParameters2& np_c)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor           face_descriptor;
  typedef std::pair<face_descriptor, face_descriptor>                           Pair_of_faces;

#ifdef CGAL_DEBUG_CLIPPING
  std::cout << "Clip single CC of size " << vertices(cc).size() << " nv " << faces(cc).size() << " nf" << std::endl;
#endif

  std::vector<Pair_of_faces> intersecting_faces;
  self_intersections(cc, std::back_inserter(intersecting_faces), np_cc);

  if(intersecting_faces.empty())
    return clip(cc, clipper, np_cc, np_c);
  else
    return clip_single_self_intersecting_cc(intersecting_faces, cc, clipper, np_cc, np_c);
}

} // namespace internal

namespace experimental {

// Try to split the mesh in multiple connected components and clip them independently
template <typename TriangleMesh,
          typename NamedParameters1, typename NamedParameters2>
bool clip_self_intersecting_mesh(TriangleMesh& tm,
                                 TriangleMesh& clipper,
                                 const NamedParameters1& np_tm,
                                 const NamedParameters2& np_c)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor           face_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters1>::type      VPM;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  typedef CGAL::dynamic_vertex_property_t<Point>                                P_property_tag;
  typedef typename boost::property_map<TriangleMesh, P_property_tag>::type      DVPM;

  typedef typename boost::graph_traits<TriangleMesh>::faces_size_type           faces_size_type;
  typedef typename boost::property_map<TriangleMesh,
                                       face_patch_id_t<faces_size_type> >::type PatchIDMap;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  // @todo keep this?
  if(!does_self_intersect(tm, np_tm))
    return clip(tm, clipper, np_tm, np_c);

  // Splitting connected components avoids having to treat a lot of "easy" self-intersections
  // (for example from two spheres intersecting)
  PatchIDMap pidmap = get(CGAL::face_patch_id_t<faces_size_type>(), tm);
  faces_size_type num_cc = connected_components(tm, pidmap);

#ifdef CGAL_DEBUG_CLIPPING
  std::cout << num_cc << " connected component(s)" << std::endl;
#endif

  // order ccs by size
  std::vector<int> number_of_faces_in_cc(num_cc, 0);
  for(face_descriptor f : faces(tm))
    number_of_faces_in_cc[get(pidmap, f)] += 1;

  std::vector<int> ordered_cc_ids(num_cc);
  std::iota(std::begin(ordered_cc_ids), std::end(ordered_cc_ids), 0);
  std::sort(std::begin(ordered_cc_ids), std::end(ordered_cc_ids),
            [&number_of_faces_in_cc](const int i, const int j) -> bool {
              return (number_of_faces_in_cc[i] > number_of_faces_in_cc[j]);
            });

  for(std::size_t i=0; i<ordered_cc_ids.size(); ++i)
    std::cout << "CC #" << i << " " << number_of_faces_in_cc[i] << " nf" << std::endl;

  std::vector<TriangleMesh> ccs(num_cc - 1); // largest CC stays in 'tm'
  std::vector<DVPM> cc_vpms(num_cc - 1);

  // Extract all but the largest CC
  for(faces_size_type i=1; i<num_cc; ++i)
  {
    cc_vpms[i] = get(P_property_tag(), ccs[i]);
    CGAL::Face_filtered_graph<TriangleMesh> tm_cc(tm, ordered_cc_ids[i], pidmap);
    CGAL::copy_face_graph(tm_cc, ccs[ordered_cc_ids[i]],
                          np_tm, parameters::vertex_point_map(cc_vpms[i]));
  }

  bool res = true;
  const CGAL::Bbox_3 clipper_bbox = CGAL::Polygon_mesh_processing::bbox(clipper, np_c);

  // Clip CC by CC
  keep_connected_components(tm, std::vector<faces_size_type>(1, ordered_cc_ids[0]), pidmap);

  const CGAL::Bbox_3 tm_bbox = CGAL::Polygon_mesh_processing::bbox(tm, np_tm);
  if(CGAL::do_overlap(tm_bbox, clipper_bbox))
  {
    if(!internal::clip_single_cc(tm, clipper, np_tm, np_c))
      res = false;
  }

  for(faces_size_type i=1; i<num_cc; ++i)
  {
    TriangleMesh& tm_cc = ccs[ordered_cc_ids[i]];
    DVPM cc_vpm = cc_vpms[ordered_cc_ids[i]];

    std::cout << "CC w/ " << num_vertices(tm_cc) << " nv " << num_faces(tm_cc) << " nf" << std::endl;

    const CGAL::Bbox_3 cc_bbox = CGAL::Polygon_mesh_processing::bbox(tm_cc, parameters::vertex_point_map(cc_vpm));

    if(CGAL::do_overlap(cc_bbox, clipper_bbox))
      if(!internal::clip_single_cc(tm_cc, clipper,  parameters::vertex_point_map(cc_vpm), np_c)) // @todo np forwarding
        res = false;

    CGAL::copy_face_graph(tm_cc, tm, parameters::vertex_point_map(cc_vpm), np_tm);
  }

  return res;
}

template <typename TriangleMesh,
          typename NamedParameters1, typename NamedParameters2>
bool generic_clip(TriangleMesh& tm,
                  TriangleMesh& clipper,
                  const NamedParameters1& np_tm,
                  const NamedParameters2& np_c)
{
  return clip_self_intersecting_mesh(tm, clipper, np_tm, np_c);
}

template <typename TriangleMesh, typename NamedParameters>
bool generic_clip(TriangleMesh& tm, TriangleMesh& clipper, const NamedParameters& np_tm)
{
  return generic_clip(tm, clipper, np_tm, parameters::all_default());
}

template <typename TriangleMesh>
bool generic_clip(TriangleMesh& tm, TriangleMesh& clipper)
{
  return generic_clip(tm, clipper, parameters::all_default(), parameters::all_default());
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL


#endif //CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CLIP_EXPERIMENTAL_H
