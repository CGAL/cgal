// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot, Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_CLIP_SELF_INTERSECTING_H
#define CGAL_POLYGON_MESH_PROCESSING_CLIP_SELF_INTERSECTING_H

#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/assertions.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <iostream>
#include <fstream>

namespace CGAL {
namespace Polygon_mesh_processing {

template <typename TriangleMesh, typename IEVPM,
          typename NamedParametersTM, typename NamedParametersC>
void clip_mesh_exactly_with_clipper_copy(TriangleMesh& tm,
                                         IEVPM tm_evpm,
                                         TriangleMesh clipper, // INTENTIONAL COPY
                                         const NamedParametersTM& np_tm,
                                         const NamedParametersC& np_c)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor        vertex_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, NamedParametersTM>::type        Geom_traits;
  typedef CGAL::Exact_predicates_exact_constructions_kernel                    EPECK;
  typedef EPECK::Point_3                                                       EPoint_3;

  typedef CGAL::dynamic_vertex_property_t<EPoint_3>                            VEPP;
  typedef typename boost::template property_map<TriangleMesh, VEPP>::type      EVPM;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParametersC>::type     C_VPM;
  C_VPM clipper_vpm = choose_parameter(get_parameter(np_c, internal_np::vertex_point),
                                       get_property_map(vertex_point, clipper));

  CGAL::Cartesian_converter<Geom_traits, EPECK> to_exact;

  EVPM clipper_evpm = get(VEPP(), clipper);
  for(vertex_descriptor vd : vertices(clipper))
    put(clipper_evpm, vd, to_exact(get(clipper_vpm, vd)));

#ifdef CGAL_DEBUG_CLIPPING
  const bool valid_input = tm.is_valid() && clipper.is_valid() &&
                           !does_self_intersect(tm, CGAL::parameters::vertex_point_map(tm_evpm)) &&
                           !does_self_intersect(clipper, CGAL::parameters::vertex_point_map(clipper_evpm)) &&
                           does_bound_a_volume(clipper);
  if(!valid_input)
  {
    std::cerr << "Invalid input for clip()" << std::endl;
    std::cerr << "is tm valid: " << tm.is_valid() << std::endl;
    std::cerr << "is clipper valid: " << clipper.is_valid() << std::endl;
    std::cerr << "does part self intersect? " << does_self_intersect(tm, CGAL::parameters::vertex_point_map(tm_evpm)) << std::endl;
    std::cerr << "does clipper self intersect? " << does_self_intersect(clipper, CGAL::parameters::vertex_point_map(clipper_evpm)) << std::endl;
    std::cerr << "clipper bounds a volume? " << does_bound_a_volume(clipper) << std::endl;
    write_polygon_mesh("results/bad_cc.off", tm, parameters::stream_precision(17));
    std::exit(1);
  }
#endif

  clip(tm, clipper,
       CGAL::parameters::vertex_point_map(tm_evpm).throw_on_self_intersection(true),
       CGAL::parameters::vertex_point_map(clipper_evpm));

  CGAL_postcondition(CGAL::is_valid_polygon_mesh(tm));
  CGAL_postcondition(tm.is_valid()); // @fixme SM
  CGAL_postcondition(clipper.is_valid());
}

template <typename TriangleMesh, typename VPM, typename EVPM>
void fill_triangle_mesh(const TriangleMesh& tm,
                        VPM tm_vpm,
                        EVPM tm_evpm,
                        typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
                        TriangleMesh& si_face)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor          face_descriptor;

  typedef CGAL::Exact_predicates_exact_constructions_kernel                    EPECK;
  typedef EPECK::Point_3                                                       EPoint_3;

  // @fixme SM
  auto si_face_evpm = si_face.template add_property_map<vertex_descriptor, EPoint_3>("v:evpm", EPoint_3()).first;

#if 1 // whether to use add_face() or make_triangle()
  VPM si_face_vpm = get_property_map(CGAL::vertex_point, si_face);

  vertex_descriptor vd0 = add_vertex(si_face);
  vertex_descriptor vd1 = add_vertex(si_face);
  vertex_descriptor vd2 = add_vertex(si_face);
  put(si_face_vpm, vd0, get(tm_vpm, source(hd, tm)));
  put(si_face_vpm, vd1, get(tm_vpm, target(hd, tm)));
  put(si_face_vpm, vd2, get(tm_vpm, target(next(hd, tm), tm)));
  put(si_face_evpm, vd0, get(tm_evpm, source(hd, tm)));
  put(si_face_evpm, vd1, get(tm_evpm, target(hd, tm)));
  put(si_face_evpm, vd2, get(tm_evpm, target(next(hd, tm), tm)));
  CGAL::Euler::add_face(std::initializer_list<vertex_descriptor>{vd0, vd1, vd2}, si_face);
#else
  Point_3 p0 = get(tm_vpm, source(hd, tm)),
          p1 = get(tm_vpm, target(hd, tm)),
          p2 = get(tm_vpm, target(next(hd, tm), tm)); // @fixme make sure that this correspondency is correct

  make_triangle(p0, p1, p2, si_face);
  face_descriptor si_fd = *(faces(si_face).begin());
  halfedge_descriptor si_hd = halfedge(si_fd, si_face);
  vertex_descriptor vd0 = source(si_hd, si_face);
  vertex_descriptor vd1 = target(si_hd, si_face);
  vertex_descriptor vd2 = target(next(si_hd, si_face), si_face);
  put(si_face_evpm, vd0, get(tm_evpm, source(hd, tm)));
  put(si_face_evpm, vd1, get(tm_evpm, target(hd, tm)));
  put(si_face_evpm, vd2, get(tm_evpm, target(next(hd, tm), tm)));
#endif
}

template <typename TriangleMesh, typename FacePairRange, typename Clipper,
          typename NamedParametersTM, typename NamedParametersC>
void clip_single_cc_self_intersecting_mesh(TriangleMesh& tm,
                                           const FacePairRange& self_intersecting_faces,
                                           Clipper& clipper,
                                           const NamedParametersTM& np_tm,
                                           const NamedParametersC& np_c)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor          face_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, NamedParametersTM>::type        Geom_traits;

  // @todo use a lazy exact vertex point map
  typedef CGAL::Exact_predicates_exact_constructions_kernel                    EPECK;
  typedef EPECK::Point_3                                                       EPoint_3;
  typedef CGAL::dynamic_vertex_property_t<EPoint_3>                            VEPP;
  typedef typename boost::template property_map<TriangleMesh, VEPP>::type      EVPM;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParametersTM>::type    TM_VPM;

  typedef CGAL::dynamic_edge_property_t<bool>                                  EBP;
  typedef typename boost::template property_map<TriangleMesh, EBP>::type       Ecm;

  CGAL_precondition(does_self_intersect(tm, np_tm));

  using parameters::get_parameter;
  using parameters::choose_parameter;

  TM_VPM tm_vpm = choose_parameter(get_parameter(np_tm, internal_np::vertex_point),
                                   get_property_map(vertex_point, tm));

#ifdef CGAL_DEBUG_CLIPPING
  std::cout << "CC self-intersects" << std::endl;
#endif

  // 1. Gather the problematic faces
  // 2. Copy (one-by-one) the problematic faces into independent meshes
  // 3. Remove the problematic faces from 'tm'
  // 4. Clip (tm - problematic_faces)
  // 2bis. Clip the problematic faces
  // 5. Re-add the clipped problematic faces from step 2.

  // 1. ------------------------------
  std::set<face_descriptor> si_faces;
  for(const auto& fp : self_intersecting_faces)
  {
    si_faces.insert(fp.first);
    si_faces.insert(fp.second);
  }

  std::size_t independent_faces_n = si_faces.size();
#ifdef CGAL_DEBUG_CLIPPING
  std::cout << independent_faces_n << " independent faces" << std::endl;
#endif

  // 2. ------------------------------
  // Switch to exact
  CGAL::Cartesian_converter<Geom_traits, EPECK> to_exact;
  CGAL::Cartesian_converter<EPECK, Geom_traits> to_input;

  EVPM tm_evpm = get(VEPP(), tm);
  for(vertex_descriptor vd : vertices(tm))
    put(tm_evpm, vd, to_exact(get(tm_vpm, vd)));

  // forced to do that because clip() wants both meshes to be of the same type (rather than just requiring
  // that the two VPMs have matching types)
  std::vector<TriangleMesh> clipped_si_faces;
  clipped_si_faces.reserve(independent_faces_n);

  std::size_t counter = 0;
  for(face_descriptor fd : si_faces)
  {
#ifdef CGAL_DEBUG_CLIPPING
    std::cout << "  Build single mesh for face: " << fd << " (" << counter+1 << "/" << independent_faces_n << ")" << std::endl;
#endif

    const halfedge_descriptor hd = halfedge(fd, tm);
    if(!is_degenerate_triangle_face(face(hd, tm), tm, np_tm))
    {
      clipped_si_faces.resize(counter + 1);
      TriangleMesh& si_face = clipped_si_faces.back();
      fill_triangle_mesh(tm, tm_vpm, tm_evpm, hd, si_face);
      ++counter;
    }
  }
#ifdef CGAL_DEBUG_CLIPPING
  std::cout << clipped_si_faces.size() << " problematic faces to clip" << std::endl;
#endif

  // 3. ------------------------------
  for(const face_descriptor fd : si_faces)
    CGAL::Euler::remove_face(halfedge(fd, tm), tm);

  // 3bis. iteratively remove faces incident to non-manifold vertices

  std::vector<halfedge_descriptor> non_manifold_cones;
  CGAL::Polygon_mesh_processing::non_manifold_vertices(tm, std::back_inserter(non_manifold_cones));
  std::set<vertex_descriptor> nm_vertices;

  std::unordered_set<face_descriptor> nm_faces;
  for(halfedge_descriptor nm_hd : non_manifold_cones)
  {
    // keep the faces incident to the nm vertex for only one (the first one) of the cones
    if(nm_vertices.insert(target(nm_hd, tm)).second)
      continue;

    for(halfedge_descriptor hd : CGAL::halfedges_around_target(nm_hd, tm))
      if(tm.face(hd) != face_descriptor())
        nm_faces.insert(tm.face(hd));
  }

  clipped_si_faces.reserve(clipped_si_faces.size() + nm_faces.size());
  for(face_descriptor fd : nm_faces)
  {
    clipped_si_faces.resize(clipped_si_faces.size() + 1);
    TriangleMesh& si_face = clipped_si_faces.back();
    halfedge_descriptor hd = tm.halfedge(fd);
    fill_triangle_mesh(tm, tm_vpm, tm_evpm, hd, si_face);
    CGAL::Euler::remove_face(hd, tm);
  }

  CGAL_postcondition(tm.is_valid());
  CGAL_postcondition(!does_self_intersect(tm, np_tm));

#ifdef CGAL_DEBUG_CLIPPING
  std::cout << "Pruned mesh: " << vertices(tm).size() << " nv "  << faces(tm).size() << " nf " << std::endl;
  write_polygon_mesh("pruned_mesh.off", tm, parameters::stream_precision(17));
#endif

  // 4. ------------------------------
#ifdef CGAL_DEBUG_CLIPPING
  std::cout << "Clipping CC's main part (w/o self-intersecting faces)" << std::endl;
#endif

  clip_mesh_exactly_with_clipper_copy(tm, tm_evpm, clipper, np_tm, np_c);

#ifdef CGAL_DEBUG_CLIPPING
  write_polygon_mesh("pruned_mesh_clipped.off", tm, parameters::stream_precision(17));
#endif

  // 5. ------------------------------
  counter = 0;
  for(TriangleMesh& clipped_si_face : clipped_si_faces)
  {
#ifdef CGAL_DEBUG_CLIPPING
    std::cout << "  Clipping face: (" << counter++ << "/" << independent_faces_n << ")" << std::endl;
#endif

    auto si_face_evpm = clipped_si_face.template property_map<vertex_descriptor, EPoint_3>("v:evpm").first;
    CGAL_assertion((clipped_si_face.template property_map<vertex_descriptor, EPoint_3>("v:evpm").second));

    clip_mesh_exactly_with_clipper_copy(clipped_si_face, si_face_evpm, clipper, np_tm, np_c);
    if(is_empty(clipped_si_face))
      continue;

#ifdef CGAL_DEBUG_CLIPPING
    std::cout << "  Re-attaching part!" << std::endl;
    std::cout << "    " << num_vertices(clipped_si_face) << " nv " << num_faces(clipped_si_face) << " nf" << std::endl;
    std::cout << "    copying..." << std::endl;
#endif

    CGAL::copy_face_graph(clipped_si_face, tm,
                          CGAL::parameters::vertex_point_map(si_face_evpm),
                          CGAL::parameters::vertex_point_map(tm_evpm));

#ifdef CGAL_DEBUG_CLIPPING
    std::cout << "\t" << vertices(tm).size() << " nv " << faces(tm).size() << " nf" << std::endl;
    std::cout << "    stitching..." << std::endl;
#endif

    stitch_borders(tm, CGAL::parameters::vertex_point_map(tm_evpm)); // @todo do something local

#ifdef CGAL_DEBUG_CLIPPING
    std::cout << "\t" << vertices(tm).size() << " nv " << faces(tm).size() << " nf" << std::endl;
#endif
  }

  // 5bis. (Bring everyone back to the real world)
  for(vertex_descriptor vd : vertices(tm))
    put(tm_vpm, vd, to_input(get(tm_evpm, vd)));

#ifdef CGAL_DEBUG_CLIPPING
  std::cout << "Done with this CC!" << std::endl;
#endif
}

template <typename TriangleMesh, typename Clipper,
          typename NamedParametersTM, typename NamedParametersC>
void clip_single_cc(TriangleMesh& tm,
                    Clipper clipper, // INTENTIONAL COPY
                    const NamedParametersTM& np_tm,
                    const NamedParametersC& np_c)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;
  typedef std::pair<face_descriptor, face_descriptor>                     Pair_of_faces;

#ifdef CGAL_DEBUG_CLIPPING
  std::cout << "Clip single CC of size " << vertices(tm).size() << " nv " << faces(tm).size() << " nf" << std::endl;
#endif
  std::vector<Pair_of_faces> intersecting_faces;
  self_intersections(tm, std::back_inserter(intersecting_faces), np_tm);

  if(intersecting_faces.empty())
    clip(tm, clipper, np_tm, np_c);
  else
    clip_single_cc_self_intersecting_mesh(tm, intersecting_faces, clipper, np_tm, np_c);
}

// Try to split the mesh in multiple connected components and clip them independently
template <typename TriangleMesh, typename Clipper,
          typename NamedParametersTM, typename NamedParametersC>
void clip_self_intersecting_mesh(TriangleMesh& tm,
                                 Clipper& clipper,
                                 const NamedParametersTM& np_tm,
                                 const NamedParametersC& np_c)
{
  typedef typename boost::property_map<TriangleMesh, CGAL::face_patch_id_t<int> >::type PatchIDMap;
  PatchIDMap pidmap = get(CGAL::face_patch_id_t<int>(), tm);

  int num_cc = connected_components(tm, pidmap, np_tm);

#ifdef CGAL_DEBUG_CLIPPING
  std::cout << num_cc << " connected component(s)" << std::endl;
#endif

  std::vector<TriangleMesh> ccs(num_cc);
  for(int i=0; i<num_cc; ++i)
  {
    CGAL::Face_filtered_graph<TriangleMesh> filtered_tm(tm, i, pidmap, np_tm);
    CGAL::copy_face_graph(filtered_tm, ccs[i]);
  }

  // @fixme need to build the 'np_tm_cc' named parameter from the ground up (not to include the VPM)

  clear(tm);
  for(TriangleMesh& tm_cc : ccs)
  {
#ifndef CGAL_USE_CLIPPER_PLANE
    CGAL::Bbox_3 b1 = CGAL::Polygon_mesh_processing::bbox(tm_cc),
                 b2 = CGAL::Polygon_mesh_processing::bbox(clipper, np_c);
    if(CGAL::do_overlap(b1, b2))
      clip_single_cc(tm_cc, clipper, parameters::all_default(), np_c);
#endif

    CGAL::copy_face_graph(tm_cc, tm, parameters::all_default(), np_tm);
  }
}

template <typename TriangleMesh, typename Clipper,
          typename NamedParametersTM, typename NamedParametersC>
void generic_clip(TriangleMesh& tm,
                  Clipper& clipper,
                  const NamedParametersTM& np_tm,
                  const NamedParametersC& np_c)
{
  // @todo do the self-intersection test only once
  if(!does_self_intersect(tm, np_tm))
    clip(tm, clipper, np_tm, np_c); // @fixme some clip take a single np argument
  else
    clip_self_intersecting_mesh(tm, clipper, np_tm, np_c);
}

template <typename TriangleMesh, typename Clipper>
void generic_clip(TriangleMesh& tm,
                  Clipper& clipper)
{
  return generic_clip(tm, clipper, parameters::all_default(), parameters::all_default());
}


} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_CLIP_SELF_INTERSECTING_H
