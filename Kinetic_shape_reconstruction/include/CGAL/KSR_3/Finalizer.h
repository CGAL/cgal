// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot, Dmitry Anisimov

#ifndef CGAL_KSR_3_FINALIZER_H
#define CGAL_KSR_3_FINALIZER_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR/parameters.h>
#include <CGAL/KSR/conversions.h>

#include <CGAL/KSR_3/Data_structure.h>

namespace CGAL {
namespace KSR_3 {

template<typename GeomTraits>
class Finalizer {

public:
  using Kernel = GeomTraits;

private:
  using FT          = typename Kernel::FT;
  using Point_2     = typename Kernel::Point_2;
  using Point_3     = typename Kernel::Point_3;
  using Vector_2    = typename Kernel::Vector_2;
  using Vector_3    = typename Kernel::Vector_3;
  using Segment_3   = typename Kernel::Segment_3;
  using Line_3      = typename Kernel::Line_3;
  using Plane_3     = typename Kernel::Plane_3;
  using Direction_2 = typename Kernel::Direction_2;

  using Data_structure = KSR_3::Data_structure<Kernel>;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  using PVertex = typename Data_structure::PVertex;
  using PEdge   = typename Data_structure::PEdge;
  using PFace   = typename Data_structure::PFace;

  using Support_plane      = typename Data_structure::Support_plane;
  using Intersection_graph = typename Data_structure::Intersection_graph;
  using Volume_cell        = typename Data_structure::Volume_cell;

  using Mesh           = typename Data_structure::Mesh;
  using Vertex_index   = typename Data_structure::Vertex_index;
  using Face_index     = typename Data_structure::Face_index;
  using Edge_index     = typename Data_structure::Edge_index;
  using Halfedge_index = typename Data_structure::Halfedge_index;

  struct Vertex_info {
    bool tagged;
    PVertex pvertex;
    IVertex ivertex;
    Vertex_info() :
    tagged(false),
    pvertex(Data_structure::null_pvertex()),
    ivertex(Data_structure::null_ivertex())
    { }
  };

  struct Face_info {
    std::size_t index;
    std::size_t input;
    Face_info() :
    index(KSR::uninitialized()),
    input(KSR::uninitialized())
    { }
  };

  using VBI = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Kernel>;
  using FBI = CGAL::Triangulation_face_base_with_info_2<Face_info, Kernel>;
  using CFB = CGAL::Constrained_triangulation_face_base_2<Kernel, FBI>;
  using TDS = CGAL::Triangulation_data_structure_2<VBI, CFB>;
  using TAG = CGAL::Exact_intersections_tag;
  using EDT = CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, TAG>;
  using CDT = CGAL::Constrained_triangulation_plus_2<EDT>;
  using CID = typename CDT::Constraint_id;

  using Vertex_handle = typename CDT::Vertex_handle;
  using Face_handle   = typename CDT::Face_handle;
  using Edge          = typename CDT::Edge;

  using Parameters     = KSR::Parameters_3<FT>;
  using Kinetic_traits = KSR::Kinetic_traits_3<Kernel>;

public:
  Finalizer(Data_structure& data, const Parameters& parameters) :
  m_data(data), m_parameters(parameters), m_kinetic_traits(parameters.use_hybrid_mode)
  { }

  void clean() {

    std::size_t stop_value = 1;
    const bool should_be_removed = false;
    std::size_t num_hanging_pfaces = detect_hanging_pfaces(should_be_removed);

    if (num_hanging_pfaces >= stop_value) {
      if (m_parameters.verbose) {
        std::cout << "* number of hanging pfaces: " << num_hanging_pfaces << std::endl;
      }
      if (should_be_removed) return;
      const std::size_t num_added_pfaces = fill_holes(should_be_removed);
      CGAL_assertion(num_added_pfaces > 0);
      if (m_parameters.verbose) {
        std::cout << "* number of added pfaces: " << num_added_pfaces << std::endl;
      }
      num_hanging_pfaces = detect_hanging_pfaces(should_be_removed);
    }
    CGAL_assertion_msg(num_hanging_pfaces < stop_value,
      "ERROR: DO WE STILL HAVE HANGING PFACES?");
  }

  void create_polyhedra() {

    std::cout.precision(20);
    // for (std::size_t i = 0; i < number_of_support_planes(); ++i)
    //   std::cout << "num pfaces sp " << i << ": " << pfaces(i).size() << std::endl;

    CGAL_assertion(m_data.check_bbox());
    CGAL_assertion(m_data.check_interior());
    CGAL_assertion(m_data.check_vertices());
    CGAL_assertion(m_data.check_edges());
    create_volumes();
    CGAL_assertion(m_data.check_faces());
  }

  void clear() {
    // to be done
  }

private:
  Data_structure& m_data;
  const Parameters& m_parameters;
  Kinetic_traits m_kinetic_traits;

  /*******************************
  **          CLEANING          **
  ********************************/

  std::size_t detect_hanging_pfaces(const bool should_be_removed) {

    bool quit = true;
    std::size_t num_removed_pfaces = 0;
    do {
      quit = true;
      const auto iedges = m_data.iedges();
      for (const auto iedge : iedges) {
        const std::size_t num_pfaces =
          initialize_pface_removal(iedge, should_be_removed);
        if (num_pfaces != 0) {
          num_removed_pfaces += num_pfaces;
          if (should_be_removed) {
            quit = false; break;
          }
        }
      }
    } while (!quit);
    return num_removed_pfaces;
  }

  std::size_t initialize_pface_removal(
    const IEdge& iedge, const bool should_be_removed) {

    std::vector<PFace> pfaces;
    std::size_t num_removed_pfaces = 0;
    m_data.incident_faces(iedge, pfaces);
    if (pfaces.size() == 1) {
      return remove_pfaces(iedge, pfaces[0], should_be_removed);
    }
    if (pfaces.size() == 2) {
      const auto& pface0 = pfaces[0];
      const auto& pface1 = pfaces[1];
      if (pface0.first >= 6 && pface1.first >= 6 && pface0.first != pface1.first) {
        return remove_pfaces(iedge, pface0, should_be_removed);
      }
    }
    return num_removed_pfaces;
  }

  std::size_t remove_pfaces(
    const IEdge& init_iedge, const PFace& init_pface,
    const bool should_be_removed, const bool debug = false) {

    if (!should_be_removed) {
      if (debug) dump_pface(m_data, init_pface, "hang-" + m_data.str(init_pface));
      return 1; // use this to just count!
    }

    std::set<PFace> unique;
    std::vector< std::pair<Halfedge_index, PFace> > nfaces;
    const Halfedge_index init_he = find_crossing_he(init_iedge, init_pface);
    collect_connected_pfaces(init_he, init_pface, unique, nfaces);

    if (debug) {
      dump_pface(m_data, init_pface, "hang-" + m_data.str(init_pface));
      std::cout << "* found faces to remove: " << nfaces.size() << std::endl;
    }

    std::size_t num_removed_pfaces = 0;
    for (const auto& item : nfaces) {
      const auto& he = item.first;
      const auto& nface = item.second;
      const bool success = remove_pface(he, nface);
      if (success) ++num_removed_pfaces;
    }
    CGAL_assertion(num_removed_pfaces == nfaces.size());
    return num_removed_pfaces;
  }

  const Halfedge_index find_crossing_he(
    const IEdge& iedge, const PFace& pface) {

    const auto& mesh = m_data.mesh(pface.first);
    const auto pedges = m_data.pedges_of_pface(pface);
    bool found_pedge = false;
    for (const auto pedge : pedges) {
      CGAL_assertion(m_data.has_iedge(pedge));
      if (m_data.iedge(pedge) == iedge) {
        found_pedge = true;

        const auto he = mesh.halfedge(pedge.second);
        const auto op = mesh.opposite(he);
        const auto face1 = mesh.face(he);
        const auto face2 = mesh.face(op);
        const bool has_face1 = (face1 != Support_plane::Mesh::null_face());
        const bool has_face2 = (face2 != Support_plane::Mesh::null_face());
        if (!has_face1) {
          return op;
        } else if (!has_face2) {
          return he;
        } else {
          CGAL_assertion_msg(false, "ERROR: CROSSING HE IS NOT FOUND!");
        }
      }
    }
    CGAL_assertion(found_pedge);
    return Halfedge_index();
  }

  void collect_connected_pfaces(
    const Halfedge_index crossing_he,
    const PFace& pface,
    std::set<PFace>& unique,
    std::vector< std::pair<Halfedge_index, PFace> >& nfaces) {

    const auto pair = unique.insert(pface);
    if (!pair.second) return;

    CGAL_assertion(crossing_he != Halfedge_index());
    CGAL_assertion(pface != m_data.null_pface());
    CGAL_assertion(pface.second != Support_plane::Mesh::null_face());
    nfaces.push_back(std::make_pair(crossing_he, pface));

    const auto& mesh = m_data.mesh(pface.first);
    const auto pedges = m_data.pedges_of_pface(pface);
    for (const auto pedge : pedges) {
      CGAL_assertion(m_data.has_iedge(pedge));

      const PVertex pvertex(pface.first, 0);
      bool is_occupied_edge, bbox_reached;
      // std::tie(is_occupied_edge, bbox_reached) = m_data.is_occupied(pvertex, ivertex, m_data.iedge(pedge));
      std::tie(is_occupied_edge, bbox_reached) = m_data.is_occupied(pvertex, m_data.iedge(pedge));
      if (is_occupied_edge || bbox_reached) continue;

      const auto he = mesh.halfedge(pedge.second);
      const auto op = mesh.opposite(he);
      const auto face1 = mesh.face(he);
      const auto face2 = mesh.face(op);

      const auto nface1 = PFace(pface.first, face1);
      const auto nface2 = PFace(pface.first, face2);
      const bool has_nface1 = (face1 != Support_plane::Mesh::null_face());
      const bool has_nface2 = (face2 != Support_plane::Mesh::null_face());

      if (nface1 == pface) {
        if (has_nface2) {
          // std::cout << "adding nface2" << std::endl;
          collect_connected_pfaces(op, nface2, unique, nfaces);
        }
        continue;
      }
      if (nface2 == pface) {
        if (has_nface1) {
          // std::cout << "adding nface1" << std::endl;
          collect_connected_pfaces(he, nface1, unique, nfaces);
        }
        continue;
      }
      CGAL_assertion_msg(false, "ERROR: NO PFACE FOUND!");
    }
  }

  bool remove_pface(const Halfedge_index he, const PFace& pface) {

    const std::string plane_idx = std::to_string(pface.first);
    const std::string face_idx  = std::to_string(pface.second);

    // std::cout << "removing " << m_data.str(pface) << std::endl;
    // dump_pface(m_data, pface, "removed-pface-" + plane_idx + "-" + face_idx);

    auto& mesh = m_data.mesh(pface.first);
    CGAL::Euler::remove_face(he, mesh);
    return true;
  }

  std::size_t fill_holes(const bool already_removed) {

    // TODO: REIMPLEMENT IN A BETTER WAY:
    // First, sort all hanging pfaces by the number of potentially added pfaces;
    // then, start from the one that has the minimum such pfaces;
    // then, check again, because, after this insertion, other pfaces can be
    // reclassified into normal pfaces and there is no need to handle them.
    // I should also precompute CDT since I may have separate holes in the same plane.
    // See real-data-test -> test-15-polygons for k = 1.
    // If the hanging face is alone that is all its edges, which do not hang, are occupied,
    // I think it is better to remove it instead of adding a lot of new pfaces. Otherwise,
    // one small pface may lead to a lot of new pfaces. Should I do the same for two pfaces as well?

    bool quit = true;
    std::size_t num_added_pfaces = 0;
    CGAL_assertion(!already_removed);
    if (already_removed) return num_added_pfaces;

    do {
      quit = true;
      const auto iedges = m_data.iedges();
      for (const auto iedge : iedges) {
        const std::size_t num_pfaces =
          initialize_pface_insertion(iedge);
        if (num_pfaces != 0) {
          num_added_pfaces += num_pfaces;
          quit = false; break;
        }
      }
    } while (!quit);
    return num_added_pfaces;
  }

  std::size_t initialize_pface_insertion(const IEdge& iedge) {

    std::vector<PFace> pfaces;
    m_data.incident_faces(iedge, pfaces);
    if (pfaces.size() == 1) {
      // std::cout << "- hang iedge: " << m_data.segment_3(iedge) << std::endl;
      // std::cout << "- working out hanging: " << m_data.str(pfaces[0]) << std::endl;
      // dump_pface(m_data, pfaces[0], "hang-" + m_data.str(pfaces[0]));
      // CGAL_assertion_msg(false,
      // "TODO: IMPLEMENT CASE WITH ONE HANGING PFACE!");
      return create_pfaces(iedge, pfaces[0]);
    }

    std::size_t num_added_pfaces = 0;
    if (pfaces.size() == 2) {
      const auto& pface0 = pfaces[0];
      const auto& pface1 = pfaces[1];
      if (pface0.first >= 6 && pface1.first >= 6 && pface0.first != pface1.first) {
        std::cout << "- hang iedge: " << m_data.segment_3(iedge) << std::endl;
        std::cout << "- working out hanging: " << m_data.str(pface0) << "/" << m_data.str(pface1) << std::endl;
        dump_pface(m_data, pface0, "hang0-" + m_data.str(pface0));
        dump_pface(m_data, pface1, "hang1-" + m_data.str(pface1));
        CGAL_assertion_msg(false,
        "TODO: CAN WE HAVE TWO HANGING PFACES FROM DIFFERENT PLANES?");
        // num_added_pfaces += create_pfaces(iedge, pface0);
        // num_added_pfaces += create_pfaces(iedge, pface1);
      }
    }
    return num_added_pfaces;
  }

  std::size_t create_pfaces(
    const IEdge& init_iedge, const PFace& init_pface, const bool debug = false) {

    CDT cdt;
    std::map<CID, IEdge> map_intersections;
    const std::size_t support_plane_idx = init_pface.first;
    initialize_cdt(support_plane_idx, cdt, map_intersections);

    if (debug) {
      dump_2d_surface_mesh(m_data, support_plane_idx,
        "iter-10000-surface-mesh-before-" + std::to_string(support_plane_idx));
      dump_cdt(m_data, support_plane_idx, cdt,
        "/Users/monet/Documents/gf/kinetic/logs/volumes/initial-");
    }

    // CGAL_assertion_msg(false, "TODO: DEBUG THIS PART!");

    tag_cdt_exterior_faces(support_plane_idx, cdt, map_intersections);
    if (debug) {
      dump_cdt(m_data, support_plane_idx, cdt,
        "/Users/monet/Documents/gf/kinetic/logs/volumes/exterior-");
    }

    // CGAL_assertion_msg(false, "TODO: DEBUG THIS PART!");

    const auto num_original_pfaces = tag_cdt_interior_faces(cdt);
    if (debug) {
      std::cout << "- num original pfaces: " << num_original_pfaces << std::endl;
      dump_cdt(m_data, support_plane_idx, cdt,
        "/Users/monet/Documents/gf/kinetic/logs/volumes/interior-");
    }

    // CGAL_assertion_msg(false, "TODO: DEBUG THIS PART!");

    const Face_handle init_fh = find_initial_face(cdt, init_iedge);
    CGAL_assertion(init_fh != Face_handle());
    const auto num_detected_pfaces = tag_cdt_potential_faces(
      support_plane_idx, cdt, init_fh, num_original_pfaces);

    if (debug) {
      std::cout << "- num detected pfaces: " << num_detected_pfaces << std::endl;
      dump_cdt(m_data, support_plane_idx, cdt,
        "/Users/monet/Documents/gf/kinetic/logs/volumes/potential-");
    }

    // CGAL_assertion_msg(false, "TODO: DEBUG THIS PART!");

    const auto num_created_pfaces = insert_pfaces(support_plane_idx, cdt);
    if (debug) {
      std::cout << "- num created pfaces: " << num_created_pfaces << std::endl;
      dump_2d_surface_mesh(m_data, support_plane_idx,
        "iter-10000-surface-mesh-after-" + std::to_string(support_plane_idx));
    }

    CGAL_assertion(num_created_pfaces == num_detected_pfaces);
    reconnect_pvertices_to_ivertices(cdt);
    reconnect_pedges_to_iedges(cdt, map_intersections);

    // CGAL_assertion_msg(false, "TODO: CREATE MISSING PFACES!");
    return num_created_pfaces;
  }

  void initialize_cdt(
    const std::size_t sp_idx, CDT& cdt,
    std::map<CID, IEdge>& map_intersections) const {

    // Create unique ivertices.
    std::set<IVertex> ivertices;
    const auto& iedges = m_data.iedges(sp_idx);
    for (const auto& iedge : iedges) {
      ivertices.insert(m_data.source(iedge));
      ivertices.insert(m_data.target(iedge));
    }
    CGAL_assertion(ivertices.size() > 0);

    // Insert ivertices.
    std::map<IVertex, Vertex_handle> vhs_map;
    for (const auto& ivertex : ivertices) {
      const auto point = m_data.to_2d(sp_idx, ivertex);
      const auto vh = cdt.insert(point);
      vh->info().ivertex = ivertex;
      vhs_map[ivertex] = vh;
    }

    // Connect pvertices to ivertices.
    const auto all_pvertices = m_data.pvertices(sp_idx);
    for (const auto pvertex : all_pvertices) {
      CGAL_assertion(m_data.has_ivertex(pvertex));
      const auto ivertex = m_data.ivertex(pvertex);
      CGAL_assertion(vhs_map.find(ivertex) != vhs_map.end());
      const auto& vh = vhs_map.at(ivertex);
      vh->info().pvertex = pvertex;
    }

    // Insert iedges.
    for (const auto& iedge : iedges) {
      const auto isource = m_data.source(iedge);
      const auto itarget = m_data.target(iedge);
      CGAL_assertion(vhs_map.find(isource) != vhs_map.end());
      CGAL_assertion(vhs_map.find(itarget) != vhs_map.end());

      const auto& vh_source = vhs_map.at(isource);
      const auto& vh_target = vhs_map.at(itarget);
      const auto cid = cdt.insert_constraint(vh_source, vh_target);
      map_intersections.insert(std::make_pair(cid, iedge));
    }
  }

  void tag_cdt_exterior_faces(
    const std::size_t sp_idx, const CDT& cdt,
    const std::map<CID, IEdge>& map_intersections) const {

    std::queue<Face_handle> todo;
    todo.push(cdt.incident_faces(cdt.infinite_vertex()));
    while (!todo.empty()) {
      const auto fh = todo.front();
      todo.pop();
      if (fh->info().index != KSR::uninitialized()) {
        continue;
      }
      fh->info().index = KSR::no_element();

      for (std::size_t i = 0; i < 3; ++i) {
        const auto next = fh->neighbor(i);
        const auto edge = std::make_pair(fh, i);
        const bool is_border_edge =
          is_border(sp_idx, cdt, edge, map_intersections);
        if (!is_border_edge) {
          todo.push(next);
        }
      }
    }
    CGAL_assertion(todo.size() == 0);
  }

  bool is_border(
    const std::size_t sp_idx, const CDT& cdt, const Edge& edge,
    const std::map<CID, IEdge>& map_intersections) const {

    if (!cdt.is_constrained(edge))
      return false;

    const auto fh = edge.first;
    const auto id = edge.second;
    const std::size_t im = (id + 1) % 3;
    const std::size_t ip = (id + 2) % 3;

    const auto vh1 = fh->vertex(im);
    const auto vh2 = fh->vertex(ip);

    const auto ctx_begin = cdt.contexts_begin(vh1, vh2);
    const auto ctx_end   = cdt.contexts_end(vh1, vh2);

    for (auto cit = ctx_begin; cit != ctx_end; ++cit) {
      const auto iter = map_intersections.find(cit->id());
      if (iter == map_intersections.end()) continue;
      const auto& iedge = iter->second;
      CGAL_assertion(iedge != m_data.null_iedge());
      if (m_data.has_pedge(sp_idx, iedge)) return true;
    }
    return false;
  }

  std::size_t tag_cdt_interior_faces(const CDT& cdt) const {

    std::size_t face_index = 0;
    std::queue<Face_handle> todo;
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
      CGAL_assertion(todo.size() == 0);
      if (fit->info().index != KSR::uninitialized()) {
        continue;
      }

      todo.push(fit);
      std::size_t num_faces = 0;
      while (!todo.empty()) {
        const auto fh = todo.front();
        todo.pop();
        if (fh->info().index != KSR::uninitialized()) {
          continue;
        }
        fh->info().index = face_index;
        ++num_faces;

        for (std::size_t i = 0; i < 3; ++i) {
          const auto next = fh->neighbor(i);
          const auto edge = std::make_pair(fh, i);
          const bool is_constrained_edge = cdt.is_constrained(edge);
          if (!is_constrained_edge) {
            todo.push(next);
          }
        }
      }
      ++face_index;
      CGAL_assertion(todo.size() == 0);
    }
    return face_index;
  }

  Face_handle find_initial_face(
    const CDT& cdt, const IEdge& init_iedge) const {

    CGAL_assertion(init_iedge != m_data.null_iedge());
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
      if (fit->info().index != KSR::no_element()) continue;

      for (std::size_t i = 0; i < 3; ++i) {
        const auto edge = std::make_pair(fit, i);
        const auto iedge = find_iedge(edge);
        if (iedge == m_data.null_iedge()) {
          CGAL_assertion(!cdt.is_constrained(edge));
          continue;
        }
        if (iedge == init_iedge) {
          return static_cast<Face_handle>(fit);
        }
      }
    }
    CGAL_assertion_msg(false, "ERROR: NO INITIAL FACE FOUND!");
    return Face_handle();
  }

  IEdge find_iedge(const Edge& edge) const {

    const auto& fh = edge.first;
    const auto& id = edge.second;

    const auto im = (id + 1) % 3;
    const auto ip = (id + 2) % 3;

    const auto& vh1 = fh->vertex(im);
    const auto& vh2 = fh->vertex(ip);

    const auto& iv1 = vh1->info().ivertex;
    const auto& iv2 = vh2->info().ivertex;

    // std::cout << "iv1: " << point_3(iv1) << std::endl;
    // std::cout << "iv2: " << point_3(iv2) << std::endl;
    CGAL_assertion(iv1 != m_data.null_ivertex()); // if cdt has extra vertices with no ivertex,
    CGAL_assertion(iv2 != m_data.null_ivertex()); // just comment out these assertions

    IEdge iedge = m_data.null_iedge();
    if (m_data.igraph().is_edge(iv1, iv2)) {
      iedge = m_data.igraph().edge(iv1, iv2);
    } else if (m_data.igraph().is_edge(iv2, iv1)) {
      iedge = m_data.igraph().edge(iv2, iv1);
    }
    return iedge;
  }

  std::size_t tag_cdt_potential_faces(
    const std::size_t sp_idx,
    const CDT& cdt,
    const Face_handle& init_fh,
    const std::size_t num_faces) const {

    CGAL_assertion(init_fh != Face_handle());
    CGAL_assertion(init_fh->info().index == KSR::no_element());
    if (init_fh == Face_handle()) return 0;

    std::size_t face_index = num_faces;
    std::queue<Face_handle> todo_ext, todo_int;

    todo_ext.push(init_fh);
    while (!todo_ext.empty()) {
      const auto first = todo_ext.front();
      todo_ext.pop();

      bool is_new_face_detected = false;
      CGAL_assertion(todo_int.size() == 0);
      todo_int.push(first);
      while (!todo_int.empty()) {
        const auto fh = todo_int.front();
        todo_int.pop();
        if (fh->info().index != KSR::no_element())    continue;
        if (fh->info().input != KSR::uninitialized()) continue;

        is_new_face_detected = true;
        fh->info().index = face_index;
        fh->info().input = face_index;
        for (std::size_t i = 0; i < 3; ++i) {
          const auto next = fh->neighbor(i);
          const auto edge = std::make_pair(fh, i);
          bool is_exterior = false, is_interior = false;
          std::tie(is_exterior, is_interior) = is_crossing(sp_idx, cdt, edge);
          if (is_exterior) {
            CGAL_assertion(!is_interior);
            todo_ext.push(next);
          }
          if (is_interior) {
            CGAL_assertion(!is_exterior);
            todo_int.push(next);
          }
        }
      }
      if (is_new_face_detected) ++face_index;
      CGAL_assertion(todo_int.size() == 0);
    }

    const auto num_detected_pfaces = face_index - num_faces;
    CGAL_assertion(num_detected_pfaces > 0);
    return num_detected_pfaces;
  }

  std::pair<bool, bool> is_crossing(
    const std::size_t sp_idx, const CDT& cdt, const Edge& edge) const {

    const auto& init_fh = edge.first;
    const auto& init_id = edge.second;
    const auto fh = init_fh->neighbor(init_id);

    if (fh->info().index != KSR::no_element())    return std::make_pair(false, false);
    if (fh->info().input != KSR::uninitialized()) return std::make_pair(false, false);

    const auto iedge = find_iedge(edge);
    if (iedge == m_data.null_iedge()) {
      CGAL_assertion(!cdt.is_constrained(edge));
      return std::make_pair(false, true);
    }

    auto pvertex = m_data.null_pvertex();
    pvertex.first = sp_idx;
    bool is_occupied_edge = false, is_bbox_reached = false;
    std::tie(is_occupied_edge, is_bbox_reached) = m_data.is_occupied(pvertex, iedge);
    if (is_occupied_edge || is_bbox_reached) return std::make_pair(false, false);
    return std::make_pair(true, false);
  }

  std::size_t insert_pfaces(
    const std::size_t sp_idx, const CDT& cdt) {

    std::set<std::size_t> done;
    std::size_t num_created_pfaces = 0;

    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
      CGAL_assertion(fit->info().index != KSR::uninitialized());
      if (fit->info().input == KSR::uninitialized()) { // skip all faces with no input
        continue;
      }

      // Search for a constrained edge.
      Edge edge;
      for (std::size_t i = 0; i < 3; ++i) {
        edge = std::make_pair(fit, i);
        if (cdt.is_constrained(edge)) {
          break;
        }
      }

      // Skip pure interior faces.
      if (!cdt.is_constrained(edge)) {
        continue;
      }

      // If face index is already a part of the set, skip.
      const auto fh = edge.first;
      if (!done.insert(fh->info().index).second) {
        continue;
      }

      // Start from the constrained edge and traverse all constrained edges / boundary
      // of the triangulation part that is tagged with the same face index.
      // While traversing, add all missing pvertices.
      auto curr = edge;
      std::vector<PVertex> new_pvertices;
      do {
        const auto curr_face = curr.first;
        const int idx = curr.second;

        const auto source = curr_face->vertex(cdt.ccw(idx));
        const auto target = curr_face->vertex(cdt.cw (idx));
        if (source->info().pvertex == m_data.null_pvertex()) {
          source->info().pvertex = m_data.add_pvertex(sp_idx, source->point());
        }
        source->info().tagged = true;
        new_pvertices.push_back(source->info().pvertex);

        // Search for the next constrained edge.
        auto next = std::make_pair(curr_face, cdt.ccw(idx));
        while (!cdt.is_constrained(next)) {

          const auto next_face = next.first->neighbor(next.second);
          // Should be the same original polygon.
          CGAL_assertion(next_face->info().index == edge.first->info().index);

          const int next_idx = cdt.ccw(next_face->index(next.first));
          next = std::make_pair(next_face, next_idx);
        }
        // Check wether next source == previous target.
        CGAL_assertion(next.first->vertex(cdt.ccw(next.second)) == target);
        curr = next;

      } while (curr != edge);
      CGAL_assertion(curr == edge);

      // Add a new pface.
      const auto pface = m_data.add_pface(new_pvertices);
      ++num_created_pfaces;
      CGAL_assertion(pface != PFace());
    }

    // CGAL_assertion_msg(false, "TODO: INSERT DETECTED PFACES!");
    return num_created_pfaces;
  }

  void reconnect_pvertices_to_ivertices(const CDT& cdt) {

    // Reconnect only those, which have already been connected.
    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
      if (!vit->info().tagged) continue;
      if (vit->info().pvertex != m_data.null_pvertex() &&
          vit->info().ivertex != m_data.null_ivertex()) {
        m_data.connect(vit->info().pvertex, vit->info().ivertex);
      }
    }
  }

  void reconnect_pedges_to_iedges(
    const CDT& cdt,
    const std::map<CID, IEdge>& map_intersections) {

    // Reconnect only those, which have already been connected.
    for (const auto& item : map_intersections) {
      const auto& cid   = item.first;
      const auto& iedge = item.second;

      if (iedge == m_data.null_iedge()) {
        continue;
      }
      CGAL_assertion(iedge != m_data.null_iedge());

      auto vit = cdt.vertices_in_constraint_begin(cid);
      while (true) {
        auto next = vit; ++next;
        if (next == cdt.vertices_in_constraint_end(cid)) { break; }
        const auto a = *vit;
        const auto b = *next;
        vit = next;

        if (
          a->info().pvertex == m_data.null_pvertex() ||
          b->info().pvertex == m_data.null_pvertex()) {
          continue;
        }

        if (!a->info().tagged || !b->info().tagged) {
          continue;
        }

        CGAL_assertion(a->info().pvertex != m_data.null_pvertex());
        CGAL_assertion(b->info().pvertex != m_data.null_pvertex());
        // std::cout << "a: " << point_3(a->info().pvertex) << std::endl;
        // std::cout << "b: " << point_3(b->info().pvertex) << std::endl;
        // std::cout << "e: " << segment_3(iedge) << std::endl;
        m_data.connect(a->info().pvertex, b->info().pvertex, iedge);
      }
    }
  }

  /*******************************
  **     EXTRACTING VOLUMES     **
  ********************************/

  void create_volumes() {

    // Initialize an empty volume map.
    auto& volumes = m_data.volumes();
    auto& map_volumes = m_data.pface_neighbors();
    auto& volume_level_map = m_data.volume_level_map();

    volumes.clear();
    std::map<int, Point_3> centroids;
    map_volumes.clear();
    for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      const auto pfaces = m_data.pfaces(i);
      for (const auto pface : pfaces)
        map_volumes[pface] = std::make_pair(-1, -1);
    }

    // First, traverse only boundary volumes.
    bool is_found_new_volume = false;
    std::size_t volume_size = 0;
    int num_volumes  = 0;
    int volume_index = 0;
    int volume_level = 0;
    for (std::size_t i = 0; i < 6; ++i) {
      const auto pfaces = m_data.pfaces(i);
      for (const auto pface : pfaces) {
        CGAL_assertion(pface.first < 6);
        std::tie(is_found_new_volume, volume_size) = traverse_boundary_volume(
          pface, volume_index, num_volumes, map_volumes, centroids);
        if (is_found_new_volume) {
          CGAL_assertion(m_data.check_volume(volume_index, volume_size, map_volumes));
          ++volume_index;
        }
      }
    }
    if (m_parameters.verbose) {
      std::cout << "* found boundary volumes: "<< volume_index << std::endl;
    }
    num_volumes = volume_index;
    CGAL_assertion(num_volumes > 0);
    volume_level_map[volume_level] = static_cast<std::size_t>(num_volumes);
    ++volume_level;

    // Then traverse all other volumes if any.
    std::vector<PFace> other_pfaces;
    for (std::size_t i = 6; i < m_data.number_of_support_planes(); ++i) {
      const auto pfaces = m_data.pfaces(i);
      for (const auto pface : pfaces) {
        CGAL_assertion(pface.first >= 6);
        other_pfaces.push_back(pface);
      }
    }

    std::sort(other_pfaces.begin(), other_pfaces.end(),
      [&](const PFace& pface1, const PFace& pface2) -> bool {
        const auto pedges1 = m_data.pedges_of_pface(pface1);
        const auto pedges2 = m_data.pedges_of_pface(pface2);
        return pedges1.size() > pedges2.size();
      }
    );

    bool quit = true;
    do {
      quit = true;
      const int before = volume_index;
      for (const auto& other_pface : other_pfaces) {
        std::tie(is_found_new_volume, volume_size) = traverse_interior_volume(
          other_pface, volume_index, num_volumes, map_volumes, centroids);
        if (is_found_new_volume) {
          quit = false;
          CGAL_assertion(m_data.check_volume(volume_index, volume_size, map_volumes));
          ++volume_index;
        }
      }
      const int after = volume_index;
      if (m_parameters.verbose) {
        std::cout << "* found interior volumes: "<< after - before << std::endl;
      }
      num_volumes = volume_index;
      CGAL_assertion(after >= before);
      if (after > before) {
        volume_level_map[volume_level] = static_cast<std::size_t>(after - before);
        ++volume_level;
      }

    } while (!quit);

    // Now, set final volumes and their neighbors.
    for (const auto& item : map_volumes) {
      const auto& pface = item.first;
      const auto& pair  = item.second;

      if (pair.first == -1) {
        dump_pface(m_data, pface, "face-debug");
        std::cout << "DEBUG face: " << m_data.str(pface) << " "   << std::endl;
        std::cout << "DEBUG  map: " << pair.first << " : " << pair.second << std::endl;
      }

      CGAL_assertion(pair.first != -1);
      if (volumes.size() <= static_cast<std::size_t>(pair.first))
        volumes.resize(pair.first + 1);
      volumes[pair.first].add_pface(pface, pair.second);

      if (pface.first < 6 && pair.second == -1) continue;
      CGAL_assertion(pair.second != -1);
      if (volumes.size() <= static_cast<std::size_t>(pair.second))
        volumes.resize(pair.second + 1);
      volumes[pair.second].add_pface(pface, pair.first);
    }
    for (auto& volume : volumes)
      create_cell_pvertices(volume);

    if (m_parameters.verbose) {
      std::cout << "* created volumes: " << volumes.size() << std::endl;
      if (m_parameters.export_all) dump_volumes(m_data, "final");
      for (std::size_t i = 0; i < volumes.size(); ++i) {
        const auto& volume = volumes[i];
        CGAL_assertion(volume.pfaces.size() > 3);
        if (m_parameters.debug) {
          std::cout <<
          " VOLUME "     << std::to_string(i) << ": "
          " pvertices: " << volume.pvertices.size() <<
          " pfaces: "    << volume.pfaces.size()    << std::endl;
        }
      }
    }

    CGAL_assertion(volumes.size() == centroids.size());
    for (std::size_t i = 0; i < volumes.size(); ++i) {
      auto& volume = volumes[i];
      volume.set_index(i);
      volume.set_centroid(centroids.at(i));
    }
  }

  const std::pair<bool, std::size_t> traverse_boundary_volume(
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::map<int, Point_3>& centroids) const {

    CGAL_assertion(num_volumes  == 0);
    CGAL_assertion(volume_index >= 0);
    if (pface.first >= 6) return std::make_pair(false, 0);
    CGAL_assertion(pface.first < 6);
    const auto& pair = map_volumes.at(pface);
    CGAL_assertion(pair.second == -1);
    if (pair.first != -1) return std::make_pair(false, 0);
    CGAL_assertion(pair.first == -1);

    std::deque<PFace> queue;
    queue.push_front(pface);

    Point_3 volume_centroid;
    std::size_t volume_size = 0;

    while (!queue.empty()) {
      // print_queue(volume_index, queue);
      const auto query = queue.front();
      queue.pop_front();
      propagate_pface(
        false, query, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    }

    if (m_parameters.debug) {
      std::cout << "- FOUND VOLUME " << volume_index << ", (SIZE/BARYCENTER): "
      << volume_size << " / " << volume_centroid << std::endl;
    }
    centroids[volume_index] = volume_centroid;
    return std::make_pair(true, volume_size);
  }

  const std::pair<bool, std::size_t> traverse_interior_volume(
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::map<int, Point_3>& centroids) const {

    CGAL_assertion(volume_index > 0);
    CGAL_assertion(volume_index >= num_volumes);

    if (pface.first < 6) return std::make_pair(false, 0);
    CGAL_assertion(pface.first >= 6);
    const auto& pair = map_volumes.at(pface);
    if (pair.second != -1) {
      CGAL_assertion(pair.first != -1);
      return std::make_pair(false, 0);
    }
    CGAL_assertion(pair.second == -1);
    if (pair.first == -1) {
      CGAL_assertion(pair.second == -1);
      return std::make_pair(false, 0);
    }
    CGAL_assertion(pair.first != -1);
    if (pair.first >= num_volumes) return std::make_pair(false, 0);
    CGAL_assertion(pair.first < num_volumes);

    std::deque<PFace> queue;
    queue.push_front(pface);

    Point_3 volume_centroid;
    std::size_t volume_size = 0;

    while (!queue.empty()) {
      // print_queue(volume_index, queue);
      const auto query = queue.front();
      queue.pop_front();
      propagate_pface(
        false, query, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    }
    if (m_parameters.debug) {
      std::cout << "- FOUND VOLUME " << volume_index << ", (SIZE/BARYCENTER): "
      << volume_size << " / " << volume_centroid << std::endl;
    }
    centroids[volume_index] = volume_centroid;
    return std::make_pair(true, volume_size);
  }

  void print_queue(
    const int volume_index,
    const std::deque<PFace>& queue) const {

    // if (volume_index != -1) return;
    std::cout << "QUEUE: " << std::endl;
    for (const auto& pface : queue) {
      std::cout << volume_index << " "
      << pface.first << " " << pface.second << std::endl;
    }
  }

  void propagate_pface(
    const bool verbose,
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<int, Point_3>& centroids,
    std::size_t& volume_size,
    Point_3& volume_centroid,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::deque<PFace>& queue) const {

    const bool is_boundary = is_boundary_pface(
      pface, volume_index, num_volumes, map_volumes);
    if (is_boundary) {
      propagate_boundary_pface(
        verbose, pface, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    } else {
      propagate_interior_pface(
        verbose, pface, volume_index, num_volumes, centroids,
        volume_size, volume_centroid, map_volumes, queue);
    }
  }

  bool is_boundary_pface(
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<PFace, std::pair<int, int> >& map_volumes) const {

    CGAL_assertion(volume_index >= 0);
    if (pface.first < 6) return true;
    CGAL_assertion(pface.first >= 6);
    if (num_volumes == 0) return false;
    CGAL_assertion(num_volumes  > 0);
    CGAL_assertion(volume_index > 0);
    CGAL_assertion(volume_index >= num_volumes);

    const auto& pair = map_volumes.at(pface);
    if (pair.first == -1) {
      CGAL_assertion(pair.second == -1);
      return false;
    }
    CGAL_assertion(pair.first != -1);
    if (pair.first < num_volumes) return true;
    CGAL_assertion(pair.first >= num_volumes);
    return false;
  }

  void propagate_boundary_pface(
    const bool verbose,
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<int, Point_3>& centroids,
    std::size_t& volume_size,
    Point_3& volume_centroid,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::deque<PFace>& queue) const {

    auto& pair = map_volumes.at(pface);
    if (pair.first >= num_volumes) return;
    CGAL_assertion(pair.first < num_volumes);
    if (pair.second != -1) {
      CGAL_assertion(pair.first != -1);
      return;
    }
    CGAL_assertion(pair.second == -1);

    if (pair.first == -1) {
      pair.first  = volume_index;
    } else {
      pair.second = volume_index;
    }

    Point_3 centroid = m_data.centroid_of_pface(pface);
    if (num_volumes > 0) {
      // std::cout << "SHIFTING CENTROID" << std::endl;

      CGAL_assertion(pair.first < num_volumes);
      CGAL_assertion(centroids.find(pair.first) != centroids.end());
      const auto& other_centroid = centroids.at(pair.first);
      const auto plane = m_data.plane_of_pface(pface);
      auto vec1 = plane.orthogonal_vector();
      vec1 = KSR::normalize(vec1);
      auto vec2 = Vector_3(centroid, other_centroid);
      vec2 = KSR::normalize(vec2);

      // TODO: CAN WE AVOID THIS VALUE?
      const FT tol = KSR::tolerance<FT>();
      const FT dot_product = vec1 * vec2;

      if (dot_product < FT(0)) {
        centroid += tol * vec1;
      } else {
        centroid -= tol * vec1;
      }
      volume_centroid = CGAL::barycenter(
        volume_centroid, static_cast<FT>(volume_size), centroid, FT(1));

    } else {
      volume_centroid = CGAL::barycenter(
        volume_centroid, static_cast<FT>(volume_size), centroid, FT(1));
    }

    // std::cout << "volume centroid: " << volume_centroid << std::endl;
    ++volume_size;

    if (verbose) {
      // std::cout << "BND PFACE MAP: (" <<
      // pair.first << ", " << pair.second << ")" << std::endl;
      std::cout << "DUMPING BND PFACE: " <<
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second) << std::endl;
      dump_pface(m_data, pface, "bnd-pface-" +
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second));
    }

    std::vector<PFace> nfaces, bnd_nfaces, int_nfaces, all_nfaces;
    const auto pedges = m_data.pedges_of_pface(pface);
    for (const auto pedge : pedges) {
      CGAL_assertion(m_data.has_iedge(pedge));
      m_data.incident_faces(m_data.iedge(pedge), nfaces);
      split_pfaces(
        pface, volume_index, num_volumes, map_volumes, nfaces,
        bnd_nfaces, int_nfaces, all_nfaces);

      if (num_volumes == 0) {
        CGAL_assertion(bnd_nfaces.size() == 1);
        CGAL_assertion(int_nfaces.size() == 0 || int_nfaces.size() == 1);
      }

      if (int_nfaces.size() == 1) {
        queue.push_back(int_nfaces[0]);
        continue;
      }

      if (int_nfaces.size() == 0 && bnd_nfaces.size() == 1) {
        queue.push_front(bnd_nfaces[0]);
        continue;
      }

      if (all_nfaces.size() == 0) {
        dump_info(m_data, pface, pedge, nfaces);
        std::cout << "DEBUG: num nfaces: " << nfaces.size() << std::endl;
      }
      CGAL_assertion(all_nfaces.size() > 0);

      const auto found_nface = find_using_2d_directions(
      volume_index, volume_centroid, pface, pedge, all_nfaces);
      if (found_nface == m_data.null_pface()) continue;

      if (is_boundary_pface(
        found_nface, volume_index, num_volumes, map_volumes)) {
        queue.push_front(found_nface);
      } else {
        queue.push_back(found_nface);
      }
    }
  }

  void propagate_interior_pface(
    const bool verbose,
    const PFace& pface,
    const int volume_index,
    const int num_volumes,
    const std::map<int, Point_3>& /* centroids */,
    std::size_t& volume_size,
    Point_3& volume_centroid,
    std::map<PFace, std::pair<int, int> >& map_volumes,
    std::deque<PFace>& queue) const {

    CGAL_assertion(num_volumes >= 0);
    auto& pair = map_volumes.at(pface);
    if (pair.first != -1 && pair.second != -1) return;
    CGAL_assertion(pair.second == -1);
    if (pair.first == volume_index) return;
    CGAL_assertion(pair.first != volume_index);
    if (pair.first != -1) {
      pair.second = volume_index;
    } else {
      pair.first  = volume_index;
    }

    if (verbose) {
      std::cout << "pface: " << m_data.str(pface) << std::endl;
      std::cout << "pair: " <<
      std::to_string(pair.first) << "/" << std::to_string(pair.second) << std::endl;
    }

    const Point_3 centroid = m_data.centroid_of_pface(pface);
    volume_centroid = CGAL::barycenter(
      volume_centroid, static_cast<FT>(volume_size), centroid, FT(1));
    // std::cout << "volume centroid: " << volume_centroid << std::endl;
    ++volume_size;

    if (verbose) {
      // std::cout << "INT PFACE MAP: (" <<
      // pair.first << ", " << pair.second << ")" << std::endl;
      std::cout << "DUMPING INT PFACE: " <<
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second) << std::endl;
      dump_pface(m_data, pface, "int-pface-" +
        std::to_string(volume_index) + "-" +
        std::to_string(pface.first) + "-" +
        std::to_string(pface.second));
    }

    std::vector<PFace> nfaces, bnd_nfaces, int_nfaces, all_nfaces;
    const auto pedges = m_data.pedges_of_pface(pface);
    for (const auto pedge : pedges) {
      CGAL_assertion(m_data.has_iedge(pedge));
      m_data.incident_faces(m_data.iedge(pedge), nfaces);
      split_pfaces(
        pface, volume_index, num_volumes, map_volumes, nfaces,
        bnd_nfaces, int_nfaces, all_nfaces);

      if (all_nfaces.size() == 0) {
        dump_info(m_data, pface, pedge, nfaces);
        std::cout << "DEBUG: num nfaces: " << nfaces.size() << std::endl;
      }
      CGAL_assertion(all_nfaces.size() > 0);

      const auto found_nface = find_using_2d_directions(
      volume_index, volume_centroid, pface, pedge, all_nfaces);
      if (found_nface == m_data.null_pface()) continue;

      if (is_boundary_pface(
        found_nface, volume_index, num_volumes, map_volumes)) {
        queue.push_front(found_nface);
      } else {
        queue.push_back(found_nface);
      }
    }
  }

  void split_pfaces(
    const PFace& current,
    const int volume_index,
    const int num_volumes,
    const std::map<PFace, std::pair<int, int> >& map_volumes,
    const std::vector<PFace>& pfaces,
    std::vector<PFace>& bnd_pfaces,
    std::vector<PFace>& int_pfaces,
    std::vector<PFace>& all_pfaces) const {

    bnd_pfaces.clear();
    int_pfaces.clear();
    all_pfaces.clear();
    for (const auto& pface : pfaces) {
      if (pface == current) continue;
      CGAL_assertion(pface != current);
      all_pfaces.push_back(pface);

      const auto& pair = map_volumes.at(pface);
      if (num_volumes > 0 && pair.first != -1) {
        if (pair.first < num_volumes && pair.second != -1) {
          if (pair.second < num_volumes) {
            continue;
          }
          CGAL_assertion(pair.second >= num_volumes);
        }
      }
      if (is_boundary_pface(
        pface, volume_index, num_volumes, map_volumes))  {
        bnd_pfaces.push_back(pface);
      } else {
        int_pfaces.push_back(pface);
      }
    }
  }

  const PFace find_using_2d_directions(
    const int /* volume_index */,
    const Point_3& volume_centroid,
    const PFace& pface,
    const PEdge& pedge,
    const std::vector<PFace>& nfaces) const {

    CGAL_assertion(nfaces.size() > 0);
    if (nfaces.size() == 1) return nfaces[0];
    const bool debug = false;
      // ( volume_index == 31 &&
      //   pface.first == 8 &&
      //   static_cast<std::size_t>(pface.second) == 7);

    if (debug) {
      dump_info(m_data, pface, pedge, nfaces);
    }
    CGAL_assertion(nfaces.size() > 1);

    Point_3 center = m_data.centroid_of_pface(pface);
    const Segment_3 segment = m_data.segment_3(pedge);
    const Line_3 line(segment.source(), segment.target());
    Point_3 midp = CGAL::midpoint(segment.source(), segment.target());
    // std::cout << "midp: " << midp << std::endl;
    Vector_3 norm(segment.source(), segment.target());
    norm = KSR::normalize(norm);
    const Plane_3 plane(midp, norm);

    std::vector<Point_3> points;
    points.reserve(nfaces.size() + 2);

    points.push_back(midp);
    points.push_back(center);
    for (const auto& nface : nfaces) {
      center = m_data.centroid_of_pface(nface);
      points.push_back(center);
    }
    CGAL_assertion(points.size() >= 3);

    for (auto& point : points) {
      point = plane.projection(point);
    }

    if (debug) {
      dump_frame(points, "volumes/directions-init");
    }

    const FT cx = volume_centroid.x();
    const FT cy = volume_centroid.y();
    const FT cz = volume_centroid.z();
    FT d = (norm.x() * cx + norm.y() * cy + norm.z() * cz + plane.d());

    // std::cout << "1 d: " << d << std::endl;
    // std::cout << "1 norm: " << norm << std::endl;
    const Plane_3 tr_plane(midp + norm * d, norm);
    Point_3 inter;
    const bool is_intersection_found =
      m_kinetic_traits.intersection(line, tr_plane, inter);
    if (!is_intersection_found) {
      std::cout << "d = " << d << std::endl;
    }
    CGAL_assertion(is_intersection_found);
    // std::cout << "inter: " << inter << std::endl;

    d = KSR::distance(midp, inter);
    norm = Vector_3(midp, inter);
    // std::cout << "2 d: " << d << std::endl;
    // std::cout << "2 norm: " << norm << std::endl;

    if (d != FT(0)) {
      CGAL_assertion(norm != Vector_3(FT(0), FT(0), FT(0)));
      norm = KSR::normalize(norm);
      for (auto& point : points) {
        point += norm * d;
      }
    }

    if (debug) {
      auto extended = points;
      extended.push_back(volume_centroid);
      dump_frame(extended, "volumes/directions");
    }

    std::vector< std::pair<Direction_2, PFace> > dir_edges;
    dir_edges.reserve(nfaces.size() + 1);

    const Point_2 proj_0 = plane.to_2d(points[0]);
    for (std::size_t i = 1; i < points.size(); ++i) {
      const Point_2 proj_i = plane.to_2d(points[i]);
      const Vector_2 vec(proj_0, proj_i);
      if (i == 1) {
        dir_edges.push_back(std::make_pair(Direction_2(vec), pface));
      } else {
        dir_edges.push_back(std::make_pair(Direction_2(vec), nfaces[i - 2]));
      }
    }
    CGAL_assertion(dir_edges.size() == nfaces.size() + 1);

    const Point_2 proj_vc = plane.to_2d(volume_centroid);
    const Vector_2 vec(proj_0, proj_vc);
    const Direction_2 ref_dir(vec);

    std::sort(dir_edges.begin(), dir_edges.end(), [&](
      const std::pair<Direction_2, PFace>& p,
      const std::pair<Direction_2, PFace>& q) -> bool {
        return p.first < q.first;
      }
    );

    const std::size_t n = dir_edges.size();
    for (std::size_t i = 0; i < n; ++i) {
      if (dir_edges[i].second == pface) {

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        const auto& dir_prev = dir_edges[im].first;
        const auto& dir_curr = dir_edges[i].first;
        const auto& dir_next = dir_edges[ip].first;

        if (debug) {
          dump_pface(m_data, dir_edges[im].second, "prev");
          dump_pface(m_data, dir_edges[ip].second, "next");
        }

        if (ref_dir.counterclockwise_in_between(dir_prev, dir_curr)) {
          if (debug) {
            std::cout << "found prev" << std::endl;
            exit(EXIT_SUCCESS);
          }
          return dir_edges[im].second;
        } else if (ref_dir.counterclockwise_in_between(dir_curr, dir_next)) {
          if (debug) {
            std::cout << "found next" << std::endl;
            exit(EXIT_SUCCESS);
          }
          return dir_edges[ip].second;
        } else {
          // return m_data.null_pface();
          std::cout << "ERROR: WRONG ORIENTATIONS!" << std::endl;
          dump_info(m_data, pface, pedge, nfaces);
          dump_frame(points, "volumes/directions-init");
          auto extended = points;
          extended.push_back(volume_centroid);
          dump_frame(extended, "volumes/directions");
          CGAL_assertion_msg(false, "ERROR: WRONG ORIENTATIONS!");
        }
      }
    }

    CGAL_assertion_msg(false, "ERROR: NEXT PFACE IS NOT FOUND!");
    return m_data.null_pface();
  }

  void create_cell_pvertices(Volume_cell& cell) {
    cell.pvertices.clear();
    for (const auto& pface : cell.pfaces) {
      for (const auto pvertex : m_data.pvertices_of_pface(pface)) {
        cell.pvertices.insert(pvertex);
      }
    }
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_FINALIZER_H
