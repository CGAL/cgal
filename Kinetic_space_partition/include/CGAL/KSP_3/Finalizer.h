// Copyright (c) 2023 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Florent Lafarge, Dmitry Anisimov, Simon Giraudot

#ifndef CGAL_KSP_3_FINALIZER_H
#define CGAL_KSP_3_FINALIZER_H

#include <CGAL/license/Kinetic_space_partition.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

// Internal includes.
#include <CGAL/KSP/utils.h>
#include <CGAL/KSP/debug.h>
#include <CGAL/KSP/parameters.h>

#include <CGAL/KSP_3/Data_structure.h>

namespace CGAL {
namespace KSP_3 {
namespace internal {

#ifdef DOXYGEN_RUNNING
#else

template<typename GeomTraits, typename IntersectionKernel>
class Finalizer {

public:
  using Kernel = GeomTraits;

private:
  using Intersection_kernel = IntersectionKernel;
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Vector_2 = typename Kernel::Vector_2;
  using Vector_3 = typename Kernel::Vector_3;
  using Segment_3 = typename Kernel::Segment_3;
  using Line_3 = typename Kernel::Line_3;
  using Plane_3 = typename Kernel::Plane_3;
  using Direction_2 = typename Kernel::Direction_2;
  using Tetrahedron_3 = typename Kernel::Tetrahedron_3;

  using From_exact = CGAL::Cartesian_converter<IntersectionKernel, Kernel>;

  using Data_structure = CGAL::KSP_3::internal::Data_structure<Kernel, IntersectionKernel>;

  using IVertex = typename Data_structure::IVertex;
  using IEdge = typename Data_structure::IEdge;

  using PVertex = typename Data_structure::PVertex;
  using PEdge = typename Data_structure::PEdge;
  using PFace = typename Data_structure::PFace;

  using Support_plane = typename Data_structure::Support_plane;
  using Intersection_graph = typename Data_structure::Intersection_graph;
  using Volume_cell = typename Data_structure::Volume_cell;

  using Mesh = typename Data_structure::Mesh;
  using Vertex_index = typename Data_structure::Vertex_index;
  using Face_index = typename Data_structure::Face_index;
  using Edge_index = typename Data_structure::Edge_index;
  using Halfedge_index = typename Data_structure::Halfedge_index;

  using F_component_map = typename Mesh::template Property_map<typename Support_plane::Face_index, typename boost::graph_traits<typename Support_plane::Mesh>::faces_size_type>;
  using E_constraint_map = typename Mesh::template Property_map<typename Support_plane::Edge_index, bool>;

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
      index(-1),
      input(-1)
    { }
  };

  using Parameters = KSP::internal::Parameters_3<FT>;

public:
  Finalizer(Data_structure& data, const Parameters& parameters) :
    m_data(data), m_parameters(parameters)
  { }

  void create_polyhedra() {
    if (m_parameters.debug) {
      for (std::size_t sp = 0; sp < m_data.number_of_support_planes(); sp++) {
        dump_2d_surface_mesh(m_data, sp, m_data.prefix() + "after-partition-sp" + std::to_string(sp));
      }
    }

    std::cout.precision(20);
    CGAL_assertion(m_data.check_bbox());
    CGAL_assertion(m_data.check_interior());
    CGAL_assertion(m_data.check_vertices());
    CGAL_assertion(m_data.check_edges());

    merge_facets_connected_components();

    create_volumes();

    CGAL_assertion(m_data.check_faces());
  }

private:
  Data_structure& m_data;
  const Parameters& m_parameters;

  /*******************************
  **     EXTRACTING VOLUMES     **
  ********************************/

  void calculate_centroid(Volume_cell& volume) {
    // First find a point in the interior of the volume cell.
    FT x = 0, y = 0, z = 0;
    for (const PVertex& v : volume.pvertices) {
      Point_3 p = m_data.point_3(v);
      x += p.x();
      y += p.y();
      z += p.z();
    }
    Point_3 inside(x / volume.pvertices.size(), y / volume.pvertices.size(), z / volume.pvertices.size());

    // Now create a vector of tetrahedrons.
    std::vector<Tetrahedron_3> tets;
    tets.reserve(volume.pvertices.size());

    for (const PFace& f : volume.pfaces) {
      // Orientation of the tetrahedron depends on the orientation of the support plane of the polygon.
      Support_plane sp = m_data.support_plane(f.first);
      Vector_3 n = sp.plane().orthogonal_vector();
      const Mesh& m = sp.mesh();
      Halfedge_index first = m.halfedge(f.second);
      Point_3 p = m_data.point_3(PVertex(f.first, m.target(first)));

      bool positive_side = (inside - p) * n > 0;

      Halfedge_index h = m.next(first);
      Point_3 a = m_data.point_3(PVertex(f.first, m.target(h)));
      h = m.next(h);
      do {
        Point_3 b = m_data.point_3(PVertex(f.first, m.target(h)));

        if (positive_side)
          tets.push_back(Tetrahedron_3(p, a, b, inside));
        else
          tets.push_back(Tetrahedron_3(p, b, a, inside));

        a = b;
        h = m.next(h);
      } while (h != first);
    }

    volume.centroid = CGAL::centroid(tets.begin(), tets.end(), CGAL::Dimension_tag<3>());
  }

  void create_volumes() {
    // Initialize an empty volume map.
    auto& volumes = m_data.volumes();//std::vector<Volume_cell>
    auto& map_volumes = m_data.pface_neighbors();//std::map<PFace, std::pair<int, int>

    volumes.clear();
    map_volumes.clear();
    for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      const auto pfaces = m_data.pfaces(i);
      for (const auto pface : pfaces)
        map_volumes[pface] = std::make_pair(-1, -1);
    }

    for (std::size_t sp = 0; sp < m_data.number_of_support_planes(); sp++) {
      for (auto pface : m_data.pfaces(sp)) {
        segment_adjacent_volumes(pface, volumes, map_volumes);
      }
    }

    std::map<PFace, std::size_t>& face2index = m_data.face_to_index();
    std::vector<std::pair<int, int> >& face2volumes = m_data.face_to_volumes();
    std::size_t num_faces = 0;

    // Adjust neighbor information in volumes
    for (std::size_t i = 0; i < volumes.size(); i++) {
      auto& v = volumes[i];
      v.index = i;
      v.faces.resize(v.pfaces.size());
      for (std::size_t f = 0; f < volumes[i].pfaces.size(); f++) {
        auto& pf = volumes[i].pfaces[f];
        auto it = face2index.find(pf);
        if (it == face2index.end()) {
          face2index[pf] = num_faces;
          v.faces[f] = num_faces++;
        }
        else
          v.faces[f] = it->second;
      }

      if (face2volumes.size() < num_faces)
        face2volumes.resize(num_faces, std::pair<int, int>(-1, -1));

      for (std::size_t j = 0; j < volumes[i].neighbors.size(); j++) {
        auto& pair = map_volumes.at(volumes[i].pfaces[j]);
        if (pair.second == -1)
          pair.second = -static_cast<int>(volumes[i].pfaces[j].first + 1);
        volumes[i].neighbors[j] = (pair.first == static_cast<int>(i)) ? pair.second : pair.first;
        face2volumes[v.faces[j]] = pair;
      }
    }

    m_data.face_is_part_of_input_polygon().resize(num_faces);
    for (const auto& p : face2index)
      m_data.face_is_part_of_input_polygon()[p.second] = m_data.support_plane(p.first.first).is_initial(p.first.second);

    m_data.face_to_vertices().resize(num_faces);
    m_data.face_to_support_plane().resize(num_faces);

    for (auto& volume : volumes) {
      create_cell_pvertices(volume);
      calculate_centroid(volume);
      volume.pface_oriented_outwards.resize(volume.pfaces.size());
      for (std::size_t i = 0; i < volume.pfaces.size(); i++) {
        PVertex vtx = *m_data.pvertices_of_pface(volume.pfaces[i]).begin();
        volume.pface_oriented_outwards[i] = ((m_data.point_3(vtx) - volume.centroid) * m_data.support_plane(volume.pfaces[i]).plane().orthogonal_vector() < 0);
      }
    }

    remove_collinear_vertices();
  }

  void segment_adjacent_volumes(const PFace& pface, std::vector<Volume_cell>& volumes, std::map<PFace, std::pair<int, int> >& map_volumes) {

    // Check whether this face is already part of one or two volumes
    auto& pair = map_volumes.at(pface);
    if (pair.first != -1 && pair.second != -1) return;

    std::size_t volume_indices[] = { static_cast<std::size_t>(-1), static_cast<std::size_t>(-1) };
    std::size_t other[] = { static_cast<std::size_t>(-1), static_cast<std::size_t>(-1) };

    // Start new volume cell
    // First of pair is positive side, second is negative
    if (pair.first == -1) {
      volume_indices[0] = volumes.size();
      pair.first = static_cast<int>(volumes.size());
      volumes.push_back(Volume_cell());
      volumes.back().add_pface(pface, pair.second);
    }
    else {
      other[0] = pair.first;
      if (pface.first < 6)
        // Thus for a bbox pair.second is always -1. Thus if pair.first is already set, there is nothing to do.
        return;
    }

    if (pair.second == -1 && pface.first >= 6) {
      volume_indices[1] = volumes.size();
      pair.second = static_cast<int>(volumes.size());
      volumes.push_back(Volume_cell());
      volumes.back().add_pface(pface, pair.first);
    }
    else other[1] = pair.second;

    // 0 is queue on positive side, 1 is queue on negative side
    std::queue<std::pair<PFace, Oriented_side> > queue[2];

    // Get neighbors with smallest dihedral angle
    const auto pedges = m_data.pedges_of_pface(pface);
    std::vector<PFace> neighbor_faces, adjacent_faces;
    for (const auto pedge : pedges) {
      CGAL_assertion(m_data.has_iedge(pedge));

      m_data.incident_faces(m_data.iedge(pedge), neighbor_faces);

      if (neighbor_faces.size() == 2) {
        // If there is only one neighbor, the edge is on the corner of the bbox.
        // Thus the only neighbor needs to be a bbox face.
        PFace neighbor = (neighbor_faces[0] == pface) ? neighbor_faces[1] : neighbor_faces[0];
        CGAL_assertion(neighbor.first < 6 && pface.first < 6);
        Oriented_side side = oriented_side(pface, neighbor);
        Oriented_side inverse_side = oriented_side(neighbor, pface);
        CGAL_assertion(side != COPLANAR && inverse_side != COPLANAR);

        if (side == ON_POSITIVE_SIDE && volume_indices[0] != static_cast<std::size_t>(-1)) {
          if (associate(neighbor, volume_indices[0], inverse_side, volumes, map_volumes))
            queue[0].push(std::make_pair(neighbor, inverse_side));
        }
        else if (side == ON_NEGATIVE_SIDE && volume_indices[1] != static_cast<std::size_t>(-1))
          if (associate(neighbor, volume_indices[1], inverse_side, volumes, map_volumes))
            queue[1].push(std::make_pair(neighbor, inverse_side));

        continue;
      }

      PFace positive_side, negative_side;

      find_adjacent_faces(pface, pedge, neighbor_faces, positive_side, negative_side);
      CGAL_assertion(positive_side != negative_side);

      if (volume_indices[0] != -1) {
        Oriented_side inverse_side = (positive_side.first == pface.first) ? ON_POSITIVE_SIDE : oriented_side(positive_side, pface);
        if (associate(positive_side, volume_indices[0], inverse_side, volumes, map_volumes))
          queue[0].push(std::make_pair(positive_side, inverse_side));
      }

      if (volume_indices[1] != -1) {
        Oriented_side inverse_side = (negative_side.first == pface.first) ? ON_NEGATIVE_SIDE : oriented_side(negative_side, pface);
        if (associate(negative_side, volume_indices[1], inverse_side, volumes, map_volumes))
          queue[1].push(std::make_pair(negative_side, inverse_side));
      }
    }

    // Propagate both queues if volumes on either side of the pface are not segmented.
    for (std::size_t i = 0; i < 2; i++) {
      if (volume_indices[i] != -1) {
        while (!queue[i].empty()) {
          propagate_volume(queue[i], volume_indices[i], volumes, map_volumes);
        }
      }
    }
  }

  void propagate_volume(std::queue<std::pair<PFace, Oriented_side> >& queue, std::size_t volume_index, std::vector<Volume_cell>& volumes, std::map<PFace, std::pair<int, int> >& map_volumes) {
    PFace pface;
    Oriented_side seed_side;
    std::tie(pface, seed_side) = queue.front();
    queue.pop();

    // Get neighbors with smallest dihedral angle
    const auto pedges = m_data.pedges_of_pface(pface);
    std::vector<PFace> neighbor_faces, adjacent_faces;
    for (const auto pedge : pedges) {
      CGAL_assertion(m_data.has_iedge(pedge));
      m_data.incident_faces(m_data.iedge(pedge), neighbor_faces);

      if (neighbor_faces.size() == 2) {
        // If there is only one neighbor, the edge is on the corner of the bbox.
        // Thus the only neighbor needs to be a bbox face.
        PFace neighbor = (neighbor_faces[0] == pface) ? neighbor_faces[1] : neighbor_faces[0];
        CGAL_assertion(neighbor.first < 6 && pface.first < 6);
        CGAL_assertion(oriented_side(pface, neighbor) == seed_side);

        Oriented_side inverse_side = oriented_side(neighbor, pface);

        CGAL_assertion(inverse_side == ON_POSITIVE_SIDE);

        if (associate(neighbor, volume_index, inverse_side, volumes, map_volumes))
          queue.push(std::make_pair(neighbor, inverse_side));
        continue;
      }

      PFace positive_side, negative_side;

      find_adjacent_faces(pface, pedge, neighbor_faces, positive_side, negative_side);
      CGAL_assertion(positive_side != negative_side);

      if (seed_side == ON_POSITIVE_SIDE) {
        Oriented_side inverse_side = (pface.first == positive_side.first) ? seed_side : oriented_side(positive_side, pface);
        if (associate(positive_side, volume_index, inverse_side, volumes, map_volumes))
          queue.push(std::make_pair(positive_side, inverse_side));
      }
      else {
        Oriented_side inverse_side = (pface.first == negative_side.first) ? seed_side : oriented_side(negative_side, pface);
        if (associate(negative_side, volume_index, inverse_side, volumes, map_volumes))
          queue.push(std::make_pair(negative_side, inverse_side));
      }
    }
  }

  bool associate(const PFace& pface, std::size_t volume_index, Oriented_side side, std::vector<Volume_cell>& volumes, std::map<PFace, std::pair<int, int> >& map_volumes) {
    auto& pair = map_volumes.at(pface);

    CGAL_assertion(side != COPLANAR);

    if (side == ON_POSITIVE_SIDE) {
      if (pair.first == static_cast<int>(volume_index))
        return false;

      pair.first = static_cast<int>(volume_index);
      volumes[volume_index].add_pface(pface, pair.second);
      return true;
    }
    else if (side == ON_NEGATIVE_SIDE) {
      if (pair.second == static_cast<int>(volume_index))
        return false;

      pair.second = static_cast<int>(volume_index);
      volumes[volume_index].add_pface(pface, pair.first);
      return true;
    }

    CGAL_assertion(false);
    return false;
  }

  void find_adjacent_faces(const PFace& pface, const PEdge& pedge, const std::vector<PFace>& neighbor_faces, PFace& positive_side, PFace& negative_side) const {
    CGAL_assertion(neighbor_faces.size() > 2);

    // for each face, find vertex that is not collinear with the edge
    // take 2d directions orthogonal to edge
    // sort and take neighbors to current face

    const Segment_3 segment = m_data.segment_3(pedge);
    Vector_3 norm(segment.source(), segment.target());
    norm = KSP::internal::normalize(norm);
    const Plane_3 plane(segment.source(), norm);
    Point_2 source2d = plane.to_2d(segment.source());

    std::vector< std::pair<Direction_2, PFace> > dir_edges;
    // Get orientation towards edge of current face
    Point_2 v2d = plane.to_2d(m_data.centroid_of_pface(pface));
    dir_edges.push_back(std::make_pair(Direction_2(source2d - v2d), pface));

    // Get orientation towards edge of other faces
    for (const PFace& face : neighbor_faces) {
      if (face == pface)
        continue;

      // Taking just the direction of the line instead of the point? (Still need to take care of the sign)
      auto& sp = m_data.support_plane(face.first);
      auto& mesh = sp.mesh();
      auto h = mesh.halfedge(face.second);
      auto first = h;

      Point_3 point;
      FT dist = 0;
      do {
        Point_3 p = sp.to_3d(mesh.point(mesh.target(h)));
        Vector_3 dist_in_plane = (p - segment.source());
        dist_in_plane -= norm * (dist_in_plane * norm);
        FT d = dist_in_plane.squared_length();
        if (d > dist) {
          dist = d;
          point = p;
        }
        h = mesh.next(h);
      } while (first != h);

      Point_2 p = plane.to_2d(point);

      dir_edges.push_back(std::make_pair(Direction_2(source2d - p), face));
    }

    CGAL_assertion(dir_edges.size() == neighbor_faces.size());

    // Sort directions
    std::sort(dir_edges.begin(), dir_edges.end(), [&](
      const std::pair<Direction_2, PFace>& p,
      const std::pair<Direction_2, PFace>& q) -> bool {
        return p.first < q.first;
      }
    );

    std::size_t n = dir_edges.size();
    for (std::size_t i = 0; i < n; ++i) {
      if (dir_edges[i].second == pface) {

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        // Check which side the faces are on
        Orientation im_side = oriented_side(pface, dir_edges[im].second);
        Orientation ip_side = oriented_side(pface, dir_edges[ip].second);

        //The angles decide where to push them, not necessarily the side.
        //At least in case of a bbox corner intersected with a third plane breaks this.

        // Both on the negative side is not possible
        CGAL_assertion(im_side != ON_NEGATIVE_SIDE || ip_side != ON_NEGATIVE_SIDE);
        CGAL_assertion(im_side != COPLANAR || ip_side != COPLANAR);

        // If two are positive it has to be a bbox corner situation.
        if (im_side == ON_POSITIVE_SIDE && ip_side == ON_POSITIVE_SIDE) {
          CGAL_assertion(pface.first < 6);
          if (dir_edges[im].second.first < 6) {
            positive_side = dir_edges[ip].second;
            negative_side = dir_edges[im].second;
          }
          else if (dir_edges[ip].second.first < 6) {
            positive_side = dir_edges[im].second;
            negative_side = dir_edges[ip].second;
          }
          else CGAL_assertion(false);
        }
        else if (ip_side == ON_POSITIVE_SIDE || im_side == ON_NEGATIVE_SIDE) {
          positive_side = dir_edges[ip].second;
          negative_side = dir_edges[im].second;
        }
        else {
          positive_side = dir_edges[im].second;
          negative_side = dir_edges[ip].second;
        }

        return;
      }
    }

    CGAL_assertion_msg(false, "ERROR: NEXT PFACE IS NOT FOUND!");
  }

  Oriented_side oriented_side(const PFace& a, const PFace& b) const {
    FT max_dist = 0;
    if (a.first == b.first)
      return COPLANAR;
    Oriented_side side;
    const Plane_3& p = m_data.support_plane(a.first).plane();
    for (auto v : m_data.pvertices_of_pface(b)) {
      Point_3 pt = m_data.point_3(v);
      FT dist = (p.point() - pt) * p.orthogonal_vector();
      if (CGAL::abs(dist) > max_dist) {
        side = p.oriented_side(m_data.point_3(v));
        max_dist = CGAL::abs(dist);
      }
    }

    return side;
  }

  void merge_facets_connected_components() {
    // Purpose: merge facets between the same volumes. Every pair of volumes can have at most one contact polygon (which also has to be convex)
    // Precondition: all volumes are convex, the contact area between each pair of volumes is empty or convex
    std::vector<E_constraint_map> edge_constraint_maps(m_data.number_of_support_planes());

    for (std::size_t sp = 0; sp < m_data.number_of_support_planes(); sp++) {
      //dump_2d_surface_mesh(m_data, sp, "face_merge/" + m_data.prefix() + std::to_string(sp) + "-before");
      typename Support_plane::Mesh& mesh = m_data.support_plane(sp).mesh();

      edge_constraint_maps[sp] = mesh.template add_property_map<typename Support_plane::Edge_index, bool>("e:keep", true).first;
      F_component_map fcm = mesh.template add_property_map<typename Support_plane::Face_index, typename boost::graph_traits<typename Support_plane::Mesh>::faces_size_type>("f:component", 0).first;

      std::size_t num = 0;

      for (auto e : mesh.edges()) {
        IEdge iedge = m_data.iedge(PEdge(sp, e));

        if (is_occupied(iedge, sp)) {
          edge_constraint_maps[sp][e] = true;
          num++;
        }
        else
          edge_constraint_maps[sp][e] = false;
      }

      CGAL::Polygon_mesh_processing::connected_components(mesh, fcm, CGAL::parameters::edge_is_constrained_map(edge_constraint_maps[sp]));

      merge_connected_components(sp, mesh, fcm, edge_constraint_maps[sp]);

      mesh.collect_garbage();
    }
  }

  void merge_connected_components(std::size_t sp_idx, typename Support_plane::Mesh& mesh, F_component_map& fcm, E_constraint_map ecm) {
    using Halfedge = typename Support_plane::Halfedge_index;
    using Face = typename Support_plane::Face_index;

    std::vector<bool> initial_component;
    std::size_t num_components = 0;

    for (const auto& f : mesh.faces())
      num_components = (std::max<std::size_t>)(num_components, fcm[f]);

    initial_component.resize(num_components + 1, false);
    Support_plane& sp = m_data.support_plane(sp_idx);
    for (const auto& f : mesh.faces()) {
      if (sp.is_initial(f))
        initial_component[fcm[f]] = true;
    }

    //std::vector<Halfedge> remove_edges;
    std::vector<bool> remove_vertices(mesh.vertices().size(), true);
    std::vector<bool> remove_faces(mesh.faces().size(), true);

    std::vector<bool> visited_halfedges(mesh.halfedges().size(), false);
    for (auto h : mesh.halfedges()) {
      Point_3 s = sp.to_3d(mesh.point(mesh.source(h)));
      Point_3 t = sp.to_3d(mesh.point(mesh.target(h)));
      std::vector<Point_3> pts;
      pts.push_back(s);
      pts.push_back(t);
      if (visited_halfedges[h])
        continue;

      visited_halfedges[h] = true;
      // Skip the outside edges.
      if (mesh.face(h) == mesh.null_face())
        continue;

      Face f0 = mesh.face(h);

      typename boost::graph_traits<typename Support_plane::Mesh>::faces_size_type c0 = fcm[f0], c_other;
      Face f_other = mesh.face(mesh.opposite(h));

      // Check whether the edge is between different components.
      if (f_other != mesh.null_face()) {
        if (c0 == fcm[f_other]) {
          continue;
        }
      }

      set_halfedge(f0, h, mesh);
      remove_faces[f0] = false;
      if (initial_component[c0])
        sp.set_initial(f0);

      // Find halfedge loop around component.
      std::vector<Halfedge> loop;
      loop.push_back(h);
      remove_vertices[mesh.target(h)] = false;
      Halfedge first = h;
      do {
        Halfedge n = h;
        //Point_3 tn = m_data.support_plane(sp).to_3d(mesh.point(mesh.target(h)));

        do {
          if (n == h)
            n = mesh.next(n);
          else
            n = mesh.next(mesh.opposite(n));

          //Point_3 tn2 = m_data.support_plane(sp).to_3d(mesh.point(mesh.target(h)));
          visited_halfedges[n] = true;

          Face_index fn = mesh.face(n);

          f_other = mesh.face(mesh.opposite(n));
          if (f_other == mesh.null_face())
            break;
          c_other = fcm[f_other];
          if (c0 == c_other && ecm[Edge_index(n >> 1)])
            std::cout << "edge and face constraint map inconsistent1" << std::endl;

          if (c0 != c_other && !ecm[Edge_index(n >> 1)])
            std::cout << "edge and face constraint map inconsistent2" << std::endl;
        } while (c0 == c_other && n != h);

        if (n == h) {
          // Should not happen.
          std::cout << "Searching for next edge of connected component failed" << std::endl;
        }
        loop.push_back(n);
        set_face(n, f0, mesh);
        set_next(h, n, mesh);
        set_halfedge(mesh.target(h), h, mesh);
        remove_vertices[mesh.target(n)] = false;
        h = n;
      } while (h != first);
    }

    // Remove all vertices in remove_vertices and all edges marked in constrained list
    for (auto f : mesh.faces()) {
      if (remove_faces[f])
        remove_face(f, mesh);
    }

    for (auto e : mesh.edges()) {
      // If edge is not marked as constrained it can be removed.
      if (!ecm[e]) {
        remove_edge(e, mesh);
      }
    }

    for (auto v : mesh.vertices()) {
      if (remove_vertices[v])
        remove_vertex(v, mesh);
    }

    if (!mesh.is_valid(true)) {
      std::cout << "mesh is not valid after merging faces of sp " << sp_idx << std::endl;
    }
  }

  bool is_occupied(IEdge iedge, std::size_t sp) {
    const std::set<std::size_t> planes = m_data.intersected_planes(iedge);
    for (std::size_t j : planes) {
      if (sp == j)
        continue;

      m_data.support_plane(j).mesh().is_valid(true);

      for (auto e2 : m_data.support_plane(j).mesh().edges()) {
        if (iedge == m_data.iedge(PEdge(j, e2))) {
          return true;
        }
      }
    }

    return false;
  }

  void create_cell_pvertices(Volume_cell& cell) {
    From_exact from_exact;
    std::vector<int>& ivertex2vertex = m_data.ivertex_to_index();
    std::vector<Point_3>& vertices = m_data.vertices();
    std::vector<typename Intersection_kernel::Point_3>& exact_vertices = m_data.exact_vertices();
    std::vector<std::vector<std::size_t> >& face2vertices = m_data.face_to_vertices();
    std::vector<std::size_t>& face2sp = m_data.face_to_support_plane();
    cell.pvertices.clear();
    for (std::size_t f = 0; f < cell.pfaces.size(); f++) {
      const auto& pface = cell.pfaces[f];
      bool face_filled = !face2vertices[cell.faces[f]].empty();

      if (!face_filled) {
        face2vertices[cell.faces[f]].reserve(m_data.pvertices_of_pface(pface).size());
        face2sp[cell.faces[f]] = pface.first;
      }

      for (const auto pvertex : m_data.pvertices_of_pface(pface)) {
        CGAL_assertion(m_data.has_ivertex(pvertex));
        cell.pvertices.insert(pvertex);

        IVertex ivertex = m_data.ivertex(pvertex);
        if (ivertex2vertex[ivertex] == -1) {
          ivertex2vertex[ivertex] = static_cast<int>(vertices.size());
          if (!face_filled)
            face2vertices[cell.faces[f]].push_back(vertices.size());
          else
            std::cout << "Should not happen" << std::endl;
          vertices.push_back(from_exact(m_data.point_3(ivertex)));
          exact_vertices.push_back(m_data.point_3(ivertex));
        }
        else if (!face_filled)
          face2vertices[cell.faces[f]].push_back(ivertex2vertex[ivertex]);
      }
    }
  }

  void remove_collinear_vertices() {
    auto& vertices = m_data.face_to_vertices();
    std::vector<bool> coll(m_data.exact_vertices().size(), true);
    std::unordered_map<std::size_t, std::vector<std::size_t> > vtx2face;

    for (std::size_t f = 0; f < vertices.size(); f++) {
      for (std::size_t i = 0; i < vertices[f].size(); i++) {
        if (!coll[vertices[f][i]])
          continue;
        const typename Intersection_kernel::Point_3& a = m_data.exact_vertices()[vertices[f][(i - 1 + vertices[f].size()) % vertices[f].size()]];
        const typename Intersection_kernel::Point_3& b = m_data.exact_vertices()[vertices[f][i]];
        const typename Intersection_kernel::Point_3& c = m_data.exact_vertices()[vertices[f][(i + 1) % vertices[f].size()]];
        if (!CGAL::collinear(a, b, c))
          coll[vertices[f][i]] = false;
        else
          vtx2face[vertices[f][i]].push_back(f);
      }
    }

    for (std::size_t i = 0; i < coll.size(); i++) {
      if (!coll[i])
        continue;
      const auto& f = vtx2face[i];
      for (std::size_t j = 0; j < f.size(); j++) {
        for (std::size_t v = 0; v < vertices[f[j]].size(); v++) {
          if (vertices[f[j]][v] == i) {
            vertices[f[j]].erase(vertices[f[j]].begin() + v);
            break;
          }
        }
      }
    }
  }
};

#endif //DOXYGEN_RUNNING

} // namespace internal
} // namespace KSP_3
} // namespace CGAL

#endif // CGAL_KSP_3_FINALIZER_H
