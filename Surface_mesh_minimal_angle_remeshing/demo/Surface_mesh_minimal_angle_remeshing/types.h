#ifndef CGAL_TYPES_H
#define CGAL_TYPES_H

// CGAL
#include <CGAL/Simple_cartesian.h>
// Local
#include "Minangle_remesh.h"

typedef CGAL::Color Color;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef PMP::Minangle_remesh<Kernel> Minangle_remesh;
typedef PMP::internal::Minangle_remesher<Kernel> Minangle_remesher;
typedef Minangle_remesher::Mesh_properties_ Mesh_properties;
typedef Minangle_remesher::Bbox Bbox;         // in Minangle_remesher
typedef Minangle_remesher::FT FT;
typedef Minangle_remesher::Point Point;
typedef Minangle_remesher::Vector Vector;
typedef Minangle_remesher::Normal Normal;
typedef Minangle_remesher::Mesh Mesh;
typedef Minangle_remesher::Face_tree Face_tree;
typedef Minangle_remesher::halfedge_descriptor halfedge_descriptor;
typedef Minangle_remesher::edge_descriptor edge_descriptor;
typedef Minangle_remesher::vertex_descriptor vertex_descriptor;
typedef Minangle_remesher::face_descriptor face_descriptor;
typedef Minangle_remesher::Halfedge_list Halfedge_list;
typedef Minangle_remesher::Edge_list Edge_list;
typedef Minangle_remesher::Vertex_list Vertex_list;
typedef Minangle_remesher::Face_list Face_list;
typedef Minangle_remesher::Point_list Point_list;
typedef Minangle_remesher::Point_iter Point_iter;
typedef Minangle_remesher::Point_const_iter Point_const_iter;
typedef Minangle_remesher::Color_list Color_list;
typedef Minangle_remesher::Color_iter Color_iter;
typedef Minangle_remesher::Color_const_iter Color_const_iter;
typedef Minangle_remesher::Point_Comp Point_Comp;
typedef Minangle_remesher::Point_pair Point_pair;
typedef Minangle_remesher::Link Link;
typedef Minangle_remesher::Link_list Link_list;
typedef Minangle_remesher::Link_list_iter Link_list_iter;
typedef Minangle_remesher::Link_list_const_iter Link_list_const_iter;
typedef Minangle_remesher::Link_iter_list_iter Link_iter_list_iter;
typedef Minangle_remesher::Link_iter_list_const_iter Link_iter_list_const_iter;
typedef Minangle_remesher::Link_pointer_list Link_pointer_list;
typedef Minangle_remesher::Link_pointer_iter Link_pointer_iter;
typedef Minangle_remesher::Link_pointer_const_iter Link_pointer_const_iter;
typedef Mesh_properties::Bvd Bvd;
typedef Mesh_properties::Halfedge_around_target_circulator Halfedge_around_target_circulator;

enum DrawType {
  k_mesh = 0,
  k_all_voronoi,
  k_vertex_voronoi,
  k_edge_voronoi,
  k_face_voronoi
};

enum RenderType {
  k_plain_faces = 0, // faces with plain color
  k_ifi_faces,       // faces with interpolated feature intensity colors
  k_mr_faces,
  k_classifications,  // vertex type (feature, crease, smooth)
  k_gaussian_curvature,
  k_maximal_halfedge_dihedral,
  k_normal_dihedral,
  k_feature_intensity,
  k_capacity,
  k_weight
};

struct halfedge2edge {
  halfedge2edge(const Mesh &mesh, std::vector<edge_descriptor> &edges)
  : m_mesh(mesh), m_edges(edges) {
  }

  void operator() (const halfedge_descriptor &h) const {
    m_edges.push_back(edge(h, m_mesh));
  }

  const Mesh &m_mesh;
  std::vector<edge_descriptor> &m_edges;
};

#endif // CGAL_TYPES_H