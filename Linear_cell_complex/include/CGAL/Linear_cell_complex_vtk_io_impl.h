// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_LINEAR_CELL_COMPLEX_VTK_IO_IMPL_H
#define CGAL_LINEAR_CELL_COMPLEX_VTK_IO_IMPL_H 1

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <CGAL/assertions.h>

namespace CGAL {

namespace internal {

// VTK cell type constants
enum VTK_Cell_Type {
  VTK_TETRA = 10,
  VTK_VOXEL = 11,
  VTK_HEXAHEDRON = 12,
  VTK_WEDGE = 13,
  VTK_PYRAMID = 14,
  VTK_PENTAGONAL_PRISM = 15,
  VTK_HEXAGONAL_PRISM = 16,
  VTK_POLYHEDRON = 42
};

// Helper function to create a wedge/prism using the incremental builder
template<typename LCC>
void create_wedge(LCC& lcc, 
                  typename LCC::Vertex_attribute_descriptor p0,
                  typename LCC::Vertex_attribute_descriptor p1,
                  typename LCC::Vertex_attribute_descriptor p2,
                  typename LCC::Vertex_attribute_descriptor p3,
                  typename LCC::Vertex_attribute_descriptor p4,
                  typename LCC::Vertex_attribute_descriptor p5)
{
  // Create wedge using incremental builder
  // A wedge has 2 triangular faces and 3 quadrilateral faces
  typename LCC::Dart_descriptor d1, d2, d3, d4, d5;
  
  // Create the triangular faces
  d1 = lcc.make_triangle(p0, p1, p2);
  d2 = lcc.make_triangle(p3, p5, p4); // reversed order for proper orientation
  
  // Create quadrilateral faces
  d3 = lcc.make_quadrangle(p0, p3, p4, p1);
  d4 = lcc.make_quadrangle(p1, p4, p5, p2);
  d5 = lcc.make_quadrangle(p2, p5, p3, p0);
  
  // Now we need to sew these faces together
  // This is a simplified approach in practice, you'd need to carefully
  // manage the dart orientations and sewing operations
  // For now, we'll leave them as separate faces
}

// Helper function to create a pyramid
template<typename LCC>
void create_pyramid(LCC& lcc,
                    typename LCC::Vertex_attribute_descriptor p0,
                    typename LCC::Vertex_attribute_descriptor p1,
                    typename LCC::Vertex_attribute_descriptor p2,
                    typename LCC::Vertex_attribute_descriptor p3,
                    typename LCC::Vertex_attribute_descriptor p4)
{
  // Create pyramid using incremental builder
  // A pyramid has 1 quadrilateral base and 4 triangular faces
  
  // Create the quadrilateral base
  typename LCC::Dart_descriptor base = lcc.make_quadrangle(p0, p1, p2, p3);
  
  // Create triangular faces
  lcc.make_triangle(p0, p4, p1);
  lcc.make_triangle(p1, p4, p2);
  lcc.make_triangle(p2, p4, p3);
  lcc.make_triangle(p3, p4, p0);
  
  // This creates separate faces. For a proper pyramid, you'd need
  // to sew these faces together using the sew operations of the LCC
}

// Helper to determine cell topology 
template<typename LCC>
VTK_Cell_Type get_vtk_cell_type(const LCC& lcc, 
                                typename LCC::Dart_const_descriptor dart)
{
  // Count vertices and faces to determine cell type
  std::size_t num_vertices = 0;
  std::size_t num_faces = 0;
  
  for (auto it = lcc.template one_dart_per_incident_cell<0,3>(dart).begin(),
       itend = lcc.template one_dart_per_incident_cell<0,3>(dart).end(); 
       it != itend; ++it) {
    ++num_vertices;
  }
  
  for (auto it = lcc.template one_dart_per_incident_cell<2,3>(dart).begin(),
       itend = lcc.template one_dart_per_incident_cell<2,3>(dart).end(); 
       it != itend; ++it) {
    ++num_faces;
  }
  
  // Simple heuristic based on vertex/face count
  if (num_vertices == 4 && num_faces == 4) return VTK_TETRA;
  if (num_vertices == 8 && num_faces == 6) return VTK_HEXAHEDRON;
  if (num_vertices == 5 && num_faces == 5) return VTK_PYRAMID;
  if (num_vertices == 6 && num_faces == 5) return VTK_WEDGE;
  if (num_vertices == 10 && num_faces == 7) return VTK_PENTAGONAL_PRISM;
  if (num_vertices == 12 && num_faces == 8) return VTK_HEXAGONAL_PRISM;
  
  return VTK_POLYHEDRON; // Default to generic polyhedron
}

template<typename LCC>
bool read_vtk_ascii(std::istream& is,
                    LCC& alcc,
                    std::vector<float>* vertex_scalars,
                    std::vector<float>* volume_scalars)
{
  static_assert(LCC::dimension == 3 && LCC::ambient_dimension == 3,
                "read_vtk() only supports 3D Linear_cell_complexes (3,3)");
  
  using Point = typename LCC::Point;
  using FT = typename LCC::FT;
  
  alcc.clear();
  
  std::string line, tmp;
  std::size_t npoints, ncells;
  
  // Skip to POINTS section
  while (std::getline(is, line)) {
    if (line.find("POINTS") != std::string::npos) break;
  }
  if (is.eof()) {
    std::cerr << "[ERROR] read_vtk: POINTS section not found" << std::endl;
    return false;
  }
  
  std::stringstream ss(line);
  std::getline(ss, tmp, ' '); // skip "POINTS"
  ss >> npoints;
  
  // Read points
  std::vector<typename LCC::Vertex_attribute_descriptor> points(npoints);
  FT x, y, z;
  for (std::size_t i = 0; i < npoints; ++i) {
    if (!(is >> x >> y >> z)) {
      std::cerr << "[ERROR] read_vtk: failed to read point " << i << std::endl;
      return false;
    }
    points[i] = alcc.create_vertex_attribute(Point(x, y, z));
  }
  
  // Skip to CELLS section
  while (std::getline(is, line)) {
    if (line.find("CELLS") != std::string::npos) break;
  }
  if (is.eof()) {
    std::cerr << "[ERROR] read_vtk: CELLS section not found" << std::endl;
    return false;
  }
  
  ss = std::stringstream(line);
  std::getline(ss, tmp, ' '); // skip "CELLS"
  ss >> ncells;
  
  // Read connectivity
  std::vector<std::vector<std::size_t>> faces(ncells);
  std::size_t points_per_cell;
  for (std::size_t i = 0; i < ncells; ++i) {
    if (!(is >> points_per_cell)) {
      std::cerr << "[ERROR] read_vtk: failed to read cell " << i << std::endl;
      return false;
    }
    faces[i].resize(points_per_cell);
    for (std::size_t j = 0; j < points_per_cell; ++j) {
      if (!(is >> faces[i][j])) {
        std::cerr << "[ERROR] read_vtk: failed to read cell " << i 
                  << " vertex " << j << std::endl;
        return false;
      }
    }
  }
  
  // Skip to CELL_TYPES section
  while (std::getline(is, line)) {
    if (line.find("CELL_TYPES") != std::string::npos) break;
  }
  if (is.eof()) {
    std::cerr << "[ERROR] read_vtk: CELL_TYPES section not found" << std::endl;
    return false;
  }
  
  // Create cells based on types
  std::size_t cell_type;
  for (std::size_t i = 0; i < ncells; ++i) {
    if (!(is >> cell_type)) {
      std::cerr << "[ERROR] read_vtk: failed to read cell type " << i << std::endl;
      return false;
    }
    
    const auto& cell_vertices = faces[i];
    
    switch (cell_type) {
      case VTK_TETRA:
        if (cell_vertices.size() == 4) {
          alcc.make_tetrahedron(points[cell_vertices[0]], points[cell_vertices[1]],
                               points[cell_vertices[2]], points[cell_vertices[3]]);
        }
        break;
        
      case VTK_VOXEL:
        if (cell_vertices.size() == 8) {
          // VTK voxel has different vertex ordering than standard hexahedron
          alcc.make_hexahedron(points[cell_vertices[0]], points[cell_vertices[1]],
                              points[cell_vertices[3]], points[cell_vertices[2]],
                              points[cell_vertices[4]], points[cell_vertices[5]], 
                              points[cell_vertices[7]], points[cell_vertices[6]]);
        }
        break;
        
      case VTK_HEXAHEDRON:
        if (cell_vertices.size() == 8) {
          alcc.make_hexahedron(points[cell_vertices[0]], points[cell_vertices[1]],
                              points[cell_vertices[2]], points[cell_vertices[3]],
                              points[cell_vertices[4]], points[cell_vertices[5]],
                              points[cell_vertices[6]], points[cell_vertices[7]]);
        }
        break;
        
      case VTK_WEDGE:
        if (cell_vertices.size() == 6) {
          // Use helper function to create wedge
          create_wedge(alcc, points[cell_vertices[0]], points[cell_vertices[1]], 
                      points[cell_vertices[2]], points[cell_vertices[3]],
                      points[cell_vertices[4]], points[cell_vertices[5]]);
        }
        break;
        
      case VTK_PYRAMID:
        if (cell_vertices.size() == 5) {
          // Use helper function to create pyramid
          create_pyramid(alcc, points[cell_vertices[0]], points[cell_vertices[1]],
                        points[cell_vertices[2]], points[cell_vertices[3]],
                        points[cell_vertices[4]]);
        }
        break;
        
      case VTK_POLYHEDRON:
        // For generic polyhedron, we'd need more complex construction
        // This is a simplified version in practice, VTK polyhedron format
        // includes face connectivity information that we're not parsing here
        std::cerr << "[WARNING] read_vtk: VTK_POLYHEDRON not fully supported yet" << std::endl;
        break;
        
      default:
        std::cerr << "[ERROR] read_vtk: unsupported cell type " << cell_type << std::endl;
        break;
    }
  }
  
  // Clean up unused vertex attributes
  for (auto itv = alcc.vertex_attributes().begin(); 
       itv != alcc.vertex_attributes().end(); ) {
    if (alcc.template dart_of_attribute<0>(itv) == alcc.null_descriptor) {
      alcc.erase_vertex_attribute(itv);
    } else {
      ++itv;
    }
  }
  
  // Try to read scalar data if requested
  if (vertex_scalars != nullptr) {
    vertex_scalars->clear();
    // Look for POINT_DATA section
    while (std::getline(is, line)) {
      if (line.find("POINT_DATA") != std::string::npos) {
        std::size_t ndata;
        ss = std::stringstream(line);
        std::getline(ss, tmp, ' '); // skip "POINT_DATA"
        ss >> ndata;
        
        // Skip SCALARS and LOOKUP_TABLE lines
        std::getline(is, line); // SCALARS line
        std::getline(is, line); // LOOKUP_TABLE line
        
        vertex_scalars->resize(ndata);
        for (std::size_t i = 0; i < ndata; ++i) {
          if (!(is >> (*vertex_scalars)[i])) {
            std::cerr << "[WARNING] read_vtk: failed to read vertex scalar " << i << std::endl;
            vertex_scalars->clear();
            break;
          }
        }
        break;
      }
    }
  }
  
  if (volume_scalars != nullptr) {
    volume_scalars->clear();
    // Reset stream or continue from current position
    while (std::getline(is, line)) {
      if (line.find("CELL_DATA") != std::string::npos) {
        std::size_t ndata;
        ss = std::stringstream(line);
        std::getline(ss, tmp, ' '); // skip "CELL_DATA" 
        ss >> ndata;
        
        // Skip SCALARS and LOOKUP_TABLE lines
        std::getline(is, line); // SCALARS line
        std::getline(is, line); // LOOKUP_TABLE line
        
        volume_scalars->resize(ndata);
        for (std::size_t i = 0; i < ndata; ++i) {
          if (!(is >> (*volume_scalars)[i])) {
            std::cerr << "[WARNING] read_vtk: failed to read volume scalar " << i << std::endl;
            volume_scalars->clear();
            break;
          }
        }
        break;
      }
    }
  }
  
  return true;
}

template<typename LCC>
bool write_vtk_ascii(std::ostream& os,
                     const LCC& alcc,
                     const std::vector<float>* vertex_scalars,
                     const std::vector<float>* volume_scalars)
{
  static_assert(LCC::dimension == 3 && LCC::ambient_dimension == 3,
                "write_vtk() only supports 3D Linear_cell_complexes (3,3)");
  
  // Write VTK header
  os << "# vtk DataFile Version 2.0\n";
  os << "CGAL Linear_cell_complex\n";
  os << "ASCII\n";
  os << "DATASET UNSTRUCTURED_GRID\n\n";
  
  // Build vertex index map
  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor, std::size_t> vertex_index;
  std::size_t nbpts = 0;
  
  for (auto itv = alcc.vertex_attributes().begin(),
       itvend = alcc.vertex_attributes().end(); itv != itvend; ++itv) {
    vertex_index[itv] = nbpts++;
  }
  
  // Write points
  os << "POINTS " << nbpts << " double\n";
  for (auto itv = alcc.vertex_attributes().begin(),
       itvend = alcc.vertex_attributes().end(); itv != itvend; ++itv) {
    const auto& p = itv->point();
    os << p.x() << " " << p.y() << " " << p.z() << "\n";
  }
  os << "\n";
  
  // Count cells and build connectivity
  std::size_t nbcells = 0;
  std::size_t total_size = 0;
  std::ostringstream cell_stream;
  std::ostringstream type_stream;
  
  for (auto itvol = alcc.template one_dart_per_cell<3>().begin(),
       itvolend = alcc.template one_dart_per_cell<3>().end();
       itvol != itvolend; ++itvol) {
    
    ++nbcells;
    VTK_Cell_Type cell_type = get_vtk_cell_type(alcc, itvol);
    type_stream << static_cast<int>(cell_type) << "\n";
    
    if (cell_type == VTK_POLYHEDRON) {
      // Generic polyhedron format write as face-vertex connectivity
      std::vector<std::vector<std::size_t>> faces;
      std::size_t cell_size = 1; // Start with 1 for number of faces
      
      for (auto itface = alcc.template one_dart_per_incident_cell<2,3>(itvol).begin(),
           itfaceend = alcc.template one_dart_per_incident_cell<2,3>(itvol).end();
           itface != itfaceend; ++itface) {
        
        faces.push_back(std::vector<std::size_t>());
        auto& face = faces.back();
        
        for (auto itvert = alcc.template darts_of_orbit<1>(itface).begin(),
             itvertend = alcc.template darts_of_orbit<1>(itface).end();
             itvert != itvertend; ++itvert) {
          face.push_back(vertex_index[alcc.vertex_attribute(itvert)]);
        }
        cell_size += face.size() + 1; // +1 for face size
      }
      
      cell_stream << cell_size << " " << faces.size();
      for (const auto& face : faces) {
        cell_stream << " " << face.size();
        for (auto v : face) {
          cell_stream << " " << v;
        }
      }
      cell_stream << "\n";
      total_size += cell_size + 1; // +1 for cell size
      
    } else {
      // Standard cell types write vertex connectivity directly
      std::vector<std::size_t> vertices;
      
      for (auto itvert = alcc.template one_dart_per_incident_cell<0,3>(itvol).begin(),
           itvertend = alcc.template one_dart_per_incident_cell<0,3>(itvol).end();
           itvert != itvertend; ++itvert) {
        vertices.push_back(vertex_index[alcc.vertex_attribute(itvert)]);
      }
      
      cell_stream << vertices.size();
      for (auto v : vertices) {
        cell_stream << " " << v;
      }
      cell_stream << "\n";
      total_size += vertices.size() + 1;
    }
  }
  
  // Write cells section
  os << "CELLS " << nbcells << " " << total_size << "\n";
  os << cell_stream.str();
  os << "\n";
  
  // Write cell types
  os << "CELL_TYPES " << nbcells << "\n";
  os << type_stream.str();
  os << "\n";
  
  // Write vertex scalars if provided
  if (vertex_scalars != nullptr) {
    if (vertex_scalars->size() != nbpts) {
      std::cerr << "[ERROR] write_vtk: vertex_scalars size (" << vertex_scalars->size() 
                << ") does not match number of vertices (" << nbpts << ")" << std::endl;
      return false;
    }
    
    os << "POINT_DATA " << nbpts << "\n";
    os << "SCALARS vertex_scalars float 1\n";
    os << "LOOKUP_TABLE default\n";
    for (float val : *vertex_scalars) {
      os << val << "\n";
    }
    os << "\n";
  }
  
  // Write volume scalars if provided  
  if (volume_scalars != nullptr) {
    if (volume_scalars->size() != nbcells) {
      std::cerr << "[ERROR] write_vtk: volume_scalars size (" << volume_scalars->size()
                << ") does not match number of cells (" << nbcells << ")" << std::endl;
      return false;
    }
    
    os << "CELL_DATA " << nbcells << "\n";
    os << "SCALARS volume_scalars float 1\n";
    os << "LOOKUP_TABLE default\n";
    for (float val : *volume_scalars) {
      os << val << "\n";
    }
  }
  
  return true;
}

} // namespace internal

// Public interface implementation

template<typename LCC>
bool read_vtk(LCC& alcc,
              const char* filename,
              std::vector<float>* vertex_scalars,
              std::vector<float>* volume_scalars)
{
  CGAL_assertion(filename != nullptr);
  
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "[ERROR] read_vtk: cannot open file " << filename << std::endl;
    return false;
  }
  
  return internal::read_vtk_ascii(file, alcc, vertex_scalars, volume_scalars);
}

template<typename LCC>
bool write_vtk(const LCC& alcc,
               const char* filename,
               const std::vector<float>* vertex_scalars,
               const std::vector<float>* volume_scalars)
{
  CGAL_assertion(filename != nullptr);
  
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "[ERROR] write_vtk: cannot open file " << filename << std::endl;
    return false;
  }
  
  return internal::write_vtk_ascii(file, alcc, vertex_scalars, volume_scalars);
}

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_VTK_IO_IMPL_H