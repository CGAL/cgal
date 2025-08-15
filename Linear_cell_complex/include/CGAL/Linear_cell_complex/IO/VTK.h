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

#ifndef CGAL_LCC_IO_VTK_H
#define CGAL_LCC_IO_VTK_H 1
#include <CGAL/Linear_cell_complex_incremental_builder_3.h>
#include <CGAL/assertions.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace CGAL {
namespace IO {

/** \file VTK.h
 * Functions to import/export 3D Linear_cell_complex from/to VTK legacy ASCII format.
 *
 * Only supports:
 * - Linear_cell_complex_for_combinatorial_map<3,3>
 * - VTK legacy ASCII format (.vtk files)
 * - Optional scalar fields for vertices and volumes
 *
 * Supported VTK cell types:
 * - VTK_TETRA (10): Tetrahedron
 * - VTK_VOXEL (11): Voxel (special hexahedron ordering)
 * - VTK_HEXAHEDRON (12): Hexahedron
 * - VTK_WEDGE (13): Prism/Wedge
 * - VTK_PYRAMID (14): Pyramid
 * - VTK_PENTAGONAL_PRISM (15): Pentagonal prism
 * - VTK_HEXAGONAL_PRISM (16): Hexagonal prism
 * - VTK_POLYHEDRON (42): Generic polyhedron
 */

// ============================================================================
//                              Declarations
// ============================================================================

/**
 * \brief Read a VTK legacy ASCII file and load it into a 3D Linear_cell_complex.
 * \ingroup PkgLinearCellComplexExamples
 *
 * \tparam LCC must be a Linear_cell_complex_for_combinatorial_map<3,3>
 * \tparam VertexScalarType Type for vertex scalar data (default: float)
 * \tparam VolumeScalarType Type for volume scalar data (default: float)
 * \param alcc The Linear_cell_complex to populate (will be cleared first)
 * \param filename Path to the VTK file
 * \param vertex_scalars Optional output vector to store per-vertex scalar values.
 *                      If provided, will be resized to match number of vertices.
 * \param volume_scalars Optional output vector to store per-volume scalar values.
 *                      If provided, will be resized to match number of volumes.
 * \return `true` if loading was successful, `false` otherwise
 */
template <typename LCC, typename VertexScalarType = float, typename VolumeScalarType = float>
bool read_VTK(LCC& alcc,
              const char* filename,
              std::vector<VertexScalarType>* vertex_scalars = nullptr,
              std::vector<VolumeScalarType>* volume_scalars = nullptr);

/**
 * \brief Write a 3D Linear_cell_complex to a VTK legacy ASCII file.
 * \ingroup PkgLinearCellComplexExamples
 *
 * \tparam LCC must be a Linear_cell_complex_for_combinatorial_map<3,3>
 * \tparam VertexScalarType Type for vertex scalar data (default: float)
 * \tparam VolumeScalarType Type for volume scalar data (default: float)
 * \param alcc The Linear_cell_complex to export
 * \param filename Path to the output VTK file
 * \param vertex_scalars Optional per-vertex scalar data. If provided, must have
 *                      same size as number of vertex attributes in the LCC.
 * \param volume_scalars Optional per-volume scalar data. If provided, must have
 *                      same size as number of 3-cells in the LCC.
 * \return `true` if writing was successful, `false` otherwise
 */
template <typename LCC, typename VertexScalarType = float, typename VolumeScalarType = float>
bool write_VTK(const LCC& alcc,
               const char* filename,
               const std::vector<VertexScalarType>* vertex_scalars = nullptr,
               const std::vector<VolumeScalarType>* volume_scalars = nullptr);

// Advanced versions with functors
template <typename LCC, typename VertexScalarReader, typename CellScalarReader>
bool read_VTK(LCC& alcc, const char* filename, VertexScalarReader vertex_reader, CellScalarReader cell_reader);

template <typename LCC, typename PointFunctor, typename CellFunctor>
bool write_VTK(const LCC& alcc, const char* filename, PointFunctor ptval, CellFunctor cellval);

// ============================================================================
//                          Implementation details
// ============================================================================

namespace internal {

// Helper: read a scalar of the right type from stream
template <typename Functor>
bool read_scalar_by_vtk_type(std::istream& is, const std::string& vtk_type, std::size_t n, Functor f) {
  if(vtk_type == "float") {
    for(std::size_t i = 0; i < n; ++i) {
      float v;
      if(!(is >> v))
        return false;
      f(i, v);
    }
  } else if(vtk_type == "double") {
    for(std::size_t i = 0; i < n; ++i) {
      double v;
      if(!(is >> v))
        return false;
      f(i, v);
    }
  } else if(vtk_type == "int") {
    for(std::size_t i = 0; i < n; ++i) {
      int v;
      if(!(is >> v))
        return false;
      f(i, v);
    }
  } else if(vtk_type == "unsigned_int") {
    for(std::size_t i = 0; i < n; ++i) {
      unsigned int v;
      if(!(is >> v))
        return false;
      f(i, v);
    }
  } else if(vtk_type == "short") {
    for(std::size_t i = 0; i < n; ++i) {
      short int v;
      if(!(is >> v))
        return false;
      f(i, v);
    }
  } else if(vtk_type == "unsigned_short") {
    for(std::size_t i = 0; i < n; ++i) {
      unsigned short int v;
      if(!(is >> v))
        return false;
      f(i, v);
    }
  } else if(vtk_type == "char") {
    for(std::size_t i = 0; i < n; ++i) {
      char v;
      int tmp;
      if(!(is >> tmp))
        return false;
      v = static_cast<char>(tmp);
      f(i, v);
    }
  } else if(vtk_type == "unsigned_char") {
    for(std::size_t i = 0; i < n; ++i) {
      unsigned char v;
      int tmp;
      if(!(is >> tmp))
        return false;
      v = static_cast<unsigned char>(tmp);
      f(i, v);
    }
  } else if(vtk_type == "long") {
    for(std::size_t i = 0; i < n; ++i) {
      long int v;
      if(!(is >> v))
        return false;
      f(i, v);
    }
  } else if(vtk_type == "unsigned_long") {
    for(std::size_t i = 0; i < n; ++i) {
      unsigned long int v;
      if(!(is >> v))
        return false;
      f(i, v);
    }
  } else {
    std::cerr << "[ERROR] read_VTK: unsupported scalar type: " << vtk_type << std::endl;
    return false;
  }
  return true;
}

// VTK type name mapping
template <typename T> struct gettype
{
  static std::string name() { return "unknown"; }
};
template <> struct gettype<bool>
{
  static std::string name() { return "bit"; }
};
template <> struct gettype<unsigned char>
{
  static std::string name() { return "unsigned_char"; }
};
template <> struct gettype<char>
{
  static std::string name() { return "char"; }
};
template <> struct gettype<unsigned short int>
{
  static std::string name() { return "unsigned_short"; }
};
template <> struct gettype<short int>
{
  static std::string name() { return "short"; }
};
template <> struct gettype<unsigned int>
{
  static std::string name() { return "unsigned_int"; }
};
template <> struct gettype<int>
{
  static std::string name() { return "int"; }
};
template <> struct gettype<unsigned long int>
{
  static std::string name() { return "unsigned_long"; }
};
template <> struct gettype<long int>
{
  static std::string name() { return "long"; }
};
template <> struct gettype<float>
{
  static std::string name() { return "float"; }
};
template <> struct gettype<double>
{
  static std::string name() { return "double"; }
};

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

// Helper: detect VTK cell type from a 3-cell
template <typename LCC, typename Dart> inline VTK_Cell_Type get_vtk_cell_type(const LCC& lcc, Dart itvol) {
  // Heuristic: count number of vertices and faces
  std::size_t nbv = 0, nbf = 0;
  for(auto itv = lcc.template one_dart_per_incident_cell<0, 3>(itvol).begin(),
           itvend = lcc.template one_dart_per_incident_cell<0, 3>(itvol).end();
      itv != itvend; ++itv)
    ++nbv;
  for(auto itf = lcc.template one_dart_per_incident_cell<2, 3>(itvol).begin(),
           itfend = lcc.template one_dart_per_incident_cell<2, 3>(itvol).end();
      itf != itfend; ++itf)
    ++nbf;

  if(nbv == 4 && nbf == 4)
    return VTK_TETRA;
  if(nbv == 5 && nbf == 5)
    return VTK_PYRAMID;
  if(nbv == 6 && nbf == 5)
    return VTK_WEDGE;
  if(nbv == 8 && nbf == 6)
    return VTK_HEXAHEDRON;
  if(nbv == 10 && nbf == 7)
    return VTK_PENTAGONAL_PRISM;
  if(nbv == 12 && nbf == 8)
    return VTK_HEXAGONAL_PRISM;
  return VTK_POLYHEDRON;
}

template <typename LCC, typename VertexScalarReader, typename CellScalarReader>
bool read_lcc_from_vtk_ascii(std::istream& is, LCC& alcc, VertexScalarReader vertex_reader, CellScalarReader cell_reader) {
  static_assert(LCC::dimension == 3 && LCC::ambient_dimension == 3,
                "read_VTK() only supports 3D Linear_cell_complexes (3,3)");

  using Point = typename LCC::Point;
  using FT = typename LCC::FT;

  alcc.clear();

  Linear_cell_complex_incremental_builder_3<LCC> ib(alcc);

  std::string line, tmp;
  std::size_t npoints, ncells;

  // Skip to POINTS section
  while(std::getline(is, line)) {
    if(line.find("POINTS") != std::string::npos)
      break;
  }
  if(is.eof()) {
    std::cerr << "[ERROR] read_VTK: POINTS section not found" << std::endl;
    return false;
  }

  std::stringstream ss(line);
  std::getline(ss, tmp, ' '); // skip "POINTS"
  ss >> npoints;

  // Read points
  std::vector<typename LCC::Vertex_attribute_descriptor> points(npoints);
  for(std::size_t i = 0; i < npoints; ++i) {
    FT x, y, z;
    if(!(is >> x >> y >> z)) {
      std::cerr << "[ERROR] read_VTK: failed to read point " << i << std::endl;
      return false;
    }
    points[i] = ib.add_vertex(Point(x, y, z));
  }

  // Skip to CELLS section
  while(std::getline(is, line)) {
    if(line.find("CELLS") != std::string::npos)
      break;
  }
  if(is.eof()) {
    std::cerr << "[ERROR] read_VTK: CELLS section not found" << std::endl;
    return false;
  }

  ss = std::stringstream(line);
  std::getline(ss, tmp, ' '); // skip "CELLS"
  ss >> ncells;

  // Read connectivity
  std::vector<std::vector<std::size_t>> faces(ncells);
  std::size_t points_per_cell;
  for(std::size_t i = 0; i < ncells; ++i) {
    if(!(is >> points_per_cell)) {
      std::cerr << "[ERROR] read_VTK: failed to read cell " << i << std::endl;
      return false;
    }
    faces[i].resize(points_per_cell);
    for(std::size_t j = 0; j < points_per_cell; ++j) {
      if(!(is >> faces[i][j])) {
        std::cerr << "[ERROR] read_VTK: failed to read cell " << i << " vertex " << j << std::endl;
        return false;
      }
    }
  }

  // Skip to CELL_TYPES section
  while(std::getline(is, line)) {
    if(line.find("CELL_TYPES") != std::string::npos)
      break;
  }
  if(is.eof()) {
    std::cerr << "[ERROR] read_VTK: CELL_TYPES section not found" << std::endl;
    return false;
  }

  // Create cells based on types
  std::size_t cell_type;
  for(std::size_t i = 0; i < ncells; ++i) {
    if(!(is >> cell_type)) {
      std::cerr << "[ERROR] read_VTK: failed to read cell type " << i << std::endl;
      return false;
    }
    const auto& v = faces[i];
    switch(cell_type) {
    case 10: // TETRA
      if(v.size() == 4)
        make_tetrahedron_with_builder(ib, v[0], v[1], v[2], v[3]);
      break;
    case 11: // VOXEL
      if(v.size() == 8)
        make_hexahedron_with_builder(ib, v[0], v[1], v[3], v[2], v[4], v[5], v[7], v[6]);
      break;
    case 12: // HEXAHEDRON
      if(v.size() == 8)
        make_hexahedron_with_builder(ib, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
      break;
    case 13: // PRISM (WEDGE)
      if(v.size() == 6)
        make_prism_with_builder(ib, v[0], v[1], v[2], v[3], v[4], v[5]);
      break;
    case 14: // PYRAMID
      if(v.size() == 5)
        make_pyramid_with_builder(ib, v[0], v[1], v[2], v[3], v[4]);
      break;
    case 15: // PENTAGONAL_PRISM
      if(v.size() == 10)
        make_pentagonal_prism_with_builder(ib, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9]);
      break;
    case 16: // HEXAGONAL_PRISM
      if(v.size() == 12)
        make_hexagonal_prism_with_builder(ib, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11]);
      break;
    case 42: // GENERIC CELL
      make_generic_cell_with_builder(ib, v);
      break;
    default:
      std::cerr << "[ERROR] read_VTK: type " << cell_type << " unknown." << std::endl;
    }
  }

  // Clean up unused vertex attributes
  for(auto itv = alcc.vertex_attributes().begin(); itv != alcc.vertex_attributes().end();) {
    if(alcc.template dart_of_attribute<0>(itv) == alcc.null_descriptor) {
      alcc.erase_vertex_attribute(itv);
    } else {
      ++itv;
    }
  }

  // Read POINT_DATA scalars if present
  while(std::getline(is, line)) {
    if(line.find("POINT_DATA") != std::string::npos) {
      std::size_t ndata;
      ss = std::stringstream(line);
      std::getline(ss, tmp, ' '); // skip "POINT_DATA"
      ss >> ndata;

      std::getline(is, line); // SCALARS line
      std::stringstream scalars_line(line);
      std::string scalars, name, vtk_type;
      scalars_line >> scalars >> name >> vtk_type;

      std::getline(is, line); // LOOKUP_TABLE line

      if(!read_scalar_by_vtk_type(is, vtk_type, ndata, vertex_reader)) {
        std::cerr << "[ERROR] read_VTK: failed to read POINT_DATA" << std::endl;
      }
      break;
    }
  }

  // Read CELL_DATA scalars if present
  while(std::getline(is, line)) {
    if(line.find("CELL_DATA") != std::string::npos) {
      std::size_t ndata;
      ss = std::stringstream(line);
      std::getline(ss, tmp, ' '); // skip "CELL_DATA"
      ss >> ndata;

      std::getline(is, line); // SCALARS line
      std::stringstream scalars_line(line);
      std::string scalars, name, vtk_type;
      scalars_line >> scalars >> name >> vtk_type;

      std::getline(is, line); // LOOKUP_TABLE line

      if(!read_scalar_by_vtk_type(is, vtk_type, ndata, cell_reader)) {
        std::cerr << "[ERROR] read_VTK: failed to read CELL_DATA" << std::endl;
      }
      break;
    }
  }

  return true;
}

template <typename LCC, typename PointFunctor, typename CellFunctor>
bool write_lcc_to_vtk_ascii(std::ostream& os, const LCC& alcc, PointFunctor ptval, CellFunctor cellval) {
  static_assert(LCC::dimension == 3 && LCC::ambient_dimension == 3,
                "write_VTK() only supports 3D Linear_cell_complexes (3,3)");

  // Write VTK header
  os << "# vtk DataFile Version 2.0\n";
  os << "CGAL Linear_cell_complex\n";
  os << "ASCII\n";
  os << "DATASET UNSTRUCTURED_GRID\n\n";

  // Build vertex index map
  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor, std::size_t> vertex_index;
  std::size_t nbpts = 0;

  for(auto itv = alcc.vertex_attributes().begin(), itvend = alcc.vertex_attributes().end(); itv != itvend; ++itv) {
    vertex_index[itv] = nbpts++;
  }

  // Write points
  os << "POINTS " << nbpts << " double\n";
  for(auto itv = alcc.vertex_attributes().begin(), itvend = alcc.vertex_attributes().end(); itv != itvend; ++itv) {
    const auto& p = itv->point();
    os << p.x() << " " << p.y() << " " << p.z() << "\n";
  }
  os << "\n";

  // Count cells and build connectivity
  std::size_t nbcells = 0;
  std::size_t total_size = 0;
  std::ostringstream cell_stream;
  std::ostringstream type_stream;

  for(auto itvol = alcc.template one_dart_per_cell<3>().begin(), itvolend = alcc.template one_dart_per_cell<3>().end();
      itvol != itvolend; ++itvol)
  {
    ++nbcells;
    VTK_Cell_Type cell_type = get_vtk_cell_type(alcc, itvol);
    type_stream << static_cast<int>(cell_type) << "\n";

    if(cell_type == VTK_POLYHEDRON) {
      // Generic polyhedron format write as face-vertex connectivity
      std::vector<std::vector<std::size_t>> faces;
      std::size_t cell_size = 1; // Start with 1 for number of faces

      for(auto itface = alcc.template one_dart_per_incident_cell<2, 3>(itvol).begin(),
               itfaceend = alcc.template one_dart_per_incident_cell<2, 3>(itvol).end();
          itface != itfaceend; ++itface)
      {
        faces.push_back(std::vector<std::size_t>());
        auto& face = faces.back();

        for(auto itvert = alcc.template darts_of_orbit<1>(itface).begin(),
                 itvertend = alcc.template darts_of_orbit<1>(itface).end();
            itvert != itvertend; ++itvert)
        {
          face.push_back(vertex_index[alcc.vertex_attribute(itvert)]);
        }
        cell_size += face.size() + 1; // +1 for face size
      }

      cell_stream << cell_size << " " << faces.size();
      for(const auto& face : faces) {
        cell_stream << " " << face.size();
        for(auto v : face) {
          cell_stream << " " << v;
        }
      }
      cell_stream << "\n";
      total_size += cell_size + 1; // +1 for cell size

    } else {
      // Standard cell types write vertex connectivity directly
      std::vector<std::size_t> vertices;

      for(auto itvert = alcc.template one_dart_per_incident_cell<0, 3>(itvol).begin(),
               itvertend = alcc.template one_dart_per_incident_cell<0, 3>(itvol).end();
          itvert != itvertend; ++itvert)
      {
        vertices.push_back(vertex_index[alcc.vertex_attribute(itvert)]);
      }

      cell_stream << vertices.size();
      for(auto v : vertices) {
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

  // Write vertex scalars if ptval is not nullptr
  if constexpr(!std::is_same_v<PointFunctor, std::nullptr_t>) {
    os << "POINT_DATA " << nbpts << "\n";
    os << "SCALARS vertex_scalars " << gettype<decltype(ptval(alcc, alcc.null_dart_descriptor))>::name() << " 1\n";
    os << "LOOKUP_TABLE default\n";
    for(auto itv = alcc.vertex_attributes().begin(), itvend = alcc.vertex_attributes().end(); itv != itvend; ++itv) {
      auto dart = alcc.template dart_of_attribute<0>(itv);
      CGAL_assertion(dart != alcc.null_dart_descriptor);
      os << ptval(alcc, dart) << "\n";
    }
  }

  // Write cell scalars if cellval is not nullptr/pointer
  if constexpr(!std::is_same_v<CellFunctor, std::nullptr_t> && !std::is_pointer_v<CellFunctor>) {
    os << "CELL_DATA " << nbcells << "\n";
    os << "SCALARS volume_scalars " << gettype<decltype(cellval(alcc, alcc.null_dart_descriptor))>::name() << " 1\n";
    os << "LOOKUP_TABLE default\n";
    for(auto itvol = alcc.template one_dart_per_cell<3>().begin(),
             itvolend = alcc.template one_dart_per_cell<3>().end();
        itvol != itvolend; ++itvol)
    {
      os << cellval(alcc, itvol) << "\n";
    }
  }

  return true;
}

} // namespace internal

// ============================================================================
//                        Public interface implementation
// ============================================================================

// Functor-based versions
template <typename LCC, typename VertexScalarReader, typename CellScalarReader>
inline bool read_VTK(LCC& alcc, const char* filename, VertexScalarReader vertex_reader, CellScalarReader cell_reader) {
  CGAL_assertion(filename != nullptr);
  std::ifstream file(filename);
  if(!file.is_open()) {
    std::cerr << "[ERROR] read_VTK: cannot open file " << filename << std::endl;
    return false;
  }
  return internal::read_lcc_from_vtk_ascii(file, alcc, vertex_reader, cell_reader);
}

template <typename LCC, typename PointFunctor, typename CellFunctor>
inline bool write_VTK(const LCC& alcc, const char* filename, PointFunctor ptval, CellFunctor cellval) {
  CGAL_assertion(filename != nullptr);
  std::ofstream file(filename);
  if(!file.is_open()) {
    std::cerr << "[ERROR] write_VTK: cannot open file " << filename << std::endl;
    return false;
  }
  return internal::write_lcc_to_vtk_ascii(file, alcc, ptval, cellval);
}

// Vector-based versions (convenience wrappers)
template <typename LCC, typename VertexScalarType, typename VolumeScalarType>
inline bool read_VTK(LCC& alcc,
              const char* filename,
              std::vector<VertexScalarType>* vertex_scalars,
              std::vector<VolumeScalarType>* volume_scalars) {
  auto v_writer = [&](std::size_t i, auto val) {
    if(vertex_scalars) {
      if(vertex_scalars->size() <= i)
        vertex_scalars->resize(i + 1);
      (*vertex_scalars)[i] = static_cast<VertexScalarType>(val);
    }
  };
  auto c_writer = [&](std::size_t i, auto val) {
    if(volume_scalars) {
      if(volume_scalars->size() <= i)
        volume_scalars->resize(i + 1);
      (*volume_scalars)[i] = static_cast<VolumeScalarType>(val);
    }
  };
  return read_VTK(alcc, filename, v_writer, c_writer);
}

template <typename LCC, typename VertexScalarType, typename VolumeScalarType>
inline bool write_VTK(const LCC& alcc,
               const char* filename,
               const std::vector<VertexScalarType>* vertex_scalars,
               const std::vector<VolumeScalarType>* volume_scalars) {
  // Build index maps
  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor, std::size_t> vertex_indices;
  std::size_t idx = 0;
  for(auto itv = alcc.vertex_attributes().begin(), itvend = alcc.vertex_attributes().end(); itv != itvend; ++itv)
    vertex_indices[itv] = idx++;

  std::unordered_map<typename LCC::Dart_const_descriptor, std::size_t> volume_indices;
  idx = 0;
  for(auto itvol = alcc.template one_dart_per_cell<3>().begin(), itvolend = alcc.template one_dart_per_cell<3>().end();
      itvol != itvolend; ++itvol)
    volume_indices[itvol] = idx++;

  return write_VTK(
      alcc, filename,
      [vertex_scalars, &vertex_indices](const LCC& lcc, typename LCC::Dart_const_descriptor d) -> VertexScalarType {
        if(vertex_scalars)
          return (*vertex_scalars)[vertex_indices.at(lcc.template attribute<0>(d))];
        return VertexScalarType();
      },
      [volume_scalars, &volume_indices](const LCC& lcc, typename LCC::Dart_const_descriptor d) -> VolumeScalarType {
        if(volume_scalars)
          return (*volume_scalars)[volume_indices.at(d)];
        return VolumeScalarType();
      });
}

} // namespace IO
} // namespace CGAL

#endif // CGAL_LCC_IO_VTK_H