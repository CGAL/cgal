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
#define CGAL_LCC_IO_VTK_H

#include <CGAL/Linear_cell_complex_incremental_builder_3.h>
#include <CGAL/assertions.h>
#include <CGAL/Element_topo.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace CGAL {
namespace IO {

/*
 * Functions to import/export 3D Linear_cell_complex from/to VTK legacy ASCII
 * format.
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

/*
 * Read a VTK legacy ASCII file and load it into a 3D Linear_cell_complex.
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
template <typename LCC, typename VertexScalarType,
          typename VolumeScalarType>
bool read_VTK(const char* filename,
              LCC& alcc,
              std::vector<VertexScalarType>* vertex_scalars,
              std::vector<VolumeScalarType>* volume_scalars);

/*
 * Write a 3D Linear_cell_complex to a VTK legacy ASCII file.
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
template <typename LCC, typename VertexScalarType,
          typename VolumeScalarType>
bool write_VTK(const char* filename,
               const LCC& alcc,
               const std::vector<VertexScalarType>* vertex_scalars,
               const std::vector<VolumeScalarType>* volume_scalars);

// "Advanced" versions with functors
template <typename LCC, typename PointFunctor, typename CellFunctor>
bool write_VTK_with_fct(const char* filename, const LCC& alcc,
                        PointFunctor ptval, CellFunctor cellval);

// ============================================================================
//                          Implementation details
// ============================================================================

namespace internal
{
  /////////////////////////////////////////////////////////////////////////////
  // VTK type name mapping
  // bit, unsigned_char, char, unsigned_short, short, unsigned_int, int,
  // unsigned_long, long, float, double.
  template<typename T>
  struct gettype
  { static std::string name() { return "unknown"; }};
  template<>
  struct gettype<bool>
  { static std::string name() { return "bit"; }};
  template<>
  struct gettype<unsigned char>
  { static std::string name() { return "unsigned_char"; }};
  template<>
  struct gettype<char>
  { static std::string name() { return "char"; }};
  template<>
  struct gettype<unsigned short int>
  { static std::string name() { return "unsigned_short"; }};
  template<>
  struct gettype<short int>
  { static std::string name() { return "short"; }};
  template<>
  struct gettype<unsigned int>
  { static std::string name() { return "unsigned_int"; }};
  template<>
  struct gettype<int>
  { static std::string name() { return "int"; }};
  template<>
  struct gettype<unsigned long int>
  { static std::string name() { return "unsigned_long"; }};
  template<>
  struct gettype<long int>
  { static std::string name() { return "long"; }};
  template<>
  struct gettype<float>
  { static std::string name() { return "float"; }};
  template<>
  struct gettype<double>
  { static std::string name() { return "double"; }};
  /////////////////////////////////////////////////////////////////////////////
  // VTK cell type constants
  enum VTK_Cell_Type
  {
    VTK_TETRA = 10,
    VTK_VOXEL = 11,
    VTK_HEXAHEDRON = 12,
    VTK_WEDGE = 13, // Prism
    VTK_PYRAMID = 14,
    VTK_PENTAGONAL_PRISM = 15,
    VTK_HEXAGONAL_PRISM = 16,
    VTK_POLYHEDRON = 42 // Generic cell
  };
  /////////////////////////////////////////////////////////////////////////////
  /// Write cell_data.
  template<typename FCT>
  struct Write_cell_data
  {
    /// nb is the number of cells,
    /// fct is a function having 3 parameters: a lcc,  a dart_descriptor,
    ///    an the index of the cell.
    template<typename LCC>
    static void run(std::ofstream& fo, LCC& lcc, std::size_t nb, FCT fct)
    {
      fo<<"CELL_DATA "<<nb<<std::endl;
      fo<<"SCALARS cell_scalars "
         <<gettype<decltype(fct(lcc, lcc.null_dart_descriptor, 0))>::name()
         <<" 1"<<std::endl;
      fo<<"LOOKUP_TABLE default"<<std::endl;
      std::size_t i=0;
      for(auto itvol=lcc.template one_dart_per_cell<3>().begin(),
           itvolend=lcc.template one_dart_per_cell<3>().end();
           itvol!=itvolend; ++itvol, ++i)
      { fo<<fct(lcc, itvol, i)<<std::endl; }
      fo<<std::endl;
    }
  };
  template<>
  struct Write_cell_data<std::nullptr_t>
  {
    template<typename LCC>
    static void run(std::ofstream&, LCC&, std::size_t, std::nullptr_t)
    {}
  };
  /////////////////////////////////////////////////////////////////////////////
  /// Write point_data.
  template<typename FCT>
  struct Write_point_data
  {
    /// nb is the number of cells,
    /// fct is a function having 3 parameters: a lcc,  a dart_descriptor,
    ///    an the index of the cell.
    template<typename LCC>
    static void run(std::ofstream& fo, LCC& lcc, std::size_t nb, FCT fct)
    {
      fo<<"POINT_DATA "<<nb<<std::endl;
      fo<<"SCALARS point_scalars "
         <<gettype<decltype(fct(lcc, lcc.null_dart_descriptor, 0))>::name()
         <<" 1"<<std::endl;
      fo<<"LOOKUP_TABLE default"<<std::endl;
      std::size_t i=0;
      for(auto itv=lcc.vertex_attributes().begin(),
           itvend=lcc.vertex_attributes().end(); itv!=itvend; ++itv, ++i)
      { fo<<fct(lcc, lcc.template dart_of_attribute<0>(itv), i)<<std::endl; }
      fo<<std::endl;
    }
  };
  /////////////////////////////////////////////////////////////////////////////
  template<>
  struct Write_point_data<std::nullptr_t>
  {
    template<typename LCC>
    static void run(std::ofstream&, LCC&, std::size_t, std::nullptr_t)
    {}
  };
  /////////////////////////////////////////////////////////////////////////////
  // Read data, stored values as T.
  template<typename T>
  bool read_data(std::istream& fi, std::string& line, std::vector<T>& data)
  {
    std::string txt, data_type;
    std::size_t nb;
    std::istringstream inputline(line);
    inputline>>txt>>nb;  // "CELL_DATA xxx"
    fi>>txt>>txt; // "SCALARS cell_scalars "
    fi>>data_type>>txt; // type for data
    fi>>txt>>txt; // "LOOKUP_TABLE default"
    if(!fi.good())
    { return false; }
    data.clear();
    data.reserve(nb);
    for(std::size_t i=0; i<nb; ++i)
    {
      if(!(fi>>txt))
      { return false; }

      std::stringstream ss{txt};
      T t;
      ss>>t;

      data.push_back(t);
    }
    return true;
  }
  /////////////////////////////////////////////////////////////////////////////
  // Helper: detect VTK cell type from a 3-cell
  template<typename LCC>
  VTK_Cell_Type get_vtk_cell_type(const LCC& lcc,
                                  typename LCC::Dart_const_descriptor itvol,
                                  typename LCC::Dart_const_descriptor& sd)
  {
    using namespace CGAL::CMap::Element_topo;
    cell_topo vol_type=get_cell_topo<3>(lcc, itvol, sd);
    switch(vol_type)
    {
      case TETRAHEDRON: return VTK_TETRA;
      case PYRAMID: return VTK_PYRAMID;
      case PRISM: return VTK_WEDGE;
      case HEXAHEDRON: return VTK_HEXAHEDRON;
      // case PENTAGONAL_PRISM: return VTK_PENTAGONAL_PRISM;
      // case HEXAGONAL_PRISM: return VTK_HEXAGONAL_PRISM;
        //       24 QUADRATIC_TETRA
        //       25 QUADRATIC_HEXAHEDRON
        //       26 QUADRATIC_WEDGE
        //       27 QUADRATIC_PYRAMID
      default: break;
    }
    return VTK_POLYHEDRON;
  }
  /////////////////////////////////////////////////////////////////////////////
template <typename LCC, typename VertexScalarType=float,
          typename CellScalarType=float>
bool read_lcc_from_vtk_ascii(std::istream& is, LCC& alcc,
                             std::vector<VertexScalarType>* vertex_scalars=nullptr,
                             std::vector<CellScalarType>* cell_scalars=nullptr)
{
  static_assert(LCC::dimension==3 && LCC::ambient_dimension==3,
                "read_VTK() only supports 3D Linear_cell_complexes (3,3)");

  using Point=typename LCC::Point;
  using FT=typename LCC::FT;

  Linear_cell_complex_incremental_builder_3<LCC> ib(alcc);

  std::string line, tmp;
  std::size_t npoints, ncells;

  // Skip to POINTS section
  while(std::getline(is, line) && line.find("POINTS")==std::string::npos)
  {}
  if(is.eof())
  {
    std::cerr<<"[ERROR] read_VTK: POINTS section not found"<<std::endl;
    return false;
  }

  std::stringstream ss(line);
  std::getline(ss, tmp, ' '); // skip "POINTS"
  ss>>npoints;

  // Read points
  std::vector<typename LCC::Vertex_attribute_descriptor> points(npoints);
  for(std::size_t i=0; i<npoints; ++i)
  {
    FT x, y, z;
    if(!(is>>x>>y>>z))
    {
      std::cerr<<"[ERROR] read_VTK: failed to read point "<<i<<std::endl;
      return false;
    }
    points[i]=ib.add_vertex(Point(x, y, z));
  }

  // Skip to CELLS section
  while(std::getline(is, line) && line.find("CELLS")==std::string::npos)
  {}
  if(is.eof())
  {
    std::cerr<<"[ERROR] read_VTK: CELLS section not found"<<std::endl;
    return false;
  }

  ss=std::stringstream(line);
  std::getline(ss, tmp, ' '); // skip "CELLS"
  ss>>ncells;

  // Read connectivity
  std::vector<std::vector<std::size_t>> faces(ncells);
  std::size_t points_per_cell;
  for(std::size_t i=0; i<ncells; ++i)
  {
    if(!(is>>points_per_cell))
    {
      std::cerr<<"[ERROR] read_VTK: failed to read cell "<<i<<std::endl;
      return false;
    }
    faces[i].resize(points_per_cell);
    for(std::size_t j=0; j<points_per_cell; ++j)
    {
      if(!(is>>faces[i][j]))
      {
        std::cerr<<"[ERROR] read_VTK: failed to read cell "<<i<<" vertex "<<j<< std::endl;
        return false;
      }
    }
  }

  // Skip to CELL_TYPES section
  while(std::getline(is, line) && line.find("CELL_TYPES")==std::string::npos)
  {}
  if(is.eof())
  {
    std::cerr<<"[ERROR] read_VTK: CELL_TYPES section not found"<<std::endl;
    return false;
  }

  // Create cells based on types
  std::size_t cell_type;
  for(std::size_t i = 0; i<ncells; ++i)
  {
    if(!(is>>cell_type))
    {
      std::cerr<<"[ERROR] read_VTK: failed to read cell type "<<i<< std::endl;
      return false;
    }
    const auto& v=faces[i];
    switch(cell_type)
    {
    case VTK_TETRA:
      if(v.size()==4)
      { make_tetrahedron_with_builder(ib, v[0], v[1], v[2], v[3]); }
      break;
    case VTK_VOXEL:
      if(v.size()==8)
      { make_hexahedron_with_builder(ib, v[0], v[1], v[3], v[2], v[4], v[5],
                                     v[7], v[6]); }
      break;
    case VTK_HEXAHEDRON:
      if(v.size()==8)
      { make_hexahedron_with_builder(ib, v[0], v[1], v[2], v[3], v[4], v[5],
                                     v[6], v[7]); }
      break;
    case VTK_WEDGE: // PRISM
      if(v.size()==6)
      { make_prism_with_builder(ib, v[0], v[1], v[2], v[3], v[4], v[5]); }
      break;
    case VTK_PYRAMID:
      if(v.size()==5)
      { make_pyramid_with_builder(ib, v[0], v[1], v[2], v[3], v[4]); }
      break;
    case VTK_PENTAGONAL_PRISM:
      if(v.size()==10)
      { make_pentagonal_prism_with_builder(ib, v[0], v[1], v[2], v[3], v[4],
                                           v[5], v[6], v[7], v[8], v[9]); }
      break;
    case VTK_HEXAGONAL_PRISM:
      if(v.size()==12)
      { make_hexagonal_prism_with_builder(ib, v[0], v[1], v[2], v[3], v[4],
                                          v[5], v[6], v[7], v[8], v[9],
                                          v[10], v[11]); }
      break;
    case VTK_POLYHEDRON: // GENERIC CELL
      make_generic_cell_with_builder(ib, v);
      break;
    default:
      std::cerr<<"[ERROR] read_VTK: type "<<cell_type<<" unknown."<<std::endl;
    }
  }

  // Clean up unused vertex attributes
  for(auto itv=alcc.vertex_attributes().begin();
       itv!=alcc.vertex_attributes().end(); ++itv)
  {
    if(alcc.template dart_of_attribute<0>(itv)==alcc.null_descriptor)
    { alcc.erase_vertex_attribute(itv); }
  }

  if(vertex_scalars!=nullptr)
  { vertex_scalars->clear(); }
  if(cell_scalars!=nullptr)
  { cell_scalars->clear(); }

  while(std::getline(is, line))
  {
    // Read POINT_DATA scalars if present
    if(vertex_scalars!=nullptr && line.find("POINT_DATA")!=std::string::npos)
    {
      if(!read_data(is, line, *vertex_scalars))
      {
        std::cerr<<"[ERROR] read_VTK: error when reading POINT_DATA."
                 <<std::endl;
      }
    }
    // Read CELL_DATA scalars if present
    else if(cell_scalars!=nullptr && line.find("CELL_DATA")!=std::string::npos)
    {
      if(!read_data(is, line, *cell_scalars))
      {
        std::cerr<<"[ERROR] read_VTK: error when reading CELL_DATA."
                  <<std::endl;
      }
    }
  }
  return true;
}
/////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool write_lcc_topo_to_vtk_ascii(std::ostream& os, const LCC& alcc,
                                 std::size_t& nbpts, std::size_t& nbcells)
{
  static_assert(LCC::dimension==3 && LCC::ambient_dimension==3,
                "write_VTK() only supports 3D Linear_cell_complexes (3,3)");

  // Write VTK header
  os<<"# vtk DataFile Version 2.0\n";
  os<<"CGAL Linear_cell_complex\n";
  os<<"ASCII\n";
  os<<"DATASET UNSTRUCTURED_GRID\n\n";

  // Build vertex index map and write points
  std::unordered_map<typename LCC::Vertex_attribute_const_descriptor, std::size_t>
      index;
  nbpts=0;
  os<<"POINTS "<<alcc.vertex_attributes().size()<<" double"<<std::endl;
  for(auto itv=alcc.vertex_attributes().begin(),
       itvend=alcc.vertex_attributes().end(); itv!=itvend; ++itv)
  {
    os<<" "<<itv->point()<<std::endl;
    index[itv]=nbpts++;
  }
  os<<std::endl;

  // Count cells and build connectivity
  nbcells=0;
  std::size_t total_size=0;
  std::ostringstream cell_stream, type_stream;
  typename LCC::Dart_const_descriptor sd;

  // Write cells section
  for(typename LCC::template One_dart_per_cell_range<3>::const_iterator
           itvol=alcc.template one_dart_per_cell<3>().begin(),
       itvolend=alcc.template one_dart_per_cell<3>().end();
      itvol!=itvolend; ++itvol)
  {
    ++nbcells;
    ++total_size; // for the number of vertices
    VTK_Cell_Type cell_type=get_vtk_cell_type(alcc, itvol, sd);
    type_stream<<static_cast<int>(cell_type)<<std::endl;

    if(cell_type==VTK_TETRA)
    {
      cell_stream<<" 4 "
         <<index[alcc.vertex_attribute(sd)]<<" "
         <<index[alcc.vertex_attribute(alcc.template beta<1>(sd))]<<" "
         <<index[alcc.vertex_attribute(alcc.template beta<0>(sd))]<<" "
         <<index[alcc.vertex_attribute(alcc.template beta<2, 0>(sd))]<<std::endl;
      total_size+=4;
    }
    else if(cell_type==VTK_PYRAMID)
    {
      cell_stream<<" 5 "
         <<index[alcc.vertex_attribute(sd)]<<" "
         <<index[alcc.vertex_attribute(alcc.template beta<1>(sd))]<<" "
         <<index[alcc.vertex_attribute(alcc.template beta<1,1>(sd))]<<" "
         <<index[alcc.vertex_attribute(alcc.template beta<0>(sd))]<<" "
         <<index[alcc.vertex_attribute(alcc.template beta<2,0>(sd))]<<std::endl;
      total_size+=5;
    }
    else if(cell_type==VTK_WEDGE)
    {
      cell_stream<<" 6 "
         <<index[alcc.vertex_attribute(sd)]<<" "
         <<index[alcc.vertex_attribute(alcc.template beta<1>(sd))]<<" "
         <<index[alcc.vertex_attribute(alcc.template beta<0>(sd))]<<" ";
      // Move to the up face
      typename LCC::Dart_const_descriptor d2=alcc.template beta<2, 1, 1, 2>(sd);
      cell_stream<<index[alcc.vertex_attribute(alcc.template beta<1>(d2))]<<" "
         <<index[alcc.vertex_attribute(d2)]<<" "
         <<index[alcc.vertex_attribute(alcc.template beta<0>(d2))]<<std::endl;
      total_size+=6;
    }
    else if(cell_type==VTK_HEXAHEDRON)
    {
      cell_stream<<" 8 ";
      for(unsigned int i=0; i<4; ++i)
      {
        cell_stream<<index[alcc.vertex_attribute(sd)]<<" ";
        sd=alcc.template beta<1>(sd);
      }
      typename LCC::Dart_const_descriptor d2=alcc.template beta<2, 1, 1, 2, 1>(sd);
      // Darts associated with particles 4, 5, 6, 7
      for(unsigned int i = 0; i < 4; i++)
      {
        cell_stream<<index[alcc.vertex_attribute(d2)]<<" ";
        d2 = alcc.template beta<0>(d2);
      }
      cell_stream<<std::endl;
      total_size+=8;
    }
    // TODO: 15 PENTAGONAL_PRISM
    //       16 HEXAGONAL_PRISM
    //       24 QUADRATIC_TETRA
    //       25 QUADRATIC_HEXAHEDRON
    //       26 QUADRATIC_WEDGE
    //       27 QUADRATIC_PYRAMID
    else
    {
      // Generic polyhedron format write as face-vertex connectivity
      std::vector<std::vector<std::size_t>> faces;
      std::size_t cell_size=1; // Start with 1 for number of faces
      ++total_size; // for the same reason
      for(auto itface=alcc.template one_dart_per_incident_cell<2, 3, 2>(itvol).begin(),
               itfaceend=alcc.template one_dart_per_incident_cell<2, 3, 2>(itvol).end();
          itface!=itfaceend; ++itface)
      {
        faces.push_back(std::vector<std::size_t>());
        typename LCC::Dart_const_descriptor curdh=itface;
        do
        {
          faces.back().push_back(index[alcc.vertex_attribute(curdh)]);
          curdh=alcc.template beta<1>(curdh);
        }
        while(curdh!=itface);
        cell_size+=faces.back().size()+1; // +1 for the number of vertices in the face
      }
      cell_stream<<cell_size<<" "<<faces.size();
      for(const auto& face : faces)
      {
        cell_stream<<" "<<face.size();
        for(auto v : face)
        { cell_stream<<" "<<v; }
        total_size+=face.size()+1; // +1 for the number of vertices in the face
      }
      cell_stream<<std::endl;
    }
  }
  os<<"CELLS "<<nbcells<<" "<<total_size<<std::endl;
  os<<cell_stream.str()<<std::endl;

  // Write cell types
  os<<"CELL_TYPES "<<nbcells<<std::endl;
  os<<type_stream.str()<<std::endl;

  return true;
}

} // namespace internal

// ============================================================================
//                        Public interface implementation
// ============================================================================

////////////////////////////////////////////////////////////////////////////////////
template <typename LCC, typename VertexScalarType, typename VolumeScalarType>
bool read_VTK(const char* filename, LCC& alcc,
              std::vector<VertexScalarType>* vertex_scalars,
              std::vector<VolumeScalarType>* volume_scalars)
{
  CGAL_assertion(filename!=nullptr);
  std::ifstream file(filename);
  if(!file.is_open())
  {
    std::cerr<<"[ERROR] read_VTK: cannot open file "<<filename<<std::endl;
    return false;
  }
  return internal::read_lcc_from_vtk_ascii(file, alcc,
                                           vertex_scalars, volume_scalars);
}

template <typename LCC>
bool read_VTK(const char* filename, LCC& alcc)
{ return read_VTK<LCC, float, float>(filename, alcc, nullptr, nullptr); }

template <typename LCC, typename VertexScalarType>
bool read_VTK(const char* filename, LCC& alcc,
              std::vector<VertexScalarType>* vertex_scalars)
{ return read_VTK<LCC, VertexScalarType, float>
      (filename, alcc, vertex_scalars, nullptr); }

template <typename LCC, typename VolumeScalarType>
bool read_VTK(const char* filename, LCC& alcc,
              std::nullptr_t,
              std::vector<VolumeScalarType>* volume_scalars)
{ return read_VTK<LCC, float, VolumeScalarType>
      (filename, alcc, nullptr, volume_scalars); }
////////////////////////////////////////////////////////////////////////////////////
template <typename LCC, typename PointFunctor, typename CellFunctor>
inline bool write_VTK_with_fct(const char* filename, const LCC& alcc,
                               PointFunctor pointfct, CellFunctor cellfct)
{
  CGAL_assertion(filename!=nullptr);
  std::ofstream file(filename);
  if(!file.good())
  {
    std::cerr<<"[ERROR] write_VTK: cannot open file "<<filename<<std::endl;
    return false;
  }
  std::size_t nbpts=0, nbcells=0;
  bool res=internal::write_lcc_topo_to_vtk_ascii(file, alcc, nbpts, nbcells);
  if(res)
  {
    if(pointfct)
    { internal::Write_point_data<PointFunctor>::
          run(file, alcc, nbpts, pointfct); }
    if(cellfct)
    { internal::Write_cell_data<CellFunctor>::
          run(file, alcc, nbcells, cellfct); }
  }
  file.close();
  return true;
}
////////////////////////////////////////////////////////////////////////////////////
template <typename LCC, typename VertexScalarType, typename VolumeScalarType>
bool write_VTK(const char* filename, const LCC& alcc,
               const std::vector<VertexScalarType>* vertex_scalars,
               const std::vector<VolumeScalarType>* volume_scalars)
{
  std::function<VertexScalarType(const LCC&,
                                 typename LCC::Dart_const_descriptor,
                                 std::size_t i)> vertexfct;
  std::function<VolumeScalarType(const LCC&,
                                 typename LCC::Dart_const_descriptor,
                                 std::size_t i)> cellfct;
  if(vertex_scalars!=nullptr)
  {
    vertexfct=[&vertex_scalars](const LCC&, typename LCC::Dart_const_descriptor,
                               std::size_t i) -> VertexScalarType
    { return (*vertex_scalars)[i]; };
  }

  if(volume_scalars!=nullptr)
  {
    cellfct=[&volume_scalars](const LCC&, typename LCC::Dart_const_descriptor,
                             std::size_t i) -> VolumeScalarType
    { return (*volume_scalars)[i]; };
  }

  return write_VTK_with_fct(filename, alcc, vertexfct, cellfct);
}

template <typename LCC>
bool write_VTK(const char* filename, const LCC& alcc)
{
  return write_VTK<LCC, float, float>(filename, alcc, nullptr, nullptr);
}

template <typename LCC, typename VertexScalarType>
bool write_VTK(const char* filename, const LCC& alcc,
               const std::vector<VertexScalarType>* vertex_scalars)
{
  return write_VTK<LCC, VertexScalarType, float>(filename, alcc, vertex_scalars,
                                                 nullptr);
}

template <typename LCC, typename VolumeScalarType>
bool write_VTK(const char* filename, const LCC& alcc,
               std::nullptr_t,
               const std::vector<VolumeScalarType>* volume_scalars)
{
  return write_VTK<LCC, float, VolumeScalarType>(filename, alcc, nullptr,
                                                 volume_scalars);
}
////////////////////////////////////////////////////////////////////////////////////

} // namespace IO
} // namespace CGAL

#endif // CGAL_LCC_IO_VTK_H
