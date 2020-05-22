// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//                 Mael Rouxel-Labbé

#ifndef CGAL_BGL_IO_VTK_H
#define CGAL_BGL_IO_VTK_H

#include <CGAL/IO/VTK.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#ifdef CGAL_USE_VTK
#include <vtkCell.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#endif

#if defined(CGAL_USE_VTK) || defined(DOXYGEN_RUNNING)

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

namespace IO {
namespace internal {

template <typename FaceGraph, typename NameddParameters>
bool vtkPointSet_to_polygon_mesh(vtkPointSet* poly_data,
                                 FaceGraph& g,
                                 const NameddParameters& np)
{
  typedef typename CGAL::GetVertexPointMap<FaceGraph, NameddParameters>::type      VPM;
  typedef typename boost::property_traits<VPM>::value_type                         Point;
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor               vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor                 face_descriptor;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(CGAL::vertex_point, g));

  vtkIdType nb_points = poly_data->GetNumberOfPoints();
  vtkIdType nb_cells = poly_data->GetNumberOfCells();

  // extract points
  std::vector<vertex_descriptor> vertex_map(nb_points);
  for(vtkIdType i=0; i<nb_points; ++i)
  {
    double coords[3];
    poly_data->GetPoint(i, coords);

    vertex_descriptor v = add_vertex(g);
    put(vpm, v, Point(coords[0], coords[1], coords[2]));
    vertex_map[i] = v;
  }

  // extract cells
  for(vtkIdType i=0; i<nb_cells; ++i)
  {
    int cell_type = poly_data->GetCellType(i);
    if(cell_type != 5
       && cell_type != 7
       && cell_type != 9) // only supported cells are triangles, quads and polygons
      continue;

    vtkCell* cell_ptr = poly_data->GetCell(i);

    vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
    if(nb_vertices < 3)
      return false;

    std::vector<vertex_descriptor> vr(nb_vertices);
    for(vtkIdType k=0; k<nb_vertices; ++k)
    {
      vtkIdType id = cell_ptr->GetPointId(k);
      vr[k] = vertex_map[id];
    }

    face_descriptor f = CGAL::Euler::add_face(vr, g);
    if(f == boost::graph_traits<FaceGraph>::null_face())
      return false;
  }

  return true;
}

} // namespace internal
} // namespace IO

/*!
 * \ingroup PkgBGLIOFct
 *
 * \brief reads a PolyData in the VTP format into a triangulated surface mesh.
 *
 * \tparam FaceGraph a model of `FaceListGraph`.
 *
 * \param fname the path to the file that will be read.
 * \param g the output mesh.
 *
 * \pre \cgal needs to be configured with the VTK Libraries for this function to be available.
 * \cgalNamedParamsBegin
 *  \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
 *    If this parameter is omitted, an internal property map for
 *    `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \pre The data must represent a 2-manifold
 *
 * \see \ref IOStreamVTK
*/
template<typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_VTP(const char* fname, FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  vtkSmartPointer<vtkPointSet> data;
  vtkSmartPointer<IO::internal::ErrorObserverVtk> obs =
    vtkSmartPointer<IO::internal::ErrorObserverVtk>::New();

  data = vtkPolyData::SafeDownCast(IO::internal::read_vtk_file<vtkXMLPolyDataReader>(fname, obs)->GetOutput());

  return IO::internal::vtkPointSet_to_polygon_mesh(data, g, np);
}

template<typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_VTP(const std::string& fname, FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  return read_VTP(fname.c_str(), g, np);
}
template<typename FaceGraph>
bool read_VTP(const char* fname, FaceGraph& g) { return read_VTP(fname, g, parameters::all_default()); }
template<typename FaceGraph>
bool read_VTP(const std::string& fname, FaceGraph& g) { return read_VTP(fname, g, parameters::all_default()); }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

namespace IO {
namespace internal {

// writes the polys appended data at the end of the .vtp file
template <typename FaceGraph, typename NamedParameters>
void write_polys(std::ostream& os,
                 const FaceGraph& g,
                 const NamedParameters& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor           face_descriptor;

  typedef typename CGAL::GetInitializedVertexIndexMap<FaceGraph, NamedParameters>::const_type Vimap;
  Vimap V = CGAL::get_initialized_vertex_index_map(g, np);

  std::vector<std::size_t> connectivity_table;
  std::vector<std::size_t> offsets;
  std::vector<unsigned char> cell_type(num_faces(g), 5);  // triangle == 5

  std::size_t off = 0;

  for(const face_descriptor f : faces(g))
  {
    off += 3;
    offsets.push_back(off);
    for(const vertex_descriptor v : vertices_around_face(halfedge(f, g), g))
      connectivity_table.push_back(get(V, v));
  }

  write_vector<std::size_t>(os, connectivity_table);
  write_vector<std::size_t>(os, offsets);
  write_vector<unsigned char>(os, cell_type);
}

//todo use named params for maps
template <typename FaceGraph, typename NamedParameters>
void write_polys_tag(std::ostream& os,
                     const FaceGraph& g,
                     bool binary,
                     std::size_t& offset,
                     const NamedParameters& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor           face_descriptor;

  typedef typename CGAL::GetInitializedVertexIndexMap<FaceGraph, NamedParameters>::const_type Vimap;
  Vimap V = CGAL::get_initialized_vertex_index_map(g, np);

  std::string formatattribute = binary ? " format=\"appended\"" : " format=\"ascii\"";

  std::string typeattribute;
  switch(sizeof(std::size_t))
  {
    case 8: typeattribute = " type=\"UInt64\""; break;
    case 4: typeattribute = " type=\"UInt32\""; break;
    default: CGAL_error_msg("Unknown size of std::size_t");
  }

  // Write connectivity table
  os << "    <Polys>\n"
     << "      <DataArray Name=\"connectivity\""
     << formatattribute << typeattribute;

  // if binary output, just write the xml tag
  if(binary)
  {
    os << " offset=\"" << offset << "\"/>\n";
    offset += (3 * num_faces(g)+ 1) * sizeof(std::size_t);
    // 3 indices (size_t) per triangle + length of the encoded data (size_t)
  }
  else
  {
    os << ">\n";

    for(face_descriptor f : faces(g))
    {
      for(vertex_descriptor v : vertices_around_face(halfedge(f, g), g))
        os << get(V, v) << " ";
    }
    os << "      </DataArray>\n";
  }

  // Write offsets
  os << "      <DataArray Name=\"offsets\""
     << formatattribute << typeattribute;

  if(binary) {  // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (num_faces(g) + 1) * sizeof(std::size_t);
    // 1 offset (size_t) per triangle + length of the encoded data (size_t)
  }
  else
  {
    os << ">\n";
    std::size_t polys_offset = 0;

    for(face_descriptor f : faces(g))
    {
      polys_offset += 3;
      os << polys_offset << " ";
    }
    os << "      </DataArray>\n";
  }

  // Write cell type (triangle == 5)
  os << "      <DataArray Name=\"types\""
     << formatattribute << " type=\"UInt8\"";

  if(binary)
  {
    os << " offset=\"" << offset << "\"/>\n";
    offset += num_faces(g) + sizeof(std::size_t);
    // 1 unsigned char per cell + length of the encoded data (size_t)
  }
  else
  {
    os << ">\n";
    for(std::size_t i = 0; i< num_faces(g); ++i)
      os << "5 ";
    os << "      </DataArray>\n";
  }
  os << "    </Polys>\n";
}

//todo : use namedparams for points and ids
//overload for facegraph
template <typename FaceGraph, typename NamedParameters>
void write_points_tag(std::ostream& os,
                      const FaceGraph& g,
                      bool binary,
                      std::size_t& offset,
                      const NamedParameters& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor               vertex_descriptor;

  typedef typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::const_type VPM;
  typedef typename boost::property_traits<VPM>::value_type                         Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel                              Gt;
  typedef typename Gt::FT                                                          FT;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  std::string format = binary ? "appended" : "ascii";
  std::string type = (sizeof(FT) == 8) ? "Float64" : "Float32";

  os << "    <Points>\n"
     << "      <DataArray type =\"" << type << "\" NumberOfComponents=\"3\" format=\""
     << format;

  if(binary)
  {
    os << "\" offset=\"" << offset << "\"/>\n";
    offset += 3 * num_vertices(g) * sizeof(FT) + sizeof(std::size_t);
    // 3 coords per points + length of the encoded data (size_t)
  }
  else
  {
    os << "\">\n";
    for(const vertex_descriptor v : vertices(g))
      os << get(vpm, v).x() << " " << get(vpm, v).y() << " " << get(vpm, v).z() << " ";
    os << "      </DataArray>\n";
  }
  os << "    </Points>\n";
}

// writes the points appended data at the end of the .vtp file
template <typename FaceGraph, typename NamedParameters>
void write_polys_points(std::ostream& os,
                        const FaceGraph& g,
                        const NamedParameters& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor               vertex_descriptor;

  typedef typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::const_type VPM;
  typedef typename boost::property_traits<VPM>::value_type                         Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel                              Gt;
  typedef typename Gt::FT                                                          FT;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  std::vector<FT> coordinates;

  for(vertex_descriptor v : vertices(g))
  {
    coordinates.push_back(get(vpm, v).x());
    coordinates.push_back(get(vpm, v).y());
    coordinates.push_back(get(vpm, v).z());
  }

  write_vector<FT>(os, coordinates);
}

} // namespace internal
} // namespace IO

/*! \ingroup PkgBGLIOFct
 *
 * \brief  writes a triangulated surface mesh in the `PolyData` XML format.
 *
 * \tparam FaceGraph a model of `FaceListGraph` with only triangle faces.
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param os the stream used for writing.
 * \param g the triangle mesh to be written.
 * \param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the
 * ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{use_binary_mode} a Boolean indicating if the
 *    data should be written in binary (`true`, the default) or in ASCII (`false`).
 *     \cgalParamEnd
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to
 * the vertices of `g`. If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_Point` must be available in `FaceGraph`.
 *     \cgalParamEnd
 *    \cgalParamBegin{vertex_index_map} the property map with the indices associated to
 * the vertices of `g`.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \see \ref IOStreamVTK
 */
template<typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_VTP(std::ostream& os,
               const FaceGraph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  os << "<?xml version=\"1.0\"?>\n"
     << "<VTKFile type=\"PolyData\" version=\"0.1\"";

#ifdef CGAL_LITTLE_ENDIAN
  os << " byte_order=\"LittleEndian\"";
#else // CGAL_BIG_ENDIAN
  os << " byte_order=\"BigEndian\"";
#endif

  switch(sizeof(std::size_t))
  {
    case 4: os << " header_type=\"UInt32\""; break;
    case 8: os << " header_type=\"UInt64\""; break;
  default: CGAL_error_msg("Unknown size of std::size_t"); return false;
  }

  os << ">\n"
     << "  <PolyData>" << "\n";
  os << "  <Piece NumberOfPoints=\"" << num_vertices(g)
     << "\" NumberOfPolys=\"" << num_faces(g) << "\">\n";

  std::size_t offset = 0;
  const bool binary = choose_parameter(get_parameter(np, internal_np::use_binary_mode), true);

  IO::internal::write_points_tag(os, g, binary, offset, np);
  IO::internal::write_polys_tag(os, g, binary, offset, np);

  os << "   </Piece>\n"
     << "  </PolyData>\n";
  if(binary)
  {
    os << "<AppendedData encoding=\"raw\">\n_";
    IO::internal::write_polys_points(os, g, np);
    IO::internal::write_polys(os, g, np);
  }
  os << "</VTKFile>\n";

  return true;
}

/*! \ingroup PkgBGLIOFct
 *
 * \brief  writes a triangulated surface mesh the file `fname`, in the `PolyData` XML format.
 *
 * \tparam FaceGraph a model of `FaceListGraph` with only triangle faces.
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param os the stream used for writing.
 * \param g the triangle mesh to be written.
 * \param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the
 * ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{use_binary_mode} a Boolean indicating if the
 *    data should be written in binary (`true`, the default) or in ASCII (`false`).
 *     \cgalParamEnd
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to
 * the vertices of `g`. If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_Point` must be available in `FaceGraph`.
 *     \cgalParamEnd
 *    \cgalParamBegin{vertex_index_map} the property map with the indices associated to
 * the vertices of `g`.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \see \ref IOStreamVTK
 */
template<typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_VTP(const char* fname, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream os(fname);
  return write_VTP(os, g, np);
}

template<typename FaceGraph>
bool write_VTP(std::ostream& os, const FaceGraph& g) { return write_VTP(os, g, CGAL::parameters::all_default()); }
template<typename FaceGraph>
bool write_VTP(const char* fname, const FaceGraph& g) { return write_VTP(fname, g, parameters::all_default()); }
template<typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_VTP(const std::string& fname, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np) {
  return write_VTP(fname.c_str(), g, np);
}
template<typename FaceGraph>
bool write_VTP(const std::string& fname, const FaceGraph& g) {
  return write_VTP(fname, g, parameters::all_default());
}

} // namespace CGAL

#endif // defined(CGAL_USE_VTK) || defined(DOXYGEN_RUNNING)

#endif // CGAL_BGL_IO_VTK_H

