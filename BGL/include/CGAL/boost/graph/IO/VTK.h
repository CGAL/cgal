// Copyright (c) 2015-2020  GeometryFactory (France).  All rights reserved.
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

#include <fstream>
#include <string>
#include <vector>

#ifdef CGAL_USE_VTK
#include <vtkCell.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#endif

#if defined(CGAL_USE_VTK) || defined(DOXYGEN_RUNNING)

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#define CGAL_DEPRECATED
#endif

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {
namespace internal {

template <typename Graph, typename NameddParameters>
bool vtkPointSet_to_polygon_mesh(vtkPointSet* poly_data,
                                 Graph& g,
                                 const NameddParameters& np)
{
  typedef typename CGAL::GetVertexPointMap<Graph, NameddParameters>::type      VPM;
  typedef typename boost::property_traits<VPM>::value_type                     Point;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor               vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor                 face_descriptor;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(CGAL::vertex_point, g));

  vtkIdType nb_points = poly_data->GetNumberOfPoints();
  vtkIdType nb_cells = poly_data->GetNumberOfCells();

  reserve(g, num_vertices(g) + static_cast<typename boost::graph_traits<Graph>::vertices_size_type>(nb_points),
             num_edges(g) + static_cast<typename boost::graph_traits<Graph>::edges_size_type>(nb_points + nb_cells),
             num_faces(g) + static_cast<typename boost::graph_traits<Graph>::faces_size_type>(nb_cells));

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
    if(f == boost::graph_traits<Graph>::null_face())
      return false;
  }

  return true;
}

} // namespace internal

/*!
 * \ingroup PkgBGLIoFuncsVTP
 *
 * \brief reads a PolyData in the \ref IOStreamVTK into a triangulated surface mesh.
 *
 * The data is expected to represent a 2-manifold (possibly with borders).
 *
 * \attention The graph `g` is not cleared, and the data from the file are appended.
 *
 * \tparam Graph a model of `MutableFaceGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the path to the file that will be read
 * \param g the output mesh
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `g`}
 *     \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `Graph`.}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 *
 * \returns `true` if reading was successful, `false` otherwise.
*/
template<typename Graph,
         typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_VTP(const std::string& fname,
              Graph& g,
              const CGAL_BGL_NP_CLASS& np)
{
  std::ifstream test(fname);
  if(!test.good())
  {
    std::cerr<<"File doesn't exist."<<std::endl;
    return false;
  }
  test.close();
  vtkSmartPointer<vtkPointSet> data;
  vtkSmartPointer<internal::ErrorObserverVtk> obs =
    vtkSmartPointer<internal::ErrorObserverVtk>::New();

  data = vtkPolyData::SafeDownCast(internal::read_vtk_file<vtkXMLPolyDataReader>(fname, obs)->GetOutput());
  if (obs->GetError())
    return false;
  return internal::vtkPointSet_to_polygon_mesh(data, g, np);
}

/// \cond SKIP_IN_MANUAL

template<typename Graph>
bool read_VTP(const std::string& fname, Graph& g) { return read_VTP(fname, g, parameters::all_default()); }

/// \endcond

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace internal {

// writes the polys appended data at the end of the .vtp file
template <typename Graph, typename NamedParameters>
void write_polys(std::ostream& os,
                 const Graph& g,
                 const NamedParameters& np)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor           face_descriptor;

  typedef typename CGAL::GetInitializedVertexIndexMap<Graph, NamedParameters>::const_type Vimap;
  Vimap V = CGAL::get_initialized_vertex_index_map(g, np);

  std::vector<std::size_t> connectivity_table;
  std::vector<std::size_t> offsets;
  std::vector<unsigned char> cell_type(faces(g).size(), 5);  // triangle == 5

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
template <typename Graph, typename NamedParameters>
void write_polys_tag(std::ostream& os,
                     const Graph& g,
                     bool binary,
                     std::size_t& offset,
                     const NamedParameters& np)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor           face_descriptor;

  typedef typename CGAL::GetInitializedVertexIndexMap<Graph, NamedParameters>::const_type Vimap;
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
    offset += (3 * faces(g).size() + 1) * sizeof(std::size_t);
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
    offset += (faces(g).size() + 1) * sizeof(std::size_t);
    // 1 offset (size_t) per triangle + length of the encoded data (size_t)
  }
  else
  {
    os << ">\n";
    std::size_t polys_offset = 0;

    for(std::size_t i = 0; i< faces(g).size(); ++i)
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
    offset += faces(g).size() + sizeof(std::size_t);
    // 1 unsigned char per cell + length of the encoded data (size_t)
  }
  else
  {
    os << ">\n";
    for(std::size_t i=0; i<faces(g).size(); ++i)
      os << "5 ";
    os << "      </DataArray>\n";
  }
  os << "    </Polys>\n";
}

//todo : use namedparams for points and ids
//overload for facegraph
template <typename Graph, typename NamedParameters>
void write_points_tag(std::ostream& os,
                      const Graph& g,
                      bool binary,
                      std::size_t& offset,
                      const NamedParameters& np)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor               vertex_descriptor;

  typedef typename CGAL::GetVertexPointMap<Graph, NamedParameters>::const_type     VPM;
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
    offset += 3 * vertices(g).size() * sizeof(FT) + sizeof(std::size_t);
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
template <typename Graph, typename NamedParameters>
void write_polys_points(std::ostream& os,
                        const Graph& g,
                        const NamedParameters& np)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor                   vertex_descriptor;

  typedef typename CGAL::GetVertexPointMap<Graph, NamedParameters>::const_type     VPM;
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

/*! \ingroup PkgBGLIoFuncsVTP
 *
 * \brief writes a triangulated surface mesh in the `PolyData` XML format (\ref IOStreamVTK).
 *
 * \tparam Graph a model of `FaceListGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param os the output stream
 * \param g the triangle mesh to be output
 * \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{use_binary_mode}
 *     \cgalParamDescription{indicates whether data should be written in binary (`true`) or in ASCII (`false`)}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `g`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `Graph`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_index_map}
 *     \cgalParamDescription{a property map associating to each vertex of `graph` a unique index between `0` and `num_vertices(graph) - 1`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{no vertex indices in the output}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{`the precision of the stream `os``}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre `g` contains only triangular faces
 *
 * \returns `true` if writing was successful, `false` otherwise.
 */
template<typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_VTP(std::ostream& os,
               const Graph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  if(!os.good())
    return false;

  set_stream_precision_from_NP(os, np);

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
  os << "  <Piece NumberOfPoints=\"" << vertices(g).size()
     << "\" NumberOfPolys=\"" << faces(g).size() << "\">\n";

  std::size_t offset = 0;
  const bool binary = choose_parameter(get_parameter(np, internal_np::use_binary_mode), true);

  internal::write_points_tag(os, g, binary, offset, np);
  internal::write_polys_tag(os, g, binary, offset, np);

  os << "   </Piece>\n"
     << "  </PolyData>\n";
  if(binary)
  {
    os << "<AppendedData encoding=\"raw\">\n_";
    internal::write_polys_points(os, g, np);
    internal::write_polys(os, g, np);
  }
  os << "</VTKFile>" << std::endl;

  return os.good();
}

/*! \ingroup PkgBGLIoFuncsVTP
 *
 * \brief writes a triangulated surface mesh the file `fname`, in the `PolyData` XML format (\ref IOStreamVTK).
 *
 * \tparam Graph a model of `FaceListGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the name of the output file
 * \param g the triangle mesh to be output
 * \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the
 * ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{use_binary_mode}
 *     \cgalParamDescription{indicates whether data should be written in binary (`true`) or in ASCII (`false`)}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `g`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `Graph`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_index_map}
 *     \cgalParamDescription{a property map associating to each vertex of `graph` a unique index between `0` and `num_vertices(graph) - 1`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{no vertex indices in the output}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{`6`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre `g` contains only triangular faces
 *
 * \returns `true` if writing was successful, `false` otherwise.
 */
template<typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_VTP(const std::string& fname, const Graph& g, const CGAL_BGL_NP_CLASS& np)
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  std::ofstream os;
  if(binary){
    os.open(fname, std::ios::binary);
    CGAL::IO::set_mode(os, CGAL::IO::BINARY);
  }
  else
    os.open(fname);

  return write_VTP(os, g, np);
}

/// \cond SKIP_IN_MANUAL

template<typename Graph>
bool write_VTP(std::ostream& os, const Graph& g) { return write_VTP(os, g, CGAL::parameters::all_default()); }
template<typename Graph>
bool write_VTP(const std::string& fname, const Graph& g) { return write_VTP(fname, g, parameters::all_default()); }

/// \endcond

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
 \ingroup PkgBGLIOFctDeprecated

 \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::write_VTP()` should be used instead.
*/
template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool write_vtp(std::ostream& os, const Graph& g, const CGAL_BGL_NP_CLASS& np)
{
  return IO::write_VTP(os, g, np);
}

template <typename Graph>
CGAL_DEPRECATED bool write_vtp(std::ostream& os, const Graph& g)
{
  return write_vtp(os, g, parameters::all_default());
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // defined(CGAL_USE_VTK) || defined(DOXYGEN_RUNNING)

#endif // CGAL_BGL_IO_VTK_H
