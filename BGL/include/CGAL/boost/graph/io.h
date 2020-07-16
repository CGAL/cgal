// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_IO_H
#define CGAL_BOOST_GRAPH_IO_H

#include <boost/container/flat_map.hpp>

#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/IO/write_vtk.h>

namespace CGAL {

/*!
 \ingroup PkgBGLIOFct

  writes the graph `g` in the wrl format (VRML 2.0).

  \param os the output stream
  \param g the graph to be written
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `FaceGraph`.}
    \cgalParamNEnd
  \cgalNamedParamsEnd
*/
template <typename FaceGraph, typename NamedParameters>
bool write_wrl(std::ostream& os,
               const FaceGraph& g,
               const NamedParameters& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  boost::container::flat_map<vertex_descriptor,vertices_size_type> reindex;
  int n = 0;

  os << "#VRML V2.0 utf8\n"
    "Group {\n"
    "children [\n"
    "Shape {\n"
    "appearance DEF A1 Appearance {\n"
    "material Material {\n"
    "diffuseColor .6 .5 .9\n"
    "}\n"
    "}\n"
    "appearance\n"
    "Appearance {\n"
    "material DEF Material Material {}\n"
    "}\n"
    "}\n"
    "Group {\n"
    "children [\n"
    "Shape {\n"
    "appearance Appearance { material USE Material }\n"
    "geometry IndexedFaceSet {\n"
    "convex FALSE\n"
    "solid  FALSE\n"
    "coord  Coordinate {\n"
    "point [\n";

  for(vertex_descriptor v : vertices(g)){
    os <<  get(vpm,v) << ",\n";
      reindex[v]=n++;
  }
  os << "] #point\n"
    "} #coord Coordinate\n"
    "coordIndex  [\n";
   for(face_descriptor f : faces(g)){
    for(vertex_descriptor v : vertices_around_face(halfedge(f,g),g)){
      os << reindex[v] << ",";
    }
    os << "-1,\n";
   }

  os << "] #coordIndex\n"
    "} #geometry\n"
    "} #Shape\n"
    "] #children\n"
    "} #group\n"
    "]\n"
    "}\n";

  return os.good();
}

template <typename FaceGraph>
bool write_wrl(std::ostream& os,
               const FaceGraph& g)
{
  return write_wrl(os, g,
                   parameters::all_default());
}

/*!
 \ingroup PkgBGLIOFct

  writes the graph `g` in the OFF format.

  \param os the output stream
  \param g the graph to be written
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `FaceGraph`.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
*/
template <typename FaceGraph, typename NamedParameters>
bool write_off(std::ostream& os,
               const FaceGraph& g,
               const NamedParameters& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;
  typedef typename boost::graph_traits<FaceGraph>::faces_size_type faces_size_type;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));
  vertices_size_type nv = static_cast<vertices_size_type>(std::distance(vertices(g).first, vertices(g).second));
  faces_size_type nf = static_cast<faces_size_type>(std::distance(faces(g).first, faces(g).second));

  os << "OFF\n" << nv << " " << nf << " 0\n";
  boost::container::flat_map<vertex_descriptor,vertices_size_type> reindex;
  int n = 0;
  for(vertex_descriptor v : vertices(g)){
    os << get(vpm,v) << '\n';
    reindex[v]=n++;
  }

  for(face_descriptor f : faces(g)){
    os << degree(f,g);
    for(vertex_descriptor v : vertices_around_face(halfedge(f,g),g)){
      os << " " << reindex[v];
    }
    os << '\n';
  }
  return os.good();
}


/*!
   \ingroup PkgBGLIOFct
    writes the graph `g` in the OFF format into a file named `fname`.
    \sa Overloads of this function for specific models of the concept `FaceGraph`.

  */
template <typename FaceGraph, typename NamedParameters>
bool write_off(const char* fname,
               const FaceGraph& g,
               const NamedParameters& np)
{
  std::ofstream out(fname);
  if(out.good()){
    return write_off(out,g, np);
  }
  return false;
}

template <typename FaceGraph, typename NamedParameters>
bool write_off(const std::string& fname,
               const FaceGraph& g,
               const NamedParameters& np)
{ return write_off(fname.c_str(), g, np); }


template <typename FaceGraph>
bool write_off(std::ostream& os,
               const FaceGraph& g)
{
  return write_off(os, g,
                   parameters::all_default());
}
template <typename FaceGraph>
bool write_off(const char* fname,
               const FaceGraph& g)
{
  return write_off(fname,g,
                   parameters::all_default());
}

template <typename FaceGraph>
bool write_off(const std::string& fname,
               const FaceGraph& g)
{ return write_off(fname, g,
                   parameters::all_default()); }

  namespace internal { namespace read_off_tools {

  inline bool is_whitespace(const std::string& s)
  {
    for(unsigned int i=0; i < s.size(); i++){
      if(s[i] != ' ' && s[i] != '\t'){
        return false;
      }
    }
    return true;
  }

inline std::string next_non_comment(std::istream& is)
{
  std::string line;
  do {
    std::getline(is, line);
  }while(line[0] == '#' || is_whitespace(line));
  return line;
}

    }
  } // namespace internal


/*!
 \ingroup PkgBGLIOFct

 reads the graph `g` from data in the OFF format. Ignores comment lines which start with a hash, and lines with whitespace.

 \param is the input stream
 \param g the graph to be read
 \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

 \cgalNamedParamsBegin
   \cgalParamNBegin{vertex_point_map}
     \cgalParamDescription{a property map associating points to the vertices of `g`}
     \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
                    as key type and `%Point_3` as value type}
     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                     must be available in `FaceGraph`.}
   \cgalParamNEnd
 \cgalNamedParamsEnd

 \sa Overloads of this function for specific models of the concept `FaceGraph`.
 \pre The data must represent a 2-manifold

 \attention The graph `g` is not cleared, and the data from the stream are added.
*/
template <typename FaceGraph, typename NamedParameters>
bool read_off(std::istream& is,
              FaceGraph& g,
              NamedParameters np)
{
  using namespace internal::read_off_tools;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;
  typedef typename boost::graph_traits<FaceGraph>::faces_size_type faces_size_type;

  typedef typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::type Vpm;
  typedef  typename boost::property_traits<Vpm>::value_type Point_3;

  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(CGAL::vertex_point, g));
  vertices_size_type nv, nvf;
  faces_size_type nf;
  int ignore;

  std::string line = next_non_comment(is);
  {
    std::istringstream iss(line);
    std::string off;
    iss >> off;
    CGAL_assertion( off == "OFF" || off == "COFF");
  }
  line = next_non_comment(is);
  {
    std::istringstream iss(line);
    iss >> nv >> nf >> ignore;
  }

  std::vector<vertex_descriptor> vertices(nv);
  Point_3 p;
  for(vertices_size_type i=0; i < nv; i++){
    line = next_non_comment(is);
    std::istringstream iss(line);
    iss >> p;
    vertices[i] = add_vertex(g);
    put(vpm,vertices[i],p);
  }

  for(faces_size_type i=0; i < nf; i++){
    line = next_non_comment(is);
    std::istringstream iss(line);
    iss >> nvf;
    std::vector<vertex_descriptor> face(nvf);
    for(vertices_size_type j = 0; j < nvf; j++){
      faces_size_type fvi;
      iss >> fvi;
      face[j] = vertices[fvi];
    }
    Euler::add_face(face,g);
  }
  return true;
}

template <typename FaceGraph>
bool read_off(std::istream& is,
              FaceGraph& g)
{
  return read_off(is, g, parameters::all_default());
}

/*!
   \ingroup PkgBGLIOFct
    reads the graph `g` from data in the OFF format. Ignores comment lines which start with a hash, and lines with whitespace.
    \sa Overloads of this function for specific models of the concept `FaceGraph`.
    \pre The data must represent a 2-manifold
    \attention The graph `g` is not cleared, and the data from the stream are added.

  */
template <typename FaceGraph, typename NamedParameters>
bool read_off(const char* fname,
              FaceGraph& g,
              NamedParameters np)
{
  std::ifstream in(fname);
  if(in.good()){
    return read_off(in, g, np);
  }
  return false;
}

template <typename FaceGraph>
bool read_off(const char* fname,
              FaceGraph& g)
{
  return read_off(fname, g, parameters::all_default());
}

template <typename FaceGraph, typename NamedParameters>
bool read_off(const std::string& fname,
              FaceGraph& g,
              NamedParameters np)
{ return read_off(fname.c_str(), g, np); }

template <typename FaceGraph>
bool read_off(const std::string& fname,
              FaceGraph& g)
{ return read_off(fname, g, parameters::all_default()); }

template <typename FaceGraph, typename NamedParameters>
bool write_inp(std::ostream& os,
               const FaceGraph& g,
               std::string name,
               std::string type,
               const NamedParameters& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;

  typedef typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::const_type VPM;
  typedef typename boost::property_traits<VPM>::value_type Point_3;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  os << "*Part, name=" << name << "\n*Node\n";
  boost::container::flat_map<vertex_descriptor,vertices_size_type> reindex;
  int n = 1;
  for(vertex_descriptor v : vertices(g)){
    Point_3 p =  get(vpm,v);
    os << n << ", " << p.x() << ", " << p.y() << ", " << p.z() << '\n';
    reindex[v]=n++;
  }
  n = 1;
  os << "*Element, type=" << type << std::endl;
  for(face_descriptor f : faces(g)){
    os << n++;
    for(vertex_descriptor v : vertices_around_face(halfedge(f,g),g)){
      os << ", " << reindex[v];
    }
    os << '\n';
  }
  os << "*End Part"<< std::endl;
  return os.good();
}
// conveniance overload
template <typename FaceGraph>
bool write_inp(std::ostream& os,
               const FaceGraph& g,
               std::string name,
               std::string type)
{
  return write_inp(os, g, name, type, parameters::all_default());
}

namespace internal {
  namespace write_vtp {

// writes the polys appended data at the end of the .vtp file
template <class Mesh,
          typename NamedParameters>
void
write_polys(std::ostream& os,
            const Mesh & mesh,
            const NamedParameters& np)
{
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_iterator face_iterator;

  typedef typename CGAL::GetInitializedVertexIndexMap<Mesh, NamedParameters>::const_type Vimap;
  Vimap V = CGAL::get_initialized_vertex_index_map(mesh, np);

  std::vector<std::size_t> connectivity_table;
  std::vector<std::size_t> offsets;
  std::vector<unsigned char> cell_type(num_faces(mesh),5);  // triangle == 5

  std::size_t off = 0;
  for( face_iterator fit = faces(mesh).begin() ;
       fit != faces(mesh).end() ;
       ++fit )
  {
    off += 3;
    offsets.push_back(off);
    for(vertex_descriptor v :
                  vertices_around_face(halfedge(*fit, mesh), mesh))
        connectivity_table.push_back(V[v]);
  }
  write_vector<std::size_t>(os,connectivity_table);
  write_vector<std::size_t>(os,offsets);
  write_vector<unsigned char>(os,cell_type);
}
//todo use named params for maps
template <class Mesh,
          typename NamedParameters>
void
write_polys_tag(std::ostream& os,
                const Mesh & mesh,
                bool binary,
                std::size_t& offset,
                const NamedParameters& np)
{
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_iterator face_iterator;

  typedef typename CGAL::GetInitializedVertexIndexMap<Mesh, NamedParameters>::const_type Vimap;
  Vimap V = CGAL::get_initialized_vertex_index_map(mesh, np);

  std::string formatattribute =
    binary ? " format=\"appended\"" : " format=\"ascii\"";

  std::string typeattribute;
  switch(sizeof(std::size_t)) {
  case 8: typeattribute = " type=\"UInt64\""; break;
  case 4: typeattribute = " type=\"UInt32\""; break;
  default: CGAL_error_msg("Unknown size of std::size_t");
  }

  // Write connectivity table
  os << "    <Polys>\n"
     << "      <DataArray Name=\"connectivity\""
     << formatattribute << typeattribute;

  if (binary) { // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (3 * num_faces(mesh)+ 1) * sizeof(std::size_t);
    // 3 indices (size_t) per triangle + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";
    for( face_iterator fit = faces(mesh).begin() ;
         fit != faces(mesh).end() ;
         ++fit )
    {
      for(vertex_descriptor v :
                    vertices_around_face(halfedge(*fit, mesh), mesh))
          os << V[v] << " ";
    }
    os << "      </DataArray>\n";
  }

  // Write offsets
  os   << "      <DataArray Name=\"offsets\""
       << formatattribute << typeattribute;

  if (binary) {  // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (num_faces(mesh) + 1) * sizeof(std::size_t);
    // 1 offset (size_t) per triangle + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";
    std::size_t polys_offset = 0;
    for( face_iterator fit = faces(mesh).begin() ;
         fit != faces(mesh).end() ;
         ++fit )
    {
      polys_offset += 3;
      os << polys_offset << " ";
    }
    os << "      </DataArray>\n";
  }

  // Write cell type (triangle == 5)
  os   << "      <DataArray Name=\"types\""
       << formatattribute << " type=\"UInt8\"";

  if (binary) {
    os << " offset=\"" << offset << "\"/>\n";
    offset += num_faces(mesh) + sizeof(std::size_t);
    // 1 unsigned char per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";
    for(std::size_t i = 0; i< num_faces(mesh); ++i)
      os << "5 ";
    os << "      </DataArray>\n";
  }
  os << "    </Polys>\n";
}

//todo : use namedparams for points and ids
//overload for facegraph
template <class Mesh,
          typename NamedParameters>
void
write_points_tag(std::ostream& os,
                 const Mesh & mesh,
                 bool binary,
                 std::size_t& offset,
                 const NamedParameters& np)
{
  typedef typename boost::graph_traits<Mesh>::vertex_iterator vertex_iterator;
  typedef typename CGAL::GetVertexPointMap<Mesh, NamedParameters>::const_type Vpmap;
  using parameters::get_parameter;
  using parameters::choose_parameter;
  Vpmap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, mesh));
  typedef typename boost::property_traits<Vpmap>::value_type Point_t;
  typedef typename CGAL::Kernel_traits<Point_t>::Kernel Gt;
  typedef typename Gt::FT FT;

  std::string format = binary ? "appended" : "ascii";
  std::string type = (sizeof(FT) == 8) ? "Float64" : "Float32";

  os << "    <Points>\n"
     << "      <DataArray type =\"" << type << "\" NumberOfComponents=\"3\" format=\""
     << format;

  if (binary) {
    os << "\" offset=\"" << offset << "\"/>\n";
    offset += 3 * num_vertices(mesh) * sizeof(FT) + sizeof(std::size_t);
    // 3 coords per points + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";
    for( vertex_iterator vit = vertices(mesh).begin();
         vit != vertices(mesh).end();
         ++vit)
      {
      os << get(vpm, *vit).x() << " " << get(vpm, *vit).y() << " "
         << get(vpm, *vit).z() << " ";
      }
    os << "      </DataArray>\n";
  }
  os << "    </Points>\n";
}


// writes the points appended data at the end of the .vtp file
template <class Mesh,
          class NamedParameters>
void
write_polys_points(std::ostream& os,
                   const Mesh & mesh,
                   const NamedParameters& np)
{
  typedef typename boost::graph_traits<Mesh>::vertex_iterator vertex_iterator;
  typedef typename CGAL::GetVertexPointMap<Mesh, NamedParameters>::const_type Vpmap;
  using parameters::get_parameter;
  using parameters::choose_parameter;
  Vpmap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, mesh));
  typedef typename boost::property_traits<Vpmap>::value_type Point_t;
  typedef typename CGAL::Kernel_traits<Point_t>::Kernel Gt;
  typedef typename Gt::FT FT;
  std::vector<FT> coordinates;
  for( vertex_iterator vit = vertices(mesh).begin();
       vit != vertices(mesh).end();
       ++vit)
    {
      coordinates.push_back(get(vpm, *vit).x());
      coordinates.push_back(get(vpm, *vit).y());
      coordinates.push_back(get(vpm, *vit).z());
    }
  write_vector<FT>(os,coordinates);
}

  } // end namespace CGAL::internal::write_vtp
} // end namespace CGAL::internal

/*!\ingroup PkgBGLIOFct
 *
 * \brief  writes a triangulated surface mesh in the `PolyData` XML format.
 *
 * \tparam TriangleMesh a model of `FaceListGraph` with only triangle faces.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param os the stream used for writing.
 * \param mesh the triangle mesh to be written.
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{use_binary_mode}
 *     \cgalParamDescription{Boolean indicating if the data should be written in binary (`true`) or in ASCII (`false`)}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `mesh`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_index_map}
 *     \cgalParamDescription{a property map associating to each vertex of `mesh` a unique index between `0` and `num_vertices(mesh) - 1`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{an automatically indexed internal map}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 */
template<class TriangleMesh,
         class NamedParameters>
void write_vtp(std::ostream& os,
               const TriangleMesh& mesh,
               const NamedParameters& np)
{
  os << "<?xml version=\"1.0\"?>\n"
     << "<VTKFile type=\"PolyData\" version=\"0.1\"";
#ifdef CGAL_LITTLE_ENDIAN
  os << " byte_order=\"LittleEndian\"";
#else // CGAL_BIG_ENDIAN
  os << " byte_order=\"BigEndian\"";
#endif
  switch(sizeof(std::size_t)) {
  case 4: os << " header_type=\"UInt32\""; break;
  case 8: os << " header_type=\"UInt64\""; break;
  default: CGAL_error_msg("Unknown size of std::size_t");
  }
  os << ">\n"
     << "  <PolyData>" << "\n";

  os << "  <Piece NumberOfPoints=\"" << num_vertices(mesh)
     << "\" NumberOfPolys=\"" << num_faces(mesh) << "\">\n";
  std::size_t offset = 0;
  const bool binary = parameters::choose_parameter(parameters::get_parameter(np, internal_np::use_binary_mode), true);
  internal::write_vtp::write_points_tag(os,mesh,binary,offset, np);
  internal::write_vtp::write_polys_tag(os,mesh,binary,offset, np);
  os << "   </Piece>\n"
     << "  </PolyData>\n";
  if (binary) {
    os << "<AppendedData encoding=\"raw\">\n_";
    internal::write_vtp::write_polys_points(os,mesh, np);
    internal::write_vtp::write_polys(os,mesh, np);
  }
  os << "</VTKFile>\n";
}

template<class TriangleMesh>
void write_vtp(std::ostream& os,
               const TriangleMesh& mesh)
{
  write_vtp(os, mesh, CGAL::parameters::all_default());
}

} // namespace CGAL

#endif // CGAL_BOOST_GRAPH_IO_H
