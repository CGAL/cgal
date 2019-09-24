// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
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
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Kernel_traits.h>
#ifdef CGAL_USE_VTK
#include <CGAL/IO/VTK/vtk_internals.h>
#endif
#include <CGAL/IO/write_vtk.h>
#include <CGAL/internal/Generic_facegraph_builder.h>
#include <CGAL/IO/STL/STL_reader.h>
#include <CGAL/IO/OBJ/OBJ_reader.h>
#include <CGAL/IO/OBJ/File_writer_wavefront.h>

namespace CGAL {
  /*!
   \ingroup PkgBGLIOFct
    writes the graph `g` in the wrl format (VRML 2.0).
    
    \cgalNamedParamsBegin
    *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
    *       If this parameter is omitted, an internal property map for
    *       `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
    * \cgalNamedParamsEnd
    */ 
template <typename FaceGraph, typename NamedParameters>
bool write_wrl(std::ostream& os,
               const FaceGraph& g,
               const NamedParameters& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;
  
  typename Polygon_mesh_processing::GetVertexPointMap<FaceGraph, NamedParameters>::const_type
      vpm = choose_param(get_param(np, internal_np::vertex_point),
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
    
    \cgalNamedParamsBegin
    *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
    *       If this parameter is omitted, an internal property map for
    *       `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
    * \cgalNamedParamsEnd
    
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

  typename Polygon_mesh_processing::GetVertexPointMap<FaceGraph, NamedParameters>::const_type
      vpm = choose_param(get_param(np, internal_np::vertex_point),
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
    
    \cgalNamedParamsBegin
    *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
    *       If this parameter is omitted, an internal property map for
    *       `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
    * \cgalNamedParamsEnd
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

  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;
  typedef typename boost::graph_traits<FaceGraph>::faces_size_type faces_size_type;

  typedef typename Polygon_mesh_processing::GetVertexPointMap<FaceGraph, NamedParameters>::type Vpm;
  typedef  typename boost::property_traits<Vpm>::value_type Point_3;
  
  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
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

  typedef typename Polygon_mesh_processing::GetVertexPointMap<FaceGraph, NamedParameters>::const_type VPM;
  typedef typename boost::property_traits<VPM>::value_type Point_3;

  VPM vpm = choose_param(get_param(np, internal_np::vertex_point),
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

namespace GOCAD_internal{
//Use CRTP to gain access to the protected members without getters/setters.
template <class Facegraph, class P>
class GOCAD_builder : public CGAL::internal::IO::Generic_facegraph_builder<Facegraph, P, GOCAD_builder<Facegraph, P> >
{
  typedef GOCAD_builder<Facegraph, P> Self;
  typedef CGAL::internal::IO::Generic_facegraph_builder<Facegraph, P, Self> Base;
  typedef typename Base::Point_3 Point_3;
  typedef typename Base::Points_3 Points_3;
  typedef typename Base::Facet Facet;
  typedef typename Base::Surface Surface;
public:
  GOCAD_builder(std::istream& is_)
    :Base(is_){}
  void do_construct(Facegraph& graph)
  {
    typedef typename boost::graph_traits<Facegraph>::vertex_descriptor
        vertex_descriptor;

    std::vector<vertex_descriptor> vertices(this->meshPoints.size());
    for(std::size_t id = 0; id < this->meshPoints.size(); ++id)
    {
      vertices[id] = add_vertex( this->meshPoints[id], graph);
    }
    //    graph.begin_surface( meshPoints.size(), mesh.size());
    typedef typename Points_3::size_type size_type;

    for(size_type i=0; i < this->mesh.size(); i++){
      std::array<vertex_descriptor, 3> face;
      face[0] = vertices[this->mesh[i][0]];
      face[1] = vertices[this->mesh[i][1]];
      face[2] = vertices[this->mesh[i][2]];

      CGAL::Euler::add_face(face, graph);
    }
  }

  void
  read(std::istream& input, Points_3& points, Surface& surface)
  {
    int offset = 0;
    char c;
    std::string s, tface("TFACE");
    int i,j,k;
    Point_3 p;
    bool vertices_read = false;
    while(input >> s){
      if(s == tface){
        break;
      }
      std::string::size_type idx;

      if((idx = s.find("name")) != std::string::npos){
        std::istringstream str(s.substr(idx+5));
        str >> this->name;
      }
      if((idx = s.find("color")) != std::string::npos){
        std::istringstream str(s.substr(idx+6));
        str >> this->color;
      }
    }
    std::getline(input, s);

    while(input.get(c)){
      if((c == 'V')||(c == 'P')){
        input >> s >> i >> p;
        if(! vertices_read){
          vertices_read = true;
          offset -= i; // Some files start with index 0 others with 1
        }

        points.push_back(p);

      } else if(vertices_read && (c == 'T')){
        input >> c >> c >> c >>  i >> j >> k;
        typename Base::Facet new_face(3);
        new_face[0] = offset+i;
        new_face[1] = offset+j;
        new_face[2] = offset+k;
        surface.push_back(new_face);
      } else if(c == 'E'){
        break;
      }
      std::getline(input, s);
    }
  }

};
}//end GOCAD_internal

/*!
   \ingroup PkgBGLIOFct
    reads the graph `face_graph` from data in the TS format.
    `name` and `color` will be filled according to the values contained in the file.

    \pre The data must represent a 2-manifold
    \attention The graph `face_graph` is not cleared, and the data from the stream are added.
    \see \ref IOStreamGocad
  */
template <typename FaceGraph>
bool
read_gocad(FaceGraph& face_graph, std::istream& in, std::string& name, std::string& color)
{
  //typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename boost::property_traits<typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type>::value_type Point_3;

  GOCAD_internal::GOCAD_builder<FaceGraph, Point_3> builder(in);
  builder(face_graph);
  name=builder.name;
  color=builder.color;

  return in.good() && face_graph.is_valid();
}

/*!
   \ingroup PkgBGLIOFct
    writes the graph `face_graph` in the TS format into `os`. `name` is the
    mandatory name that will be assigned to `face_graph`in the file.
  */
template <typename FaceGraph>
bool
write_gocad(FaceGraph& face_graph, std::ostream& os, const std::string& name)
{
  os << "GOCAD TSurf 1\n"
    "HEADER {\n"
    "name:";
  os << name << std::endl;
  os << "*border:on\n"
    "*border*bstone:on\n"
    "}\n"
    "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n"
    "NAME Default\n"
    "AXIS_NAME \"X\" \"Y\" \"Z\"\n"
    "AXIS_UNIT \"m\" \"m\" \"m\"\n"
    "ZPOSITIVE Elevation\n"
    "END_ORIGINAL_COORDINATE_SYSTEM\n"
    "TFACE\n";

  os.precision(16);
  typedef typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type VPMap;
  VPMap vpmap = get(CGAL::vertex_point, face_graph);
  std::map<typename boost::graph_traits<FaceGraph>::vertex_descriptor, int> id_map;
  {
    typename boost::graph_traits<FaceGraph>::vertex_iterator it, end;
    it = vertices(face_graph).begin();
    end = vertices(face_graph).end();
    int i=0;
    for(; it != end; ++it){
      id_map[*it] = i;
      os << "VRTX " << i << " " << get(vpmap, *it) << "\n";
      ++i;
    }
  }

  {
    typename boost::graph_traits<FaceGraph>::face_iterator it, end;
    it = faces(face_graph).begin();
    end = faces(face_graph).end();
    for(; it != end; ++it){
      os << "TRGL " << id_map[target(prev(halfedge(*it, face_graph), face_graph), face_graph)] << " "
         << id_map[target(halfedge(*it, face_graph), face_graph)] << " "
         << id_map[target(next(halfedge(*it, face_graph), face_graph), face_graph)] << "\n";
    }
  }

  os << "END" << std::endl;

  return true;
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
  typedef typename CGAL::Polygon_mesh_processing::GetVertexIndexMap<Mesh, NamedParameters>::type Vimap;
  Vimap V = choose_param(get_param(np, CGAL::internal_np::vertex_index),
                           get_const_property_map(CGAL::internal_np::vertex_index, mesh));

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
        connectivity_table.push_back(get(V, v));
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
  typedef typename CGAL::Polygon_mesh_processing::GetVertexIndexMap<Mesh, NamedParameters>::type Vimap;
  Vimap V = choose_param(get_param(np, CGAL::internal_np::vertex_index),
                           get_const_property_map(CGAL::internal_np::vertex_index, mesh));

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
          os << get(V, v) << " ";
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
  typedef typename CGAL::Polygon_mesh_processing::GetVertexPointMap<Mesh, NamedParameters>::const_type Vpmap;
  Vpmap vpm = choose_param(get_param(np, CGAL::vertex_point),
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
  typedef typename CGAL::Polygon_mesh_processing::GetVertexPointMap<Mesh, NamedParameters>::const_type Vpmap;
  Vpmap vpm = choose_param(get_param(np, CGAL::vertex_point),
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

/*! \ingroup PkgBGLIOFct
 *
 * \brief  writes a triangulated surface mesh in the `PolyData` XML format.
 *
 * \tparam TriangleMesh a model of `FaceListGraph` with only triangle faces.
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param os the stream used for writing.
 * \param mesh the triangle mesh to be written.
 * \param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the
 * ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{use_binary_mode} a Boolean indicating if the
 *    data should be written in binary (`true`, the default) or in ASCII (`false`).
 *     \cgalParamEnd
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to
 * the vertices of `mesh`. If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `TriangleMesh`.
 *     \cgalParamEnd
 *    \cgalParamBegin{vertex_index_map} the property map with the indices associated to
 * the vertices of `mesh`. If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_index_t` must be available in `TriangleMesh`.
 *     \cgalParamEnd
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
  const bool binary = boost::choose_param(boost::get_param(np, internal_np::use_binary_mode), true);
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

#ifdef CGAL_USE_VTK
namespace VTK_internal{

template <typename FaceGraph>
bool vtkPointSet_to_polygon_mesh(vtkPointSet* poly_data,
                                 FaceGraph& face_graph)
{
  typedef typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type VPMap;
  typedef typename boost::property_map_value<FaceGraph, CGAL::vertex_point_t>::type Point_3;
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;

  VPMap vpmap = get(CGAL::vertex_point, face_graph);

  // get nb of points and cells
  vtkIdType nb_points = poly_data->GetNumberOfPoints();
  vtkIdType nb_cells = poly_data->GetNumberOfCells();
  //extract points
  std::vector<vertex_descriptor> vertex_map(nb_points);
  for (vtkIdType i = 0; i<nb_points; ++i)
  {
    double coords[3];
    poly_data->GetPoint(i, coords);

    vertex_descriptor v = add_vertex(face_graph);
    put(vpmap, v, Point_3(coords[0], coords[1], coords[2]));
    vertex_map[i]=v;
  }

  //extract cells
  for (vtkIdType i = 0; i<nb_cells; ++i)
  {
    int cell_type = poly_data->GetCellType(i);
    if(cell_type != 5
       && cell_type != 7
       && cell_type != 9) //only supported cells are triangles, quads and polygons
      continue;
    vtkCell* cell_ptr = poly_data->GetCell(i);

    vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
    if (nb_vertices < 3)
      return false;
    std::vector<vertex_descriptor> vr(nb_vertices);
    for (vtkIdType k=0; k<nb_vertices; ++k){
      vtkIdType id = cell_ptr->GetPointId(k);
      vr[k]=vertex_map[id];
    }

    CGAL::Euler::add_face(vr, face_graph);
  }
  return true;
}
} //end VTK_internal

template<class FaceGraph>
bool read_vtp(const char* filename, FaceGraph& face_graph)
{
  vtkSmartPointer<vtkPointSet> data;
  if(!CGAL::read_vtp_file(filename, data))
  {
    return false;
  }
  return VTK_internal::vtkPointSet_to_polygon_mesh(data, face_graph);
}
#endif //CGAL_USE_VTK

#ifdef DOXYGEN_RUNNING
/*! \ingroup PkgBGLIOFct
 * \brief  reads a PolyData in the VTP format into a triangulated surface mesh.
 *
 * \tparam FaceGraph a model of `FaceListGraph`.
 *
 * \param filename the path to the file that will be read.
 * \param face_graph the output mesh.
 *
 * \pre \cgal needs to be configured with the VTK Libraries for this function to be available.
 */
template<class FaceGraph>
bool read_vtp(const char* filename, FaceGraph& face_graph);

#endif

/*!
  \ingroup PkgBGLIOFct
  writes the graph `tm` in the stream `out` in the STL format.
  \pre The graph must contain only triangle faces.
  */
template <class TriangleMesh>
std::ostream&
write_STL(const TriangleMesh& tm, std::ostream& out)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::const_type Vpm;
  typedef typename boost::property_traits<Vpm>::reference Point_3_ref;
  typedef typename boost::property_traits<Vpm>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel::Vector_3 Vector_3;

  Vpm vpm = get(boost::vertex_point, tm);

  if (get_mode(out) == IO::BINARY)
  {
    out << "FileType: Binary                                                                ";
    const boost::uint32_t N32 = static_cast<boost::uint32_t>(faces(tm).size());
    out.write(reinterpret_cast<const char *>(&N32), sizeof(N32));

    for(face_descriptor f : faces(tm))
    {
      halfedge_descriptor h = halfedge(f, tm);
      Point_3_ref p = get(vpm, target(h, tm));
      Point_3_ref q = get(vpm, target(next(h, tm), tm));
      Point_3_ref r = get(vpm, source(h, tm));

      Vector_3 n = collinear(p,q,r) ? Vector_3(1,0,0):
                                      unit_normal(p,q,r);

      const float coords[12]={
        static_cast<float>(n.x()), static_cast<float>(n.y()), static_cast<float>(n.z()),
        static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()),
        static_cast<float>(q.x()), static_cast<float>(q.y()), static_cast<float>(q.z()),
        static_cast<float>(r.x()), static_cast<float>(r.y()), static_cast<float>(r.z()) };

      for (int i=0; i<12; ++i)
        out.write(reinterpret_cast<const char *>(&coords[i]), sizeof(coords[i]));
      out << "  ";
    }
  }
  else
  {
    out << "solid\n";
    for(face_descriptor f : faces(tm))
    {
      halfedge_descriptor h = halfedge(f, tm);
      Point_3_ref p = get(vpm, target(h, tm));
      Point_3_ref q = get(vpm, target(next(h, tm), tm));
      Point_3_ref r = get(vpm, source(h, tm));

      Vector_3 n = collinear(p,q,r) ? Vector_3(1,0,0):
                                      unit_normal(p,q,r);
      out << "facet normal " << n << "\nouter loop\n";
      out << "vertex " << p << "\n";
      out << "vertex " << q << "\n";
      out << "vertex " << r << "\n";
      out << "endloop\nendfacet\n";
    }
    out << "endsolid\n";
  }
  return out;
}


namespace STL_internal
{
//Use CRTP to gain access to the protected members without getters/setters.
template <class Facegraph, class P>
class STL_builder : public CGAL::internal::IO::Generic_facegraph_builder<Facegraph, P, STL_builder<Facegraph, P> >
{
  typedef STL_builder<Facegraph, P> Self;
  typedef CGAL::internal::IO::Generic_facegraph_builder<Facegraph, P, Self> Base;
  typedef typename Base::Point_3 Point_3;
  typedef typename Base::Points_3 Points_3;
  typedef typename Base::Facet Facet;
  typedef typename Base::Surface Surface;
public:
  STL_builder(std::istream& is_)
    :Base(is_){}
  void do_construct(Facegraph& graph)
  {
    typedef typename boost::graph_traits<Facegraph>::vertex_descriptor
        vertex_descriptor;

    std::vector<vertex_descriptor> vertices(this->meshPoints.size());
    for(std::size_t id = 0; id < this->meshPoints.size(); ++id)
    {
      vertices[id] = add_vertex( this->meshPoints[id], graph);
    }
    //    graph.begin_surface( meshPoints.size(), mesh.size());
    typedef typename Points_3::size_type size_type;

    for(size_type i=0; i < this->mesh.size(); i++){
      std::array<vertex_descriptor, 3> face;
      face[0] = vertices[this->mesh[i][0]];
      face[1] = vertices[this->mesh[i][1]];
      face[2] = vertices[this->mesh[i][2]];

      CGAL::Euler::add_face(face, graph);
    }
  }

  void
  read(std::istream& input, Points_3& points, Surface& surface)
  {
    read_STL(input, points, surface);
  }

};
} // end STL_internal

/*!
  \ingroup PkgBGLIOFct
  reads the graph `tm` from the stream `in` in the STL format.
  \pre The data must represent a 2-manifold
  */
template <class TriangleMesh>
bool
read_STL(TriangleMesh& tm, std::istream& in)
{
  //typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename boost::property_traits<typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type>::value_type Point_3;

  STL_internal::STL_builder<TriangleMesh, Point_3> builder(in);
  builder(tm);
  bool ok = in.good() || in.eof();
  ok &= tm.is_valid();
  return ok;
}







namespace OBJ_internal
{
//Use CRTP to gain access to the protected members without getters/setters.
template <class Facegraph, class P>
class OBJ_builder : public CGAL::internal::IO::Generic_facegraph_builder<Facegraph, P, OBJ_builder<Facegraph, P> >
{
  typedef OBJ_builder<Facegraph, P> Self;
  typedef CGAL::internal::IO::Generic_facegraph_builder<Facegraph, P, Self> Base;
  typedef typename Base::Point_3 Point_3;
  typedef typename Base::Points_3 Points_3;
  typedef typename Base::Facet Facet;
  typedef typename Base::Surface Surface;
public:
  OBJ_builder(std::istream& is_)
    :Base(is_){}
  void do_construct(Facegraph& graph)
  {
    typedef typename boost::graph_traits<Facegraph>::vertex_descriptor
        vertex_descriptor;

    std::vector<vertex_descriptor> vertices(this->meshPoints.size());
    for(std::size_t id = 0; id < this->meshPoints.size(); ++id)
    {
      vertices[id] = add_vertex( this->meshPoints[id], graph);
    }
    typedef typename Points_3::size_type size_type;

    for(size_type i=0; i < this->mesh.size(); i++){
      std::vector<vertex_descriptor> face(this->mesh[i].size());
      for(std::size_t j=0; j< face.size(); ++j)
        face[j] = vertices[this->mesh[i][j]];

      CGAL::Euler::add_face(face, graph);
    }
  }

  void
  read(std::istream& input, Points_3& points, Surface& surface)
  {
    read_OBJ(input, points, surface);
  }

};
} // end STL_internal

/*!
  \ingroup PkgBGLIOFct
  reads the graph `tm` from the stream `in` in the OBJ format.
  \returns `true` if the resulting mesh is valid.
  \pre The data must represent a 2-manifold
  */
template <class TriangleMesh>
bool
read_OBJ(TriangleMesh& tm, std::istream& in)
{
  //typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename boost::property_traits<typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type>::value_type Point_3;

  OBJ_internal::OBJ_builder<TriangleMesh, Point_3> builder(in);
  builder(tm);
  bool ok = in.good() || in.eof();
  ok &= tm.is_valid();
  return ok;
}

template <typename FaceGraph>
bool
write_OBJ(const FaceGraph& face_graph, std::ostream& os)
{
  // writes M to `out' in the format provided by `writer'.
  CGAL::File_writer_wavefront writer;
  typedef typename boost::graph_traits<FaceGraph >::vertex_iterator VCI;
  typedef typename boost::graph_traits<FaceGraph >::face_iterator   FCI;
  typedef typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type VPmap;
  VPmap map = get(CGAL::vertex_point, face_graph);
  // Print header.
  writer.write_header( os,
                       num_vertices(face_graph),
                       num_halfedges(face_graph),
                       num_faces(face_graph));

  std::map<typename boost::graph_traits<FaceGraph>::vertex_descriptor, std::size_t> index_map;
  auto hint = index_map.begin();
  std::size_t id = 0;

  for( VCI vi = vertices(face_graph).begin(); vi != vertices(face_graph).end(); ++vi) {
    writer.write_vertex( ::CGAL::to_double( get(map, *vi).x()),
                         ::CGAL::to_double( get(map, *vi).y()),
                         ::CGAL::to_double( get(map, *vi).z()));

    hint = index_map.insert(hint, std::make_pair(*vi, id++));
  }

  writer.write_facet_header();
  for( FCI fi = faces(face_graph).begin(); fi != faces(face_graph).end(); ++fi) {
    CGAL::Halfedge_around_face_circulator<FaceGraph> hc(halfedge(*fi, face_graph), face_graph);
    auto hc_end = hc;
    std::size_t n = circulator_size( hc);
    CGAL_assertion( n >= 3);
    writer.write_facet_begin( n);
    do {
      writer.write_facet_vertex_index(index_map[target(*hc, face_graph)]);
      ++hc;
    } while( hc != hc_end);
    writer.write_facet_end();
  }
  writer.write_footer();
  return os.good();
}
} // namespace CGAL


#endif // CGAL_BOOST_GRAPH_IO_H
