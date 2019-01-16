// Copyright (c) 2018 GeometryFactory (France).
// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// Author(s)     : Laurent RINEAU, Stephane Tayeb, Maxime Gimeno

#ifndef CGAL_VTP_IO_H
#define CGAL_VTP_IO_H

#include <CGAL/license/Polyhedron.h>


#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/IO/write_vtk.h>


//todo try to factorize with functors
namespace CGAL{
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
    BOOST_FOREACH(vertex_descriptor v,
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
      BOOST_FOREACH(vertex_descriptor v,
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

//public API

template<class TriangleMesh, 
         class NamedParameters>
void write_VTP(std::ostream& os,
                    const TriangleMesh& mesh,
                    bool binary,
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
  write_points_tag(os,mesh,binary,offset, np);
  write_polys_tag(os,mesh,binary,offset, np);
  os << "   </Piece>\n"
     << "  </PolyData>\n";
  if (binary) {
    os << "<AppendedData encoding=\"raw\">\n_"; 
    write_polys_points(os,mesh, np);
    write_polys(os,mesh, np);
  }
  os << "</VTKFile>\n";
}

template<class TriangleMesh>
void write_VTP(std::ostream& os,
                      const TriangleMesh& mesh,
                      bool binary = true)
{
  write_VTP(os, mesh, binary, CGAL::parameters::all_default());
}

} //end CGAL
#endif // CGAL_VTP_IO_H
