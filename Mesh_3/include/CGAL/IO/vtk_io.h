// Copyright (c) 2018 GeometryFactory (France).
// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Laurent RINEAU, Stephane Tayeb, Maxime Gimeno

#ifndef CGAL_VTK_IO_H
#define CGAL_VTK_IO_H

#include <CGAL/license/Mesh_3.h>


#include <fstream>
#include <vector>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>


//todo try to factorize with functors
namespace CGAL{
// writes the appended data into the .vtu file
template <class FT> 
void
write_vector(std::ostream& os, 
             const std::vector<FT>& vect) 
{
  const char* buffer = reinterpret_cast<const char*>(&(vect[0]));
  std::size_t size = vect.size()*sizeof(FT);
  
  os.write(reinterpret_cast<const char *>(&size), sizeof(std::size_t)); // number of bytes encoded
  os.write(buffer, vect.size()*sizeof(FT));                     // encoded data
}

template <class C3T3>
void 
write_cells_tag(std::ostream& os,
                const C3T3 & c3t3,
                std::map<typename C3T3::Triangulation::Vertex_handle, std::size_t> & V,
                bool binary,
                std::size_t& offset)
{
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
  std::string formatattribute =
    binary ? " format=\"appended\"" : " format=\"ascii\"";

  std::string typeattribute;
  switch(sizeof(std::size_t)) {
  case 8: typeattribute = " type=\"UInt64\""; break;
  case 4: typeattribute = " type=\"UInt32\""; break;
  default: CGAL_error_msg("Unknown size of std::size_t");
  }

  // Write connectivity table
  os << "    <Cells>\n"
     << "      <DataArray Name=\"connectivity\""
     << formatattribute << typeattribute;
  
  if (binary) { // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (4 * c3t3.number_of_cells() + 1) * sizeof(std::size_t);
    // 4 indices (size_t) per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";   
    for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
         cit != c3t3.cells_in_complex_end() ;
         ++cit )
    {
      for (int i=0; i<4; i++)
        os << V[cit->vertex(i)] << " ";
    }
    os << "      </DataArray>\n";
  }
  
  // Write offsets
  os   << "      <DataArray Name=\"offsets\""
       << formatattribute << typeattribute;
  
  if (binary) {  // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (c3t3.number_of_cells() + 1) * sizeof(std::size_t);
    // 1 offset (size_t) per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    std::size_t cells_offset = 0;
    for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
	 cit != c3t3.cells_in_complex_end() ;
	 ++cit )
      {
	cells_offset += 4;
	os << cells_offset << " ";
      }  
    os << "      </DataArray>\n";
  }

  // Write cell type (tetrahedra == 10)
  os   << "      <DataArray Name=\"types\""
       << formatattribute << " type=\"UInt8\"";

  if (binary) {
    os << " offset=\"" << offset << "\"/>\n";
    offset += c3t3.number_of_cells() + sizeof(std::size_t);
    // 1 unsigned char per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
	 cit != c3t3.cells_in_complex_end() ;
	 ++cit )
      os << "10 ";
    os << "      </DataArray>\n";
  }
  os << "    </Cells>\n";
}


template <class C3T3>
void
write_cells(std::ostream& os,
            const C3T3 & c3t3,
            std::map<typename C3T3::Triangulation::Vertex_handle, std::size_t> & V,
            std::vector<float>& mids)
{
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
  std::vector<std::size_t> connectivity_table;
  std::vector<std::size_t> offsets;
  std::vector<unsigned char> cell_type(c3t3.number_of_cells(),10);  // tetrahedra == 10
  
  std::size_t off = 0;
  for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
       cit != c3t3.cells_in_complex_end() ;
       ++cit )
    {
      off += 4;
      offsets.push_back(off);
      for (int i=0; i<4; i++)
	connectivity_table.push_back(V[cit->vertex(i)]);
      mids.push_back(cit->subdomain_index());
    }
  write_vector<std::size_t>(os,connectivity_table);
  write_vector<std::size_t>(os,offsets);
  write_vector<unsigned char>(os,cell_type);
}

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
  typedef typename CGAL::Polygon_mesh_processing::GetVertexPointMap<Mesh, NamedParameters>::const_type Vpmap;
  typedef typename CGAL::Polygon_mesh_processing::GetVertexIndexMap<Mesh, NamedParameters>::type Vimap;
  Vimap V = choose_param(get_param(np, CGAL::internal_np::vertex_index),
                           get_const_property_map(CGAL::internal_np::vertex_index, mesh));
  
  typedef typename boost::property_traits<Vpmap>::value_type Point_t;
  typedef typename CGAL::Kernel_traits<Point_t>::Kernel Gt;
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
  typedef typename CGAL::Polygon_mesh_processing::GetVertexPointMap<Mesh, NamedParameters>::const_type Vpmap;
  typedef typename CGAL::Polygon_mesh_processing::GetVertexIndexMap<Mesh, NamedParameters>::type Vimap;
  Vimap V = choose_param(get_param(np, CGAL::internal_np::vertex_index),
                           get_const_property_map(CGAL::internal_np::vertex_index, mesh));
  
  
  typedef typename boost::property_traits<Vpmap>::value_type Point_t;
  typedef typename CGAL::Kernel_traits<Point_t>::Kernel Gt;
  
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
// writes the points tags before binary data is appended
template <class Tr>
void 
write_points_tag(std::ostream& os,
                 const Tr & tr,
                 std::map<typename Tr::Vertex_handle, std::size_t> & V,
                 bool binary,
                 std::size_t& offset) 
{
  std::size_t dim = 3;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::FT FT;

  std::size_t inum = 0;
  std::string format = binary ? "appended" : "ascii";
  std::string type = (sizeof(FT) == 8) ? "Float64" : "Float32";

  os << "    <Points>\n"
     << "      <DataArray type =\"" << type << "\" NumberOfComponents=\"3\" format=\""
     << format;

  if (binary) {
    os << "\" offset=\"" << offset << "\"/>\n";    
    offset += 3 * tr.number_of_vertices() * sizeof(FT) + sizeof(std::size_t);
    // dim coords per points + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end();
         ++vit)
    {
      V[vit] = inum++;
        os << vit->point()[0] << " ";
        os << vit->point()[1] << " ";
        if(dim == 3) 
          os << vit->point()[2] << " ";
        else
          os << 0.0 << " ";
      }
    os << "      </DataArray>\n";
  }
  os << "    </Points>\n";
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

// writes the points appended data at the end of the .vtu file 
template <class Tr>
void
write_points(std::ostream& os,
	     const Tr & tr,
	     std::map<typename Tr::Vertex_handle, 
             std::size_t> & V)
{
  std::size_t dim = 3;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::FT FT;

  std::size_t inum = 0;
  std::vector<FT> coordinates;
  for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
       vit != tr.finite_vertices_end();
       ++vit)
    {
      V[vit] = inum++;  // binary output => the map has not been filled yet
      coordinates.push_back(vit->point()[0]);
      coordinates.push_back(vit->point()[1]);
      coordinates.push_back(dim == 3 ? vit->point()[2] : 0.0);
    }
  write_vector<FT>(os,coordinates);
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
// writes the attribute tags before binary data is appended
template <class T>
void 
write_attribute_tag(std::ostream& os,
		    const std::string& attr_name,
		    const std::vector<T>& attribute,
		    bool binary,
		    std::size_t& offset)
{
  std::string format = binary ? "appended" : "ascii";
  std::string type = (sizeof(T) == 8) ? "Float64" : "Float32";
  os << "      <DataArray type=\"" << type << "\" Name=\"" << attr_name << "\" format=\"" << format; 

  if (binary) {
    os << "\" offset=\"" << offset << "\"/>\n";    
    offset += attribute.size() * sizeof(T) + sizeof(std::size_t);
  }
  else {
    typedef typename std::vector<T>::const_iterator Iterator;
    os << "\">\n";   
    for (Iterator it = attribute.begin();
	 it != attribute.end();
	 ++it )
      os << *it << " ";
    os << "      </DataArray>\n";
  }
}

// writes the attributes appended data at the end of the .vtu file 
template <typename FT>
void
write_attributes(std::ostream& os,
		 const std::vector<FT>& att)
{
  write_vector(os,att);
}

//public API

template<class TriangleMesh, 
         class NamedParameters>
void write_polydata(std::ostream& os,
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
void write_polydata(std::ostream& os,
                      const TriangleMesh& mesh,
                      bool binary = true)
{
  write_polydata(os, mesh, binary, CGAL::parameters::all_default());
} 


template <class C3T3>
void write_unstructured_grid_3(std::ostream& os,
                             const C3T3& c3t3)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename Tr::Vertex_handle Vertex_handle;
  const Tr& tr = c3t3.triangulation();
  std::map<Vertex_handle, std::size_t> V;
  //write header
  os << "<?xml version=\"1.0\"?>\n"
     << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
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
     << "  <UnstructuredGrid>" << "\n";
  
  os << "  <Piece NumberOfPoints=\"" << tr.number_of_vertices() 
     << "\" NumberOfCells=\"" << c3t3.number_of_cells() << "\">\n";
  bool binary = true;
  std::size_t offset = 0;
  write_points_tag(os,tr,V,binary,offset);
  write_cells_tag(os,c3t3,V,binary,offset);
  std::vector<float> mids;      
    os << "   <CellData Domain=\"MeshDomain";
    os << "\">\n";
    write_attribute_tag(os,"MeshDomain",mids,binary,offset);
    os << "    </CellData>\n";
  os << "   </Piece>\n"
     << "  </UnstructuredGrid>\n";
  if (binary) {
    os << "<AppendedData encoding=\"raw\">\n_"; 
    write_points(os,tr,V);  // write points before cells to fill the std::map V
    write_cells(os,c3t3,V, mids);//todo mids should be filled by write_attribute_tag
    write_attributes(os,mids);
  }
  os << "</VTKFile>\n";
}

} //end CGAL
#endif // CGAL_VTK_IO_H
