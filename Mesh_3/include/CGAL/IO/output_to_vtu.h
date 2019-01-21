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

#ifndef CGAL_OUTPUT_TO_VTU_H
#define CGAL_OUTPUT_TO_VTU_H

#include <CGAL/license/Mesh_3.h>


#include <iostream>
#include <vector>
#include <map>
#include <CGAL/assertions.h>
#include <CGAL/IO/io.h>
#include <CGAL/IO/write_vtk.h>

//todo try to factorize with functors
namespace CGAL{

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


template <class Tr>
void 
write_c3t3_points_tag(std::ostream& os,
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

// writes the points appended data at the end of the .vtu file 
template <class Tr>
void
write_c3t3_points(std::ostream& os,
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
template <class C3T3>
void output_to_vtu(std::ostream& os,
                   const C3T3& c3t3,
                   IO::Mode mode = IO::BINARY)
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
  std::size_t offset = 0;
  const bool binary = (mode == IO::BINARY);
  write_c3t3_points_tag(os,tr,V,binary,offset);
  write_cells_tag(os,c3t3,V,binary,offset); // fills V if the mode is ASCII
  std::vector<float> mids;      
  os << "   <CellData Domain=\"MeshDomain";
  os << "\">\n";
  write_attribute_tag(os,"MeshDomain",mids,binary,offset);
  os << "    </CellData>\n";
  os << "   </Piece>\n"
     << "  </UnstructuredGrid>\n";
  if (binary) {
    os << "<AppendedData encoding=\"raw\">\n_"; 
    write_c3t3_points(os,tr,V); // fills V if the mode is BINARY
    write_cells(os,c3t3,V, mids);
    write_attributes(os,mids);
  }
  os << "</VTKFile>\n";
}

} //end CGAL
#endif // CGAL_VTK_IO_H
