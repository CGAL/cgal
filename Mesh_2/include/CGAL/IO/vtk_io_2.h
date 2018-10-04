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
//
//
// Author(s)     : Laurent RINEAU, Stephane Tayeb, Maxime Gimeno

#ifndef CGAL_VTK_IO_H
#define CGAL_VTK_IO_H

#include <fstream>
#include <vector>

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

// writes the cells tags before binary data is appended

template <class CDT>
void 
write_cells_tag_2(std::ostream& os,
                  const CDT & tr,
                  std::size_t number_of_triangles,
                  std::map<typename CDT::Vertex_handle, std::size_t> & V,
                  bool binary,
                  std::size_t& offset)
{
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
    offset += (3 * number_of_triangles + 1) * sizeof(std::size_t); 
    // 3 indices (size_t) per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";   
    for(typename CDT::Finite_faces_iterator 
            fit = tr.finite_faces_begin(),
            end = tr.finite_faces_end();
          fit != end; ++fit)
      {
      if(fit->is_in_domain())
      {
        os << V[fit->vertex(0)] << " ";
        os << V[fit->vertex(2)] << " ";
        os << V[fit->vertex(1)] << " ";
      }
    }
    os << "      </DataArray>\n";
  }
  
  // Write offsets
  os   << "      <DataArray Name=\"offsets\""
       << formatattribute << typeattribute;
  
  if (binary) {  // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (number_of_triangles + 1) * sizeof(std::size_t);
    // 1 offset (size_t) per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    std::size_t cells_offset = 0;
    for(typename CDT::Finite_faces_iterator fit = 
        tr.finite_faces_begin() ;
        fit != tr.finite_faces_end() ;
        ++fit )
    {
      if(fit->is_in_domain())
      {
        cells_offset += 3;
        os << cells_offset << " ";
      }
    }  
    os << "      </DataArray>\n";
  }

  // Write cell type (triangles == 5)
  os   << "      <DataArray Name=\"types\""
       << formatattribute << " type=\"UInt8\"";

  if (binary) {
    os << " offset=\"" << offset << "\"/>\n";
    offset += number_of_triangles + sizeof(std::size_t);
    // 1 unsigned char per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    for(typename CDT::Finite_faces_iterator fit = 
        tr.finite_faces_begin() ;
        fit != tr.finite_faces_end() ;
        ++fit )
    {
      if(fit->is_in_domain())
      {
        os << "5 ";
      }
    }
    os << "      </DataArray>\n";
  }
  os << "    </Cells>\n";
}

// writes the cells appended data at the end of the .vtu file 
template <class CDT>
void
write_cells_2(std::ostream& os,
              const CDT & tr,
              std::size_t number_of_triangles,
              std::map<typename CDT::Vertex_handle, std::size_t> & V)
{
  std::vector<std::size_t> connectivity_table;
  std::vector<std::size_t> offsets;
  std::vector<unsigned char> cell_type(number_of_triangles,5);  // triangles == 5
  
  std::size_t off = 0;
  for(typename CDT::Finite_faces_iterator 
      fit = tr.finite_faces_begin(),
      end = tr.finite_faces_end();
      fit != end; ++fit)
  {
    if(fit->is_in_domain())
    {
      off += 3;
      offsets.push_back(off);
      connectivity_table.push_back(V[fit->vertex(0)]);
      connectivity_table.push_back(V[fit->vertex(2)]);
      connectivity_table.push_back(V[fit->vertex(1)]);
    }
  }
  write_vector<std::size_t>(os,connectivity_table);
  write_vector<std::size_t>(os,offsets);
  write_vector<unsigned char>(os,cell_type);
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
  std::size_t dim = 2;
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
write_points(std::ostream& os,
	     const Tr & tr,
	     std::map<typename Tr::Vertex_handle, 
             std::size_t> & V)
{
  std::size_t dim = 2;
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


template <class CDT>
void write_unstructured_grid_2(std::ostream& os,
                               const CDT& tr,
                               bool binary = true)
{
  typedef typename CDT::Vertex_handle Vertex_handle;
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
  
  int number_of_triangles = 0;
  for(typename CDT::Finite_faces_iterator 
      fit = tr.finite_faces_begin(),
      end = tr.finite_faces_end();
      fit != end; ++fit)
  {
    if(fit->is_in_domain()) ++number_of_triangles;
  }
  os << "  <Piece NumberOfPoints=\"" << tr.number_of_vertices() 
     << "\" NumberOfCells=\"" << number_of_triangles << "\">\n";
  std::size_t offset = 0;
  write_points_tag(os,tr,V,binary,offset);
  write_cells_tag_2(os,tr,number_of_triangles, V,binary,offset);
  os << "   </Piece>\n"
     << "  </UnstructuredGrid>\n";
  if (binary) {
    os << "<AppendedData encoding=\"raw\">\n_"; 
    write_points(os,tr,V);  // write points before cells to fill the std::map V
    write_cells_2(os,tr, number_of_triangles, V);
  }
  os << "</VTKFile>\n";
}

} //end CGAL
#endif // CGAL_VTK_IO_H
