// Copyright (c) 2018 GeometryFactory (France).
// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
#include <boost/variant.hpp>

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
    os << ">\n";
    for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
         cit != c3t3.cells_in_complex_end() ;
         ++cit )
    {
      for (int i=0; i<4; i++)
        os << V[cit->vertex(i)] << " ";
    }
    os << "\n      </DataArray>\n";
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
    os << ">\n";
    std::size_t cells_offset = 0;
    for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
         cit != c3t3.cells_in_complex_end() ;
         ++cit )
      {
        cells_offset += 4;
        os << cells_offset << " ";
      }
    os << "\n      </DataArray>\n";
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
    os << ">\n";
    for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
         cit != c3t3.cells_in_complex_end() ;
         ++cit )
      os << "10 ";
    os << "\n      </DataArray>\n";
  }
  os << "    </Cells>\n";
}


template <class C3T3>
void
write_cells(std::ostream& os,
            const C3T3 & c3t3,
            std::map<typename C3T3::Triangulation::Vertex_handle, std::size_t> & V)
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
    }
  write_vector<std::size_t>(os,connectivity_table);
  write_vector<std::size_t>(os,offsets);
  write_vector<unsigned char>(os,cell_type);
}


template <class Tr>
void
write_c3t3_points_tag(std::ostream& os,
                      const Tr & tr,
                      std::size_t size_of_vertices,
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
    offset += 3 * size_of_vertices * sizeof(FT) + sizeof(std::size_t);
    // dim coords per points + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";
    for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end();
         ++vit)
    {
      if(vit->in_dimension() <= -1) continue;
      V[vit] = inum++;
      os << vit->point()[0] << " ";
      os << vit->point()[1] << " ";
      if(dim == 3)
        os << vit->point()[2] << " ";
      else
        os << 0.0 << " ";
    }
    os << "\n      </DataArray>\n";
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
      if(vit->in_dimension() <= -1) continue;
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
  std::string type = "";
  if(std::is_floating_point<T>::value)
  {
    type = (sizeof(T) == 8) ? "Float64" : "Float32";
  }
  else
  {
    if(sizeof(T) == 1)
      type =  "UInt8";
    else if(sizeof(T) == 4)
      type = "UInt32";
    else
      type = "UInt64";
  }
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
    os << "\n      </DataArray>\n";
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

enum VTU_ATTRIBUTE_TYPE{
  DOUBLE=0,
  UNIT_8,
  SIZE_TYPE
};

typedef boost::variant<const std::vector<double>*, const std::vector<uint8_t>*, const std::vector<std::size_t>* > Vtu_attributes;

template <class C3T3>
void output_to_vtu_with_attributes(std::ostream& os,
                                   const C3T3& c3t3,
                                   std::vector<std::pair<const char*, Vtu_attributes> >&attributes,
                                   IO::Mode mode = IO::BINARY)
{
  //CGAL_assertion(attributes.size() == attribute_types.size());
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
  const std::size_t number_of_vertices =
    tr.number_of_vertices() -  c3t3.number_of_far_points();
  os << "  <Piece NumberOfPoints=\"" << number_of_vertices
     << "\" NumberOfCells=\"" << c3t3.number_of_cells() << "\">\n";
  std::size_t offset = 0;

  const bool binary = (mode == IO::BINARY);
  write_c3t3_points_tag(os,tr,number_of_vertices,V,binary,offset);
  write_cells_tag(os,c3t3,V,binary,offset); // fills V if the mode is ASCII
  os << "    <CellData Scalars=\""<<attributes.front().first<<"\">\n";
  for(std::size_t i = 0; i< attributes.size(); ++i)
  {
    switch(attributes[i].second.which()){
    case 0:
      write_attribute_tag(os,attributes[i].first, *boost::get<const std::vector<double>* >(attributes[i].second), binary,offset);
      break;
    case 1:
      write_attribute_tag(os,attributes[i].first, *boost::get<const std::vector<uint8_t>* >(attributes[i].second), binary,offset);
      break;
    default:
      write_attribute_tag(os,attributes[i].first, *boost::get<const std::vector<std::size_t>* >(attributes[i].second), binary,offset);
      break;
    }
  }
  os << "    </CellData>\n";
  os << "   </Piece>\n"
     << "  </UnstructuredGrid>\n";
  if (binary) {
    os << "<AppendedData encoding=\"raw\">\n_";
    write_c3t3_points(os,tr,V); // fills V if the mode is BINARY
    write_cells(os,c3t3,V);
    for(std::size_t i = 0; i< attributes.size(); ++i)
      switch(attributes[i].second.which()){
      case 0:
        write_attributes(os, *boost::get<const std::vector<double>* >(attributes[i].second));
        break;
      case 1:
        write_attributes(os, *boost::get<const std::vector<uint8_t>* >(attributes[i].second));
        break;
      default:
        write_attributes(os, *boost::get<const std::vector<std::size_t>* >(attributes[i].second));
        break;
      }
  }
  os << "</VTKFile>\n";
}



//public API
template <class C3T3>
void output_to_vtu(std::ostream& os,
               const C3T3& c3t3,
               IO::Mode mode = IO::BINARY)
{
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
  std::vector<double> mids;
  for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
       cit != c3t3.cells_in_complex_end() ;
       ++cit )
  {
    double v = cit->subdomain_index();
    mids.push_back(v);
  }

  std::vector<std::pair<const char*, Vtu_attributes > > atts;
  Vtu_attributes v = &mids;
  atts.push_back(std::make_pair("MeshDomain", v));
  output_to_vtu_with_attributes(os, c3t3, atts, mode);
}

} //end CGAL
#endif // CGAL_VTK_IO_H
