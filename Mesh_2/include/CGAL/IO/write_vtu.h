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
//
// Author(s)     : Laurent RINEAU, Stephane Tayeb, Maxime Gimeno

#ifndef CGAL_WRITE_VTU_H
#define CGAL_WRITE_VTU_H

#include <CGAL/license/Mesh_2.h>

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <CGAL/assertions.h>
#include <CGAL/IO/io.h>

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
    // 3 indices (size_t) per triangle + length of the encoded data (size_t)
    offset += (3 * number_of_triangles + 1) * sizeof(std::size_t);
    // 2 indices (size_t) per edge (size_t)
    offset += (2 * std::distance(tr.constrained_edges_begin(),
                                 tr.constrained_edges_end())) * sizeof(std::size_t);
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
    offset += (number_of_triangles +std::distance(tr.constrained_edges_begin(),
                                                  tr.constrained_edges_end()) + 1)
        * sizeof(std::size_t);
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
    offset += number_of_triangles
        + std::distance(tr.constrained_edges_begin(),
                        tr.constrained_edges_end())
        + sizeof(std::size_t);
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
  cell_type.resize(cell_type.size() + std::distance(tr.constrained_edges_begin(),
                                                     tr.constrained_edges_end()), 3);  // line == 3

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
  for(typename CDT::Constrained_edges_iterator
      cei = tr.constrained_edges_begin(),
      end = tr.constrained_edges_end();
      cei != end; ++cei)
  {
    off += 2;
    offsets.push_back(off);
    for(int i=0; i<3; ++i)
    {
      if(i != cei->second)
        connectivity_table.push_back(V[cei->first->vertex(i)]);
    }
  }
  write_vector<std::size_t>(os,connectivity_table);
  write_vector<std::size_t>(os,offsets);
  write_vector<unsigned char>(os,cell_type);
}

// writes the points tags before binary data is appended
template <class Tr>
void
write_cdt_points_tag(std::ostream& os,
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
write_cdt_points(std::ostream& os,
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

// writes the attribute tags before binary data is appended
template <class T>
void
write_attribute_tag_2 (std::ostream& os,
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
write_attributes_2(std::ostream& os,
                 const std::vector<FT>& att)
{
  write_vector(os,att);
}

template <class CDT>
void write_vtu_with_attributes(std::ostream& os,
               const CDT& tr,
               std::vector<std::pair<const char*, const std::vector<double>*> >& attributes,
               IO::Mode mode = IO::BINARY)
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
     << "\" NumberOfCells=\"" << number_of_triangles + std::distance(tr.constrained_edges_begin(), tr.constrained_edges_end()) << "\">\n";
  std::size_t offset = 0;
  const bool binary = (mode == IO::BINARY);
  write_cdt_points_tag(os,tr,V,binary,offset);
  write_cells_tag_2(os,tr,number_of_triangles, V,binary,offset);
  if(attributes.empty())
      os << "    <CellData >\n";
  else
    os << "    <CellData Scalars=\""<<attributes.front().first<<"\">\n";
  for(std::size_t i = 0; i< attributes.size(); ++i)
  {
    write_attribute_tag_2(os,attributes[i].first, *attributes[i].second, binary,offset);
  }
  os << "    </CellData>\n";
  os << "   </Piece>\n"
     << "  </UnstructuredGrid>\n";
  if (binary) {
    os << "<AppendedData encoding=\"raw\">\n_";
    write_cdt_points(os,tr,V);  // write points before cells to fill the std::map V
    write_cells_2(os,tr, number_of_triangles, V);
    for(std::size_t i = 0; i< attributes.size(); ++i)
      write_attributes_2(os, *attributes[i].second);
  }
  os << "</VTKFile>\n";
}


template <class CDT>
void write_vtu(std::ostream& os,
               const CDT& tr,
               IO::Mode mode = IO::BINARY)
{
  std::vector<std::pair<const char*, const std::vector<double>*> > dummy_atts;
  write_vtu_with_attributes(os, tr, dummy_atts, mode);
}

} //end CGAL
#endif // CGAL_WRITE_VTU_H
