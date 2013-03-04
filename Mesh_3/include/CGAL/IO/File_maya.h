// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
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
// $URL:  $
// $Id:  $
//
//
// Author(s)     : Clement Jamin

#ifndef CGAL_IO_FILE_MAYA_H
#define CGAL_IO_FILE_MAYA_H

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <CGAL/utility.h>

namespace CGAL {

//-------------------------------------------------------
// IO functions
//-------------------------------------------------------
  
template <class C3T3>
void
output_to_maya(std::ostream& os, 
               const C3T3& c3t3, 
               bool surfaceOnly = true)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Facets_in_complex_iterator Facet_iterator;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point_3;

#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cerr << "Output to maya:\n";
#endif
  
  const Tr& tr = c3t3.triangulation();

  //-------------------------------------------------------
  // File output
  //-------------------------------------------------------

  //-------------------------------------------------------
  // Header
  //-------------------------------------------------------
  os << std::setprecision(20);

  os << "//Maya ASCII 2011 scene" << std::endl;
  os << "//Name: testMaya3.ma" << std::endl;
  os << "//Last modified: Wed, Jan 25, 2012 05:54:26 PM" << std::endl;
  os << "//Codeset: 1252" << std::endl;
  os << "requires maya \"2011\";" << std::endl;
  os << "currentUnit -l centimeter -a degree -t film;" << std::endl;
  os << "fileInfo \"application\" \"maya\";" << std::endl;
  os << "fileInfo \"product\" \"Maya 2011\";" << std::endl;
  os << "fileInfo \"version\" \"2011\";" << std::endl;
  os << "fileInfo \"cutIdentifier\" \"201003190014-771504\";" << std::endl;
  os << "fileInfo \"license\" \"education\";" << std::endl;

  std::string name = "Mesh_3";
  os << "createNode mesh -n \"" << name << "Shape\" -p \"" << name << "\";" << std::endl;
  os << "  setAttr -k off \".v\";" << std::endl;
  os << "  setAttr \".uvst[0].uvsn\" -type \"string\" \"map1\";" << std::endl;
  os << "  setAttr \".cuvs\" -type \"string\" \"map1\";" << std::endl;
  os << "  setAttr \".dcol\" yes;" << std::endl;
  os << "  setAttr \".dcc\" -type \"string\" \"Ambient+Diffuse\";" << std::endl;
  
  os << "  connectAttr \"" << name << "Shape.iog\" \":initialShadingGroup.dsm\" -na;\n\n";
    
  //-------------------------------------------------------
  // Colors
  //------------------------------------------------------
  
  /*os << "  setAttr \".ccls\" -type \"string\" \"colorSet\";\n";
  os << "  setAttr \".clst[0].clsn\" -type \"string\" \"colorSet\";\n";
  os << "  setAttr \".clst[0].rprt\" 3;\n";
  os << "  setAttr -s " << 3 << " \".clst[0].clsp[0:" << 3-1 << "]\"" << std::endl;
  os << "    10 50 250" << std::endl;
  os << "    100 250 50" << std::endl;
  os << "    0 200 200" << std::endl;
  os << "  ;\n";*/
  
  //-------------------------------------------------------
  // Vertices
  //------------------------------------------------------
  
  std::map<Vertex_handle, int> V;
  std::stringstream vertices_sstr;
  int num_vertices = 0;
  for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
       vit != tr.finite_vertices_end();
       ++vit)
  {
    if ( (surfaceOnly  && c3t3.in_dimension(vit) <= 2)
      || !surfaceOnly)
    {
      V[vit] = num_vertices++;
      Point_3 p = vit->point();
      vertices_sstr << "    " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " " << CGAL::to_double(p.z()) << std::endl;
    }
  }
    
  os << "  setAttr -s " << num_vertices << " \".vt[0:" << num_vertices-1 << "]\"" << std::endl;
  os << vertices_sstr.str();
  os << ";\n";

  /*
  // Triangles
  os << "setAttr -s " << QString().setNum(number_of_triangles) << " \".fc[0:" << QString().setNum(number_of_triangles-1) << "]\"  -type \"polyFaces\" \n";
  for (int i=0;i<triangle.size();i=i+3){
    int a=int(triangle[i+0]);
    int b=int(triangle[i+1]);
    int c=int(triangle[i+2]);
    os << "f 3 " << QString().setNum(a) << " " << QString().setNum(b) << " " << QString().setNum(c) << " ";
    if (mat.materialNode==VertexColor){
      os << "mc 0 3 " << QString().setNum(index[i*2]) << " " << QString().setNum(index[i*2+2]) << " " << QString().setNum(index[i*2+4]) << " ";
    }
  }
  os << ";\n\n";*/

  //-------------------------------------------------------
  // Get edges and facets
  //-------------------------------------------------------

  typename C3T3::size_type number_of_triangles = c3t3.number_of_facets_in_complex();
  
  std::stringstream facets_sstr;
  //std::stringstream normals_sstr;

  //normals_sstr << "  setAttr -s " << number_of_triangles*3 << " \".n[0:" << number_of_triangles*3-1 << "]\"  -type \"float3\" \n";

  // Save edges
  typedef std::vector<std::pair<int, int> > EdgeList;
  EdgeList edges;

  // Surface only
  if (surfaceOnly)
  {
    facets_sstr <<  "  setAttr -s " << number_of_triangles 
      << " \".fc[0:" << number_of_triangles-1 << "]\"  -type \"polyFaces\" \n";
    int c = 0;
    for( Facet_iterator fit = c3t3.facets_in_complex_begin();
         fit != c3t3.facets_in_complex_end();
         ++fit, ++c)
    {
      int indices[3];
      //Point_3 points[3];
      facets_sstr << "    f 3 ";
      for (int j = 0, i = (fit->second + 1) % 4 ; j < 3 ; i = (i+1)%4, ++j)
      {
        const Vertex_handle& vh = fit->first->vertex(i);
        indices[j] = V[vh];
        //points[j] = vh->point();
      }
    
      // Reverse triangle orientation?
      bool reverse_triangle = 
           (fit->second % 2 == 0 && !c3t3.is_in_complex(fit->first))
        || (fit->second % 2 != 0 && c3t3.is_in_complex(fit->first));
      if (reverse_triangle)
      {
        std::swap(indices[1], indices[2]);
        //std::swap(points[1], points[2]);
      }
      //Kernel::Vector_3 n = cross_product(points[1] - points[0], points[2] - points[0]);
      //n = n / CGAL::sqrt(n*n);
      // Add the normal 3 times
      //normals_sstr << "    " << n.x() << " " << n.y() << " " << n.z() << std::endl;
      //normals_sstr << "    " << n.x() << " " << n.y() << " " << n.z() << std::endl;
      //normals_sstr << "    " << n.x() << " " << n.y() << " " << n.z() << std::endl;

      // 3 edges
      for (int i = 0 ; i < 3 ; ++i)
      {
        std::pair<int, int> edge = std::make_pair(
          (std::min)(indices[i], indices[(i+1)%3]),
          (std::max)(indices[i], indices[(i+1)%3]));
        size_t pos = std::find(edges.begin(), edges.end(), edge) - edges.begin();
        if (pos == edges.size()) // Not found?
        {
          edges.push_back(edge);
        }
        // ith edge of triangle
        facets_sstr << pos << " ";
      }
    
      // 1 triangles
      facets_sstr << std::endl;
      // Colors
      //facets_sstr << "    mc 0 3 " << rand()%3 << " " << rand()%3 << " " << rand()%3 << std::endl;
    }
  }
  // Tetrahedra = 4 facets for each
  else
  {
    facets_sstr <<  "  setAttr -s " << 4*c3t3.number_of_cells_in_complex() 
      << " \".fc[0:" << 4*c3t3.number_of_cells_in_complex()-1 << "]\"  -type \"polyFaces\" \n";
    int c = 0;
    for( Cell_iterator cit = c3t3.cells_in_complex_begin();
         cit != c3t3.cells_in_complex_end();
         ++cit, ++c)
    {
      for (int facet_i = 0 ; facet_i < 4 ; ++facet_i)
      {
        int indices[3];
        //Point_3 points[3];
        facets_sstr << "    f 3 ";
        for (int j = 0, i = (facet_i + 1) % 4 ; j < 3 ; i = (i+1)%4, ++j)
        {
          const Vertex_handle& vh = cit->vertex(i);
          indices[j] = V[vh];
          //points[j] = vh->point();
        }
    
        // Reverse triangle orientation?
        bool reverse_triangle = (facet_i % 2 != 0 && c3t3.is_in_complex(cit, facet_i));
        if (reverse_triangle)
        {
          std::swap(indices[1], indices[2]);
          //std::swap(points[1], points[2]);
        }
        //Kernel::Vector_3 n = cross_product(points[1] - points[0], points[2] - points[0]);
        //n = n / CGAL::sqrt(n*n);
        // Add the normal 3 times
        //normals_sstr << "    " << n.x() << " " << n.y() << " " << n.z() << std::endl;
        //normals_sstr << "    " << n.x() << " " << n.y() << " " << n.z() << std::endl;
        //normals_sstr << "    " << n.x() << " " << n.y() << " " << n.z() << std::endl;

        // 3 edges
        for (int i = 0 ; i < 3 ; ++i)
        {
          std::pair<int, int> edge = std::make_pair(
            (std::min)(indices[i], indices[(i+1)%3]),
            (std::max)(indices[i], indices[(i+1)%3]));
          size_t pos = std::find(edges.begin(), edges.end(), edge) - edges.begin();
          if (pos == edges.size()) // Not found?
          {
            edges.push_back(edge);
          }
          // ith edge of triangle
          facets_sstr << pos << " ";
        }
    
        // 1 triangles
        facets_sstr << std::endl;
        // Colors
        //facets_sstr << "    mc 0 3 " << rand()%3 << " " << rand()%3 << " " << rand()%3 << std::endl;
      }
    }
  }
  facets_sstr << ";\n\n";
  //normals_sstr << ";\n\n";
  
  //-------------------------------------------------------
  // Edges
  //-------------------------------------------------------
  os << "  setAttr -s " << edges.size() << " \".ed[0:" 
     << edges.size() - 1 << "]\"" << std::endl;

  for (EdgeList::const_iterator it = edges.begin(), it_end = edges.end() ; it != it_end ; ++it)
    os << "    " << it->first << " " << it->second << " " << 0 << std::endl;

  os << ";\n";
  
  //-------------------------------------------------------
  // Normals
  //-------------------------------------------------------
  
  //os << normals_sstr.str();

  //-------------------------------------------------------
  // Facets
  //-------------------------------------------------------
  
  os << facets_sstr.str();
  
  //-------------------------------------------------------
  // Tetrahedra
  //-------------------------------------------------------
  /*os << "Tetrahedra" << std::endl
     << c3t3.number_of_cells_in_complex() << std::endl;

  for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
       cit != c3t3.cells_in_complex_end() ;
       ++cit )
  {
    for (int i=0; i<4; i++)
      os << V[cit->vertex(i)] << " ";

    os << get(cell_pmap, cit) << std::endl;
  }*/

  //-------------------------------------------------------
  // End
  //-------------------------------------------------------

#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cerr << "done.\n";
#endif
} // end output_to_maya(...)

} // end namespace CGAL

#endif // CGAL_IO_FILE_MAYA_H
