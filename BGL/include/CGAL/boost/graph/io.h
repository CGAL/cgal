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
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_IO_H
#define CGAL_BOOST_GRAPH_IO_H

#include <boost/container/flat_map.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL {
/*!
   \ingroup PkgBGLIOFct
    writes the graph `g` in the OFF format.
    \sa Overloads of this function for specific models of the concept `FaceGraph`.

  */ 
template <typename FaceGraph>
bool write_off(std::ostream& os,
               const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;
  typedef typename boost::graph_traits<FaceGraph>::faces_size_type faces_size_type;

  typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::const_type vpm = get(CGAL::vertex_point,g);
  vertices_size_type nv = static_cast<vertices_size_type>(std::distance(vertices(g).first, vertices(g).second));
  faces_size_type nf = static_cast<faces_size_type>(std::distance(faces(g).first, faces(g).second));

  os << "OFF\n" << nv << " " << nf << " 0\n";
  boost::container::flat_map<vertex_descriptor,vertices_size_type> reindex;
  int n = 0;
  BOOST_FOREACH(vertex_descriptor v, vertices(g)){
    os << get(vpm,v) << '\n';
    reindex[v]=n++;
  }
  
  BOOST_FOREACH(face_descriptor f, faces(g)){
    os << degree(f,g);
    BOOST_FOREACH(vertex_descriptor v, vertices_around_face(halfedge(f,g),g)){
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
template <typename FaceGraph>
bool write_off(const char* fname,
               const FaceGraph& g)
{
  std::ofstream out(fname);
  if(out.good()){
    return write_off(out,g);
  }
  return false;
}


  namespace internal { namespace read_off_tools {
  
  bool is_whitespace(const std::string& s)
  {
    for(unsigned int i=0; i < s.size(); i++){
      if(s[i] != ' ' && s[i] != '\t'){
        return false;
      }
    }
    return true;
  }
  
std::string next_non_comment(std::istream& is)
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
    \sa Overloads of this function for specific models of the concept `FaceGraph`.
    \pre The data must represent a 2-manifold
    \attention The graph `g` is not cleared, and the data from the stream are added.

  */ 
template <typename FaceGraph>
bool read_off(std::istream& is,
              FaceGraph& g)
{
  using namespace internal::read_off_tools;

  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;
  typedef typename boost::graph_traits<FaceGraph>::faces_size_type faces_size_type;

  typedef typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type Vpm;
  typedef  typename boost::property_traits<Vpm>::value_type Point_3;
  Vpm vpm = get(CGAL::vertex_point,g);
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


/*!
   \ingroup PkgBGLIOFct
    reads the graph `g` from data in the OFF format. Ignores comment lines which start with a hash, and lines with whitespace.
    \sa Overloads of this function for specific models of the concept `FaceGraph`.
    \pre The data must represent a 2-manifold
    \attention The graph `g` is not cleared, and the data from the stream are added.

  */ 
template <typename FaceGraph>
bool read_off(const char* fname,
              FaceGraph& g)
{
  std::ifstream in(fname);
  if(in.good()){
    return read_off(in, g);
  }
  return false;
}


} // namespace CGAL

#endif // CGAL_BOOST_GRAPH_IO_H
