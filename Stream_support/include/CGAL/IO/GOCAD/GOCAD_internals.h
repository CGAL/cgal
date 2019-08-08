// Copyright (c) 2019 GeometryFactory
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
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_GOCAD_GOCAD_INTERNALS_H
#define CGAL_IO_GOCAD_GOCAD_INTERNALS_H

#include <deque>
#include <iostream>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <CGAL/boost/graph/Euler_operations.h>
#include <iostream>

namespace CGAL{

namespace GOCAD_internal {
template <class Facegraph, class P>
class Build_from_gocad_BGL
{
  typedef P Point_3;
  typedef std::deque<Point_3> Points_3;
  typedef boost::tuple<int,int,int> Facet;
  typedef std::deque<Facet> Surface;

  std::istream& is;
  int counter;
  Points_3 meshPoints;
  Surface mesh;

public:

  std::string name, color;


  Build_from_gocad_BGL(std::istream& is_)
    : is(is_), counter(0)
  {}

  void operator()( Facegraph& graph) {
    read(is, meshPoints, mesh);

    typedef typename boost::graph_traits<Facegraph>::vertex_descriptor
        vertex_descriptor;

    std::vector<vertex_descriptor> vertices(meshPoints.size());
    for(std::size_t id = 0; id < meshPoints.size(); ++id)
    {
      vertices[id] = add_vertex( meshPoints[id], graph);
    }
//    graph.begin_surface( meshPoints.size(), mesh.size());
    typedef typename Points_3::size_type size_type;

    for(size_type i=0; i < mesh.size(); i++){
      std::array<vertex_descriptor, 3> face;
      face[0] = vertices[mesh[i].template get<0>()];
      face[1] = vertices[mesh[i].template get<1>()];
      face[2] = vertices[mesh[i].template get<2>()];

      CGAL::Euler::add_face(face, graph);
    }
  }

  void
  read(std::istream& input, Points_3& points, Surface& surface, int offset = 0)
  {
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
        str >> name;
      }
      if((idx = s.find("color")) != std::string::npos){
        std::istringstream str(s.substr(idx+6));
        str >> color;
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
        surface.push_back(boost::make_tuple(offset+i, offset+j, offset+k));
      } else if(c == 'E'){
        break;
      }
      std::getline(input, s);
    }
  }

};

template <typename FaceGraph>
bool
read_gocad(FaceGraph& polyhedron, std::istream& in, std::string& name, std::string& color)
{
  //typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename boost::property_traits<typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type>::value_type Point_3;

  Build_from_gocad_BGL<FaceGraph, Point_3> builder(in);
  builder(polyhedron);
  name=builder.name;
  color=builder.color;

  return in.good() && polyhedron.is_valid();
}

template <typename FaceGraph>
bool
write_gocad(FaceGraph& polyhedron, std::ostream& os, const std::string& name)
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
  VPMap vpmap = get(CGAL::vertex_point, polyhedron);
  std::map<typename boost::graph_traits<FaceGraph>::vertex_descriptor, int> id_map;
  {
    typename boost::graph_traits<FaceGraph>::vertex_iterator it, end;
    it = vertices(polyhedron).begin();
    end = vertices(polyhedron).end();
    int i=0;
    for(; it != end; ++it){
      id_map[*it] = i;
      os << "VRTX " << i << " " << get(vpmap, *it) << "\n";
      ++i;
    }
  }

  {
    typename boost::graph_traits<FaceGraph>::face_iterator it, end;
    it = faces(polyhedron).begin();
    end = faces(polyhedron).end();
    for(; it != end; ++it){
      os << "TRGL " << id_map[target(prev(halfedge(*it, polyhedron), polyhedron), polyhedron)] << " "
         << id_map[target(halfedge(*it, polyhedron), polyhedron)] << " "
         << id_map[target(next(halfedge(*it, polyhedron), polyhedron), polyhedron)] << "\n";
    }
  }

  os << "END" << std::endl;

  return true;
}

}//end GOCAD_internal
}//end CGAL
#endif // CGAL_IO_GOCAD_GOCAD_INTERNALS_H
 
