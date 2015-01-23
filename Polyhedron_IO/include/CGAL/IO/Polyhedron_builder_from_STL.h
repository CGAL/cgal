// Copyright (c) 2015 GeometryFactory
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_IO_POLYHEDRON_STL_BUILDER_H
#define CGAL_IO_POLYHEDRON_STL_BUILDER_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/tuple.h>

#include <iostream>

namespace CGAL{

template <class HDS>
class Polyhedron_builder_from_STL : public CGAL::Modifier_base<HDS> {
  typedef typename HDS::Vertex::Point Point_3;
  typedef std::vector<Point_3> Points_3;
  typedef cpp11::tuple<int,int,int> Facet;
  typedef std::vector<Facet> Surface;

  std::istream& is;
  int counter;
  Points_3 meshPoints;
  Surface mesh;

  bool
  read(std::istream& input, Points_3& points, Surface& surface, int /*offset*/ = 0)
  {
    std::string s, solid("solid"), facet("facet"), outer("outer"), loop("loop"), vertex("vertex"), endloop("endloop"), endsolid("endsolid");

    std::map<Point_3, int> vertex_index;
    int index = 0;
    int ijk[3];
    Point_3 p;

    input >> s;
    if(s == solid){
      std::getline(input, s);
    } else {
      std::cerr << "We expect keyword 'solid'" << std::endl;
      return false;
    }

    while(input >> s){
      if(s == endsolid){
        //std::cerr << "found endsolid" << std::endl;
      } else if(s == facet){
        //std::cerr << "found facet" << std::endl;
        std::getline(input, s); // ignore the normal
        input >> s;
        if(s != outer){
          std::cerr << "Expect 'outer' and got " << s << std::endl;
          return false;
        }
        input >> s;
        if(s != loop){
          std::cerr << "Expect 'loop' and got " << s << std::endl;
          return false;
       }
        int count = 0;
        do {
          input >> s;
          if(s == vertex){
            //      std::cerr << "found vertex" << std::endl;
            if(count < 3){
              input >> p;
              if(vertex_index.find(p) == vertex_index.end()){
                ijk[count] = index;
                vertex_index[p] = index++;
                points.push_back(p);
              } else {
                ijk[count] = vertex_index[p];
              }
              ++count;
            } else {
              std::cerr << "We can only read triangulated surfaces" << std::endl;
              return false;
            }
          }
        }while(s != endloop);

        surface.push_back(cpp11::make_tuple(ijk[0], ijk[1], ijk[2]));
      }
    }
    return true;
  }


public:

  std::string name, color;

  Polyhedron_builder_from_STL(std::istream& is_)
    : is(is_), counter(0)
  {}

  void operator()( HDS& hds) {
    if(!read(is, meshPoints, mesh)) return;

    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds);
    B.begin_surface( meshPoints.size(), mesh.size());
    typedef typename Points_3::size_type size_type;

    for(size_type i=0; i < meshPoints.size(); i++){
      B.add_vertex( meshPoints[i]);
    }
    for(size_type i=0; i < mesh.size(); i++){
      B.begin_facet();
      B.add_vertex_to_facet( mesh[i].template get<0>());
      B.add_vertex_to_facet( mesh[i].template get<1>());
      B.add_vertex_to_facet( mesh[i].template get<2>());
      B.end_facet();
    }
    if(B.error())
      {
        std::cerr << "An error occured while creating a Polyhedron" << std::endl;
        B.rollback();
      }

    B.end_surface();
  }
};

} //end of CGAL namespace

#endif // CGAL_IO_POLYHEDRON_STL_BUILDER_H
