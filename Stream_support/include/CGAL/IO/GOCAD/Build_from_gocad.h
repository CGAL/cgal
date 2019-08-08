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

#ifndef CGAL_IO_GOCAD_BUILD_FROM_GOCAD_H
#define CGAL_IO_GOCAD_BUILD_FROM_GOCAD_H

#include <deque>
#include <iostream>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <iostream>

namespace CGAL{

namespace GOCAD_internal {
template <class Facegraph, class P, class Derived>
class Build_from_gocad
{
protected:
  typedef P Point_3;
  typedef std::deque<Point_3> Points_3;
  typedef boost::tuple<int,int,int> Facet;
  typedef std::deque<Facet> Surface;

public:
  std::string name, color;
  Build_from_gocad(std::istream& is_)
    : is(is_), counter(0)
  {}

  void do_construct(){} //specific to Facegraph (declared in BGL)
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

  void operator()( Facegraph& graph)
  {
    read(this->is, this->meshPoints, this->mesh);
    static_cast<Derived*>(this)->do_construct(graph);
  }

protected:
  std::istream& is;
  int counter;
  Points_3 meshPoints;
  Surface mesh;
};



}//end GOCAD_internal
}//end CGAL
#endif // CGAL_IO_GOCAD_BUILD_FROM_GOCAD_H
 
