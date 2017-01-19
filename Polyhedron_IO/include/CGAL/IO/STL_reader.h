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

#ifndef CGAL_IO_STL_READER_H
#define CGAL_IO_STL_READER_H

#include <CGAL/license/Polyhedron.h>


#include <CGAL/array.h>
#include <boost/cstdint.hpp> 
#include <vector>
#include <map>
#include <iostream>
#include <string>

namespace CGAL{

  bool
  read_STL( std::istream& input,
            std::vector< cpp11::array<double,3> >& points,
            std::vector< cpp11::array<int,3> >& facets,
            bool verbose = false)
  {
    std::string s, solid("solid");
    std::map<cpp11::array<double,3>, int> pmap;
    int index = 0;
    cpp11::array<int,3> ijk;
    cpp11::array<double,3> p;

    char line[80];
    for(int i=0;i < 80; i++){
      boost::uint8_t c;
      input.read(reinterpret_cast<char*>(&c), sizeof(c));
      line[i]=c;
      if(i==5){
        s = std::string(line,5);
        if(s == solid){
          break;
        }
      }
    }

    if(s!= solid){
      boost::uint32_t N32;
      input.read(reinterpret_cast<char*>(&N32), sizeof(N32));
      unsigned int N = N32;

      for(unsigned int i=0; i < N; i++){
        float normal[3];
        input.read(reinterpret_cast<char*>(&normal[0]), sizeof(normal[0]));
        input.read(reinterpret_cast<char*>(&normal[1]), sizeof(normal[1]));
        input.read(reinterpret_cast<char*>(&normal[2]), sizeof(normal[2]));

        for(int j=0; j < 3; j++){
          float x,y,z;
          input.read(reinterpret_cast<char*>(&x), sizeof(x));
          input.read(reinterpret_cast<char*>(&y), sizeof(y));
          input.read(reinterpret_cast<char*>(&z), sizeof(z));
          p[0]=x; p[1]=y; p[2]=z;
          std::map<cpp11::array<double,3>, int>::iterator iti=
            pmap.insert(std::make_pair(p,-1)).first;
          if(iti->second==-1){
            ijk[j] = index;
            iti->second = index++;
            points.push_back(p);
          } else {
            ijk[j] = iti->second;
          }
        }
        if((ijk[0] != ijk[1]) && 
           (ijk[0] != ijk[2]) &&
           (ijk[1] != ijk[2])){
          facets.push_back(ijk);
        }else{
          if(verbose){
            std::cerr << "ignore degenerate face" << std::endl;
          }
        }
        char c;
        input.read(reinterpret_cast<char*>(&c), sizeof(c));
        input.read(reinterpret_cast<char*>(&c), sizeof(c));
      }
      return true;
    } else {
      std::string facet("facet"),
        outer("outer"),
        loop("loop"),
        vertex("vertex"),
        endloop("endloop"),
        endsolid("endsolid");

      while(input >> s){
        if(s == endsolid){
          //std::cerr << "found endsolid" << std::endl;
        } else if(s == facet){
          //std::cerr << "found facet" << std::endl;
          std::getline(input, s); // ignore the normal
          input >> s;
          if(s != outer){
            if (verbose)
              std::cerr << "Expect 'outer' and got " << s << std::endl;
            return false;
          }
          input >> s;
          if(s != loop){
            if (verbose)
              std::cerr << "Expect 'loop' and got " << s << std::endl;
            return false;
          }
          int count = 0;
          do {
            input >> s;
            if(s == vertex){
              //      std::cerr << "found vertex" << std::endl;
              if(count < 3){
                input >> p[0] >> p[1] >> p[2];
                std::map<cpp11::array<double,3>, int>::iterator iti=
                  pmap.insert(std::make_pair(p,-1)).first;
                if(iti->second==-1){
                  ijk[count] = index;
                  iti->second = index++;
                  points.push_back(p);
                } else {
                  ijk[count] = iti->second;
                }
                ++count;
              } else {
                if (verbose)
                  std::cerr << "We can only read triangulated surfaces" << std::endl;
                return false;
              }
            }
          }while(s != endloop);
          
          facets.push_back(ijk);
        }
      }
      return true;
    }
  }
} // namespace CGAL

#endif // CGAL_IO_STL_READER_H
