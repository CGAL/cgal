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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_IO_STL_READER_H
#define CGAL_IO_STL_READER_H

#include <CGAL/array.h>
#include <boost/cstdint.hpp> 
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <cctype>

namespace CGAL{

  bool
  read_STL( std::istream& input,
            std::vector< cpp11::array<double,3> >& points,
            std::vector< cpp11::array<int,3> >& facets,
            bool verbose = false)
  {
    bool is_binary_file = false;
    
    std::string s, solid("solid"), facet("facet");
    std::map<cpp11::array<double,3>, int> pmap;
    int index = 0;
    cpp11::array<int,3> ijk;
    cpp11::array<double,3> p;

    char line[80];
    int i = 0, ni = 0;

      // read the 5 first characters to check if the first word is "solid"
      boost::uint8_t c;
      for(; i < 5; i++){
        input.read(reinterpret_cast<char*>(&c), sizeof(c));
        line[i]=c;
      }

      s = std::string(line,5);
      if(s == solid){
        // we found the keyword "solid" which is supposed to indicate the file is Ascii
        // But it might still be binary, so we have to find out if it is followed
        // by an (optional) name and then the keyword "facet"
        // When we find "facet" we conclude that it is really Ascii

        // first skip whitespaces after solid
        do {
          input.read(reinterpret_cast<char*>(&c), sizeof(c));
          line[i++]=c;
        }while(isspace(c) && ( i < 80));
        if(i==80){
          is_binary_file = true;
          goto done;
        }
        // now c is not whitespace
        ni = i-1; // here starts either the name or the keyword "facet"
        do {
          input.read(reinterpret_cast<char*>(&c), sizeof(c));
          line[i++]=c;
        }while(! isspace(c) && ( i < 80));
        s = std::string(line+ni, (i-1) - ni);
#       ifdef CGAL_DEBUG_BINARY_HEADER
          std::cout << "|" << s  << "|" << std::endl;
#       endif        
        if(s == facet){
          goto done;
        } else if(i == 80){
          // the entire header is a name
          is_binary_file = true;
          goto done;
        }
        
        // we continue to read what comes after the name
        
        // now c is whitespace, skip other whitespaces
        do {
          input.read(reinterpret_cast<char*>(&c), sizeof(c));
          line[i++]=c;
        }while(isspace(c) && ( i < 80));
        if(i==80){
          is_binary_file = true;
          goto done;
        }
        
        // now c is not whitespace
        ni = i-1; // here starts either "facet", or it is really binary
        do {
          input.read(reinterpret_cast<char*>(&c), sizeof(c));
          line[i++]=c;
        }while(! isspace(c) && ( i < 80));
          s = std::string(line+ni, (i-1) - ni);
#       ifdef CGAL_DEBUG_BINARY_HEADER
          std::cout << "|" << s  << "|" << std::endl;
#       endif
        if(s == facet){
          goto done;
        } else {
          for(; i < 80; i++){
            input.read(reinterpret_cast<char*>(&c), sizeof(c));
          }
          is_binary_file = true;
          goto done;
        }
      }else{
        // we read the other 75 characters of the header
        for(; i < 80; i++){
          input.read(reinterpret_cast<char*>(&c), sizeof(c));
        }
        is_binary_file = true;
      }

  done:

    if(is_binary_file){
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

      // It is an Ascii file but the first occurence of "facet" has already be parsed
      bool first_facet = true;
      std::string outer("outer"),
        loop("loop"),
        vertex("vertex"),
        endloop("endloop"),
        endsolid("endsolid");
      s = facet;
      while(first_facet || (input >> s)){
        first_facet = false;
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
