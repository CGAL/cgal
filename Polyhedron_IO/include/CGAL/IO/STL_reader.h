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
// Author(s)     : Andreas Fabri,
//                 Mael Rouxel-Labbé

#ifndef CGAL_IO_STL_READER_H
#define CGAL_IO_STL_READER_H

#include <CGAL/array.h>
#include <CGAL/IO/io.h>

#include <boost/cstdint.hpp> 

#include <cctype>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace CGAL {

bool read_ASCII_facet(std::istream& input,
                      std::vector<cpp11::array<double,3> >& points,
                      std::vector<cpp11::array<int,3> >& facets,
                      int& index,
                      std::map<cpp11::array<double, 3>, int >& index_map,
                      bool verbose = false)
{
  // Here, we have already read the word 'facet' and are looking to read till 'endfacet'

  std::string s;
  std::string vertex("vertex"),
              endfacet("endfacet");

  int count = 0;
  cpp11::array<double, 3> p;
  cpp11::array<int, 3> ijk;

  while(input >> s)
  {
    if(s == endfacet)
    {
      if(count != 3)
      {
        if(verbose)
          std::cerr << "Error: only triangulated surfaces are supported" << std::endl;

        return false;
      }

      facets.push_back(ijk);
      return true;
    }
    else if(s == vertex)
    {
      if(count >= 3)
      {
        if(verbose)
          std::cerr << "Error: only triangulated surfaces are supported" << std::endl;

        return false;
      }

      if(!(input >> iformat(p[0]) >> iformat(p[1]) >> iformat(p[2])))
      {
        if(verbose)
          std::cerr << "Error while reading point coordinates (premature end of file)" << std::endl;

        return false;
      }
      else
      {
        std::map<cpp11::array<double, 3>, int>::iterator iti=
            index_map.insert(std::make_pair(p, -1)).first;

        if(iti->second == -1)
        {
          ijk[count] = index;
          iti->second = index++;
          points.push_back(p);
        }
        else
        {
          ijk[count] = iti->second;
        }
      }

      ++count;
    }
  }

  if(verbose)
    std::cerr << "Error while reading facet (premature end of file)" << std::endl;

  return false;
}

bool parse_ASCII_STL(std::istream& input,
                     std::vector<cpp11::array<double,3> >& points,
                     std::vector<cpp11::array<int,3> >& facets,
                     bool verbose = false)
{
  if(verbose)
    std::cout << "Parsing ASCII file..." << std::endl;

  if(!input.good())
    return true;

  // Here, we have already read the word 'solid'

  int index = 0;
  std::map<cpp11::array<double, 3>, int> index_map;

  std::string s;
  std::string facet("facet"),
              endsolid("endsolid");

  while(input >> s)
  {
    if(s == facet)
    {
      if(!read_ASCII_facet(input, points, facets, index, index_map, verbose))
        return false;
    }
    else if(s == endsolid)
    {
      return true;
    }
  }

  if(verbose)
    std::cerr << "Error while parsing ASCII file" << std::endl;

  return false;
}

bool parse_binary_STL(std::istream& input,
                      std::vector<cpp11::array<double,3> >& points,
                      std::vector<cpp11::array<int,3> >& facets,
                      bool verbose = false)
{
  if(verbose)
    std::cout << "Parsing binary file..." << std::endl;

  // Start from the beginning again to simplify things
  input.clear();
  input.seekg(0, std::ios::beg);

  if(!input.good())
    return true;

  // Discard the first 80 chars (unused header)
  int pos = 0;
  char c;

  if(verbose)
    std::cout << "header: ";

  while(pos < 80)
  {
    input.read(reinterpret_cast<char*>(&c), sizeof(c));
    if(!input.good())
      break;

    if(verbose)
      std::cout << c;

    ++pos;
  }

  if(verbose)
    std::cout << std::endl;

  if(pos != 80)
    return true; // empty file

  int index = 0;
  std::map<cpp11::array<double, 3>, int> index_map;

  boost::uint32_t N32;
  if(!(input.read(reinterpret_cast<char*>(&N32), sizeof(N32))))
  {
    if(verbose)
      std::cerr << "Error while reading number of facets" << std::endl;

    return false;
  }

  unsigned int N = N32;
  if(verbose)
    std::cout << N << " facets to read" << std::endl;

  for(unsigned int i=0; i<N; ++i)
  {
    float normal[3];
    if(!(input.read(reinterpret_cast<char*>(&normal[0]), sizeof(normal[0]))) ||
       !(input.read(reinterpret_cast<char*>(&normal[1]), sizeof(normal[1]))) ||
       !(input.read(reinterpret_cast<char*>(&normal[2]), sizeof(normal[2]))))
    {
      if(verbose)
        std::cerr << "Error while reading normal coordinates (premature end of file)" << std::endl;

      return false;
    }

    cpp11::array<int, 3> ijk;

    for(int j=0; j<3; ++j)
    {
      float x,y,z;
      if(!(input.read(reinterpret_cast<char*>(&x), sizeof(x))) ||
         !(input.read(reinterpret_cast<char*>(&y), sizeof(y))) ||
         !(input.read(reinterpret_cast<char*>(&z), sizeof(z))))
      {
        if(verbose)
          std::cerr << "Error while reading vertex coordinates (premature end of file)" << std::endl;

        return false;
      }

      cpp11::array<double, 3> p;
      p[0] = x; p[1] = y; p[2] = z;

      std::map<cpp11::array<double, 3>, int>::iterator iti =
        index_map.insert(std::make_pair(p, -1)).first;

      if(iti->second == -1)
      {
        ijk[j] = index;
        iti->second = index++;
        points.push_back(p);
      }
      else
      {
        ijk[j] = iti->second;
      }
    }

    facets.push_back(ijk);

    // Read so-called attribute byte count and ignore it
    char c;
    if(!(input.read(reinterpret_cast<char*>(&c), sizeof(c))) ||
       !(input.read(reinterpret_cast<char*>(&c), sizeof(c))))
    {
      if(verbose)
        std::cerr << "Error while reading attribute byte count (premature end of file)" << std::endl;

      return false;
    }
  }

  return true;
}

bool read_STL(std::istream& input,
              std::vector<cpp11::array<double,3> >& points,
              std::vector<cpp11::array<int,3> >& facets,
              bool verbose = false)
{
  int pos = 0;

  // Ignore all initial whitespace
  char c;

  while(input.read(reinterpret_cast<char*>(&c), sizeof(c)))
  {
    if(!isspace(c))
    {
      input.unget(); // move back to the first interesting char
      break;
    }
    ++pos;
  }

  if(!input.good()) // reached the end
    return true;

  // If we have gone beyond 80 characters and have not read anything yet,
  // then this must be an ASCII file.
  if(pos > 80)
    return parse_ASCII_STL(input, points, facets, verbose);

  // We are within the first 80 characters, both ASCII and binary are possible

  // Read the 5 first characters to check if the first word is "solid"
  std::string s, solid("solid");

  char word[5];
  if(input.read(reinterpret_cast<char*>(&word[0]), sizeof(c)) &&
     input.read(reinterpret_cast<char*>(&word[1]), sizeof(c)) &&
     input.read(reinterpret_cast<char*>(&word[2]), sizeof(c)) &&
     input.read(reinterpret_cast<char*>(&word[3]), sizeof(c)) &&
     input.read(reinterpret_cast<char*>(&word[4]), sizeof(c)))
  {
    s = std::string(word, 5);
    pos += 5;
  }
  else
    return true; // empty file

  // If the first word is not 'solid', the file must be binary
  if(s != solid)
  {
    if(parse_binary_STL(input, points, facets, verbose))
    {
      return true;
    }
    else
    {
      // If we failed to read it as a binary, try as ASCII just in case...
      // The file does not start with 'solid' anyway, so it's fine to reset it.
      input.clear();
      input.seekg(0, std::ios::beg);
      return parse_ASCII_STL(input, points, facets, verbose);
    }
  }

  // Now, we have found the keyword "solid" which is supposed to indicate that the file is ASCII
  if(parse_ASCII_STL(input, points, facets, verbose))
  {
    // correctly read the input as an ASCII file
    return true;
  }
  else // Failed to read the ASCII file
  {
    // It might have actually have been a binary file... ?
    return parse_binary_STL(input, points, facets, verbose);
  }
}

} // namespace CGAL

#endif // CGAL_IO_STL_READER_H
