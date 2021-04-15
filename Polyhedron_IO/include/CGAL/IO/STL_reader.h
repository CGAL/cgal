// Copyright (c) 2015 GeometryFactory
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri,
//                 Mael Rouxel-Labbé

#ifndef CGAL_IO_STL_READER_H
#define CGAL_IO_STL_READER_H

#include <CGAL/IO/io.h>
#include <CGAL/IO/reader_helpers.h>

#include <CGAL/Container_helper.h>

#include <boost/cstdint.hpp>

#include <cctype>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace CGAL {

template <class Point, class Triangle>
bool read_ASCII_facet(std::istream& input,
                      std::vector<Point>& points,
                      std::vector<Triangle>& facets,
                      int& index,
                      std::map<Point, int >& index_map,
                      bool verbose = false)
{
  // Here, we have already read the word 'facet' and are looking to read till 'endfacet'

  std::string s;
  std::string vertex("vertex"),
              endfacet("endfacet");

  int count = 0;
  double x,y,z;
  Point p;
  Triangle ijk;
  CGAL::internal::resize(ijk, 3);

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

      if(!(input >> iformat(x) >> iformat(y) >> iformat(z)))
      {
        if(verbose)
          std::cerr << "Error while reading point coordinates (premature end of file)" << std::endl;

        return false;
      }
      else
      {
        IO::internal::fill_point(x, y, z, p);
        typename std::map<Point, int>::iterator iti = index_map.insert(std::make_pair(p, -1)).first;

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

template <class Point, class Triangle>
bool parse_ASCII_STL(std::istream& input,
                     std::vector<Point>& points,
                     std::vector<Triangle>& facets,
                     bool verbose = false)
{
  if(verbose)
    std::cout << "Parsing ASCII file..." << std::endl;

  if(!input.good())
    return true;

  // Here, we have already read the word 'solid'

  int index = 0;
  std::map<Point, int> index_map;

  std::string s, facet("facet"), endsolid("endsolid");

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

template <typename Point, typename Triangle>
bool parse_binary_STL(std::istream& input,
                      std::vector<Point>& points,
                      std::vector<Triangle>& facets,
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
  std::map<Point, int> index_map;

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

    Triangle ijk;
    CGAL::internal::resize(ijk, 3);

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

      Point p;
      IO::internal::fill_point(x, y, z, p);

      typename std::map<Point, int>::iterator iti = index_map.insert(std::make_pair(p, -1)).first;

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

//
// Read a file with `.stl` format.
//
// \tparam Point must be a model of the concept `RandomAccessContainer` or a %CGAL point type
// \tparam Triangle must be a model of the concept `RandomAccessContainer`
//
// \param input the input stream
// \param points a container that will contain the points used in the .stl file
// \param polygons a container that will contain the triangles used in the .stl file
// \param verbose whether to enable or not a sanity log
//
// \returns `true` if the reading process went well, `false` otherwise
//
// \warning `points` and `facets` are not cleared: new points and triangles are added to the back
//          of the containers.
//
// Although the STL file format uses triangles, it is convenient to be able to use vectors
// and other models of the `SequenceContainer` (instead of arrays) for the face type,
// to avoid having to convert the to apply polygon soup reparation algorithms.
template <class Point, class Triangle>
bool read_STL(std::istream& input,
              std::vector<Point>& points,
              std::vector<Triangle>& facets,
              bool verbose = false)
{
  int pos = 0;

  // Ignore all initial whitespace
  unsigned char c;

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
