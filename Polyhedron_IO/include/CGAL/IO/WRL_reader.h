// Copyright (c) 2018 GeometryFactory (France).
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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Sébastien Loriot,
//                 Mael Rouxel-Labbé
//

#ifndef CGAL_IO_WRL_READER_H
#define CGAL_IO_WRL_READER_H

#include <CGAL/assertions.h>
#include <CGAL/array.h>
#include <CGAL/IO/io.h>
#include <CGAL/use.h>

#include <boost/foreach.hpp>
#include <boost/cstdint.hpp>

#include <cctype>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace CGAL {

namespace IO {

namespace internal {

void discard_whitespace(std::istream& input)
{
  if(!input.good())
    return;

  char peeked_c = input.peek();
  while(std::isspace(static_cast<unsigned char>(peeked_c)))
  {
    peeked_c = input.get();
    if(!input.good())
      return;
    peeked_c = input.peek();
  }
}

// Checks if we are encountering 'c' next, ignoring whitespace
bool is_next_character(std::istream& input,
                       char c,
                       const bool verbose = false)
{
  input >> std::ws; // discards potential whitespace

  char peeked_c = input.peek();
  if(verbose)
    std::cout << "Peeking... " << peeked_c << ", c was: " << c << std::endl;

  return (peeked_c == c);
}

// reading "1324]" as an int goes fine, but reading "[12548" as an int is no-no,
// so we gotta discard the bracket safely first incase it touches the int or w/e is after
void read_opening_bracket(std::istream& input)
{
  CGAL_assertion(is_next_character(input, '['));

  discard_whitespace(input);
  char bracket = input.get();
  CGAL_USE(bracket);
  CGAL_assertion(bracket == '[');
}

template <typename Point, typename Face>
void read_WRL_single_shape(std::istream& input,
                           std::vector<Point>& points,
                           std::vector<Face>& facets,
                           const bool verbose = false)
{
  std::string s;
  double x, y, z;
  int i;

  while(input >> s)
  {
    if(verbose)
      std::cout << "Reading... " << s << std::endl;

    if(s == "point")
    {
      read_opening_bracket(input);

      if(is_next_character(input, ']'))
        continue;

      bool met_closing_bracket = false;
      while(input >> x >> y >> z >> s)
      {
        CGAL_assertion(s == ",");
        if(verbose)
          std::cout << "Read point: " << x << " " << y << " " << z << std::endl;
        points.push_back(CGAL::make_array(x, y, z));

        if(is_next_character(input, ']'))
        {
          met_closing_bracket = true;
          break;
        }
      }
      CGAL_assertion(met_closing_bracket);
    }
    else if(s == "coordIndex")
    {
      read_opening_bracket(input);

      if(is_next_character(input, ']'))
        continue;

      Face facet;
      bool met_closing_bracket = false;
      while(input >> i)
      {
        if(i == -1)
        {
          if(verbose)
            std::cout << "Add facet of size: " << facet.size() << std::endl;

          facets.push_back(facet);
          facet.clear();
        }
        else
        {
          facet.push_back(i);
        }

        // ignore the potential spaces / commas
        if(input >> std::ws && input.peek() == ',')
          input.ignore();

        if(is_next_character(input, ']'))
        {
          met_closing_bracket = true;
          break;
        }
      }

      CGAL_assertion(met_closing_bracket);
      break; // finished reading facets means we are done reading the shape
    }
  }
}

  int digit_counter(std::size_t i)
{
  CGAL_precondition(i >= 0);
  int n = 1;
  while(i/=10)
    ++n;
  return n;
}

// Merge the polygon soup given by 'points' and 'faces' into the cumulated 'all_...'. Identical
// points are NOT merged along the way.
template <typename Point, typename Face>
void merge_polygon_soups(const std::vector<Point>& points,
                         const std::vector<Face>& faces,
                         std::vector<Point>& all_points,
                         std::vector<Face>& all_faces)
{
  std::size_t nv = points.size(), nf = faces.size(), apn = all_points.size();

  all_points.reserve(apn + nv);
  all_faces.reserve(all_faces.size() + nf);

  std::copy(points.begin(), points.end(), std::back_inserter(all_points));

  for(std::size_t fi=0; fi<nf; ++fi)
  {
    const std::size_t facet_size = faces[fi].size();
    Face face;
    face.reserve(facet_size);

    for(std::size_t sfi=0; sfi<facet_size; ++sfi)
    {
      int new_id = apn + faces.at(fi).at(sfi);
      face.push_back(new_id);
    }

    all_faces.push_back(face);
  }
}

// The function below does NOT merge the shapes together
template <typename Point, typename Face>
void write_shapes_to_OFF(const std::string filename,
                         const std::vector<std::vector<Point> >& points,
                         const std::vector<std::vector<Face> >& meshes,
                         const bool verbose = false)
{
  std::size_t number_of_shapes = points.size();

  CGAL_precondition(points.size() == meshes.size());
  const int d = digit_counter(number_of_shapes);
  CGAL_USE(d);

  for(std::size_t i=0; i<number_of_shapes; ++i)
  {
    const std::vector<Point>& single_shape_points = points[i];
    const std::vector<Face>& single_shape_facets = meshes[i];

    std::stringstream oss;
    if(meshes.size() > 1)
      oss << filename << "_" << std::setfill('0') << std::setw(d) << i+1 << ".off" << std::ends;
    else
      oss << filename << ".off" << std::ends;

    if(verbose)
      std::cout << "Writing (OFF) to: " << oss.str() << std::endl;

    std::ofstream out(oss.str().c_str());

    out << "OFF\n" << single_shape_points.size() << " " << single_shape_facets.size() << " 0" << '\n';

    for(std::size_t pi=0; pi<single_shape_points.size(); ++pi)
      out << single_shape_points[pi][0] << " "
          << single_shape_points[pi][1] << " "
          << single_shape_points[pi][2] << '\n';

    for(std::size_t fi=0; fi<single_shape_facets.size(); ++fi)
    {
      std::size_t facet_size = single_shape_facets[fi].size();
      out << facet_size << " ";
      for(std::size_t sfi=0; sfi<facet_size; ++sfi)
        out << single_shape_facets[fi][sfi] << " ";
      out << "\n";
    }
  }
}

// @todo could probably be merged with the function below
template <typename Point, typename Face>
bool read_WRL_and_merge_shapes(std::istream& input,
                               std::vector<Point>& all_points,
                               std::vector<Face>& all_faces,
                               const bool verbose = false)
{
  std::string s;

  while(input >> s)
  {
    // Ignore everything until 'Shape' is read
    if(s == "Shape")
    {
      std::vector<Point> single_shape_points;
      std::vector<Face> single_shape_facets;
      read_WRL_single_shape(input, single_shape_points, single_shape_facets, verbose);
      merge_polygon_soups(single_shape_points, single_shape_facets, all_points, all_faces);
      if(verbose)
        std::cout << "Read... " << single_shape_facets.size() << " facets" << std::endl;
    }
  }

  if(verbose)
    std::cout << "Read " << all_faces.size() << " facets (total)" << std::endl;

  return true;
}

// Very rough, nowhere near covering the full standard
template <typename Point, typename Face>
bool read_WRL(std::istream& input,
              std::vector<std::vector<Point> >& points,
              std::vector<std::vector<Face> >& meshes,
              const bool verbose = false)
{
  std::string s;

  while(input >> s)
  {
    // Ignore everything until 'Shape' is read
    if(s == "Shape")
    {
      std::vector<Point> single_shape_points;
      std::vector<Face> single_shape_facets;
      read_WRL_single_shape(input, single_shape_points, single_shape_facets, verbose);
      points.push_back(single_shape_points);
      meshes.push_back(single_shape_facets);
      if(verbose)
        std::cout << "Read... " << single_shape_facets.size() << " facets" << std::endl;
    }
  }

  if(verbose)
    std::cout << "Read " << meshes.size() << " shapes" << std::endl;

  return true;
}

bool convert_WRL_to_OFF(const char* filename)
{
  std::cout << "convert: " << filename << std::endl;

  std::vector<std::vector<CGAL::cpp11::array<double, 3> > > points;
  std::vector<std::vector<std::vector<int> > > shapes;

  std::ifstream in(filename);
  if(!in.good() || !read_WRL(in, points, shapes, false))
  {
    std::cerr << "Could not read file: " << filename << std::endl;
    return false;
  }

  std::string filename_core(filename);
  filename_core = filename_core.substr(0, filename_core.find_last_of('.'));
  write_shapes_to_OFF(filename_core, points, shapes);

  return true;
}

} // end namespace internal

} // end namespace IO

} // end namespace CGAL

#endif // CGAL_IO_WRL_READER_H
