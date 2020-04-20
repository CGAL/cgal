// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_STL_H
#define CGAL_IO_STL_H

#include <CGAL/IO/STL/STL_reader.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <boost/range/value_type.hpp>

#include <iostream>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

/*!
 * \ingroup IOstreamFunctions
 *
 * reads the content of `is` into `points` and `facets`, in the STL format.
 *
 * \see \ref IOStreamSTL
 */
template <typename PointRange, typename TriangleRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_STL(std::istream& is,
              PointRange& points,
              TriangleRange& facets,
              const CGAL_BGL_NP_CLASS& /*np*/) // might become useful one day for face normals
{
  const bool verbose = false;
  int pos = 0;

  // Ignore all initial whitespace
  unsigned char c;

  while(is.read(reinterpret_cast<char*>(&c), sizeof(c)))
  {
    if(!isspace(c))
    {
      is.unget(); // move back to the first interesting char
      break;
    }
    ++pos;
  }

  if(!is.good()) // reached the end
    return true;

  // If we have gone beyond 80 characters and have not read anything yet,
  // then this must be an ASCII file.
  if(pos > 80)
    return IO::internal::parse_ASCII_STL(is, points, facets, verbose);

  // We are within the first 80 characters, both ASCII and binary are possible

  // Read the 5 first characters to check if the first word is "solid"
  std::string s;

  char word[6];
  if(is.read(reinterpret_cast<char*>(&word[0]), sizeof(c)) &&
     is.read(reinterpret_cast<char*>(&word[1]), sizeof(c)) &&
     is.read(reinterpret_cast<char*>(&word[2]), sizeof(c)) &&
     is.read(reinterpret_cast<char*>(&word[3]), sizeof(c)) &&
     is.read(reinterpret_cast<char*>(&word[4]), sizeof(c)) &&
     is.read(reinterpret_cast<char*>(&word[5]), sizeof(c)))
  {
    s = std::string(word, 5);
    pos += 5;
  }
  else
  {
    return true; // empty file
  }

  // If the first word is not 'solid', the file must be binary
  if(s != "solid" || (word[5] !='\n' && word[5] != ' '))
  {
    if(IO::internal::parse_binary_STL(is, points, facets, verbose))
    {
      return true;
    }
    else
    {
      // If we failed to read it as a binary, try as ASCII just in case...
      // The file does not start with 'solid' anyway, so it's fine to reset it.
      is.clear();
      is.seekg(0, std::ios::beg);
      return IO::internal::parse_ASCII_STL(is, points, facets, verbose);
    }
  }

  // Now, we have found the keyword "solid" which is supposed to indicate that the file is ASCII
  is.clear();
  is.seekg(0, std::ios::beg); // the parser needs to read all "solid" to work correctly.
  if(IO::internal::parse_ASCII_STL(is, points, facets, verbose))
  {
    // correctly read the input as an ASCII file
    return true;
  }
  else // Failed to read the ASCII file
  {
    // It might have actually have been a binary file... ?
    return IO::internal::parse_binary_STL(is, points, facets, verbose);
  }
}

//todo doc this too
template <typename PointRange, typename TriangleRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_STL(const char* fname,
              PointRange& points,
              TriangleRange& facets,
              const CGAL_BGL_NP_CLASS& np)
{
  std::ifstream in(fname);
  return read_STL(in, points, facets, np);
}

template <typename PointRange, typename TriangleRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_STL(const std::string& fname, PointRange& points, TriangleRange& facets, const CGAL_BGL_NP_CLASS& np)
{
  return read_STL(fname.c_str(), points, facets, np);
}

template <typename PointRange, typename TriangleRange>
bool read_STL(std::istream& is, PointRange& points, TriangleRange& facets)
{
  return read_STL(is, points, facets, parameters::all_default());
}

template <typename PointRange, typename TriangleRange>
bool read_STL(const char* fname, PointRange& points, TriangleRange& facets)
{
  return read_STL(fname, points, facets, parameters::all_default());
}

template <typename PointRange, typename TriangleRange>
bool read_STL(const std::string& fname, PointRange& points, TriangleRange& facets)
{
  return read_STL(fname, points, facets, parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
 * \ingroup IOstreamFunctions
 *
 * writes the content of `points` and `facets` in `out`, in the STL format.
 *
 * \see \ref IOStreamSTL
 */
template <class PointRange, class TriangleRange>
std::ostream& write_STL(std::ostream& out,
                        const PointRange& points,
                        const TriangleRange& facets)
{
  typedef typename boost::range_value<TriangleRange>::type            Triangle;
  typedef typename boost::range_value<PointRange>::type               Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel                 K;
  typedef typename K::Vector_3                                        Vector_3;

  if (get_mode(out) == IO::BINARY)
  {
    out << "FileType: Binary                                                                ";
    const boost::uint32_t N32 = static_cast<boost::uint32_t>(facets.size());
    out.write(reinterpret_cast<const char *>(&N32), sizeof(N32));

    for(const Triangle& face : facets)
    {
      const Point& p = points[face[0]];
      const Point& q = points[face[1]];
      const Point& r = points[face[2]];

      const Vector_3 n = collinear(p,q,r) ? Vector_3(1,0,0) : unit_normal(p,q,r);

      const float coords[12] = { static_cast<float>(n.x()), static_cast<float>(n.y()), static_cast<float>(n.z()),
                                 static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()),
                                 static_cast<float>(q.x()), static_cast<float>(q.y()), static_cast<float>(q.z()),
                                 static_cast<float>(r.x()), static_cast<float>(r.y()), static_cast<float>(r.z()) };

      for (int i=0; i<12; ++i)
        out.write(reinterpret_cast<const char *>(&coords[i]), sizeof(coords[i]));
      out << "  ";
    }
  }
  else
  {
    out << "solid\n";
    for(const Triangle& face : facets)
    {
      const Point& p = points[face[0]];
      const Point& q = points[face[1]];
      const Point& r = points[face[2]];

      const Vector_3 n = collinear(p,q,r) ? Vector_3(1,0,0) : unit_normal(p,q,r);
      out << "facet normal " << n << "\nouter loop\n";
      out << "vertex " << p << "\n";
      out << "vertex " << q << "\n";
      out << "vertex " << r << "\n";
      out << "endloop\nendfacet\n";
    }
    out << "endsolid\n";
  }

  return out;
}

template <typename PointRange, typename TriangleRange>
bool write_STL(const char* fname, const PointRange& points, const TriangleRange& facets)
{
  std::ofstream out(fname);
  return write_STL(out, points, facets);
}

template <typename PointRange, typename TriangleRange>
bool write_STL(const std::string& fname, const PointRange& points, const TriangleRange& facets)
{
  return write_STL(fname.c_str(), points, facets);
}

} // namespace CGAL

#endif // CGAL_IO_STL_H
