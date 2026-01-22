// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_IO_PLY_READ_XYZ_POINTS_H
#define CGAL_IO_PLY_READ_XYZ_POINTS_H

#include <CGAL/IO/XYZ.h>
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/type_traits/is_iterator.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace CGAL {

namespace IO {

// doxygen in ../XYZ.h
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
bool read_XYZ(std::istream& is,
              OutputIterator output,
              const CGAL_NP_CLASS& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3::Fake_point_range<OutputIteratorValueType> PointRange;

  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, CGAL_NP_CLASS> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;
  typedef typename NP_helper::Geom_traits Kernel;

  PointMap point_map = NP_helper::get_point_map(np);
  NormalMap normal_map = NP_helper::get_normal_map(np);

  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  //typedef typename value_type_traits<OutputIterator>::type Enriched_point;
  typedef OutputIteratorValueType Enriched_point;

  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;

  if(!is)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // scan points
  long pointsCount; // number of points in file
  int lineNumber = 0; // line counter
  std::string line; // line buffer
  std::istringstream iss;

  while(getline(is,line))
  {
    // position + normal
    FT x,y,z;
    FT nx,ny,nz;

    ++lineNumber;

    // Trims line buffer
    line.erase(line.find_last_not_of (" ")+1);
    line.erase(0, line.find_first_not_of (" "));

    // Skips comment or empty line...
    if (line.length() == 0 || line[0] == '#')
    {
      continue;
    }
    // ...or reads position...
    else
    {
      iss.clear();
      iss.str(line);
      if (iss >> IO::iformat(x) >> IO::iformat(y) >> IO::iformat(z))
      {
        Point point(x,y,z);
        Vector normal = CGAL::NULL_VECTOR;
        // ... + normal...
        if (iss >> IO::iformat(nx))
        {
          // In case we could read one number, we expect that there are two more
          if(iss >> IO::iformat(ny) >> IO::iformat(nz)){
            normal = Vector(nx,ny,nz);
          } else {
            std::cerr << "Error line " << lineNumber << " of file (incomplete normal coordinates)" << std::endl;
            return false;
          }
        }

        Enriched_point pwn;
        put(point_map,  pwn, point);  // point_map[pwn] = point

        put(normal_map, pwn, normal); // normal_map[pwn] = normal

        *output++ = pwn;
        continue;
      }
    }

    // ...or skips number of points on first line (optional)
    if (lineNumber == 1 && std::istringstream(line) >> pointsCount)
    {
      continue;
    }
    else // if wrong file format
    {
      std::cerr << "Error line " << lineNumber << " of file (expected point coordinates)" << std::endl;
      return false;
    }
  }

  if(is.eof())
    is.clear(is.rdstate() & ~std::ios_base::failbit); // set by getline

  return true;
}

// documentation in ../XYZ.h
template <typename OutputIteratorValueType,
          typename OutputIterator,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
bool read_XYZ(const std::string& fname,
              OutputIterator output,
              const CGAL_NP_CLASS& np)
{
  std::ifstream is(fname);
  return read_XYZ<OutputIteratorValueType>(is, output, np);
}

/// \cond SKIP_IN_MANUAL

// variants with default output iterator value type
template <typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(std::istream& is, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values(),
              std::enable_if_t<CGAL::is_iterator<OutputIterator>::value>* = nullptr)
{
  return read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

template <typename OutputIterator,typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(const std::string& fname, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values())
{
  std::ifstream is(fname);
  return read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

/// \endcond

} // namespace IO



} // namespace CGAL

#endif // CGAL_IO_PLY_READ_XYZ_POINTS_H
