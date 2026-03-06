// Copyright (c) 2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :

#ifndef CGAL_IO_FILE_C3T3_IO_H
#define CGAL_IO_FILE_C3T3_IO_H

#include <CGAL/license/SMDS_3.h>

#include <iostream>
#include <string>
#include <limits>
#include <vector>

#include <CGAL/SMDS_3/io_signature.h>

namespace CGAL {

namespace IO {

namespace internal {

  template <typename IOStream>
  void set_binary(IOStream& ios, const bool binary)
  {
    if(binary)
      CGAL::IO::set_binary_mode(ios);
    else
    {
      CGAL::IO::set_ascii_mode(ios);
      ios.precision(std::numeric_limits<double>::max_digits10);
    }
  }

} // end namespace internal


std::vector<std::string> default_c3t3_properties()
{
  return std::vector<std::string>{"cell:subdomain_index",
                                  "facet:surface_patch_index",
                                  "edge:curve_index",
                                  "vertex:corner_index"};
}

  /**
   * @ingroup PkgSMDS3IOFunctions
   * @brief outputs a mesh complex
   *
   * @tparam C3T3 Type of mesh complex, model of `MeshComplex_3InTriangulation_3`
   *
   * @param os the output stream
   * @param c3t3 the mesh complex
   *
   * @sa `CGAL::IO::load_c3t3()`
   *
   * @todo move parameters to named parameters
   */
template <class C3T3>
bool
save_c3t3(std::ostream& os
        , const C3T3& c3t3
        , bool binary = false
        , const std::vector<std::string>& properties = default_c3t3_properties())
{
  os << "CGAL c3t3 " << CGAL::Get_io_signature<C3T3>()() << "\n";
  internal::set_binary(os, binary);

  // write TDS
  // write points
  //// #vertices
  //// for each vertex: point coordinates

  // write cells
  //// #cells
  //// for each cell: 4 vertex indices (use tds.write_cells())

  // write cells properties, matched by names in `properties`
  //// for each property "cell:ppty" in `properties`:
  ////// property name "cell:ppty"
  ////// property type (e.g. int, double, std::string)
  ////// for each cell: index + value of "cell::ppty",
  //////                in the same order as cells are written

  // write facets properties, matched by names in `properties`
  //// # facets
  //// for each property "facet:ppty" in `properties`:
  ////// property name "facet:ppty"
  ////// property type (e.g. int, double, std::string)
  ////// for each facet: value of "facet::ppty",
  //////                in the same order as facets from the iterator

  // write edges properties, matched by names in `properties`
  //// # edges
  //// for each property "edge:ppty" in `properties`:
  ////// property name "edge:ppty"
  ////// property type (e.g. int, double, std::string)
  ////// for each edge: index + value of "edge::ppty",
  //////                in the same order as edges from the iterator

  // write vertices properties, matched by names in `properties`
  //// # vertices
  //// for each property "vertex:ppty" in `properties`:
  ////// property name "vertex:ppty"
  ////// property type (e.g. int, double, std::string)
  ////// for each vertex: index + value of "vertex::ppty",
  //////                in the same order as vertices from the iterator

  return !os.fail();
}

/**
 * @ingroup PkgSMDS3IOFunctions
 * @brief loads a mesh complex
 *
 * @tparam C3T3 Type of mesh complex, model of `MeshComplex_3InTriangulation_3`
 *
 * @param is the input stream
 * @param c3t3 the mesh complex
 *
 * @sa `CGAL::IO::save_c3t3()`
 */
template <class C3T3>
bool load_c3t3(std::istream& is
             , C3T3& c3t3
             , const bool binary = false
             , const std::vector<std::string>& properties = default_c3t3_properties())
{
  internal::set_binary(is, binary);

  return !is.fail();
}

} // end namespace IO

} // end namespace CGAL

#endif // CGAL_IO_FILE_C3T3_IO_H
