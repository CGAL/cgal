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

#include <CGAL/Unique_hash_map.h>
#include <CGAL/IO/internal/C3t3_io_helpers.h>
#include <CGAL/IO/io.h>

#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <variant>

#include <CGAL/SMDS_3/io_signature.h>


namespace CGAL {

namespace IO {

  /**
   * @ingroup PkgSMDS3IOFunctions
   * @brief outputs a mesh complex
   *
   * @tparam C3t3 Type of mesh complex, model of `MeshComplex_3InTriangulation_3`
   * @tparam Property //todo
   *
   * @param os the output stream
   * @param c3t3 the mesh complex
   *
   * @sa `CGAL::IO::load_c3t3()`
   *
   * @todo move parameters to named parameters
   * @todo use cell indices instead of pairs/triples/quadruples of vertex indices to
   *  match properties with cells/facets/edges
   */
template <class C3t3, typename Properties>
bool
save_c3t3(std::ostream& os
        , const C3t3& c3t3
        , bool binary
        , const Properties& properties)
{
  // usage example of the named parameters
  // CGAL::parameters::vertex_properties(Circum_getter())
  using Property = typename Properties::value_type;
  using Tr = typename C3t3::Triangulation;
  using Vertex_handle = typename Tr::Vertex_handle;

  const Tr& tr = c3t3.triangulation();

  os << "CGAL c3t3 " << CGAL::Get_io_signature<C3t3>()() << "\n";
  if(binary)
    os << "binary";
  internal::set_binary(os, binary);

  // write points
  //// #vertices
  os << "Vertices" << "\n";
  os << tr.number_of_vertices() << "\n";
  //// points coordinates, in the same order as vertices are written
  CGAL::Unique_hash_map<Vertex_handle, std::size_t> vertices;
  std::size_t i = 0;
  vertices[tr.infinite_vertex()] = i++;
  for(const auto v : tr.finite_vertex_handles())
  {
    vertices[v] = i++;
    os << tr.point(v) << "\n";
  }

  // write cells
  //// #cells
  os << "Cells" << "\n";
  //// write the nb of cells, and
  //// for each cell: 4 vertex indices
  tr.tds().print_cells(os, vertices);

  // write cells properties, matched by names in `properties`
  //// for each property "cell:ppty" in `properties`:
  for(const auto& pp : properties)
  {
    std::visit([&](auto& ppty)
    {
      os << "# " << ppty.type() << " "
                 << std::to_string(ppty.number_of_simplices(c3t3)) << " "
                 << std::to_string(ppty.number_of_components()) << " "
                 << ppty.name() << "\n";

      using S = typename std::decay_t<decltype(ppty)>::key_type;
      for(const S& c : internal::get_simplex_range<S>(c3t3))
      {
        if(!ppty.has_value(c, c3t3))
          continue;
        internal::write_property(os, c, c3t3, vertices);
        os << ppty(c, c3t3);
        if(!binary)
          os << "\n";
      }
      if(!binary)
        os << "\n";
    }, pp);
  }

  ////// property type (e.g. int, double, std::string)
  ////// property name "cell:ppty"
  ////// for each cell: index + value of "cell::ppty",
  //////                in the same order as cells are written

  // write facets properties, matched by names in `properties`
  //// # facets
  //// for each property "facet:ppty" in `properties`:
  ////// property type (e.g. int, double, std::string)
  ////// property name "facet:ppty"
  ////// for each facet: value of "facet::ppty",
  //////                in the same order as facets from the iterator

  // write edges properties, matched by names in `properties`
  //// # edges
  //// for each property "edge:ppty" in `properties`:
  ////// property type (e.g. int, double, std::string)
  ////// property name "edge:ppty"
  ////// for each edge: index + value of "edge::ppty",
  //////                in the same order as edges from the iterator

  // write vertices properties, matched by names in `properties`
  //// # vertices
  //// for each property "vertex:ppty" in `properties`:
  ////// property type (e.g. int, double, std::string)
  ////// property name "vertex:ppty"
  ////// for each vertex: index + value of "vertex::ppty",
  //////                in the same order as vertices from the iterator

  return !os.fail();
}

template <class C3t3>
bool
save_c3t3(std::ostream& os
        , const C3t3& c3t3
        , bool binary = false)
{
  using P = std::variant<internal::Subdomain_index_property<C3t3>,
                         internal::Surface_patch_index_property<C3t3>,
                         internal::Curve_index_property<C3t3>,
                         internal::Corner_index_property<C3t3>>;
  std::vector<P> c3t3_properties{internal::Subdomain_index_property<C3t3>(),
                                 internal::Surface_patch_index_property<C3t3>(),
                                 internal::Curve_index_property<C3t3>(),
                                 internal::Corner_index_property<C3t3>()};
  return save_c3t3(os, c3t3, binary, c3t3_properties);
}

/**
 * @ingroup PkgSMDS3IOFunctions
 * @brief loads a mesh complex
 *
 * @tparam C3t3 Type of mesh complex, model of `MeshComplex_3InTriangulation_3`
 *
 * @param is the input stream
 * @param c3t3 the mesh complex
 *
 * @sa `CGAL::IO::save_c3t3()`
 */
template <class C3t3> //, typename Properties>
bool load_c3t3(std::istream& is
             , C3t3& c3t3
             , const bool binary)
{
  using Tr = typename C3t3::Triangulation;
  using Cell_handle = typename Tr::Cell_handle;
  using Facet = typename Tr::Facet;
  using Edge = typename Tr::Edge;
  using Vertex_handle = typename Tr::Vertex_handle;
  using FT = typename Tr::Geom_traits::FT;
  using Point = typename Tr::Point;

  std::string s;
  if(!(is >> s))
    return false;

  bool bin = (s == "binary");
  CGAL_assertion(bin == binary);
  if(bin) {
    if(!(is >> s))
      return false;
  }
  if(s != "CGAL" || !(is >> s) || s != "c3t3")
    return false;

  internal::set_binary(is, binary);

  c3t3.clear();
  Tr& tr = c3t3.triangulation();

  tr.tds().clear();
  tr.infinite_vertex() = tr.tds().create_vertex();

//  // todo : fix reading header
//  std::size_t n;
//  int d;
//  if(IO::is_ascii(is)) {
//    is >> d >> n;
//  } else {
//    read(is, d);
//    read(is, n);
//  }
//  if(!is)
//    return false;

  std::string str;
  do {
    is >> str;
  }
  while(str.compare("Vertices") != 0);
  std::size_t n;
  is >> n;

  // read points
  int d = 3;//todo : save/load dimension
  std::vector<Vertex_handle> vertices;
//  if(d > 3 || d < -2 || (n + 1) > vertices.max_size()) {
//    is.setstate(std::ios_base::failbit);
//    return false;
//  }
  tr.tds().set_dimension(d);
  vertices.resize(n + 1);
  vertices[0] = tr.infinite_vertex(); // the infinite vertex is numbered 0

  for(std::size_t i = 1; i <= n; i++) {
    vertices[i] = tr.tds().create_vertex();
    Point p;
    if(!(is >> p))
      return false;
    tr.set_point(vertices[i], p);
  }

  // read cells and fill the TDS
  std::vector<Cell_handle> cells;
  do {
    is >> str;
  } while(str.compare("Cells") != 0);

  tr.tds().read_cells(is, vertices, n, cells);

  // go to properties
  bool dimension_is_a_property = false;
  do
  {
    s.clear();
    if(!bin)
    {
      while(s.empty() && !is.eof()) // skip empty lines
        std::getline(is, s);
      if(is.eof())
      {
        c3t3.rescan_after_load_of_triangulation();
        return true;
      }
    }

    dimension_is_a_property = dimension_is_a_property
                            || (s.find("vertex:dimension") != std::string::npos);

    std::string ppty_value_type;
    std::size_t nb_simplices;
    std::size_t nb_components_per_simplex;
    std::string ppty_name;
    std::istringstream iss(s);
    if(!internal::read_property_header(iss,
          ppty_value_type, nb_simplices, nb_components_per_simplex, ppty_name))
    {
      std::cerr << "load_c3t3(): expected property header, got \"" << s << "\"\n";
      return false;
    }

    const std::string simplex_type_str = ppty_name.substr(0, ppty_name.find(':'));
    const auto simplex_type = internal::parse_simplex_type(simplex_type_str);
    const auto value_type = internal::parse_value_type(ppty_value_type);

    // read all simplices
    for(int i = 0; i < nb_simplices; ++i)
    {
      const auto [simplex_indices, ppty_values]
        = internal::read_simplex_and_values(simplex_type,
                                            value_type,
                                            nb_components_per_simplex,
                                            is);
      internal::apply_property(simplex_indices, ppty_name, ppty_values,
                               vertices, c3t3);
    }
  }
  while(!is.eof());

  c3t3.rescan_after_load_of_triangulation();

  // set dimension for each vertex
  if(!dimension_is_a_property)
    internal::set_vertex_dimensions(c3t3);

  // set Index for each vertex
  internal::set_vertex_indices(c3t3);

  return true;
}

//template <class C3t3>
//bool load_c3t3(std::istream& is
//             , C3t3& c3t3
//             , const bool binary = false)
//{
//  using P = std::variant<internal::Subdomain_index_property<C3t3>,
//                         internal::Surface_patch_index_property<C3t3>,
//                         internal::Curve_index_property<C3t3>,
//                         internal::Corner_index_property<C3t3>>;
//  std::vector<P> c3t3_properties{internal::Subdomain_index_property<C3t3>(),
//                                 internal::Surface_patch_index_property<C3t3>(),
//                                 internal::Curve_index_property<C3t3>(),
//                                 internal::Corner_index_property<C3t3>()};
//  return load_c3t3(is, c3t3, binary, c3t3_properties);
//}



///** * @ingroup PkgSMDS3IOFunctions
//  @todo implement and document */
//bool convert_old_to_new_c3t3(std::istream& is, std::ostream& os)
//{
//  // read old c3t3 from `is`
//  // write new c3t3 to `os` using `save_c3t3()`
//  return !is.fail() && !os.fail();
//}

} // end namespace IO

} // end namespace CGAL

#endif // CGAL_IO_FILE_C3T3_IO_H
