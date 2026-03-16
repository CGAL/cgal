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
    os << v->point() << "\n";
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
template <class C3t3,
          typename Properties>
bool load_c3t3(std::istream& is
             , C3t3& c3t3
             , const bool binary
             , const Properties& properties)
{
  using Tr = typename C3t3::Triangulation;
  using Cell_handle = typename Tr::Cell_handle;
  using Facet = typename Tr::Facet;
  using Edge = typename Tr::Edge;
  using Vertex_handle = typename Tr::Vertex_handle;

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

    std::size_t n;
  int d;
  if(IO::is_ascii(is)) {
    is >> d >> n;
  } else {
    read(is, d);
    read(is, n);
  }
  if(!is)
    return false;

  // read points
  std::vector<Vertex_handle> vertices;
  if(d > 3 || d < -2 || (n + 1) > vertices.max_size()) {
    is.setstate(std::ios_base::failbit);
    return false;
  }
  tr.tds().set_dimension(d);
  vertices.resize(n + 1);
  vertices[0] = tr.infinite_vertex(); // the infinite vertex is numbered 0

  for(std::size_t i = 1; i <= n; i++) {
    vertices[i] = tr.tds().create_vertex();
    if(!(is >> *vertices[i]))
      return false;
  }

  // read cells and fill the TDS
  std::vector<Cell_handle> cells;
  std::size_t m;
  tr.tds().read_cells(is, vertices, m, cells);

  // go to properties
  s.clear();
  do
  {
    if(!bin)
    {
      while(s.empty() && !is.eof()) // skip empty lines
        std::getline(is, s);
    }

    std::string ppty_name;
    int nb_elements;
    int nb_values_per_element;
    std::string ppty_type;
    if(IO::is_ascii(is))
    {
      std::istringstream iss(s);
      std::string hash;
      if(!(iss >> hash >> ppty_name >> nb_elements >> nb_values_per_element >> ppty_type)
        || hash != "#")
        return false;
    }
    else
    {
      if(!(is >> ppty_name >> nb_elements >> nb_values_per_element >> ppty_type))
        return false;
    }

  //  // todo : read the type of the property,
  //  //        and use it to read the value of the property
    if(internal::is_cell_property(ppty_name))
    {
      ;//internal::read_cell_property(is, tr, ppty.substr(6));
    }
    else if(internal::is_facet_property(ppty_name))
    {
      ;//internal::read_facet_property(is, tr, ppty.substr(6));
    }
    else if(internal::is_edge_property(ppty_name))
    {
      ;//internal::read_edge_property(is, tr, ppty.substr(5));
    }
    else if(internal::is_vertex_property(ppty_name))
    {
      ;//internal::read_vertex_property(is, tr, ppty.substr(7));
    }
    else
    {
      std::cerr << "load_c3t3(): unknown property name \"" << ppty_name << "\"\n";
      return false;
    }
  }
  while(!is.eof());

  c3t3.rescan_after_load_of_triangulation();

  return !is.fail();
}

template <class C3t3>
bool load_c3t3(std::istream& is
             , C3t3& c3t3
             , const bool binary = false)
{
  using P = std::variant<internal::Subdomain_index_property<C3t3>,
                         internal::Surface_patch_index_property<C3t3>,
                         internal::Curve_index_property<C3t3>,
                         internal::Corner_index_property<C3t3>>;
  std::vector<P> c3t3_properties{internal::Subdomain_index_property<C3t3>(),
                                 internal::Surface_patch_index_property<C3t3>(),
                                 internal::Curve_index_property<C3t3>(),
                                 internal::Corner_index_property<C3t3>()};
  return load_c3t3(is, c3t3, binary, c3t3_properties);
}



/** * @ingroup PkgSMDS3IOFunctions
  @todo implement and document */
bool convert_old_to_new_c3t3(std::istream& is, std::ostream& os)
{
  // read old c3t3 from `is`
  // write new c3t3 to `os` using `save_c3t3()`
  return !is.fail() && !os.fail();
}

} // end namespace IO

} // end namespace CGAL

#endif // CGAL_IO_FILE_C3T3_IO_H
