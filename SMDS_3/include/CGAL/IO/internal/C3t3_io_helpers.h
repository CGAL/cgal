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

#ifndef CGAL_IO_INTERNAL_C3T3_IO_HELPERS_H
#define CGAL_IO_INTERNAL_C3T3_IO_HELPERS_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/IO/C3t3_io.h>

#include <CGAL/SMDS_3/io_signature.h>

#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <memory>
#include <type_traits>


namespace CGAL {

namespace IO {

template <typename C3t3,
          typename Simplex,
          typename Value_type>
class Property
{
public:
  virtual ~Property() = default;
  virtual std::string name() const = 0;
  // nb elements per simplex (3 for coordinates, 1 for index)
  virtual int number_of_components() const = 0;
  virtual std::string type() const = 0;
  // nb simplices (cells, facets, edges, vertices) for which the property is defined
  virtual int number_of_simplices(const C3t3&) const = 0;
  virtual Value_type operator()(const Simplex&, const C3t3&) const = 0;
  virtual bool has_value(const Simplex&, const C3t3&) const = 0;
};


namespace internal {

//  template <class Base, class... Derived>
//  std::vector<std::unique_ptr<Base>>
//  make_properties_vector()
//  {
//    std::vector<std::unique_ptr<Base>> v;
//    v.reserve(sizeof...(Derived));
//    (v.emplace_back(std::make_unique<Derived>()), ...);
//    return v;
//  }
//
  bool is_cell_property(const std::string& ppty_name)
  {
    return ppty_name.compare(0, 5, "cell:") == 0;
  }
  bool is_facet_property(const std::string& ppty_name)
  {
    return ppty_name.compare(0, 6, "facet:") == 0;
  }
  bool is_edge_property(const std::string& ppty_name)
  {
    return ppty_name.compare(0, 5, "edge:") == 0;
  }
  bool is_vertex_property(const std::string& ppty_name)
  {
    return ppty_name.compare(0, 7, "vertex:") == 0;
  }

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

  //////////////////////
  // WRITE PROPERTIES //
  //////////////////////
  template <typename T>
  std::string type_to_string(){ return "unknown"; }
  template <>
  std::string type_to_string<double>() { return "double";}
  template <>
  std::string type_to_string<float>() { return "float"; }
  template <>
  std::string type_to_string<int>() { return "int"; }


  template <typename C3t3>
  class Subdomain_index_property
    : public Property<C3t3, typename C3t3::Cell_handle, typename C3t3::Subdomain_index>
  {
  public:
    using key_type   = typename C3t3::Cell_handle;
    using value_type = typename C3t3::Subdomain_index;
    std::string name()         const { return "cell:subdomain_index"; }
    int number_of_components() const { return 1; }
    std::string type()         const { return type_to_string<value_type>(); }
    int number_of_simplices(const C3t3& c3t3) const
    {
      return c3t3.number_of_cells_in_complex();
    }
    bool has_value(const key_type& c, const C3t3& c3t3) const
    {
      return c3t3.is_in_complex(c);
    }
    value_type operator()(const key_type& c, const C3t3& c3t3) const
    {
      return c3t3.subdomain_index(c);
    }
  };

  template <typename C3t3>
  class Surface_patch_index_property
    : public Property<C3t3, typename C3t3::Facet, typename C3t3::Surface_patch_index>
  {
  public:
    using key_type = typename C3t3::Facet;
    using value_type = typename C3t3::Surface_patch_index;
    std::string name()         const { return "facet:surface_patch_index"; }
    int number_of_components() const { return 1; }
    std::string type()         const { return type_to_string<value_type>(); }
    int number_of_simplices(const C3t3& c3t3) const
    {
      return c3t3.number_of_facets_in_complex();
    }
    bool has_value(const typename C3t3::Facet& f, const C3t3& c3t3) const
    {
        return c3t3.is_in_complex(f);
    }
    value_type operator()(const typename C3t3::Facet& f, const C3t3& c3t3) const
    {
      return c3t3.surface_patch_index(f);
    }
  };

  template <typename C3t3>
  class Curve_index_property
    : public Property<C3t3, typename C3t3::Edge, typename C3t3::Curve_index>
  {
  public:
    using key_type = typename C3t3::Edge;
    using value_type = typename C3t3::Curve_index;
    std::string name()         const { return "edge:curve_index"; }
    int number_of_components() const { return 1; }
    std::string type()         const { return type_to_string<value_type>(); }
    int number_of_simplices(const C3t3& c3t3) const
    {
      return c3t3.number_of_edges_in_complex();
    }
    bool has_value(const typename C3t3::Edge& e, const C3t3& c3t3) const
    {
      return c3t3.is_in_complex(e);
    }
    value_type operator()(const typename C3t3::Edge& e, const C3t3& c3t3) const
    {
      return c3t3.curve_index(e);
    }
  };

  template <typename C3t3>
  class Corner_index_property
    : public Property<C3t3, typename C3t3::Vertex_handle, typename C3t3::Corner_index>
  {
  public:
    using key_type = typename C3t3::Vertex_handle;
    using value_type = typename C3t3::Corner_index;
    std::string name()         const { return "vertex:corner_index"; }
    int number_of_components() const { return 1; }
    std::string type()         const { return type_to_string<value_type>(); }
    int number_of_simplices(const C3t3& c3t3) const
    {
      return c3t3.number_of_vertices_in_complex();
    }
    bool has_value(const typename C3t3::Vertex_handle& v, const C3t3& c3t3) const
    {
      return c3t3.is_in_complex(v);
    }
    value_type operator()(const typename C3t3::Vertex_handle& v, const C3t3& c3t3) const
    {
      return c3t3.corner_index(v);
    }
  };



  //////////////////////
  // READ PROPERTIES ///
  //////////////////////
  template <typename Tr>
  void read_cell_property(std::istream& is,
                          Tr& tr,
                          const std::string& ppty_name)
  {
    std::cout << "read_cell_property(): property name = " << ppty_name << "\n";
    //todo
  }

  template <typename Tr>
  void read_facet_property(std::istream& is,
                           Tr& tr,
                           const std::string& ppty_name)
  {
    std::cout << "read_facet_property(): property name = " << ppty_name << "\n";
    // todo
  }

  template <typename Tr>
  void read_edge_property(std::istream& is,
                          Tr& tr,
                          const std::string& ppty_name)
  {
    std::cout << "read_edge_property(): property name = " << ppty_name << "\n";
    // todo
  }

  template <typename Tr>
  void read_vertex_property(std::istream& is,
                            Tr& tr,
                            const std::string& ppty_name)
  {
    std::cout << "read_vertex_property(): property name = " << ppty_name << "\n";
    // todo
  }

  template <typename Simplex, typename C3t3>
  auto get_simplex_range(const C3t3& c3t3)
  {
    if constexpr(std::is_same<Simplex, typename C3t3::Cell_handle>::value)
      return c3t3.triangulation().all_cell_handles();
    else if constexpr(std::is_same<Simplex, typename C3t3::Facet>::value)
      return c3t3.triangulation().all_facets();
    else if constexpr(std::is_same<Simplex, typename C3t3::Edge>::value)
      return c3t3.triangulation().all_edges();
    else if constexpr(std::is_same<Simplex, typename C3t3::Vertex_handle>::value)
      return c3t3.triangulation().all_vertex_handles();
    else
      CGAL_assertion(false);
  };

  template <typename Simplex, typename C3t3, typename V2I>
  void write_property(std::ostream& os,
                      const Simplex& s,
                      const C3t3& c3t3,
                      const V2I& vertices)
  {
    const auto& v = c3t3.triangulation().vertices(s);
    for(const auto& v_i : v)
      os << vertices[v_i] << " ";
  }

} // end namespace internal
} // end namespace IO
} // end namespace CGAL

#endif // CGAL_IO_INTERNAL_C3T3_IO_HELPERS_H
