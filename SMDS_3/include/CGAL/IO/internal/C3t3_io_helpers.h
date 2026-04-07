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

#include <boost/container/small_vector.hpp>

#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <memory>
#include <type_traits>
#include <cstdint>


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
  // portable type names for properties
  template <typename T>
  std::string type_to_string()           { return "unknown"; }
  template <>
  std::string type_to_string<int8_t>()   { return "int8_t"; }
  template <>
  std::string type_to_string<uint8_t>()  { return "uint8_t"; }
  template <>
  std::string type_to_string<int16_t>()  { return "int16_t"; }
  template <>
  std::string type_to_string<uint16_t>() { return "uint16_t"; }
  template <>
  std::string type_to_string<int32_t>()  { return "int32_t"; }
  template <>
  std::string type_to_string<uint32_t>() { return "uint32_t"; }
  template <>
  std::string type_to_string<int64_t>()  { return "int64_t"; }
  template <>
  std::string type_to_string<uint64_t>() { return "uint64_t"; }

  // in c++23, <stdfloat> will define
  // std::float16_t
  // std::float32_t
  // std::float64_t
  // std::float128_t
  template <>
  std::string type_to_string<float>() { return "float32_t"; }
  template <>
  std::string type_to_string<double>() { return "float64_t"; }

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

  bool read_property_header(std::istream& is,
                            std::string& type,
                            std::size_t& nb_simplices,
                            std::size_t& nb_components_per_simplex,
                            std::string& ppty_name)
  {
    std::string hash;
    if(IO::is_ascii(is))
    {
      if(!(is >> hash >> type >> nb_simplices >> nb_components_per_simplex >> ppty_name)
        || hash != "#")
        return false;
    }
    else
    { //todo : deal with end of line in binary mode
      if(!(is >> hash >> type >> nb_simplices >> nb_components_per_simplex >> ppty_name))
        return false;
    }
    return true;
  }

  enum class NumberType { int8, int16, int32, int64, uint8, uint16, uint32, uint64, float32, float64, unknown = 100 };
  using NumberTypeVariant = std::variant<int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t, float, double>;

  enum class SimplexType { vertex = 1, edge, facet, cell, unknown = 100};
  std::size_t number_of_vertices(const SimplexType& st)
  {
    return static_cast<std::size_t>(st);
  }

  NumberType parse_value_type(const std::string& s)
  {
    if(s == "int8_t")         return NumberType::int8;
    else if(s == "int16_t")   return NumberType::int16;
    else if(s == "int32_t")   return NumberType::int32;
    else if(s == "int64_t")   return NumberType::int64;
    else if(s == "uint8_t")   return NumberType::uint8;
    else if(s == "uint16_t")  return NumberType::uint16;
    else if(s == "uint32_t")  return NumberType::uint32;
    else if(s == "uint64_t")  return NumberType::uint64;
    else if(s == "float32_t") return NumberType::float32;
    else if(s == "float64_t") return NumberType::float64;

    CGAL_assertion(false);
    return NumberType::unknown;
  }

  SimplexType parse_simplex_type(const std::string& s)
  {
    if(s == "vertex")     return SimplexType::vertex;
    else if(s == "edge")  return SimplexType::edge;
    else if(s == "facet") return SimplexType::facet;
    else if(s == "cell")  return SimplexType::cell;

    CGAL_assertion(false);
    return SimplexType::unknown;
  }

  template<typename V>
  NumberTypeVariant read_property_value(std::istream& is)
  {
    V value;
    is >> value;
    return value;
  }

  inline NumberTypeVariant read_by_number_type(const NumberType& nt,
                                               std::istream& is)
  {
    switch(nt)
    {
      case NumberType::int8:
        return read_property_value<int8_t>(is);
      case NumberType::int16:
        return read_property_value<int16_t>(is);
      case NumberType::int32:
        return read_property_value<int32_t>(is);
      case NumberType::int64:
        return read_property_value<int64_t>(is);
      case NumberType::uint8:
        return read_property_value<uint8_t>(is);
      case NumberType::uint16:
        return read_property_value<uint16_t>(is);
      case NumberType::uint32:
        return read_property_value<uint32_t>(is);
      case NumberType::uint64:
        return read_property_value<uint64_t>(is);
      case NumberType::float32:
        return read_property_value<float>(is);
      case NumberType::float64:
        return read_property_value<double>(is);
    };
    std::cerr << "Unknown number type " << std::endl;
    CGAL_assertion(false);
    return NumberTypeVariant();
  }

  inline auto
  read_simplex_and_values(const SimplexType& st,
                          const NumberType& nt,
                          const std::size_t& nb_components,
                          std::istream& is)
  {
    const std::size_t nb_indices = number_of_vertices(st);
    std::vector<std::size_t> indices(nb_indices);
    for(std::size_t i = 0; i < nb_indices; ++i)
    {
      std::size_t index;
      is >> index;
      indices[i] = index;
    }

    std::vector<NumberTypeVariant> values(nb_components);
    for(std::size_t i = 0; i < nb_components; ++i)
      values[i] = read_by_number_type(nt, is);

    return std::make_pair(indices, values);
  }

  template<typename C3t3>
  void get_vertex_and_apply_property(const std::size_t& i,
                                     const NumberTypeVariant& value,
                                     const std::string& ppty_name,
                                     const std::vector<typename C3t3::Vertex_handle>& vertices,
                                     C3t3 &c3t3)
  {
    using Corner_index = typename C3t3::Corner_index;
    typename C3t3::Vertex_handle v = vertices[i];

    if(ppty_name == "vertex:corner_index" && !c3t3.is_in_complex(v))
      c3t3.add_to_complex(v, std::get<Corner_index>(value));
  }
  template <typename C3t3>
  void get_edge_and_apply_property(const std::vector<std::size_t>& indices,
                                     const NumberTypeVariant& value,
                                     const std::string& ppty_name,
                                     const std::vector<typename C3t3::Vertex_handle>& vertices,
                                     C3t3& c3t3)
  {
    CGAL_assertion(indices.size() == 2);
    using Curve_index = typename C3t3::Curve_index;
    const typename C3t3::Vertex_handle v0 = vertices[indices[0]],
                                       v1 = vertices[indices[1]];

    if(ppty_name == "edge:curve_index" && !c3t3.is_in_complex(v0, v1))
      c3t3.add_to_complex(v0, v1, std::get<Curve_index>(value));
  }

  template <typename C3t3>
  void get_facet_and_apply_property(const std::vector<std::size_t>& indices,
                                     const NumberTypeVariant& value,
                                     const std::string& ppty_name,
                                     const std::vector<typename C3t3::Vertex_handle>& vertices,
                                     C3t3& c3t3)
  {
    CGAL_assertion(indices.size() == 3);
    using Surface_patch_index = typename C3t3::Surface_patch_index;
    using Cell_handle = typename C3t3::Cell_handle;

    const typename C3t3::Vertex_handle v0 = vertices[indices[0]],
                                       v1 = vertices[indices[1]],
                                       v2 = vertices[indices[2]];

    Cell_handle ch;
    int i, j, k;
    CGAL_assertion_code(bool ok =)
    c3t3.triangulation().is_facet(vertices[indices[0]],
                                  vertices[indices[1]],
                                  vertices[indices[2]],
                                  ch, i, j, k);
    const int findex = 6 - i - j - k;
    CGAL_assertion(ok);

    if(ppty_name == "facet:surface_patch_index" && !c3t3.is_in_complex(ch, findex))
      c3t3.add_to_complex(ch, findex, std::get<Surface_patch_index>(value));
  }
  template <typename C3t3>
  void get_cell_and_apply_property(const std::vector<std::size_t>& indices,
                                     const NumberTypeVariant& value,
                                     const std::string& ppty_name,
                                     const std::vector<typename C3t3::Vertex_handle>& vertices,
                                     C3t3& c3t3)
  {
    using Subdomain_index = typename C3t3::Subdomain_index;
    using Cell_handle = typename C3t3::Cell_handle;

    Cell_handle ch;
    CGAL_assertion_code(bool ok =)
    c3t3.triangulation().is_cell(vertices[indices[0]],
                                 vertices[indices[1]],
                                 vertices[indices[2]],
                                 vertices[indices[3]],
                                 ch);
    CGAL_assertion(ok);

    if(ppty_name == "cell:subdomain_index" && !c3t3.is_in_complex(ch))
      c3t3.add_to_complex(ch, std::get<Subdomain_index>(value));
  }

  template<typename C3t3>
  void apply_property(const std::vector<std::size_t>& simplex_indices,
                      const std::string& ppty_name,
                      const std::vector<NumberTypeVariant>& ppty_values,
                      const std::vector<typename C3t3::Vertex_handle>& vertices,
                      C3t3& c3t3)
  {
    const NumberTypeVariant ppty_value = ppty_values[0];
    //todo : deal with more than 1 ppty value
    switch(simplex_indices.size())
    {
    case 1:
      get_vertex_and_apply_property(simplex_indices[0], ppty_value, ppty_name, vertices, c3t3);
      break;
    case 2:
      get_edge_and_apply_property(simplex_indices, ppty_value, ppty_name, vertices, c3t3);
      break;
    case 3:
      get_facet_and_apply_property(simplex_indices, ppty_value, ppty_name, vertices, c3t3);
      break;
    case 4:
      get_cell_and_apply_property(simplex_indices, ppty_value, ppty_name, vertices, c3t3);
      break;
    default:
      std::cerr << "Error : simplex_indices has size "
                << simplex_indices.size() << std::endl;
    };
  }

  template <typename Simplex, typename C3t3>
  auto get_simplex_range(const C3t3& c3t3)
  {
    if constexpr(std::is_same_v<Simplex, typename C3t3::Cell_handle>)
      return c3t3.triangulation().all_cell_handles();
    else if constexpr(std::is_same_v<Simplex, typename C3t3::Facet>)
      return c3t3.triangulation().all_facets();
    else if constexpr(std::is_same_v<Simplex, typename C3t3::Edge>)
      return c3t3.triangulation().all_edges();
    else if constexpr(std::is_same_v<Simplex, typename C3t3::Vertex_handle>)
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
