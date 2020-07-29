// Copyright (c) 2020 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Thien Hoang <thienvhoang99@gmail.com>
//
#ifndef CGAL_GENERIC_MAP_SELECTOR_H
#define CGAL_GENERIC_MAP_SELECTOR_H

#include <CGAL/license/Surface_mesh_topology.h>

// #include <CGAL/Polyhedron_3.h>
// #include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Polygonal_schema.h>

#include <tuple>

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

  struct Items_for_shortest_noncontractible_cycle
  {
    template <class Map>
    struct Dart_wrapper
    {
      using Vertex_attribute = CGAL::Cell_attribute<Map, int>;
      using Attributes       = std::tuple<Vertex_attribute>;
    };
  };

  template <class Mesh_, class Items_>
  struct SNC_for_generalized_map
  {
    using Mesh_original              = Mesh_;
    using Generic_map                = CGAL::Generalized_map<2, Items_>;
    using Dart_const_handle_original = typename Mesh_original::Dart_const_handle;
    using Copy_to_origin_map         = boost::unordered_map<typename Generic_map::Dart_handle,
                                                            Dart_const_handle_original>;
    using Origin_to_copy_map         = boost::unordered_map<Dart_const_handle_original,
                                                            typename Generic_map::Dart_handle>;

    static void copy(Generic_map& target, const Mesh_original& source,
                     Origin_to_copy_map& origin_to_copy,
                     Copy_to_origin_map& copy_to_origin,
                     typename Generic_map::size_type mark_perforated)
    {
      target.copy_from_const(source, &origin_to_copy, &copy_to_origin, true, mark_perforated);
    }
  };

  template <class Mesh_, class Items_>
  struct SNC_for_combinatorial_map
  {
    using Mesh_original              = Mesh_;
    using Generic_map                = CGAL::Combinatorial_map<2, Items_>;
    using Dart_const_handle_original = typename Mesh_original::Dart_const_handle;
    using Copy_to_origin_map         = boost::unordered_map<typename Generic_map::Dart_handle,
                                                            Dart_const_handle_original>;
    using Origin_to_copy_map         = boost::unordered_map<Dart_const_handle_original,
                                                            typename Generic_map::Dart_handle>;

    static void copy(Generic_map& target, const Mesh_original& source,
                     Origin_to_copy_map& origin_to_copy, Copy_to_origin_map& copy_to_origin,
                     typename Generic_map::size_type mark_perforated)
    {
      target.copy_from_const(source, &origin_to_copy, &copy_to_origin, true, mark_perforated);
    }
  };

  template <class Mesh_, class Items_>
  struct Generic_map_selector
  {
    using Mesh_original              = Mesh_;
    using Generic_map                = CGAL::Combinatorial_map<2, Items_>;
    using Dart_const_handle_original = typename boost::graph_traits<Mesh_original>::halfedge_descriptor;
    using Copy_to_origin_map         = boost::unordered_map<typename Generic_map::Dart_handle,
                                                            Dart_const_handle_original>;
    using Origin_to_copy_map         = boost::unordered_map<Dart_const_handle_original,
                                                            typename Generic_map::Dart_handle>;

    static void copy(Generic_map& target, const Mesh_original& source,
                     Origin_to_copy_map& origin_to_copy, Copy_to_origin_map& copy_to_origin,
                     typename Generic_map::size_type mark_perforated)
    {
      target.import_from_halfedge_graph(source, &origin_to_copy, &copy_to_origin, true, mark_perforated);
    }
  };

  template <unsigned int d, class Refs, class Items, class Alloc, class Storage, class Items2>
  struct Generic_map_selector<CGAL::Generalized_map_base<d, Refs, Items, Alloc, Storage>, Items2 >
    : SNC_for_generalized_map< CGAL::Generalized_map_base<d, Refs, Items, Alloc, Storage>, Items2 > {};

  template <unsigned int d, class Items, class Alloc, class Storage, class Items2>
  struct Generic_map_selector< CGAL::Generalized_map<d, Items, Alloc, Storage>, Items2 >
    : SNC_for_generalized_map< CGAL::Generalized_map<d, Items, Alloc, Storage>, Items2 > {};

  template <unsigned int d, unsigned int d2, class Traits, class Items,
            class Alloc, template<unsigned int,class,class,class,class> class Map,
            class Storage, class Items2>
  struct Generic_map_selector< CGAL::Linear_cell_complex_for_generalized_map
      <d, d2, Traits, Items, Alloc, Map, Storage>, Items2 >
      : SNC_for_generalized_map< CGAL::Linear_cell_complex_for_generalized_map
      <d, d2, Traits, Items, Alloc, Map, Storage>, Items2 > {};

  template <class Items_, class Alloc_, class Storage_, class Items2>
  struct Generic_map_selector< CGAL::Surface_mesh_topology::Polygonal_schema_with_generalized_map
                               <Items_, Alloc_, Storage_>, Items2 > : SNC_for_generalized_map
  <CGAL::Surface_mesh_topology::Polygonal_schema_with_generalized_map
      <Items_, Alloc_, Storage_>, Items2 > {};

  template <unsigned int d, class Refs, class Items, class Alloc, class Storage, class Items2>
  struct Generic_map_selector< CGAL::Combinatorial_map_base<d, Refs, Items, Alloc, Storage>, Items2 >
  : SNC_for_combinatorial_map< CGAL::Combinatorial_map_base<d, Refs, Items, Alloc, Storage>, Items2 > {};

  template <unsigned int d, class Items, class Alloc, class Storage, class Items2>
  struct Generic_map_selector< CGAL::Combinatorial_map<d, Items, Alloc, Storage>, Items2 >
  : SNC_for_combinatorial_map< CGAL::Combinatorial_map<d, Items, Alloc, Storage>, Items2 > {};

  template <unsigned int d, unsigned int d2, class Traits, class Items,
            class Alloc, template<unsigned int,class,class,class,class> class Map,
            class Storage, class Items2>
  struct Generic_map_selector< CGAL::Linear_cell_complex_for_combinatorial_map
      <d, d2, Traits, Items, Alloc, Map, Storage>, Items2 >
      : SNC_for_combinatorial_map< CGAL::Linear_cell_complex_for_combinatorial_map
      <d, d2, Traits, Items, Alloc, Map, Storage>, Items2 > {};

  template <class Items_, class Alloc_, class Storage_, class Items2>
  struct Generic_map_selector< CGAL::Surface_mesh_topology::Polygonal_schema_with_combinatorial_map
                               <Items_, Alloc_, Storage_>, Items2 > :
      SNC_for_combinatorial_map<CGAL::Surface_mesh_topology::
      Polygonal_schema_with_combinatorial_map<Items_, Alloc_, Storage_>, Items2 > {};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif
