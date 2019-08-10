#ifndef CGAL_GENERIC_MAP_SELECTOR_H
#define CGAL_GENERIC_MAP_SELECTOR_H

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

  struct Attributes {
    template <class GenericMap>
    struct Dart_wrapper {
      using Vertex_attribute = CGAL::Cell_attribute<GenericMap, int>;
      using Edge_attribute   = CGAL::Cell_attribute<GenericMap, int>;
      using Face_attribute   = void;
      using Attributes       = CGAL::cpp11::tuple<Vertex_attribute, Edge_attribute, Face_attribute>;
    };
  };

  template <class Mesh_>
  struct SNC_for_generalized_map {
    using Mesh_original              = Mesh_;
    using Generic_map                = CGAL::Generalized_map<2, Attributes>;
    using Dart_handle_original       = typename Mesh_original::Dart_handle;
    using Dart_const_handle_original = typename Mesh_original::Dart_const_handle;
    using Copy_to_origin_map         = boost::unordered_map<typename Generic_map::Dart_handle, 
                                                            Dart_handle_original>;
    using Origin_to_copy_map         = boost::unordered_map<Dart_handle_original, 
                                                            typename Generic_map::Dart_handle>;

    static void copy(Generic_map& target, Mesh_original& source,
                     Origin_to_copy_map& origin_to_copy, Copy_to_origin_map& copy_to_origin)
    {
      target.copy(source, &origin_to_copy, &copy_to_origin);
    }
  };

  template <class Mesh_>
  struct SNC_for_combinatorial_map {
    using Mesh_original              = Mesh_;
    using Generic_map                = CGAL::Combinatorial_map<2, Attributes>;
    using Dart_handle_original       = typename Mesh_original::Dart_handle;
    using Dart_const_handle_original = typename Mesh_original::Dart_const_handle;
    using Copy_to_origin_map         = boost::unordered_map<typename Generic_map::Dart_handle, 
                                                            Dart_handle_original>;
    using Origin_to_copy_map         = boost::unordered_map<Dart_handle_original, 
                                                            typename Generic_map::Dart_handle>;

    static void copy(Generic_map& target, Mesh_original& source,
                     Origin_to_copy_map& origin_to_copy, Copy_to_origin_map& copy_to_origin)
    {
      target.copy(source, &origin_to_copy, &copy_to_origin);
    }
  };

  template <class Mesh_>
  struct Generic_map_selector {
    using Mesh_original              = Mesh_;
    using Generic_map                = CGAL::Combinatorial_map<2, Attributes>;
    using Dart_handle_original       = typename boost::graph_traits<Mesh_original>::halfedge_descriptor;
    using Dart_const_handle_original = typename boost::graph_traits<Mesh_original>::halfedge_descriptor;
    using Copy_to_origin_map         = boost::unordered_map<typename Generic_map::Dart_handle, 
                                                            Dart_handle_original>;
    using Origin_to_copy_map         = boost::unordered_map<Dart_handle_original, 
                                                            typename Generic_map::Dart_handle>;

    static void copy(Generic_map& target, Mesh_original& source,
                     Origin_to_copy_map& origin_to_copy, Copy_to_origin_map& copy_to_origin)
    {
      target.import_from_halfedge_graph(source, &origin_to_copy, &copy_to_origin);
    }
  };

  template <unsigned int d, class Items, class Alloc, class Storage>
  struct Generic_map_selector< CGAL::Generalized_map<d, Items, Alloc, Storage> > 
    : SNC_for_generalized_map< CGAL::Generalized_map<d, Items, Alloc, Storage> > {};

  template <unsigned int d, class Refs, class Items, class Alloc, class Storage>
  struct Generic_map_selector< CGAL::Generalized_map_base
                                <d, Refs, Items, Alloc, Storage> > 
    : SNC_for_generalized_map< CGAL::Generalized_map_base
                                <d, Refs, Items, Alloc, Storage> > {};

  template <unsigned int d, unsigned int d2, class Traits, class Items,
            class Alloc, template<unsigned int,class,class,class,class> class Map, class Storage>
  struct Generic_map_selector< CGAL::Linear_cell_complex_for_generalized_map
                              <d, d2, Traits, Items, Alloc, Map, Storage> >
    : SNC_for_generalized_map< CGAL::Linear_cell_complex_for_generalized_map
                              <d, d2, Traits, Items, Alloc, Map, Storage> > {};

  template <unsigned int d, class Items, class Alloc, class Storage>
  struct Generic_map_selector< CGAL::Combinatorial_map<d, Items, Alloc, Storage> >
  : SNC_for_combinatorial_map< CGAL::Combinatorial_map<d, Items, Alloc, Storage> > {};

  template <unsigned int d, class Refs, class Items, class Alloc, class Storage>
  struct Generic_map_selector< CGAL::Combinatorial_map_base
                                <d, Refs, Items, Alloc, Storage> >
  : SNC_for_combinatorial_map< CGAL::Combinatorial_map_base
                                <d, Refs, Items, Alloc, Storage> > {};

  template <unsigned int d, unsigned int d2, class Traits, class Items,
            class Alloc, template<unsigned int,class,class,class,class> class Map, class Storage>
  struct Generic_map_selector< CGAL::Linear_cell_complex_for_combinatorial_map
                              <d, d2, Traits, Items, Alloc, Map, Storage> >
  : SNC_for_combinatorial_map< CGAL::Linear_cell_complex_for_combinatorial_map
                              <d, d2, Traits, Items, Alloc, Map, Storage> > {};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif