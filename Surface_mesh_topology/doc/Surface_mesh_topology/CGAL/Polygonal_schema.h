
namespace CGAL {
namespace Surface_mesh_topology {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses

    The class `Polygonal_schema_with_combinatorial_map` is a model of `PolygonalSchema` using `Combinatorial_map` as underlying model.

    \tparam Items a model of the `PolygonalSchemaItems` concept. Equal to `Polygonal_schema_min_items` by default.

    \tparam Alloc has to match the standard allocator requirements. The `rebind` mechanism  `Alloc` will be used to create appropriate allocators internally with value type `Dart`. Equal to `CGAL_ALLOCATOR(int)` from `<CGAL/memory.h>` by default.
  */
  template<typename Items, typename Alloc>
  class Polygonal_schema_with_combinatorial_map: public CGAL::Combinatorial_map<2, Items, Alloc>
  {
  };

  /*!
    \ingroup PkgSurfaceMeshTopologyClasses

    The class `Polygonal_schema_with_generalized_map` is a model of `PolygonalSchema` using `Generalized_map` as underlying model.

    \tparam Items a model of the `PolygonalSchemaItems` concept. Equal to `Polygonal_schema_min_items` by default.

    \tparam Alloc has to match the standard allocator requirements. The `rebind` mechanism  `Alloc` will be used to create appropriate allocators internally with value type `Dart`. Equal to `CGAL_ALLOCATOR(int)` from `<CGAL/memory.h>` by default.
  */
  template<typename Items, typename Alloc>
  class Polygonal_schema_with_generalized_map: public CGAL::Generalized_map<2, Items, Alloc>
  {
  };
}
}
