
namespace CGAL {
namespace Surface_mesh_topology {
/*!
\ingroup PkgSurfaceMeshTopologyClasses

The class `Polygonal_schema_min_items` defines a struct with a `std::string` as the information associated with darts, and no attribute is enabled.

\cgalModels{PolygonalSchemaItems}

\cgalHeading{Example}

The following example shows one implementation of the `Polygonal_schema_min_items` class.

\code{.cpp}

struct Polygonal_schema_min_items
{
  template <class PS>
  struct Dart_wrapper
  {
    struct Info_for_darts
    { std::string m_label; };
    typedef Info_for_darts Dart_info;
  };
};

\endcode

\sa `CGAL::Surface_mesh_topology::Polygonal_schema_with_combinatorial_map<Items,Alloc>`
\sa `CGAL::Surface_mesh_topology::Polygonal_schema_with_generalized_map<Items,Alloc>`

*/
struct Polygonal_schema_min_items {

}; /* end Linear_cell_complex_min_items */

} /* end namespace Surface_mesh_topology */
} /* end namespace CGAL */
