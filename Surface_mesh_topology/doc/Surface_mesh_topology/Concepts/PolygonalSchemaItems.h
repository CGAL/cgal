/*!
  \ingroup PkgSurfaceMeshTopologyConcepts
  \cgalConcept

  The concept `PolygonalSchemaItems` allows to customize a PolygonalSchema  by choosing the information associated with darts, and by enabling and disabling some attributes.

  \cgalRefines GenericMapItems

  \cgalHasModel \link CGAL::Surface_mesh_topology::Polygonal_schema_min_items `CGAL::Surface_mesh_topology::Polygonal_schema_min_items`\endlink
  */

class PolygonalSchemaItems
{
    /*!
  - `%Dart_wrapper<Map>::%Dart_info`, should be a class having a public data member char* m_label.
    */
};
