
namespace CGAL {

/*!
\ingroup PkgTDS3Classes

The class `Triangulation_ds_vertex_base_3` can be used as the base vertex
for a 3D-triangulation data structure, it is a model of the concept
`TriangulationDSVertexBase_3`.

Note that if the triangulation data structure is used as a parameter of a
geometric triangulation (Section \ref TDS3secdesign and
Chapter \ref chapterTriangulation3 "3D Triangulations"), then the vertex base class has to
fulfill additional geometric requirements, i.e.\ it has to be a model of the
concept `TriangulationVertexBase_3`.

This base class can be used directly or can serve as a base to derive
other base classes with some additional attributes (a color for
example) tuned for a specific application.

\cgalModels{TriangulationDSVertexBase_3}

\tparam TDS should not be specified (see Section \ref tds3cyclic and examples)

\sa `CGAL::Triangulation_vertex_base_3`
\sa `CGAL::Triangulation_ds_cell_base_3`

*/
template< typename TDS = void >
class Triangulation_ds_vertex_base_3 {

}; /* end Triangulation_ds_vertex_base_3 */
} /* end namespace CGAL */
