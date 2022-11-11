namespace CGAL {

/*!
\ingroup PkgTriangulationsVertexCellClasses

A `Triangulation_face` is a model of the concept `TriangulationDSFace`.

\tparam TriangulationDataStructure_ must be a model of the concept
`TriangulationDataStructure`.
Actually, `Triangulation_face` needs only that this concept defines the types
`Full_cell_handle`,
`Vertex_handle`, and
`Maximal_dimension`.

\cgalModels `TriangulationDSFace`

\sa `TriangulationDSFace`
\sa `TriangulationDataStructure`

*/
template< typename TriangulationDataStructure_ >
class Triangulation_face {
}; /* end Triangulation_face */

} /* end namespace CGAL */
