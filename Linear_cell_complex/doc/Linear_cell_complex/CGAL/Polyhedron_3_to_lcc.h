namespace CGAL{
/*!
\ingroup PkgLinearCellComplexConstructions

Imports `apoly`into `lcc`. Objects are added in `lcc`, existing darts are not modified. Returns a dart created during the import.
\pre \link GenericMap::dimension `LCC::dimension`\endlink \f$ \geq\f$ 2 and \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3.

\tparam LCC a model of `LinearCellComplex`
\tparam PolygonMesh a model of `FaceGraph`

\sa `CGAL::read_plane_graph_in_lcc<LCC>`
\sa `CGAL::triangulation_3_to_lcc<LCC,Triangulation>`
*/
template<class LCC,class PolygonMesh>
typename LCC::Dart_descriptor polyhedron_3_to_lcc(LCC& lcc,
                                                  const PolygonMesh &apoly);
}
