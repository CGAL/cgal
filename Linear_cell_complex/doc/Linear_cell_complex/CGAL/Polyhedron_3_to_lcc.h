namespace CGAL{
/*!
\ingroup PkgLinearCellComplexConstructions

Imports `apoly` (a `Polyhedron_3`) into `lcc`, a model of the `LinearCellComplex` concept. Objects are added in `lcc`, existing darts are not modified. Returns a dart created during the import.
\pre \link GenericMap::dimension `LCC::dimension`\endlink \f$ \geq\f$ 2 and \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3.

\sa `CGAL::plane_graph_to_lcc<LCC>` (formerly `import_from_plane_graph`)
\sa `CGAL::triangulation_3_to_lcc<LCC,Triangulation>` (formerly `import_from_triangulation_3`)
*/
template<class LCC,class Polyhedron>
typename LCC::Dart_descriptor polyhedron_3_to_lcc(LCC& lcc,
const Polyhedron &apoly);
}
