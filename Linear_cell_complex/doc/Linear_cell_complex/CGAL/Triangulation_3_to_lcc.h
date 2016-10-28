namespace CGAL{

/*!
\ingroup PkgLinearCellComplexConstructions

Imports `atr` (a `Triangulation_3`) into `lcc`. 
Objects are added in `lcc`, existing darts are not modified.
Returns a dart created during the import.
\pre \ref CombinatorialMap::dimension "LCC::dimension"\f$ \geq\f$ 3 and
  \ref Linear_cell_complex::ambient_dimension "LCC::ambient_dimension"==3.

\sa `CGAL::import_from_plane_graph<LCC>`
\sa `CGAL::import_from_polyhedron_3<LCC,Polyhedron>`
*/
template <class LCC,class Triangulation_>
typename LCC::Dart_handle import_from_triangulation_3(LCC& lcc,
const Triangulation_&atr);
}
