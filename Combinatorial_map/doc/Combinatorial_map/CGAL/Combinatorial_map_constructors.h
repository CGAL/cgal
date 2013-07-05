namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsConstructions

Creates a combinatorial hexahedron (six combinatorial quadrangles linked
together by \f$ \beta_2\f$), and adds it in `cm`.
Returns a handle on one dart of this combinatorial hexahedron.
\pre `CMap::dimension` \f$\geq\f$ 2.

\sa `CGAL::make_edge<CMap>`
\sa `CGAL::make_combinatorial_polygon<CMap>`
\sa `CGAL::make_combinatorial_tetrahedron<CMap>`

*/
template < class CMap >
typename CMap::Dart_handle make_combinatorial_hexahedron(CMap& cm);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsConstructions

Creates a combinatorial polygon of length `lg` (`lg` darts
linked by \f$ \beta_1\f$), and adds it in `cm`.
Returns a handle on one dart of this combinatorial polygon.
\pre `CMap::dimension`\f$ \geq\f$ 1 and `lg`\f$ >\f$ 0.


\sa `CGAL::make_edge<CMap>`
\sa `CGAL::make_combinatorial_tetrahedron<CMap>`
\sa `CGAL::make_combinatorial_hexahedron<CMap>`
*/
template < class CMap > typename CMap::Dart_handle
make_combinatorial_polygon(CMap& cm, unsigned int lg);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsConstructions

Creates a combinatorial tetrahedron (four combinatorial triangles linked
together by \f$ \beta_2\f$), and adds it in `cm`.
Returns a handle on one dart of this combinatorial tetrahedron.
\pre `CMap::dimension`\f$ \geq\f$ 2.

\sa `CGAL::make_edge<CMap>`
\sa `CGAL::make_combinatorial_polygon<CMap>`
\sa `CGAL::make_combinatorial_hexahedron<CMap>`
*/
template < class CMap >
typename CMap::Dart_handle make_combinatorial_tetrahedron(CMap& cm);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsConstructions

Creates an isolated edge (two darts linked by \f$ \beta_2\f$) and adds it in `cm`.
Returns a handle on one dart of this edge.
\pre `CMap::dimension`\f$ \geq\f$ 2.

\sa `CGAL::make_combinatorial_polygon<CMap>`
\sa `CGAL::make_combinatorial_tetrahedron<CMap>`
\sa `CGAL::make_combinatorial_hexahedron<CMap>`
*/
template < class CMap >
typename CMap::Dart_handle make_edge(CMap& cm);

} /* namespace CGAL */

