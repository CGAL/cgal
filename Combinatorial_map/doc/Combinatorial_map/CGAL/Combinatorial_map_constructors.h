namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsConstructions

\deprecated Creates a combinatorial hexahedron. Deprecated. Use `cm.make_combinatorial_hexahedron()` instead.
*/
template < class CMap >
typename CMap::Dart_handle make_combinatorial_hexahedron(CMap& cm);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsConstructions

\deprecated Creates a combinatorial polygon of length `lg`.  Deprecated. Use `cm.make_combinatorial_polygon()` instead.
*/
template < class CMap > typename CMap::Dart_handle
make_combinatorial_polygon(CMap& cm, unsigned int lg);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsConstructions

\deprecated Creates a combinatorial tetrahedron.  Deprecated. Use `cm.make_combinatorial_tetrahedron()` instead.
*/
template < class CMap >
typename CMap::Dart_handle make_combinatorial_tetrahedron(CMap& cm);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsConstructions

\deprecated Creates an isolated edge.  Deprecated. Use `cm.make_edge()` instead.
*/
template < class CMap >
typename CMap::Dart_handle make_edge(CMap& cm);

} /* namespace CGAL */

