namespace CGAL {

/*!
\ingroup PkgLinearCellComplexOperations

Returns the scaled Jacobian of the hexahedron containing `d`, a quality measure
of hexahedra between -1 and 1. At each corner of the hexahedron, the scaled
Jacobian is the determinant of the three unit vectors along its incident edges.
The value returned is the minimum over the eight corners: 1 for a cube, 0 or
less for a degenerate or inverted hexahedron. A corner with an edge of length 0
has a value of 0.

The sign depends on the orientation of the darts. For example, the hexahedron
created by `make_hexahedron()` with the points (0,0,0), (1,0,0), (1,1,0), (0,1,0),
(0,1,1), (0,0,1), (1,0,1), (1,1,1) has a scaled Jacobian of 1, and the same
hexahedron with mirrored points has a scaled Jacobian of -1.

\pre \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3,
`d` \f$ \in \f$ \link GenericMap::darts `lcc.darts()`\endlink, and
`lcc.is_volume_combinatorial_hexahedron(d)`.
*/
template <class LCC>
double scaled_jacobian_of_hexahedron(const LCC& lcc,
                                     typename LCC::Dart_const_descriptor d);

} /* namespace CGAL */
