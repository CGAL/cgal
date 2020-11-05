namespace CGAL {

/*!
\ingroup PkgNef3IOFunctions

This function creates a 3D Nef polyhedron from an OFF file which
is read from input stream `in`. The purpose of `OFF_to_nef_3`
is to create a Nef polyhedron from an OFF file that cannot be handled
by the `Nef_polyhedron_3` constructors. It handles double
coordinates while using a homogenous kernel, non-coplanar facets,
surfaces with boundaries, self-intersecting surfaces, and single
facets. Every closed volume gets marked. The function returns the
number of facets it could not handle.

\sa `CGAL::Nef_polyhedron_3<Traits>`

*/
template<class Nef_polyhedron_3>
std::size_t OFF_to_nef_3(std::istream& in, Nef_polyhedron_3& N);

} /* namespace CGAL */
