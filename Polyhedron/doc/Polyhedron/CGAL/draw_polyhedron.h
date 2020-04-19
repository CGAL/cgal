namespace CGAL {

/*!
\ingroup PkgDrawPolyhedron

opens a new window and draws `apoly`, an instance of the `CGAL::Polyhedron_3` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam POLY an instance of the `CGAL::Polyhedron_3` class.
\param apoly the polyhedron to draw.

*/
template<class POLY>
void draw(const POLY& apoly);

} /* namespace CGAL */

