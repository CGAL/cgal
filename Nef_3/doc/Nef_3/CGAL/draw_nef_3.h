namespace CGAL {

/*!
\ingroup PkgDrawNef3

Open a new window and draws `anef3`, the `Nef_polyhedron_3`. A call to this function is blocking, that is the program continues as soon as the user closes the window.
 This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam Nef3 a model of the `Nef_polyhedron_3` concept.
\param anef3 the nef polyhedron to draw.

*/
template<class Nef3>
void draw(const Nef3& anef3);

} /* namespace CGAL */