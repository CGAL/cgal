namespace CGAL {

/*!
\ingroup PkgDrawTriangulation2

opens a new window and draws `at2`, a model of the `TriangulationDataStructure_2` concept. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam T2 a model of the `TriangulationDataStructure_2` concept.
\param at2 the triangulation to draw.

*/
template<class T2>
void draw(const T2& at2);

} /* namespace CGAL */

