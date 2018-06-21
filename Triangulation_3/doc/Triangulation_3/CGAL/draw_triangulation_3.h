namespace CGAL {
  
/*!
\ingroup PkgDrawTriangulation3

Open a new window and draw `at3`, a model of the `TriangulationDataStructure_3` concept. The function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam T3 a model of the `TriangulationDataStructure_3` concept.
\param at3 the triangulation to draw.

*/
template<class T3>
void draw(const T3& at3);

} /* namespace CGAL */

