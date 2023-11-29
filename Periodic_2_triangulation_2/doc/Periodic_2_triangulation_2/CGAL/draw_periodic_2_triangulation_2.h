namespace CGAL {

/*!
\ingroup PkgDrawPeriodic2Triangulation2

opens a new window and draws `ap2t2`, the `Periodic_2_Triangulation_2`. A call to this function is blocking, that is the program continues as soon as the user closes the window.
This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam P2T2 a model of the `Periodic_2TriangulationTraits_2` concept.
\param ap2t2 the 2D periodic trinagulation to draw.

*/
template<class P2T2>
void draw(const P2T2& ap2t2);

} /* namespace CGAL */
