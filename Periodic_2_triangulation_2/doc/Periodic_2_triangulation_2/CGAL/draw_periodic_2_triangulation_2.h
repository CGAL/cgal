namespace CGAL {

/*!
\ingroup PkgDrawPeriodic2Triangulation2

opens a new window and draws `ap2t2`, the `Periodic_2_Triangulation_2`. A call to this function is blocking, that is the program continues as soon as the user closes the window.
 This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam P2T2 a model of the `Periodic_2TriangulationTraits_2` concept.
\param ap2t2 the 2D periodic trinagulation to draw.

*/
template<class P2T2>
void draw(const P2T2& ap2t2);

} /* namespace CGAL */
