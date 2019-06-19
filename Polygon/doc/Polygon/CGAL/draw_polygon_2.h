namespace CGAL {
  
/*!
\ingroup PkgDrawPolygon2

opens a new window and draws `ap`, an instance of the `CGAL::Polygon_2` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam P an instance of the `CGAL::Polygon_2` class.
\param ap the polygon to draw.

*/
template<class P>
void draw(const P& ap);

} /* namespace CGAL */
