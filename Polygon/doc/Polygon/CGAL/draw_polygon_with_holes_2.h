namespace CGAL {
  
/*!
\ingroup PkgDrawPolygonWithHoles2

opens a new window and draws `aph`, an instance of the `CGAL::Polygon_with_holes_2` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam PH an instance of the `CGAL::Polygon_with_holes_2` class.
\param aph the polygon with holes to draw.

*/
template<class PH>
void draw(const PH& aph);

} /* namespace CGAL */
