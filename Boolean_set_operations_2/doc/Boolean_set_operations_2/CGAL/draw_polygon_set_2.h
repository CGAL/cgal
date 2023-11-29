namespace CGAL {

/*!
\ingroup PkgDrawPolygonSet2

opens a new window and draws `aps`, an instance of the `CGAL::Polygon_set_2` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam PS an instance of the `CGAL::Polygon_set_2` class.
\param aps the polygon set to draw.

*/
template<class PS>
void draw(const PS& aps);

} /* end namespace CGAL */
