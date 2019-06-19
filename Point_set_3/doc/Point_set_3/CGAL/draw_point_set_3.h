namespace CGAL {
  
/*!
\ingroup PkgDrawPointSet3D

opens a new window and draws `aps`, an instance of the `CGAL::Point_set_3` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam PS an instance of the `CGAL::Point_set_3` class.
\param aps the point set to draw.

*/
template<class PS>
void draw(const PS& aps);

} /* namespace CGAL */
