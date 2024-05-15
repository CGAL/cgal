namespace CGAL {

/*!
\ingroup PkgDrawLinearCellComplex

opens a new window and draws `alcc`, a model of the `LinearCellComplex` concept. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam LCC a model of the `LinearCellComplex` concept.
\param alcc the linear cell complex to draw.

*/
template<class LCC>
void draw(const LCC& alcc);

} /* namespace CGAL */

