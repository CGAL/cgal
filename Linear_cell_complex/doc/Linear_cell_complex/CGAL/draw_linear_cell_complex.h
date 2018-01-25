namespace CGAL {
  
/*!
\ingroup PkgDrawLinearCellComplex

\cgalModifBegin
Open a new window and draw `alcc`, a model of the `LinearCellComplex` concept. The function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam LCC a model of the `LinearCellComplex` concept.
\param alcc The linear cell complex to draw.

\cgalModifEnd 
*/
template<class LCC>
void draw(const LCC& alcc);

} /* namespace CGAL */

