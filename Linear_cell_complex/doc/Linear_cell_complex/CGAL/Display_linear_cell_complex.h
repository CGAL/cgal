namespace CGAL {
  
/*!
\ingroup PkgDisplayLinearCellComplex

\cgalModifBegin
Open a window and display `alcc`, a model of the `LinearCellComplex` concept. The running of the program will continue when the window will be closed. This function requires that CGAL_QT5 was enabled, and is only available if flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\cgalModifEnd 
*/
template<class LCC>
void display(const LCC& alcc);

} /* namespace CGAL */

