namespace CGAL
{

/*!
\ingroup PkgFeatureGraphRef

The enum `Surface_element_type` describe the different types that an element on the surface can have.

\sa `SurfaceElement_3`
*/
enum Surface_element_type
{
  POINT_ELEMENT = 0, //!< Descibes a 0D point element
  LINE_ELEMENT = 1, //!< Descibes a 1D line element
  SURFACE_ELEMENT = 2 //!< Descibes a 2D surface element
};

} /* namespace CGAL */