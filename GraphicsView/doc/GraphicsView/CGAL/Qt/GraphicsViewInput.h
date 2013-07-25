namespace CGAL {
namespace Qt {
/*!
\ingroup PkgGraphicsViewInputClasses

An object of type `GraphicsViewInput` can emit a signal with `CGAL::Object` as argument. 

*/

class GraphicsViewInput {
public:

/// \name Signals 
/// @{

/*!
A signal that emits a `CGAL::Object`. 
*/ 
void generate(Object); 

/// @}

}; /* end GraphicsViewInput */
} /* end namespace Qt */
} /* end namespace CGAL */
