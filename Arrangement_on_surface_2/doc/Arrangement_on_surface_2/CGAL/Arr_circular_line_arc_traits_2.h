
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2TraitsClasses

This class is a traits class for \cgal arrangements, built on top of a
model of concept `CircularKernel`. It provides curves that can be
of both types
`CGAL::Line_arc_2<CircularKernel>` or
`CGAL::Circular_arc_2<CircularKernel>`.

It uses the <A HREF="https://www.boost.org/doc/html/variant.html">boost::variant</A>.

\cgalModels `ArrangementTraits_2`

*/
template< typename CircularKernel >
class Arr_circular_line_arc_traits_2 {
public:

}; /* end Arr_circular_line_arc_traits_2 */
} /* end namespace CGAL */
