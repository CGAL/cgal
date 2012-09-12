
namespace CGAL {

/*!
\ingroup PkgArrangement2

The class `Arr_circle_segment_traits_2` is a model of the `ArrangementTraits_2` concept 
and can be used to construct and maintain arrangements of circular arcs 
and line segments. 

The traits class must be instantiated with a geometric kernel, such that 
the supporting circles of the circular arcs are of type `Kernel::Circle_2` 
and the supporting lines of the line segments are of type `Kernel::Line_2`. 
Thus, the coordinates of the center of supporting circles, and its squared 
radius are of type `Kernel::FT`, which should be an exact rational 
number-type; similarly, the coefficients of each supporting line 
\f$ ax + by + c = 0\f$ are also of type `Kernel::FT`. Note however that 
the intersection point between two such arcs do not have rational coordinates 
in general. For this reason, we do not require the endpoints of the input arcs 
and segments to have rational coordinates. 

The nested `Point_2` type defined by the traits class is therefore 
<I>different</I> than the `Kernel::Point_2` type. Its coordinates are 
of type `CoordNT`, which an instantiation of 
`CGAL::Sqrt_extension<NT,ROOT>` where `NT = ROOT = Kernel::FT`. 
Moreover, the third and fourth (hidden) template parameters of 
`CGAL::Sqrt_extension<NT,ROOT>` are set to `CGAL::Tag_true`, which 
enables efficient comparison among different extensions. 

For more details see the documentation of \ref ::CGAL::Sqrt_extension<NT,ROOT> 

While `Arr_circle_segment_traits_2` models the concept 
`ArrangementDirectionalXMonotoneTraits_2`, the implementation of 
the `Arr_mergeable_2` operation does not enforce the input curves 
to have the same direction as a precondition. Moreover, `Arr_circle_segment_traits_2` 
supports the merging of curves of opposite directions. 

\models ::ArrangementTraits_2 
\models ::ArrangementDirectionalXMonotoneTraits_2 

CONVERROR 3 nested classes missing 

*/
template< typename Kernel >
class Arr_circle_segment_traits_2 {

}; /* end Arr_circle_segment_traits_2 */
} /* end namespace CGAL */
