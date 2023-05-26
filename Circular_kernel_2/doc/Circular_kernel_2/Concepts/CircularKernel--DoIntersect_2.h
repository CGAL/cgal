
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

Testing whether two curves intersect.

\cgalRefines{Kernel::DoIntersect_2}

\sa \link do_intersect_grp `CGAL::do_intersect()` \endlink

*/

class CircularKernel::DoIntersect_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
determines if two geometric objects of type `Type1` and `Type2` intersect or not,
for all pairs `Type1` and `Type2`, where the types `Type1` and `Type2` can be any of the
following:
- `CircularKernel::Line_2`
- `CircularKernel::Line_arc_2`
- `CircularKernel::Circle_2`
- `CircularKernel::Circular_arc_2`
*/
bool operator()
(const Type1 & obj1, const Type2 & obj2);

/// @}

}; /* end CircularKernel::DoIntersect_2 */

