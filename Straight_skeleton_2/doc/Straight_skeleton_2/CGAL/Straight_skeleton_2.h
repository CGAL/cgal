namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Classes

\cgalModels `StraightSkeleton_2`

The class `Straight_skeleton_2` provides a model for the
`StraightSkeleton_2` concept which is the class
type used to represent a straight skeleton.

The only purpose of this class is to protect all the modifying
operations in a `HalfedgeDS`. Normal users should not modify a
straight skeleton. If an advanced user needs to get access to the
modifying operations, it must call the required methods through the
`Base` class.

\tparam Traits must be a model of `Kernel`

\sa `StraightSkeletonVertex_2`
\sa `StraightSkeletonHalfedge_2`
*/
template< typename Traits>
class Straight_skeleton_2
  : public CGAL::HalfedgeDS_default<Traits>
{
public:

  /// Access to the base (HDS) class
  typedef unspecified_type Base;

}; /* end Straight_skeleton_2 */

} /* end namespace CGAL */
