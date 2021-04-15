
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2Tags

The categories `Left_side_category`, `Right_side_category`,
`Bottom_side_category`, and `Top_side_category`,
nested in any model of the `ArrangementBasicTraits_2`, must be
convertible to `Arr_oblivious_side_tag`. `Arr_oblivious_side_tag` is an empty construct used
for dispatching functions based on type of curves that induce the
arrangement.

\sa `Arr_open_side_tag`
\sa `ArrangementBasicTraits_2`

*/

struct Arr_oblivious_side_tag {

}; /* end Arr_oblivious_side_tag */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2Tags

All the four types `Left_side_category`,
`Right_side_category`, `Bottom_side_category`,
and `Top_side_category` nested in any model of the
concept `ArrangementOpenBoundaryTraits` must be convertible
to `Arr_open_side_tag`, which derives from `Arr_oblivious_side_tag`. It
implies that some curves are expected to approach the left, right,
bottom, or top sides of the open boundary of the parameter
space. `Arr_open_side_tag` is an empty construct used for dispatching
functions based on type of curves that induce the arrangement.

\sa `Arr_oblivious_side_tag`
\sa `ArrangementOpenBoundaryTraits_2`

*/

struct Arr_open_side_tag {

}; /* end Arr_open_side_tag */
} /* end namespace CGAL */
