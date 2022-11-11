namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Tags
 *
 * This type tag is used to indicate that the condition of a boundary side of
 * the parameter space is irrelevant. More precisely, when curves inserted into
 * the arrangement are never expected to reach a boundary side, either left,
 * right, bottom, or top, then the corresponding categories, namely,
 * `Left_side_category`, `Right_side_category`, `Bottom_side_category`, and
 * `Top_side_category`, nested in every geometry traits class, must be
 * convertible to the type `Arr_oblivious_side_tag`.
 *
 * `Arr_oblivious_side_tag` is an empty construct used for dispatching functions
 * based on type of curves that induce the arrangement.
 *
 * \sa `Arr_open_side_tag`
 * \sa `Arr_closed_side_tag`
 * \sa `Arr_contracted_side_tag`
 * \sa `Arr_identified_side_tag`
 * \sa `ArrangementBasicTraits_2`
 */
struct Arr_oblivious_side_tag {};

}

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Tags
 *
 * This type tag is used to indicate that a side of the parameter space, either
 * left, right, bottom, or top, is open, and curves that approach this side
 * might be inserted into the arrangement. More precisely, when this condition
 * occurs at the left, right, bottom, or top side, then the corresponding
 * categories, namely, `Left_side_category`, `Right_side_category`,
 * `Bottom_side_category`, and `Top_side_category`, nested in every geometry
 * traits class, must be convertible to the type `Arr_open_side_tag`. For
 * example, all categories above, nested in every model of the
 * `ArrangementOpenBoundaryTraits_2` concept, must be convertible to
 * `Arr_open_side_tag`, as curves are expected to approach all the four boundary
 * sides of the parameter space (i.e., left, right, bottom, and top).
 *
 * `Arr_oblivious_side_tag` is an empty construct used for dispatching functions
 * based on type of curves that induce the arrangement.
 *
 * \sa `Arr_oblivious_side_tag`
 * \sa `Arr_closed_side_tag`
 * \sa `Arr_contracted_side_tag`
 * \sa `Arr_identified_side_tag`
 * \sa `ArrangementOpenBoundaryTraits_2`
 */
struct Arr_open_side_tag : {};

}

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Tags
 *
 * This type tag is used to indicate that a side of the parameter space, either
 * left, right, bottom, or top, is closed, and curves that reach this side might
 * be inserted into the arrangement. More precisely, when this condition occurs
 * at the left, right, bottom, or top side, then the corresponding categories,
 * namely, `Left_side_category`, `Right_side_category`, `Bottom_side_category`,
 * and `Top_side_category`, nested in every geometry traits class, must be
 * convertible to the type `Arr_closed_side_tag`. At this point none of the
 * traits provided with the \ref PkgArrangementOnSurface2 package supports this
 * condition
 *
 * `Arr_oblivious_side_tag` is an empty construct used for dispatching functions
 * based on type of curves that induce the arrangement.
 *
 * \sa `Arr_oblivious_side_tag`
 * \sa `Arr_open_side_tag`
 * \sa `Arr_contracted_side_tag`
 * \sa `Arr_identified_side_tag`
 * \sa `ArrangementOpenBoundaryTraits_2`
 */
struct Arr_closed_side_tag {};

}

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Tags
 *
 * This type tag is used to indicate that a side of the parameter space, either
 * left, right, bottom, or top, is contracted, and curves that approach this
 * side might be inserted into the arrangement. More precisely, when this
 * condition occurs at the left, right, bottom, or top side, then the
 * corresponding categories, namely, `Left_side_category`,
 * `Right_side_category`, `Bottom_side_category`, and `Top_side_category`,
 * nested in every geometry traits class, must be convertible to the type
 * `Arr_contracted_side_tag`. For example, the `Bottom_side_category` and
 * `Top_side_category` category types, nested in every model of the
 * `ArrangementSphericalBoundaryTraits_2 concept` (such as any instance of the
 * `Arr_geodesic_arc_on_sphere_traits_2` class template) must be convertible to
 * `Arr_contracted_side_tag`
 *
 * `Arr_oblivious_side_tag` is an empty construct used for dispatching functions
 * based on type of curves that induce the arrangement.
 *
 * \sa `Arr_oblivious_side_tag`
 * \sa `Arr_open_side_tag`
 * \sa `Arr_closed_side_tag`
 * \sa `Arr_identified_side_tag`
 * \sa `ArrangementOpenBoundaryTraits_2`
 */
struct Arr_contracted_side_tag {};

}

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Tags
 *
 * This type tag is used to indicate that a side of the parameter space, either
 * left, right, bottom, or top, is identified, and curves that approach this
 * side might be inserted into the arrangement. More precisely, when this
 * condition occurs at the left, right, bottom, or top side, then the
 * corresponding categories, namely, `Left_side_category`,
 * `Right_side_category`, `Bottom_side_category`, and `Top_side_category`,
 * nested in every geometry traits class, must be convertible to the type
 * `Arr_identified_side_tag`. For example, the `Left_side_category` and
 * `Right_side_category` category types, nested in every model of the
 * `ArrangementSphericalBoundaryTraits_2 concept` (such as any instance of the
 * `Arr_geodesic_arc_on_sphere_traits_2` class template) must be convertible to
 * `Arr_identified_side_tag`
 *
 * `Arr_oblivious_side_tag` is an empty construct used for dispatching functions
 * based on type of curves that induce the arrangement.
 *
 * \sa `Arr_oblivious_side_tag`
 * \sa `Arr_open_side_tag`
 * \sa `Arr_closed_side_tag`
 * \sa `Arr_contracted_side_tag`
 * \sa `ArrangementOpenBoundaryTraits_2`
 */
struct Arr_identified_side_tag {};

}
