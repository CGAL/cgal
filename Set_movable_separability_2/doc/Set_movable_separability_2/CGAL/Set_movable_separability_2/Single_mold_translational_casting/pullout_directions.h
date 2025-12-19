namespace CGAL {
namespace Set_movable_separability_2 {
namespace Single_mold_translational_casting {

/*! \ingroup pullout_directions_grp
 *
 * Given a simple polygon and an edge of the polygon, this function determines
 * whether a cavity (of a mold in the plane) that has the shape of the polygon
 * can be used so that the polygon could be casted in the mold using the input
 * edge as the top edge and then pulled out of the mold without colliding
 * into the mold (but possibly sliding along the mold surface). If the polygon
 * is <em>castable</em> this way, the function computes the closed range of
 * pullout directions.
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param polygon the input polygon.
 * \param it an iterator to an edge in polygon.
 * \param traits the traits to use.
 * \return a pair of elements, where the first is a Boolean that indicates
 *         whether the input edge is a valid top edge, and the second
 *         is a closed range of pullout directions represented as a pair
 *         of the extreme directions in the range. If the input edge is not
 *         a valid top edge, the range is nondeterministic.
 *
 * \pre `polygon` must be non-degenerate (has at least 3 vertices) and simple,
 * and it does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
pullout_directions
(const CGAL::Polygon_2<CastingTraits_2>& polygon,
 const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& it,
 CastingTraits_2& traits = CastingTraits_2());

/*! \ingroup pullout_directions_grp
 *
 * Same as above with the additional `orientation` argument.
 * If the orientation of the polygon is known upon invocation, specify it.
 * Otherwise, it has to be computed.  Note that finding the orientation of a
 * polygon requires time linear in the number of edges.
 *
 * \param polygon the input polygon.
 * \param it an iterator to an edge in polygon.
 * \param orientation the orientation of `polygon`.
 * \param traits the traits to use.
 * \return a pair of elements, where the first is a Boolean that indicates
 *         whether the input edge is a valid top edge, and the second
 *         is a closed range of pullout directions represented as a pair
 *         of the extreme directions in the range. If the input edge is not
 *         a valid top edge, the range is nondeterministic.
 *
 * \pre `polygon` must be non-degenerate (has at least 3 vertices) and simple,
 * and it does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
pullout_directions
(const CGAL::Polygon_2<CastingTraits_2>& polygon,
 const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& it,
 CGAL::Orientation orientation,
 CastingTraits_2& traits = CastingTraits_2());

} // namespace Single_mold_translational_casting
} // namespace Set_movable_separability_2
} // namespace CGAL
