namespace CGAL {
namespace Set_movable_separability_2 {

/*! \ingroup PkgSetMovableSeparability2Funcs
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
 * \param[in] pgn the input polygon.
 * \param[in] i the index of an edge in pgn.
 * \param[in] traits the traits to use.
 * \return a pair of elements, where the first is a Boolean that indicates
 *         whether the input edge is a valid top edge, and the second
 *         is a closed range of pullout directions represented as a pair
 *         of the extreme directions in the range. If the input edge is not
 *         a valid top edge, the range is nondeterministic.
 *
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
pullout_directions_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i,
 CastingTraits_2& traits = CastingTraits_2());

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Same as above with the additional `orientation` argument.
 * If the orientation of the polygon is known upon invocation, specify it.
 * Otherwise, it has to be computed.  Note that finding the orientation of a
 * polygon requires time linear in the number of edges.
 *
 * \param[in] pgn the input polygon.
 * \param[in] i the index of an edge in pgn.
 * \param[in] orientation the orientation of `pgn`.
 * \param[in] traits the traits to use.
 * \return a pair of elements, where the first is a Boolean that indicates
 *         whether the input edge is a valid top edge, and the second
 *         is a closed range of pullout directions represented as a pair
 *         of the extreme directions in the range. If the input edge is not
 *         a valid top edge, the range is nondeterministic.
 *
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
pullout_directions_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i,
 CGAL::Orientation orientation,
 CastingTraits_2& traits = CastingTraits_2());

} /* end namesapce Set_movable_separability_2 */
} /* end namesapce CGAL */
