namespace CGAL {
namespace Set_movable_separability_2 {

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Given a simple polygon and a direction, this function determines whether a
 * cavity (of a mold in the plane) that has the shape of the polygon could be
 * casted in the mold and then pulled out of the mold in the given direction
 * without colliding into the mold (but possibly sliding along the mold
 * surface). If the polygon is not castable at all, the function returns `false`
 * whatsoever.
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param[in] pgn the input polygon.
 * \param[in] d the inspected direction.
 * \param[in] traits the traits to use.
 * \return a pair of elements, where the first is a Boolean that indicates
 *         whether `pgn` can be pulled out in the `d` direction, and the
 *         second is the index of the corresponding top edge in `pgn`.
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
std::pair<bool, size_t>
is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 typename CastingTraits_2::Direction_2& d,
 const CastingTraits_2& traits = CastingTraits_2());

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Same as above with the additional `orientation` argument.
 * If the orientation of the polygon is known upon invocation, specify it.
 * Otherwise, it has to be computed.  Note that finding the orientation of a
 * polygon requires time linear in the number of edges.
 *
 * \param[in] pgn the input polygon.
 * \param[in] d the inspected direction.
 * \param[in] orientation the orientation of `pgn`.
 * \param[in] traits the traits to use.
 * \return a pair of elements, where the first is a Boolean that indicates
 *         whether `pgn` can be pulled out in the `d` direction, and the
 *         second is the index of the corresponding top edge in `pgn`.
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
std::pair<bool, size_t>
is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 typename CastingTraits_2::Direction_2& d,
 CGAL::Orientation orientation,
 const CastingTraits_2& traits = CastingTraits_2());

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Given a simple polygon, an edge of the polygon, and a direction, this
 * function determines whether a cavity (of a mold in the plane) that has the
 * shape of the polygon can be used so that the polygon could be casted in the
 * mold and then pulled out of the mold in the given direction such that the
 * given edge is used as the top edge without colliding into the mold (but
 * possibly sliding along the mold surface). Observe, that if polygon can be
 * pulled out in the given direction, but with a top edge different than the
 * given one, the function returns `false`.
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param[in] pgn the input polygon.
 * \param[in] i the index of an edge in pgn.
 * \param[in] d the tested direction.
 * \param[in] traits the traits to use.
 * \return true if `pgn` can be pulled out in the `d` direction with the
 *         edge identified by `i` being the top edge, and `false` otherwise.
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
bool is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i,
 typename CastingTraits_2::Direction_2& d,
 const CastingTraits_2& traits = CastingTraits_2());

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Same as above with the additional `orientation` argument.
 * If the orientation of the polygon is known upon invocation, specify it.
 * Otherwise, it has to be computed.  Note that finding the orientation of a
 * polygon requires time linear in the number of edges.
 *
 * \param[in] pgn the input polygon.
 * \param[in] i the index of an edge in pgn.
 * \param[in] d the tested direction.
 * \param[in] orientation the orientation of `pgn`.
 * \param[in] traits the traits to use.
 * \return true if `pgn` can be pulled out in the `d` direction with the
 *         edge identified by `i` being the top edge, and `false` otherwise.
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
bool is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i,
 typename CastingTraits_2::Direction_2& d,
 CGAL::Orientation orientation,
 const CastingTraits_2& traits = CastingTraits_2());

} /* end namesapce Set_movable_separability_2 */
} /* end namesapce CGAL */
