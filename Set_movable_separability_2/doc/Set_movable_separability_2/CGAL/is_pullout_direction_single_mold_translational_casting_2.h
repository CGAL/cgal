namespace CGAL {
namespace Set_movable_separability_2 {

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Given a simple polygon and a direction, this function determines whether a
 * cavity (of a mold in the plane) that has the shape of the polygon could be
 * casted in the mold and then pulled out of the mold in the given direction
 * without colliding into the mold (but possibly sliding along the mold
 * surface). Naturally, if the polygon is not castable at all, the function
 * returns `false` whatsoever.
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param[in] pgn the input polygon.
 * \param[in] d the tested direction.
 * \return true if `pgn` can be pulled out in the `d` direction, when rotated
 *         such that the normal to the edge identified by `i` is parallel to
 *         the \f$y\f$-axis.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 typename CastingTraits_2::Direction_2& d);

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Same as above with the additional traits argument.
 * \param[in] pgn the input polygon.
 * \param[in] d the tested direction.
 * \param[in] traits the traits to use.
 * \return true if `pgn` can be pulled out in the `d` direction, when rotated
 *         such that the normal to the edge identified by `i` is parallel to
 *         the \f$y\f$-axis.
 * \pre pgn is castable and the edge identified by `i` is a top edge.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 typename CastingTraits_2::Direction_2& d, CastingTraits_2& traits);

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Given a simple polygon, an edge of the polygon, and a direction, this
 * function determines whether a cavity (of a mold in the plane) that has the
 * shape of the polygon can be used so that the polygon could be casted in the
 * mold and then pulled out of the mold in the given direction without colliding
 * into the mold (but possibly sliding along the mold surface). Observe, that if
 * polygon can be pulled out in the given direction, but with a top edge
 * different than the given one, the function returns `false`.
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param[in] pgn the input polygon.
 * \param[in] i the index of an edge in pgn.
 * \param[in] d the tested direction.
 * \return true if `pgn` can be pulled out in the `d` direction, when rotated
 *         such that the normal to the edge identified by `i` is parallel to
 *         the \f$y\f$-axis.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i,
 typename CastingTraits_2::Direction_2& d);

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Same as above with the additional traits argument.
 * \param[in] pgn the input polygon.
 * \param[in] i the index of an edge in pgn.
 * \param[in] d the tested direction.
 * \param[in] traits the traits to use.
 * \return true if `pgn` can be pulled out in the `d` direction, when rotated
 *         such that the normal to the edge identified by `i` is parallel to
 *         the \f$y\f$-axis.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i,
 typename CastingTraits_2::Direction_2& d, CastingTraits_2& traits);

} /* end namesapce Set_movable_separability_2 */
} /* end namesapce CGAL */
