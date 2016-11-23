namespace CGAL {
namespace Set_movable_separability_2 {

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Given a simple polygon, an edge of the polygon, and a direction, this
 * function determines whether a cavity (of a mold in the plane) that has the
 * shape of the polygon rotated, such that the normal to the given edge is
 * parallel to the \f$y\f$-axis, can be used so that the polygon could be casted
 * in the mold and then pulled out of the mold in the given direction without
 * colliding into the mold (but possibly sliding along the mold surface). It is
 * required that the polygon is <em>castable</em> this way.
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
 * \pre pgn is castable and the edge identified by `i` is a top edge.
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
 * \pre pgn is castable and the edge identified by `i` is a top edge.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
is_pullout_direction_single_mold_translational_casting_2
(const CGAL::Polygon_2<CastingTraits_2>& pgn, size_t i,
 typename CastingTraits_2::Direction_2& d, CastingTraits_2& traits);

} /* end namesapce Set_movable_separability_2 */
} /* end namesapce CGAL */
