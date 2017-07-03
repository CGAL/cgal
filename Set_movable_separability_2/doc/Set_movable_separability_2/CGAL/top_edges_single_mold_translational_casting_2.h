namespace CGAL {
namespace Set_movable_separability_2 {
namespace Single_mold_translational_casting {

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Given a simple polygon, this function determines whether a cavity (of a mold
 * in the plane) that has the shape of the polygon can be used so that the
 * polygon could be casted in the mold and then pulled out of the mold without
 * colliding into the mold (but possibly sliding along the mold surface). If the
 * polygon is <em>castable</em>, the function computes the set of top edges of
 * such cavities and the corresponding closed ranges of pullout directions.  Let
 * \f$n\f$ denote the normal to a top \f$e\f$, and let \f$d\f$ denote a pullout
 * direction of \f$e\f$. Naturally, the angle between \f$n\f$ and \f$d\f$ must
 * be in the open range (-90&deg;, 90&deg;); that is, \f$n \cdot d > 0\f$. Each
 * top edge and corresponding range is added to a container referred to by a
 * given output iterator.
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param[in] pgn the input polygon.
 * \param[out] oi the output iterator. Its value type is a pair, where
 *             (i) the first element in the pair identifies a valid top edge
 *                 represented by its index the type of which is convertible to
                   `size_t`, and
 *             (ii) the second element is a closed range of pullout directions
 *                  represented as a pair of the extreme directions in the
 *                  range of type `Kernel::Direction_2`.
 * \param[in] traits the traits to use.
 * \return the past-the-end iterator of the output container.
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2, typename OutputIterator>
OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits>& pgn,
                         OutputIterator oi,
                         CastingTraits_2& traits = CastingTraits_2());

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Same as above with the additional `orientation` argument.
 * If the orientation of the polygon is known upon invocation, specify it.
 * Otherwise, it has to be computed.  Note that finding the orientation of a
 * polygon requires time linear in the number of edges.
 *
 * \param[in] pgn the input polygon.
 * \param[out] oi the output iterator. Its value type is a pair, where
 *             (i) the first element in the pair identifies a valid top edge
 *                 represented by its index the type of which is convertible to
                   `size_t`, and
 *             (ii) the second element is a closed range of pullout directions
 *                  represented as a pair of the extreme directions in the
 *                  range of type `Kernel::Direction_2`.
 * \param[in] orientation the orientation of `pgn`.
 * \param[in] traits the traits to use.
 * \return the past-the-end iterator of the output container.
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2, typename OutputIterator>
OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits>& pgn,
                         OutputIterator oi,
                         CGAL::Orientation orientation,
                         CastingTraits_2& traits = CastingTraits_2());

} // namespace Single_mold_translational_casting
} // namesapce Set_movable_separability_2
} // namesapce CGAL
