namespace CGAL {
namespace Set_movable_separability_2 {
namespace Single_mold_translational_casting {

/*! \ingroup top_edges_grp
 *
 * Given a simple polygon, this function determines whether a cavity (of a mold
 * in the plane) that has the shape of the polygon can be used so that the
 * polygon could be casted in the mold and then pulled out of the mold without
 * colliding into the mold (but possibly sliding along the mold surface). If the
 * polygon is <em>castable</em>, the function computes the set of top edges of
 * such cavities and the corresponding closed ranges of pullout directions.  Let
 * \f$n\f$ denote the normal to a top edge \f$e\f$, and let \f$d\f$ denote a
 * pullout direction of \f$e\f$. Naturally, the angle between \f$n\f$ and
 * \f$d\f$ must be in the open range (-90&deg;, 90&deg;); that is, \f$n \cdot d
 * > 0\f$. Each top edge and corresponding range is added to a container
 * referred to by a given output iterator.
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param polygon the input polygon.
 * \param oi the output iterator. Its dereference type is a pair, where
 *             (i) the first element in the pair is an iterator to a top edge,
 *                 and
 *             (ii) the second element is a closed range of pullout directions
 *                  represented as a pair of the extreme directions in the
 *                  range of type `CastingTraits_2::Direction_2`.
 * \param traits the traits to use.
 * \return the past-the-end iterator of the output container.
 * \pre `polygon` must be non-degenerate (has at least 3 vertices) and simple,
 * and it does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2, typename OutputIterator>
OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits>& polygon,
                         OutputIterator oi,
                         CastingTraits_2& traits = CastingTraits_2());

/*! \ingroup top_edges_grp
 *
 * Same as above with the additional `orientation` argument.
 * If the orientation of the polygon is known upon invocation, specify it.
 * Otherwise, it has to be computed.  Note that finding the orientation of a
 * polygon requires time linear in the number of edges.
 *
 * \param polygon the input polygon.
 * \param oi the output iterator. Its dereference type is a pair, where
 *             (i) the first element in the pair is an iterator to a top edge,
 *                 and
 *             (ii) the second element is a closed range of pullout directions
 *                  represented as a pair of the extreme directions in the
 *                  range of type `CastingTraits_2::Direction_2`.
 * \param orientation the orientation of `polygon`.
 * \param traits the traits to use.
 * \return the past-the-end iterator of the output container.
 * \pre `polygon` must be non-degenerate (has at least 3 vertices) and simple,
 * and it does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2, typename OutputIterator>
OutputIterator top_edges(const CGAL::Polygon_2<CastingTraits>& polygon,
                         OutputIterator oi,
                         CGAL::Orientation orientation,
                         CastingTraits_2& traits = CastingTraits_2());

} // namespace Single_mold_translational_casting
} // namespace Set_movable_separability_2
} // namespace CGAL
