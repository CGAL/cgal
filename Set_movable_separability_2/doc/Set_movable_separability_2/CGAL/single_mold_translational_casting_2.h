namespace CGAL {
namespace Set_movable_separability_2 {

/*! \ingroup PkgSetMovableSeparability2Funcs
 *
 * Given a simple polygon, this function determines whether a cavity (of a mold
 * in the plane) that has the shape of the polygon can be used so that the
 * polygon could be casted in the mold and then pulled out of the mold without
 * colliding into the mold (but possibly sliding along the mold surface). If the
 * polygon is <em>castable</em>, the function computes the set of top edges of
 * such cavities and the corresponding closed ranges of pull directions.
 * Assuming the top edge normal is parallel to the \f$y\f$-axis, every direction
 * in a range must have a positive component in the positive
 * \f$y\f$-direction. Each top edge and corresponding range is added to a
 * container referred to by a given output iterator.
 *
 * \param[in] pgn the input polygon.
 * \param[out] oi the output iterator. Its value type is a pair, where
 *             (i) the first element in the pair identifies a valid top edge
 *                 represented by its index the type of which is convertible to
                   `size_t`, and
 *             (ii) the second element is a closed range of pull-out directions
 *                  represented as a pair of the extreme directions in the
 *                  range of type `Kernel::Direction_2`.
 * \return the past-the-end iterator of the output container.
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename Kernel, typename OutputIterator>
OutputIterator
single_mold_translational_casting_2(const CGAL::Polygon_2<Kernel>& pgn,
                                    OutputIterator oi);

} /* end namesapce Set_movable_separability_2 */
} /* end namesapce CGAL */
