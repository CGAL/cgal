namespace CGAL {

/*!
 * \ingroup PkgCasting2Funcs
 * This function computes all top edges of a simple polygon.
 * \param pgn[in] the input polygon.
 * \param oi[out] the output iterator. Its value type is a pair.
 *        a closed range of pull-out directions represented as a pair of the
 *        extreme directions in the
 * \return the past-the-end iterator of the output container.
 * \pre the polygon must be non-degenerate (has at least 3 vertices), simple,
 *      and does not have three consecutive collinear vertices.
 */
template <typename Kernel, typename OutputIterator>
OutputIterator
single_mold_translational_casting_2(const CGAL::Polygon_2<Kernel>& pgn,
                                    OutputIterator oi);

} /* end namesapce CGAL */
