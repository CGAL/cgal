namespace CGAL {
namespace Regularized_boolean_set_operations_2 {

/*! \addtogroup boolean_complement Complement Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_complement
 *
 * There are several overloaded function templates called `complement()` that
 * compute the \e complement of a given polygon `pgn`. Depending on the type of
 * the polygon `pgn` the complement is either a single (general) polygon with
 * holes, or several (general) poylgons with holes. In the latter case the
 * `complement()` function template inserts the resulting poylgons with holes
 * into a container via an output iterator.
 *
 * A function template in this group has one of the following signatures:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>1.1.</b></td><td>`void complement(const Type1& pgn, Type2& res, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>1.2.</b></td><td>`void complement(const Type1& pgn, Type2& res);`</td></tr>
 * <tr><td align="right"><b>2.</b></td><td>`void complement(const Type1& pgn, Type2& res, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * \tparam UsePolylines determines whether the boundary of the input polygon is
 * treated as a cyclic sequence of single (\f$x\f$-monotone) segments or as a
 * cyclic sequence of (\f$x\f$-monotone) polylines. If substituted with
 * `CGAL::Tag_true`, which is the default, the input polygon is converted to a
 * general polygon bounded by polylines before the operation is actually
 * performed. Then, the resulting general polygon with holes is converted back
 * to a standard polygon. If substituted with `CGAL::Tag_false`, the input
 * polygon is used as is. Refer to \ref bso_ssectraits_sel for more information.
 *
 *   - The types `Type` and `Type2` of the parameters must be convertible to the
 * types specified in a row in the table below, respectively.
 *   - The types that apply to signature (<b>1.1.</b>) above are restricted to those
 * listed in rows <b>1</b> and <b>2</b> in the table below.
 *   - The types that apply to signature (<b>1.2.</b>) above are restricted to those
 * listed in rows <b>3</b> and <b>4</b> in the table below.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>&nbsp;</th><th>`Type1`</th>                                  <th>`Type2`</th></tr>
 * <tr><td><b>1</b></td><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td><b>2</b></td><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td><b>3</b></td><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td></tr>
 * <tr><td><b>4</b></td><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
 *
 * \param pgn the input polygon.
 * \param res the resulting polygon.
 * \param traits an optional traits object.
 *
 * \sa \link boolean_do_intersect `CGAL::Regularized_boolean_operations_2::do_intersect()` \endlink
 * \sa \link boolean_intersection `CGAL::Regularized_boolean_operations_2::intersection()` \endlink
 * \sa \link boolean_join `CGAL::Regularized_boolean_operations_2::join()` \endlink
 * \sa \link boolean_difference `CGAL::Regularized_boolean_operations_2::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::Regularized_boolean_operations_2::symmetric_difference()` \endlink
 */

/// @{

//////// Traits-less

/*! computes the complement of a polygon.
 * \param pgn the input polygon.
 * \param res the resulting complement of \p pgn.
 */
template <typename Kernel, typename Container>
void complement(const Polygon_2<Kernel, Container>& pgn,
                Polygon_with_holes_2<Kernel, Container>& res);

/*! computes the complement of a polygon.
 * \tparam UsePolylines determines whether the boundary of `pgn` is treated as a
 *         cyclic sequence of single (\f$x\f$-monotone) segments or as a cyclic
 *         sequence of (\f$x\f$-monotone) polylines. If substituted with
 *         `CGAL::Tag_true`, which is the default, `pgn` is converted to a
 *         general polygon bounded by polylines before the operation is actually
 *         performed. Then, the resulting general polygon with holes is
 *         converted back to a standard polygon. If substituted with
 *         `CGAL::Tag_false`, `pgn` is used as is. Refer to \ref
 *         bso_ssectraits_sel for more information.
 * \param pgn the input polygon.
 * \param res the resulting complement of \p pgn.
 */
template <typename Kernel, typename Container, typename UsePolylines>
void complement(const Polygon_2<Kernel, Container>& pgn,
                Polygon_with_holes_2<Kernel, Container>& res,
                UsePolylines = Tag_true());

/*! computes the complement of a general polygon.
 * \param pgn the input polygon.
 * \param res the complement of \p pgn.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
void
complement(const General_polygon_2<ArrTraits>& pgn,
           General_polygon_with_holes_2<General_polygon_2<Arr_traits>>& res);

/*! computes the complement of a polygon with holes.
 * \param pgn the input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                          OutputIterator oi);

/*! computes the complement of a polygon with holes.
 * \tparam UsePolylines determines whether the boundaries of `pgn` is treated as
 *         cyclic sequences of single (\f$x\f$-monotone) segments or as cyclic
 *         sequences of (\f$x\f$-monotone) polylines. If substituted with
 *         `CGAL::Tag_true`, which is the default, `pgn` is converted to a
 *         general polygon with holes bounded by polylines before the operation
 *         is actually performed. Then, the resulting general polygons with
 *         holes are converted back to a standard polygon with holes. If
 *         substituted with `CGAL::Tag_false`, `pgn` is used as is. Refer to
 *         \ref bso_ssectraits_sel for more information.
 * \param pgn the input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Traits, typename OutputIterator, typename UsePolylines>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                          OutputIterator oi,
                          UsePolylines = Tag_true());

/*! computes the complement of a general polygon with holes.
 * \param pgn the input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<<General_polygon_2<ArrTraits>>`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
complement(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn,
           OutputIterator oi);

//////// With Traits

/*! computes the complement of a polygon.
 * \param pgn the input polygon.
 * \param res the resulting complement of \p pgn
 * \param traits a traits object.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
void complement(const Polygon_2<Kernel, Container>& pgn,
                Polygon_with_holes_2<Kernel, Container>& res,
                const GpsTraits& traits);

/*! computes the complement of a general polygon.
 * \param pgn the input polygon.
 * \param res the resulting complement of \p pgn
 * \param traits a traits object.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename GpsTraits>
void
complement(const General_polygon_2<ArrTraits>& pgn,
           General_polygon_with_holes_2<General_polygon_2<Arr_traits>>& res,
           const GpsTraits& traits);

/*! computes the complement of a polygon with holes.
 * \param pgn the input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename GpsTraits>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the complement of the general polygon with holes.
 * \param pgn the input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<<General_polygon_2<ArrTraits>>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Polygon, typename OutputIterato, typename GpsTraitsr>
OutputIterator complement(const General_polygon_with_holes_2<Polygon>& pgn,
                          OutputIterator oi,
                          const GpsTraits& traits);

/// @}

} } /* namespace CGAL::Regularized_boolean_set_operations_2 */
