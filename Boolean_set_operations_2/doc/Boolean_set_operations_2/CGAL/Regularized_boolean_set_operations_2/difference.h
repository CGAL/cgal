namespace CGAL {
namespace Regularized_boolean_set_operations_2 {

/*! \addtogroup boolean_difference Difference Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_difference
 *
 * There are several overloaded function templates called `difference()` that
 * compute the \e difference between two polygons and insert the resulting
 * polygons with holes into a container via an output iterator.
 *
 * A function template in this group has one of the two following signatures:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>1.1.</b></td><td>`OutputIterator difference(const Type1& pgn1, const Type2& pgn2, OutputIterator oi, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>1.2.</b></td><td>`OutputIterator difference(const Type1& pgn1, const Type2& pgn2, OutputIterator oi);`</td></tr>
 * <tr><td align="right"><b>2.</b></td><td>`OutputIterator difference(const Type1& pgn1, const Type2& pgn2, OutputIterator oi, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * \tparam UsePolylines determines whether the boundaries of the input polygons
 * are treated as cyclic sequences of single (\f$x\f$-monotone) segments or as
 * cyclic sequences of (\f$x\f$-monotone) polylines. If substituted with
 * `CGAL::Tag_true`, which is the default, the input polygons are converted to
 * general polygons bounded by polylines before the operation is actually
 * performed. Then, the resulting general polygons with holes are converted back
 * to standard polygons. If substituted with `CGAL::Tag_false`, the input
 * polygons are used as is. Refer to \ref bso_ssectraits_sel for more information.
 *
 *   - The types `Type1` and `Type2` of the parameters must be convertible to the
 * types specified in a row in the table below, respectively.  The 3rd column
 * specifies the corresponding dereference type of the output iterator.
 *   - The types that apply to signature (<b>1.1.</b>) above are restricted to those
 * listed in rows <b>1&ndash;4</b> in the table below.
 *   - The types that apply to signature (<b>1.2.</b>) above are restricted to those
 * listed in rows <b>5&ndash;8</b> in the table below.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>&nbsp;</th><th>`Type1`</th>                        <th>`Type2`</th>                                       <th>%Dereference Type of `oi`</th></tr>
 * <tr><td><b>1</b></td><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_2`</td>                   <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td><b>2</b></td><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td>        <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td><b>3</b></td><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_2`</td>                   <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td><b>4</b></td><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td>        <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td><b>5</b></td><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_2`</td>           <td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td><b>6</b></td><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td><td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td><b>7</b></td><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_2`</td>           <td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td><b>8</b></td><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td><td>`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
 *
 * \param pgn1,pgn2 the input polygons.
 * \param oi the output iterator for the result.
 * \param traits an optional traits object.
 *
 * \sa \link boolean_complement `CGAL::Regularized_boolean_operations_2::complement()` \endlink
 * \sa \link boolean_do_intersect `CGAL::Regularized_boolean_operations_2::do_intersect()` \endlink
 * \sa \link boolean_intersection `CGAL::Regularized_boolean_operations_2::intersection()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::Regularized_boolean_operations_2::symmetric_difference()` \endlink
 */

/// @{

//////// Traits-less

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                          const Polygon_2<Kernel, Container>& pgn2,
                          OutputIterator oi);

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted to a general polygons bounded by polylines
 *         before the operation is actually performed. Then, the resulting
 *         general polygons with holes are converted back to standard
 *         polygons with holes. If substituted with `CGAL::Tag_false`, `pgn1`
 *         and `pgn2` are used as is. Refer to \ref bso_ssectraits_sel for more
 *         information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename UsePolylines>
OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                          const Polygon_2<Kernel, Container>& pgn2,
                          OutputIterator oi,
                          UsePolylines = Tag_true());

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                          const Polygon_with_holes_2<Kernel,Container>& pgn2,
                          OutputIterator oi);

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted to a general polygon and a general polygon
 *         with holes, respectively, bounded by polylines before the operation
 *         is actually performed. Then, the resulting general polygons with
 *         holes are converted back to standard polygons with holes. If
 *         substituted with `CGAL::Tag_false`, `pgn1` and `pgn2` are used as
 *         is. Refer to \ref bso_ssectraits_sel for more information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename UsePolylines>
OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                          const Polygon_with_holes_2<Kernel,Container>& pgn2,
                          OutputIterator oi,
                          UsePolylines = Tag_true());

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                          const Polygon_2<Kernel, Container>& pgn2,
                          OutputIterator oi);

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted to a general polygon with holes and a general
 *         polygon, respectively, bounded by polylines before the operation
 *         is actually performed. Then, the resulting general polygons with
 *         holes are converted back to standard polygons with holes. If
 *         substituted with `CGAL::Tag_false`, `pgn1` and `pgn2` are used as
 *         is. Refer to \ref bso_ssectraits_sel for more information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename UsePolylines>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                          const Polygon_2<Kernel, Container>& pgn2,
                          OutputIterator oi,
                          UsePolylines = Tag_true());

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                          const Polygon_with_holes_2<Kernel, Container>& pgn2,
                          OutputIterator oi);

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted to general polygons with holes, bounded by
 *         polylines before the operation is actually performed. Then, the
 *         resulting general polygons with holes are converted back to standard
 *         polygons with holes. If substituted with `CGAL::Tag_false`, `pgn1`
 *         and `pgn2` are used as is. Refer to \ref bso_ssectraits_sel for more
 *         information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename UsePolylines>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                          const Polygon_with_holes_2<Kernel, Container>& pgn2,
                          OutputIterator oi,
                          UsePolylines = Tag_true());

/*! computes the difference of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator difference(const General_polygon_2<ArrTraits>& pgn1,
                          const General_polygon_2<ArrTraits>& pgn2,
                          OutputIterator oi);

/*! computes the difference of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
           const General_polygon_2<ArrTraits>& pgn2,
           OutputIterator oi);

/*! computes the difference of two polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
difference(const General_polygon_2<ArrTraits>& pgn1,
           const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
           OutputIterator oi);

/*! computes the difference of two general polygons with holes and inserts the
 * resulting general polygons with holes into a container via an output
 * iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<Polygon>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Polygon, typename OutputIterator>
OutputIterator difference(const General_polygon_with_holes_2<Polygon>& pgn1,
                          const General_polygon_with_holes_2<Polygon>& pgn2,
                          OutputIterator oi);

//////// With Traits

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename GpsTraits>
OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                          const Polygon_2<Kernel, Container>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename GpsTraits>
OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                          const Polygon_with_holes_2<Kernel,Container>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename GpsTraits>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                          const Polygon_2<Kernel, Container>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of two polygons with holes and inserts the resulting
 * polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename GpsTraits>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                          const Polygon_with_holes_2<Kernel, Container>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator difference(const General_polygon_2<ArrTraits>& pgn1,
                          const General_polygon_2<ArrTraits>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator
difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
           const General_polygon_2<ArrTraits>& pgn2,
           OutputIterator oi,
           const GpsTraits& traits);

/*! computes the difference of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator
difference(const General_polygon_2<ArrTraits>& pgn1,
           const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
           OutputIterator oi,
           const GpsTraits& traits);

/*! computes the difference of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<Polygon>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename Polygon, typename OutputIterator, typename GpsTraits>
OutputIterator difference(const General_polygon_with_holes_2<Polygon>& pgn1,
                          const General_polygon_with_holes_2<Polygon>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);
/// @}

} } /* namespace CGAL::Regularized_boolean_set_operations_2 */
