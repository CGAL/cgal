namespace CGAL {
namespace Regularized_boolean_set_operations_2 {

/*! \addtogroup boolean_intersection Intersection Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_intersection
 *
 * There are several overloaded function templates called `intersection()` that
 * compute the \e intersection of two or more polygons and insert the resulting
 * polygons with holes into a container via an output iterator.
 *
 * A function template in this group that accepts two input polygons has one of
 * the following signatures:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>1.1.</b></td><td>`OutputIterator %intersection(const Type1& pgn1, const Type2& pgn2, OutputIterator oi, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>1.2.</b></td><td>`OutputIterator %intersection(const Type1& pgn1, const Type2& pgn2, OutputIterator oi);`</td></tr>
 * <tr><td align="right"><b>2.</b></td><td>`OutputIterator %intersection(const Type1& pgn1, const Type2& pgn2, OutputIterator oi, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * There are also function templates that accept one or two ranges of input polygons:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>3.1.</b></td><td>`OutputIterator intersection(InputIterator begin, InputIterator end, OutputIterator oi, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>3.2.</b></td><td>`OutputIterator intersection(InputIterator begin, InputIterator end, OutputIterator oi);`</td></tr>
 * <tr><td align="right"><b>4.</b></td><td>`OutputIterator intersection(InputIterator begin, InputIterator end, OutputIterator oi, const GpsTraits& traits);`</td></tr>
 * <tr><td align="right"><b>5.1.</b></td><td>`OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, OutputIterator oi, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>5.2.</b></td><td>`OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, OutputIterator oi);`</td></tr>
 * <tr><td align="right"><b>6.</b></td><td>`OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, OutputIterator oi, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * \tparam UsePolylines determines whether the boundaries of the input polygons
 * are treated as cyclic sequences of single (\f$x\f$-monotone) segments or as
 * cyclic sequences of (\f$x\f$-monotone) polylines. If substituted with
 * `CGAL::Tag_true`, which is the default, the input polygons are converted to
 * general polygons bounded by polylines before the operation is actually
 * performed. Then, the resulting general polygons with holes are converted back
 * to standard polygons with holes. If substituted with `CGAL::Tag_false`, the input
 * polygons are used as is. Refer to \ref bso_ssectraits_sel for more information.
 *
 * \param oi the output iterator for the result.
 *
 *   - The types `Type1` and `Type2` of the parameters and of
 *     `InputIterator1::value_type` and `InputIterator2::value_type` must be
 *     convertible to the types specified in a row in the table below,
 *     respectively.  The 3rd column specifies the corresponding dereference
 *     type of the output iterator.
 *
 *   - The types that apply to signatures (<b>1.1.</b>) and (<b>5.1.</b>) above
 *     are restricted to those listed in rows <b>1&ndash;4</b> in the table
 *     below.
 *
 *   - The types that apply to signatures (<b>1.2.</b>) and (<b>5.2.</b>) above
 *     are restricted to those listed in rows <b>5&ndash;8</b> in the table
 *     below.
 *
 *   - The type of `InputIterator::value_type` in (<b>3.1.</b>) above
 *     must be convertible to either `Polygon_2` or `Polygon_with_holes_2`.
 *
 *   - The type of `InputIterator::value_type` in (<b>3.2.</b>) above must be
 *     convertible to either `General_polygon_2` or
 *     `General_polygon_with_holes_2`.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>&nbsp;</th><th>`Type1`</th>                                         <th>`Type2`</th>                                       <th>%Dereference Type of `oi`</th></tr>
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
 * \param (end,begin) a range of input polygons.
 * \param (end1,begin1) the first range of input polygons.
 * \param (end2,begin2) the second range of input polygons.
 * \param oi the output iterator for the result.
 * \param traits an optional traits object.
 *
 * \sa \link boolean_complement `CGAL::Regularized_boolean_operations_2::complement()` \endlink
 * \sa \link boolean_do_intersect `CGAL::Regularized_boolean_operations_2::do_intersect()` \endlink
 * \sa \link boolean_join `CGAL::Regularized_boolean_operations_2::join()` \endlink
 * \sa \link boolean_difference `CGAL::Regularized_boolean_operations_2::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::Regularized_boolean_operations_2::symmetric_difference()` \endlink
 */

/// @{

/*! computes the intersection of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \tparam UsePolylines determines whether the boundary of `pgn` is treated as a
 *         cyclic sequence of single (\f$x\f$-monotone) segments or as a cyclic
 *         sequence of (\f$x\f$-monotone) polylines. If substituted with
 *         `CGAL::Tag_true`, which is the default, `pgn` is converted to a
 *         general polygon bounded by polylines before the operation is actually
 *         performed. Then, the resulting general polygon with holes is
 *         converted back to a standard polygon. If substituted with
 *         `CGAL::Tag_false`, `pgn` is used as is. Refer to \ref
 *         bso_ssectraits_sel for more information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename UsePolylines>
OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            OutputIterator oi,
                            UsePolylines = Tag_true());

/*! computes the intersection of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \tparam UsePolylines determines whether the boundaries of `pgn` is treated as
 *         cyclic sequences of single (\f$x\f$-monotone) segments or as cyclic
 *         sequences of (\f$x\f$-monotone) polylines. If substituted with
 *         `CGAL::Tag_true`, which is the default, `pgn` is converted to a
 *         general polygon with holes bounded by polylines before the operation
 *         is actually performed. Then, the resulting general polygons with
 *         holes are converted back to a standard polygon with holes. If
 *         substituted with `CGAL::Tag_false`, `pgn` is used as is. Refer to
 *         \ref bso_ssectraits_sel for more information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename UsePolylines>
OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            OutputIterator oi,
                            UsePolylines = Tag_true());

/*! computes the intersection of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of two polygons and inserts the resulting polygons
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
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <ctypename Kernel, ctypename Container, ctypename OutputIterator,
          typename UsePolylines>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            OutputIterator oi,
                            UsePolylines = Tag_true());

/*! computes the intersection of two polygons and inserts the resulting polygons
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of two polygons and inserts the resulting polygons
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
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename UsePolylines>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            OutputIterator oi,
                            UsePolylines = Tag_true());

/*! computes the intersection of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator intersection(const General_polygon_2<ArrTraits>& pgn1,
                            const General_polygon_2<ArrTraits>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
intersection(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
             const General_polygon_2<ArrTraits>& pgn2,
             OutputIterator oi);

/*! computes the intersection of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
intersection(const General_polygon_2<ArrTraits>& pgn1,
             const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
             OutputIterator oi);

/*! computes the intersection of two general polygons with holes and inserts the
 * resulting general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<Polygon>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Polygon, typename OutputIterator>
OutputIterator intersection(const General_polygon_with_holes_2<Polygon>& pgn1,
                            const General_polygon_with_holes_2<Polygon>& pgn2,
                            OutputIterator oi);


/*! Given a range of polygons (resp. general polygons) or a range of
 * polygons with holes (resp. general polygons with holes) computes the
 * intersection of all polygons in the range and inserts the resulting polygons
 * with holes (resp. general polygons with holes) into a container via an output
 * iterator.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or `Polygon_with_holes_2`
 *        (resp. `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or `Polygon_with_holes_2`
 *        (resp. `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (resp. `General_polygons_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator, typename OutputIterator>
OutputIterator intersection(InputIterator begin, InputIterator end,
                            OutputIterator oi);

/*! Given a range of polygons (resp. general polygons) or a range of
 * polygons with holes (resp. general polygons with holes) computes the
 * intersection of all polygons in the range and inserts the resulting polygons
 * with holes (resp. general polygons with holes) into a container via an output
 * iterator.
 * \tparam UsePolylines determines whether the boundaries of the polygons in the
 *         input range are treated as cyclic sequences of single
 *         (\f$x\f$-monotone) segments or as a cyclic sequences of
 *         (\f$x\f$-monotone) polylines. If substituted with `CGAL::Tag_true`,
 *         which is the default, the input polygons are converted to general
 *         polygon with holes , bounded by polylines before the operation is
 *         actually performed. Then, the resulting general polygons with holes
 *         are converted back to standard polygons with holes.  If substituted
 *         with `CGAL::Tag_false`, `pgn1` and `pgn2` are used as is. Refer to
 *         \ref bso_ssectraits_sel for more information.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or
 *        `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or
 *        `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (resp. `General_polygons_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator, typename OutputIterator,
          typename UsePolylines>
OutputIterator intersection(InputIterator begin, InputIterator end,
                            OutputIterator oi,
                            UsePolylines = Tag_true());

/*! Given a range of polygons (resp. general polygons) and a range of
 * polygons with holes (resp. general polygons with holes) computes the
 * intersection of all polygons in the two ranges and inserts the resulting
 * polygons with holes (resp. general polygons with holes) into a container via
 * an output iterator.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (resp. `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (resp. `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (resp. `General_polygons_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator1, typename InputIterator2,
typename OutputIterator>
OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1,
                            InputIterator2 begin2, InputIterator2 end2,
                            OutputIterator oi);

/*! Given a range of polygons (resp. general polygons) and a range of
 * polygons with holes (resp. general polygons with holes) computes the
 * intersection of all polygons in the two ranges and inserts the resulting
 * polygons with holes (resp. general polygons with holes) into a container via
 * an output iterator.
 * \tparam UsePolylines determines whether the boundaries of the polygons in the
 *         input ranges are treated as cyclic sequences of single
 *         (\f$x\f$-monotone) segments or as a cyclic sequences of
 *         (\f$x\f$-monotone) polylines. If substituted with `CGAL::Tag_true`,
 *         which is the default, the input polygons are converted to general
 *         polygon with holes , bounded by polylines before the operation is
 *         actually performed. Then, the resulting general polygons with holes
 *         are converted back to standard polygons with holes.  If substituted
 *         with `CGAL::Tag_false`, `pgn1` and `pgn2` are used as is. Refer to
 *         \ref bso_ssectraits_sel for more information.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (resp. `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (resp. `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (resp. `General_polygons_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename UsePolylines>
OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1,
                            InputIterator2 begin2, InputIterator2 end2,
                            OutputIterator oi,
                            UsePolylines = Tag_true());

//////// With Traits

/*! computes the intersection between two polygons and inserts the
 * resulting polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename GpsTraits>
OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            OutputIterator oi,
                            const GpsTraits& traits);

/*! computes the intersection between two polygons and inserts the
 * resulting polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename GpsTraits>
OutputIterator
intersection(const Polygon_2<Kernel, Container>& pgn1,
             const Polygon_with_holes_2<Kernel, Container>& pgn2,
             OutputIterator oi,
             const GpsTraits& traits);

/*! computes the intersection between two polygons and inserts the
 * resulting polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename GpsTraits>
OutputIterator
intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
             const Polygon_2<Kernel, Container>& pgn2,
             OutputIterator oi,
             const GpsTraits& traits);


/*! computes the intersection between two polygons with holes and
 * inserts the resulting polygons with holes into a container via an output
 * iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename GpsTraits>
OutputIterator
intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
             const Polygon_with_holes_2<Kernel, Container>& pgn2,
             OutputIterator oi,
             const GpsTraits& traits);

/*! computes the intersection between two general polygons and inserts
 * the resulting general  polygons with holes into a container via an output
 * iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator intersection(const General_polygon_2<ArrTraits>& pgn1,
                            const General_polygon_2<ArrTraits>& pgn2,
                            OutputIterator oi,
                            const GpsTraits& traits);


/*! computes the intersection between two general polygons and inserts
 * the resulting general  polygons with holes into a container via an output
 * iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator
intersection(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
             const General_polygon_2<ArrTraits>& pgn2,
             OutputIterator oi,
             const GpsTraits& traits);


/*! computes the intersection between two general polygons and inserts
 * the resulting general  polygons with holes into a container via an output
 * iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits>>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator
intersection(const General_polygon_2<ArrTraits>& pgn1,
             const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
             OutputIterator oi,
             const GpsTraits& traits);

/*! computes the intersection between two general polygons and inserts
 * the resulting general  polygons with holes into a container via an output
 * iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<Polygon>`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename Polygon, typename OutputIterator, typename GpsTraits>
OutputIterator
intersection(const General_polygon_with_holes_2<Polygon>& pgn1,
             const General_polygon_with_holes_2<Polygon>& pgn2,
             OutputIterator oi,
             const GpsTraits& traits);

/*! Given a range of polygons (resp. general polygons) or a range of
 * polygons with holes (resp. general polygons with holes) computes the
 * intersection  of all polygons in the range and inserts the resulting
 * polygons with holes (resp. general polygons with holes) into a container via
 * an output iterator.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or
 *        `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or
 *        `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (resp. `General_polygons_with_holes_2`).
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator, typename OutputIterator, typename GpsTraits>
OutputIterator intersection(InputIterator begin, InputIterator end,
                            OutputIterator oi,
                            const GpsTraits& traits);

/*! Given a range of polygons (resp. general polygons) and a range of
 * polygons with holes (resp. general polygons with holes) computes the
 * intersection between all polygons in the two ranges and inserts the
 * resulting polygons with holes (resp. general polygons with holes) into a
 * container via an output iterator.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (resp. `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (resp. `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (resp. `General_polygons_with_holes_2`).
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename GpsTraits>
OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1,
                                    InputIterator2 begin2, InputIterator2 end2,
                                    OutputIterator oi,
                                    const GpsTraits& traits);

/// @}

} } /* namespace CGAL::Regularized_boolean_set_operations_2 */
