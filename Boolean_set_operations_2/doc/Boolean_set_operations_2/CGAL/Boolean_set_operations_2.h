// CGAL/Boolean_set_operations_2/complement.h
namespace CGAL {

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
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 * \sa \link boolean_intersection `CGAL::intersection()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_difference `CGAL::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
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

} /* namespace CGAL */


// // CGAL/Boolean_set_operations_2/difference.h
namespace CGAL {

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
 * \sa \link boolean_complement `CGAL::complement()` \endlink
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 * \sa \link boolean_intersection `CGAL::intersection()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
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

} /* namespace CGAL */

// CGAL/Boolean_set_operations_2/do_intersect.h
namespace CGAL {

/*! \addtogroup boolean_do_intersect Polygon Intersection Testing Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_do_intersect
 *
 * There are several overloaded function templates called `do_intersect()`
 * that determine whether the interior of two or more polygons intersect.
 *
 * A function template in this group that accepts two input polygons has one of
 * the following signatures:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>1.1.</b></td><td>`bool do_intersect(const Type1& pgn1, const Type2& pgn2, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>1.2.</b></td><td>`bool do_intersect(const Type1& pgn1, const Type2& pgn2);`</td></tr>
 * <tr><td align="right"><b>  2.</b></td><td>`bool do_intersect(const Type1& pgn1, const Type2& pgn2, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * There are also function templates that accept one or two ranges of input polygons:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>3.1.</b></td><td>`bool do_intersect(InputIterator begin, InputIterator end, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>3.2.</b></td><td>`bool do_intersect(InputIterator begin, InputIterator end);`</td></tr>
 * <tr><td align="right"><b>  4.</b></td><td>`bool do_intersect(InputIterator begin, InputIterator end, const GpsTraits& traits);`</td></tr>
 * <tr><td align="right"><b>5.1.</b></td><td>`bool do_intersect(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>5.2.</b></td><td>`bool do_intersect(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2);`</td></tr>
 * <tr><td align="right"><b>  6.</b></td><td>`bool do_intersect(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * \tparam UsePolylines determines whether the boundary of the input polygons
 * are treated as a cyclic sequence of single (\f$x\f$-monotone) segments or as
 * a cyclic sequence of (\f$x\f$-monotone) polylines. If substituted with
 * `CGAL::Tag_true`, which is the default, the input polygons are converted to
 * general polygons bounded by polylines before the operation is actually
 * performed. If substituted with `CGAL::Tag_false`, the input polygons are used
 * as is. Refer to \ref bso_ssectraits_sel for more information.
 *
 *   - The types `Type1` and `Type2` of the parameters of
 *     `InputIterator1::value_type` and `InputIterator2::value_type` must be
 *     convertible to the types specified in a row in the table below,
 *     respectively.
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
 * <tr><th>&nbsp;</th><th>`Type1`</th>                                         <th>`Type2`</th></tr>
 * <tr><td><b>1</b></td><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_2`</td></tr>
 * <tr><td><b>2</b></td><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td><b>3</b></td><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_2`</td></tr>
 * <tr><td><b>4</b></td><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td><b>5</b></td><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_2`</td></tr>
 * <tr><td><b>6</b></td><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td></tr>
 * <tr><td><b>7</b></td><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_2`</td></tr>
 * <tr><td><b>8</b></td><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
 *
 * \param pgn1,pgn2 the input polygons.
 * \param (end,begin) a range of input polygons.
 * \param (end1,begin1) the first range of input polygons.
 * \param (end2,begin2) the second range of input polygons.
 * \param traits an optional traits object.
 *
 * \sa \link boolean_complement `CGAL::complement()` \endlink
 * \sa \link boolean_intersection `CGAL::intersection()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_difference `CGAL::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
 */

/// @{

//////// Traits-less

/*! determines whether two polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2);

/*! determines whether two polygons intersect in their interior.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted to general polygons, bounded by polylines
 *         before the operation is actually performed. If substituted with
 *         `CGAL::Tag_false`, `pgn1` and `pgn2` are used as is. Refer to \ref
 *         bso_ssectraits_sel for more information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2,
                  UsePolylines = Tag_true());

/*! determines whether two polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2);

/*! determines whether two polygons intersect in their interior.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted to a general polygon and a general polygon
 *         with holes, respectively, bounded by polylines before the operation
 *         is actually performed. If substituted with `CGAL::Tag_false`, `pgn1`
 *         and `pgn2` are used as is. Refer to \ref bso_ssectraits_sel for more
 *         information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2,
                  UsePolylines = Tag_true());

/*! determines whether two polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2);

/*! determines whether two polygons intersect in their interior.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted to a general polygon with holes and a general
 *         polygon, respectively, bounded by polylines before the operation
 *         is actually performed. If substituted with `CGAL::Tag_false`, `pgn1`
 *         and `pgn2` are used as is. Refer to \ref bso_ssectraits_sel for more
 *         information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2,
                  UsePolylines = Tag_true());

/*! determines whether two polygons with holes intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2);

/*! determines whether two polygons with holes intersect in their interior.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted to general polygon with holes , bounded by
 *         polylines before the operation is actually performed. If substituted
 *         with `CGAL::Tag_false`, `pgn1` and `pgn2` are used as is. Refer to
 *         \ref bso_ssectraits_sel for more information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2,
                  UsePolylines = Tag_true());

/*! determines whether two general polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
bool do_intersect(const General_polygon_2<ArrTraits>& pgn1,
                  const General_polygon_2<ArrTraits>& pgn2);

/*! determines whether two general polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
bool
do_intersect(const General_polygon_2<ArrTraits>& pgn1,
             const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2);

/*! determines whether two general polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
bool do_intersect(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
                  const General_polygon_2<ArrTraits>& pgn2);

/*! determines whether two general polygons with holes intersect in their
 * interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 */
template <typename Polygon>
bool do_intersect(const General_polygon_with_holes_2<Polygon>& pgn1,
                  const General_polygon_with_holes_2<Polygon>& pgn2);

/*! Given a range of polygons or a range of polygons with holes (respectively a range
 * of general polygons or a range of general polygons with holes) determines
 * whether the open polygons (respectively general polygons) in the range have a common
 * point.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return `true` if the pairwise intersections of all open polygons or polygons
 *         with holes (respectively general polygons or general polygons with holes) in
 *         the range [*begin,*end) overlap, and `false` otherwise.
 */
template <typename InputIterator>
bool do_intersect(InputIterator begin, InputIterator end);

/*! Given a range of polygons or a range of polygons with holes (respectively a range
 * of general polygons or a range of general polygons with holes) determines
 * whether the open polygons (respectively general polygons) in the range have a common
 * point.
 * \tparam UsePolylines determines whether the boundaries of the polygons in the
 *         input range are treated as cyclic sequences of single
 *         (\f$x\f$-monotone) segments or as a cyclic sequences of
 *         (\f$x\f$-monotone) polylines. If substituted with `CGAL::Tag_true`,
 *         which is the default, the input polygons are converted to general
 *         polygon with holes , bounded by polylines before the operation is
 *         actually performed. If substituted with `CGAL::Tag_false`, `pgn1` and
 *         `pgn2` are used as is. Refer to \ref bso_ssectraits_sel for more
 *         information.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return `true` if the pairwise intersections of all open polygons or polygons
 *         with holes (respectively general polygons or general polygons with holes) in
 *         the range [*begin,*end) overlap, and `false` otherwise.
 */
template <typename InputIterator, typename UsePolylines>
bool do_intersect(InputIterator begin, InputIterator end,
                  UsePolylines = Tag_true());

/*! Given a range of polygons (respectively general polygons) and a range of polygons
 * with holes (respectively general polygons with holes) determines whether the open
 * polygons (respectively general polygons) in the two ranges have a common point.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return `true` if the pairwise intersections of all open polygons (respectively
 *        general polygons) and polygons with holes (respectively general polygons with
 *        holes) in the ranges [*begin1,*end1) and [*begin2,*end2),
 *        respectively, overlap, and `false` otherwise.
 */
template <typename InputIterator1, typename InputIterator2>
bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                  InputIterator2 begin2, InputIterator2 end2);

/*! Given a range of polygons (respectively general polygons) and a range of polygons
 * with holes (respectively general polygons with holes) determines whether the open
 * polygons (respectively general polygons) in the two ranges have a common point.
 * \tparam UsePolylines determines whether the boundaries of the polygons in the
 *         input ranges are treated as cyclic sequences of single
 *         (\f$x\f$-monotone) segments or as a cyclic sequences of
 *         (\f$x\f$-monotone) polylines. If substituted with `CGAL::Tag_true`,
 *         which is the default, the input polygons are converted to general
 *         polygon with holes , bounded by polylines before the operation is
 *         actually performed. If substituted with `CGAL::Tag_false`, `pgn1` and
 *         `pgn2` are used as is. Refer to \ref bso_ssectraits_sel for more
 *         information.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return `true` if the pairwise intersections of all open polygons (respectively
 *        general polygons) and polygons with holes (respectively general polygons with
 *        holes) in the ranges [*begin1,*end1) and [*begin2,*end2),
 *        respectively, overlap, and `false` otherwise.
 */
template <typename InputIterator1, typename InputIterator2,
          typename UsePolylines>
bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                  InputIterator2 begin2, InputIterator2 end2,
                  UsePolylines = Tag_true());

//////// With Traits

/*! determines whether two polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2,
                  const GpsTraits& traits);

/*! determines whether two polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2,
                  const GpsTraits& traits,
                  const GpsTraits& traits);

/*! determines whether two polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2,
                  const GpsTraits& traits);

/*! determines whether two polygons with holes intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2,
                  const GpsTraits& traits);

/*! determines whether two general polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename GpsTraits>
bool do_intersect(const General_polygon_2<ArrTraits>& pgn1,
                  const General_polygon_2<ArrTraits>& pgn2,
                  const GpsTraits& traits);

/*! determines whether two general polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename GpsTraits>
bool
do_intersect(const General_polygon_2<ArrTraits>& pgn1,
             const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
             const GpsTraits& traits);

/*! determines whether two general polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename GpsTraits>
bool
do_intersect(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
             const General_polygon_2<ArrTraits>& pgn2,
             const GpsTraits& traits);

/*! determines whether two general polygons with holes intersect in their
 * interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interior and `false`
 *         otherwise.
 *
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Polygon, typename GpsTraits>
bool do_intersect(const General_polygon_with_holes_2<Polygon>& pgn1,
                  const General_polygon_with_holes_2<Polygon>& pgn2,
                  const GpsTraits& traits);

/*! Given a range of polygons or a range of polygons with holes (respectively a range
 * of general polygons or a range of general polygons with holes) determines
 * whether the open polygons (respectively general polygons) in the range have a common
 * point.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param traits a traits object.
 * \return `true` if the pairwise intersections of all open polygons or polygons
 *         with holes (respectively general polygons or general polygons with holes) in
 *         the range [*begin,*end) overlap, and `false` otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator, typename GpsTraits>
bool do_intersect(InputIterator begin, InputIterator end,
                  const GpsTraits& traits);

/*! Given a range of polygons (respectively general polygons) and a range of polygons
 * with holes (respectively general polygons with holes) determines whether the open
 * polygons (respectively general polygons) in the two ranges have a common point.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param traits a traits object.
 * \return `true` if the pairwise intersections of all open polygons (respectively
 *        general polygons) and polygons with holes (respectively general polygons with
 *        holes) in the ranges [*begin1,*end1) and [*begin2,*end2),
 *        respectively, overlap, and `false` otherwise.
 */
template <typename InputIterator1, typename InputIterator2, typename GpsTraits>
bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                  InputIterator2 begin2, InputIterator2 end2,
                  const GpsTraits& traits);

/// @}
} /* namespace CGAL */

// CGAL/Boolean_set_operations_2/intersection.h
namespace CGAL {

/*! \addtogroup boolean_intersection Polygon Intersection Functions
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
 * \sa \link boolean_complement `CGAL::complement()` \endlink
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_difference `CGAL::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
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


/*! Given a range of polygons (respectively general polygons) or a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * intersection of all polygons in the range and inserts the resulting polygons
 * with holes (respectively general polygons with holes) into a container via an output
 * iterator.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or `Polygon_with_holes_2`
 *        (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or `Polygon_with_holes_2`
 *        (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator, typename OutputIterator>
OutputIterator intersection(InputIterator begin, InputIterator end,
                            OutputIterator oi);

/*! Given a range of polygons (respectively general polygons) or a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * intersection of all polygons in the range and inserts the resulting polygons
 * with holes (respectively general polygons with holes) into a container via an output
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
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator, typename OutputIterator,
          typename UsePolylines>
OutputIterator intersection(InputIterator begin, InputIterator end,
                            OutputIterator oi,
                            UsePolylines = Tag_true());

/*! Given a range of polygons (respectively general polygons) and a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * intersection of all polygons in the two ranges and inserts the resulting
 * polygons with holes (respectively general polygons with holes) into a container via
 * an output iterator.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator1, typename InputIterator2,
typename OutputIterator>
OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1,
                            InputIterator2 begin2, InputIterator2 end2,
                            OutputIterator oi);

/*! Given a range of polygons (respectively general polygons) and a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * intersection of all polygons in the two ranges and inserts the resulting
 * polygons with holes (respectively general polygons with holes) into a container via
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
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
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

/*! Given a range of polygons (respectively general polygons) or a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * intersection  of all polygons in the range and inserts the resulting
 * polygons with holes (respectively general polygons with holes) into a container via
 * an output iterator.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator, typename OutputIterator, typename GpsTraits>
OutputIterator intersection(InputIterator begin, InputIterator end,
                            OutputIterator oi,
                            const GpsTraits& traits);

/*! Given a range of polygons (respectively general polygons) and a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * intersection between all polygons in the two ranges and inserts the
 * resulting polygons with holes (respectively general polygons with holes) into a
 * container via an output iterator.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
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

} /* namespace CGAL */

// CGAL/Boolean_set_operations_2/join.h
namespace CGAL {

/*! \addtogroup boolean_join Union Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_union
 *
 * There are several overloaded function templates called `join()` that
 * compute the \e union of two polygons.
 *
 * A function template in this group has one of the following signatures:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>1.1.</b></td><td>`bool join(const Type1& pgn1, const Type2& pgn2, Type3& res, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>1.2.</b></td><td>`bool join(const Type1& pgn1, const Type2& pgn2, Type3& res);`</td></tr>
 * <tr><td align="right"><b>2.</b></td><td>`bool join(const Type1& pgn1, const Type2& pgn2, Type3& res, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * There are also function templates that accept one or two ranges of input polygons:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>3.1.</b></td><td>`OutputIterator join(InputIterator begin, InputIterator end, OutputIterator oi, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>3.2.</b></td><td>`OutputIterator join(InputIterator begin, InputIterator end, OutputIterator oi);`</td></tr>
 * <tr><td align="right"><b>4.</b></td><td>`OutputIterator join(InputIterator begin, InputIterator end, OutputIterator oi, const GpsTraits& traits);`</td></tr>
 * <tr><td align="right"><b>5.1.</b></td><td>`OutputIterator join(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, OutputIterator oi, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>5.2.</b></td><td>`OutputIterator join(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, OutputIterator oi);`</td></tr>
 * <tr><td align="right"><b>6.</b></td><td>`OutputIterator join(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, OutputIterator oi, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * \tparam UsePolylines determines whether the boundaries of the input polygons
 * are treated as cyclic sequences of single (\f$x\f$-monotone) segments or as
 * cyclic sequences of (\f$x\f$-monotone) polylines. If substituted with
 * `CGAL::Tag_true`, which is the default, the input polygons are converted to
 * general polygons bounded by polylines before the operation is actually
 * performed. Then, the resulting general polygon with holes is converted back
 * to a standard polygon. If substituted with `CGAL::Tag_false`, the input
 * polygons are used as is. Refer to \ref bso_ssectraits_sel for more information.
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
 * <tr><th>&nbsp;</th><th>`Type1`</th>                                       <th>`Type2`</th>                                       <th>`Type3`</th></tr>
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
 * \param res the resulting polygon.
 * \param oi the output iterator for the result.
 * \param traits an optional traits object.
 *
 * \sa \link boolean_complement `CGAL::complement()` \endlink
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 * \sa \link boolean_intersection `CGAL::intersection()` \endlink
 * \sa \link boolean_difference `CGAL::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
 */

/// @{

//////// Traits-less

/*! computes the union of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_2<Kernel, Container>& pgn1,
          const Polygon_2<Kernel, Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res);

/*! computes the union of two polygons.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted into a general polygons bounded by polylines
 *         before the operation is actually performed. Then, the resulting
 *         general polygons with holes are converted back to standard
 *         polygons with holes. If substituted with `CGAL::Tag_false`, `pgn1`
 *         and `pgn2` are used as is. Refer to \ref bso_ssectraits_sel for more
 *         information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool join(const Polygon_2<Kernel, Container>& pgn1,
          const Polygon_2<Kernel, Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res,
          UsePolylines = Tag_true());

/*! computes the union of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_2<Kernel, Container>& pgn1,
          const Polygon_with_holes_2<Kernel,Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res);

/*! computes the union of two polygons.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted into a general polygon and a general polygon
 *         with holes, respectively, bounded by polylines before the operation
 *         is actually performed. Then, the resulting general polygons with
 *         holes are converted back to standard polygons with holes. If
 *         substituted with `CGAL::Tag_false`, `pgn1` and `pgn2` are used as
 *         is. Refer to \ref bso_ssectraits_sel for more information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool join(const Polygon_2<Kernel, Container>& pgn1,
          const Polygon_with_holes_2<Kernel,Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res,
          UsePolylines = Tag_true());

/*! computes the union of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_with_holes_2<Kernel, Container>& pgn1,
          const Polygon_2<Kernel, Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res);

/*! computes the union of two polygons.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted into a general polygon with holes and a general
 *         polygon, respectively, bounded by polylines before the operation
 *         is actually performed. Then, the resulting general polygons with
 *         holes are converted back to standard polygons with holes. If
 *         substituted with `CGAL::Tag_false`, `pgn1` and `pgn2` are used as
 *         is. Refer to \ref bso_ssectraits_sel for more information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool join(const Polygon_with_holes_2<Kernel, Container>& pgn1,
          const Polygon_2<Kernel, Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res,
          UsePolylines = Tag_true());

/*! computes the union of two polygons with holes.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_with_holes_2<Kernel, Container>& pgn1,
          const Polygon_with_holes_2<Kernel, Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res);

/*! computes the union of two polygons with holes.
 * \tparam UsePolylines determines whether the boundaries of `pgn1` and `pgn2`
 *         are treated as cyclic sequences of single (\f$x\f$-monotone) segments
 *         or as a cyclic sequences of (\f$x\f$-monotone) polylines. If
 *         substituted with `CGAL::Tag_true`, which is the default, `pgn1` and
 *         `pgn2` are converted into general polygons with holes, bounded by
 *         polylines before the operation is actually performed. Then, the
 *         resulting general polygons with holes are converted back to standard
 *         polygons with holes. If substituted with `CGAL::Tag_false`, `pgn1`
 *         and `pgn2` are used as is. Refer to \ref bso_ssectraits_sel for more
 *         information.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool join(const Polygon_with_holes_2<Kernel, Container>& pgn1,
          const Polygon_with_holes_2<Kernel, Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res,
          UsePolylines = Tag_true());

/*! computes the union of two general polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
bool join(const General_polygon_2<ArrTraits>& pgn1,
          const General_polygon_2<ArrTraits>& pgn2,
          General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& res);

/*! computes the union of two general polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
bool
join(const General_polygon_2<ArrTraits>& pgn1,
     const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& res);

/*! computes the union of two general polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
bool
join(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
     const General_polygon_2<ArrTraits>& pgn1,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& res);

/*! computes the union of two general polygons with holes.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Polygon>
bool join(const General_polygon_with_holes_2<Polygon>& pgn1,
          const General_polygon_with_holes_2<Polygon>& pgn2,
          General_polygon_with_holes_2<Polygon>& res);

/*! Given a range of polygons (respectively general polygons) or a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * union of all polygons in the range and inserts the resulting polygons
 * with holes (respectively general polygons with holes) into a container via an output
 * iterator.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator, typename OutputIterator>
OutputIterator join(InputIterator begin, InputIterator end,
                    OutputIterator oi);

/*! Given a range of polygons (respectively general polygons) or a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * union of all polygons in the range and inserts the resulting polygons
 * with holes (respectively general polygons with holes) into a container via an output
 * iterator.
 * \tparam UsePolylines determines whether the boundaries of the polygons in the
 *         input range are treated as cyclic sequences of single
 *         (\f$x\f$-monotone) segments or as a cyclic sequences of
 *         (\f$x\f$-monotone) polylines. If substituted with `CGAL::Tag_true`,
 *         which is the default, the input polygons are converted into general
 *         polygon with holes , bounded by polylines before the operation is
 *         actually performed. Then, the resulting general polygons with holes
 *         are converted back to standard polygons with holes.  If substituted
 *         with `CGAL::Tag_false`, `pgn1` and `pgn2` are used as is. Refer to
 *         \ref bso_ssectraits_sel for more information.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator, typename OutputIterator,
          typename UsePolylines>
OutputIterator join(InputIterator begin, InputIterator end,
                    OutputIterator oi,
                    UsePolylines = Tag_true());

/*! Given a range of polygons (respectively general polygons) and a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * union of all polygons in the two ranges and inserts the resulting
 * polygons with holes (respectively general polygons with holes) into a container via
 * an output iterator.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                    InputIterator2 begin2, InputIterator2 end2,
                    OutputIterator oi);

//////// With Traits

/*! computes the union of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \param traits a traits object.
 * \return `true` if the two input polygons overlap.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
bool join(const Polygon_2<Kernel, Container>& pgn1,
          const Polygon_2<Kernel, Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res,
          const GpsTraits& traits);


/*! computes the union of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \param traits a traits object.
 * \return `true` if the two input polygons overlap.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
bool join(const Polygon_2<Kernel, Container>& pgn1,
          const Polygon_with_holes_2<Kernel,Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res,
          const GpsTraits& traits);


/*! computes the union of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \param traits a traits object.
 * \return `true` if the two input polygons overlap.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
bool join(const Polygon_with_holes_2<Kernel, Container>& pgn2,
          const Polygon_2<Kernel, Container>& pgn1,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res,
          const GpsTraits& traits);


/*! computes the union of two polygons with holes.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \param traits a traits object.
 * \return `true` if the two input polygons overlap.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
bool join(const Polygon_with_holes_2<Kernel, Container>& pgn2,
          const Polygon_with_holes_2<Kernel, Container>& pgn1,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container>>& res,
          const GpsTraits& traits);


/*! computes the union of two general polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \param traits a traits object.
 * \return `true` if the two input polygons overlap.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename GpsTraits>
bool join(const General_polygon_2<ArrTraits>& pgn1,
          const General_polygon_2<ArrTraits>& pgn2,
          General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& res,
          const GpsTraits& traits);


/*! computes the union of two general polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \param traits a traits object.
 * \return `true` if the two input polygons overlap.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename GpsTraits>
bool
join(const General_polygon_2<ArrTraits>& pgn1,
     const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& res,
     const GpsTraits& traits);


/*! computes the union of two general polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \param traits a traits object.
 * \return `true` if the two input polygons overlap.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename GpsTraits>
bool join(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
          const General_polygon_2<ArrTraits>& pgn1,
          General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& res,
          const GpsTraits& traits);


/*! computes the union of two general polygons with holes.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \param traits a traits object.
 * \return `true` if the two input polygons overlap.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Polygon, typename GpsTraits>
bool join(const General_polygon_with_holes_2<Polygon>& pgn1,
          const General_polygon_with_holes_2<Polygon>& pgn2,
          General_polygon_with_holes_2<Polygon>& res,
          const GpsTraits& traits);

/*! Given a range of polygons (respectively general polygons) or a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * union of all polygons in the range and inserts the resulting polygons
 * with holes (respectively general polygons with holes) into a container via an output
 * iterator.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator, typename OutputIterator, typename GpsTraits>
OutputIterator join(InputIterator begin, InputIterator end,
                    OutputIterator oi,
                    const GpsTraits& traits);

/*! Given a range of polygons (respectively general polygons) and a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * union of all polygons in the two ranges and inserts the resulting
 * polygons with holes (respectively general polygons with holes) into a container via
 * an output iterator.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename GpsTraits>
OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                    InputIterator2 begin2, InputIterator2 end2,
                    OutputIterator oi,
                    const GpsTraits& traits);

/*! Given a range of polygons (respectively general polygons) and a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * union of all polygons in the two ranges and inserts the resulting
 * polygons with holes (respectively general polygons with holes) into a container via
 * an output iterator.
 * \tparam UsePolylines determines whether the boundaries of the polygons in the
 *         input ranges are treated as cyclic sequences of single
 *         (\f$x\f$-monotone) segments or as a cyclic sequences of
 *         (\f$x\f$-monotone) polylines. If substituted with `CGAL::Tag_true`,
 *         which is the default, the input polygons are converted into general
 *         polygon with holes , bounded by polylines before the operation is
 *         actually performed. Then, the resulting general polygons with holes
 *         are converted back to standard polygons with holes.  If substituted
 *         with `CGAL::Tag_false`, `pgn1` and `pgn2` are used as is. Refer to
 *         \ref bso_ssectraits_sel for more information.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator1, typename InputIterator2,
class OutputIterator, typename UsePolylines>
OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                    InputIterator2 begin2, InputIterator2 end2,
                    OutputIterator oi,
                    UsePolylines = Tag_true());

/// @}

} /* namespace CGAL */

// CGAL/Boolean_set_operations_2/oriented_side.h
namespace CGAL {

/*! \addtogroup boolean_oriented_side Oriented Side Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_oriented_side
 *
 * There are several overloaded function templates called `Oriented_side()`
 * that compute the relative position of either (i) two polygons or (ii) a
 * point and a polygon. This group of function templates is divided into two
 * subgroups.
 *
 * \cgalHeading{Oriented Side of two Polygons}
 *
 * Every function template in the first subgroup accepts two polygons `pgn1` and
 * `pgn2`.  It returns `ON_POSITIVE_SIDE` if the two given polygons `pgn1` and
 * `pgn2` intersect in their interiors, `ON_NEGATIVE_SIDE` if `pgn1` and `pgn2`
 * do not intersect at all, and `ON_ORIENTED_BOUNDARY` if `pgn1` and `pgn2`
 * intersect only in their boundaries.
 *
 * A function template in this subgroup has one of the two following signatures:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>1.1.</b></td><td>`Oriented_side oriented_side(const Type1& pgn1, const Type2& pgn2, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>1.2.</b></td><td>`Oriented_side oriented_side(const Type1& pgn1, const Type2& pgn2);`</td></tr>
 * <tr><td align="right"><b>2.</b></td><td>`Oriented_side oriented_side(const Type1& pgn1, const Type2& pgn2, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * \cgalHeading{Oriented Side of a Point and a Polygon}
 *
 * The functions in the second group accept a point `p` and a polygon `pgn`.
 * Each function in this group returns `ON_POSITIVE_SIDE` if the point `p`
 * is in the interior of `pgn`, `ON_NEGATIVE_SIDE` if `p` is in the exterior
 * of `pgn`, and `ON_ORIENTED_BOUNDARY` if `p` is on the boundary of `pgn`.
 *
 * A function in this subgroup has one of the following signatures:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>3.1.</b></td><td>`Oriented_side oriented_side(const Point_2& p, const Type& pgn, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>3.2.</b></td><td>`Oriented_side oriented_side(const Point_2& p, const Type& pgn);`</td></tr>
 * <tr><td align="right"><b>  4.</b></td><td>`Oriented_side oriented_side(const Point_2& p, const Type& pgn, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * \tparam UsePolylines determines whether the boundaries of the input polygons
 * are treated as cyclic sequences of single (\f$x\f$-monotone) segments or as
 * cyclic sequences of (\f$x\f$-monotone) polylines. If substituted with
 * `CGAL::Tag_true`, which is the default, the input polygons are converted to
 * general polygons bounded by polylines before the operation is actually
 * performed. If substituted with `CGAL::Tag_false`, the input polygons are used
 * as is. Refer to \ref bso_ssectraits_sel for more information.
 *
 *   - The types `Type1` and `Type2` of the parameters must be convertible to the
 * types specified in a row in the following table, respectively.
 *   - The types that apply to signature (<b>1.1.</b>) above are restricted to those
 * listed in rows <b>1&ndash;4</b> in the table below.
 *   - The types that apply to signature (<b>1.2.</b>) above are restricted to those
 * listed in rows <b>5&ndash;8</b> in the table below.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>&nbsp;</th><th>Type1</th>                                         <th>Type2</th></tr>
 * <tr><td><b>1</b></td><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_2`</td></tr>
 * <tr><td><b>2</b></td><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td><b>3</b></td><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_2`</td></tr>
 * <tr><td><b>4</b></td><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td><b>5</b></td><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_2`</td></tr>
 * <tr><td><b>6</b></td><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td></tr>
 * <tr><td><b>7</b></td><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_2`</td></tr>
 * <tr><td><b>8</b></td><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
 *
 *   - In (<b>3.1.</b>) `Type` must be convertible to either `Polygon_2` or `Polygon_with_holes_2`.
 *   - In (<b>3.2.</b>) `Type` must be convertible to either `General_polygon_2` or `General_polygon_with_holes_2`.
 *   - In (<b>4.</b>) `Type` must be convertible to any of the four types above.
 *
 * \param p the input point.
 * \param pgn the input polygon.
 * \param pgn1,pgn2 the input polygons.
 * \param traits an optional traits object.
 *
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 */

/// @{

// Polygon--Polygon Traits-less

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 */
template <typename Kernel, typename Container>
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2);

/*! computes the relative position of two polygons.
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
 */
template <typename Kernel, typename Container, typename UsePolylines>
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            UsePolylines = Tag_true());

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 */
template <typename Kernel, typename Container>
Oriented_side
oriented_side(const Polygon_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2);

/*! computes the relative position of two polygons.
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
 */
template <typename Kernel, typename Container, typename UsePolylines>
Oriented_side
oriented_side(const Polygon_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2,
              UsePolylines = Tag_true());

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 */
template <typename Kernel, typename Container>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2);

/*! computes the relative position of two polygons.
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
 */
template <typename Kernel, typename Container, typename UsePolylines>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            UsePolylines = Tag_true());

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 */
template <typename Kernel, typename Container>
Oriented_side
oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2);

/*! computes the relative position of two polygons.
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
 */
template <typename Kernel, typename Container, typename UsePolylines>
Oriented_side
oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2,
              UsePolylines = Tag_true());

/*! computes the relative position of two polygons.
 * \param pgn1 1st the input polygon.
 * \param pgn2 the 2nd input polygon.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
Oriented_side oriented_side(const General_polygon_2<ArrTraits>& pgn1,
                            const General_polygon_2<ArrTraits>& pgn2);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
Oriented_side
oriented_side(const General_polygon_2<ArrTraits>& pgn1,
              const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
Oriented_side
oriented_side(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
              const General_polygon_2<ArrTraits>& pgn2);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 */
template <typename Polygon>
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon>& pgn1,
                            const General_polygon_with_holes_2<Polygon>& pgn2);

//////// Polygon--Polygon With Traits

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            const GpsTraits& traits);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            const GpsTraits& traits);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            const GpsTraits& traits);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            const GpsTraits& traits);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 *
 */
template <typename ArrTraits, typename GpsTraits>
Oriented_side oriented_side(const General_polygon_2<ArrTraits>& pgn1,
                            const General_polygon_2<ArrTraits>& pgn2,
                            const GpsTraits& traits);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 *
 */
template <typename ArrTraits, typename GpsTraits>
Oriented_side
oriented_side(const General_polygon_2<ArrTraits>& pgn1,
              const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
              const GpsTraits& traits);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 *
 */
template <typename ArrTraits, typename GpsTraits>
Oriented_side
oriented_side(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
              const General_polygon_2<ArrTraits>& pgn2,
              const GpsTraits& traits);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Polygon, typename GpsTraits>
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon>& pgn1,
                            const General_polygon_with_holes_2<Polygon>& pgn2,
                            const GpsTraits& traits);

// Point--Polygon Traits-less

/*! computes the relative position of a point and a polygon.
 * \param p the input point.
 * \param pgn the input polygon.
 */
template <typename Kernel, typename Container>
Oriented_side oriented_side(const Point_2& p,
                            const Polygon_2<Kernel, Container>& pgn);

/*! computes the relative position of a point and a polygon.
 * \param p the input point.
 * \param pgn the input polygon.
 */
template <typename Kernel, typename Container>
Oriented_side oriented_side(const Point_2& p,
                            const Polygon_with_holes_2<Kernel, Container>& pgn);

/*! computes the relative position of a point and a polygon.
 * \param p the input point.
 * \param pgn the input polygon.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
Oriented_side oriented_side(const Point_2& p,
                            const General_polygon_2<ArrTraits>& pgn);

/*! computes the relative position of a point and a polygon.
 * \param p the input point.
 * \param pgn the input polygon.
 */
template <typename Polygon>
Oriented_side oriented_side(const Point_2& p,
                            const General_polygon_with_holes_2<Polygon>& pgn);

//////// Point--Polygon With Traits

/*! computes the relative position of a point and a polygon.
 * \param p the input point.
 * \param pgn the input polygon.
 * \param traits a traits object.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Point_2& p,
                            const Polygon_2<Kernel, Container>& pgn,
                            const GpsTraits& traits);

/*! computes the relative position of a point and a polygon.
 * \param p the input point.
 * \param pgn the input polygon.
 * \param traits a traits object.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Point_2& p,
                            const Polygon_with_holes_2<Kernel, Container>& pgn,
                            const GpsTraits& traits);

/*! computes the relative position of a point and a polygon.
 * \param p the input point.
 * \param pgn the input polygon.
 * \param traits a traits object.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename GpsTraits>
Oriented_side oriented_side(const Point_2& p,
                            const General_polygon_2<ArrTraits>& pgn,
                            const GpsTraits& traits);

/*! computes the relative position of a point and a polygon.
 * \param p the input point.
 * \param pgn the input polygon.
 * \param traits a traits object.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Polygon, typename GpsTraits>
Oriented_side oriented_side(const Point_2& p,
                            const General_polygon_with_holes_2<Polygon>& pgn,
                            const GpsTraits& traits);

/// @}

} /* namespace CGAL */

// CGAL/Boolean_set_operations_2/symmetric_difference.h
namespace CGAL {

/*! \addtogroup boolean_symmetric_difference Symmetric Difference Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_symmetric_difference
 *
 * There are several overloaded function templates called
 * `symmetric_difference()` that compute the \e symmetric difference
 * between two or more input polygons and insert the resulting
 * polygons with holes into a container via an output iterator.
 *
 * A function template in this group that accepts two input polygons has one of
 * the following signatures:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>1.1.</b></td><td>`OutputIterator symmetric_difference(const Type1& pgn1, const Type2& pgn2, OutputIterator oi, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>1.2.</b></td><td>`OutputIterator symmetric_difference(const Type1& pgn1, const Type2& pgn2, OutputIterator oi);`</td></tr>
 * <tr><td align="right"><b>2.</b></td><td>`OutputIterator symmetric_difference(const Type1& pgn1, const Type2& pgn2, OutputIterator oi, const GpsTraits& traits);`</td></tr>
 * </table>
 *
 * There are also function templates that accept one or two ranges of input polygons:
 * <table cellpadding=3 border="0">
 * <tr><td align="right"><b>3.1.</b></td><td>`OutputIterator symmetric_difference(InputIterator begin, InputIterator end, OutputIterator oi, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>3.2.</b></td><td>`OutputIterator symmetric_difference(InputIterator begin, InputIterator end, OutputIterator oi);`</td></tr>
 * <tr><td align="right"><b>4.</b></td><td>`OutputIterator symmetric_difference(InputIterator begin, InputIterator end, OutputIterator oi, const GpsTraits& traits);`</td></tr>
 * <tr><td align="right"><b>5.1.</b></td><td>`OutputIterator symmetric_difference(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, OutputIterator oi, UsePolylines = Tag_true());`</td></tr>
 * <tr><td align="right"><b>5.2.</b></td><td>`OutputIterator symmetric_difference(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, OutputIterator oi);`</td></tr>
 * <tr><td align="right"><b>6.</b></td><td>`OutputIterator symmetric_difference(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, InputIterator2 end2, OutputIterator oi, const GpsTraits& traits);`</td></tr>
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
 * \sa \link boolean_complement `CGAL::complement()` \endlink
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 * \sa \link boolean_intersection `CGAL::intersection()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_difference `CGAL::difference()` \endlink
 */

/// @{

//////// Traits-less

/*! computes the symmetric difference between two polygons and inserts the
 * resulting polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                                    const Polygon_2<Kernel, Container>& pgn2,
                                    OutputIterator oi);

/*! computes the symmetric difference between two polygons and inserts the
 * resulting polygons with holes into a container via an output iterator.
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
template <typename Kernel, typename Container, typename OutputIterator,
          typename UsePolylines>
OutputIterator symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                                    const Polygon_2<Kernel, Container>& pgn2,
                                    OutputIterator oi,
                                    UsePolylines = Tag_true());

/*! computes the symmetric difference between two polygons and inserts the
 * resulting polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator
symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi);

/*! computes the symmetric difference between two polygons and inserts the
 * resulting polygons with holes into a container via an output iterator.
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
OutputIterator
symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi,
                     UsePolylines = Tag_true());

/*! computes the symmetric difference between two polygons and inserts the
 * resulting polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_2<Kernel, Container>& pgn2,
                     OutputIterator oi);

/*! computes the symmetric difference between two polygons and inserts the
 * resulting polygons with holes into a container via an output iterator.
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
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename UsePolylines>
OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_2<Kernel, Container>& pgn2,
                     OutputIterator oi,
                     UsePolylines = Tag_true());

/*! computes the symmetric difference between two polygons with holes and
 * inserts the resulting polygons with holes into a container via an output
 * iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi);

/*! computes the symmetric difference between two polygons with holes and
 * inserts the resulting polygons with holes into a container via an output
 * iterator.
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
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename UsePolylines>
OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi,
                     UsePolylines = Tag_true());

/*! computes the symmetric difference between two general polygons and inserts
 * the resulting general  polygons with holes into a container via an output
 * iterator.
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
OutputIterator symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                                    const General_polygon_2<ArrTraits>& pgn2,
                                    OutputIterator oi);


/*! computes the symmetric difference between two general polygons and inserts
 * the resulting general  polygons with holes into a container via an output
 * iterator.
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
symmetric_difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
                     const General_polygon_2<ArrTraits>& pgn2,
                     OutputIterator oi);


/*! computes the symmetric difference between two general polygons and inserts
 * the resulting general  polygons with holes into a container via an output
 * iterator.
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
symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                     const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
                     OutputIterator oi);

/*! computes the symmetric difference between two general polygons and inserts
 * the resulting general  polygons with holes into a container via an output
 * iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<Polygon>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Polygon, typename OutputIterator>
OutputIterator
symmetric_difference(const General_polygon_with_holes_2<Polygon>& pgn1,
                     const General_polygon_with_holes_2<Polygon>& pgn2,
                     OutputIterator oi);

/*! Given a range of polygons (respectively general polygons) or a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * symmetric difference  of all polygons in the range and inserts the resulting
 * polygons with holes (respectively general polygons with holes) into a container via
 * an output iterator. A point is contained in the symmetric difference, if and
 * only if it is contained in an odd number of input polygons.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator, typename OutputIterator>
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
                                    OutputIterator oi);

/*! Given a range of polygons (respectively general polygons) or a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * symmetric difference  of all polygons in the range and inserts the resulting
 * polygons with holes (respectively general polygons with holes) into a container via
 * an output iterator. A point is contained in the symmetric difference, if and
 * only if it is contained in an odd number of input polygons.
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
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator, typename OutputIterator,
          typename UsePolylines>
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
                                    OutputIterator oi,
                                    UsePolylines = Tag_true());

/*! Given a range of polygons (respectively general polygons) and a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * symmetric difference between all polygons in the two ranges and inserts the
 * resulting polygons with holes (respectively general polygons with holes) into a
 * container via an output iterator. A point is contained in the symmetric
 * difference, if and only if it is contained in an odd number of input
 * polygons.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
OutputIterator symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                                    InputIterator2 begin2, InputIterator2 end2,
                                    OutputIterator oi);

/*! Given a range of polygons (respectively general polygons) and a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * symmetric difference between all polygons in the two ranges and inserts the
 * resulting polygons with holes (respectively general polygons with holes) into a
 * container via an output iterator. A point is contained in the symmetric
 * difference, if and only if it is contained in an odd number of input
 * polygons.
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
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \return the past-the-end iterator of the output container.
 */
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator,
          typename UsePolylines>
OutputIterator symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                                    InputIterator2 begin2, InputIterator2 end2,
                                    OutputIterator oi,
                                    UsePolylines = Tag_true());

//////// With Traits

/*! computes the symmetric difference between two polygons and inserts the
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
OutputIterator symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                                    const Polygon_2<Kernel, Container>& pgn2,
                                    OutputIterator oi,
                                    const GpsTraits& traits);

/*! computes the symmetric difference between two polygons and inserts the
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
symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi,
                     const GpsTraits& traits);

/*! computes the symmetric difference between two polygons and inserts the
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
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_2<Kernel, Container>& pgn2,
                     OutputIterator oi,
                     const GpsTraits& traits);


/*! computes the symmetric difference between two polygons with holes and
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
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi,
                     const GpsTraits& traits);

/*! computes the symmetric difference between two general polygons and inserts
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
OutputIterator symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                                    const General_polygon_2<ArrTraits>& pgn2,
                                    OutputIterator oi,
                                    const GpsTraits& traits);


/*! computes the symmetric difference between two general polygons and inserts
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
symmetric_difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
                     const General_polygon_2<ArrTraits>& pgn2,
                     OutputIterator oi,
                     const GpsTraits& traits);


/*! computes the symmetric difference between two general polygons and inserts
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
symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                     const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
                     OutputIterator oi,
                     const GpsTraits& traits);

/*! computes the symmetric difference between two general polygons and inserts
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
symmetric_difference(const General_polygon_with_holes_2<Polygon>& pgn1,
                     const General_polygon_with_holes_2<Polygon>& pgn2,
                     OutputIterator oi,
                     const GpsTraits& traits);

/*! Given a range of polygons (respectively general polygons) or a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * symmetric difference  of all polygons in the range and inserts the resulting
 * polygons with holes (respectively general polygons with holes) into a container via
 * an output iterator. A point is contained in the symmetric difference, if and
 * only if it is contained in an odd number of input polygons.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (respectively `General_polygon_2`) or
 *        `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator, typename OutputIterator, typename GpsTraits>
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
                                    OutputIterator oi,
                                    const GpsTraits& traits);

/*! Given a range of polygons (respectively general polygons) and a range of
 * polygons with holes (respectively general polygons with holes) computes the
 * symmetric difference between all polygons in the two ranges and inserts the
 * resulting polygons with holes (respectively general polygons with holes) into a
 * container via an output iterator. A point is contained in the symmetric
 * difference, if and only if it is contained in an odd number of input
 * polygons.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (respectively `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (respectively `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (respectively `General_polygon_with_holes_2`).
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename GpsTraits>
OutputIterator symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                                    InputIterator2 begin2, InputIterator2 end2,
                                    OutputIterator oi,
                                    const GpsTraits& traits);

/// @}

} /* namespace CGAL */
