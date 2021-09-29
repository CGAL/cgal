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
 * A function template in this group has one of the two following signatures:
 *   - `void complement(const Type1& pgn, Type2& res);`
 *   - `void complement(const Type1& pgn, Type2& res, const GpsTraits& traits);`
 *
 * \cgalHeading{Parameters}
 *
 * The types `Type` and `Type2` of the parameters must be convertible to the
 * types specified in a row in the table below, respectively.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>`Type1`</th>                                       <th>`Type2`</th></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
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

/*! computes the complement of a general polygon.
 * \param pgn the input polygon.
 * \param res the complement of \p pgn.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
void
complement(const General_polygon_2<ArrTraits>& pgn,
           General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& res);

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

/*! computes the complement of a general polygon with holes.
 * \param pgn the input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
complement(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn,
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
           General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& res,
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
 *             `General_polygon_with_holes_2<<General_polygon_2<ArrTraits> >`.
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
 *   - `OutputIterator difference(const Type1& pgn1, const Type2& pgn2,
 *                                OutputIterator oi);`
 *   - `OutputIterator difference(const Type1& pgn1, const Type2& pgn2,
 *                                OutputIterator oi, const GpsTraits& traits);`
 *
 * \param oi the output iterator for the result.
 *
 * The types `Type1` and `Type2` of the parameters must be convertible to the
 * types specified in a row in the table below, respectively.  The 3rd column
 * specifies the corresponding dereference type of the output iterator.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>`Type1`</th>                                       <th>`Type2`</th>                                       <th>%Dereference Type of `oi`</th></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_2`</td>                   <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td>        <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_2`</td>                   <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td>        <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_2`</td>           <td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td><td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_2`</td>           <td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td><td>`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
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

/*! computes the difference of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
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
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
           const General_polygon_2<ArrTraits>& pgn2,
           OutputIterator oi);

/*! computes the difference of two polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
difference(const General_polygon_2<ArrTraits>& pgn1,
           const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
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
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
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
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator
difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
           const General_polygon_2<ArrTraits>& pgn2,
           OutputIterator oi,
           const GpsTraits& traits);

/*! computes the difference of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator
difference(const General_polygon_2<ArrTraits>& pgn1,
           const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
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

namespace CGAL {

/*! \addtogroup boolean_do_intersect Intersection Testing Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_do_intersect
 *
 * There are several overloaded function templates called `do_intersect()`
 * that determine whether the interior of two or more polygons intersect.
 *
 * A function template in this group that accepts two input polygons has one of
 * the two following signatures:
 *   - `bool do_intersect(const Type1& pgn1, const Type2& pgn2);`
 *   - `bool do_intersect(const Type1& pgn1, const Type2& pgn2,
                          const GpsTraits& traits);`
 *
 * \cgalHeading{Parameters}
 *
 * The types `Type1` and `Type2` of the parameters must be convertible to the
 * types specified in a row in the table below, respectively.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>`Type1`</th>                                       <th>`Type2`</th></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_2`</td></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2);

/*! determines whether two polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2);

/*! determines whether two polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2);

/*! determines whether two polygons with holes intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2);

/*! determines whether two general polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
bool
do_intersect(const General_polygon_2<ArrTraits>& pgn1,
             const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2);

/*! determines whether two general polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
bool do_intersect(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
                  const General_polygon_2<ArrTraits>& pgn2);

/*! determines whether two general polygons with holes intersect in their
 * interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 */
template <typename Polygon>
bool do_intersect(const General_polygon_with_holes_2<Polygon>& pgn1,
                  const General_polygon_with_holes_2<Polygon>& pgn2);

/*! Given a range of polygons or a range of polygons with holes (resp. a range
 * of general polygons or a range of general polygons with holes) determines
 * whether the open polygons (resp. general polygons) in the range have a common
 * point.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or `Polygon_with_holes_2`
 *        (resp. `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or `Polygon_with_holes_2`
 *        (resp. `General_polygon_with_holes_2`).
 * \return `true` if the pairwise intersections of all open polygons or polygons
 *         with holes (resp. general polygons or general polygons with holes) in
 *         the range [*begin,*end) overlap, and `false` otherwise.
 */
template <typename InputIterator>
bool do_intersect(InputIterator begin, InputIterator end);

/*! Given a range of polygons (resp. general polygons) and a range of polygons
 * with holes (resp. general polygons with holes) determines whether the open
 * polygons (resp. general polygons) in the two ranges have a common point.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (resp. `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (resp. `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \return `true` if the pairwise intersections of all open polygons (resp.
 *        general polygons) and polygons with holes (resp. general polygons with
 *        holes) in the ranges [*begin1,*end1) and [*begin2,*end2),
 *        respectively, overlap, and `false` otherwise.
 */
template <typename InputIterator1, typename InputIterator2>
bool do_intersect(InputIterator1 begin1,
                  InputIterator1 end1,
                  InputIterator2 begin2,
                  InputIterator2 end2);

//////// With Traits

/*! determines whether two polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename GpsTraits>
bool
do_intersect(const General_polygon_2<ArrTraits>& pgn1,
             const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
             const GpsTraits& traits);

/*! determines whether two general polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename GpsTraits>
bool
do_intersect(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
             const General_polygon_2<ArrTraits>& pgn2,
             const GpsTraits& traits);

/*! determines whether two general polygons with holes intersect in their
 * interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Polygon, typename GpsTraits>
bool do_intersect(const General_polygon_with_holes_2<Polygon>& pgn1,
                  const General_polygon_with_holes_2<Polygon>& pgn2,
                  const GpsTraits& traits);

/*! Given a range of polygons or a range of polygons with holes (resp. a range
 * of general polygons or a range of general polygons with holes) determines
 * whether the open polygons (resp. general polygons) in the range have a common
 * point.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or `Polygon_with_holes_2`
 *        (resp. `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or `Polygon_with_holes_2`
 *        (resp. `General_polygon_with_holes_2`).
 * \param traits a traits object.
 * \return `true` if the pairwise intersections of all open polygons or polygons
 *         with holes (resp. general polygons or general polygons with holes) in
 *         the range [*begin,*end) overlap, and `false` otherwise.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator, typename GpsTraits>
bool do_intersect(InputIterator begin, InputIterator end,
                  const GpsTraits& traits);

/*! Given a range of polygons (resp. general polygons) and a range of polygons
 * with holes (resp. general polygons with holes) determines whether the open
 * polygons (resp. general polygons) in the two ranges have a common point.
 * \param begin1 the first iterator of the 1st input range. Its value type is
 *        `Polygon_2` (resp. `General_polygon_2`).
 * \param end1 the past-the-end iterator of the 1st input range. Its value
 *        type is `Polygon_2` (resp. `General_polygon_2`).
 * \param begin2 the first iterator of the 2nd input range. Its value type
 *        is `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param end2 the past-the-end iterator of the 2nd input range. Its value
 *        type is `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param traits a traits object.
 * \return `true` if the pairwise intersections of all open polygons (resp.
 *        general polygons) and polygons with holes (resp. general polygons with
 *        holes) in the ranges [*begin1,*end1) and [*begin2,*end2),
 *        respectively, overlap, and `false` otherwise.
 */
template <typename InputIterator1, typename InputIterator2, typename GpsTraits>
bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                  InputIterator2 begin2, InputIterator2 end2,
                  const GpsTraits& traits);

/// @}
} /* namespace CGAL */

namespace CGAL {

/*! \addtogroup boolean_intersection Intersection Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_intersection
 *
 * There are several overloaded function templates called `intersection()` that
 * compute the \e intersection of two or more polygons and insert the resulting
 * polygons with holes into a container via an output iterator.
 *
 * A function template in this group that accepts two input polygons has one of
 * the two following signatures:
 *   - `OutputIterator %intersection(const Type1& pgn1, const Type2& pgn2,
 *                                   OutputIterator oi);`
 *   - `OutputIterator %intersection(const Type1& pgn1, const Type2& pgn2,
 *                                   OutputIterator oi, const GpsTraits& traits);`
 *
 * \param oi the output iterator for the result.
 *
 * The types `Type1` and `Type2` of the parameters must be convertible to the
 * types specified in a row in the table below, respectively.  The 3rd column
 * specifies the corresponding dereference type of the output iterator.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>`Type1`</th>                                       <th>`Type2`</th>                                       <th>%Dereference Type of `oi`</th></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_2`</td>                   <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td>        <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_2`</td>                   <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td>        <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_2`</td>           <td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td><td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_2`</td>           <td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td><td>`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
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

/*! computes the intersection of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
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
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
intersection(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
             const General_polygon_2<ArrTraits>& pgn2,
             OutputIterator oi);

/*! computes the intersection of two general polygons and inserts the resulting
 * general polygons with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
intersection(const General_polygon_2<ArrTraits>& pgn1,
             const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
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


/*! Given a range of polygons (resp. general polygons) or a range of general
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

/*! Given a range of polygons (resp. general polygons) and a range of general
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
OutputIterator intersection(InputIterator1 begin1,
                            InputIterator1 end1,
                            InputIterator2 begin2,
                            InputIterator2 end2,
                            OutputIterator oi);

/// @}

} /* namespace CGAL */

namespace CGAL {

/*! \addtogroup boolean_join Union Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_union
 *
 * There are several overloaded function templates called `join()` that
 * compute the \e union of two polygons.
 *
 * A function template in this group has one of the two following signatures:
 *   - `bool join(const Type1& pgn1, const Type2& pgn2, Type3& res);`
 *   - `bool join(const Type1& pgn1, const Type2& pgn2, Type3& res,
 *                const GpsTraits& traits);`
 *
 * \cgalHeading{Parameters}
 *
 * The types `Type1`, `Type2`, and `Type3, of the parameters must be convertible
 * to the types specified in a row in the table below, respectively.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>`Type1`</th>                                       <th>`Type2`</th>                                       <th>`Type3`</th></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_2`</td>                   <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td>        <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_2`</td>                   <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td>        <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_2`</td>           <td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td><td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_2`</td>           <td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td><td>`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
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
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! computes the union of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_2<Kernel, Container>& pgn1,
          const Polygon_with_holes_2<Kernel,Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! computes the union of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_with_holes_2<Kernel, Container>& pgn2,
          const Polygon_2<Kernel, Container>& pgn1,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! computes the union of two polygons with holes.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param res the resulting union of \p pgn1 and \p pgn2.
 * \return `true` if the two input polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_with_holes_2<Kernel, Container>& pgn2,
          const Polygon_with_holes_2<Kernel, Container>& pgn1,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


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
          General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res);


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
     const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res);


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
join(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
     const General_polygon_2<ArrTraits>& pgn1,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res);


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

/*! Given a range of polygons (resp. general polygons) or a range of general
 * polygons with holes (resp. general polygons with holes) computes the
 * union of all polygons in the range and inserts the resulting polygons
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
OutputIterator join(InputIterator begin, InputIterator end,
                    OutputIterator oi);

/*! Given a range of polygons (resp. general polygons) and a range of general
 * polygons with holes (resp. general polygons with holes) computes the
 * union of all polygons in the two ranges and inserts the resulting
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
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res,
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
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res,
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
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res,
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
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res,
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
          General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res,
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
     const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res,
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
bool join(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
          const General_polygon_2<ArrTraits>& pgn1,
          General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res,
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

/*! Given a range of polygons (resp. general polygons) or a range of general
 * polygons with holes (resp. general polygons with holes) computes the
 * union of all polygons in the range and inserts the resulting polygons
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
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator, typename OutputIterator, typename GpsTraits>
OutputIterator join(InputIterator begin, InputIterator end,
                    OutputIterator oi,
                    const GpsTraits& traits);

/*! Given a range of polygons (resp. general polygons) and a range of general
 * polygons with holes (resp. general polygons with holes) computes the
 * union of all polygons in the two ranges and inserts the resulting
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
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator1, typename InputIterator2, typename OutputIterator, typename GpsTraits>
OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                    InputIterator2 begin2, InputIterator2 end2,
                    OutputIterator oi,
                    const GpsTraits& traits);

/// @}
} /* namespace CGAL */

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
 *   - `Oriented_side oriented_side(const Type1& pgn1, const Type2& pgn2);`
 *   - `Oriented_side oriented_side(const Type1& pgn1, const Type2& pgn2,
 *                                  const GpsTraits& traits);`
 *
 * \cgalHeading{Parameters}
 *
 * The types `Type1` and `Type2` of the parameters must be convertible to the types specified in a row in the following table, respectively.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>Type1</th>                                         <th>Type2</th></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_2`</td></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
 *
 * \cgalHeading{Oriented Side of a Point and a Polygon}
 *
 * The functions in the second group accept a point `p` and a polygon `pgn`.
 * Each function in this group returns `ON_POSITIVE_SIDE` if the point `p`
 * is in the interior of `pgn`, `ON_NEGATIVE_SIDE` if `p` is in the exterior
 * of `pgn`, and `ON_ORIENTED_BOUNDARY` if `p` is on the boundary of `pgn`.
 *
 * A function in this subgroup has one of the two following signatures:
 *   - `Oriented_side oriented_side(const Point_2& p, const Type& pgn);`
 *   - `Oriented_side oriented_side(const Point_2& p, const Type& pgn,
 *                                  const GpsTraits& traits);`
 *
 * \cgalHeading{Parameters}
 *
 * `Type` must be convertible to one of
 *   `Polygon_2`,
 *   `Polygon_with_holes_2`,
 *   `General_polygon_2`, or
 *   `General_polygon_with_holes_2`.
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
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 */
template <typename Kernel, typename Container>
Oriented_side
oriented_side(const Polygon_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 */
template <typename Kernel, typename Container>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 */
template <typename Kernel, typename Container>
Oriented_side
oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2);

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
              const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
Oriented_side
oriented_side(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
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
              const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
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
oriented_side(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
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
 * the two following signatures:
 *   - `OutputIterator symmetric_difference(const Type1& pgn1, const Type2& pgn2,
 *                                          OutputIterator oi);`
 *   - `OutputIterator symmetric_difference(const Type1& pgn1, const Type2& pgn2,
 *                                          OutputIterator oi, const GpsTraits& traits);`
 *
 * \param oi the output iterator for the result.
 *
 * The types `Type1` and `Type2` of the parameters must be convertible to the
 * types specified in a row in the table below, respectively.  The 3rd column
 * specifies the corresponding dereference type of the output iterator.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>`Type1`</th>                                       <th>`Type2`</th>                                       <th>%Dereference Type of `oi`</th></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_2`</td>                   <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_2`</td>                   <td valign="center">`Polygon_with_holes_2`</td>        <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_2`</td>                   <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`Polygon_with_holes_2`</td>        <td valign="center">`Polygon_with_holes_2`</td>        <td>`Polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_2`</td>           <td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_2`</td>           <td valign="center">`General_polygon_with_holes_2`</td><td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_2`</td>           <td>`General_polygon_with_holes_2`</td></tr>
 * <tr><td valign="center">`General_polygon_with_holes_2`</td><td valign="center">`General_polygon_with_holes_2`</td><td>`General_polygon_with_holes_2`</td></tr>
 * </table>
 * </div>
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

/*! computes the symmetric difference between two general polygons and inserts
 * the resulting general  polygons with holes into a container via an output
 * iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
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
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
symmetric_difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
                     const General_polygon_2<ArrTraits>& pgn2,
                     OutputIterator oi);


/*! computes the symmetric difference between two general polygons and inserts
 * the resulting general  polygons with holes into a container via an output
 * iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator
symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                     const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
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

/*! Given a range of polygons (resp. general polygons) or a range of general
 * polygons with holes (resp. general polygons with holes) computes the
 * symmetric difference  of all polygons in the range and inserts the resulting
 * polygons with holes (resp. general polygons with holes) into a container via
 * an output iterator. A point is contained in the symmetric difference, if and
 * only if it is contained in an odd number of input polygons.
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
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
                                    OutputIterator oi);

/*! Given a range of polygons (resp. general polygons) and a range of general
 * polygons with holes (resp. general polygons with holes) computes the
 * intersection of all polygons in the two ranges and inserts the resulting
 * polygons with holes (resp. general polygons with holes) into a container via
 * an output iterator. A point is contained in the symmetric difference, if and
 * only if it is contained in an odd number of input polygons.
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
OutputIterator symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                                    InputIterator2 begin2, InputIterator2 end2,
                                    OutputIterator oi);

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
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
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
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator
symmetric_difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
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
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept
 *      `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator
symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                     const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
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

/*! Given a range of polygons (resp. general polygons) or a range of general
 * polygons with holes (resp. general polygons with holes) computes the
 * symmetric difference  of all polygons in the range and inserts the resulting
 * polygons with holes (resp. general polygons with holes) into a container via
 * an output iterator. A point is contained in the symmetric difference, if and
 * only if it is contained in an odd number of input polygons.
 * \param begin the first iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or `Polygon_with_holes_2`
 *        (resp. `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or `Polygon_with_holes_2`
 *        (resp. `General_polygon_with_holes_2`).
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2` (resp. `General_polygons_with_holes_2`).
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename InputIterator, typename OutputIterator, typename GpsTraits>
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
                                    OutputIterator oi,
                                    const GpsTraits& traits);

/*! Given a range of polygons (resp. general polygons) and a range of general
 * polygons with holes (resp. general polygons with holes) computes the
 * intersection of all polygons in the two ranges and inserts the resulting
 * polygons with holes (resp. general polygons with holes) into a container via
 * an output iterator. A point is contained in the symmetric difference, if and
 * only if it is contained in an odd number of input polygons.
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
OutputIterator symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                                    InputIterator2 begin2, InputIterator2 end2,
                                    OutputIterator oi,
                                    const GpsTraits& traits);

/// @}

} /* namespace CGAL */
