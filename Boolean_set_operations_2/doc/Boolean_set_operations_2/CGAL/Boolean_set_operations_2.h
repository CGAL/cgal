namespace CGAL {

/*! \addtogroup boolean_complement Complement Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_complement
 *
 * There are several overloaded function templates called `complement()` that
 * compute the \e complement of a given polygon `pgn`. Depending on the type of
 * the polygon `pgn` the complement is either a single (general) polygon with
 * holes, or several (general) poylgons with holes. In the latter case the
 * `complement()` function template inserts the resulting poylgons with holes into
 * a container via an output iterator.
 *
 * \param pgn The input polygon. Its type must be convertible to one of the
 *        types `Polygon_2`, `General_polygon_2`, `Polygon_with_holes_2`, or
 *        `General_polygon_with_holes_2`.
 *
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 * \sa \link boolean_intersection `CGAL::intersection()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_difference `CGAL::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
 */

/// @{

// Traits-less

/*! Computes the complement of a polygon.
 * \param pgn the input polygon.
 * \param res the resulting complement of \p pgn.
 */
template <typename Kernel, typename Container>
void complement(const Polygon_2<Kernel, Container>& pgn,
                Polygon_with_holes_2<Kernel, Container>& res);

/*! Computes the complement of a general polygon.
 * \param pgn the input polygon.
 * \param res the complement of \p pgn.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
void complement(const General_polygon_2<ArrTraits>& pgn,
                General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& res);

/*! Computes the complement of a polygon with holes.
 * \param pgn the input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                          OutputIterator oi);

/*! Computes the complement of a general polygon with holes.
 * \param pgn the input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator complement(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn,
                          OutputIterator oi);

// With Traits

/*! Computes the complement of a polygon.
 * \param pgn the input polygon.
 * \param res the resulting complement of \p pgn
 * \param traits a traits object.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 */
template <typename Kernel, typename Container, typename GpsTraits>
void complement(const Polygon_2<Kernel, Container>& pgn,
                Polygon_with_holes_2<Kernel, Container>& res,
                const GpsTraits& traits);

/*! Computes the complement of a general polygon.
 * \param pgn the input polygon.
 * \param res the resulting complement of \p pgn
 * \param traits a traits object.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
 */
template <typename ArrTraits, typename GpsTraits>
void complement(const General_polygon_2<ArrTraits>& pgn,
                General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& res,
                const GpsTraits& traits);

/*! Computes the complement of a polygon with holes.
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

/*! Computes the complement of the general polygon with holes.
 * \param pgn the input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertible to
 *             `General_polygon_with_holes_2<<General_polygon_2<ArrTraits> >`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
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
 * computes the \e difference between two polygons `pgn1` and `pgn2` and inserts
 * the resulting polygons with holes into a container via the output iterator `oi`.
 *
 * \param oi the output iterator for the result.
 *
 * The types `Type1` and `Type2` of the parameters must be convertible to the types specified in a row in the table below, respectively.
 * The 3rd column specifies the corresponding dereference type of the output iterator.
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
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 * \sa \link boolean_intersection `CGAL::intersection()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
 */

/// @{

//////// Traits-less

/*! computes the difference of the polygons and inserts the resulting polygon
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

/*! computes the difference of the polygons and inserts the resulting polygon
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

/*! computes the difference of the polygons and inserts the resulting polygon
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

/*! computes the difference of the polygons and inserts the resulting polygon
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

/*! computes the difference of the polygons and inserts the resulting polygon
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator difference(const General_polygon_2<ArrTraits>& pgn1,
                          const General_polygon_2<ArrTraits>& pgn2,
                          OutputIterator oi);

/*! computes the difference of the polygons and inserts the resulting polygon
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
                          const General_polygon_2<ArrTraits>& pgn2,
                          OutputIterator oi);

/*! computes the difference of the polygons and inserts the resulting polygon
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator difference(const General_polygon_2<ArrTraits>& pgn1,
                          const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
                          OutputIterator oi);

/*! computes the difference of the polygons and inserts the resulting polygon
 * with holes into a container via an output iterator.
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

/*! computes the difference of the polygons and inserts the resulting polygon
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
  template <typename Kernel, typename Container, typename OutputIterator, typename GpsTraits>
OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                          const Polygon_2<Kernel, Container>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of the polygons and inserts the resulting polygon
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
template <typename Kernel, typename Container, typename OutputIterator, typename GpsTraits>
OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                          const Polygon_with_holes_2<Kernel,Container>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of the polygons and inserts the resulting polygon
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
template <typename Kernel, typename Container, typename OutputIterator, typename GpsTraits>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                          const Polygon_2<Kernel, Container>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of the polygons and inserts the resulting polygon
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
template <typename Kernel, typename Container, typename OutputIterator, typename GpsTraits>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                          const Polygon_with_holes_2<Kernel, Container>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of the polygons and inserts the resulting polygon
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator difference(const General_polygon_2<ArrTraits>& pgn1,
                          const General_polygon_2<ArrTraits>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of the polygons and inserts the resulting polygon
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
                          const General_polygon_2<ArrTraits>& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of the polygons and inserts the resulting polygon
 * with holes into a container via an output iterator.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param oi the output iterator for the result.
 *           Its dereference type must be convertibe to
 *             `General_polygon_with_holes_2<General_polygon_2<ArrTraits> >`.
 * \param traits a traits object.
 * \return the past-the-end iterator of the output container.
 * \pre GpsTraits must be a model of `GeneralPolygonSetTraits_2`.
 */
template <typename ArrTraits, typename OutputIterator, typename GpsTraits>
OutputIterator difference(const General_polygon_2<ArrTraits>& pgn1,
                          const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! computes the difference of the polygons and inserts the resulting polygon
 * with holes into a container via an output iterator.
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
 * Each one of these functions computes if the interior of two given
 * polygons `pgn1` and `pgn2` intersect.
 *
 * The signature of the function is:
 *   - `bool do_intersect(const Type1& pgn1, const Type2& pgn2);`
 *
 * \cgalHeading{Parameters}
 *
 * The types of the parameters of the \link ref_bso_do_intersect
 * `do_intersect()` \endlink function are any of the following combinations.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th> Type1 </th><th>Type2</th></tr>
 * <tr><td valign="center">Polygon_2</td><td valign="center">Polygon_2</td></tr>
 * <tr><td valign="center">Polygon_2</td><td valign="center">Polygon_with_holes_2</td></tr>
 * <tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_2</td></tr>
 * <tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_with_holes_2</td></tr>
 * <tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_2</td></tr>
 * <tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_with_holes_2</td></tr>
 * <tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr>
 * <tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr>
 * </table>
 * </div>
 *
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
 *
 * \sa \link boolean_intersection `CGAL::intersection()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_difference `CGAL::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
 */

/// @{

/*! returns `true` if the polygons `pgn1` and `pgn2` intersect in their interior.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2);

/*! returns `true` if the polygons `pgn1` and `pgn2` intersect in their interior.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2);

/*! returns `true` if the polygons `pgn1` and `pgn2` intersect in their interior.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2);

/*! returns `true` if the polygons `pgn1` and `pgn2` intersect in their interior.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2);

/*! returns `true` if the general polygons `pgn1` and `pgn2` intersect in their
 * interior.
 */
template <typename ArrTraits>
bool do_intersect(const General_polygon_2<ArrTraits>& pgn1,
                  const General_polygon_2<ArrTraits>& pgn2);

/*! returns `true` if the general polygons `pgn1` and `pgn2` intersect in their
 * interior.
 */
template <typename ArrTraits>
bool do_intersect(const General_polygon_2<ArrTraits>& pgn1,
                  const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2);

/*! returns `true` if the general polygons `pgn1` and `pgn2` intersect in their
 * interior.
 */
template <typename ArrTraits>
bool do_intersect(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
                  const General_polygon_2<ArrTraits>& pgn2);


/*! returns `true` if the general polygons `pgn1` and `pgn2` intersect in their
 * interior.
 */
template <typename Polygon>
bool do_intersect(const General_polygon_with_holes_2<Polygon>& pgn1,
                  const General_polygon_with_holes_2<Polygon>& pgn2);

/*! returns `true`, if the set of general polygons (or general polygons with
 * holes) in the given range intersect in their interior, and `false`
 * otherwise. (The value type of the input iterator is used to distinguish
 * between the two).
 */
template <typename InputIterator>
bool do_intersect(InputIterator begin, InputIterator end);

/*! returns `true`, if the set of general polygons and general polygons with
 * holes in the given two ranges respectively intersect in their interior, and
 * `false` otherwise.
 */
template <typename InputIterator1, typename InputIterator2>
bool do_intersect(InputIterator1 pgn_begin1,
                  InputIterator1 pgn_end1,
                  InputIterator2 pgn_begin2,
                  InputIterator2 pgn_end2);
/// @}
} /* namespace CGAL */

namespace CGAL {

/*! \addtogroup boolean_intersection Intersection Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_intersection
 *
 * Each one of these functions computes the intersection of two given
 * polygons `pgn1` and `pgn2`, inserts the resulting polygons with
 * holes into an output container through a given output iterator
 * `oi`, and returns the output iterator. The value type of the
 * `OutputIterator` is either `Polygon_with_holes_2` or
 * `General_polygon_with_holes_2`.
 *
 * The signature of the function is:
 *   - `%OutputIterator %intersection(const Type1& pgn1, const Type2& pgn2,
 *                                    %OutputIterator oi);`
 *
 * \cgalHeading{Parameters}
 *
 * The types of the parameters of the `intersection()` function are any of the
 * following combinations.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th> Type1</th><th> Type2</th></tr>
 * <tr><td valign="center">Polygon_2</td><td valign="center">Polygon_2</td></tr>
 * <tr><td valign="center">Polygon_2</td><td valign="center">Polygon_with_holes_2</td></tr>
 * <tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_2</td></tr>
 * <tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_with_holes_2</td></tr>
 * <tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_2</td></tr>
 * <tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_with_holes_2</td></tr>
 * <tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr>
 * <tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr>
 * </table>
 * </div>
 *
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
 *
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_difference `CGAL::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
 *
 */

/// @{

/*! computes the intersection of the polygons `pgn1` and `pgn2` into the output iterator `oi`.
 *  The value type of `oi` is `Polygon_with_holes_2`.
 */
OutputIterator intersection(const Type1& pgn1, const Type2& pgn2,
                            OutputIterator oi);

/*! computes the intersection of the polygons `pgn1` and `pgn2` into the output iterator `oi`.
 * The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of the polygons `pgn1` and `pgn2` into the output iterator `oi`.
 * The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of the polygons `pgn1` and `pgn2` into the output iterator `oi`.
 * The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of the polygons `pgn1` and `pgn2` into the output iterator `oi`.
 * The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of the general polygons `pgn1` and `pgn2` into the output iterator `oi`.
 * The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator intersection(const General_polygon_2<ArrTraits>& pgn1,
                            const General_polygon_2<ArrTraits>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of the general polygons `pgn1` and `pgn2` into the output iterator `oi`.
 * The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator intersection(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
                            const General_polygon_2<ArrTraits>& pgn2,
                            OutputIterator oi);

/*! computes the intersection of the general polygons `pgn1` and `pgn2` into the output iterator `oi`.
 * The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator intersection(const General_polygon_2<ArrTraits>& pgn1,
                            const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
                            OutputIterator oi);

/*! computes the intersection of the general polygons `pgn1` and `pgn2` into the output iterator `oi`.
 * The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename Polygon, typename OutputIterator>
OutputIterator intersection(const General_polygon_with_holes_2<Polygon>& pgn1,
                            const General_polygon_with_holes_2<Polygon>& pgn2,
                            OutputIterator oi);


/*! computes the intersection of the general polygons (or general polygons with
 * holes) in the given range. (The value type of the input iterator is used to
 * distinguish between the two.) The result, represented by a set of general
 * polygon with holes, is written into the output iterator `oi`.  The output
 * iterator is returned. The value type of the `OutputIterator` is
 * `Traits::Polygon_with_holes_2`.
 */
template <typename InputIterator, typename OutputIterator>
OutputIterator intersection(InputIterator begin, InputIterator end,
                            OutputIterator oi);

/*! computes the intersection of the general polygons and general polygons with
 * holes in the given two ranges. The result, represented by a set of general
 * polygon with holes, is written into the output iterator `oi`.  The output
 * iterator is returned. The value type of the `OutputIterator` is
 * `Traits::Polygon_with_holes_2`.
 */
template <typename InputIterator1, typename InputIterator2,
typename OutputIterator>
OutputIterator intersection(InputIterator1 pgn_begin1,
                            InputIterator1 pgn_end1,
                            InputIterator2 pgn_begin2,
                            InputIterator2 pgn_end2,
                            OutputIterator oi);

/// @}

} /* namespace CGAL */

namespace CGAL {

/*! \addtogroup boolean_join Union Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_union
 *
 * Each one of these functions computes the union of two given polygons `pgn1` and
 * `pgn2`. If the two given polygons overlap, it returns `true`, and places the
 * resulting polygon in `p`. Otherwise, it returns `false`.
 *
 * The signature of the function is:
 *   - `bool join(const Type1& pgn1, const Type2& pgn2, General_polygon_with_holes_2& res);`
 *
 * \cgalHeading{Parameters}
 *
 * The types of the parameters of the `join()` function are any of the following combinations.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th> Type1</th><th> Type2</th></tr>
 * <tr><td valign="center">Polygon_2</td><td valign="center">Polygon_2</td></tr>
 * <tr><td valign="center">Polygon_2</td><td valign="center">polygon_with_holes_2</td></tr>
 * <tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_2</td></tr>
 * <tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_with_holes_2</td></tr>
 * <tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_2</td></tr>
 * <tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_with_holes_2</td></tr>
 * <tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr>
 * <tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr>
 * </table>
 * </div>
 *
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
 *
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 * \sa \link boolean_intersection `CGAL::intersection()` \endlink
 * \sa \link boolean_difference `CGAL::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
 */

/// @{

/*! computes the union of the polygons `pgn1` and `pgn2` into the polygon with holes `res`.
 *  Returns `true` if the two given polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_2<Kernel, Container>& pgn1,
          const Polygon_2<Kernel, Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! computes the union of the polygons `pgn1` and `pgn2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_2<Kernel, Container>& pgn1,
          const Polygon_with_holes_2<Kernel,Container>& pgn2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! computes the union of the polygons `pgn1` and `pgn2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_with_holes_2<Kernel, Container>& pgn2,
          const Polygon_2<Kernel, Container>& pgn1,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! computes the union of the polygons `pgn1` and `pgn2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_with_holes_2<Kernel, Container>& pgn2,
          const Polygon_with_holes_2<Kernel, Container>& pgn1,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! computes the union of the general polygons `pgn1` and `pgn2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename ArrTraits>
bool join(const General_polygon_2<ArrTraits>& pgn1,
          const General_polygon_2<ArrTraits>& pgn2,
          General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res);


/*! computes the union of the polygons `pgn1` and `pgn2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename ArrTraits>
bool join(const General_polygon_2<ArrTraits>& pgn1,
          const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
          General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res);


/*! computes the union of the general polygons `pgn1` and `pgn2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename ArrTraits>
bool join(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
          const General_polygon_2<ArrTraits>& pgn1,
          General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res);


/*! computes the union of the general polygons `pgn1` and `pgn2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename Polygon>
bool join(const General_polygon_with_holes_2<Polygon>& pgn1,
          const General_polygon_with_holes_2<Polygon>& pgn2,
          Traits::Polygon_with_holes_2& res);

/*! computes the union of the general polygons (or general polygons with holes)
 * in the given range. (The value type of the input iterator is used to
 * distinguish between the two.) The result, represented by a set of general
 * polygon with holes, is written into the output iterator `oi`.  The output
 * iterator is returned. The value type of the `OutputIterator` is
 * `Traits::Polygon_with_holes_2`.
 */
template <typename InputIterator, typename OutputIterator>
OutputIterator join(InputIterator begin, InputIterator end,
                    OutputIterator oi);

/*! computes the union of the general polygons and general polygons with holes
 * in the given two ranges. The result, represented by a set of general polygon
 * with holes, is written into the output iterator `oi`.  The output iterator is
 * returned. The value type of the `OutputIterator` is
 * `Traits::Polygon_with_holes_2`.
 */
template <typename InputIterator1, typename InputIterator2, typename OutputIterator>
OutputIterator join(InputIterator1 pgn_begin1, InputIterator1 pgn_end1,
                    InputIterator2 pgn_begin2, InputIterator2 pgn_end2,
                    OutputIterator oi);

/// @}
} /* namespace CGAL */

namespace CGAL {
/*! \addtogroup boolean_oriented_side Oriented Side Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_oriented_side
 *
 * There are several overloaded function templatess called `Oriented_side()`
 * that computes the relative position of either (i) two polygons or (ii) a
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
 *   - `Oriented_side oriented_side(const Type1& pgn1, const Type2& pgn2, const GpsTraits& traits);`
 *
 * The types `Type1` and `Type2` of the parameters must be convertible to the types specified in a row in the following table, respectively.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th>Type1</th>                                       <th>Type2</th></tr>
 * <tr><td valign="center">Polygon_2</td>                   <td valign="center">Polygon_2</td></tr>
 * <tr><td valign="center">Polygon_2</td>                   <td valign="center">Polygon_with_holes_2</td></tr>
 * <tr><td valign="center">Polygon_with_holes_2</td>        <td valign="center">Polygon_2</td></tr>
 * <tr><td valign="center">Polygon_with_holes_2</td>        <td valign="center">Polygon_with_holes_2</td></tr>
 * <tr><td valign="center">General_polygon_2</td>           <td valign="center">General_polygon_2</td></tr>
 * <tr><td valign="center">General_polygon_2</td>           <td valign="center">General_polygon_with_holes_2</td></tr>
 * <tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr>
 * <tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr>
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
 *   - `Oriented_side oriented_side(const Point_2& p, const Type& pgn, const GpsTraits& traits);`
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
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
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
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2);

/*! computes the relative position of two polygons.
 * \param pgn1 1st the input polygon.
 * \param pgn2 the 2nd input polygon.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
Oriented_side oriented_side(const General_polygon_2<ArrTraits>& pgn1,
                            const General_polygon_2<ArrTraits>& pgn2);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
Oriented_side oriented_side(const General_polygon_2<ArrTraits>& pgn1,
                            const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 */
template <typename ArrTraits>
Oriented_side oriented_side(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
                            const General_polygon_2<ArrTraits>& pgn2);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 */
template <typename Polygon>
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon>& pgn1,
                            const General_polygon_with_holes_2<Polygon>& pgn2);

// Polygon--Polygon With Traits

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
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
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
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
 *
 */
template <typename ArrTraits, typename GpsTraits>
Oriented_side oriented_side(const General_polygon_2<ArrTraits>& pgn1,
                            const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
                            const GpsTraits& traits);

/*! computes the relative position of two polygons.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \param traits a traits object.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
 *
 */
template <typename ArrTraits, typename GpsTraits>
Oriented_side oriented_side(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
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
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
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

// Point--Polygon With Traits

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
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
 *
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
 * Each one of these functions computes the symmetric difference between
 * two given polygons `pgn1` and `pgn2`, and inserts the resulting
 * polygons with holes into an output container through the output
 * iterator `oi`. The value type of the `OutputIterator` is either
 * `Polygon_with_holes_2` or
 * `General_polygon_with_holes_2`.
 *
 * The signature of the function is:
 *   - `%OutputIterator symmetric_difference(const Type1& pgn1, const Type2& pgn2, %OutputIterator oi);`
 *
 * \cgalHeading{Parameters}
 *
 * The types of the parameters of the `symmetric_difference()` function are any of the following combinations.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th> Arg 1 type</th><th>Arg 2 type</th></tr>
 * <tr><td valign="center">Polygon_2</td><td valign="center">Polygon_2</td></tr>
 * <tr><td valign="center">Polygon_2</td><td valign="center">Polygon_with_holes_2</td></tr>
 * <tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_2</td></tr>
 * <tr><td valign="center">Polygon_with_holes_2</td><td valign="center">Polygon_with_holes_2</td></tr>
 * <tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_2</td></tr>
 * <tr><td valign="center">General_polygon_2</td><td valign="center">General_polygon_with_holes_2</td></tr>
 * <tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_2</td></tr>
 * <tr><td valign="center">General_polygon_with_holes_2</td><td valign="center">General_polygon_with_holes_2</td></tr>
 * </table>
 * </div>
 *
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
 *
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 * \sa \link boolean_intersection `CGAL::intersection()` \endlink
 * \sa \link boolean_join `CGAL::join()` \endlink
 * \sa \link boolean_difference `CGAL::difference()` \endlink
 *
 */

/// @{

OutputIterator symmetric_difference(const Type1& pgn1, const Type2& pgn2,
                                    OutputIterator oi);


template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                                    const Polygon_2<Kernel, Container>& pgn2,
                                    OutputIterator oi);


template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator
symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi);

template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_2<Kernel, Container>& pgn2,
                     OutputIterator oi);


template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi);


template <typename ArrTraits, typename OutputIterator>
OutputIterator symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                                    const General_polygon_2<ArrTraits>& pgn2,
                                    OutputIterator oi);


template <typename ArrTraits, typename OutputIterator>
OutputIterator
symmetric_difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
                     const General_polygon_2<ArrTraits>& pgn2,
                     OutputIterator oi);


template <typename ArrTraits, typename OutputIterator>
OutputIterator
symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                     const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
                     OutputIterator oi);


template <typename Polygon, typename OutputIterator>
OutputIterator
symmetric_difference(const General_polygon_with_holes_2<Polygon>& pgn1,
                     const General_polygon_with_holes_2<Polygon>& pgn2,
                     OutputIterator oi);

/*! computes the symmetric difference of the general polygons (or general
 * polygons with holes) in the given range. A point is contained in the
 * symmetric difference, if and only if it is contained in an odd number of
 * input polygons. (The value type of the input iterator is used to distinguish
 * between the two.) The result, represented by a set of general polygon with
 * holes, is inserted into an output container through a given output iterator
 * `oi`. The output iterator is returned. The value type of the `OutputIterator`
 * is `Traits::Polygon_with_holes_2`.
 */
template <typename InputIterator, typename OutputIterator>
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
                                    OutputIterator oi);

/*! computes the symmetric difference of the general polygons and general
 * polygons with holes in the given two ranges. A point is contained in the
 * symmetric difference, if and only if it is contained in an odd number of
 * input polygons. The result, represented by a set of general polygon with
 * holes, is inserted into an output container through a given output iterator
 * `oi`. The output iterator is returned. The value type of the `OutputIterator`
 * is `Traits::Polygon_with_holes_2`.
 */
template <typename InputIterator1, typename InputIterator2, typename OutputIterator>
OutputIterator symmetric_difference(InputIterator1 pgn_begin1,
                                    InputIterator1 pgn_end1,
                                    InputIterator2 pgn_begin2,
                                    InputIterator2 pgn_end2,
                                    OutputIterator oi);
/// @}

} /* namespace CGAL */
