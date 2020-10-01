namespace CGAL {

/*! \addtogroup boolean_complement Complement Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_complement
 *
 * There are several overloaded functions called `complement` that computes
 * the complement of a given polygon `pgn` . Depending on the type of polygon
 * `pgn` the complement is either a single (general) polygon with holes, or
 * several (general) poylgons with holes. In the latter case the `complement
 * function` writes them into an output iterator `oi`.
 *
 * \param pgn The input polygon for the `complement` function. Its type must
 *        be convertible to one of the types `Polygon_2`, `General_polygon_2`,
 *        `Polygon_with_holes_2`, or `General_polygon_with_holes_2`.
 *
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
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
 * \param pgn the input polygon
 * \param res the complement of \p pgn.
 */
template <typename Kernel, typename Container>
void complement(const Polygon_2<Kernel, Container>& pgn,
                Polygon_with_holes_2<Kernel, Container>& res);

/*! Computes the complement of a general polygon.
 * \param pgn the input polygon
 * \param res the complement of \p pgn
 */
template <typename ArrTraits>
void complement(const General_polygon_2<ArrTraits>& pgn,
                General_polygon_with_holes_2<ArrTraits>& res);

/*! Computes the complement of a polygon with holes.
 * \param pgn the input polygon
 * \param oi the output iterator for the result.
 *           Its dereference type is `Polygon_with_holes_2<Kernel, Container>`.
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                          OutputIterator oi);

/*! Computes the complement of a general polygon with holes.
 * \param pgn the input polygon
 * \param oi the output iterator for the result.
 *           Its dereference type is
 *             `General_polygon_with_holes_2<<General_polygon_2<ArrTraits> >`.
 * \return the past-the-end iterator of the output container.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator complement(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn,
                          OutputIterator oi);

// With Traits

/*! Computes the complement of a polygon.
 * \param pgn the input polygon
 * \param res the complement of \p pgn
 * \param traits a model of `GeneralPolygonSetTraits_2`
 */
template <typename Kernel, typename Container, typename GpsTraits>
void complement(const Polygon_2<Kernel, Container>& pgn,
                Polygon_with_holes_2<Kernel, Container>& res,
                const GpsTraits& traits);

/*! Computes the complement of a general polygon.
 * \param pgn the input polygon
 * \param res the complement of \p pgn
 * \param traits a model of `GeneralPolygonSetTraits_2`
 */
template <typename ArrTraits, typename GpsTraits>
void complement(const General_polygon_2<ArrTraits>& pgn,
                General_polygon_with_holes_2<ArrTraits>& res,
                const GpsTraits& traits);

/*! Computes the complement of a polygon with holes.
 * \param pgn the input polygon
 * \param oi the output iterator for the result.
 *           Its dereference type is `Polygon_with_holes_2<Kernel, Container>`.
 * \param traits a model of `GeneralPolygonSetTraits_2`
 * \return the past-the-end iterator of the output container.
 */
template <typename Kernel, typename Container, typename OutputIterator,
          typename GpsTraits>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                          OutputIterator oi,
                          const GpsTraits& traits);

/*! Computes the complement of the general polygon with holes.
 * \param pgn the input polygon
 * \param oi the output iterator for the result.
 *           Its dereference type is
 *             `General_polygon_with_holes_2<<General_polygon_2<ArrTraits> >`.
 * \param traits a model of `GeneralPolygonSetTraits_2`
 * \return the past-the-end iterator of the output container.
 */
template <typename ArrTraits, typename OutputIterato, typename GpsTraitsr>
OutputIterator complement(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn,
                          OutputIterator oi,
                          const GpsTraits& traits);

/// @}

} /* namespace CGAL */

namespace CGAL {

/*! \addtogroup boolean_difference Difference Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_difference
 *
 * Each one of these functions computes the difference between two given
 * polygons `p1` and `p2`, and inserts the resulting polygons
 * with holes into an output container through the output iterator `oi`.
 * The value type of the `OutputIterator` is either
 * `Polygon_with_holes_2` or
 * `General_polygon_with_holes_2`.
 *
 * The signature of the function is:
 *   - `%OutputIterator %difference(const Type1& p1, const Type2& p2, %OutputIterator oi);`
 *
 * \cgalHeading{Parameters}
 *
 * The types of the parameters of the `difference()` function are any of the
 * following combinations.
 *
 * <div align="left">
 * <table cellpadding=3 border="1">
 * <tr><th> Type1</th><th>Type2</th></tr>
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
 * \sa \link boolean_symmetric_difference `CGAL::symmetric_difference()` \endlink
 */

/// @{

/*! writes the difference of the polygons `p1` and `p2` into the output iterator
 * `oi`.
 *  The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator difference(const Polygon_2<Kernel, Container>& p1,
                          const Polygon_2<Kernel, Container>& p2,
                          OutputIterator oi);

/*! writes the difference of the polygons `p1` and `p2` into the output iterator
 * `oi`.
 *  The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator difference(const Polygon_2<Kernel, Container>& p1,
                          const Polygon_with_holes_2<Kernel,Container>& p2,
                          OutputIterator oi);

/*! writes the difference of the polygons `p1` and `p2` into the output iterator
 * `oi`.
 *  The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container>& p1,
                          const Polygon_2<Kernel, Container>& p2,
                          OutputIterator oi);

/*! writes the difference of the polygons `p1` and `p2` into the output iterator
 * `oi`.
 * The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator difference(const Polygon_with_holes_2<Kernel, Container>& p1,
                          const Polygon_with_holes_2<Kernel, Container>& p2,
                          OutputIterator oi);

/*! writes the difference of the general polygons `p1` and `p2` into the output
 * iterator `oi`.
 * The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator difference(const General_polygon_2<ArrTraits>& p1,
                          const General_polygon_2<ArrTraits>& p2,
                          OutputIterator oi);

/*! writes the difference of the general polygons `p1` and `p2` into the output
 * iterator `oi`.
 *  The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p1,
                          const General_polygon_2<ArrTraits>& p2,
                          OutputIterator oi);

/*! writes the difference of the general polygons `p1` and `p2` into the output
 * iterator `oi`.
 *  The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator difference(const General_polygon_2<ArrTraits>& p1,
                          const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p2,
                          OutputIterator oi);

/*! writes the difference of the general polygons `p1` and `p2` into the output
 * iterator `oi`.
 * The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename Polygon, typename OutputIterator>
OutputIterator difference(const General_polygon_with_holes_2<Polygon>& p1,
                          const General_polygon_with_holes_2<Polygon>& p2,
                          OutputIterator oi);
/// @}

} /* namespace CGAL */

namespace CGAL {

/*! \addtogroup boolean_do_intersect Intersection Testing Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_do_intersect
 *
 * Each one of these functions computes if the interior of two given
 * polygons `p1` and `p2` intersect.
 *
 * The signature of the function is:
 *   - `bool do_intersect(const Type1& p1, const Type2& p2);`
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

/*! returns `true` if the polygons `p1` and `p2` intersect in their interior.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_2<Kernel, Container>& p1,
                  const Polygon_2<Kernel, Container>& p2);

/*! returns `true` if the polygons `p1` and `p2` intersect in their interior.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_2<Kernel, Container>& p1,
                  const Polygon_with_holes_2<Kernel, Container>& p2);

/*! returns `true` if the polygons `p1` and `p2` intersect in their interior.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& p1,
                  const Polygon_2<Kernel, Container>& p2);

/*! returns `true` if the polygons `p1` and `p2` intersect in their interior.
 */
template <typename Kernel, typename Container>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& p1,
                  const Polygon_with_holes_2<Kernel, Container>& p2);

/*! returns `true` if the general polygons `p1` and `p2` intersect in their
 * interior.
 */
template <typename ArrTraits>
bool do_intersect(const General_polygon_2<ArrTraits>& p1,
                  const General_polygon_2<ArrTraits>& p2);

/*! returns `true` if the general polygons `p1` and `p2` intersect in their
 * interior.
 */
template <typename ArrTraits>
bool do_intersect(const General_polygon_2<ArrTraits>& p1,
                  const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p2);

/*! returns `true` if the general polygons `p1` and `p2` intersect in their
 * interior.
 */
template <typename ArrTraits>
bool do_intersect(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p1,
                  const General_polygon_2<ArrTraits>& p2);


/*! returns `true` if the general polygons `p1` and `p2` intersect in their
 * interior.
 */
template <typename Polygon>
bool do_intersect(const General_polygon_with_holes_2<Polygon>& p1,
                  const General_polygon_with_holes_2<Polygon>& p2);

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
 * polygons `p1` and `p2`, inserts the resulting polygons with
 * holes into an output container through a given output iterator
 * `oi`, and returns the output iterator. The value type of the
 * `OutputIterator` is either `Polygon_with_holes_2` or
 * `General_polygon_with_holes_2`.
 *
 * The signature of the function is:
 *   - `%OutputIterator %intersection(const Type1& p1, const Type2& p2,
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

/*! writes the intersection of the polygons `p1` and `p2` into the output iterator `oi`.
 *  The value type of `oi` is `Polygon_with_holes_2`.
 */
OutputIterator intersection(const Type1& p1, const Type2& p2,
                            OutputIterator oi);

/*! writes the intersection of the polygons `p1` and `p2` into the output iterator `oi`.
 * The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_2<Kernel, Container>& p1,
                            const Polygon_2<Kernel, Container>& p2,
                            OutputIterator oi);

/*! writes the intersection of the polygons `p1` and `p2` into the output iterator `oi`.
 * The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_2<Kernel, Container>& p1,
                            const Polygon_with_holes_2<Kernel, Container>& p2,
                            OutputIterator oi);

/*! writes the intersection of the polygons `p1` and `p2` into the output iterator `oi`.
 * The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container>& p1,
                            const Polygon_2<Kernel, Container>& p2,
                            OutputIterator oi);

/*! writes the intersection of the polygons `p1` and `p2` into the output iterator `oi`.
 * The value type of `oi` is `Polygon_with_holes_2`.
 */
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator intersection(const Polygon_with_holes_2<Kernel, Container>& p1,
                            const Polygon_with_holes_2<Kernel, Container>& p2,
                            OutputIterator oi);

/*! writes the intersection of the general polygons `p1` and `p2` into the output iterator `oi`.
 * The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator intersection(const General_polygon_2<ArrTraits>& p1,
                            const General_polygon_2<ArrTraits>& p2,
                            OutputIterator oi);

/*! writes the intersection of the general polygons `p1` and `p2` into the output iterator `oi`.
 * The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator intersection(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p1,
                            const General_polygon_2<ArrTraits>& p2,
                            OutputIterator oi);

/*! writes the intersection of the general polygons `p1` and `p2` into the output iterator `oi`.
 * The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename ArrTraits, typename OutputIterator>
OutputIterator intersection(const General_polygon_2<ArrTraits>& p1,
                            const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p2,
                            OutputIterator oi);

/*! writes the intersection of the general polygons `p1` and `p2` into the output iterator `oi`.
 * The value type of `oi` is `General_polygon_with_holes_2`.
 */
template <typename Polygon, typename OutputIterator>
OutputIterator intersection(const General_polygon_with_holes_2<Polygon>& p1,
                            const General_polygon_with_holes_2<Polygon>& p2,
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
 * Each one of these functions computes the union of two given polygons `p1` and
 * `p2`. If the two given polygons overlap, it returns `true`, and places the
 * resulting polygon in `p`. Otherwise, it returns `false`.
 *
 * The signature of the function is:
 *   - `bool join(const Type1& p1, const Type2& p2, General_polygon_with_holes_2& res);`
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

/*! writes the union of the polygons `p1` and `p2` into the polygon with holes `res`.
 *  Returns `true` if the two given polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_2<Kernel, Container>& p1,
          const Polygon_2<Kernel, Container>& p2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! writes the union of the polygons `p1` and `p2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_2<Kernel, Container>& p1,
          const Polygon_with_holes_2<Kernel,Container>& p2,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! writes the union of the polygons `p1` and `p2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_with_holes_2<Kernel, Container>& p2,
          const Polygon_2<Kernel, Container>& p1,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! writes the union of the polygons `p1` and `p2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename Kernel, typename Container>
bool join(const Polygon_with_holes_2<Kernel, Container>& p2,
          const Polygon_with_holes_2<Kernel, Container>& p1,
          General_polygon_with_holes_2<Polygon_2<Kernel, Container> >& res);


/*! writes the union of the general polygons `p1` and `p2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename ArrTraits>
bool join(const General_polygon_2<ArrTraits>& p1,
          const General_polygon_2<ArrTraits>& p2,
          General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res);


/*! writes the union of the polygons `p1` and `p2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename ArrTraits>
bool join(const General_polygon_2<ArrTraits>& p1,
          const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p2,
          General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res);


/*! writes the union of the general polygons `p1` and `p2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename ArrTraits>
bool join(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p2,
          const General_polygon_2<ArrTraits>& p1,
          General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res);


/*! writes the union of the general polygons `p1` and `p2` into the polygon with holes `res`.
 * Returns `true` if the two given polygons overlap.
 */
template <typename Polygon>
bool join(const General_polygon_with_holes_2<Polygon>& p1,
          const General_polygon_with_holes_2<Polygon>& p2,
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
 * There are several overloaded functions called `Oriented_side()` that computes
 * the relative position of two polygons or of a point and a polygon. This group
 * of functions is  divided into two subgroups.
 *
 * \cgalHeading{Oriented Side of two Polygons}
 *
 * The functions in the first subgroup accept two polygons `pgn1` and `pgn2`.
 * Each function in this group returns `ON_POSITIVE_SIDE` if the two
 * given polygons `pgn1` and `pgn2` intersect in their interiors,
 * `ON_NEGATIVE_SIDE` if `pgn1` and `pgn2` do not intersect at all, and
 * `ON_ORIENTED_BOUNDARY` if `pgn1` and `pgn2` intersect only in their
 * boundaries.
 *
 * This group is further divided into two sub-subgroups.
 * A function in the first sub-subgroup has one of the two following signatures:
 *   - `template <typename Kernel, typename Container> Oriented_side oriented_side(const Type1& pgn1, const Type2& pgn2);`
 *   - `template <typename Kernel, typename Container, typename GpsTraits> Oriented_side oriented_side(const Type1& pgn1, const Type2& pgn2, const GpsTraits& traits);`
 *
 * The types of `pgn1` and `pgn2` must be convertible to one of the two types `Polygon_2<Kernel, Container>` or Polygon_with_holes_2<Kernel, Container>`.
 *
 * A function in the second sub-subgroup has one of the two following signatures:
 *   - `template <typename %ArrTraits> Oriented_side oriented_side(const Type1& pgn1, const Type2& pgn2);`
 *   - `template <typename %ArrTraits, typename GpsTraits> Oriented_side oriented_side(const Type1& pgn1, const Type2& pgn2, const GpsTraits& traits);`
 *
 * The types of `pgn1` and `pgn2` must be convertible to one of the two types `General_polygon_2<%ArrTraits>` or `General_polygon_with_holes_2<General_polygon_2<Arr_traits> >`.
 *
 * \cgalHeading{Oriented Side of a Point and a Polygon}
 *
 * The functions in the second group accept a point `p` and a polygon `pgn`.
 * Each function in this group returns `ON_POSITIVE_SIDE` if the point `p`
 * is in the interior of `pgn`, `ON_NEGATIVE_SIDE` if `p` is in the exterior
 * of `pgn`, and `ON_ORIENTED_BOUNDARY` if `p` is on the boundary of `pgn`.
 * This group is also further divided into two sub-subgroups.
 * A function in the first sub-subgroup has one of the two following signatures:
 *   - `template <typename %ArrTraits> Oriented_side oriented_side(const Point_2& p, const Type& pgn);`
 *   - `template <typename %ArrTraits, typename GpsTraits> Oriented_side oriented_side(const Point_2& p, const Type& pgn, const GpsTraits& traits);`
 *
 * `Type` must be convertible to one of `Polygon_2<Kernel, Container>` or `Polygon_with_holes_2<Kernel, Container>`.
 *
 * A function in the second sub-subgroup has one of the two following signatures:
 *   - `template <typename %ArrTraits> Oriented_side oriented_side(const Point_2& p, const Type2& pgn);`
 *   - `template <typename %ArrTraits, typename GpsTraits> Oriented_side oriented_side(const Point_2& p, const Type2& pgn, const GpsTraits& traits);`
 * `Type` must be convertible to one of `General_polygon_2<%ArrTraits>` or `General_polygon_with_holes_2<General_polygon_2<Arr_traits> >`.
 *
 * \pre `GpsTraits` must be a model of the concept `GeneralPolygonSetTraits_2`.
 * \pre `%ArrTraits` must be a model of the concept `ArrangementDirectionalXMonotoneTraits_2`.
 * \pre A model of `GpsTraits` must derive from a type that is convertible to a model of `%ArrTraits`.
 *
 * \sa \link boolean_do_intersect `CGAL::do_intersect()` \endlink
 */

/// @{

// Polygon--Polygon Traits-less
template <typename Kernel, typename Container>
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& p1,
                            const Polygon_2<Kernel, Container>& p2);


template <typename Kernel, typename Container>
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& p1,
                            const Polygon_with_holes_2<Kernel, Container>& p2);


template <typename Kernel, typename Container>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& p1,
                            const Polygon_2<Kernel, Container>& p2);


template <typename Kernel, typename Container>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& p1,
                            const Polygon_with_holes_2<Kernel, Container>& p2);


template <typename ArrTraits>
Oriented_side oriented_side(const General_polygon_2<ArrTraits>& p1,
                            const General_polygon_2<ArrTraits>& p2);


template <typename ArrTraits>
Oriented_side oriented_side(const General_polygon_2<ArrTraits>& p1,
                            const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p2);


template <typename ArrTraits>
Oriented_side oriented_side(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p1,
                            const General_polygon_2<ArrTraits>& p2);


template <typename Polygon>
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon>& p1,
                            const General_polygon_with_holes_2<Polygon>& p2);

// Polygon--Polygon With Traits
template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& p1,
                            const Polygon_2<Kernel, Container>& p2,
                            const GpsTraits& traits);


template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& p1,
                            const Polygon_with_holes_2<Kernel, Container>& p2,
                            const GpsTraits& traits);


template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& p1,
                            const Polygon_2<Kernel, Container>& p2,
                            const GpsTraits& traits);


template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& p1,
                            const Polygon_with_holes_2<Kernel, Container>& p2,
                            const GpsTraits& traits);


template <typename ArrTraits, typename GpsTraits>
Oriented_side oriented_side(const General_polygon_2<ArrTraits>& p1,
                            const General_polygon_2<ArrTraits>& p2,
                            const GpsTraits& traits);


template <typename ArrTraits, typename GpsTraits>
Oriented_side oriented_side(const General_polygon_2<ArrTraits>& p1,
                            const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p2,
                            const GpsTraits& traits);


template <typename ArrTraits, typename GpsTraits>
Oriented_side oriented_side(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p1,
                            const General_polygon_2<ArrTraits>& p2,
                            const GpsTraits& traits);


template <typename Polygon, typename GpsTraits>
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon>& p1,
                            const General_polygon_with_holes_2<Polygon>& p2,
                            const GpsTraits& traits);

// Point--Polygon Traits-less
template <typename Kernel, typename Container>
Oriented_side oriented_side(const Point_2& p1,
                            const Polygon_2<Kernel, Container>& p2);


template <typename Kernel, typename Container>
Oriented_side oriented_side(const Point_2& p1,
                            const Polygon_with_holes_2<Kernel, Container>& p2);


template <typename ArrTraits>
Oriented_side oriented_side(const Point_2& p1,
                            const General_polygon_2<ArrTraits>& p2);


template <typename ArrTraits>
Oriented_side oriented_side(const Point_2& p1,
                            const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p2);

// Point--Polygon With Traits
template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Point_2& p1,
                            const Polygon_2<Kernel, Container>& p2,
                            const GpsTraits& traits);


template <typename Kernel, typename Container, typename GpsTraits>
Oriented_side oriented_side(const Point_2& p1,
                            const Polygon_with_holes_2<Kernel, Container>& p2,
                            const GpsTraits& traits);


template <typename ArrTraits, typename GpsTraits>
Oriented_side oriented_side(const Point_2& p1,
                            const General_polygon_2<ArrTraits>& p2,
                            const GpsTraits& traits);


template <typename ArrTraits, typename GpsTraits>
Oriented_side oriented_side(const Point_2& p1,
                            const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p2,
                            const GpsTraits& traits);

/// @}
} /* namespace CGAL */

namespace CGAL {

/*! \addtogroup boolean_symmetric_difference Symmetric Difference Functions
 * \ingroup PkgBooleanSetOperations2Ref
 * \anchor ref_bso_symmetric_difference
 *
 * Each one of these functions computes the symmetric difference between
 * two given polygons `p1` and `p2`, and inserts the resulting
 * polygons with holes into an output container through the output
 * iterator `oi`. The value type of the `OutputIterator` is either
 * `Polygon_with_holes_2` or
 * `General_polygon_with_holes_2`.
 *
 * The signature of the function is:
 *   - `%OutputIterator symmetric_difference(const Type1& p1, const Type2& p2, %OutputIterator oi);`
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

OutputIterator symmetric_difference(const Type1& p1, const Type2& p2,
                                    OutputIterator oi);


template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator symmetric_difference(const Polygon_2<Kernel, Container>& p1,
                                    const Polygon_2<Kernel, Container>& p2,
                                    OutputIterator oi);


template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator
symmetric_difference(const Polygon_2<Kernel, Container>& p1,
                     const Polygon_with_holes_2<Kernel, Container>& p2,
                     OutputIterator oi);

template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& p1,
                     const Polygon_2<Kernel, Container>& p2,
                     OutputIterator oi);


template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& p1,
                     const Polygon_with_holes_2<Kernel, Container>& p2,
                     OutputIterator oi);


template <typename ArrTraits, typename OutputIterator>
OutputIterator symmetric_difference(const General_polygon_2<ArrTraits>& p1,
                                    const General_polygon_2<ArrTraits>& p2,
                                    OutputIterator oi);


template <typename ArrTraits, typename OutputIterator>
OutputIterator
symmetric_difference(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p1,
                     const General_polygon_2<ArrTraits>& p2,
                     OutputIterator oi);


template <typename ArrTraits, typename OutputIterator>
OutputIterator
symmetric_difference(const General_polygon_2<ArrTraits>& p1,
                     const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& p2,
                     OutputIterator oi);


template <typename Polygon, typename OutputIterator>
OutputIterator
symmetric_difference(const General_polygon_with_holes_2<Polygon>& p1,
                     const General_polygon_with_holes_2<Polygon>& p2,
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
