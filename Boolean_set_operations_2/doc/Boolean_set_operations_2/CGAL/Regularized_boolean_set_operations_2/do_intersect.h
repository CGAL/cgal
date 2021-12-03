namespace CGAL {
namespace Regularized_boolean_set_operations_2 {

/*! \addtogroup boolean_do_intersect Intersection Testing Functions
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
 * \sa \link boolean_complement `CGAL::Regularized_boolean_operations_2::complement()` \endlink
 * \sa \link boolean_intersection `CGAL::Regularized_boolean_operations_2::intersection()` \endlink
 * \sa \link boolean_join `CGAL::Regularized_boolean_operations_2::join()` \endlink
 * \sa \link boolean_difference `CGAL::Regularized_boolean_operations_2::difference()` \endlink
 * \sa \link boolean_symmetric_difference `CGAL::Regularized_boolean_operations_2::symmetric_difference()` \endlink
 * \sa \sa `do_intersect_polygons_grp`
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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2,
                  UsePolylines = Tag_true());

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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2,
                  UsePolylines = Tag_true());

/*! determines whether two polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2,
                  UsePolylines = Tag_true());

/*! determines whether two polygons with holes intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
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
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
 *         otherwise.
 */
template <typename Kernel, typename Container, typename UsePolylines>
bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2,
                  UsePolylines = Tag_true());

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
             const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2);

/*! determines whether two general polygons intersect in their interior.
 * \param pgn1 the 1st input polygon.
 * \param pgn2 the 2nd input polygon.
 * \return `true` if `pgn1` and `pgn2` intersect in their interiro and `false`
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
 *        either `Polygon_2` (resp. `General_polygon_2`) or
 *        `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or
 *        `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \return `true` if the pairwise intersections of all open polygons or polygons
 *         with holes (resp. general polygons or general polygons with holes) in
 *         the range [*begin,*end) overlap, and `false` otherwise.
 */
template <typename InputIterator>
bool do_intersect(InputIterator begin, InputIterator end);

/*! Given a range of polygons or a range of polygons with holes (resp. a range
 * of general polygons or a range of general polygons with holes) determines
 * whether the open polygons (resp. general polygons) in the range have a common
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
 *        either `Polygon_2` (resp. `General_polygon_2`) or
 *        `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or
 *        `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \return `true` if the pairwise intersections of all open polygons or polygons
 *         with holes (resp. general polygons or general polygons with holes) in
 *         the range [*begin,*end) overlap, and `false` otherwise.
 */
template <typename InputIterator, typename UsePolylines>
bool do_intersect(InputIterator begin, InputIterator end,
                  UsePolylines = Tag_true());

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
bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                  InputIterator2 begin2, InputIterator2 end2);

/*! Given a range of polygons (resp. general polygons) and a range of polygons
 * with holes (resp. general polygons with holes) determines whether the open
 * polygons (resp. general polygons) in the two ranges have a common point.
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
             const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn2,
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
do_intersect(const General_polygon_with_holes_2<General_polygon_2<ArrTraits>>& pgn1,
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
 *        either `Polygon_2` (resp. `General_polygon_2`) or
 *        `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
 * \param end the past-the-end iterator of the input range. Its value type is
 *        either `Polygon_2` (resp. `General_polygon_2`) or
 *        `Polygon_with_holes_2` (resp. `General_polygon_with_holes_2`).
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
} } /* namespace CGAL::Regularized_boolean_set_operations_2 */
