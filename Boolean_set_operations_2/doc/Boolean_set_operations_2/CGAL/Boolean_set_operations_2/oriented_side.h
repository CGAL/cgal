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
