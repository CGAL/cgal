namespace CGAL{
/*!
 * \ingroup IOstreamFunctions
 * writes the content of `in` in `points` and `polygons`, in the OFF format.
 * \see \ref IOStreamOFF
 */
template <class Point_3, class Polygon_3>
bool
read_OFF( std::istream& in,
          std::vector< Point_3 >& points,
          std::vector< Polygon_3 >& polygons,
          bool /* verbose */ = false);

/*!
 * \ingroup IOstreamFunctions
 * writes the content of `points` and `polygons` in `out`, in the OFF format.
 * \see \ref IOStreamOFF
 */
template <class Point_3, class Polygon_3>
bool
write_OFF(std::ostream& out,
          std::vector< Point_3 >& points,
          std::vector< Polygon_3 >& polygons);
}
