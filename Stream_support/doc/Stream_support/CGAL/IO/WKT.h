namespace CGAL{
//! \ingroup PkgStreamSupportRef
//! \brief `read_point_WKT()` fills a `Point` from a  WKT stream. The first line starting with POINT
//! in the stream will be used.
//!
//! \tparam Point can be a `CGAL::Point_2` or `CGAL::Point_3`.
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//!
//! \see `CGAL::Point_2`
//! \see `CGAL::Point_3`
template<typename Point>
std::istream&
read_point_WKT( std::istream& in,
                Point& point );

//! \ingroup PkgStreamSupportRef
//! \brief `read_multi_point_WKT()` overwrites the content of a `MultiPoint`
//! with the first line starting with MULTIPOINT in the stream.
//!
//! \tparam MultiPoint must be a model of `RandomAccessRange` of `CGAL::Point_2` or `CGAL::Point_3`,
//! and have:
//! - a function `push_back()` that takes the same point type,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//!
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//! \see `CGAL::Point_2`
//! \see `CGAL::Point_3`
template<typename MultiPoint>
std::istream&
read_multi_point_WKT( std::istream& in,
                      MultiPoint& mp );


//! \ingroup PkgStreamSupportRef
//! \brief `read_linestring_WKT()` fills a `Linestring` from a WKT stream.
//! The first line starting with LINESTRING in the stream will be used.
//!
//! \tparam Linestring must be a model of `RandomAccessRange` of `CGAL::Point_2`,
//! and have:
//! - a function `push_back()` that takes a `CGAL::Point_2`.
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//! \see `CGAL::Point_2`
template<typename LineString>
std::istream&
read_linestring_WKT( std::istream& in,
                     LineString& polyline );

//! \ingroup PkgStreamSupportRef
//! \brief `read_multi_linestring_WKT()` overwrites the content of a `MultiLineString`
//! with the first line starting with MULTILINESTRING in the stream.
//!
//! \tparam MultiLineString must be a model of `RandomAccessRange` of `Linestring`,
//! and have:
//! - a function `push_back()` that takes a `Linestring`,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//!
//! \see `CGAL::Point_2`
template<typename MultiLineString>
std::istream&
read_multi_linestring_WKT( std::istream& in,
                           MultiLineString& mls );

//! \ingroup PkgStreamSupportRef
//! \brief `read_polygon_WKT()` fills `polygon` from a WKT stream.
//! The first line starting with POLYGON in the stream will be used.
//!
//! \tparam Polygon is a `CGAL::General_polygon_with_holes_2`.
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//!
//! \see `CGAL::General_polygon_with_holes_2`
template<typename Polygon>
std::istream&
read_polygon_WKT( std::istream& in,
                  Polygon& polygon );

//! \ingroup PkgStreamSupportRef
//! \brief `read_multi_polygon_WKT()` overwrites the content of a `MultiPolygon`
//! with the first line starting with MULTIPOLYGON in the stream.
//!
//! \tparam MultiPolygon must be a model of `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`,
//! and have:
//! - a function `push_back()` that takes a `CGAL::General_polygon_with_holes_2`,
//! - a function `clear()`,
//! - a function `resize()` that takes an `size_type`
//! - an `operator[]()` that takes a `size_type`.
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//! \see `CGAL::General_polygon_with_holes_2`

template<typename MultiPolygon>
std::istream&
read_multi_polygon_WKT( std::istream& in,
                        MultiPolygon& polygons );

//! \ingroup PkgStreamSupportRef
//! \brief `write_point_WKT()` writes `point` into a WKT stream.
//! \tparam Point is a `CGAL::Point_2`
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//! \see `CGAL::Point_2`
template<typename Point>
std::ostream&
write_point_WKT( std::ostream& out,
                 const Point& point );

//! \ingroup PkgStreamSupportRef
//! \brief `write_polygon_WKT()` writes `poly` into a WKT stream.
//! \tparam Polygon  must be a `CGAL::General_polygon_with_holes_2`
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//! \see `CGAL::General_polygon_with_holes_2`
template<typename Polygon>
std::ostream&
write_polygon_WKT( std::ostream& out,
                   const Polygon& poly );

//! \ingroup PkgStreamSupportRef
//! \brief `write_linestring_WKT()` writes the content of `ls`
//! into a WKT stream.
//! \tparam LineString must be a `RandomAccessRange` of `CGAL::Point_2`.
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//!\see `CGAL::Point_2`
template<typename LineString>
std::ostream&
write_linestring_WKT( std::ostream& out,
                      LineString ls );

//! \ingroup PkgStreamSupportRef
//! \brief `write_multi_point_WKT()` writes the content of `mp`
//! into a WKT stream.
//! \tparam MultiPoint must be a `RandomAccessRange` of `CGAL::Point_2`.
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//!\see `CGAL::Point_2`
template<typename MultiPoint>
std::ostream&
write_multi_point_WKT( std::ostream& out,
                       MultiPoint& mp );

//! \ingroup PkgStreamSupportRef
//! \brief `write_multi_polygon_WKT()` writes the content of `polygons`
//! into a WKT stream.
//! \tparam MultiPolygon must be a `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`.
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//!\see `CGAL::General_polygon_with_holes_2`
template<typename MultiPolygon>
std::ostream&
write_multi_polygon_WKT( std::ostream& out,
                         MultiPolygon& polygons );

//! \ingroup PkgStreamSupportRef
//! \brief `write_multi_linestring_WKT()` writes the content of `mls`
//! into a WKT stream.
//! \tparam MultiLineString must be a `RandomAccessRange` of `LineString`.
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//! \see `CGAL::write_linestring_WKT()`
template<typename MultiLineString>
std::ostream&
write_multi_linestring_WKT( std::ostream& out,
                            MultiLineString& mls );

//! \ingroup PkgStreamSupportRef
//!  reads the content of a WKT stream and fills
//! `points`, `polylines` and `polygons` with all the POINT, MULTIPOINT,
//! LINESTRING, MULTILINESTRING, POLYGON and MULTIPOLYGON it finds in `input`.
//! \tparam MultiPoint must be a model of `RandomAccessRange` of `CGAL::Point_2` or `CGAL::Point_3`.
//! \tparam MultiLineString must be a `RandomAccessRange` of `Linestring`.
//! \tparam MultiPolygon must be a model of `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`.
//! \attention Only Cartesian Kernels with double or float  as `FT` are supported.
//! \attention This function is only available with boost versions starting at 1.56.
//! \see `CGAL::read_linestring_WKT()`
template<typename MultiPoint,
         typename MultiLineString,
         typename MultiPolygon>
std::istream&
read_WKT( std::istream& input,
                        MultiPoint& points,
                        MultiLineString& polylines,
                        MultiPolygon& polygons);
}
