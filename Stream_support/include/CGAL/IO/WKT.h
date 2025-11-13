// Copyright (c) 2018-2020  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_WKT_H
#define CGAL_IO_WKT_H

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>

#include <CGAL/IO/WKT/traits_point.h>
#include <CGAL/IO/WKT/traits_point_3.h>
#include <CGAL/IO/WKT/traits_linestring.h>
#include <CGAL/IO/WKT/traits_polygon.h>
#include <CGAL/IO/WKT/traits_multipoint.h>
#include <CGAL/IO/WKT/traits_multilinestring.h>
#include <CGAL/IO/WKT/traits_multipolygon.h>

#include <boost/geometry/io/wkt/read.hpp>
#include <boost/geometry/io/wkt/write.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace CGAL {
namespace IO {
namespace internal {

template <typename K>
void pop_back_if_equal_to_front(CGAL::Polygon_2<K>& poly)
{
  auto last_it = std::prev(poly.end());
  if((*poly.begin()) == *last_it)
    poly.erase(last_it);
}

template <typename K>
void pop_back_if_equal_to_front(CGAL::Polygon_with_holes_2<K>& pwh)
{
  pop_back_if_equal_to_front(pwh.outer_boundary());
  for(auto& hole : pwh.holes())
    pop_back_if_equal_to_front(hole);
}

template <typename Geometry>
bool read_wkt_or_fail_stream(std::istream& in,
                             const std::string& line,
                             Geometry& geometry)
{
  try {
    boost::geometry::read_wkt(line, geometry);
  } catch(std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    in.clear(in.rdstate() | std::ios::failbit);
    return false;
  }
  return true;
}

bool get_a_new_line(std::istream& in, std::string& line)
{
  in >> std::ws; // skip whitespaces
  if(in.good()) {
    return !std::getline(in, line).fail();
  } else {
    return false;
  }
}

} // namespace internal

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief fills a `Point` from a WKT stream.
//!
//! The first line starting with POINT in the stream will be used.
//!
//! \tparam Point can be a `CGAL::Point_2` or `CGAL::Point_3`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float` as `FT` are supported.
//!
//! \see `CGAL::Point_2`
//! \see `CGAL::Point_3`
template<typename Point>
bool read_point_WKT(std::istream& in,
                    Point& point)
{
  std::string line;
  while(internal::get_a_new_line(in, line))
  {
    if(line.substr(0, 5).compare("POINT") == 0)
    {
      internal::read_wkt_or_fail_stream(in, line, point);
      break;
    }
  }

  return !in.fail();
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief overwrites the content of a `MultiPoint` with the first line starting with MULTIPOINT in the stream.
//!
//! \tparam MultiPoint must be a model of `RandomAccessRange` of `CGAL::Point_2` or `CGAL::Point_3`,
//! and have:
//! - a function `push_back()` that takes the same point type,
//! - a function `clear()`,
//! - a function `resize()` that takes a `size_type`
//! - an `operator[]()` that takes a `size_type`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//! \see `CGAL::Point_2`
//! \see `CGAL::Point_3`
template<typename MultiPoint>
bool read_multi_point_WKT(std::istream& in,
                          MultiPoint& mp)
{
  std::string line;
  bool found = false;
  while(internal::get_a_new_line(in, line))
  {
    if(line.substr(0, 10).compare("MULTIPOINT") == 0)
    {
      CGAL::internal::Geometry_container<MultiPoint, boost::geometry::multi_point_tag> gc(mp);
      found = internal::read_wkt_or_fail_stream(in, line, gc);
      break;
    }
  }
  if(! found){
    return false;
  }
  return !in.fail();
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief fills a `Linestring` from a WKT stream.
//!
//! The first line starting with LINESTRING in the stream will be used.
//!
//! \tparam Linestring must be a model of `RandomAccessRange` of `CGAL::Point_2` or `CGAL::Point_3`,
//! and have:
//! - a function `push_back()` that takes a point.
//! - a function `clear()`,
//! - a function `resize()` that takes a `size_type`
//! - an `operator[]()` that takes a `size_type`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//! \see `CGAL::Point_2`
//! \see `CGAL::Point_3`
template<typename LineString>
bool read_linestring_WKT(std::istream& in,
                         LineString& polyline)
{
  std::string line;
  bool found = false;
  while(internal::get_a_new_line(in, line))
  {
    if(line.substr(0, 10).compare("LINESTRING") == 0)
    {
      CGAL::internal::Geometry_container<LineString, boost::geometry::linestring_tag> gc(polyline);
      found = internal::read_wkt_or_fail_stream(in, line, gc);
      break;
    }
  }

  if(! found){
    return false;
  }
  return !in.fail();
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief overwrites the content of a `MultiLineString` with the first line starting with MULTILINESTRING in the stream.
//!
//! \tparam MultiLineString must be a model of `RandomAccessRange` of `Linestring`,
//! and have:
//! - a function `push_back()` that takes a `Linestring`,
//! - a function `clear()`,
//! - a function `resize()` that takes a `size_type`
//! - an `operator[]()` that takes a `size_type`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
template<typename MultiLineString>
bool read_multi_linestring_WKT(std::istream& in,
                               MultiLineString& mls)
{
  std::string line;
  bool found = false;
  while(internal::get_a_new_line(in, line))
  {
    if(line.substr(0, 15).compare("MULTILINESTRING") == 0)
    {
      using PointRange = typename MultiLineString::value_type;
      using LineString = CGAL::internal::Geometry_container<PointRange, boost::geometry::linestring_tag>;

      std::vector<LineString> pr_range;
      CGAL::internal::Geometry_container<std::vector<LineString>, boost::geometry::multi_linestring_tag> gc(pr_range);

      found = internal::read_wkt_or_fail_stream(in, line, gc);
      for(LineString& ls : gc) {
        mls.push_back(*ls.range);
      }

      break;
    }
  }
  if(! found){
    return false;
  }
  return !in.fail();
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief fills `polygon` from a WKT stream.
//!
//! The first line starting with POLYGON in the stream will be used.
//!
//! \tparam Polygon is a `CGAL::General_polygon_with_holes_2`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//! \see `CGAL::General_polygon_with_holes_2`
template<typename Polygon>
bool read_polygon_WKT(std::istream& in,
                      Polygon& polygon)
{
  std::string line;
  bool found = false;
  while(internal::get_a_new_line(in, line))
  {
    if(line.substr(0, 7).compare("POLYGON") == 0)
    {
      found = internal::read_wkt_or_fail_stream(in, line, polygon);
      internal::pop_back_if_equal_to_front(polygon);
      break;
    }
  }
  if(! found){
    return false;
  }
  return !in.fail();
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief overwrites the content of a `MultiPolygon` with the first line starting with MULTIPOLYGON in the stream.
//!
//! \tparam MultiPolygon must be a model of `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`,
//! and have:
//! - a function `push_back()` that takes a `CGAL::General_polygon_with_holes_2`,
//! - a function `clear()`,
//! - a function `resize()` that takes a `size_type`
//! - an `operator[]()` that takes a `size_type`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//! \see `CGAL::General_polygon_with_holes_2`
template<typename MultiPolygon>
bool read_multi_polygon_WKT(std::istream& in,
                            MultiPolygon& polygons)
{
  std::string line;
  bool found = false;
  while(internal::get_a_new_line(in, line))
  {
    if(line.substr(0, 12).compare("MULTIPOLYGON") == 0)
    {
      CGAL::internal::Geometry_container<MultiPolygon, boost::geometry::multi_polygon_tag> gc(polygons);
      found = internal::read_wkt_or_fail_stream(in, line, gc);

      for(auto& p : gc)
        internal::pop_back_if_equal_to_front(p);

      break;
    }
  }
  if(! found){
    return false;
  }
  return !in.fail();
}


template<typename Kernel, typename Container>
bool read_multi_polygon_WKT(std::istream& in,
                            Multipolygon_with_holes_2<Kernel,Container>& mp)
{
  return read_multi_polygon_WKT(in, mp.polygons_with_holes());
}


//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief writes `point` into a WKT stream.
//!
//! \tparam Point is a `CGAL::Point_2` or `CGAL::Point_3`
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//! \see `CGAL::Point_2`
//! \see `CGAL::Point_3`
template<typename Point>
std::ostream& write_point_WKT(std::ostream& out,
                              const Point& point)
{
  if(!out.good())
    return out;

  out << boost::geometry::wkt(point) << std::endl;
  return out;
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief writes `poly` into a WKT stream.
//!
//! \tparam Polygon must be a `CGAL::General_polygon_with_holes_2`
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//! \see `CGAL::General_polygon_with_holes_2`
template<typename Polygon>
std::ostream& write_polygon_WKT(std::ostream& out,
                                const Polygon& poly)
{
  if(!out.good())
    return out;

  out << boost::geometry::wkt(poly) << std::endl;
  return out;
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief writes the content of `ls` into a WKT stream.
//!
//! \tparam LineString must be a `RandomAccessRange` of `CGAL::Point_2` or `CGAL::Point_3`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//!\see `CGAL::Point_2`
//!\see `CGAL::Point_3`
template<typename LineString>
std::ostream& write_linestring_WKT(std::ostream& out,
                                   LineString ls)
{
  if(!out.good())
    return out;

  CGAL::internal::Geometry_container<LineString, boost::geometry::linestring_tag> gc(ls);
  out << boost::geometry::wkt(gc) << std::endl;
  return out;
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief writes the content of `mp` into a WKT stream.
//!
//! \tparam MultiPoint must be a `RandomAccessRange` of `CGAL::Point_2` or `CGAL::Point_3`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//!\see `CGAL::Point_2`
//!\see `CGAL::Point_2`
template<typename MultiPoint>
std::ostream& write_multi_point_WKT(std::ostream& out,
                                    MultiPoint& mp)
{
  if(!out.good())
    return out;

  CGAL::internal::Geometry_container<MultiPoint, boost::geometry::multi_point_tag> gc(mp);
  out << boost::geometry::wkt(gc) << std::endl;
  return out;
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief writes the content of `polygons` into a WKT stream.
//!
//! \tparam MultiPolygon must be a `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//!\see `CGAL::General_polygon_with_holes_2`
template<typename MultiPolygon>
std::ostream& write_multi_polygon_WKT(std::ostream& out,
                                      MultiPolygon& polygons)
{
  if(!out.good())
    return out;

  CGAL::internal::Geometry_container<MultiPolygon, boost::geometry::multi_polygon_tag> gc(polygons);
  out << boost::geometry::wkt(gc) << std::endl;
  return out;
}

template<typename Kernel, typename Container>
std::ostream& write_multi_polygon_WKT(std::ostream& out,
                                      Multipolygon_with_holes_2<Kernel,Container>& mp)
{
  return write_multi_polygon_WKT(out, mp.polygons_with_holes());
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! \brief writes the content of `mls` into a WKT stream.
//!
//! \tparam MultiLineString must be a `RandomAccessRange` of `LineString`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//! \see `CGAL::IO::write_linestring_WKT()`
template<typename MultiLineString>
std::ostream& write_multi_linestring_WKT(std::ostream& out,
                                         MultiLineString& mls)
{
  if(!out.good())
    return out;

  typedef typename MultiLineString::value_type                                      PointRange;
  typedef CGAL::internal::Geometry_container<PointRange, boost::geometry::linestring_tag> LineString;

  std::vector<LineString> pr_range;
  for(PointRange& pr : mls)
  {
    LineString ls(pr);
    pr_range.push_back(ls);
  }

  CGAL::internal::Geometry_container<std::vector<LineString>, boost::geometry::multi_linestring_tag> gc(pr_range);
  out << boost::geometry::wkt(gc) << std::endl;

  return out;
}

//! \ingroup PkgStreamSupportIoFuncsWKT
//!
//! reads the content of a WKT stream and fills `points`, `polylines` and `polygons`
//! with all the POINT, MULTIPOINT, LINESTRING, MULTILINESTRING, POLYGON and MULTIPOLYGON it finds in `input`.
//!
//! \tparam MultiPoint must be a model of `RandomAccessRange` of `CGAL::Point_2` or `CGAL::Point_3`.
//! \tparam MultiLineString must be a `RandomAccessRange` of `Linestring`.
//! \tparam MultiPolygon must be a model of `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`.
//!
//! \attention Only %Cartesian Kernels with `double`or `float`  as `FT` are supported.
//!
//! \see `CGAL::IO::read_linestring_WKT()`
template<typename MultiPoint,
         typename MultiLineString,
         typename MultiPolygon>
bool read_WKT(std::istream& is,
              MultiPoint& points,
              MultiLineString& polylines,
              MultiPolygon& polygons)
{
  auto fail = [&is]() { is.clear(is.rdstate() | std::ios::failbit); return false; };

  bool found = false;
  std::string line;
  while(is >> std::ws && is.good() && std::getline(is, line))
  {
    typedef typename MultiPoint::value_type Point;
    typedef typename MultiLineString::value_type LineString;
    typedef typename MultiPolygon::value_type Polygon;

    std::string::size_type header_end = line.find("("); // }
    if(header_end == std::string::npos){
      continue;
    }
    std::string type="";
    const std::string header = line.substr(0,header_end);
    const std::string types[6] = { "MULTIPOLYGON", "MULTILINESTRING", "MULTIPOINT", "POLYGON", "LINESTRING", "POINT"};
    for(int i= 0; i < 6; ++i){
      if(header.find(types[i]) != std::string::npos){
        type = types[i];
        break;
      }
    }
    if(type == ""){
      continue;
    }
    std::istringstream iss(line);

    if(type == "POINT")
    {
      Point p;
      if(!IO::read_point_WKT(iss, p) ) return fail();
      found = true;
      points.push_back(p);
    }
    else if(type == "LINESTRING")
    {
      LineString l;
      if(!IO::read_linestring_WKT(iss, l)) return fail();
      found = true;
      polylines.push_back(std::move(l));
    }
    else if(type == "POLYGON")
    {
      Polygon p;
      if(!IO::read_polygon_WKT(iss, p)) return fail();
      if(!p.outer_boundary().is_empty()){
        found = true;
        polygons.push_back(std::move(p));
      }
    }
    else if(type == "MULTIPOINT")
    {
      MultiPoint mp;
      if(!IO::read_multi_point_WKT(iss, mp)) return fail();
      for(const Point& point : mp){
        points.push_back(point);
        found = true;
      }
    }
    else if(type == "MULTILINESTRING")
    {
      MultiLineString mls;
      if(!IO::read_multi_linestring_WKT(iss, mls)) return fail();
      for(LineString& ls : mls){
        polylines.push_back(std::move(ls));
        found = true;
      }
    }
    else if(type == "MULTIPOLYGON")
    {
      MultiPolygon mp;
      if(!IO::read_multi_polygon_WKT(iss, mp)) return fail();
      for(Polygon& poly : mp){
        polygons.push_back(std::move(poly));
        found = true;
      }
    }
  }

  if(!found){
    return false;
  }

  return !is.fail();
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE
using IO::read_linestring_WKT;
using IO::read_multi_linestring_WKT;
using IO::read_multi_point_WKT;
using IO::read_multi_polygon_WKT;
using IO::read_point_WKT;
using IO::read_polygon_WKT;
using IO::read_WKT;
using IO::write_linestring_WKT;
using IO::write_multi_linestring_WKT;
using IO::write_multi_point_WKT;
using IO::write_multi_polygon_WKT;
using IO::write_point_WKT;
using IO::write_polygon_WKT;
#endif

} // namespace CGAL

#endif // CGAL_IO_WKT_H
