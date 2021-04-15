// Copyright (c) 2018  GeometryFactory Sarl (France).
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

#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <iostream>
#include <sstream>
#include <string>


#include <boost/geometry/io/wkt/read.hpp>
#include <boost/geometry/io/wkt/write.hpp>

#include <CGAL/IO/traits_point.h>
#include <CGAL/IO/traits_point_3.h>
#include <CGAL/IO/traits_linestring.h>
#include <CGAL/IO/traits_polygon.h>
#include <CGAL/IO/traits_multipoint.h>
#include <CGAL/IO/traits_multilinestring.h>
#include <CGAL/IO/traits_multipolygon.h>


namespace CGAL{

namespace internal {

  template <typename K>
  void pop_back_if_equal_to_front(CGAL::Polygon_2<K>& poly)
  {
    typename CGAL::Polygon_2<K>::iterator it = poly.end();
    --it;
    if( (*poly.begin()) == *it){
      poly.erase(it);
    }
  }

  template <typename K>
  void pop_back_if_equal_to_front(CGAL::Polygon_with_holes_2<K>& pwh)
  {
    pop_back_if_equal_to_front(pwh.outer_boundary());
    for(auto i = pwh.holes_begin(); i!= pwh.holes_end(); ++i){
      pop_back_if_equal_to_front(*i);
    }
  }
}

template<typename Point>
std::istream&
read_point_WKT( std::istream& in,
                Point& point)
{
  if(!in)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return in;
  }

  std::string line;
  while(std::getline(in, line))
  {
    std::istringstream iss(line);
    std::string type;
    iss >> type;

    if(type.substr(0, 5).compare("POINT")==0)
    {
      boost::geometry::read_wkt(line, point);
      break;
    }
  }
  return in;
}

template<typename MultiPoint>
std::istream&
read_multi_point_WKT( std::istream& in,
                     MultiPoint& mp)
{
  if(!in)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return in;
  }
  internal::Geometry_container<MultiPoint, boost::geometry::multi_point_tag> gc(mp);
  std::string line;
  while(std::getline(in, line))
  {
    std::istringstream iss(line);
    std::string type;
    iss >> type;

    if(type.substr(0, 10).compare("MULTIPOINT")==0)
    {
      boost::geometry::read_wkt(line, gc);
      break;
    }
  }
  return in;
}


template<typename LineString>
std::istream&
read_linestring_WKT( std::istream& in,
                     LineString& polyline)
{
  if(!in)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return in;
  }
  internal::Geometry_container<LineString, boost::geometry::linestring_tag> gc(polyline);
  std::string line;
  while(std::getline(in, line))
  {
    std::istringstream iss(line);
    std::string type;
    iss >> type;

    if(type.substr(0, 10).compare("LINESTRING")==0)
    {
      boost::geometry::read_wkt(line, gc);
      break;
    }
  }
  return in;
}

template<typename MultiLineString>
std::istream&
read_multi_linestring_WKT( std::istream& in,
                          MultiLineString& mls)
{
  if(!in)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return in;
  }
  typedef typename MultiLineString::value_type PointRange;
  typedef internal::Geometry_container<PointRange, boost::geometry::linestring_tag> LineString;

  std::vector<LineString> pr_range;
  internal::Geometry_container<std::vector<LineString>, boost::geometry::multi_linestring_tag> gc(pr_range);
  std::string line;
  while(std::getline(in, line))
  {
    std::istringstream iss(line);
    std::string type;
    iss >> type;

    if(type.substr(0, 15).compare("MULTILINESTRING")==0)
    {
      boost::geometry::read_wkt(line, gc);
      break;
    }
  }
  for(LineString& ls : gc)
  {
    mls.push_back(*ls.range);
  }
  return in;
}

template<typename Polygon>
std::istream&
read_polygon_WKT( std::istream& in,
                  Polygon& polygon
                  )
{
  if(!in)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return in;
  }

  std::string line;
  while(std::getline(in, line))
  {
    std::istringstream iss(line);
    std::string type;
    iss >> type;

    if(type.substr(0, 7).compare("POLYGON")==0)
    {
      try {
        boost::geometry::read_wkt(line, polygon);
      } catch( ...){
        in.setstate(std::ios::failbit);
        return in;
      };
      internal::pop_back_if_equal_to_front(polygon);
      break;
    }
  }
  return in;
}

template<typename MultiPolygon>
std::istream&
read_multi_polygon_WKT( std::istream& in,
                        MultiPolygon& polygons
                        )
{
  if(!in)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return in;
  }
  internal::Geometry_container<MultiPolygon, boost::geometry::multi_polygon_tag> gc(polygons);
  std::string line;
  while(std::getline(in, line))
  {
    std::istringstream iss(line);
    std::string type;
    iss >> type;

    if(type.substr(0, 12).compare("MULTIPOLYGON")==0)
    {
      try {
              boost::geometry::read_wkt(line, gc);
            } catch( ...){
              in.setstate(std::ios::failbit);
              return in;
            };
      for( typename
           internal::Geometry_container<MultiPolygon, boost::geometry::multi_polygon_tag>::iterator it
          = gc.begin(); it != gc.end(); ++it)
      {
        internal::pop_back_if_equal_to_front(*it);
      }
      break;
    }
  }
  return in;
}

template<typename Point>
std::ostream&
write_point_WKT( std::ostream& out,
           const Point& point
           )
{
  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return out;
  }
  out<<boost::geometry::wkt(point)<<std::endl;
  return out;
}

template<typename Polygon>
std::ostream&
write_polygon_WKT( std::ostream& out,
           const Polygon& poly
           )
{
  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return out;
  }
  out<<boost::geometry::wkt(poly)<<std::endl;
  return out;
}

template<typename LineString>
std::ostream&
write_linestring_WKT( std::ostream& out,
            LineString ls
           )
{
  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return out;
  }
  internal::Geometry_container<LineString, boost::geometry::linestring_tag> gc(ls);
  out<<boost::geometry::wkt(gc)<<std::endl;
  return out;
}

template<typename MultiPoint>
std::ostream&
write_multi_point_WKT( std::ostream& out,
            MultiPoint& mp
           )
{
  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return out;
  }
  internal::Geometry_container<MultiPoint, boost::geometry::multi_point_tag> gc(mp);
  out<<boost::geometry::wkt(gc)<<std::endl;
  return out;
}

template<typename MultiPolygon>
std::ostream&
write_multi_polygon_WKT( std::ostream& out,
            MultiPolygon& polygons
           )
{
  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return out;
  }
  internal::Geometry_container<MultiPolygon, boost::geometry::multi_polygon_tag> gc(polygons);
  out<<boost::geometry::wkt(gc)<<std::endl;
  return out;
}

template<typename MultiLineString>
std::ostream&
write_multi_linestring_WKT( std::ostream& out,
            MultiLineString& mls
           )
{
  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return out;
  }
  typedef typename MultiLineString::value_type PointRange;
  typedef internal::Geometry_container<PointRange, boost::geometry::linestring_tag> LineString;
  std::vector<LineString> pr_range;
  for(PointRange& pr : mls)
  {
    LineString ls(pr);
    pr_range.push_back(ls);
  }
  internal::Geometry_container<std::vector<LineString>, boost::geometry::multi_linestring_tag> gc(pr_range);
  out<<boost::geometry::wkt(gc)<<std::endl;
  return out;
}

template<typename MultiPoint,
         typename MultiLineString,
         typename MultiPolygon>
std::istream&
read_WKT( std::istream& input,
          MultiPoint& points,
          MultiLineString& polylines,
          MultiPolygon& polygons)
{
  do
  {
    typedef typename MultiPoint::value_type Point;
    typedef typename MultiLineString::value_type LineString;
    typedef typename MultiPolygon::value_type Polygon;
    std::string line;
    std::streampos input_pos = input.tellg();
    std::getline(input, line);
    std::istringstream iss(line);
    std::string t;
    std::string type="";
    iss >> t;
    for(std::size_t pos=0; pos < t.length(); ++pos)
    {
      char c=t[pos];
      if(c=='(')
        break;
      type.push_back(c);
    }
    input.seekg(input_pos);
    if(type == "POINT")
    {
      Point p;
      CGAL::read_point_WKT(input, p);
      points.push_back(p);
    }
    else if(type == "LINESTRING")
    {
      LineString l;
      CGAL::read_linestring_WKT(input, l);
      polylines.push_back(l);
    }
    else if(type == "POLYGON")
    {
      Polygon p;
      CGAL::read_polygon_WKT(input, p);
      if(!p.outer_boundary().is_empty())
        polygons.push_back(p);
    }
    else if(type == "MULTIPOINT")
    {
      MultiPoint mp;
      CGAL::read_multi_point_WKT(input, mp);
      for(const Point& point : mp)
          points.push_back(point);
    }
    else if(type == "MULTILINESTRING")
    {
      MultiLineString mls;
      CGAL::read_multi_linestring_WKT(input, mls);
      for(const LineString& ls : mls)
          polylines.push_back(ls);
    }
    else if(type == "MULTIPOLYGON")
    {
      MultiPolygon mp;
      CGAL::read_multi_polygon_WKT(input, mp);
      for(const Polygon& poly : mp)
          polygons.push_back(poly);
    }
  }while(input.good() && !input.eof());
  return input;
}
}//end CGAL
#endif
#endif
