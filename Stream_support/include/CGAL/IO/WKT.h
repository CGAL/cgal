// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_WKT_H
#define CGAL_WKT_H

#include <iostream>
#include <string>

#include <boost/geometry.hpp>
#include <boost/foreach.hpp>

#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <CGAL/IO/traits_point.h>
#include <CGAL/IO/traits_linestring.h>
#include <CGAL/IO/traits_polygon.h>
#include <CGAL/IO/traits_multipoint.h>
#include <CGAL/IO/traits_multilinestring.h>
#include <CGAL/IO/traits_multipolygon.h>


//! \todo Should they be in a namespace CGAL::IO or just CGAL ?
namespace CGAL{
//! \ingroup PkgIOstreams
//! \brief read_point_WKT reads the content of a .wkt file into a `Point` if possible.
//! 
//! A `Point` must inherit `CGAL::Point_2`.
//! 
//! \see CGAL::Point_2
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

//! \ingroup PkgIOstreams
//! \brief read_multipoint_WKT reads the content of a .wkt file into a `MultiPoint` if possible.
//! 
//! A `MultiPoint` must be a model of `RandomAccessRange` of `CGAL::Point_2`.
//! 
//! \see CGAL::Point_2
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


//! \ingroup PkgIOstreams
//! \brief read_linestring_WKT reads the content of a .wkt file into 
//! a `Linestring` if possible.
//!
//! A `Linestring` must be a model of `RandomAccessRange` of `CGAL::Point_2`.
//! \see CGAL::Point_2
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

//! \ingroup PkgIOstreams
//! \brief read_linestring_WKT reads the content of a .wkt file into 
//! a `MultiLineString` if possible.
//!
//! A `MultiLineString` must be a model of `RandomAccessRange` of `Linestring`.
//! 
//! \see CGAL::Point_2
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
  BOOST_FOREACH(LineString& ls, gc)
  {
    mls.push_back(*ls.range);
  }
  return in;  
}

//! \ingroup PkgIOstreams
//! \brief read_polygon_WKT reads the content of a .wkt file into a `Polygon` if possible.
//! 
//! A `Polygon` must inherit `CGAL::General_polygon_with_holes_2`.
//! 
//! \see CGAL::General_polygon_with_holes_2
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
      boost::geometry::read_wkt(line, polygon);
      break;
    }
  }
  return in;  
}

//! \ingroup PkgIOstreams
//! \brief read_multipolygon_WKT reads the content of a .wkt file into 
//! a `Multipolygon` if possible.
//!
//! A `Multipolygon` must be a model of `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`.
//! \see CGAL::General_polygon_with_holes_2

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
      boost::geometry::read_wkt(line, gc);
      break;
    }
  }
  return in;  
}

//! \ingroup PkgIOstreams
//! \brief write_point_WKT writes `point` into a .WKT file.
//! `Point` must be a `CGAL::Point_2`
//! \see CGAL::Point_2
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

//! \ingroup PkgIOstreams
//! \brief write_point_WKT writes `poly` into a .WKT file.
//! `Polygon` must be a `CGAL::General_polygon_with_holes_2`
//! \see `CGAL::General_polygon_with_holes_2`
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

//! \ingroup PkgIOstreams
//! \brief write_linestring_WKT writes the content of `ls` 
//! into a WKT file.
//! `LineString` must be a `RandomAccessRange` of `CGAL::Point_2`. 
//!\see `CGAL::Point_2`
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

//! \ingroup PkgIOstreams
//! \brief write_multi_point_WKT writes the content of `mp` 
//! into a WKT file.
//! `MultiPoint` must be a `RandomAccessRange` of `CGAL::Point_2`. 
//!\see `CGAL::Point_2`
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

//! \ingroup PkgIOstreams
//! \brief write_multi_polygon_WKT writes the content of `polygons` 
//! into a WKT file.
//! `multiPolygon` must be a `RandomAccessRange` of `CGAL::General_polygon_with_holes_2`. 
//!\see `CGAL::General_polygon_with_holes_2`
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

//! \ingroup PkgIOstreams
//! \brief write_multi_linestring_WKT writes the content of `mls` 
//! into a WKT file.
//! `MultiLineString` must be a `RandomAccessRange` of `LineString`. 
//! \see `CGAL::write_linestring_WKT`
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
  BOOST_FOREACH(PointRange& pr, mls)
  {
    LineString ls(pr);
    pr_range.push_back(ls);
  }
  internal::Geometry_container<std::vector<LineString>, boost::geometry::multi_linestring_tag> gc(pr_range);
  out<<boost::geometry::wkt(gc)<<std::endl;
  return out;
}
}//end CGAL
#endif
