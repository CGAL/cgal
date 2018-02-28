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

#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <CGAL/IO/traits_point.h>
#include <CGAL/IO/traits_linestring.h>
#include <CGAL/IO/traits_polygon.h>
#include <CGAL/IO/traits_multipoint.h>
#include <CGAL/IO/traits_multilinestring.h>
#include <CGAL/IO/traits_multipolygon.h>



namespace CGAL{

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
    
    if(type.compare("POINT")==0)
    {
      boost::geometry::read_wkt(line, point);
      break;
    }
  }
  return in;  
}

template<typename Multipoint>
std::istream&
read_multipoint_WKT( std::istream& in,
                     Multipoint& multipoint)
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
    
    if(type.compare("MULTIPOINT")==0)
    {
      boost::geometry::read_wkt(line, multipoint);
      break;
    }
  }
  return in;  
}


//! \ingroup PkgPolygonIO
//! \brief read_linestring_WKT reads the content of a .wkt file into 
//! an `std::vector<CGAL::Point_2>`.
//! 
//! \relates CGAL::Point_2
template<typename Linestring>
std::istream&
read_linestring_WKT( std::istream& in,
                     Linestring& polyline)
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
    
    if(type.compare("LINESTRING")==0)
    {
      boost::geometry::read_wkt(line, polyline);
      break;
    }
  }
  return in;  
}

template<typename Multilinestring>
std::istream&
read_multilinestring_WKT( std::istream& in,
                          Multilinestring& mls)
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
    
    if(type.compare("MULTILINESTRING")==0)
    {
      boost::geometry::read_wkt(line, mls);
      break;
    }
  }
  return in;  
}

//! \ingroup PkgPolygonIO
//! \brief read_polygon_WKT reads the content of a .wkt file into a `Polygon`.
//! 
//! A `Polygon` must inherit `CGAL::General_polygon_with_holes_2`.
//! 
//! \relates CGAL::General_polygon_with_holes_2
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
    
    if(type.compare("POLYGON")==0)
    {
      boost::geometry::read_wkt(line, polygon);
      break;
    }
  }
  return in;  
}

template<typename MultiPolygon>
std::istream&
read_multipolygon_WKT( std::istream& in,
                       MultiPolygon& polygons
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
    
    if(type.compare("MULTIPOLYGON")==0)
    {
      boost::geometry::read_wkt(line, polygons);
      break;
    }
  }
  return in;  
}
//! \ingroup PkgPolygonIO
//! \brief write_WKT writes the content of `type` into a .WKT file.
//! `Type` must have appropriate boost::geometry::traits available.
//! Such traits are already available in CGAL for the following type :
//! - `CGAL::Point_2` as WKT Point
//! - `std::list<CGAL::Point_2>` as WKT Linestring
//! - `CGAL::Polygon_with_hole_2` as WKT Polygon
//! - `std::vector<CGAL::Point_2>` as WKT Multipoint
//! - `std::vector<std::list<CGAL::Point_2>>` as WKT Multilinestring
//! - `std::vector<CGAL::Polygon_with_hole_2>` as WKT Multipolygon
//! 
//! \relates CGAL::General_polygon_with_holes_2
//! \todo We should find a way to generalize the range types instead of 
//! forcing list and vector. Maybe through a pair of <Range, Tag> to avoid 
//! conflicts between traits.
template<typename Type>
std::ostream&
write_WKT( std::ostream& out,
           const Type& type
           )
{
  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return out;
  }
  out<<boost::geometry::wkt(type)<<std::endl;
  return out;
}


}//end CGAL
#endif