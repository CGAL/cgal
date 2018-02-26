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

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Point_2.h>

#include <CGAL/IO/traits_point.h>
#include <CGAL/IO/traits_polygon.h>



namespace CGAL{

template<typename K,
         typename OutputPointIterator,
         typename OutputPolygonIterator>
bool
read_WKT( std::istream& in,
          OutputPointIterator out_point,
          OutputPolygonIterator out_polygons
          )
{
  typedef CGAL::Point_2<K> Point;
  typedef CGAL::Polygon_with_holes_2<K> Polygon;
  if(!in)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }
  std::string line;
  while(std::getline(in, line))
  {
    std::istringstream iss(line);
    std::string type;
    iss >> type;
    
    if(type.compare("POINT")==0)
    {
      Point p;
      boost::geometry::read_wkt(line, p);
      *out_point++ = p;
    }
    else if(type.compare("POLYGON")==0)
    {
      Polygon p;
      boost::geometry::read_wkt(line, p);
      *out_polygons = p;
    }
    else if(type.length()>1){
      std::cout<<"unknown WKT type: "<<type<<std::endl;
    }
  }
  return in.good();  
}


//!
//! \brief read_polygon_WKT
//! \tparam Polygon 
//! 
//! \param in
//! \param polygon
//! \return 
//! 
//! \relates CGAL::Polygon_with_holes_2
//!
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

template<typename PolygonRange>
bool
write_polygons_WKT( std::ostream& out,
          const PolygonRange& polygons
          )
{
  typedef typename boost::range_value<PolygonRange>::type Polygon;
  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }
  for(typename PolygonRange::const_iterator it = polygons.begin();
      it != polygons.end();
      ++it)
  out<<boost::geometry::wkt(*it)<<std::endl;
  return out.good();
}

template<typename Polygon>
std::ostream&
write_polygon_WKT( std::ostream& out,
          const Polygon& polygon
          )
{
  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return out;
  }
  out<<boost::geometry::wkt(polygon)<<std::endl;
  return out;
}

}//end CGAL
#endif