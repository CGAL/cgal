// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion
//                 Andreas Fabri, Ron Wein, Julien Hazebrouck

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

// Descriptions of the file format can be found at
// http://www.autodesk.com/techpubs/autocad/acad2000/dxf/
// http://www.tnt.uni-hannover.de/soft/compgraph/fileformats/docs/DXF.ascii

#ifndef CGAL_IO_DXF_VARIANT_READER_H
#define CGAL_IO_DXF_VARIANT_READER_H

#include <CGAL/IO/Dxf_reader_doubles.h>
#include <iostream>
#include <string>
#include <list>
#include <boost/variant.hpp>
#include <CGAL/array.h>


namespace CGAL {

template<class CK,class Circular_arc_2, class Line_arc_2, class OutputIterator>
  OutputIterator variant_load(std::istream& is, OutputIterator res)
{

  typedef cpp11::array<double, 3> Triplet;
  typedef typename CK::FT FT;
  typedef typename CK::Circular_arc_point_2 Circular_arc_point_2;    
  typedef typename CK::Root_of_2 Root_of_2;
  typedef typename CK::Root_for_circles_2_2 Root_for_circles_2_2;
  typedef typename CK::Line_2  Line_2;
  typedef typename CK::Point_2 Point_2;
  typedef typename CK::Circle_2 Circle_2;
  typedef typename boost::variant< Circular_arc_2, Line_arc_2 >        Arc;
  typedef std::list<Triplet> Polygon;
  typedef std::list<Polygon> Polygons;
  typedef std::list<Triplet> Circles;

  Polygons polygons;
  Circles circles;
  CGAL::Dxf_reader_doubles reader;
  
  reader(is, polygons, circles);

  std::cout << "Read " << polygons.size() << " polygons, and " 
	    << circles.size() << " circles" << std::endl;


  for(typename Circles::iterator it = circles.begin(); it != circles.end(); it++){
    Arc arc = typename CK::Construct_circular_arc_2()(typename CK::Construct_circle_2()(typename CK::Construct_point_2()((*it)[0], (*it)[1]), FT((*it)[2])));
    *res++ = arc;
  }
  
  std::map<std::pair<double,double>, Circular_arc_point_2> points;
  typename std::map<std::pair<double,double>, Circular_arc_point_2>::iterator p_cap_it;

  double bulge;

  Circular_arc_point_2 caps, capt;
  Arc arc;

  for(typename Polygons::iterator it = polygons.begin(); it != polygons.end(); it++){
    typename Polygon::iterator pit = it->begin();

    std::pair<double,double> xyfirst = std::make_pair((*pit)[0], (*pit)[1]);    
    std::pair<double,double> xyps, xypt = std::make_pair((*pit)[0], (*pit)[1]);
    Point_2 ps, pt = typename CK::Construct_point_2()(xypt.first, xypt.second);
    Point_2 first = pt;

      while(true){
      xyps = xypt;
      ps = pt;
      bulge = (*pit)[2];
      pit++;

      if(pit ==it->end()){
	break;
      }
      xypt = std::make_pair((*pit)[0], (*pit)[1]);
      pt = typename CK::Construct_point_2()(xypt.first, xypt.second);

      p_cap_it = points.find(xyps);	  
      if(p_cap_it == points.end()){
	caps = typename CK::Construct_circular_arc_point_2()(ps);
	points.insert(std::make_pair(xyps, caps));
      }else{
	caps = p_cap_it->second;
      }
      p_cap_it = points.find(xypt);	  
      if(p_cap_it == points.end()){
	capt = typename CK::Construct_circular_arc_point_2()(pt);
	points.insert(std::make_pair(xypt, capt));
      } else {
	capt = p_cap_it->second;
      } 
      
      if(bulge == 0){
	
	if(xyps != xypt){
	  typename CK::Line_2 l(ps,pt);
	  Line_arc_2 la  = typename CK::Construct_line_arc_2()(l,caps, capt);
	  arc = la;
	  *res++ = arc;
	}
      } else {
	Circular_arc_2 carc = typename CK::Construct_circular_arc_2()(typename CK::Construct_circle_2()(ps, pt, bulge),
								      caps, capt);
	arc = carc;
	*res++ = arc;
      }
    }

    
    if(bulge == 0){
      if(xypt != xyfirst){
	arc = typename CK::Construct_line_arc_2()(typename CK::Construct_line_2()(pt, first),capt, points.find(xyfirst)->second);      
	*res++ = arc;
      }
    } else {
      Circular_arc_2 carc = typename CK::Construct_circular_arc_2()(typename CK::Construct_circle_2()(pt, first, bulge),
								    capt, points.find(xyfirst)->second);
      arc = carc;
      *res++ = arc;
    }
    
  }
  std::cout << " Loaded" << std::endl;
  
  return res;
}

}// namespace CGAL

#endif // CGAL_IO_DXF_VARIANT_READER_H
