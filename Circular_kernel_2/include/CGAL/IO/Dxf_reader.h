// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Andreas Fabri

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

// Description of the file format can be found at the following address:
// http://www.autodesk.com/techpubs/autocad/acad2000/dxf/

#ifndef CGAL_IO_DXF_READER_H
#define CGAL_IO_DXF_READER_H

#include <CGAL/license/Circular_kernel_2.h>


#include <CGAL/basic.h>
#include <iostream>
#include <string>
#include <list>

namespace CGAL {

template <typename K>
class Dxf_reader {

public:
  typedef typename K::FT       FT;
  typedef typename K::Point_2  Point_2;
  typedef typename K::Circle_2 Circle_2;


  typedef std::list<std::pair<Point_2, double> > Polygon;
  typedef std::list<Polygon> Polygons;
  typedef std::list<Circle_2> Circles;
  typedef std::list<std::pair<Point_2, FT> > Centers_and_radii;

private:

  void
  header(std::istream& is)
  {
    int n;
    double xmin, ymin;
    double xmax, ymax;
    is >> n;
    CGAL_assertion(n == 9);
    char c;
    is >> c;
    CGAL_assertion(c == '$');
    std::string str;
    is >> str;
    if(str == std::string("EXTMIN")){
      is >> n;
      CGAL_assertion(n == 10);
    is >> xmin;
    is >> n;
    CGAL_assertion(n == 20);
    is >> ymin;
    }
    is >> n;
    CGAL_assertion(n == 9);
    is >> c;
    CGAL_assertion(c == '$');
    is >> str;
    if(str == "EXTMAX"){
      is >> n;
      CGAL_assertion(n == 10);
      is >> xmax;
      is >> n;
      CGAL_assertion(n == 20);
      is >> ymax;
    }
  }


  void
  skip_header(std::istream& is)
  {
    int n;
    is >> n;
    CGAL_assertion(n == 0);
    std::string str;
    is >> str;
    CGAL_assertion(str == "SECTION");
    is >> n;
    CGAL_assertion(n == 2);
    is >> str;
    if(str == "HEADER"){
      header(is);
    }
    is >> n;
    CGAL_assertion(n == 0);
    is >> str;
    CGAL_assertion(str == "ENDSEC");
  }



  void
  read_circle(std::istream& is, Circle_2& circ)
  {
    int n;
    double cx, cy, r;
    std::string str;
    is >> n;
    CGAL_assertion(n == 8);
    is >> n;
    CGAL_assertion(n == 0);

  is >> n;
  CGAL_assertion(n == 10);
  is >> IO::iformat(cx);
  is >> n;
  CGAL_assertion(n == 20);
  is >> IO::iformat(cy);
  is >> n;
  CGAL_assertion(n == 40);
  is >> IO::iformat(r);
  FT sqr_ft(r*r);
  circ = typename K::Construct_circle_2()(Point_2(cx,cy), sqr_ft);
}

  void
  read_center_and_radius(std::istream& is, Point_2& center, FT& rft)
  {
    int n;
    double cx, cy, r;
    std::string str;
    is >> n;
    CGAL_assertion(n == 8);
    is >> n;
    CGAL_assertion(n == 0);

  is >> n;
  CGAL_assertion(n == 10);
  is >> IO::iformat(cx);
  is >> n;
  CGAL_assertion(n == 20);
  is >> IO::iformat(cy);
  is >> n;
  CGAL_assertion(n == 40);
  is >> IO::iformat(r);

  center = typename K::Construct_point_2()(cx,cy);
  rft = FT(r); // intentionally not squared
}


void
read_polygon(std::istream& is, Polygon& poly)
{
  int n = 0;
  do {
    is >> n;
    if(n != 0){
      int m;
      is >> m;
    }
  } while(n != 0);

  std::string str;
  do {
    double len;
    double x, y;
    is >> str;
    if(str == "VERTEX"){
      is >> n;
      CGAL_assertion(n == 8);
      is >> n;
      CGAL_assertion(n == 0);
      is >> n;
      CGAL_assertion(n == 10);
      is >> IO::iformat(x);
      is >> n;
      CGAL_assertion(n == 20);
      is >> IO::iformat(y);
      is >> n;
      len = 0;
      if(n == 42){
        is >> len;
      } else {
        CGAL_assertion(n == 0);
      }
      poly.push_back(std::make_pair(typename K::Construct_point_2()(x,y), len));
    }

  } while (str != "SEQEND");
  is >> n;
  CGAL_assertion(n == 8);
  is >> n;
  CGAL_assertion(n == 0);


}


void
read_entities(std::istream& is, Polygons& polys, Circles& circles)
{
  int n;
  //double x, y;
  std::string str;
  is >> n;
  CGAL_assertion(n == 0);
  is >> str;
  CGAL_assertion(str == "SECTION");
  is >> n;
  is >> str;
  CGAL_assertion(str == "ENTITIES");
  do {
    is >> n;
    CGAL_assertion(n == 0);
    is >> str;
    if(str == "POLYLINE"){
      Polygon p;
      polys.push_back(p);
      read_polygon(is, polys.back());
    } else if(str == "CIRCLE"){
      Circle_2 c;
      read_circle(is,c);
      circles.push_back(c);
    } else if(str == "ENDSEC"){

    } else {
      CGAL_error_msg( "unknown entity" );
    }
  } while(str != "ENDSEC");
  is >> n;
  CGAL_assertion(n == 0);
  is >> str;
  CGAL_assertion(str == "EOF");
}

void
read_entities(std::istream& is, Polygons& polys, Centers_and_radii& car)
{
  int n;
  std::string str;
  is >> n;
  CGAL_assertion(n == 0);
  is >> str;
  CGAL_assertion(str == "SECTION");
  is >> n;
  is >> str;
  CGAL_assertion(str == "ENTITIES");
  do {
    is >> n;
    CGAL_assertion(n == 0);
    is >> str;
    if(str == "POLYLINE"){
      Polygon p;
      polys.push_back(p);
      read_polygon(is, polys.back());
    } else if(str == "CIRCLE"){
      Point_2 center;
      FT radius;
      read_center_and_radius(is,center, radius);
      car.push_back(std::make_pair(center, radius));
    } else if(str == "ENDSEC"){

    } else {
      CGAL_error_msg( "unknown entity" );
    }
  } while(str != "ENDSEC");
  is >> n;
  CGAL_assertion(n == 0);
  is >> str;
  CGAL_assertion(str == "EOF");
}

public:

void operator()(std::istream& is, Polygons& polygons, Circles& circles)
{
  skip_header(is);
  read_entities(is, polygons, circles);
}

void operator()(std::istream& is, Polygons& polygons, Centers_and_radii& car)
{
  skip_header(is);
  read_entities(is, polygons, car);
}
};

} //namespace CGAL

#endif // CGAL_IO_DXF_READER_H
