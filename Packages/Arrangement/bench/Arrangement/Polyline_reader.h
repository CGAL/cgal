#ifndef POLYLINE_READER_H
#define POLYLINE_READER_H

#include <CGAL/basic.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <list>

#include "numberType.h"

template <class Traits>
class Polyline_reader
{

public:
  typedef typename Traits::Point_2    Point_2;
  typedef typename Traits::Curve_2    Curve_2;
  typedef typename Traits::X_monotone_curve_2  X_monotone_curve_2;
  typedef std::list<Curve_2> CurveList;

  int ReadData(const char * filename, CurveList & curves,
               CGAL::Bench_parse_args::FormatId format,
               CGAL::Bbox_2 & bbox)
  {
    std::ifstream file(filename);
    curves.clear();

    int num_polylines, num_segments;
    int ix, iy;
    std::vector<Point_2> points;
    int i, j;

    file >> num_polylines;
    for (i = 0; i < num_polylines; i++) 
    {
      file >> num_segments;
      points.clear();
      for (j = 0; j < num_segments; j++)
      {
        file >> ix >> iy;
        points.push_back (Point_2(NT(ix),NT(iy)));
      }

      Curve_2   polyline(points);
      curves.push_back(polyline);
      if (i == 0)
        bbox = polyline.bbox();
      else
        bbox = bbox + polyline.bbox();
    }
    
    return 0;
  }
};

#endif
