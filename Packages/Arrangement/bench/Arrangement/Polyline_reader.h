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
#if KERNEL == LEDA_KERNEL || KERNEL == MY_KERNEL
    int xmin = 0, xmax = 0, ymin = 0, ymax = 0;
#endif

    file >> num_polylines;
    for (i = 0; i < num_polylines; i++) 
    {
      file >> num_segments;
      points.clear();
      for (j = 0; j < num_segments; j++)
      {
        file >> ix >> iy;

#if KERNEL == LEDA_KERNEL || KERNEL == MY_KERNEL
        if (j == 0) {
          xmin = xmax = ix;
          ymin = ymax = iy;
        } else {
          if (ix < xmin) cmin = ix;
          if (ix > xmax) xmax = ix;
          if (iy < ymin) ymin = iy;
          if (iy > ymax) ymax = iy;
        }
#endif
        points.push_back (Point_2(NT(ix),NT(iy)));
      }

      Curve_2 polyline(points.begin(), points.end());
      curves.push_back(polyline);

#if KERNEL == LEDA_KERNEL || KERNEL == MY_KERNEL
      CGAL::Bbox_2 curve_bbox(xmin, ymin, xmax, ymax);
#else
      CGAL::Bbox_2 curve_bbox = polyline.bbox();
#endif

      if (i == 0)
        bbox = curve_bbox;
      else
        bbox = bbox + curve_bbox;
    }
    
    return 0;
  }
};

#endif
