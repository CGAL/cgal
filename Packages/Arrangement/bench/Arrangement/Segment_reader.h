#ifndef SEGMENT_READER_H
#define SEGMENT_READER_H

#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <list>

#include "numberType.h"

template <class Traits>
class Segment_reader
{

public:
  typedef typename Traits::Point_2    Point_2;
  typedef typename Traits::Curve_2    Curve_2;
  typedef typename Traits::X_curve_2  X_curve_2;
  typedef std::list<Curve_2> CurveList;

  int ReadData(const char * filename, CurveList & curveList,
               CGAL::Bench_parse_args::FormatId format,
               CGAL::Bbox_2 & bbox)
  {
    std::ifstream inp(filename);
    if (!inp.is_open()) {
      std::cerr << "Cannot open file " << filename << "!" << std::endl;
      return -1;
    }
    int count;
    inp >> count;
    
    int i;
    for (i = 0; i < count; i++) {
      NT x0, y0, x1, y1;
      if (format == CGAL::Bench_parse_args::FORMAT_RAT) {
        inp >> x0 >> y0 >> x1 >> y1;
      } else if (format == CGAL::Bench_parse_args::FORMAT_INT) {
        int ix0, iy0, ix1, iy1;
        inp >> ix0 >> iy0 >> ix1 >> iy1;
        x0 = ix0; y0 = iy0; x1 = ix1; y1 = iy1;
      } else if (format == CGAL::Bench_parse_args::FORMAT_FLT) {
        float ix0, iy0, ix1, iy1;
        inp >> ix0 >> iy0 >> ix1 >> iy1;
        x0 = ix0; y0 = iy0; x1 = ix1; y1 = iy1;
      } else {
        std::cerr << "Illegal format!" << std::endl;
        return -1;
      }

#if defined(USE_LAZY_RAT) || defined(USE_LAZY_QUOTIENT)
      WNT lazy_exact_x0(x0);
      WNT lazy_exact_x1(x1);
      WNT lazy_exact_y0(y0);
      WNT lazy_exact_y1(y1);
      Point_2 p1(lazy_exact_x0, lazy_exact_y0);
      Point_2 p2(lazy_exact_x1, lazy_exact_y1);
#else
      Point_2 p1(x0, y0);
      Point_2 p2(x1, y1);
#endif
      // if (p1 == p2) continue;
      Curve_2 curve(p1, p2);
      curveList.push_back(curve);

      // Update the bounding box of the arrangement.
#if defined(USE_LEDA_KERNEL) || defined(USE_MY_KERNEL)
      double xmin, ymin, xmax, ymax;
      if (p1.xcoord() < p2.xcoord()) {
        xmin = CGAL::to_double(p1.xcoord());
        xmax = CGAL::to_double(p2.xcoord());
      } else {
        xmin = CGAL::to_double(p2.xcoord());
        xmax = CGAL::to_double(p1.xcoord());
      }
      if (p1.ycoord() < p2.ycoord()) {
        ymin = CGAL::to_double(p1.ycoord());
        ymax = CGAL::to_double(p2.ycoord());
      } else {
        ymin = CGAL::to_double(p2.ycoord());
        ymax = CGAL::to_double(p1.ycoord());
      }
      
      CGAL::Bbox_2 curve_bbox(xmin, ymin, xmax, ymax);
#else
      CGAL::Bbox_2 curve_bbox = curve.bbox();
#endif
      if (i == 0) bbox = curve_bbox;
      else bbox = bbox + curve_bbox;
    }
    inp.close();
    return 0;
  }
};

#endif
