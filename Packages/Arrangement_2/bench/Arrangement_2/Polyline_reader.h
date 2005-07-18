#ifndef POLYLINE_READER_H
#define POLYLINE_READER_H

#include <CGAL/basic.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <list>

#include "number_type.h"
#include "Input_traits.h"

template <class Traits>
class Polyline_reader {
public:
  typedef typename Traits::Point_2    Point_2;
  typedef typename Traits::Curve_2    Curve_2;
  typedef typename Traits::X_monotone_curve_2  X_monotone_curve_2;

  template<class OutputIterator>
  int read_data(const char * filename, OutputIterator curves_out,
               CGAL::Bench_parse_args::FormatId format,
               CGAL::Bbox_2 & bbox)
  {
    std::ifstream file(filename);

#if KERNEL == LEDA_KERNEL || KERNEL == MY_KERNEL
    int xmin = 0, xmax = 0, ymin = 0, ymax = 0;
#endif

    unsigned int num_polylines;
    file >> num_polylines;
    for (unsigned int i = 0; i < num_polylines; i++) {
      unsigned int num_points;
      file >> num_points;
      std::vector<Point_2> points;
      points.clear();
      for (unsigned int j = 0; j < num_points; j++) {
        WNT x, y;
        if (format == CGAL::Bench_parse_args::FORMAT_RAT) {
          Input_traits<WNT>::Input_rat_type ix, iy;
          file >> ix >> iy;
          x = ix; y = iy;
        } else if (format == CGAL::Bench_parse_args::FORMAT_INT) {
          Input_traits<WNT>::Input_int_type ix, iy;
          file >> ix >> iy;
          x = (WNT) ix; y = (WNT) iy;
        } else if (format == CGAL::Bench_parse_args::FORMAT_FLT) {
          Input_traits<WNT>::Input_float_type ix, iy;
          file >> ix >> iy;
          x = (WNT) ix; y = (WNT) iy;
        } else {
          std::cerr << "Illegal format!" << std::endl;
          return -1;
        }
        Point_2 p(x, y);
        points.push_back(p);

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
      }

      Curve_2 polyline(points.begin(), points.end());
      ++curves_out = polyline;

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
