#ifndef SEGMENT_READER_H
#define SEGMENT_READER_H

#include <CGAL/basic.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <list>

#include "numberType.h"
#include "Input_traits.h"

template <class Traits>
class Segment_reader
{

public:
  typedef typename Traits::Point_2    Point_2;
  typedef typename Traits::Curve_2    Curve_2;
  typedef typename Traits::X_monotone_curve_2  X_monotone_curve_2;

  template<class OutputIterator>
  int read_data(const char * filename, OutputIterator curves_out,
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
      WNT x0, y0, x1, y1;
      if (format == CGAL::Bench_parse_args::FORMAT_RAT) {
        Input_traits<WNT>::Input_rat_type ix0, iy0, ix1, iy1;
        inp >> ix0 >> iy0 >> ix1 >> iy1;
        x0 = ix0; y0 = iy0; x1 = ix1; y1 = iy1;
      } else if (format == CGAL::Bench_parse_args::FORMAT_INT) {
        Input_traits<WNT>::Input_int_type ix0, iy0, ix1, iy1;
        inp >> ix0 >> iy0 >> ix1 >> iy1;
        x0 = (WNT) ix0; y0 = (WNT) iy0; x1 = (WNT) ix1; y1 = (WNT) iy1;
      } else if (format == CGAL::Bench_parse_args::FORMAT_FLT) {
        Input_traits<WNT>::Input_float_type ix0, iy0, ix1, iy1;
        inp >> ix0 >> iy0 >> ix1 >> iy1;
        x0 = (WNT) ix0; y0 = (WNT) iy0; x1 = (WNT) ix1; y1 = (WNT) iy1;
      } else {
        std::cerr << "Illegal format!" << std::endl;
        return -1;
      }
      
      Point_2 p1(x0, y0);
      Point_2 p2(x1, y1);

      // if (p1 == p2) continue;
      Curve_2 curve(p1, p2);
      ++curves_out = curve;

      // Update the bounding box of the arrangement.
#if KERNEL == LEDA_KERNEL || KERNEL == MY_KERNEL
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
