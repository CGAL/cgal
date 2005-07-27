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

/*! Convert WNT to int. Used to calculate the bounding box */
template <class T_NT>
class Toint {
public:
  int operator()(T_NT & n) { return static_cast<int>(n); }
};

/*! Convert leda_rational to int. Used to calculate the bounding box */
#if BENCH_NT == LEDA_RAT_NT
template <>
class Toint<leda_rational> {
public:
  int operator()(leda_rational & n)
  {
    if (n.denominator() == n.numerator()) return 1;
    leda_integer lint = n.denominator() / n.denominator();
    return lint.to_long();
  }
};
#endif

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

#if BENCH_KERNEL == LEDA_KERNEL || BENCH_KERNEL == MY_KERNEL
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

#if BENCH_KERNEL == LEDA_KERNEL || BENCH_KERNEL == MY_KERNEL
        Toint<WNT> toint;
        if (j == 0) {
          xmin = xmax = toint(x);
          ymin = ymax = toint(y);
        } else {
          if (x < xmin) xmin = toint(x);
          if (x > xmax) xmax = toint(x);
          if (y < ymin) ymin = toint(y);
          if (y > ymax) ymax = toint(y);
        }
#endif
      }

      Curve_2 polyline(points.begin(), points.end());
      ++curves_out = polyline;

#if BENCH_KERNEL == LEDA_KERNEL || BENCH_KERNEL == MY_KERNEL
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
