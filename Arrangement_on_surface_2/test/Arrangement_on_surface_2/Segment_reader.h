#ifndef SEGMENT_READER_H
#define SEGMENT_READER_H

#include <CGAL/Bbox_2.h>
#include <iostream>
#include <fstream>
#include <list>

template <class Traits>
class Segment_reader {
public:
  typedef typename Traits::Kernel     Kernel;
  typedef typename Kernel::FT         NT;
  typedef typename Traits::Point_2    Point_2;
  typedef typename Traits::Curve_2    Curve_2;
  typedef typename Traits::X_monotone_curve_2  X_monotone_curve_2;

  template<class OutputIterator>
  int read_data(const char* filename, OutputIterator curves_out,
                CGAL::Bbox_2& bbox)
  {
    std::ifstream inp(filename);
    if (!inp.is_open()) {
      std::cerr << "Cannot open file " << filename << "!" << std::endl;
      return -1;
    }
    int count = 0;
    inp >> count;
        //std::cout << "count ="<<count <<std::endl;

    NT    x1, y1, x2, y2;
    int   i;

    for (i = 0; i < count; i++) {
      inp >> x1 >> y1 >> x2 >> y2;
          //std::cout << "x1 ="<<x1<< "x2 ="<<x2<< "y1 ="<<y1<< "y2 ="<<y2<<std::endl;

      Point_2 p1(x1, y1);
      Point_2 p2(x2, y2);

      Curve_2 curve(p1, p2);
      ++curves_out = curve;

      // Update the bounding box of the arrangement.
      CGAL::Bbox_2 curve_bbox = curve.bbox();

      if (i == 0)
        bbox = curve_bbox;
      else
        bbox = bbox + curve_bbox;
    }
    inp.close();
    return 0;
  }
};

#endif
