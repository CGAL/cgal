#ifndef CK_CONIC_READER_H
#define CK_CONIC_READER_H

#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <list>

#include "numberType.h"

template <class Traits>
class Ck_conic_reader {
public:
  typedef typename Traits::Kernel                       Curved_k;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::Point_2                      Point_2;
  
  template<class OutputIterator>
  int read_data(const char* filename, OutputIterator curves_out,
                CGAL::Bench_parse_args::FormatId format,
                CGAL::Bbox_2 & bbox)
  {
    std::ifstream inp(filename);
    if (!inp.is_open()) {
      std::cerr << "Cannot open file " << filename << "!" << std::endl;
      return -1;
    }

    std::vector<int> C;
    std::copy(std::istream_iterator<int>(inp),
              std::istream_iterator<int>(), std::back_inserter(C));  
    for(size_t i = 0; i < C.size(); i += 6) {
      Curve_2 cv;
      add_conic(cv, C[i], C[i+1], C[i+2], C[i+3], C[i+4], C[i+5]);
      ++curves_out = cv;
#if 0
      CGAL::Bbox_2 curve_bbox = cv.bbox();
      if (i == 0) bbox = curve_bbox;
      else bbox = bbox + curve_bbox;
#else
      bbox = CGAL::Bbox_2(-7, -7, 28, 29);
#endif
    }
    inp.close();
    return 0;
  }

  /*! */
  void add_conic(Curve_2 & cv, int r=0,int s=0,int t=0,int u=0,int v=0,int w=0)
  {
    FT gc[6] = {w,v,u,t,s,r};
    cv =
      ECG::Conic_arc_2<Curved_k>
      (ECG::Conic_2<Curved_k>(gc[5],gc[4],gc[3],gc[2],gc[1],gc[0]));
  }
};

#endif
