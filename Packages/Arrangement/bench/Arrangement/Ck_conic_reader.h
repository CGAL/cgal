#ifndef CK_CONIC_READER_H
#define CK_CONIC_READER_H

#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>

#include "numberType.h"

template <class Traits>
class Ck_conix_reader {
public:
  typedef typename Traits::Point_2               Point_2;
  typedef typename Traits::Curve_2               Curve_2;
  typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;
  typedef std::list<Curve_2>                     Curve_list;
    
  int read_data(const char* filename, Curve_list& curves,
                CGAL::Bench_parse_args::FormatId format,
                CGAL::Bbox_2 & bbox)
  {
    return 0;
  }
};

#endif
