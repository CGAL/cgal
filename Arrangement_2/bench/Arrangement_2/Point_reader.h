#ifndef POINT_READER_H
#define POINT_READER_H

#include <CGAL/basic.h>
#include <CGAL/Bbox_2.h>
#include <iostream>
#include <fstream>
#include <list>

#include "number_type.h"
#include "Input_traits.h"
#include "Option_parser.h"

template <class Traits>
class Point_reader {
public:
  typedef typename Traits::Point_2    Point_2;

  template<class OutputIterator>
  int read_data(const char * filename, OutputIterator points_out,
                Option_parser::Format_id format)  {
    std::ifstream inp(filename);
    if (!inp.is_open()) {
      std::cerr << "Cannot open file " << filename << "!" << std::endl;
      return -1;
    }
    int count;
    inp >> count;
    
    int i;
    for (i = 0; i < count; ++i) {
      WNT x, y;
      if (format == Option_parser::FORMAT_RATIONAL) {
        Input_traits<WNT>::Input_rat_type ix, iy;
        inp >> ix >> iy;
        x = ix; y = iy;
      } else if (format == Option_parser::FORMAT_INT) {
        Input_traits<WNT>::Input_int_type ix, iy;
        inp >> ix >> iy;
        x = (WNT) ix; y = (WNT) iy; 
      } else if (format == Option_parser::FORMAT_FLOAT) {
        Input_traits<WNT>::Input_float_type ix, iy;
        inp >> ix >> iy;
        x = (WNT) ix; y = (WNT) iy; 
      } else {
        std::cerr << "Illegal format!" << std::endl;
        return -1;
      }
      
      Point_2 pnt(x, y);
      ++points_out = pnt;
    }

    inp.close();
    return 0;
  }
};

#endif
