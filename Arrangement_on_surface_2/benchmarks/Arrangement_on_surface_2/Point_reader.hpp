#ifndef POINT_READER_HPP
#define POINT_READER_HPP

#include <CGAL/basic.h>
#include <CGAL/Bbox_2.h>
#include <iostream>
#include <fstream>
#include <list>

#include "number_type.hpp"
#include "Number_type_traits.hpp"
#include "Option_parser.hpp"

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
      Number_type x, y;
      if (format == Option_parser::FORMAT_RATIONAL) {
        Input_traits<Number_type>::Input_rat_type ix, iy;
        inp >> ix >> iy;
        x = ix; y = iy;
      } else if (format == Option_parser::FORMAT_INT) {
        Input_traits<Number_type>::Input_int_type ix, iy;
        inp >> ix >> iy;
        x = (Number_type) ix; y = (Number_type) iy; 
      } else if (format == Option_parser::FORMAT_FLOAT) {
        Input_traits<Number_type>::Input_float_type ix, iy;
        inp >> ix >> iy;
        x = (Number_type) ix; y = (Number_type) iy; 
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
