#ifndef EXACUS_CONIC_READER_H
#define EXACUS_CONIC_READER_H

#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>

#include "numberType.h"

#include <LiS/file_io.h>


template <class Traits>
class Exacus_conic_reader {
    
public:
    typedef typename Traits::Point_2               Point_2;
    typedef typename Traits::Curve_2               Curve_2;
    typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;
    
    template<class OutputIterator>
    int read_data(const char* filename, OutputIterator curves_out,
                  CGAL::Bench_parse_args::FormatId format,
                  CGAL::Bbox_2 & bbox) {
        
        Curve_2 cv;

        std::list< Curve_2 > curves;
        bool success = LiS::read_file(filename, curves);
        if (!success) {
            return 0;
        }
        std::copy(curves.begin(), curves.end(), curves_out);
        
#if 0
        // TODO set these boxes
        for (typename std::list< Curve_2 >::iterator it = curves.begin();
             it != curves.end();
             it++) {
            CGAL::Bbox_2 curve_bbox = cv.bounding_box();
            if (i == 0) {
                bbox = curve_bbox;
            } else {
                bbox = bbox + curve_bbox;
            }
        }
#else
        bbox = CGAL::Bbox_2(-1000, -1000, 1000, 1000);
#endif
        return static_cast< int >(curves.size());
    }
};

#endif // EXACUS_CONIC_READER
