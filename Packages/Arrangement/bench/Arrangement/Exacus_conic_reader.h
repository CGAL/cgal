#ifndef EXACUS_CONIX_READER_H
#define EXACUS_CONIX_READER_H

#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>

#include "numberType.h"

template <class Traits>
class Exacus_conix_reader {
    
public:
    typedef typename Traits::Point_2               Point_2;
    typedef typename Traits::Curve_2               Curve_2;
    typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;
    typedef std::list<Curve_2>                     CurveList;
    
    int ReadData(const char* filename, CurveList& curves,
                 CGAL::Bench_parse_args::FormatId format,
                 CGAL::Bbox_2 & bbox) {
        
        Curve_2 cv;
        
        std::ifstream inp(filename);
        if (!inp.is_open()) {
            std::cerr << "Cannot open file " << filename << "!" << std::endl;
            return -1;
        }
        int count;
        inp >> count;
        for (int i = 0; i < count; i++) {
            inp >> cv;
            curves.push_back(cv);
            // TODO bounding box
#if 0
            CGAL::Bbox_2 curve_bbox = cv.bounding_box();
            if (i == 0) bbox = curve_bbox;
            else bbox = bbox + curve_bbox;      
#else
            // TODO set these boxes
            bbox = CGAL::Bbox_2(-10, -10, 10, 10);
#endif
        }
        inp.close();
        return 0;
    }
    
    void skip_comments(std::ifstream& is, char* one_line) {
        while( !is.eof() ){
            is.getline( one_line, 128 );
            if( one_line[0] != '#' ){
                break;
            }
        }  
    }
};

#endif // EXACUS_CONIX_READER
