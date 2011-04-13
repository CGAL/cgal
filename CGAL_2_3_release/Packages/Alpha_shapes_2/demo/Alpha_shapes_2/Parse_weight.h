#ifndef PARSE_H
#define PARSE_H

//#include <CGAL/bool.h>

class Options {
public:
    Options()
        :  file_input(true), file_output(false),
	   Delaunay(false), contour(false),
           regularized(false), init(false), number_of_points(5000),
           min(0.0), max(612.0), winx(512), winy(512)
  { 
    CGAL_CLIB_STD::strcpy(finname,"./Data/fin");
    CGAL_CLIB_STD::strcpy(foutname,"./Data/fout");
  }

    char program[100];
    char finname[100];
    char foutname[100];
    bool file_input;
    bool file_output;
    bool Delaunay;
    bool contour;
    bool regularized;
    bool weight;
    bool init;
    int  number_of_points;
    double min;
    double max;
    int winx;
    int winy;
};


void usage(char* program);

bool
parse(int argc, char* argv[], Options &opt);

#endif
