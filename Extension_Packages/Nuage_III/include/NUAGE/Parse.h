#ifndef PARSE_H
#define PARSE_H


#define AF_CGAL_CLIB_STD 
//#define AF_CGAL_CLIB_STD std

class Options {
public:
    Options()
        :  file_input(true), file_output(false),
	   Delaunay(false), contour(false), Section_file(false),
           number_of_points(-1), max_connected_comp(-1),
	   DELTA(.86), K_init(1.1), K_step(.1), K(5), out_format(0),
	   NB_BORDER_MAX(15)
  { 
    AF_CGAL_CLIB_STD::strcpy(finname,"finput"); // af: does not compile with bcc -O2 
    //std::strcpy(finname,"finput");

    AF_CGAL_CLIB_STD::strcpy(foutname,"foutput");// af: does not compile with bcc -O2
    //strcpy(foutname,"foutput"); 
  }

    char program[100];
    char finname[100];
    char foutname[100];
    bool file_input;
    bool file_output;
    bool Delaunay;
    bool contour;
    bool Section_file;
    int  number_of_points;
    int  max_connected_comp;
    double DELTA;
    double K_init;
    double K_step;
    double K;
    int out_format;
    int NB_BORDER_MAX;
};


void usage(char* program);

bool parse(int argc, char* argv[], Options &opt);

#endif
