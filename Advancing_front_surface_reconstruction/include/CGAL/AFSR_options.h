#ifndef CGAL_AFSR_OPTIONS_H
#define CGAL_AFSR_OPTIONS_H


namespace CGAL {

class AFSR_options {
public:
    AFSR_options()
        :  file_input(true), file_output(false),
	   Delaunay(false), contour(false), binary(false), xyz(false), 
           Section_file(false), max_connected_comp(-1),
	   delta(.86), K_init(1.1), K_step(.1), K(5), out_format(0),
	   NB_BORDER_MAX(15), red(0), green(0), blue(0), no_header(false), area(0), perimeter(0), 
	   abs_area(0), abs_perimeter(0)
  { 
    std::strcpy(finname,"finput");
    std::strcpy(foutname,"foutput"); 
  }

  char program[100];
  char finname[100];
  char foutname[100];
  bool file_input;
  bool file_output;
  bool Delaunay;
  bool contour;
  bool binary;
  bool xyz;
  bool Section_file;
  int  max_connected_comp;
  double delta;
  double K_init;
  double K_step;
  double K;
  int out_format;
  int NB_BORDER_MAX;
  double red, green, blue;
  bool no_header;
  double area, perimeter, abs_area, abs_perimeter;
};


} // namespace CGAL

#endif  // CGAL_AFSR_OPTIONS_H
