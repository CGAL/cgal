#include <map>
#include <string>
#include <CGAL/Gray_level_image_3.h>

// type "pointer to implicit function"
typedef double (Implicit_function) (double, double, double);

extern Implicit_function generic_inrimage_function;

typedef std::map<std::string, double> Double_options;
typedef std::map<std::string, Implicit_function*> Implicit_function_map;
typedef std::map<std::string, std::string> String_options;

extern Double_options double_options;
extern Implicit_function_map functions;
extern String_options string_options;
typedef CGAL::Gray_level_image_3<double> Gray_image;
extern Gray_image* isosurface;

void init_parameters();
