#include <map>
#include <set>
#include <string>

// type "pointer to implicit function"
typedef double (Implicit_function) (double, double, double);

typedef std::map<std::string, double> Double_options;
typedef std::map<std::string, Implicit_function*> Implicit_function_map;
typedef std::map<std::string, std::string> String_options;
typedef std::set<std::string> Check_strings;

extern String_options string_options;
extern Double_options double_options;
extern Implicit_function_map functions;
extern Check_strings check_strings;

typedef void (Usage)(std::string = "");
extern Usage* usage_ptr;

void init_parameters();

#ifdef CGAL_SURFACE_MESHER_TEST_OPTIONS
void check_all_options_have_been_used(); // debug function
#endif

std::string get_string_option(std::string);

double get_double_option(std::string);

Implicit_function* get_function(std::string);
