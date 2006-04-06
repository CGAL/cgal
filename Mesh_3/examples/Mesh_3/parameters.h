#include <map>
#include <string>

// type "pointer to implicit function"
typedef double (Implicit_function) (double, double, double);

typedef std::map<std::string, double> Double_options;
typedef std::map<std::string, Implicit_function*> Implicit_function_map;
typedef std::map<std::string, std::string> String_options;

extern Double_options double_options;
extern Implicit_function_map functions;
extern String_options string_options;

void init_parameters();
