#include <map>
#include <string>

extern int number_of_initial_points;
extern bool bipolar_oracle;

typedef std::map<std::string, double> Double_options;
extern Double_options double_options;

void init_parameters();
