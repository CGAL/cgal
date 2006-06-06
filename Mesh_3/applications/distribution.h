#include <vector>
#include "Distribution_displayer.h"

typedef std::vector<double> Qualities;
typedef std::vector<int> Distribution;

void compute_distribution(const Qualities& qualities,
                          const double max_quality,
                          Distribution& distribution);

void display_distribution(Distribution_displayer* display,
                          const Distribution& distribution,
                          const double echelle);
