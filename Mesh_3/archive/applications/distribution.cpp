#include "distribution.h"

void compute_distribution(const Qualities& qualities,
                          const double max_quality,
                          Distribution& distribution)
{
  const int number_of_classes = distribution.size();

  const int qualities_size = qualities.size();

  for(int j = 0; j < qualities_size; ++j)
  {
    if(qualities[j] < max_quality)
    {
      ++distribution[static_cast<int>((qualities[j]/max_quality)
                                      *number_of_classes)];
    }
  }
}

void display_distribution(Distribution_displayer* display,
                          const Distribution& distribution,
                          const double echelle)
{
  const int number_of_classes = distribution.size();

  if( number_of_classes == 0 ) return;
  const double width = 1.0 / number_of_classes;

  display->fill_rectangle(0., 0., 1., 1., CGAL::IO::Color(200, 200, 200));
  for(int k = 0; k < number_of_classes; ++k)
    if(distribution[k]>0)
    {
      const double height = ( distribution[k]+0. ) * echelle;
      display->fill_rectangle(k    * width, 0,
                              (k+1)* width, height,
                              CGAL::black());
    }
    else
      display->segment(k     * width, 0.,
                       (k+1) * width, 0.,
                       CGAL::IO::red());
}
