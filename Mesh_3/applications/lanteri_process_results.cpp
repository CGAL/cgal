#include "distribution.h"
#include "Gd_displayer.h"

#include "lanteri_process_results.h"

#include <vector>
#include <string>

#include <iostream>

const int global_number_of_classes = 20;

bool process_aux_1(const std::vector<Qualities>& qualities,
                   const std::string filename,
                   const int number_of_classes,
                   std::ostream* out_stream = &std::cout,
                   // output stream
                   double max = 1.0
                   )
{
  const int maximum_index = qualities.size() - 1;

  Gd_displayer main_display(200*maximum_index, 200);

  main_display.set_window(-1.1, 1.1, -1.1, 1.1,
                          0, 0, 200, 200);

  std::vector<Gd_displayer*> displays(maximum_index);
  std::vector<Distribution> distributions(maximum_index);
  
  for(int i = 0; i < maximum_index; ++i)
  {
    distributions[i].resize(number_of_classes);

    compute_distribution(qualities[i+1],
                         max,
                         distributions[i]);

    displays[i] = new Gd_displayer(main_display.image());
    displays[i]->set_window(-.1, 1.1, -.1, 1.1,
                                 i * 200, 0, 200, 200);

    const int max = *(std::max_element(distributions[i].begin(),
                                       distributions[i].end()));
    
    display_distribution(displays[i],
                         distributions[i],
                         1. / max);
  }

  *out_stream << "saving " << filename.c_str() << "...\n";
  
  return main_display.save_png(filename.c_str());
}

bool process_aux_2(const std::vector<Qualities>& qualities,
                   const std::vector<double>& length_bounds,
                   const std::string filename,
                   const int number_of_classes,
                   std::ostream* out_stream = &std::cout
                   // output stream
                   )
{
  const int maximum_index = qualities.size() - 1;

  Gd_displayer main_display(200*maximum_index, 200);

  main_display.set_window(-1.1, 1.1, -1.1, 1.1,
                          0, 0, 200, 200);

  std::vector<Gd_displayer*> displays(maximum_index);
  std::vector<Distribution> distributions(maximum_index);
  
  for(int i = 0; i < maximum_index; ++i)
  {
    distributions[i].resize(number_of_classes);

    double max_quality = *(std::max_element(qualities[i+1].begin(),
                                            qualities[i+1].end()));

    max_quality = std::max(max_quality, length_bounds[i]);

    compute_distribution(qualities[i+1],
                         max_quality,
                         distributions[i]);

    displays[i] = new Gd_displayer(main_display.image());
    displays[i]->set_window(-.1, 1.1, -.1, 1.1,
                                 i * 200, 0, 200, 200);

    const int max = *(std::max_element(distributions[i].begin(),
                                               distributions[i].end()));
    
    display_distribution(displays[i],
                         distributions[i],
                         1. / max);

    double x_position_of_length_bound = length_bounds[i] / max_quality;
    
    displays[i]->segment(x_position_of_length_bound,  0.0,
                         x_position_of_length_bound, -0.05,
                         CGAL::BLUE);
  }

  *out_stream << "saving " << filename.c_str() << "...\n";

  return main_display.save_png(filename.c_str());
}

bool process_cells(const std::vector<Qualities>& volume_cells_quality,
                   const std::vector<Qualities>& volume_cells_min_angle,
                   const std::string filename_prefix,
                   std::ostream* out_stream
                   // output stream
                   )
{
  return process_aux_1(volume_cells_quality,
                       filename_prefix + "_cells_radius_radius_ratio.png",
                       global_number_of_classes,
                       out_stream)
  && process_aux_1(volume_cells_min_angle,
                   filename_prefix + "_cells_min_angles.png",
                   global_number_of_classes,
                   out_stream,
                   90);
}

bool process_volume_edges(const std::vector<Qualities>& volume_edges_length,
                          const std::vector<double>& length_bounds,
                          const std::string filename_prefix,
                          std::ostream* out_stream
                          // output stream
                          )
{
  return process_aux_2(volume_edges_length,
                       length_bounds,
                       filename_prefix + "_volume_edges_lengths.png",
                       global_number_of_classes,
                       out_stream);
}

bool process_surface_edges(const std::vector<Qualities>& surface_edges_length,
                           const std::vector<double>& length_bounds,
                           const std::string filename_prefix,
                           std::ostream* out_stream
                           // output stream
                           )
{
  return process_aux_2(surface_edges_length,
                       length_bounds,
                       filename_prefix + "_surface_edges_lengths.png",
                       global_number_of_classes,
                       out_stream);
}
