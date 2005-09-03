#include "distribution.h"
#include "Gd_displayer.h"

#include <vector>
#include <string>

const int global_number_of_classes = 20;

bool process_aux_1(const std::vector<Qualities>& qualities,
                   const std::string filename,
                   const int number_of_classes)
{
  Gd_displayer main_display(1000, 200);

  main_display.set_window(-1.1, 1.1, -1.1, 1.1,
                          0, 0, 200, 200);

  std::vector<Gd_displayer*> displays(5);
  std::vector<Distribution> distributions(5);
  
  for(int i = 0; i < 5; ++i)
  {
    distributions[i].resize(number_of_classes);

    compute_distribution(qualities[i+1],
                         1.,
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

  return main_display.save_png(filename.c_str());
}

bool process_aux_2(const std::vector<Qualities>& qualities,
                   const std::string filename,
                   const int number_of_classes)
{
  Gd_displayer main_display(1000, 200);

  main_display.set_window(-1.1, 1.1, -1.1, 1.1,
                          0, 0, 200, 200);

  std::vector<Gd_displayer*> displays(5);
  std::vector<Distribution> distributions(5);
  
  for(int i = 0; i < 5; ++i)
  {
    distributions[i].resize(number_of_classes);

    compute_distribution(qualities[i+1],
                         *(std::max_element(qualities[i+1].begin(),
                                            qualities[i+1].end())),
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

  return main_display.save_png(filename.c_str());
}

bool process_cells(const std::vector<Qualities>& volume_cells_quality)
{
  return process_aux_1(volume_cells_quality,
                       "lanteri_cells.png",
                       global_number_of_classes);
}

bool process_volume_edges(const std::vector<Qualities>& volume_edges_lenght)
{
  return process_aux_2(volume_edges_lenght,
                       "lanteri_volume_edges.png",
                       global_number_of_classes);
}

bool process_surface_edges(const std::vector<Qualities>& surface_edges_lenght)
{
  return process_aux_2(surface_edges_lenght,
                       "lanteri_surface_edges.png",
                       global_number_of_classes);
}
