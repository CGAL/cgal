#ifndef CGAL_MESH_3_LANTERI_PROCESS_RESULTS_H
#define CGAL_MESH_3_LANTERI_PROCESS_RESULTS_H

#include <iostream>
#include <vector>
#include <string>

bool process_cells(const std::vector<Qualities>& volume_cells_quality,
                   const std::vector<Qualities>& volume_cells_min_angle,
                   const std::string filename_prefix,
                   std::ostream* out_stream = &std::cout);

bool process_volume_edges(const std::vector<Qualities>& volume_edges_length,
                          const std::vector<double>& length_bounds,
                          const std::string filename_prefix,
                          std::ostream* out_stream = &std::cout);

bool process_surface_edges(const std::vector<Qualities>& surface_edges_length,
                           const std::vector<double>& length_bounds,
                           const std::string filename_prefix,
                           std::ostream* out_stream = &std::cout);
#endif // CGAL_MESH_3_LANTERI_PROCESS_RESULTS_H
