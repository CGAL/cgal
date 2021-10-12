#ifndef CGAL_KSR_ALL_PARAMETERS_EXAMPLES_H
#define CGAL_KSR_ALL_PARAMETERS_EXAMPLES_H

// STL includes.
#include <string>

namespace CGAL {
namespace KSR {

  template<typename FT>
  struct All_parameters {

    // Path to the input data file.
    std::string data;

    // Label indices defined in the ply header:
    // ground (gi),
    // building boundary (bi),
    // building interior (ii),
    // vegetation (vi).
    std::string gi, bi, ii, vi;

    // Main parameters.
    FT scale; // meters
    FT noise; // meters

    // Boolean tags.
    const bool with_normals; // do we use normals
    const bool verbose;
    const bool debug;

    // Shape detection / shape regularization.
    std::size_t k_neighbors;
    FT distance_threshold;
    FT angle_threshold;
    std::size_t min_region_size;
    bool regularize;

    // Partitioning.
    unsigned int k_intersections;
    const unsigned int n_subdivisions;
    const FT enlarge_bbox_ratio;
    const bool reorient;

    // Reconstruction.
    FT graphcut_beta;

    // Constructor.
    All_parameters() :
    data(""),
    gi("0"), bi("1"), ii("2"), vi("3"),
    // main parameters
    scale(FT(4)),
    noise(FT(2)),
    // boolean tags
    with_normals(true),
    verbose(true),
    debug(true),
    // shape detection / shape regularization
    k_neighbors(12),
    distance_threshold(noise / FT(2)),
    angle_threshold(FT(15)),
    min_region_size(50),
    regularize(false),
    // partitioning
    k_intersections(1),
    n_subdivisions(0),
    enlarge_bbox_ratio(FT(11) / FT(10)),
    reorient(false),
    // reconstruction
    graphcut_beta(FT(1) / FT(2))
    { }

    // Update all parameters, which depend on scale and noise.
    void update_dependent() {
      distance_threshold = noise / FT(2);
    }
  };

} // KSR
} // CGAL

#endif // CGAL_KSR_ALL_PARAMETERS_EXAMPLES_H
