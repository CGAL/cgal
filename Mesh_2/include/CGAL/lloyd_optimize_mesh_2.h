#ifndef CGAL_LLOYD_OPTIMIZE_MESH_2_H
#define CGAL_LLOYD_OPTIMIZE_MESH_2_H

#include <CGAL/Mesh_2/Mesh_global_optimizer_2.h>
#include <CGAL/Mesh_2/Lloyd_move_2.h>
#include <CGAL/Mesh_2/Mesh_sizing_field.h>
#include <fstream>

#include <boost/parameter.hpp>
#include <boost/parameter/name.hpp>

BOOST_PARAMETER_NAME( cdt )
BOOST_PARAMETER_NAME( (max_iteration_number, tag) max_iteration_number_ )
BOOST_PARAMETER_NAME( (convergence, tag) convergence_)
BOOST_PARAMETER_NAME( (time_limit, tag) time_limit_ )
BOOST_PARAMETER_NAME( (freeze_bound, tag) freeze_bound_)

namespace CGAL
{
  BOOST_PARAMETER_FUNCTION(
  (void),
  lloyd_optimize_mesh_2,
  tag,
  (required (in_out(cdt),*))
  (optional
    (max_iteration_number_, *, 0 )
    (convergence_, *, 0.001 )
    (time_limit_, *, 0. )
    (freeze_bound_, *, 0.001 )
    )
  )
  {
    return lloyd_optimize_mesh_2_impl(cdt,
      max_iteration_number_,
      convergence_,
      freeze_bound_,
      time_limit_);
  }

  template<typename CDT>
  void lloyd_optimize_mesh_2_impl(CDT& cdt,
                             const unsigned int max_iterations,
                             const double convergence_ratio,
                             const double freeze_bound,
                             const double time_limit)
  {
    typedef Mesh_2::Mesh_sizing_field<CDT>           Sizing;
    typedef Mesh_2::Lloyd_move_2<CDT, Sizing>        Mv;
    typedef Mesh_2::Mesh_global_optimizer_2<CDT, Mv> Optimizer;

    Optimizer lloyd(cdt,
                    convergence_ratio,
                    freeze_bound);
    lloyd.set_time_limit(time_limit);

#ifdef CGAL_MESH_2_OPTIMIZERS_DEBUG
    std::ofstream os("before_lloyd.angles.txt");
    lloyd.output_angles_histogram(os);
    os.close();
#endif

    // 1000 iteration max to avoid infinite loop
    int nb_iterations = (0 == max_iterations)
      ? 1000
      : max_iterations;

    //run optimization
    lloyd(nb_iterations);

#ifdef CGAL_MESH_2_OPTIMIZERS_DEBUG
    std::ofstream os2("after_lloyd.angles.txt");
    lloyd.output_angles_histogram(os2);
    os2.close();
#endif
  }

} //end namespace CGAL

#endif
