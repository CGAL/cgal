#include <CGAL/Polygonal_surface_reconstruction.h>

#include "Kernel_type.h"
#include "SMesh_type.h"
#include "Scene_points_with_normal_item.h"

#ifdef CGAL_USE_SCIP
#  include <CGAL/SCIP_mixed_integer_program_traits.h>
#endif

#ifdef CGAL_USE_GLPK
#  include <CGAL/GLPK_mixed_integer_program_traits.h>
#endif

typedef        CGAL::Polygonal_surface_reconstruction<Kernel> Polygonal_surface_reconstruction;

SMesh* polygonal_reconstruct (const Point_set& points,
                              double data_fitting,
                              double data_coverage,
                              double model_complexity,
                              const QString& solver_name)
{
  // Avoid warnings if no solver is available
  CGAL_USE (data_fitting);
  CGAL_USE (data_coverage);
  CGAL_USE (model_complexity);
  CGAL_USE (solver_name);

  Point_set::Property_map<int> shape_map
    = points.property_map<int>("shape").first;

        Polygonal_surface_reconstruction poly
    (points, points.point_map(), points.normal_map(), shape_map);

  SMesh* mesh = new SMesh;

#if defined(CGAL_USE_SCIP)
  if (solver_name == QString("SCIP"))
  {
    if (!poly.reconstruct<CGAL::SCIP_mixed_integer_program_traits<double> >(*mesh,
                                                                            data_fitting,
                                                                            data_coverage,
                                                                            model_complexity))
    {
      std::cerr << "Error: " << poly.error_message() << std::endl;
      delete mesh;
      return nullptr;
    }
  }
#endif

#if defined(CGAL_USE_GLPK)
  if (solver_name == QString("GLPK"))
  {
    if (!poly.reconstruct<CGAL::GLPK_mixed_integer_program_traits<double> >(*mesh,
                                                                            data_fitting,
                                                                            data_coverage,
                                                                            model_complexity))
    {
      std::cerr << "Error: " << poly.error_message() << std::endl;
      delete mesh;
      return nullptr;
    }
  }
#endif

  return mesh;
}
