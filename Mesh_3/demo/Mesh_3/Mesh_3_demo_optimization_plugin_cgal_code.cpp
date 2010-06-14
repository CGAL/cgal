//#define CGAL_MESH_3_VERBOSE
#include "C3t3_type.h"
#include "Scene_c3t3_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_segmented_image_item.h"
#include "Scene_implicit_function_item.h"

#include "implicit_functions/Implicit_function_interface.h"

#include <CGAL/optimize_mesh_3.h>
#include <CGAL/Bbox_3.h>

#include <fstream>

#include <CGAL/Timer.h>

namespace cgp = CGAL::parameters;
void treat_new_item(Scene_c3t3_item& new_item, const bool create_new_item);



// -----------------------------------
// Optimization generic function (responsible of dynamic casting)
// -----------------------------------
template <typename Function>
Scene_c3t3_item* cgal_code_optimization(Scene_c3t3_item& c3t3_item,
                                        const Function& f,
                                        const bool create_new_item,
                                        CGAL::Mesh_optimization_return_code& return_code)
{
  // Create result item
  Scene_c3t3_item* p_result_item = create_new_item ? 
    new Scene_c3t3_item(c3t3_item.c3t3()) : &c3t3_item;
  
  if ( NULL == p_result_item )
  {
    return NULL;
  }
  
  
  // Create domain using real type of c3t3_item.data_item()
  // ------------------
  
  // Image
  const Scene_segmented_image_item* image_item = 
    qobject_cast<const Scene_segmented_image_item*>(c3t3_item.data_item());
  
  if ( NULL != image_item )
  {
    // Build domain
    const Image* p_image = image_item->image();
    if ( NULL == p_image )
    {
      return NULL;
    }
    
    Image_mesh_domain domain(*p_image, 1e-6);
    
    // Launch
    return_code = f(p_result_item, domain);
    
    // Treat new item and exit
    treat_new_item(*p_result_item, create_new_item);
    return p_result_item;
  }
  
  
  // Polyhedron
  const Scene_polyhedron_item* poly_item = 
    qobject_cast<const Scene_polyhedron_item*>(c3t3_item.data_item());
  
  if ( NULL != poly_item )
  {
    // Build domain
    const Polyhedron* p_poly = poly_item->polyhedron();
    if ( NULL == p_poly )
    {
      return NULL;
    }
    
    Mesh_domain domain(*p_poly);
    
    // Launch
    return_code = f(p_result_item, domain);
    
    // Treat new item and exit
    treat_new_item(*p_result_item, create_new_item);
    return p_result_item;
  }
  
  // Function
  const Scene_implicit_function_item* function_item = 
    qobject_cast<const Scene_implicit_function_item*>(c3t3_item.data_item());
  
  if ( NULL != function_item )
  {
    // Build domain
    const Implicit_function_interface* p_function = function_item->function();
    if ( NULL == p_function ) { return NULL; }
    
    CGAL::Bbox_3 domain_bbox (p_function->bbox().xmin,
                              p_function->bbox().ymin,
                              p_function->bbox().zmin,
                              p_function->bbox().xmax,
                              p_function->bbox().ymax,
                              p_function->bbox().zmax);
    
    Function_mesh_domain domain(Function_wrapper(*p_function),
                                domain_bbox,
                                1e-7);
    
    // Launch
    return_code = f(p_result_item, domain);
    
    // Treat new item and exit
    treat_new_item(*p_result_item, create_new_item);
    return p_result_item;
  }
  
  return NULL;
}



// -----------------------------------
// Odt
// -----------------------------------
struct Odt_function
{
  double time_limit;
  double convergence_ratio;
  double freeze_ratio;
  int max_iteration_nb;
  
  template <typename Domain>
  CGAL::Mesh_optimization_return_code 
  operator()(Scene_c3t3_item* p_item, const Domain& domain) const
  {
    // Odt
    return CGAL::odt_optimize_mesh_3(p_item->c3t3(),
                                     domain,
                                     cgp::time_limit = time_limit,
                                     cgp::convergence = convergence_ratio,
                                     cgp::freeze_bound = freeze_ratio,
                                     cgp::max_iteration_number = max_iteration_nb);
  }
};


Scene_c3t3_item*
cgal_code_odt_mesh_3(Scene_c3t3_item& c3t3_item,
                     const double time_limit,
                     const double convergence_ratio,
                     const double freeze_ratio,
                     const int max_iteration_number,
                     const bool create_new_item,
                     CGAL::Mesh_optimization_return_code& return_code)
{
  Odt_function f;
  f.time_limit = time_limit;
  f.convergence_ratio = convergence_ratio;
  f.freeze_ratio = freeze_ratio;
  f.max_iteration_nb = max_iteration_number;
  
  return cgal_code_optimization(c3t3_item, f, create_new_item, return_code);
}



// -----------------------------------
// Lloyd
// -----------------------------------
struct Lloyd_function
{
  double time_limit;
  double convergence_ratio;
  double freeze_ratio;
  int max_iteration_nb;
  
  template <typename Domain>
  CGAL::Mesh_optimization_return_code  
  operator()(Scene_c3t3_item* p_item, const Domain& domain) const
  {
    // Lloyd
    return CGAL::lloyd_optimize_mesh_3(p_item->c3t3(),
                                       domain,
                                       cgp::time_limit = time_limit,
                                       cgp::convergence = convergence_ratio,
                                       cgp::freeze_bound = freeze_ratio,
                                       cgp::max_iteration_number = max_iteration_nb);
  }
};


Scene_c3t3_item*
cgal_code_lloyd_mesh_3(Scene_c3t3_item& c3t3_item,
                       const double time_limit,
                       const double convergence_ratio,
                       const double freeze_ratio,
                       const int max_iteration_number,
                       const bool create_new_item,
                       CGAL::Mesh_optimization_return_code& return_code)
{
  Lloyd_function f;
  f.time_limit = time_limit;
  f.convergence_ratio = convergence_ratio;
  f.freeze_ratio = freeze_ratio;
  f.max_iteration_nb = max_iteration_number;
  
  return cgal_code_optimization(c3t3_item, f, create_new_item, return_code);
}



// -----------------------------------
// Perturbation
// -----------------------------------
struct Perturb_function
{
  double time_limit;
  double sliver_bound;
  
  template <typename Domain>
  CGAL::Mesh_optimization_return_code  
  operator()(Scene_c3t3_item* p_item, const Domain& domain) const
  {
    // Perturbation
    return CGAL::perturb_mesh_3(p_item->c3t3(),
                                domain,
                                cgp::sliver_bound = sliver_bound,
                                cgp::time_limit = time_limit);
  }
};


Scene_c3t3_item*
cgal_code_perturb_mesh_3(Scene_c3t3_item& c3t3_item,
                         const double time_limit,
                         const double sliver_bound,
                         const bool create_new_item,
                         CGAL::Mesh_optimization_return_code& return_code)
{
  Perturb_function f;
  f.sliver_bound = sliver_bound;
  f.time_limit = time_limit;
  
  return cgal_code_optimization(c3t3_item, f, create_new_item, return_code);
}


// -----------------------------------
// Exudation
// -----------------------------------
Scene_c3t3_item*
cgal_code_exude_mesh_3(Scene_c3t3_item& c3t3_item,
                       const double time_limit,
                       const double sliver_bound,
                       const bool create_new_item,
                       CGAL::Mesh_optimization_return_code& return_code)
{
  // Create result item
  Scene_c3t3_item* p_result_item = create_new_item ? 
    new Scene_c3t3_item(c3t3_item.c3t3()) : &c3t3_item;

  if ( NULL == p_result_item )
  {
    return NULL;
  }
  
  // Exudation
  return_code = CGAL::exude_mesh_3(p_result_item->c3t3(),
                                   cgp::sliver_bound = sliver_bound,
                                   cgp::time_limit = time_limit);
  
  // Treat result and exit
  treat_new_item(*p_result_item, create_new_item);
  return p_result_item;
}



// -----------------------------------
// Helper functions
// -----------------------------------
void treat_new_item(Scene_c3t3_item& new_item, const bool create_new_item)
{
  if ( create_new_item )
  {
    const Scene_item::Bbox& bbox = new_item.bbox();
    new_item.setPosition((bbox.xmin + bbox.xmax)/2.f,
                         (bbox.ymin + bbox.ymax)/2.f,
                         (bbox.zmin + bbox.zmax)/2.f);
  }
  else
  {
    new_item.update_histogram();
  }  
}
