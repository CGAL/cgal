//#define CGAL_MESH_3_VERBOSE
#include "C3t3_type.h"
#include "Scene_c3t3_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_segmented_image_item.h"
#include "Scene_implicit_function_item.h"

#include "implicit_functions/Implicit_function_interface.h"
#include "Optimizer_thread.h"

#include <CGAL/optimize_mesh_3.h>
#include <CGAL/Bbox_3.h>

#include <fstream>

#include <CGAL/Timer.h>

namespace cgp = CGAL::parameters;


// -----------------------------------
// Optimization_function template class
// -----------------------------------
template <typename Domain, typename Function>
class Optimization_function
  : public Optimization_function_interface
{
public:
  // This class takes the responsability of d
  Optimization_function(Scene_c3t3_item* i, Domain* d, const Function& f)
    : item_(i), domain_(d), function_(f) {}
  
  virtual ~Optimization_function()
  {
    delete domain_;
  }
  
  virtual CGAL::Mesh_optimization_return_code launch() const
  {
    return function_(item_->c3t3(), *domain_);
  }
  
  virtual Scene_c3t3_item* item() const
  {
    return item_;
  }
  
private:
  Scene_c3t3_item* item_;
  Domain* domain_;
  Function function_;
};


// -----------------------------------
// Optimization generic function (responsible of dynamic casting)
// -----------------------------------
template <typename Function>
Optimizer_thread* cgal_code_optimization(Scene_c3t3_item& c3t3_item,
                                        const Function& f,
                                        const bool create_new_item)
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
    
    Image_mesh_domain* p_domain = new Image_mesh_domain(*p_image, 1e-6);
    
    // Create thread
    typedef Optimization_function<Image_mesh_domain,Function> Opt_function;
    Opt_function* p_opt_function = new Opt_function(p_result_item, p_domain, f);
    
    return new Optimizer_thread(p_opt_function);
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
    
    Mesh_domain* p_domain = new Mesh_domain(*p_poly);
    
    // Create thread
    typedef Optimization_function<Mesh_domain,Function> Opt_function;
    Opt_function* p_opt_function = new Opt_function(p_result_item, p_domain, f);
    
    return new Optimizer_thread(p_opt_function);
  }
  
  // Function
  const Scene_implicit_function_item* function_item = 
    qobject_cast<const Scene_implicit_function_item*>(c3t3_item.data_item());
  
  if ( NULL != function_item )
  {
    // Build domain
    const Implicit_function_interface* p_function = function_item->function();
    if ( NULL == p_function ) { return NULL; }
    
    CGAL::Bbox_3 dom_bbox (p_function->bbox().xmin,
                           p_function->bbox().ymin,
                           p_function->bbox().zmin,
                           p_function->bbox().xmax,
                           p_function->bbox().ymax,
                           p_function->bbox().zmax);
    
    Function_mesh_domain* p_domain =
      new Function_mesh_domain(Function_wrapper(*p_function), dom_bbox, 1e-7);
    
    // Create thread
    typedef Optimization_function<Function_mesh_domain,Function> Opt_function;
    Opt_function* p_opt_function = new Opt_function(p_result_item, p_domain, f);
    
    return new Optimizer_thread(p_opt_function);
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
  operator()(C3t3& c3t3, const Domain& domain) const
  {
    // Odt
    return CGAL::odt_optimize_mesh_3(c3t3,
                                     domain,
                                     cgp::time_limit = time_limit,
                                     cgp::convergence = convergence_ratio,
                                     cgp::freeze_bound = freeze_ratio,
                                     cgp::max_iteration_number = max_iteration_nb);
  }
};


Optimizer_thread*
cgal_code_odt_mesh_3(Scene_c3t3_item& c3t3_item,
                     const double time_limit,
                     const double convergence_ratio,
                     const double freeze_ratio,
                     const int max_iteration_number,
                     const bool create_new_item)
{
  Odt_function f;
  f.time_limit = time_limit;
  f.convergence_ratio = convergence_ratio;
  f.freeze_ratio = freeze_ratio;
  f.max_iteration_nb = max_iteration_number;
  
  return cgal_code_optimization(c3t3_item, f, create_new_item);
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
  operator()(C3t3& c3t3, const Domain& domain) const
  {
    // Lloyd
    return CGAL::lloyd_optimize_mesh_3(c3t3,
                                       domain,
                                       cgp::time_limit = time_limit,
                                       cgp::convergence = convergence_ratio,
                                       cgp::freeze_bound = freeze_ratio,
                                       cgp::max_iteration_number = max_iteration_nb);
  }
};


Optimizer_thread*
cgal_code_lloyd_mesh_3(Scene_c3t3_item& c3t3_item,
                       const double time_limit,
                       const double convergence_ratio,
                       const double freeze_ratio,
                       const int max_iteration_number,
                       const bool create_new_item)
{
  Lloyd_function f;
  f.time_limit = time_limit;
  f.convergence_ratio = convergence_ratio;
  f.freeze_ratio = freeze_ratio;
  f.max_iteration_nb = max_iteration_number;
  
  return cgal_code_optimization(c3t3_item, f, create_new_item);
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
  operator()(C3t3& c3t3, const Domain& domain) const
  {
    // Perturbation
    return CGAL::perturb_mesh_3(c3t3,
                                domain,
                                cgp::sliver_bound = sliver_bound,
                                cgp::time_limit = time_limit);
  }
};


Optimizer_thread*
cgal_code_perturb_mesh_3(Scene_c3t3_item& c3t3_item,
                         const double time_limit,
                         const double sliver_bound,
                         const bool create_new_item)
{
  Perturb_function f;
  f.sliver_bound = sliver_bound;
  f.time_limit = time_limit;
  
  return cgal_code_optimization(c3t3_item, f, create_new_item);
}


// -----------------------------------
// Exudation
// -----------------------------------
struct Exude_function
{
  double time_limit;
  double sliver_bound;
  
  CGAL::Mesh_optimization_return_code  
  operator()(C3t3& c3t3, int) const
  {
    // Perturbation
    return CGAL::exude_mesh_3(c3t3,
                              cgp::sliver_bound = sliver_bound,
                              cgp::time_limit = time_limit);
  }
};


Optimizer_thread*
cgal_code_exude_mesh_3(Scene_c3t3_item& c3t3_item,
                       const double time_limit,
                       const double sliver_bound,
                       const bool create_new_item)
{
  // Create result item
  Scene_c3t3_item* p_result_item = create_new_item ? 
    new Scene_c3t3_item(c3t3_item.c3t3()) : &c3t3_item;

  if ( NULL == p_result_item )
  {
    return NULL;
  }
  
  // Exudation
  Exude_function f;
  f.sliver_bound = sliver_bound;
  f.time_limit = time_limit;

  // Create thread
  typedef Optimization_function<int,Exude_function> Opt_function;
  Opt_function* p_opt_function = new Opt_function(p_result_item, NULL, f);
  
  return new Optimizer_thread(p_opt_function);
}