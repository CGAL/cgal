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

namespace cgp = CGAL::parameters;


// -----------------------------------
// Optimization_function_base template class
// -----------------------------------
template <typename Domain>
class Optimization_function_base
  : public Optimization_function_interface
{
public:
  /// Constructor
  /// Takes the responsability of d
  explicit
  Optimization_function_base(C3t3& c3t3, Domain* d)
    : c3t3_(c3t3), domain_(d) {}
  
  /// Destructor
  virtual ~Optimization_function_base()
  {
    delete domain_;
  }
  
  /// Launch
  virtual CGAL::Mesh_optimization_return_code launch()
  {
    return (*this)(c3t3_, *domain_);
  }
  
protected:
  /// Virtual operator() which should be overloaded
  virtual CGAL::Mesh_optimization_return_code
  operator()(C3t3& c3t3, const Domain& domain) = 0;
  
private:
  C3t3& c3t3_;
  Domain* domain_;
};

// Prototype which will be partially specialized for each Parameter class
template < typename Domain, typename Parameters >
class Optimization_function {};



// -----------------------------------
// Optimization generic function (responsible of dynamic casting)
// -----------------------------------
template <typename Parameters>
Optimizer_thread* cgal_code_optimization(Scene_c3t3_item& c3t3_item,
                                        const Parameters& param,
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
    typedef Optimization_function<Image_mesh_domain,Parameters> Opt_function;
    Opt_function* p_opt_function = new Opt_function(p_result_item->c3t3(), p_domain, param);
    
    return new Optimizer_thread(p_opt_function, p_result_item);
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
    
    Polyhedral_mesh_domain* p_domain = new Polyhedral_mesh_domain(*p_poly);
    
    // Create thread
    typedef Optimization_function<Polyhedral_mesh_domain,Parameters> Opt_function;
    Opt_function* p_opt_function = new Opt_function(p_result_item->c3t3(), p_domain, param);
    
    return new Optimizer_thread(p_opt_function, p_result_item);
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
    typedef Optimization_function<Function_mesh_domain,Parameters> Opt_function;
    Opt_function* p_opt_function = new Opt_function(p_result_item->c3t3(), p_domain, param);
    
    return new Optimizer_thread(p_opt_function, p_result_item);
  }
  
  return NULL;
}



// -----------------------------------
// Odt
// -----------------------------------
struct Odt_parameters
{
  double time_limit;
  double convergence_ratio;
  double freeze_ratio;
  int max_iteration_nb;
  
  QStringList log() const
  {
    return QStringList()
      << QString("time limit: %1").arg(time_limit)
      << QString("convergence ratio: %1").arg(convergence_ratio)
      << QString("freeze ratio: %1").arg(freeze_ratio)
      << QString("maximum iterations: %1").arg(max_iteration_nb);
  }
};


/**
 * @class Odt_function
 * Partial specialization of class Optimization_function for Odt
 * Runs odt global optimization
 */
template <typename Domain>
class Optimization_function < Domain, Odt_parameters >
  : public Optimization_function_base< Domain >
{
  // Private types
  typedef C3t3::Triangulation  Tr;
  typedef Tr::Geom_traits      Gt;
  
  typedef CGAL::Mesh_3::Mesh_sizing_field<Tr>    Sizing;
  typedef CGAL::Mesh_3::Odt_move<C3t3,Sizing>    Move;
  
  typedef typename CGAL::Mesh_3::Mesh_global_optimizer<C3t3,Domain,Move> Odt_optimizer;

  typedef Optimization_function_base< Domain > Base;  

public:
  /// Constructor
  Optimization_function(C3t3& c3t3, Domain* d, const Odt_parameters& p)
    : Base(c3t3,d)
    , odt_(NULL)
    , p_(p) {}
  
  /// Destructor
  virtual ~Optimization_function() { delete odt_; }

  /// Stops process (set time limit to 1ms)
  virtual void stop() { odt_->set_time_limit(0.001); }
  
  /// Log strings
  virtual QString name() const { return QString("Odt"); }
  virtual QStringList parameters_log() const { return p_.log(); }
  
protected:
  /// Launch odt optimization
  /// The content of this method is taken from CGAL::odt_optimize_mesh_3()
  virtual CGAL::Mesh_optimization_return_code 
  operator()(C3t3& c3t3, const Domain& domain)
  {
    if ( NULL != odt_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }

    // Create optimizer
    odt_ = new Odt_optimizer(c3t3, domain, p_.freeze_ratio, p_.convergence_ratio);
    if ( NULL == odt_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }
    
    // Set max time
    odt_->set_time_limit(p_.time_limit);
    
    // 1000 iteration max to avoid infinite loops
    int max_iteration_nb = ( 0 == p_.max_iteration_nb ) ? 1000
                                                       : p_.max_iteration_nb;
      
    // Launch optimization
    return (*odt_)(max_iteration_nb);
  }
  
private:
  Odt_optimizer* odt_;
  Odt_parameters p_;
};


/**
 * Global function cgal_code_odt_mesh_3
 */
Optimizer_thread*
cgal_code_odt_mesh_3(Scene_c3t3_item& c3t3_item,
                     const double time_limit,
                     const double convergence_ratio,
                     const double freeze_ratio,
                     const int max_iteration_number,
                     const bool create_new_item)
{
  Odt_parameters p;
  p.time_limit = time_limit;
  p.convergence_ratio = convergence_ratio;
  p.freeze_ratio = freeze_ratio;
  p.max_iteration_nb = max_iteration_number;
  
  return cgal_code_optimization(c3t3_item, p, create_new_item);
}



// -----------------------------------
// Lloyd
// -----------------------------------
struct Lloyd_parameters
{
  double time_limit;
  double convergence_ratio;
  double freeze_ratio;
  int max_iteration_nb;
  
  QStringList log() const
  {
    return QStringList()
      << QString("time limit: %1").arg(time_limit)
      << QString("convergence ratio: %1").arg(convergence_ratio)
      << QString("freeze ratio: %1").arg(freeze_ratio)
      << QString("maximum iterations: %1").arg(max_iteration_nb);
  }
};


/**
 * @class Lloyd_function
 * Partial specialization of class Optimization_function for Lloyd
 * Runs lloyd global optimization
 */
template <typename Domain>
class Optimization_function < Domain, Lloyd_parameters >
  : public Optimization_function_base< Domain >
{
  // Private types
  typedef C3t3::Triangulation  Tr;
  typedef Tr::Geom_traits      Gt;
  
  typedef CGAL::Mesh_3::Mesh_sizing_field<Tr>    Sizing;
  typedef CGAL::Mesh_3::Lloyd_move<C3t3,Sizing>  Move;
  
  typedef typename CGAL::Mesh_3::Mesh_global_optimizer<C3t3,Domain,Move> Lloyd_optimizer;
  
  typedef Optimization_function_base< Domain > Base;
  
public:
  /// Constructor
  Optimization_function(C3t3& c3t3, Domain* d, const Lloyd_parameters& p)
    : Base(c3t3,d)
    , lloyd_(NULL)
    , p_(p) {}
  
  /// Destructor
  virtual ~Optimization_function() { delete lloyd_; }

  /// Stops process (set time limit to 1ms)
  virtual void stop() { lloyd_->set_time_limit(0.001); }
  
  /// Log strings
  virtual QString name() const { return QString("Lloyd"); }
  virtual QStringList parameters_log() const { return p_.log(); }
  
protected:
  /// Launch lloyd optimization
  /// The content of this method is taken from CGAL::lloyd_optimize_mesh_3()
  virtual CGAL::Mesh_optimization_return_code 
  operator()(C3t3& c3t3, const Domain& domain)
  {
    if ( NULL != lloyd_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }
    
    // Create optimizer
    lloyd_ = new Lloyd_optimizer(c3t3, domain, p_.freeze_ratio, p_.convergence_ratio);
    if ( NULL == lloyd_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }
    
    // Set max time
    lloyd_->set_time_limit(p_.time_limit);
    
    // 1000 iteration max to avoid infinite loops
    int max_iteration_nb = ( 0 == p_.max_iteration_nb ) ? 1000
                                                        : p_.max_iteration_nb;
    
    // Launch optimization
    return (*lloyd_)(max_iteration_nb);
  }
  
private:
  Lloyd_optimizer* lloyd_;
  Lloyd_parameters p_;
};


/**
 * Global function cgal_code_lloyd_mesh_3
 */
Optimizer_thread*
cgal_code_lloyd_mesh_3(Scene_c3t3_item& c3t3_item,
                       const double time_limit,
                       const double convergence_ratio,
                       const double freeze_ratio,
                       const int max_iteration_number,
                       const bool create_new_item)
{
  Lloyd_parameters p;
  p.time_limit = time_limit;
  p.convergence_ratio = convergence_ratio;
  p.freeze_ratio = freeze_ratio;
  p.max_iteration_nb = max_iteration_number;
  
  return cgal_code_optimization(c3t3_item, p, create_new_item);
}



// -----------------------------------
// Perturbation
// -----------------------------------
struct Perturb_parameters
{
  double time_limit;
  double sliver_bound;
  
  QStringList log() const
  {
    return QStringList()
      << QString("time limit: %1").arg(time_limit)
      << QString("sliver bound: %1").arg(sliver_bound);
  }
};


/**
 * @class Perturb_function
 * Partial specialization of class Optimization_function for perturbation
 * Runs sliver perturbation
 */
template <typename Domain>
class Optimization_function < Domain, Perturb_parameters >
  : public Optimization_function_base< Domain >
{
  // Private types
  typedef C3t3::Triangulation::Geom_traits                    Gt;
  typedef CGAL::Mesh_3::Min_dihedral_angle_criterion<Gt>      Sc;
  
  typedef CGAL::Mesh_3::Sliver_perturber<C3t3,Domain,Sc>      Perturber;
  
  typedef Optimization_function_base< Domain > Base;
  
public:
  /// Constructor
  Optimization_function(C3t3& c3t3, Domain* d, const Perturb_parameters& p)
    : Base(c3t3,d)
    , perturb_(NULL)
    , p_(p) {}
  
  /// Destructor
  ~Optimization_function() { delete perturb_; }

  /// Stops process (set time limit to 1ms)
  virtual void stop() { perturb_->set_time_limit(0.001); }
  
  /// Log strings
  virtual QString name() const { return QString("Perturb"); }
  virtual QStringList parameters_log() const { return p_.log(); }
  
protected:
  /// Launch sliver perturbation
  /// The content of this method is taken from CGAL::perturb_mesh_3()
  virtual CGAL::Mesh_optimization_return_code  
  operator()(C3t3& c3t3, const Domain& domain)
  {
    if ( NULL != perturb_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }

    typedef CGAL::Mesh_3::Sq_radius_perturbation<C3t3,Domain,Sc>      Sq_radius;
    typedef CGAL::Mesh_3::Volume_perturbation<C3t3,Domain,Sc>         Volume;
    typedef CGAL::Mesh_3::Dihedral_angle_perturbation<C3t3,Domain,Sc> Dihedral_angle;
    typedef CGAL::Mesh_3::Li_random_perturbation<C3t3,Domain,Sc>      Li_random;
    
    // Build perturber
    perturb_ = new Perturber(c3t3,domain);
    if ( NULL == perturb_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }
    
    // Add perturbations
    perturb_->add_perturbation(new Sq_radius(40,0.02));
    perturb_->add_perturbation(new Volume(40,0.02));
    perturb_->add_perturbation(new Dihedral_angle(40,0.02));
    perturb_->add_perturbation(new Li_random(100,0.05));
    
    // Set max time
    perturb_->set_time_limit(p_.time_limit);
    
    // Launch perturber
    if ( p_.sliver_bound != 0 ) { return (*perturb_)(p_.sliver_bound); }
    else                       { return (*perturb_)(); }
  }

private:
  Perturber* perturb_;
  Perturb_parameters p_;
};


/**
 * Global function cgal_code_perturb_mesh_3
 */
Optimizer_thread*
cgal_code_perturb_mesh_3(Scene_c3t3_item& c3t3_item,
                         const double time_limit,
                         const double sliver_bound,
                         const bool create_new_item)
{
  Perturb_parameters p;
  p.sliver_bound = sliver_bound;
  p.time_limit = time_limit;
  
  return cgal_code_optimization(c3t3_item, p, create_new_item);
}


// -----------------------------------
// Exudation
// -----------------------------------
struct Exude_parameters
{
  double time_limit;
  double sliver_bound;
  
  QStringList log() const
  {
    return QStringList()
      << QString("time limit: %1").arg(time_limit)
      << QString("sliver bound: %1").arg(sliver_bound);
  }
};


/**
 * @class Exude_function
 * Partial specialization of class Optimization_function for exudation
 * Runs sliver exudation
 */
template <typename Domain>
class Optimization_function < Domain, Exude_parameters >
  : public Optimization_function_base< Domain >
{
  // Private types
  typedef C3t3::Triangulation::Geom_traits                  Gt;
  typedef CGAL::Mesh_3::Min_dihedral_angle_criterion<Gt>    Sc;
  typedef CGAL::Mesh_3::Slivers_exuder<C3t3, Sc>            Exuder;
  
  typedef Optimization_function_base< Domain > Base;
  
public:
  // Constructor
  Optimization_function(C3t3& c3t3, Domain* d, const Exude_parameters& p)
    : Base(c3t3,d)
    , exude_(NULL)
    , p_(p) {}
  
  /// Destructor
  ~Optimization_function() { delete exude_; }
  
  /// Stops process (set time limit to 1ms)
  virtual void stop() { exude_->set_time_limit(0.001); }
  
  // Log strings
  virtual QString name() const { return QString("Exude"); }
  virtual QStringList parameters_log() const { return p_.log(); }
  
protected:
  /// Launch sliver exudation
  /// The content of this method is taken from CGAL::exude_mesh_3()
  virtual CGAL::Mesh_optimization_return_code  
  operator()(C3t3& c3t3, const Domain&)
  {
    if ( NULL != exude_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }
    
    // Create exuder
    exude_ = new Exuder(c3t3);
    if ( NULL == exude_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }
    
    // Set time_limit
    exude_->set_time_limit(p_.time_limit);
    
    // Launch exudation
    return (*exude_)(p_.sliver_bound);
  }
  
private:
  Exuder* exude_;
  Exude_parameters p_;
};


/**
 * Global function cgal_code_exude_mesh_3
 */
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
  Exude_parameters p;
  p.sliver_bound = sliver_bound;
  p.time_limit = time_limit;

  // Create thread
  typedef Optimization_function<int,Exude_parameters> Opt_function;
  Opt_function* p_opt_function = new Opt_function(p_result_item->c3t3(), NULL, p);
  
  return new Optimizer_thread(p_opt_function, p_result_item);
}