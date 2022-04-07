#include "config.h"
#include "config_mesh_3.h"
#include "C3t3_type.h"
#include "Scene_c3t3_item.h"
#include "Scene_surface_mesh_item.h"

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
#include "Scene_image_item.h"
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
#include "Scene_implicit_function_item.h"
#include "implicit_functions/Implicit_function_interface.h"
#endif

#include "Optimizer_thread.h"

#include <CGAL/optimize_mesh_3.h>
#include <CGAL/Bbox_3.h>

#include <fstream>
#include <cstddef>

namespace cgp = CGAL::parameters;


// -----------------------------------
// Helper function
// -----------------------------------
QString translate_bool(const bool b)
{
  return b ? QString("done")
           : QString("in progress");
}


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

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
  // Image
  const Scene_image_item* image_item =
    qobject_cast<const Scene_image_item*>(c3t3_item.data_item());

  if ( NULL != image_item )
  {
    // Build domain
    const Image* p_image = image_item->image();

    if ( NULL == p_image )
    {
      return NULL;
    }

    Image_mesh_domain* p_domain =
      new Image_mesh_domain(Image_mesh_domain::create_labeled_image_mesh_domain
                              (CGAL::parameters::image = *p_image,
                               CGAL::parameters::relative_error_bound = 1e-6,
                               CGAL::parameters::construct_surface_patch_index =
                                    [](int i, int j) { return (i * 1000 + j); } ));

    // Create thread
    typedef Optimization_function<Image_mesh_domain,Parameters> Opt_function;
    Opt_function* p_opt_function = new Opt_function(p_result_item->c3t3(), p_domain, param);

    return new Optimizer_thread(p_opt_function, p_result_item);
  }
#endif

  // Surface mesh
  const Scene_surface_mesh_item* sm_item =
    qobject_cast<const Scene_surface_mesh_item*>(c3t3_item.data_item());

  if ( NULL != sm_item )
  {
    const_cast<Scene_surface_mesh_item*>(sm_item)->setItemIsMulticolor(true);
    const_cast<Scene_surface_mesh_item*>(sm_item)->computeItemColorVectorAutomatically(true);
    // Build domain
    const SMesh* smesh = sm_item->face_graph();
    if ( NULL == smesh )
    {
      return NULL;
    }
    Polyhedral_mesh_domain* sm_domain = new Polyhedral_mesh_domain(*smesh);
    if(c3t3_item.get_sharp_edges_angle() != -1 )
      sm_domain->detect_features(c3t3_item.get_sharp_edges_angle());
    else if(c3t3_item.get_detect_borders())
      sm_domain->detect_borders();

    // Create thread
    typedef Optimization_function<Polyhedral_mesh_domain,Parameters> Opt_function;
    Opt_function* p_opt_function = new Opt_function(p_result_item->c3t3(), sm_domain, param);
    return new Optimizer_thread(p_opt_function, p_result_item);
  }

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  // Function
  const Scene_implicit_function_item* function_item =
    qobject_cast<const Scene_implicit_function_item*>(c3t3_item.data_item());

  if ( NULL != function_item )
  {
    // Build domain
    const Implicit_function_interface* p_function = function_item->function();
    if ( NULL == p_function ) { return NULL; }

    CGAL::Bbox_3 dom_bbox (p_function->bbox().xmin(),
                           p_function->bbox().ymin(),
                           p_function->bbox().zmin(),
                           p_function->bbox().xmax(),
                           p_function->bbox().ymax(),
                           p_function->bbox().zmax());

    Function_mesh_domain* p_domain =
      new Function_mesh_domain(Function_wrapper(*p_function), dom_bbox, 1e-7,
                               CGAL::parameters::construct_surface_patch_index =
                                 [](int i, int j) { return (i * 1000 + j); } );

    // Create thread
    typedef Optimization_function<Function_mesh_domain,Parameters> Opt_function;
    Opt_function* p_opt_function = new Opt_function(p_result_item->c3t3(), p_domain, param);

    return new Optimizer_thread(p_opt_function, p_result_item);
  }
#endif

  return NULL;
}



// -----------------------------------
// Global optimization
// -----------------------------------
struct Global_optimization_status
{
  bool compute_moves_done;
  bool move_points_done;
  bool rebuild_restricted_delaunay_done;
  int iteration_done;

  Global_optimization_status()
    : compute_moves_done(false)
    , move_points_done(false)
    , rebuild_restricted_delaunay_done(false)
    , iteration_done(-1) {}
};

struct Global_visitor
{
  Global_visitor(Global_optimization_status* status) : p_status_(status) {}
  Global_visitor(const Global_visitor& rhs) : p_status_(rhs.p_status_) {}

  void after_compute_moves() { p_status_->compute_moves_done = true; }
  void after_move_points() { p_status_->move_points_done = true; }
  void after_rebuild_restricted_delaunay() { p_status_->rebuild_restricted_delaunay_done = true; }

  void end_of_iteration(int iteration_number)
  {
    p_status_->iteration_done = iteration_number;
    p_status_->compute_moves_done = false;
    p_status_->move_points_done = false;
    p_status_->rebuild_restricted_delaunay_done = false;
  }

private:
  Global_optimization_status* p_status_;
};


template <typename Domain>
class Global_optimization_function
  : public Optimization_function_base< Domain >
{
  typedef Global_visitor                       Visitor;
  typedef Optimization_function_base< Domain > Base;

public:
  /// Constructor
  Global_optimization_function(C3t3& c3t3, Domain* d)
    : Base(c3t3,d)
    , status_() {}

  /// Destructor
  virtual ~Global_optimization_function() {}

  // Logs
  virtual QString status(double) const
  {
    QString res = QString("Iteration %1<br /><br />"
                          "Compute moves: %2<br />")
    .arg(status_.iteration_done + 2)
    .arg(translate_bool(status_.compute_moves_done));

    if ( status_.compute_moves_done )
    {
      res += QString("Move points: %1<br />")
      .arg(translate_bool(status_.move_points_done));
    }

    if ( status_.move_points_done )
    {
      res += QString("Rebuild restricted Delaunay: %1")
      .arg(translate_bool(status_.rebuild_restricted_delaunay_done));
    }

    return res;
  }

protected:
  Global_optimization_status status_;
};


// -----------------------------------
// Odt
// -----------------------------------

#ifndef CGAL_MESH_3_DEMO_DISABLE_ODT

struct Odt_parameters
{
  double time_limit;
  double convergence_ratio;
  double freeze_ratio;
  bool do_freeze;
  int max_iteration_nb;

  QStringList log() const
  {
    return QStringList()
      << QString("time limit: %1").arg(time_limit)
      << QString("convergence ratio: %1").arg(convergence_ratio)
      << QString("freeze ratio: %1").arg(freeze_ratio)
      << QString("do freeze: %1").arg(do_freeze)
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
  : public Global_optimization_function< Domain >
{
  // Private types
  typedef C3t3::Triangulation  Tr;
  typedef CGAL::Mesh_3::Mesh_sizing_field<Tr>    Sizing;
  typedef CGAL::Mesh_3::Odt_move<C3t3,Sizing>    Move;
  typedef Global_visitor                         Visitor;

  typedef typename CGAL::Mesh_3::Mesh_global_optimizer<C3t3,Domain,Move,Visitor> Odt_optimizer;

  typedef Global_optimization_function< Domain > Base;

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
    odt_ = new Odt_optimizer(c3t3, domain, p_.freeze_ratio, p_.do_freeze, p_.convergence_ratio);
    if ( NULL == odt_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }

    // Set max time
    odt_->set_time_limit(p_.time_limit);

    // 1000 iteration max to avoid infinite loops
    int max_iteration_nb = ( 0 == p_.max_iteration_nb ) ? 1000
                                                       : p_.max_iteration_nb;

    // Launch optimization
    return (*odt_)(max_iteration_nb, Visitor(&(this->status_)));
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
#endif



// -----------------------------------
// Lloyd
// -----------------------------------

#ifndef CGAL_MESH_3_DEMO_DISABLE_LLOYD

struct Lloyd_parameters
{
  double time_limit;
  double convergence_ratio;
  double freeze_ratio;
  bool do_freeze;
  int max_iteration_nb;

  QStringList log() const
  {
    return QStringList()
      << QString("time limit: %1").arg(time_limit)
      << QString("convergence ratio: %1").arg(convergence_ratio)
      << QString("freeze ratio: %1").arg(freeze_ratio)
      << QString("do freeze: %1").arg(do_freeze)
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
  : public Global_optimization_function< Domain >
{
  // Private types
  typedef C3t3::Triangulation  Tr;
  typedef CGAL::Mesh_3::Mesh_sizing_field<Tr>    Sizing;
  typedef CGAL::Mesh_3::Lloyd_move<C3t3,Sizing>  Move;
  typedef Global_visitor                         Visitor;

  typedef typename CGAL::Mesh_3::Mesh_global_optimizer<C3t3,Domain,Move,Visitor> Lloyd_optimizer;

  typedef Global_optimization_function< Domain > Base;

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
    lloyd_ = new Lloyd_optimizer(c3t3, domain, p_.freeze_ratio, p_.do_freeze, p_.convergence_ratio);
    if ( NULL == lloyd_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }

    // Set max time
    lloyd_->set_time_limit(p_.time_limit);

    // 1000 iteration max to avoid infinite loops
    int max_iteration_nb = ( 0 == p_.max_iteration_nb ) ? 1000
                                                        : p_.max_iteration_nb;

    // Launch optimization
    return (*lloyd_)(max_iteration_nb, Visitor(&(this->status_)));
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
#endif



// -----------------------------------
// Perturbation
// -----------------------------------

#ifndef CGAL_MESH_3_DEMO_DISABLE_PERTURBER

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


struct Perturb_status
{
  double bound_reached;
  double vertices_left;

  Perturb_status() : bound_reached(0), vertices_left(0) {}
};

struct Perturb_visitor
{
  Perturb_visitor(Perturb_status* status) : p_status_(status) {}
  Perturb_visitor(const Perturb_visitor& rhs) : p_status_(rhs.p_status_) {}

  void bound_reached(const double bound) { p_status_->bound_reached = bound; }
  void end_of_perturbation_iteration(std::size_t v) { p_status_->vertices_left = v;}

private:
  Perturb_status* p_status_;
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
  typedef C3t3::Triangulation                                     Tr;
  typedef CGAL::Mesh_3::Min_dihedral_angle_criterion<Tr>          Sc;
  typedef Perturb_visitor                                         Visitor;

  typedef CGAL::Mesh_3::Sliver_perturber<C3t3,Domain,Sc,Visitor>  Perturber;

  typedef Optimization_function_base< Domain > Base;

public:
  /// Constructor
  Optimization_function(C3t3& c3t3, Domain* d, const Perturb_parameters& p)
    : Base(c3t3,d)
    , perturb_(NULL)
    , p_(p)
    , criterion_(p.sliver_bound, c3t3.triangulation()) {}

  /// Destructor
  ~Optimization_function() { delete perturb_; }

  /// Stops process (set time limit to 1ms)
  virtual void stop() { perturb_->set_time_limit(0.001); }

  /// Log strings
  virtual QString name() const { return QString("Perturb"); }
  virtual QStringList parameters_log() const { return p_.log(); }
  virtual QString status(double) const
  {
    return QString("Dihedral angle reached: %1<br /><br />"
                   "Vertices left in queue (to reach next bound): %2")
    .arg(status_.bound_reached)
    .arg(status_.vertices_left);
  }

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
    perturb_ = new Perturber(c3t3, domain, criterion_);
    if ( NULL == perturb_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }

    // Add perturbations
    perturb_->add_perturbation(new Sq_radius(40,0.05));
    perturb_->add_perturbation(new Volume(40,0.05));
    perturb_->add_perturbation(new Dihedral_angle(40,0.05));
    perturb_->add_perturbation(new Li_random(100,0.15));

    // Set max time
    perturb_->set_time_limit(p_.time_limit);

    // Set sliver bound (0 means no sliver bound)
    if ( 0 == p_.sliver_bound ) { p_.sliver_bound = criterion_.get_max_value(); }

    // Launch perturber
    return (*perturb_)(Visitor(&status_));
  }

private:
  Perturber* perturb_;
  Perturb_parameters p_;
  Perturb_status status_;
  Sc criterion_;
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
#endif

// -----------------------------------
// Exudation
// -----------------------------------

#ifndef CGAL_MESH_3_DEMO_DISABLE_EXUDER

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


struct Exude_status
{
  double cells_left_in_queue;
  double vertices_pumped;

  Exude_status() : cells_left_in_queue(0), vertices_pumped(0) {}
};

struct Exude_visitor
{
  Exude_visitor(Exude_status* status) : p_status_(status) {}
  Exude_visitor(const Exude_visitor& rhs) : p_status_(rhs.p_status_) {}

  void after_cell_pumped(std::size_t n) { p_status_->cells_left_in_queue = n; }

private:
  Exude_status* p_status_;
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
  typedef C3t3::Triangulation                               Tr;
  typedef CGAL::Mesh_3::Min_dihedral_angle_criterion<Tr>    Sc;
  typedef Exude_visitor                                     Visitor;
  typedef CGAL::Mesh_3::Slivers_exuder<C3t3,Sc,Visitor>     Exuder;

  typedef Optimization_function_base< Domain > Base;

public:
  // Constructor
  Optimization_function(C3t3& c3t3, Domain* d, const Exude_parameters& p)
    : Base(c3t3,d)
    , exude_(NULL)
    , p_(p)
    , criterion_(p.sliver_bound, c3t3.triangulation()) {}

  /// Destructor
  ~Optimization_function() { delete exude_; }

  /// Stops process (set time limit to 1ms)
  virtual void stop() { exude_->set_time_limit(0.001); }

  // Log strings
  virtual QString name() const { return QString("Exude"); }
  virtual QStringList parameters_log() const { return p_.log(); }
  virtual QString status(double) const
  {
    return QString("Cells left in queue: %1<br />")
                   //"Vertices pumped: %2")
      .arg(status_.cells_left_in_queue);
      //.arg(status_.vertices_pumped);
  }

protected:
  /// Launch sliver exudation
  /// The content of this method is taken from CGAL::exude_mesh_3()
  virtual CGAL::Mesh_optimization_return_code
  operator()(C3t3& c3t3, const Domain&)
  {
    if ( NULL != exude_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }

    // Create exuder
    exude_ = new Exuder(c3t3, criterion_);
    if ( NULL == exude_ ) { return CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR; }

    // Set time_limit
    exude_->set_time_limit(p_.time_limit);

    // Launch exudation
    return (*exude_)(Visitor(&status_));
  }

private:
  Exuder* exude_;
  Exude_parameters p_;
  Exude_status status_;
  Sc criterion_;
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
#endif
