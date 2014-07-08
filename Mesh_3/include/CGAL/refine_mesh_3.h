// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : refine_mesh_3 function declaration and implementation.
//******************************************************************************

#ifndef CGAL_REFINE_MESH_3_H
#define CGAL_REFINE_MESH_3_H

#include <CGAL/config.h>
#include <CGAL/Mesh_3/config.h>
#include <CGAL/Mesh_3/Dump_c3t3.h>
#include <CGAL/Mesh_3/global_parameters.h>
#include <CGAL/Mesh_3/Mesher_3.h>
#include <CGAL/optimize_mesh_3.h>

namespace CGAL {

  namespace details {
    
    /**
     * @class Insert_vertex_in_c3t3
     *
     * A functor designed to insert unweighted points into the triangulation
     * of a C3T3 from C3T3::Tr::Vertex , keeping the dimension and indices.
     */
    template <typename C3T3>
    class Insert_vertex_in_c3t3
    {
    private:
      typedef typename C3T3::Vertex_handle          Vertex_handle;
      typedef typename C3T3::Index                  Index;
      
      typedef typename C3T3::Triangulation          Tr;
      typedef typename Tr::Vertex                   Vertex;
      typedef typename Tr::Point                    Weighted_point;
      typedef typename Tr::Point::Weight            Weight;
      
    public:
      Insert_vertex_in_c3t3(C3T3& c3t3)
        : r_c3t3_(c3t3) {}
      
      void operator()(const Vertex& vertex) const
      {
        // Get vh properties
        int dimension = vertex.in_dimension();
        Weight w = (dimension < 2) ? vertex.point().weight() : 0;
        Weighted_point point(vertex.point().point(),w);
        Index index = vertex.index();

        // Insert point and restore handle properties
        Vertex_handle new_vertex = r_c3t3_.triangulation().insert(point);
        r_c3t3_.set_index(new_vertex, index);
        r_c3t3_.set_dimension(new_vertex, dimension);

#if defined(CGAL_LINKED_WITH_TBB)\
 && !defined(CGAL_PARALLEL_MESH_3_DO_NOT_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE)
        if (boost::is_convertible<typename C3T3::Concurrency_tag, CGAL::Parallel_tag>::value)
        {
          if (dimension == -1)
            r_c3t3_.add_far_point(new_vertex);
        }
#endif
#ifdef CGAL_SEQUENTIAL_MESH_3_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE
        if (boost::is_convertible<typename C3T3::Concurrency_tag, CGAL::Sequential_tag>::value)
        {
          if (dimension == -1)
            r_c3t3_.add_far_point(new_vertex);
        }
#endif
      }

    private:
      C3T3& r_c3t3_;
    };
  }
  
  
namespace parameters {
  
  namespace internal {
    
    const int undef_parameter = -1;
    
    // Helpers
    struct Optimization_options_base
    {
      Optimization_options_base(bool b)
      : b_(b), time_limit_(undef_parameter), bound_(undef_parameter) {}
      
      operator bool() const { return b_; }
      
      bool is_time_limit_set() const { return time_limit_ != undef_parameter; }
      void set_time_limit(double d) { time_limit_ = d; }
      double time_limit() const { return time_limit_; }
      
      bool is_bound_set() const { return bound_ != undef_parameter; }
      void set_bound(double d) { bound_ = d; }
      double bound() const { return bound_; }
      
    private:
      bool b_;
      double time_limit_;
      double bound_;
    };
    
    struct Global_optimization_options_base
    {
      Global_optimization_options_base()
      : convergence_(undef_parameter), max_it_nb_(undef_parameter) {}
      
      bool is_convergence_set() const { return convergence_ != undef_parameter; }
      void set_convergence(double d) { convergence_ = d; }
      double convergence() const { return convergence_; }
      
      bool is_max_iteration_number_set() const { return max_it_nb_ != undef_parameter; }
      void set_max_iteration_number(int i) { max_it_nb_ = i; }
      int max_iteration_number() const { return max_it_nb_; }
      
    private:
      double convergence_;
      int max_it_nb_; 
    };
    
    // Perturb
    struct Perturb_options : public Optimization_options_base
    {
      Perturb_options(bool b) : Optimization_options_base(b) {}
    };
    
    // Exude
    struct Exude_options : public Optimization_options_base
    {
      Exude_options(bool b) : Optimization_options_base(b) {}
    };
    
    // Odt
    struct Odt_options : public Optimization_options_base
    , public Global_optimization_options_base
    {
      Odt_options(bool b) : Optimization_options_base(b)
      , Global_optimization_options_base() {}
    };
    
    // Lloyd
    struct Lloyd_options : public Optimization_options_base
    , public Global_optimization_options_base
    {
      Lloyd_options(bool b) : Optimization_options_base(b)
      , Global_optimization_options_base() {}
    };

    // Various Mesh_3 option
    struct Mesh_3_options {
      Mesh_3_options() 
        : dump_after_init_prefix()
        , dump_after_refine_surface_prefix()
        , dump_after_refine_prefix()
        , dump_after_glob_opt_prefix()
        , dump_after_perturb_prefix()
        , dump_after_exude_prefix()
      {}

      std::string dump_after_init_prefix;
      std::string dump_after_refine_surface_prefix;
      std::string dump_after_refine_prefix;
      std::string dump_after_glob_opt_prefix;
      std::string dump_after_perturb_prefix;
      std::string dump_after_exude_prefix;

    }; // end struct Mesh_3_options

  } // end namespace internal
#if defined(__clang__) || defined(__GNUC__) && CGAL_GCC_VERSION >= 40600
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
// remove a warning about an unused parameter "args" in the definition of
// BOOST_PARAMETER_FUNCTION
#endif
  // -----------------------------------
  // Perturb
  // -----------------------------------
  BOOST_PARAMETER_FUNCTION((internal::Perturb_options), perturb, tag,
                           (optional (time_limit_, *, internal::undef_parameter )
                                     (sliver_bound_, *, default_values::perturb_sliver_bound )))
  { 
    internal::Perturb_options options(true);
    
    if ( internal::undef_parameter != time_limit_ )
      options.set_time_limit(time_limit_);
    
    options.set_bound(sliver_bound_);
    
    return options;
  }
  
  inline internal::Perturb_options no_perturb() { return internal::Perturb_options(false); }
  
  // -----------------------------------
  // Exude
  // -----------------------------------
  BOOST_PARAMETER_FUNCTION((internal::Exude_options), exude, tag,
                           (optional (time_limit_, *, internal::undef_parameter )
                                     (sliver_bound_, *, default_values::exude_sliver_bound )))
  { 
    internal::Exude_options options(true);
    
    if ( internal::undef_parameter != time_limit_ )
      options.set_time_limit(time_limit_);

    options.set_bound(sliver_bound_);
    
    return options;
  }
  
  inline internal::Exude_options no_exude() { return internal::Exude_options(false); }
  
  // -----------------------------------
  // Odt
  // -----------------------------------
  BOOST_PARAMETER_FUNCTION((internal::Odt_options), odt, tag,
                           (optional (time_limit_, *, 0 )
                                     (max_iteration_number_, *, 0 )
                                     (convergence_, *, default_values::odt_convergence_ratio )
                                     (freeze_bound_, *, default_values::odt_freeze_ratio )))
  { 
    internal::Odt_options options(true);

    options.set_time_limit(time_limit_);
    options.set_bound(freeze_bound_);
    options.set_convergence(convergence_);
    options.set_max_iteration_number(max_iteration_number_);
    
    return options;
  }
  
  inline internal::Odt_options no_odt() { return internal::Odt_options(false); }
  
  // -----------------------------------
  // Lloyd
  // -----------------------------------
  BOOST_PARAMETER_FUNCTION((internal::Lloyd_options), lloyd, tag,
                           (optional (time_limit_, *, 0 )
                                     (max_iteration_number_, *, 0 )
                                     (convergence_, *, default_values::lloyd_convergence_ratio )
                                     (freeze_bound_, *, default_values::lloyd_freeze_ratio )))
  { 
    internal::Lloyd_options options(true);
    
    options.set_time_limit(time_limit_);
    options.set_bound(freeze_bound_);
    options.set_convergence(convergence_);
    options.set_max_iteration_number(max_iteration_number_);
    
    return options;
  }
  
  inline internal::Lloyd_options no_lloyd() { return internal::Lloyd_options(false); }
  
  // Mesh options
  // -----------------------------------

  // Undocumented Boost parameter for refine_mesh_3 and make_mesh_3.
  // Allows to dump the mesh at given stage of the mesh generation
  // algorithm.
  BOOST_PARAMETER_FUNCTION((internal::Mesh_3_options), mesh_3_options, tag,
                           (optional 
                            (dump_after_init_prefix_, (std::string), "" )
                            (dump_after_refine_surface_prefix_, (std::string), "" )
                            (dump_after_refine_prefix_, (std::string), "" )
                            (dump_after_glob_opt_prefix_, (std::string), "" )
                            (dump_after_perturb_prefix_, (std::string), "" )
                            (dump_after_exude_prefix_, (std::string), "" )
                            )
                           )
  { 
    internal::Mesh_3_options options;

    options.dump_after_init_prefix=dump_after_init_prefix_;
    options.dump_after_refine_surface_prefix=dump_after_refine_surface_prefix_;
    options.dump_after_refine_prefix=dump_after_refine_prefix_;
    options.dump_after_glob_opt_prefix=dump_after_glob_opt_prefix_;
    options.dump_after_perturb_prefix=dump_after_perturb_prefix_;
    options.dump_after_exude_prefix=dump_after_exude_prefix_;
    
    return options;
  }
  
  // Undocumented Boost parameter for refine_mesh_3 and make_mesh_3.
  // Default Mesh_3_options: dump at every stage of the mesh generation.
  inline internal::Mesh_3_options mesh_3_dump()
  {
    internal::Mesh_3_options options;

    options.dump_after_init_prefix="mesh_dump_after_init";
    options.dump_after_refine_surface_prefix="mesh_dump_after_refine_surface";
    options.dump_after_refine_prefix="mesh_dump_after_refine";
    options.dump_after_glob_opt_prefix="mesh_dump_after_glob_opt";
    options.dump_after_perturb_prefix="mesh_dump_after_perturb";
    options.dump_after_exude_prefix="mesh_dump_after_exude";
    
    return options;
  }

#if defined(__clang__) || defined(__GNUC__) && CGAL_GCC_VERSION >= 40600
//#if defined(__GNUC__) && __GNUC__ >= 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic pop
#endif
  
  // -----------------------------------
  // Reset_c3t3 (undocumented)
  // -----------------------------------
  CGAL_MESH_BOOLEAN_PARAMETER(Reset,reset_c3t3,no_reset_c3t3)
  // CGAL_MESH_BOOLEAN_PARAMETER defined in <CGAL/Mesh_3/global_parameters.h>
  
// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/Mesh_3/config.h>
CGAL_MESH_3_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

  // -----------------------------------
  // Parameters
  // -----------------------------------
  BOOST_PARAMETER_NAME( exude_param )
  BOOST_PARAMETER_NAME( perturb_param )
  BOOST_PARAMETER_NAME( odt_param )
  BOOST_PARAMETER_NAME( lloyd_param )
  BOOST_PARAMETER_NAME( reset_param )
  BOOST_PARAMETER_NAME( mesh_options_param )

CGAL_PRAGMA_DIAG_POP
} // end namespace parameters
  
  
BOOST_PARAMETER_FUNCTION(
  (void),
  refine_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*) (criteria,*) ) // nondeduced
  (deduced 
    (optional
      (exude_param, (parameters::internal::Exude_options), parameters::exude())
      (perturb_param, (parameters::internal::Perturb_options), parameters::perturb())
      (odt_param, (parameters::internal::Odt_options), parameters::no_odt())
      (lloyd_param, (parameters::internal::Lloyd_options), parameters::no_lloyd())
      (reset_param, (parameters::Reset), parameters::reset_c3t3())
      (mesh_options_param, (parameters::internal::Mesh_3_options), 
                           parameters::internal::Mesh_3_options())
    )
  )
)
{
  return refine_mesh_3_impl(c3t3,
                            domain,
                            criteria,
                            exude_param,
                            perturb_param,
                            odt_param,
                            lloyd_param,
                            reset_param(),
                            mesh_options_param);
}
  
  
/**
 * @brief This function refines the mesh c3t3 wrt domain & criteria
 *
 * @param c3t3 the mesh to be refined.
 * @param domain the domain to be discretized
 * @param criteria the criteria
 * @param exude if \c true, an exudation step will be done at
 *   the end of the Delaunay refinement process
 * @param perturb if \c true, an explicit vertex perturbation step will be
 *   done at the end of refinement process
 * @param reset_c3t3 if \c true, a new C3T3 will be construct from param c3t3.
 *   The new c3t3 keeps only the vertices (as NON-weighted points with their
 *   dimension and Index) of the triangulation. That allows to refine a mesh
 *   which has been exuded.
 * @param mesh_3_various_options is a struct object used to pass
 * non-documented options, for debugging purpose.
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void refine_mesh_3_impl(C3T3& c3t3,
                        const MeshDomain&   domain,
                        const MeshCriteria& criteria,
                        const parameters::internal::Exude_options& exude,
                        const parameters::internal::Perturb_options& perturb,
                        const parameters::internal::Odt_options& odt,
                        const parameters::internal::Lloyd_options& lloyd,
                        bool reset_c3t3,
                        const parameters::internal::Mesh_3_options& 
                          mesh_options = parameters::internal::Mesh_3_options())
{
  typedef Mesh_3::Mesher_3<C3T3, MeshCriteria, MeshDomain> Mesher;

  // Reset c3t3 (i.e. remove weights) if needed
  if ( reset_c3t3 )
  {
    C3T3 tmp_c3t3;
    std::for_each(c3t3.triangulation().finite_vertices_begin(),
                  c3t3.triangulation().finite_vertices_end(),
                  details::Insert_vertex_in_c3t3<C3T3>(tmp_c3t3));
    // TODO: corners and edges are not restored
    c3t3.swap(tmp_c3t3);
  }
  
  dump_c3t3(c3t3, mesh_options.dump_after_init_prefix);

  // Build mesher and launch refinement process
  Mesher mesher (c3t3, domain, criteria);
  double refine_time = mesher.refine_mesh(mesh_options.dump_after_refine_surface_prefix);

  dump_c3t3(c3t3, mesh_options.dump_after_refine_prefix);

  // Odt
  if ( odt )
  {
    odt_optimize_mesh_3(c3t3,
                        domain,
                        parameters::time_limit = odt.time_limit(),
                        parameters::max_iteration_number = odt.max_iteration_number(),
                        parameters::convergence = odt.convergence(),
                        parameters::freeze_bound = odt.bound());
  }
  
  // Lloyd
  if ( lloyd )
  {
    lloyd_optimize_mesh_3(c3t3,
                          domain,
                          parameters::time_limit = lloyd.time_limit(),
                          parameters::max_iteration_number = lloyd.max_iteration_number(),
                          parameters::convergence = lloyd.convergence(),
                          parameters::freeze_bound = lloyd.bound());
  }
    
  dump_c3t3(c3t3, mesh_options.dump_after_glob_opt_prefix);

  // Perturbation
  if ( perturb )
  {
    double perturb_time_limit = refine_time;
    
    if ( perturb.is_time_limit_set() )
      perturb_time_limit = perturb.time_limit();

    perturb_mesh_3(c3t3,
                   domain,
                   parameters::time_limit = perturb_time_limit,
                   parameters::sliver_bound = perturb.bound());
  }
  
  dump_c3t3(c3t3, mesh_options.dump_after_perturb_prefix);

  // Exudation
  if ( exude )
  {
    double exude_time_limit = refine_time;
    
    if ( exude.is_time_limit_set() )
      exude_time_limit = exude.time_limit();
    
    exude_mesh_3(c3t3,
                 parameters::time_limit = exude_time_limit,
                 parameters::sliver_bound = exude.bound());
  }
  
  dump_c3t3(c3t3, mesh_options.dump_after_exude_prefix);

}

}  // end namespace CGAL


#endif // CGAL_REFINE_MESH_3_H
