// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_SIMULATION_TRAITS_H
#define CGAL_KINETIC_SIMULATION_TRAITS_H
#include <CGAL/Kinetic/basic.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Kinetic/Active_objects_vector.h>
#include <CGAL/Kinetic/Cartesian_instantaneous_kernel.h>
#include <CGAL/Kinetic/Cartesian_kinetic_kernel.h>
#include <CGAL/Kinetic/Derivitive_filter_function_kernel.h>
#include <CGAL/Kinetic/Handle_degeneracy_function_kernel.h>
#include <CGAL/Kinetic/Simulator.h>
#include <CGAL/Kinetic/Two_list_pointer_event_queue.h>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/Numeric_root_stack.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Kinetic/Heap_pointer_event_queue.h>

CGAL_KINETIC_BEGIN_NAMESPACE

template <class StaticKernel, class InstantaneousKernel,
	  class KineticKernel, class SimulatorC, class ActiveObjectsTable>
struct Simulation_traits
{
public:

  typedef ActiveObjectsTable Active_objects_table;

  typedef typename StaticKernel::FT NT;
  struct Static_kernel: public StaticKernel {};

  struct Kinetic_kernel: public KineticKernel {};

  typedef SimulatorC Simulator;

  //typedef typename Simulator::Function_kernel Function_kernel;

  typedef typename KineticKernel::Function_kernel::Root Time;

  typedef InstantaneousKernel Instantaneous_kernel;

  Simulation_traits(const Time &lb, const Time &ub): sp_(new Simulator(lb, ub)), ao_(new ActiveObjectsTable){}
  Simulation_traits(): sp_(new Simulator()), ao_(new ActiveObjectsTable) {
  }
  /*Simulation_traits(typename Simulator::Pointer sp): sp_(sp){
    }*/

  Simulator* simulator_pointer(){ return sp_.get();}
  const Simulator* simulator_pointer() const { return sp_.get();}
  Active_objects_table* active_objects_table_pointer(){ return ao_.get();}
  const Active_objects_table* active_objects_table_pointer() const { return ao_.get();}
  Static_kernel& static_kernel_object(){return k_;}
  Kinetic_kernel& kinetic_kernel_object(){return kk_;}
  //Function_kernel function_kernel_object(){return sp_->function_kernel_object();}

  const Static_kernel& static_kernel_object() const {return k_;}
  const Kinetic_kernel& kinetic_kernel_object() const {return kk_;}
  //const Function_kernel function_kernel_object() const {return sp_->function_kernel_object();}
  Instantaneous_kernel instantaneous_kernel_object() const
  {
    return Instantaneous_kernel(ao_, static_kernel_object());
  }
protected:
  Static_kernel k_;
  Kinetic_kernel kk_;
  typename Simulator::Pointer sp_;
  typename Active_objects_table::Pointer ao_;
};
CGAL_KINETIC_END_NAMESPACE

CGAL_KINETIC_BEGIN_INTERNAL_NAMESPACE
struct Sest_types
{
  typedef CGAL::Simple_cartesian<CGAL::Gmpq> Static_kernel;
  typedef Static_kernel::FT NT;
  typedef CGAL::POLYNOMIAL::Polynomial<NT> Function;
  typedef CGAL::POLYNOMIAL::Upper_bound_root_stack_Descartes_traits<Function> Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Upper_bound_root_stack<Root_stack_traits> Root_stack;
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel;
  struct Simulator_function_kernel: public CGAL::Kinetic::Handle_degeneracy_function_kernel<Function_kernel> {};
  typedef CGAL::Kinetic::Cartesian_kinetic_kernel<Simulator_function_kernel> Kinetic_kernel;
  typedef  Simulator_function_kernel::Root Time;
  typedef CGAL::Kinetic::Two_list_pointer_event_queue<Time, double> Queue_base;

  struct Event_queue: public Queue_base
  {
    Event_queue(const Time &start, const Time &end): Queue_base(start, end){}
  };

  typedef CGAL::Kinetic::Simulator<Simulator_function_kernel, Event_queue > Simulator;
};

template <class ActiveObject>
struct Suggested_exact_simulation_traits:
  public Simulation_traits<typename Sest_types::Static_kernel,
			   Cartesian_instantaneous_kernel<Active_objects_vector<ActiveObject>,
							  typename Sest_types::Static_kernel>,
			   typename Sest_types::Kinetic_kernel,
			   typename Sest_types::Simulator,
			   Active_objects_vector<ActiveObject> >
{
  typedef Simulation_traits<typename Sest_types::Static_kernel,
			    Cartesian_instantaneous_kernel<Active_objects_vector<ActiveObject>,
							   typename Sest_types::Static_kernel>,
			    typename Sest_types::Kinetic_kernel,
			    typename Sest_types::Simulator,
			    Active_objects_vector<ActiveObject> > P;
  Suggested_exact_simulation_traits(const typename P::Time &lb,
				    const typename P::Time &ub): P(lb,ub){}
};

struct Sist_types
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Static_kernel;
  typedef Static_kernel::FT NT;
  typedef CGAL::POLYNOMIAL::Polynomial<NT> Function;
  typedef CGAL::POLYNOMIAL::Root_stack_default_traits<Function> Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Numeric_root_stack<Root_stack_traits> Root_stack;
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel;
  typedef CGAL::Kinetic::Derivitive_filter_function_kernel<Function_kernel> Simulator_function_kernel;
  typedef CGAL::Kinetic::Cartesian_kinetic_kernel<Simulator_function_kernel> Kinetic_kernel;
  typedef Simulator_function_kernel::Root Time;
  typedef CGAL::Kinetic::Heap_pointer_event_queue<Time> Queue_base;

  struct Event_queue: public Queue_base
  {
    Event_queue(const Time &start, const Time &finish): Queue_base(start, finish){}
  };
  typedef CGAL::Kinetic::Simulator<Simulator_function_kernel, Event_queue > Simulator;
};

template <class ActiveObject>
struct Suggested_inexact_simulation_traits: public Simulation_traits<typename Sist_types::Static_kernel,
								     Cartesian_instantaneous_kernel<Active_objects_vector<ActiveObject>,
												    typename Sist_types::Static_kernel>,
								     typename Sist_types::Kinetic_kernel,
								     typename Sist_types::Simulator,
								     Active_objects_vector<ActiveObject> >
{
  typedef Simulation_traits<typename Sist_types::Static_kernel,
			    Cartesian_instantaneous_kernel<Active_objects_vector<ActiveObject>,
							   typename Sist_types::Static_kernel>,
			    typename Sist_types::Kinetic_kernel,
			    typename Sist_types::Simulator,
			    Active_objects_vector<ActiveObject> > P;
  Suggested_inexact_simulation_traits(const typename P::Time &lb,
				      const typename P::Time &ub): P(lb,ub){}
};

CGAL_KINETIC_END_INTERNAL_NAMESPACE
#endif
