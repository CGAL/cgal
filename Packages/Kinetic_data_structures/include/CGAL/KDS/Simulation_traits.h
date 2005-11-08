// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_SIMULATION_TRAITS_H
#define CGAL_KDS_SIMULATION_TRAITS_H
#include <CGAL/KDS/basic.h>



#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Gmpq.h>
#include <CGAL/KDS/Active_objects_vector.h>
#include <CGAL/KDS/Cartesian_instantaneous_kernel.h>
#include <CGAL/KDS/Cartesian_kinetic_kernel.h>
#include <CGAL/KDS/Derivitive_filter_function_kernel.h>
#include <CGAL/KDS/Handle_degeneracy_function_kernel.h>
#include <CGAL/KDS/Simulator.h>
#include <CGAL/KDS/Two_list_pointer_event_queue.h>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/Numeric_root_stack.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
#include <CGAL/Simple_cartesian.h>
//#include <CGAL/KDS/Heap_pointer_event_queue.h>


CGAL_KDS_BEGIN_NAMESPACE

template <class StaticKernel, class InstantaneousKernel,
	  class KineticKernel, class SimulatorC, class ActiveObjectsTable>
struct Simulation_traits {
public:
  
  typedef ActiveObjectsTable Active_objects_table;

  typedef typename StaticKernel::FT NT;
  struct Static_kernel: public StaticKernel {};

  struct Kinetic_kernel: public KineticKernel {};
  
  typedef SimulatorC Simulator;

  typedef typename Simulator::Function_kernel Function_kernel;

  typedef typename Function_kernel::Root Time;

  typedef InstantaneousKernel Instantaneous_kernel;

  Simulation_traits(const Time &lb, const Time &ub): sp_(new Simulator(lb, ub)), ao_(new ActiveObjectsTable){}
  Simulation_traits(): sp_(new Simulator()), ao_(new ActiveObjectsTable){
  }
  /*Simulation_traits(typename Simulator::Pointer sp): sp_(sp){
    }*/

  Simulator* simulator_pointer(){ return sp_.get();}
  const Simulator* simulator_pointer() const { return sp_.get();}
  Active_objects_table* active_objects_table_pointer(){ return ao_.get();}
  const Active_objects_table* active_objects_table_pointer() const { return ao_.get();}
  Static_kernel& static_kernel_object(){return k_;}
  Kinetic_kernel& kinetic_kernel_object(){return kk_;}
  Function_kernel function_kernel_object(){return sp_->function_kernel_object();}

  const Static_kernel& static_kernel_object() const {return k_;}
  const Kinetic_kernel& kinetic_kernel_object() const {return kk_;}
  const Function_kernel function_kernel_object() const {return sp_->function_kernel_object();}
  Instantaneous_kernel instantaneous_kernel_object() const {
    return Instantaneous_kernel(ao_, static_kernel_object());
  }
protected:
  Static_kernel k_;
  Kinetic_kernel kk_;
  typename Simulator::Pointer sp_;
  typename Active_objects_table::Pointer ao_;
};
CGAL_KDS_END_NAMESPACE






CGAL_KDS_BEGIN_INTERNAL_NAMESPACE
struct Sest_types {
  typedef CGAL::Simple_cartesian<CGAL::Gmpq> Static_kernel;
  typedef Static_kernel::FT NT;
  typedef CGAL::POLYNOMIAL::Polynomial<NT> Function;
  typedef CGAL::POLYNOMIAL::Upper_bound_root_stack_Descartes_traits<Function> Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Upper_bound_root_stack<Root_stack_traits> Root_stack;
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel;
  typedef CGAL::KDS::Handle_degeneracy_function_kernel<Function_kernel> Simulator_function_kernel;
  typedef CGAL::KDS::Cartesian_kinetic_kernel<Function_kernel> Kinetic_kernel;
  typedef  Simulator_function_kernel::Root Time;
  typedef CGAL::KDS::Two_list_pointer_event_queue<Time, double> Queue_base;
 
  struct Event_queue: public Queue_base{
    Event_queue(const Time &start, const Time &end): Queue_base(start, end){}
  };
  
  typedef CGAL::KDS::Simulator<Simulator_function_kernel, Event_queue > Simulator;
};

template <class ActiveObject>
struct Suggested_exact_simulation_traits: 
  public Simulation_traits<typename Sest_types::Static_kernel,
			   Cartesian_instantaneous_kernel<Active_objects_vector<ActiveObject>,
							  typename Sest_types::Static_kernel>,
			   typename Sest_types::Kinetic_kernel,
			   typename Sest_types::Simulator,
			   Active_objects_vector<ActiveObject> > {
  typedef Simulation_traits<typename Sest_types::Static_kernel,
			    Cartesian_instantaneous_kernel<Active_objects_vector<ActiveObject>,
							   typename Sest_types::Static_kernel>,
			    typename Sest_types::Kinetic_kernel,
			    typename Sest_types::Simulator,
			    Active_objects_vector<ActiveObject> > P;
  Suggested_exact_simulation_traits(const typename P::Time &lb, 
				    const typename P::Time &ub): P(lb,ub){}
};


struct Sist_types {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Static_kernel;
  typedef Static_kernel::FT NT;
  typedef CGAL::POLYNOMIAL::Polynomial<NT> Function;
  typedef CGAL::POLYNOMIAL::Root_stack_default_traits<Function> Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Numeric_root_stack<Root_stack_traits> Root_stack;
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel;
  typedef CGAL::KDS::Derivitive_filter_function_kernel<Function_kernel> Simulator_function_kernel;
  typedef CGAL::KDS::Cartesian_kinetic_kernel<Function_kernel> Kinetic_kernel;
  typedef  Simulator_function_kernel::Root Time;
  typedef CGAL::KDS::Heap_pointer_event_queue<Time> Queue_base;

  struct Event_queue: public Queue_base{
    Event_queue(const Time &start, const Time &finish): Queue_base(start, finish){}
  };
  typedef CGAL::KDS::Simulator<Simulator_function_kernel, Event_queue > Simulator;
};


template <class ActiveObject>
struct Suggested_inexact_simulation_traits: public Simulation_traits<typename Sist_types::Static_kernel,
								     Cartesian_instantaneous_kernel<Active_objects_vector<ActiveObject>,
												   typename Sist_types::Static_kernel>,
								     typename Sist_types::Kinetic_kernel,
								     typename Sist_types::Simulator,
								     Active_objects_vector<ActiveObject> > {
  typedef Simulation_traits<typename Sist_types::Static_kernel,
			    Cartesian_instantaneous_kernel<Active_objects_vector<ActiveObject>,
							   typename Sist_types::Static_kernel>,
			    typename Sist_types::Kinetic_kernel,
								   typename Sist_types::Simulator,
			    Active_objects_vector<ActiveObject> > P;
  Suggested_inexact_simulation_traits(const typename P::Time &lb, 
				      const typename P::Time &ub): P(lb,ub){}
};

CGAL_KDS_END_INTERNAL_NAMESPACE



#endif
