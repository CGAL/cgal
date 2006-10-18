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
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_TESTING_SORT_H
#define CGAL_KINETIC_TESTING_SORT_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Active_objects_listener_helper.h>
#include <CGAL/Kinetic/Cartesian_instantaneous_kernel.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/Simulator_kds_listener.h>
#include <CGAL/Kinetic/Sort_visitor_base.h>
#include <CGAL/Kinetic/Event_base.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>

CGAL_KINETIC_BEGIN_NAMESPACE

template <class KDS, class It, class RE>
class Swap_event;

//! A simple KDS which maintains objects sorted by their x coordinate
/*!  This does not use Simple_kds_base for now irrelevant
  reasons. Ditto for the lack of protection of any of the fields. The
  code is designed to be read, so read it if you want to figure out
  what is going on.
*/
template <class Traits, class Visitor=Sort_visitor_base> class Sort:
// for ref counted pointers
  public Ref_counted<Sort< Traits, Visitor> >
{
  typedef Sort<Traits, Visitor> This;
  // The way the Simulator represents time.
  typedef typename Traits::Simulator::Time Time;
  // A label for a moving primitive in the MovingObjectTable
  typedef typename Traits::Active_points_1_table::Key Object_key;
  // A label for a certificate so it can be descheduled.
  typedef typename Traits::Simulator::Event_key Event_key;
  // To shorten the names. Use the default choice for the static kernel.
  typedef typename Traits::Instantaneous_kernel Instantaneous_kernel;
  // this is used to identify pairs of objects in the list
  typedef typename std::list<Object_key>::iterator iterator;
  // The certificate generator
  typedef typename Traits::Kinetic_kernel::Is_less_x_1 Less;
  typedef Swap_event<This,iterator,typename Less::result_type> Event;
  // Redirects the Simulator notifications to function calls
  typedef typename CGAL::Kinetic::
  Simulator_kds_listener<typename Traits::Simulator::Listener,
			 This> Sim_listener;
  // Redirects the MovingObjectTable notifications to function calls
  typedef typename CGAL::Kinetic::
  Active_objects_listener_helper<typename Traits::Active_points_1_table::Listener,
				 This> MOT_listener;
public:
  // Register this KDS with the MovingObjectTable and the Simulator
  Sort(Traits tr, Visitor v=Visitor()): less_(tr.kinetic_kernel_object().is_less_x_1_object()),
					ik_(tr.instantaneous_kernel_object()), v_(v),
					tr_(tr){
    sim_listener_= Sim_listener(tr.simulator_handle(), this);
    mot_listener_= MOT_listener(tr.active_points_1_table_handle(), this);
    wrote_objects_= false;
  }

  const Visitor& visitor() const {
    return v_;
  }
  Visitor& visitor() {
    return v_;
  }

  Traits &traits() {
    return tr_;
  }
  const Traits &traits() const {
    return tr_;
  }

  typedef iterator Vertex_handle;

  /* Insert k and update the affected certificates. std::upper_bound
     returns the first place where an item can be inserted in a sorted
     list. Called by the MOT_listener.*/
  void insert(Object_key k) {
    //std::cout << "Inserting " << k <<std::endl;
    ik_.set_time(simulator()->rational_current_time());
    iterator it = std::upper_bound(sorted_.begin(), sorted_.end(),
				   k, ik_.less_x_1_object());
    /*if (it != sorted_.end()) {
      v_.remove_edge(prior(it), it);
      }*/
    sorted_.insert(it, k);
    v_.create_vertex(prior(it));
    
    rebuild_certificate(prior(it)); 
    //v_.create_edge(prior(it), it);
    if (prior(it) != sorted_.begin()) {
      rebuild_certificate(prior(prior(it)));
      //v_.create_edge(prior(prior(it)), prior(it));
    }
    //write(std::cout);
  }

  /* Rebuild the certificate for the pair of points *it and *(++it).
     If there is a previous certificate there, deschedule it.*/
  void rebuild_certificate(const iterator it) {
    CGAL_precondition(it != sorted_.end());
    if (events_.find(*it) != events_.end()) {
      simulator()->delete_event(events_[*it]); events_.erase(*it);
    }
    if (next(it)== sorted_.end()) return;
    //Less less=kk_.less_x_1_object();
    typename Less::result_type s
      = less_(object(*(it)), object(*next(it)), simulator()->current_time(),
	      simulator()->end_time());
    // the Simulator will detect if the failure time is at infinity
    Time t= s.failure_time();
    s.pop_failure_time();
    Event e(it, this,s);
    events_[*it]= simulator()->new_event(t, e);
    //} else events_[*it]= simulator()->null_event();
  }



  /* Swap the pair of objects with *it as the first element.  The old
     solver is used to compute the next root between the two points
     being swapped. This method is called by an Event object.*/
  void swap(iterator it, typename Less::result_type &s) {
    v_.before_swap(it, next(it));
    events_.erase(*it);
    iterator n= next(it);
    if (*n!= sorted_.back()) simulator()->delete_event(events_[*n]);
    events_.erase(*next(it));

    std::swap(*it, *next(it));
    if (next(it) != sorted_.end()) {
      rebuild_certificate(next(it));
      //v_.create_edge(next(it)), next(next(it));
    }
    Time t= s.failure_time(); s.pop_failure_time();
    events_[*it]= simulator()->new_event(t, Event(it, this,s));
    
    v_.after_swap(it, next(it));
    if (it != sorted_.begin()) {
      rebuild_certificate(--it);
      //v_.create_edge(it, next(it));
    }
    
    //write(std::cout);
    //std::cout << "At time " << simulator()->current_time() << std::endl;
  }

  /* Verify the structure by checking that the current coordinates are
     properly sorted for time t. This function is called by the Sim_listener.*/
  void audit() const
  {

    if (sorted_.size() <2) return;

    ik_.set_time(simulator()->rational_current_time());
    typename Instantaneous_kernel::Less_x_1 less= ik_.less_x_1_object();
    
    for (typename std::list<Object_key>::const_iterator it
	   = sorted_.begin(); *it != sorted_.back(); ++it) {
#ifdef CGAL_KINETIC_CHECK_EXACTNESS
      if (!less(*it, *next(it))) {
	std::cerr << "ERROR: objects " << *it << " and "
		  << *next(it) << " are out of order.\n";
	std::cerr << object(*it) << " and " << object(*next(it)) << std::endl;
	std::cerr << "Time is " <<simulator()->rational_current_time() << std::endl; 
	std::cerr << "ERROR: order is ";
	write(std::cerr);
	std::cerr << std::endl;
	CGAL::Kinetic::internal::fail__=true;

	if (!wrote_objects_) {
	  wrote_objects_=true;
	  std::cerr << "Objects are: ";
	  for (typename Traits::Active_objects_table::Keys_iterator kit= mot_listener_.notifier()->keys_begin();
	       kit != mot_listener_.notifier()->keys_end(); ++kit){
	    std::cerr <<  mot_listener_.notifier()->at(*kit) << std::endl;
	  }
	}
      }
#else
      if (!less(*it, *next(it))) {
	if (warned_.find(*it) == warned_.end() ||
	    warned_[*it].find(*next(it)) == warned_[*it].end()) {
	  std::cerr << "WARNING: objects " << *it << " and "
		    << *next(it) << " are out of order.\n";
	  std::cerr << object(*it) << " and " << object(*next(it)) << std::endl;
	  std::cerr << "Time is " <<simulator()->rational_current_time() << std::endl; 
	  std::cerr << "WARNING: order is ";
	  write(std::cerr);
	  if (!wrote_objects_) {
	    wrote_objects_=true;
	    std::cerr << "Objects are: ";
	    for (typename Traits::Active_points_1_table::Key_iterator kit= mot_listener_.notifier()->keys_begin();
		 kit != mot_listener_.notifier()->keys_end(); ++kit){
	      std::cerr <<  mot_listener_.notifier()->at(*kit) << std::endl;
	    }
	  }
	  warned_[*it].insert(*next(it));
	  std::cerr << std::endl;
	}
      }	
#endif
    }
  }

  /* Update the certificates adjacent to object k. This method is called by
     the MOT_listener. std::equal_range finds all items equal
     to a key in a sorted list (there can only be one).*/
  void set(Object_key k) {
    iterator it =  std::equal_range(sorted_.begin(), sorted_.end(),k).first;
    v_.modify_vertex(it);
    rebuild_certificate(it);
    if (it != sorted_.begin()) rebuild_certificate(--it);
  }

  /* Remove object k and destroy 2 certificates and create one new one.
     This function is called by the MOT_listener.*/
  void erase(Object_key k) {
    iterator it;
    for (it = sorted_.begin(); it != sorted_.end(); ++it){
      if (*it==k) break;
    }
    //iterator it =  std::equal_range(sorted_.begin(), sorted_.end(),k).first;
    CGAL_precondition(it != sorted_.end());
    CGAL_precondition(*it == k);
    iterator p= it; --p;
    v_.remove_vertex(it);
    sorted_.erase(it);
    if (p != sorted_.end()) {
      rebuild_certificate(p);
    }
    simulator()->delete_event(events_[k]);
    events_.erase(k);
  }
  template <class It> static It next(It it){ return ++it;}
  template <class It> static It prior(It it){ return --it;}
  typename Traits::Active_points_1_table::Data object(Object_key k) const
  {
    return mot_listener_.notifier()->at(k);
  }
  typename Traits::Simulator* simulator() {return sim_listener_.notifier();}
  const typename Traits::Simulator* simulator() const {return sim_listener_.notifier();}

  void write(std::ostream &out) const
  {
    for (typename std::list<Object_key>::const_iterator it
	   = sorted_.begin(); it != sorted_.end(); ++it) {
      out << *it << " ";
    }
    out << std::endl;
  }

  typedef typename std::list<Object_key>::const_iterator Key_iterator;
  Key_iterator begin() const
  {
    return sorted_.begin();
  }
  Key_iterator end() const
  {
    return sorted_.end();
  }

  Sim_listener sim_listener_;
  MOT_listener mot_listener_;
  // The points in sorted order
  std::list<Object_key> sorted_;
  // events_[k] is the certificates between k and the object after it
  std::map<Object_key, Event_key > events_;
  Less less_;
  Instantaneous_kernel ik_;
  //#ifndef NDEBUG
  mutable bool wrote_objects_;
  mutable std::map<Object_key, std::set<Object_key> > warned_;
  Visitor v_;
  Traits tr_;
  //#endif
};

template <class Traits, class Visitor>
std::ostream &operator<<(std::ostream &out, const Sort<Traits, Visitor> &s) {
  s.write(out);
  return out;
}

/* It needs to implement the time() and process() functions and
   operator<< */
template <class Sort, class Id, class Solver>
class Swap_event: public Event_base<Sort*>
{
public:
  Swap_event(Id o, Sort* sorter,
	     const Solver &s): Event_base<Sort*>(sorter),
			       left_object_(o), s_(s){}
  void process() {
    Event_base<Sort*>::kds()->swap(left_object_, s_);
  }
  void write(std::ostream &out) const {
    out << *left_object_ << "X" << *Sort::next(left_object_);
  }
  Id left_object_; Solver s_;
};


CGAL_KINETIC_END_NAMESPACE
#endif
