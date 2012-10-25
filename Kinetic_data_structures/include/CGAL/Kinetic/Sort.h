// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
#include <CGAL/Kinetic/listeners.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/internal/debug_counters.h>
#include <CGAL/Kinetic/Sort_visitor_base.h>
#include <CGAL/Kinetic/Event_base.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>

namespace CGAL { namespace Kinetic {

template <class KDS, class It, class RE>
class Swap_event;

struct Empty_data {};

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
  // for later, please ignore
  typedef typename Traits::Active_points_1_table TTable;
  typedef typename Traits::Kinetic_kernel::Compare_x_1 KLess;
  typedef typename Traits::Instantaneous_kernel::Compare_x_1 IComp;

  //typedef typename Traits::Instantaneous_kernel::Compare_x_1 ILess;

  typedef Sort<Traits, Visitor> This;
  // The way the Simulator represents time.
  typedef typename Traits::Simulator::Time Time;
  // The way the Simulator represents time.
  typedef typename Traits::Simulator::NT NT;
  // A label for a moving primitive in the MovingObjectTable
  typedef typename TTable::Key Object_key;


  // STL algorithms want less rather than compare. So we need to convert.
  struct ILess {
    ILess(IComp ic): ic_(ic){}

    bool operator()(Object_key a, Object_key b) const {
      bool ret=( ic_(a,b) == CGAL::SMALLER);
      return ret;
    }
    IComp ic_;
  };

  // A label for a certificate so it can be descheduled.
  typedef typename Traits::Simulator::Event_key Event_key;
  // To shorten the names. Use the default choice for the static kernel.
  typedef typename Traits::Instantaneous_kernel Instantaneous_kernel;
  // The table containing the points
  typedef TTable Active_objects_table;
  // the comparators, one for static and one for instantaneous
  typedef KLess Kinetic_less;
  typedef ILess Instantaneous_less;
  struct OD {
    OD(Object_key k): key_(k){}
    Object_key object() const {return key_;}
    Event_key event() const {return event_;}
    void set_event(Event_key k) {
      event_= k;
    }
    operator Object_key() const {
      return key_;
    }
    void swap(OD &o) {
      std::swap(key_, o.key_);
    }
    Object_key key_;
    Event_key event_;
  };
  // this is used to identify pairs of objects in the list
  typedef typename std::list<OD>::iterator iterator;
  // The certificate generator
  /*struct Less {
    typedef typename Traits::Kinetic_kernel::Is_less_x_1 Less_x;
    Less(Less_x x): less_(x){}
    bool operator()(const OD &a, const OD &b) const {
      return less_(a.key(), b.key());
    }
    Less less_;
    }*/

  typedef Swap_event<This,iterator, typename KLess::result_type> Event;
  friend class Swap_event<This,iterator, typename KLess::result_type>;
  // Redirects the Simulator notifications to function calls
  CGAL_KINETIC_DECLARE_LISTENERS(typename Traits::Simulator,
				 typename Active_objects_table)
public:

  // Register this KDS with the MovingObjectTable and the Simulator
  Sort(Traits tr, Visitor v=Visitor()/*,
       typedef Active_objects_table::Handle aot,
       Kinetic_less kless=tr.kinetic_kernel_object().is_less_x_1_object(),
       Instantaneous_less iless*/): compare_(tr.kinetic_kernel_object().compare_x_1_object()),
				    ik_(tr.instantaneous_kernel_object()),
				    iless_(ik_.compare_x_1_object()), v_(v),
				    aot_(tr.active_points_1_table_handle()),
				    simulator_(tr.simulator_handle()){
    CGAL_KINETIC_INITIALIZE_LISTENERS(simulator_, aot_);
    wrote_objects_= false;
  }

  const Visitor& visitor() const {
    return v_;
  }
  Visitor& visitor() {
    return v_;
  }

  /*Traits &traits() {
    return tr_;
  }
  const Traits &traits() const {
    return tr_;
    }*/

  typedef iterator Vertex_handle;

  /* Insert k and update the affected certificates. std::upper_bound
     returns the first place where an item can be inserted in a sorted
     list. Called by the MOT_listener.*/
  iterator insert(Object_key k) {
    NT nt= simulator_->next_time_representable_as_nt();
    simulator_->set_current_time(nt);
    ik_.set_time(nt);
    iterator it = std::upper_bound(sorted_.begin(), sorted_.end(),
				   k, iless_);
    CGAL_LOG(Log::LOTS, "\nInserting " << k);
    if (it != sorted_.end()) {
      CGAL_LOG(Log::LOTS, " before " << it->object() <<std::endl);
    } else {
      CGAL_LOG(Log::LOTS, " before end" <<std::endl);
    }
    /*if (it != sorted_.end()) {
      v_.remove_edge(prior(it), it);
      }*/
    v_.pre_insert_vertex(k);
    sorted_.insert(it, OD(k));


    rebuild_certificate(prior(it));
    //v_.create_edge(prior(it), it);
    if (prior(it) != sorted_.begin()) {
      rebuild_certificate(prior(prior(it)));
      //v_.create_edge(prior(prior(it)), prior(it));
    }

    v_.post_insert_vertex(prior(it));
    CGAL_LOG_WRITE(Log::LOTS, write(LOG_STREAM));
    CGAL_LOG(Log::LOTS, std::endl);
    return it;
    //write(std::cout);
  }

  /* Rebuild the certificate for the pair of points *it and *(++it).
     If there is a previous certificate there, deschedule it.*/
  void rebuild_certificate(const iterator it) {
    CGAL_LOG(Log::LOTS, "Building certifiate for " << it->object() << " and " << next(it)->object()<< std::endl);
    CGAL_precondition(it != sorted_.end());
    if (it->event() != Event_key()) {
      simulator_->delete_event(it->event());
      it->set_event(Event_key());
    }
    if (next(it) == sorted_.end()) return;
    //Less less=kk_.less_x_1_object();
    typename KLess::result_type s
      = compare_( aot_->at(next(it)->object()), aot_->at(it->object()), simulator_->current_time(),
		 simulator_->end_time());
    // the Simulator will detect if the failure time is at infinity
    if (s.will_fail()) {
      Time t= s.failure_time();
      s.pop_failure_time();
      Event e(it, this,s);
      it->set_event( simulator_->new_event(t, e));
    } else {
      it->set_event( simulator_->null_event());
    }
    //} else events_[*it]= simulator_->null_event();
  }



  /* Swap the pair of objects with *it as the first element.  The old
     solver is used to compute the next root between the two points
     being swapped. This method is called by an Event object.*/
  void swap(iterator it, typename KLess::result_type &s) {
    CGAL_LOG(Log::LOTS, "Swapping " << it->object() << " and " << next(it)->object() << std::endl);
    CGAL_LOG_WRITE(Log::LOTS, write(LOG_STREAM));
    v_.pre_swap(it, next(it));
    it->set_event(Event_key());
    iterator n= next(it);
    if (n->event() != Event_key()) {
      simulator_->delete_event(n->event());
      n->set_event(Event_key());
    }

    it->swap(*n);

    CGAL_LOG(Log::LOTS, "Updating next certificate " << std::endl);
    if (n != sorted_.end()) {
      rebuild_certificate(n);
    }
    CGAL_LOG(Log::LOTS, "Updating middle certificate " << std::endl);
    if (s.will_fail()) {
      Time t= s.failure_time(); s.pop_failure_time();
      it->set_event(simulator_->new_event(t, Event(it, this,s)));
    } else {
      it->set_event(simulator_->null_event());
    }


    CGAL_LOG(Log::LOTS, "Updating prev certificate " << std::endl);
    if (it != sorted_.begin()) {
      rebuild_certificate(prior(it));
    }
    v_.post_swap(it, n);

    CGAL_LOG_WRITE(Log::LOTS, write(LOG_STREAM));
  }

  void audit_order() const {
    //std::cout << "Auditing order at time " << ik_.time() << std::endl;
    for (typename std::list<OD>::const_iterator it
	   = sorted_.begin(); it->object() != sorted_.back().object(); ++it) {
      if (iless_(*next(it), *it)) {
#ifdef CGAL_KINETIC_CHECK_EXACTNESS
	std::cerr << "ERROR: objects " << it->object() << " and "
		  << next(it)->object() << " are out of order.\n";
	std::cerr << "Kinetic are " << aot_->at(it->object()) << " and " << aot_->at(*next(it)) << std::endl;
	std::cerr << "Time is " << ik_.time() << std::endl;
	/*std::cerr << "Static are " << ik_.current_coordinates_object()(it->object()) << " and "
	  << ik_.current_coordinates_object()(next(it)->object()) << std::endl;*/
	std::cerr << "ERROR: order is ";
#else
	if (warned_.find(*it) == warned_.end() ||
	    warned_[*it].find(*next(it)) == warned_[*it].end()) {
	  std::cerr << "NUMERIC ISSUE: objects " << it->object() << " and "
		    << next(it)->object() << " are out of order.\n";
	  std::cerr << aot_->at(it->object()) << " and " << aot_->at(next(it)->object()) << std::endl;
	  std::cerr << "Time is " << ik_.time() << std::endl;
	  std::cerr << "NUMERIC ISSUE: order is ";
	}
#endif
	write(std::cerr);
	std::cerr << std::endl;
	++internal::audit_failures__;

	if (!wrote_objects_) {
	  wrote_objects_=true;
	  std::cerr << "Objects are: ";
	  for (typename Active_objects_table::Key_iterator kit= aot_->keys_begin();
	       kit != aot_->keys_end(); ++kit){
	    std::cerr <<  aot_->at(*kit) << std::endl;
	  }
	}
      }
      if (compare_.sign_at(  aot_->at(it->object()),
                             aot_->at(next(it)->object()),
                             simulator_->current_time()) == CGAL::LARGER) {
#ifdef CGAL_KINETIC_CHECK_EXACTNESS
	std::cerr << "ERROR: kinetic objects " << it->object() << " and "
		  << next(it)->object() << " are out of order.\n";
	std::cerr << "Kinetic are " << aot_->at(it->object()) << " and " << aot_->at(*next(it)) << std::endl;
	std::cerr << "Time is " <<ik_.time() << std::endl;
	/*std::cerr << "Static are " << ik_.current_coordinates_object()(it->object()) << " and "
	  << ik_.current_coordinates_object()(next(it)->object()) << std::endl;*/
	std::cerr << "ERROR: order is ";
#else
	if (warned_.find(*it) == warned_.end() ||
	    warned_[*it].find(*next(it)) == warned_[*it].end()) {
	  std::cerr << "NUMERIC ISSUE: objects " << it->object() << " and "
		    << next(it)->object() << " are out of order.\n";
	  std::cerr << aot_->at(it->object()) << " and " << aot_->at(next(it)->object()) << std::endl;
	  std::cerr << "Time is " <<ik_.time() << std::endl;
	  std::cerr << "NUMERIC ISSUE: order is ";
	}
#endif
	write(std::cerr);
	std::cerr << std::endl;
	++internal::audit_failures__;

	if (!wrote_objects_) {
	  wrote_objects_=true;
	  std::cerr << "Objects are: ";
	  for (typename Active_objects_table::Key_iterator kit= aot_->keys_begin();
	       kit != aot_->keys_end(); ++kit){
	    std::cerr << aot_->at(*kit) << std::endl;
	  }
	}
      }
    }
  }

  /* Verify the structure by checking that the current coordinates are
     properly sorted for time t. This function is called by the Sim_listener.*/
  void audit() const
  {

    if (sorted_.size() <2) return;

    ik_.set_time(simulator_->audit_time());
    CGAL_LOG_WRITE(Log::LOTS, write(LOG_STREAM));
    CGAL_LOG(Log::LOTS, std::endl);


    //typename Instantaneous_kernel::Less_x_1 less= ik_.less_x_1_object();
    for (typename std::list<OD>::const_iterator it
	   = sorted_.begin(); it->object() != sorted_.back().object(); ++it) {
      CGAL_assertion(it->event() != Event_key());
    }
    CGAL_assertion(sorted_.back().event()==Event_key());

    audit_order();
  }

  /* Update the certificates adjacent to object k. This method is called by
     the MOT_listener. std::equal_range finds all items equal
     to a key in a sorted list (there can only be one).*/
  void set(Object_key k) {
    typename std::list<OD>::iterator it;
    for (it = sorted_.begin(); it != sorted_.end(); ++it){
      if (it->object()==k) break;
    }
    CGAL_assertion(it != sorted_.end());
    v_.change_vertex(it);
    rebuild_certificate(it);
    if (it != sorted_.begin()) rebuild_certificate(--it);
  }

  /* Remove object k and destroy 2 certificates and create one new one.
     This function is called by the MOT_listener.*/
  void erase(Object_key k) {
    iterator it;
    for (it = sorted_.begin(); it != sorted_.end(); ++it){
      if (it->object()==k) break;
    }
    //iterator it =  std::equal_range(sorted_.begin(), sorted_.end(),k).first;
    CGAL_precondition(it != sorted_.end());
    CGAL_precondition(it->object() == k);

    v_.pre_remove_vertex(it);
    if (next(it) != Iterator(end())) {
      simulator_->delete_event(it->event());
      it->set_event(Event_key());
    }
    if (it != sorted_.begin()) {
      iterator p= prior(it);
      sorted_.erase(it);
      rebuild_certificate(p);
    } else {
      sorted_.erase(it);
    }
    v_.post_remove_vertex(k);
  }
  template <class It> static It next(It it){ return ++it;}
  template <class It> static It prior(It it){ return --it;}

  void write(std::ostream &out) const {
    out << "Sort:\n";
    for (typename std::list<OD>::const_iterator it
	   = sorted_.begin(); it != sorted_.end(); ++it) {
      out << it->object() << " with event (";
      if (it->event() != Event_key()) {
        out << it->event();
      } else {
        out << "NULL";
      }
      out << ")\n";
    }
    out << std::endl << std::endl;;
  }

  typedef typename std::list<OD>::const_iterator Iterator;
  Iterator begin() const
  {
    return sorted_.begin();
  }
  Iterator end() const
  {
    return sorted_.end();
  }

  // The points in sorted order
  std::list<OD > sorted_;
  // events_[k] is the certificates between k and the object after it
  //std::map<Object_key, Event_key > events_;
  Kinetic_less compare_;
  Instantaneous_kernel ik_;
  Instantaneous_less iless_;
  //#ifndef NDEBUG
  mutable bool wrote_objects_;
  mutable std::map<Object_key, std::set<Object_key> > warned_;
  Visitor v_;
  typename Active_objects_table::Handle aot_;
  typename Traits::Simulator::Handle simulator_;
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
    out << left_object_->object() << "X" << Sort::next(left_object_)->object();
    if (s_.will_fail()) out <<  " next is " << s_.failure_time();
    else out << " out of failures";
  }
  void audit(typename Sort::Event_key
#ifndef NDEBUG
	     tk
#endif
) const {
    //std::cout << "Auditing event ";
    //write(std::cout);
    //std::cout << std::endl;
    CGAL_assertion(left_object_->event() == tk);
  }
  Id left_object_; Solver s_;
};


} } //namespace CGAL::Kinetic
#endif
