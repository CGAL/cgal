#ifndef CGAL_KDS_TESTING_SORT_H
#define CGAL_KDS_TESTING_SORT_H
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/Cartesian_instantaneous_kernel.h>
#include <algorithm>
#include <map>
#include <list>
#include <iterator>
#include <iostream>
#include <CGAL/KDS/Ref_counted.h>
#include <CGAL/KDS/Notifying_table_listener_helper.h>
#include <CGAL/KDS/Simulator_kds_listener.h>

CGAL_KDS_BEGIN_NAMESPACE

template <class KDS, class It, class RE> 
class Swap_event;

//! A simple KDS which maintains objects sorted by their x coordinate
/*!  This does not use Simple_kds_base for now irrelevant
  reasons. Ditto for the lack of protection of any of the fields. The
  code is designed to be read, so read it if you want to figure out
  what is going on.
*/
template <class Traits> class Sort:
   // for ref counted pointers
  public Ref_counted<Sort< Traits> >  {
  typedef Sort<Traits> This;
  // The way the Simulator represents time.
  typedef typename Traits::Simulator::Time Time;
  // A label for a moving primitive in the MovingObjectTable
  typedef typename Traits::Moving_point_table::Key Object_key;
  // A label for a certificate so it can be descheduled.
  typedef typename Traits::Simulator::Event_key Event_key;
  // To shorten the names. Use the default choice for the static kernel.
  typedef typename Traits::Instantaneous_kernel Instantaneous_kernel;
  // this is used to identify pairs of objects in the list
  typedef typename std::list<Object_key>::iterator iterator;
  typedef Swap_event<This,iterator,typename Traits::Simulator::Root_stack> Event;
  // Redirects the Simulator notifications to function calls
  typedef typename CGAL::KDS::
  Simulator_kds_listener<typename Traits::Simulator::Listener, 
				This> Sim_listener;
  // Redirects the MovingObjectTable notifications to function calls
  typedef typename CGAL::KDS::
  Notifying_table_listener_helper<typename Traits::Moving_point_table::Listener,
				      This> MOT_listener;
public:
   // Register this KDS with the MovingObjectTable and the Simulator
  Sort(Traits tr): sim_listener_(tr.simulator_pointer(), this),
		   mot_listener_(tr.moving_point_table_pointer(), this),
		   kk_(tr.kinetic_kernel_object()),
		   ik_(tr.instantaneous_kernel_object()){}

  /* Insert k and update the affected certificates. std::upper_bound
     returns the first place where an item can be inserted in a sorted
     list. Called by the MOT_listener.*/
  void insert(Object_key k) {
    //std::cout << "Inserting " << k <<std::endl;
    ik_.set_time(simulator()->rational_current_time());
    iterator it = std::upper_bound(sorted_.begin(), sorted_.end(),
				   k, ik_.less_x_1_object());
    sorted_.insert(it, k);
    rebuild_certificate(--it); rebuild_certificate(--it);
    //write(std::cout);
  }

  /* Rebuild the certificate for the pair of points *it and *(++it).
     If there is a previous certificate there, deschedule it.*/
  void rebuild_certificate(const iterator it) {
    if (it == sorted_.end()) return;
    if (events_.find(*it) != events_.end()) {
      simulator()->delete_event(events_[*it]); events_.erase(*it);
    }
    if (next(it)== sorted_.end()) return;     
    typename Traits::Kinetic_kernel::Less_x_1 less=kk_.less_x_1_object();
    typename Traits::Simulator::Root_stack s 
      = simulator()->root_stack_object(less(object(*(it)), 
					    object(*next(it))));
    // the Simulator will detect if the failure time is at infinity
    if (!s.empty()) {
      Time t= s.top(); 
      s.pop();
      Event e(it, this,s);
      events_[*it]= simulator()->new_event(t, e);
    } else events_[*it]= simulator()->null_event();
  }

  /* Swap the pair of objects with *it as the first element.  The old
     solver is used to compute the next root between the two points
     being swapped. This method is called by an Event object.*/
  void swap(iterator it, typename Traits::Simulator::Root_stack &s) {
    events_.erase(*it);
    iterator n= next(it);
    if (*n!= sorted_.back()) simulator()->delete_event(events_[*n]);
    events_.erase(*next(it));
    std::swap(*it, *next(it));
    rebuild_certificate(next(it));
    if (!s.empty()){
      Time t= s.top(); s.pop(); 
      events_[*it]= simulator()->new_event(t, Event(it, this,s));
    } else events_[*it]= simulator()->null_event();
    if (it != sorted_.begin()) rebuild_certificate(--it);
    //write(std::cout);
    //std::cout << "At time " << simulator()->current_time() << std::endl;
  }

  /* Verify the structure by checking that the current coordinates are
     properly sorted for time t. This function is called by the Sim_listener.*/
  void audit() const {
    if (sorted_.size() <2) return;
    ik_.set_time(simulator()->rational_current_time());
    typename Instantaneous_kernel::Less_x_1 less= ik_.less_x_1_object();
    for (typename std::list<Object_key>::const_iterator it
	   = sorted_.begin(); *it != sorted_.back(); ++it){
      if (!less(*it, *next(it))){
	std::cerr << "ERROR: objects " << *it << " and " 
		  << *next(it) << " are out of order.\n";
	std::cerr << "ERROR: order is ";
	write(std::cerr);
	std::cerr << std::endl;
      }
    }
  }

  /* Update the certificates adjacent to object k. This method is called by
     the MOT_listener. std::equal_range finds all items equal 
     to a key in a sorted list (there can only be one).*/
  void set(Object_key k) {
    iterator it =  std::equal_range(sorted_.begin(), sorted_.end(),k).first;
    rebuild_certificate(it); rebuild_certificate(--it);
  }

  /* Remove object k and destroy 2 certificates and create one new one.
     This function is called by the MOT_listener.*/
  void erase(Object_key k) {
    iterator it =  std::equal_range(sorted_.begin(), sorted_.end(),k).first;
    iterator p= it; --p;
    sorted_.erase(it);
    rebuild_certificate(p);
    simulator()->delete_event(events_[k]);
    events_.erase(k);
  }
  template <class It> static It next(It it){ return ++it;}
  typename Traits::Moving_point_table::Data object(Object_key k) const { 
    return mot_listener_.notifier()->at(k);
  }
  typename Traits::Simulator* simulator() {return sim_listener_.notifier();}
  const typename Traits::Simulator* simulator() const {return sim_listener_.notifier();}

  void write(std::ostream &out) const {
    for (typename std::list<Object_key>::const_iterator it
	   = sorted_.begin(); it != sorted_.end(); ++it){
      out << *it << " "; 
    }
    out << std::endl;
  }

  typedef typename std::list<Object_key>::const_iterator Key_iterator;
  Key_iterator begin() const {
    return sorted_.begin();
  }
  Key_iterator end() const {
    return sorted_.end();
  }

  Sim_listener sim_listener_; 
  MOT_listener mot_listener_; 
  // The points in sorted order
  std::list<Object_key> sorted_;
  // events_[k] is the certificates between k and the object after it
  std::map<Object_key, Event_key > events_;
  typename Traits::Kinetic_kernel kk_; 
  Instantaneous_kernel ik_;
};

/* It needs to implement the time() and process() functions and 
   operator<< */
template <class Sort, class Id, class Solver> 
class Swap_event {
public:
  Swap_event(Id o, Sort* sorter, 
	     const Solver &s): left_object_(o), sorter_(sorter), s_(s){}
  void process(const typename Solver::Root &){
    sorter_->swap(left_object_, s_);
  }
  Id left_object_; Sort* sorter_; Solver s_;
};
template <class S, class I, class SS>
std::ostream &operator<<(std::ostream &out,
			 const Swap_event<S,I,SS> &ev){
  return out << "swap " << *ev.left_object_ ;
}


CGAL_KDS_END_NAMESPACE
#endif
  
