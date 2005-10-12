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

#ifndef CGAL_KDS_BIN_QUEUE_H
#define CGAL_KDS_BIN_QUEUE_H
#include <CGAL/KDS/basic.h>
#include <iostream>
#include <CGAL/In_place_list.h>
#include <functional>
#include <CGAL/assertions.h>
#include <iostream>
#include <CGAL/KDS/Ref_counted.h>
#include <CGAL/KDS/internal/infinity_or_max.h>
#include <algorithm>
#include <boost/utility.hpp>

//int two_list_remaining=0;

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE



// The interface for an item stored in the ::Pointer_event_queue
template <class Priority>
class Two_list_event_queue_item: 
  public CGAL::In_place_list_base<Two_list_event_queue_item<Priority> >,
  public Ref_counted<Two_list_event_queue_item<Priority> > {

  typedef Two_list_event_queue_item<Priority> This;

  Two_list_event_queue_item(const Two_list_event_queue_item &o){bool do_not_copy;}
  void operator=(const Two_list_event_queue_item &o) {
    bool do_not_copy;
  }
public:
  Two_list_event_queue_item(): time_(infinity_or_max<Priority>()){ /*++two_list_remaining;*/}
  Two_list_event_queue_item(const Priority &t): time_(t){/*++two_list_remaining;*/}
  virtual ~Two_list_event_queue_item(){/*--two_list_remaining;*/}

  const Priority& time() const {return time_;};

  bool is_in_front_list() const {
    return front_list_;
  }
  void set_is_in_front_list(bool ft) {
    front_list_=ft;
  }

  virtual void write(std::ostream &out) const =0;
  virtual void process(const Priority &t) =0;
 
  bool operator<(const This &o) const {
    return time() < o.time();
  }
  bool operator>(const This &o) const {
    return time() > o.time();
  }
  bool operator==(const This &o) const {
    return time() == o.time();
  }
private:
  Priority time_;
  bool front_list_;
};

template <class Priority>
inline std::ostream& operator<<(std::ostream &out, const Two_list_event_queue_item<Priority> &i) {
  i.write(out);
  return out;
}


// The how a dummy item is stored in the ::Two_list_event_queue
/*
  One dummy item is used to represent events which will never happen.
*/
template <class Priority>
class Two_list_event_queue_dummy_item: public Two_list_event_queue_item<Priority> {

public:
  Two_list_event_queue_dummy_item(){}
  Two_list_event_queue_dummy_item(const Two_list_event_queue_dummy_item &):
    Two_list_event_queue_item<Priority>(){}
  virtual void process(const Priority &){
    std::cerr << "Trying to process a NULL event.\n";
    assert(0);
  }
  virtual void write(std::ostream &out) const {
    out << "Never.";
  }
  virtual ~Two_list_event_queue_dummy_item(){}
};

// The how a real item is stored in the ::Two_list_event_queue
/*
  This just stores an object of type Event and passes the virtual calls on to it. 

  The object is reference counted so you don't have to worry about the
  queue deleting it or not.
*/
template <class Priority, class Event>
class Two_list_event_queue_item_rep: public internal::Two_list_event_queue_item<Priority> {
  typedef  Two_list_event_queue_item<Priority> P;
public:
  Two_list_event_queue_item_rep(): internal::Two_list_event_queue_item<Priority>(){}
  Two_list_event_queue_item_rep(const Priority &t, Event &e): Two_list_event_queue_item<Priority>(t), 
								event_(e){}
  
  virtual void write(std::ostream &out) const {
    out << event_ << " at " << P::time();
  }
  virtual void process(const Priority &t){
    event_.process(t);
  }
  // Access the actual event
  const Event &event() const {
    return event_;
  }
  void set_event(const Event &e) {
    event_=e;
  }
  virtual ~Two_list_event_queue_item_rep(){}

protected:
  // what the fuck? Why do I need this (gcc 3.4.1)
  mutable Event event_;
};

template <class T>
struct Two_list_event_queue_item_allocator {
  typedef Two_list_event_queue_dummy_item<T> dummy_value_type;

  typedef Two_list_event_queue_item<T>* pointer;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef const pointer const_pointer;
  typedef Two_list_event_queue_item<T>& reference;
  typedef const Two_list_event_queue_item<T>& const_reference;
  typedef dummy_value_type value_type;

  Two_list_event_queue_item_allocator(){}
  pointer allocate(size_t sz, void* =0) const {
    return static_cast<pointer>(::operator new(sz*sizeof(dummy_value_type)));
  }
  void deallocate(pointer p, size_type ){
    ::operator delete(static_cast<void*>(p));
  }
  size_type max_size() const throw() {
    return std::numeric_limits<size_type>::max()/sizeof(dummy_value_type);
  }
  template <class TT>
  void construct(pointer p, const TT &){
    new(static_cast<void*>(p)) dummy_value_type();
  }
  void destroy(pointer p){
    p->~Two_list_event_queue_item<T>();
  }
};


CGAL_KDS_END_INTERNAL_NAMESPACE


CGAL_KDS_BEGIN_NAMESPACE;

template <class PriorityT, class NT>
class Two_list_pointer_event_queue;

template <class Item, class PriorityT, class NT>
struct Two_list_pointer_event_queue_key: private Item::Pointer{
  typedef Two_list_pointer_event_queue_key<Item, PriorityT, NT> This;
  typedef typename Item::Pointer P;
  Two_list_pointer_event_queue_key(){};
  Two_list_pointer_event_queue_key(Item  *p): Item::Pointer(p){};
  friend class Two_list_pointer_event_queue<PriorityT, NT>;
  void write(std::ostream &out){
    if (Item::Pointer::get()) {
      out << *Item::Pointer::get();
    } else {
      out << "null";
    }
  }
  operator bool() const {
    return Item::Pointer::get() != NULL;
  }
  bool operator!() const {
    return Item::Pointer::operator!();
  }
  bool operator==(const This &o) const {
    return Item::Pointer::get() == o.get();
  }
  bool operator!=(const This &o) const {
    return Item::Pointer::get() != o.get();
  }
  //using P::operator<; 
  bool operator<(const This& o) const {
    return P::get() < o.get();
  }

  Item* pointer() {
    return Item::Pointer::get();
  }
  //using P::operator>;
};

template  <class I, class PT, class NT>
std::ostream &operator<<(std::ostream &out, Two_list_pointer_event_queue_key<I, PT, NT> k){
  k.write(out);
  return out;
}

//! The priority queue for holding many different types of events.
/*!  This queue allows the priorities to be updated and for elements
  to be removed. The events are stored behind an interface and
  accessed through virtual functions allowing many different types of
  events to be stored in the queue at once, as long as they all share
  the same Priority.

  Currently the items in the queue are refence counted. This imposes
  some overhead, but makes accessing them much simpler since you don't
  have to worry about them being processed (and deleted or not). I am
  not sure which is better.
*/
template <class PriorityT, class NT>
class Two_list_pointer_event_queue {
  typedef Two_list_pointer_event_queue<PriorityT, NT> This;
  typedef typename internal::Two_list_event_queue_item<PriorityT> Item;
 
  typedef typename CGAL::In_place_list<Item, false,
				       internal::Two_list_event_queue_item_allocator<PriorityT> > Queue;
  typedef typename Queue::iterator Iterator;

  //static const unsigned int max_front_size=20;

public:
  typedef PriorityT Priority;
   
  typedef Two_list_pointer_event_queue_key<Item, PriorityT, NT> Key;
  
  /*  struct Key: public Item::Pointer {
    Key(){};
    Key(Item*i): Item::Pointer(i){}
    Key(typename Item::Pointer ih): Item::Pointer(ih){}
    static Key null() {
      return Key();
    }
    };*/

  //! Construct it with a suggested size of sz.
  Two_list_pointer_event_queue(Priority start_time, int =0): ub_(to_interval(start_time).first),
							     step_(1) {
    null_event_= new internal::Two_list_event_queue_dummy_item<Priority>();
  }

  //! insert value_type into the queue and return a reference to it
  /*!
    The reference can be used to update or erase it.
  */
  template <class E>
  Key insert(const Priority &t, const E & e){
    CGAL_expensive_precondition(audit());
    CGAL_precondition(t != internal::infinity_or_max<Priority>());
   
    Item *ni = make_event(t, e);
    
    //CGAL_exactness_assertion(t >= lb_);
    //lb_=std::min(t, lb_);
    
    if (t <= ub_) {
      ni->set_is_in_front_list(true);
      typename Queue::iterator iit=std::upper_bound(front_.begin(), front_.end(), *ni);
      /*if (iit == front_.begin()){
	step_= ub_- to_interval(t).first;
	}*/
      front_.insert(iit, *ni);
      CGAL_expensive_assertion(audit());
      if (front_.size() > 2*max_front_size()){
	shrink_front();
      }

    } else if (front_.empty()){
      assert(back_.empty());
      front_.push_back(*ni);
      ub_= to_interval(t).second;
      ni->set_is_in_front_list(true);
    } else {
      ni->set_is_in_front_list(false);
      back_.push_back(*ni);
    }
    CGAL_expensive_postcondition(audit());
    CGAL_expensive_postcondition(is_in_queue(Key(ni)));
    return Key(ni);
  }

  //! remove the event referenced by item from the queue
  /*!
    \todo Add check that item is actually in the list
  */
  void erase(const Key &item){
    if (item== null_event_) return;
    CGAL_expensive_precondition(is_in_queue(item));
    CGAL_expensive_precondition(audit());
    Item *i=item.get();
    if (item->is_in_front_list()){
      front_.erase(i);
      if (front_.empty() && !back_.empty()) grow_front();
    } else {
      back_.erase(i);
    }
    unmake_event(i);
    
    CGAL_expensive_postcondition(audit());
  }


  template <class E>
  const E& get(const Key &item) const {
    return reinterpret_cast<internal::Two_list_event_queue_item_rep<Priority, E>*>( item.get())->event();
  }


  //! Replace the event referenced by item with a new event referenced by ne
  /*!  They must have exactly the same times associated with
    them. This is checked when expensive checks are turned on.
  */
  template <class NE>
  Key set(const Key &item, const NE &ne) {
    CGAL_expensive_precondition(is_in_queue(item));
    CGAL_precondition(item != end_key());
    Item *oi= item.get();
    bool front= item->is_in_front_list();

    Item *ni= make_event(item->time(), ne);
    ni->set_is_in_front_list(front);

    if (front) {
      front_.insert(oi, *ni);
      front_.erase(oi);
    } else {
      back_.insert(oi, *ni);
      back_.erase(oi);
    }

    unmake_event(oi);

    CGAL_expensive_postcondition(audit());
    return Key(ni);
  }

  //! Get the time of the next event to be processed.
  /*!  It is OK to call this if the queue is empty. In that case it
    will just return infinity.
  */
  Priority front_priority() const {
    CGAL_precondition(!front_.empty());
    return front_.front().time();
  }

  //! Access the time of a particular event
  Priority priority(const Key &item) const {
    return item->time();
  }

  //! empty
  bool empty() const {
    CGAL_precondition(!front_.empty() || back_.empty());
    return front_.empty();
  }

  
  //! Remove the next event from the queue and process it.
  /*!
    Processing means that the process() method of the event object is called. 
  */
  void process_front() {
    CGAL_precondition(!empty());
    CGAL_expensive_precondition(audit());
    if (!front_.empty()){
      Item *i= &front_.front();
      CGAL_KDS_LOG(LOG_SOME, "Processing event " << *i << std::endl);
      front_.pop_front();
      CGAL_expensive_postcondition(audit());
      if (front_.empty() && !back_.empty()) grow_front();
      i->process(i->time());
      
      if (!front_.empty() && i->time() == front_.front().time()){
	CGAL_KDS_LOG(LOG_SOME, "Degeneracy at time " 
		     << i->time() << " the events are: "
		     << *i << " and " << front_.front() << std::endl);
      }

      intrusive_ptr_release(i);
    } else {
      assert(back_.empty());
    }
  }

  //! debugging
  bool print() const {
    write(std::cout);
    return true;
  }

  bool write(std::ostream &out) const {
    for (typename Queue::const_iterator it = front_.begin(); it != front_.end(); ++it){
      out << "[";
      it->write(out);
      out << "] ";
    }
    out << std::endl;
    for (typename Queue::const_iterator it = back_.begin(); it != back_.end(); ++it){
      it->write(out);
      out << " ";
    }
    out << std::endl;
    return true;
  }

  Key end_key() const {
    return null_event_;
  }
 

protected:

  template <class E>
  Item *make_event(const Priority &t, E &e){
    Item *ni= new internal::Two_list_event_queue_item_rep<Priority, E>(t, e);
    intrusive_ptr_add_ref(ni);
    return ni;
  }

  void unmake_event(Item *i){
    intrusive_ptr_release(i);
  }

  bool audit() {
    for (typename Queue::const_iterator it = front_.begin(); it != front_.end(); ++it){
      Priority t= it->time();
      CGAL_assertion(t <= ub_);
      CGAL_assertion(it->is_in_front_list());
      //CGAL_exactness_assertion(t >= lb_);
    }
    for (typename Queue::const_iterator it = back_.begin(); it != back_.end(); ++it){
      Priority t= it->time();
      CGAL_assertion(t > ub_);
      CGAL_assertion(!it->is_in_front_list());
    }
    {
      typename Queue::const_iterator it = front_.begin();
      ++it;
      for (; it != front_.end(); ++it){
	Priority tc= it->time();
	Priority tp= boost::prior(it)->time();
	CGAL_assertion(tc >= tp);
      }
    }
    return true;
  }

  bool is_in_queue(Key k){
    for (typename Queue::const_iterator it = front_.begin(); it != front_.end(); ++it){
      if (&*it == k.pointer()) return true;
    }
    for (typename Queue::const_iterator it = back_.begin(); it != back_.end(); ++it){
      if (&*it == k.pointer()) return true;
    }
    return false;
  }

  unsigned int select(Queue &source, Queue &target, NT b){
    unsigned int sz= source.size() + target.size();if (sz);
    int count=0;
    for (Iterator it= source.begin(); it != source.end(); ++it){
      // assert(it->time() >= a);
      if (it->time() <= b){
	Item *i= &*it;
	Iterator t= boost::prior(it);
	source.erase(it);
	it=t;
	target.push_back(*i);
	++count;
      }
    }
    assert(sz==source.size() + target.size());
    return count;
  }

  /*NT step() const{
    return std::max(ub_-lb_, NT(1));
    }*/

  NT av(NT a, NT b) const{
    return .5*(a+b);
  }

  template <class It>
  void set_front(It b, It e, bool val) {
    for (; b!= e; ++b){
      b->set_is_in_front_list(val);
    }
  }

  void grow_front(Queue &cand) {
    const bool dprint=false;
    CGAL_assertion(front_.empty());
    CGAL_assertion(!cand.empty());
    CGAL_assertion(step_ != 0);
    if (dprint) std::cout << "Growing front from " << ub_ << " with step " << step_;
    //lb_=ub_;
    //NT oub= ub_;
    ub_+= step_;
    unsigned int num= select(cand, front_, ub_);
    if (front_.empty()){
      //lb_=ub_;
      if (dprint) std::cout << "undershot." << std::endl;
      NT nstep = step_*2;
      if (nstep > step_) {
	step_=nstep;
	CGAL_assertion(step_!=0);
	grow_front(cand);
      }
    } else {
      //      unsigned int ncand= cand.size();
      back_.splice(back_.begin(), cand);
      if (num >  max_front_size()){
	NT nstep= step_*NT(.6+.4*max_front_size()/static_cast<double>(num));
	//else nstep = step_*.6;
	if (nstep >0) {
	  cand.swap(front_);
	  //ub_=lb_;
	  ub_-=step_;
	  if (dprint) std::cout << "...overshot" << std::endl;
	  CGAL_assertion(nstep < step_);
	  step_=nstep;
	  CGAL_assertion(step_!=0);
	  grow_front(cand);
	}
      } else {
	if (dprint) std::cout << std::endl;
      }
    }
    CGAL_postcondition(cand.empty());
  }

  void grow_front() {
    //std::cout << "Growing front from " << ub_ << " with lb " << lb_ << std::endl;
    //assert(is_valid());
    CGAL_precondition(!back_.empty());
    CGAL_precondition(front_.empty());
    CGAL_assertion_code(unsigned int sz= front_.size()+back_.size());
    Queue cand;
    cand.splice(cand.begin(), back_);
    grow_front(cand);
    set_front(front_.begin(), front_.end(), true);
    front_.sort();
    
    CGAL_assertion(sz==front_.size()+back_.size());
    CGAL_assertion(audit());
    //std::cout << "to " << ub_ << " with lb " << lb_ << std::endl;
  }

  void shrink_front() {

    typename Queue::iterator it=front_.begin();
    unsigned int mf= max_front_size();
    for (unsigned int i=0; i < mf; ++i){
      ++it;
    }

    double split= CGAL::to_interval(it->time()).second;
    while (it->time() < split && it != front_.end()) ++it;

    if (it != front_.end()) {
      //double splitubb = to_interval(it->time()).first;
      // find a simple double between split and splitubb

      set_front(it, front_.end(), false);
      back_.splice(back_.begin(), front_, it, front_.end());
      NT oub=ub_;
      //double dt=CGAL::to_interval(front_.back().time()).second;
      ub_ = NT(split);   
      step_= oub-ub_;
      CGAL_postcondition_code(if (step_<0) std::cerr << step_ << std::endl;);
      CGAL_postcondition_code(if (step_<0) std::cerr << ub_ << std::endl;);
      CGAL_postcondition_code(if (step_<0) std::cerr << oub << std::endl;);
      CGAL_postcondition_code(if (step_<0) std::cerr << front_.back().time() << std::endl;);
      CGAL_postcondition_code(if (step_==0) for (typename Queue::const_iterator it=front_.begin(); it != front_.end(); ++it) std::cout << *it << std::endl);
      CGAL_postcondition(step_>=0);
    }
  }


  unsigned int max_front_size() const {
    return std::max(10U, static_cast<unsigned int>(std::sqrt(front_.size()+back_.size())));
  }

  Queue front_, back_;
  Key null_event_;
  NT ub_, step_;
};

template <class D,class T>
std::ostream &operator<<(std::ostream &out, const Two_list_pointer_event_queue<D,T> &q){
  q.write(out);
  return out;
}

CGAL_KDS_END_NAMESPACE;


#endif
