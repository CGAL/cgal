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

#ifndef CGAL_KINETIC_BIN_QUEUE_H
#define CGAL_KINETIC_BIN_QUEUE_H
#include <CGAL/Kinetic/basic.h>
#include <iostream>
#include <CGAL/Kinetic/internal/debug_counters.h>
#include <CGAL/In_place_list.h>
#include <functional>
#include <CGAL/assertions.h>
#include <CGAL/use.h>
#include <iostream>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/internal/infinity_or_max.h>
#include <algorithm>
#include <boost/utility.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <CGAL/Kinetic/internal/debug_counters.h>

//int two_list_remaining=0;

//int growth__=0;
//int shrink__=0;
//int queue_insertions__=0;
//int queue_front_insertions__=0;

namespace CGAL { namespace Kinetic { namespace internal {


template <class Item>
struct Two_list_pointer_event_queue_key: public Item::Handle
{
  typedef Two_list_pointer_event_queue_key<Item> This;
  typedef typename Item::Handle P;
  Two_list_pointer_event_queue_key(){}
  Two_list_pointer_event_queue_key(Item  *p): Item::Handle(p){}
  std::ostream& write(std::ostream &out) {
    if (Item::Handle::get()) {
      out << "(" /*<< Item::Handle::get()*/ << ": " << *Item::Handle::get() << ")";
    }
    else {
      out << "null";
    }
    return out;
  }
  /*operator bool() const
    {
    return Item::Handle::get() != NULL;
    }*/
  bool is_valid() const {
    return  Item::Handle::get() != NULL;
  }
  /*bool operator!() const
    {
    return Item::Handle::operator!();
    }*/
  bool operator==(const This &o) const
  {
    return Item::Handle::get() == o.get();
  }
  bool operator!=(const This &o) const
  {
    return Item::Handle::get() != o.get();
  }
  //using P::operator<;
  bool operator<(const This& o) const
  {
    return P::get() < o.get();
  }

  Item* pointer() {
    return Item::Handle::get();
  }
  const Item* pointer() const
  {
    return Item::Handle::get();
  }
 
  //using P::operator>;
};

template  <class I>
std::ostream &operator<<(std::ostream &out, Two_list_pointer_event_queue_key<I> k)
{
  k.write(out);
  return out;
}

// The interface for an item stored in the ::Pointer_event_queue
template <class Priority>
class Two_list_event_queue_item:
  public CGAL::In_place_list_base<Two_list_event_queue_item<Priority> >,
  public Ref_counted<Two_list_event_queue_item<Priority> >
{

  typedef Two_list_event_queue_item<Priority> This;

  Two_list_event_queue_item(const Two_list_event_queue_item &) {}
  void operator=(const Two_list_event_queue_item &) {}
public:
  typedef Two_list_pointer_event_queue_key<This> Key;
  Two_list_event_queue_item() { /*++two_list_remaining;*/}
  Two_list_event_queue_item(const Priority &t): time_(t){/*++two_list_remaining;*/}
  virtual ~Two_list_event_queue_item(){/*--two_list_remaining;*/
  }

  enum List {FRONT, BACK, INF};

  const Priority& time() const {return time_;}

  List in_list() const
  {
    return front_list_;
  }
  void set_in_list(List lt) {
    front_list_=lt;
  }

  bool operator<(const This &o) const {
    CGAL::Comparison_result c= CGAL::compare(time(), o.time());
    if (c != CGAL::EQUAL) return c== CGAL::SMALLER;
    else {
      if (kds() < o.kds()) return true;
      else if (kds() > o.kds()) return false;
      else {
	CGAL::Comparison_result c= compare_concurrent(Key((This*) this),Key((This*) &o));
	return c==CGAL::SMALLER;
      }
    }
  }

  virtual std::ostream& write(std::ostream &out) const{
    out << "Dummy event." << std::endl;
    return out;
  }
  virtual void process() {
    CGAL_error();
    CGAL_ERROR("Writing dummy queue element.");
  }
  virtual CGAL::Comparison_result compare_concurrent(Key , Key ) const {
    CGAL_error();
    return CGAL::EQUAL;
  };
  virtual bool merge_concurrent(Key , Key ) {
    CGAL_error();
    return false;
  }
  virtual void *kds() const{return NULL;}
  virtual void audit(Key) const{};
private:
  Priority time_;
  List front_list_;
};

template <class Priority>
inline std::ostream& operator<<(std::ostream &out, const Two_list_event_queue_item<Priority> &i)
{
  i.write(out);
  return out;
}

// The how a dummy item is stored in the ::Two_list_event_queue
/*
  One dummy item is used to represent events which will never happen.
*/
/*template <class Priority>
class Two_list_event_queue_dummy_item: public Two_list_event_queue_item<Priority>
{
  typedef Two_list_event_queue_item<Priority> P;
public:
  Two_list_event_queue_dummy_item(){}
  Two_list_event_queue_dummy_item(const Two_list_event_queue_dummy_item &):
    Two_list_event_queue_item<Priority>(){}
  virtual void process() {
    std::cerr << "Trying to process a NULL event.\n";
    CGAL_error();
  }
  virtual CGAL::Comparison_result compare_concurrent(typename P::Key , typename P::Key ) const{return CGAL::EQUAL;}
  virtual bool merge_concurrent(typename P::Key, typename P::Key){
    return false;
  }
  virtual std::ostream& write(std::ostream &out) const
  {
    out << "Never";
    return out;
  }
  virtual ~Two_list_event_queue_dummy_item(){}
  };*/

// The how a real item is stored in the ::Two_list_event_queue
/*
  This just stores an object of type Event and passes the virtual calls on to it.

  The object is reference counted so you don't have to worry about the
  queue deleting it or not.
*/
template <class Priority, class Event>
class Two_list_event_queue_item_rep: public internal::Two_list_event_queue_item<Priority>
{
  typedef  Two_list_event_queue_item<Priority> P;
public:
  Two_list_event_queue_item_rep(): internal::Two_list_event_queue_item<Priority>(){}
  Two_list_event_queue_item_rep(const Priority &t, const Event &e): 
    Two_list_event_queue_item<Priority>(t),
    event_(e){}

  virtual std::ostream& write(std::ostream &out) const
  {
    event_.write(out);
    out << " at " << P::time();
    return out;
  }
  virtual void process() {
    event_.process();
  }
  virtual void audit(typename P::Key k) const {
    event_.audit(k);
  }
  virtual CGAL::Comparison_result compare_concurrent(typename P::Key a, typename P::Key b) const{
    return event_.compare_concurrent(a,b);
  }
  virtual bool merge_concurrent(typename P::Key a, typename P::Key b){
    return event_.merge_concurrent(a,b);;
  }
  virtual void *kds() const{return event_.kds();}

  // Access the actual event
  const Event &event() const
  {
    return event_;
  }

 // Access the actual event
  Event &event()
  {
    return event_;
  }
  
  virtual ~Two_list_event_queue_item_rep(){}

protected:
  Event event_;
};


/* This is needed since the list cannot allocate an element of the abstract base class. I could just make it non-abstract. Why not?*/
/*template <class T>
struct Two_list_event_queue_item_allocator
{
  typedef Two_list_event_queue_dummy_item<T> dummy_value_type;

  typedef Two_list_event_queue_item<T>* pointer;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef const pointer const_pointer;
  typedef Two_list_event_queue_item<T>& reference;
  typedef const Two_list_event_queue_item<T>& const_reference;
  typedef dummy_value_type value_type;

  Two_list_event_queue_item_allocator(){}
  pointer allocate(size_t sz, void* =0) const
  {
    return static_cast<pointer>(::operator new(sz*sizeof(dummy_value_type)));
  }
  void deallocate(pointer p, size_type ) {
    ::operator delete(static_cast<void*>(p));
  }
  size_type max_size() const throw() {
    return (std::numeric_limits<size_type>::max)()/sizeof(dummy_value_type);
  }
  template <class TT>
  void construct(pointer p, const TT &) {
    new(static_cast<void*>(p)) dummy_value_type();
  }
  void destroy(pointer p) {
    p->~Two_list_event_queue_item<T>();
  }
  };*/

} } } //namespace CGAL::Kinetic::internal

namespace CGAL { namespace Kinetic {




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
template <class FK, bool INF=false, unsigned int TARGET=8>
class Two_list_pointer_event_queue
{
  typedef typename FK::Root PriorityT;
  typedef typename FK::FT NT;
  typedef Two_list_pointer_event_queue<FK, INF, TARGET> This;
  typedef typename internal::Two_list_event_queue_item<PriorityT> Item;

  typedef typename CGAL::In_place_list<Item, false/*,
						    internal::Two_list_event_queue_item_allocator<PriorityT>*/ > Queue;
  typedef typename Queue::iterator Iterator;

public:
  typedef PriorityT Priority;

  typedef internal::Two_list_pointer_event_queue_key<Item> Key;

  //! Construct it with a suggested size of sz.
  Two_list_pointer_event_queue(Priority start_time,
                               Priority end_time,
                               FK, int =0):
    null_event_(new internal::Two_list_event_queue_item<Priority>()){
    CGAL_precondition(!INF);
    initialize(start_time, end_time);
  }

  //! Construct it with a suggested size of sz.
  Two_list_pointer_event_queue(Priority start_time,
                               FK, int =0):
    null_event_(new internal::Two_list_event_queue_item<Priority>()){
    CGAL_precondition(INF);
    initialize(start_time);
  }

  ~Two_list_pointer_event_queue() {
    // release pointers
    std::vector<Key> all;
    all.reserve(front_.size()+back_.size());
    for (typename Queue::iterator it = front_.begin();
         it != front_.end(); ++it) {
      all.push_back(Key(&*it));
    }
    for (typename Queue::iterator it = back_.begin();
         it != back_.end(); ++it) {
      all.push_back(Key(&*it));
    }
    front_.clear();
    back_.clear();
    for (unsigned int i=0; i< all.size(); ++i) {
      unmake_event(all[i].get());
    }
  }

  /*template <class E>
  Key insert_at_end(const E & e) {

  }*/

  bool is_after_end(const Priority &t) const {
    if (INF) return false; 
    else return  CGAL::compare(t,end_priority()) == CGAL::LARGER;
  }

  //! insert value_type into the queue and return a reference to it
  /*!
    The reference can be used to update or erase it.
  */
  template <class E>
  Key insert(const Priority &t, const E & e) {
    CGAL_expensive_precondition(audit());
 

    Item *ni = make_event(t, e);

    //CGAL_exactness_assertion(t >= lbs_);
    //lb_=(std::min)(t, lb_);
    if ( is_after_end(ni->time())){
      return end_key();
    } 


    if (leq_ub(ni->time())) {
      ni->set_in_list(Item::FRONT);
      typename Queue::iterator iit=std::upper_bound(front_.begin(), front_.end(), *ni);
   
      front_.insert(iit, *ni);
      CGAL_expensive_assertion(audit());
      if (front_.size() > 2*max_front_size()) {
	shrink_front();
      }
      //++queue_front_insertions__;
    }  else if (front_.empty()) {
      CGAL_assertion(back_.empty());
      CGAL_assertion(INF || CGAL::compare(t, end_priority()) != CGAL::LARGER);
      //CGAL_assertion(INF || CGAL::compare(end_priority(), std::numeric_limits<Priority>::infinity()) == CGAL::SMALLER);
      if (true){
	//++queue_front_insertions__;
	front_.push_back(*ni);
	ub_= CGAL::to_interval(t).second;
	ni->set_in_list(Item::FRONT);
      } 
    } else {
      ni->set_in_list(Item::BACK);
      back_.push_back(*ni);
    }
    CGAL_expensive_postcondition(audit());
    CGAL_expensive_postcondition(contains(Key(ni)));
    //std::cout << "Made event " << ni << std::endl;
    return Key(ni);
  }


  //! remove the event referenced by item from the queue
  /*!
    \todo Add check that item is actually in the list
  */
  void erase(const Key &item) {
    //std::cout << "Erase event " << item.pointer() << std::endl;
    if (item== end_key()) return;
#ifndef NDEBUG
    if (!contains(item)) {
      std::cerr << "Erasing event not in queue ";
      item->write(std::cerr);
      std::cerr << std::endl;
    }
#endif
    CGAL_expensive_precondition(contains(item));
    CGAL_expensive_precondition(audit());
    Item *i=item.get();
    if (item->in_list() == Item::FRONT) {
      front_.erase(i);
      if (front_.empty() && !back_.empty()) grow_front();
    }
    else if (item->in_list() == Item::BACK) {
      back_.erase(i);
    }
    if (item->in_list() == Item::FRONT || item->in_list() == Item::BACK) {
      unmake_event(i);
    }

    CGAL_expensive_postcondition(audit());
  }

  template <class E>
  const E& get(const Key &item) const
  {
    return reinterpret_cast<internal::Two_list_event_queue_item_rep<Priority, E>*>( item.get())->event();
  }

  template <class E>
  E& get(const Key &item)  {
    return reinterpret_cast<internal::Two_list_event_queue_item_rep<Priority, E>*>( item.get())->event();
  }

  //! Replace the event referenced by item with a new event referenced by ne
  /*!  They must have exactly the same times associated with
    them. This is checked when expensive checks are turned on.
  */
  template <class NE>
  Key set(const Key &item, const NE &ne) {
    CGAL_expensive_precondition(contains(item));
    CGAL_precondition(item != end_key());
    Item *oi= item.get();
    typename Item::List front= item->in_list();
    Item *ni= make_event(item->time(), ne);
    ni->set_in_list(front);
    if (front != Item::INF) {
      if (front== Item::FRONT) {
	front_.insert(oi, *ni);
	front_.erase(oi);
      }
      else {
	back_.insert(oi, *ni);
	back_.erase(oi);
      }
      unmake_event(oi);
      return Key(ni);
    }
    else {
#ifndef NDEBUG
      bool found=false;
      for (unsigned int i=0; i< inf_.size(); ++i) {
	if (inf_[i].get() == oi) {
	  inf_[i]=ni;
	  found=true;
	  break;
	}
      }
      CGAL_postcondition(found);
      CGAL_USE(found);
#endif
      unmake_event(oi);
      return Key(ni);
    }
    // unreachable

  }

  //! Get the time of the next event to be processed.
  /*!  It is OK to call this if the queue is empty. In that case it
    will just return infinity.
  */
  const Priority& front_priority() const
  {
    CGAL_precondition(!front_.empty());
    return front_.front().time();
  }


  Key front() const
  {
    CGAL_precondition(!front_.empty());
    return Key(const_cast<Item*>(&front_.front()));
  }

  //! Access the time of a particular event
  const Priority& priority(const Key &item) const
  {
    return item->time();
  }

  //! empty
  bool empty() const
  {
    CGAL_precondition(!front_.empty() || back_.empty());
    return front_.empty() 
      || CGAL::compare(front_.front().time(), end_priority()) == CGAL::LARGER;
  }

  //! Remove the next event from the queue and process it.
  /*!
    Processing means that the process() method of the event object is called.
  */
  void process_front() {
    CGAL_precondition(!empty());
    CGAL_expensive_precondition(audit());
    if (!front_.empty()) {
      Item *i= &front_.front();
      CGAL_LOG(Log::SOME, "Processing event " << *i << std::endl);
      front_.pop_front();
      CGAL_expensive_postcondition(audit());
      if (front_.empty() && !back_.empty()) grow_front();
      i->process();

      /*if (!front_.empty() && i->time() == front_.front().time()) {
	CGAL_LOG(Log::SOME, "Degeneracy at time "
	<< i->time() << " the events are: "
	<< *i << " and " << front_.front() << std::endl);
	}*/

      unmake_event(i);
    }
    else {
      CGAL_assertion(back_.empty());
    }
  }

  //! debugging
  bool print() const
  {
    write(std::cout);
    return true;
  }

  void write(const Queue &q, std::ostream& out) const {
    for (typename Queue::const_iterator it = q.begin(); it != q.end(); ++it) {
      out << "(" << &*it << ": " << *it << ")";
      out << std::endl;
    }
  }

  bool write(std::ostream &out) const
  {
    write(front_, std::cout);
    out << "--" << ub_ << "--" << std::endl;
    write(back_, std::cout);
    out << std::endl;
    return true;
  }

  Key end_key() const
  {
    return null_event_;
  }

  

  const Priority& end_priority() const
  {
    return end_time_;
  }


  void set_interval(const Priority &start_time, const Priority &end_time) {
    CGAL_precondition(!INF);
    initialize(start_time, end_time);
  }
  
  void audit_events() const {
    for (typename Queue::const_iterator it= front_.begin(); it != front_.end(); ++it) {
      it->audit(Key(const_cast<Item*>(&*it)));
    }
    for (typename Queue::const_iterator it= back_.begin(); it != back_.end(); ++it) {
      it->audit(Key(const_cast<Item*>(&*it)));
    }
  }

  void audit_event(Key k) const {
    k->audit(k);
  }

  void clear() {
    front_.clear();
    back_.clear();
    //all_in_front_=false;
  }

protected:
  void initialize(const Priority &start_time, const Priority &end_time) {
    ub_=CGAL::to_interval(start_time).second;
    // should be nextafter
    step_=1;
    //all_in_front_= false;
    end_time_=end_time;
  }

 void initialize(const Priority &start_time) {
   CGAL_precondition(INF);
   ub_=CGAL::to_interval(start_time).second;
    // should be nextafter
   step_=1;
    //all_in_front_= false;
  }
  bool leq_ub(const Priority &t) const {
    //if (all_in_front_) return true;
    //else 
    // pretty much anything fast will have a fast to_interval, so use it
    std::pair<double,double> iv= CGAL::to_interval(t);
    if (iv.first > ub_) {
      return false;
    } else if (iv.second <= ub_) {
      return true;
    } else {
      return CGAL::compare(t, Priority(ub_)) != CGAL::LARGER;
    }
  }

 
  template <class E>
  Item *make_event(const Priority &t, E &e) {
    typedef typename boost::remove_const<E>::type NCE;
    Item *ni
      = new internal::Two_list_event_queue_item_rep<Priority, NCE>(t, e);
    intrusive_ptr_add_ref(ni);
    return ni;
  }

  void unmake_event(Item *i) {
    intrusive_ptr_release(i);
  }

  bool audit() {
    for (typename Queue::const_iterator it = front_.begin(); it != front_.end(); ++it) {
      Priority t= it->time();
      CGAL_assertion(leq_ub(t));
      CGAL_assertion(it->in_list()== Item::FRONT);
      //CGAL_exactness_assertion(t >= lb_);
    }
    for (typename Queue::const_iterator it = back_.begin(); it != back_.end(); ++it) {
      Priority t= it->time();
      CGAL_assertion(!leq_ub(t));
      CGAL_assertion(it->in_list()== Item::BACK);
    }
#ifndef NDEBUG
    for (unsigned int i=0; i< inf_.size(); ++i) {
      Priority t= inf_[i]->time();
      CGAL_assertion(INF || CGAL::compare(t, end_priority())== CGAL::LARGER);
      CGAL_assertion(inf_[i]->in_list() == Item::INF);
    }
#endif
    {
      typename Queue::const_iterator it = front_.begin();
      ++it;
      for (; it != front_.end(); ++it) {
	Priority tc= it->time();
	Priority tp= boost::prior(it)->time();
#ifndef NDEBUG
	if (CGAL::compare(tc, tp) == CGAL::SMALLER) {
	  std::cout << "ERROR: Out of order " << tc << std::endl << tp << std::endl << std::endl;
	  ++internal::get_static_audit_failures();
	}
#endif
	//CGAL_assertion(tc >= tp);
      }
    }
    return true;
  }
public:
  bool contains(const Key k) const
  {
    //if (k.pointer()->time() == std::numeric_limits<Priority>::infinity()) return true;
    for (typename Queue::const_iterator it = front_.begin(); it != front_.end(); ++it) {
      if (&*it == k.pointer()) return true;
    }
    for (typename Queue::const_iterator it = back_.begin(); it != back_.end(); ++it) {
      if (&*it == k.pointer()) return true;
    }
#ifndef NDEBUG
    for (unsigned int i=0; i< inf_.size(); ++i) {
      const Key j=inf_[i];
      const Key ki= k;
      if (j==ki) return true;
    }
#else
    if (k.pointer()->in_list() == Item::INF) return true;
#endif
    return false;
  }
protected:
  unsigned int select(Queue &source, Queue &target/*, bool binf*/) {
    CGAL_assertion_code(std::size_t sz= source.size() + target.size();)
    int count=0;
    Iterator it= source.begin();
    while (it != source.end()) {
      // CGAL_assertion(it->time() >= a);
    
      if (leq_ub(it->time())) {
	Item *i= &*it;
	Iterator t= boost::next(it);
	source.erase(it);
	it=t;
	target.push_back(*i);
	++count;
      }
      else {
	++it;
      }
    }
    CGAL_assertion(sz==source.size() + target.size());
    return count;
  }

  /*NT step() const{
    return (std::max)(ub_-lb_, NT(1));
    }*/

  NT av(NT a, NT b) const
  {
    return .5*(a+b);
  }

  template <class It>
  void set_front(It b, It e, typename Item::List val) {
    for (; b!= e; ++b) {
      b->set_in_list(val);
    }
  }
  template <class C, class It>
  void make_inf(C &c, It b, It e) {
    for (It cit = b; cit != e; ++cit) {
      CGAL_assertion(INF || CGAL::compare(cit->time(), end_priority()) == CGAL::LARGER);
      //std::cout << "Dropping inf event " << &*cit << std::endl;
#ifndef NDEBUG
      inf_.push_back(&*cit);
#endif
      cit->set_in_list(Item::INF);
      It p= boost::prior(cit);
      c.erase(cit);
      unmake_event(&*cit);
      cit=p;
    }
    //c.erase(b, e);
  }

  void grow_front(Queue &cand, int recursive_count=0) {
    const bool dprint=false;
    CGAL_assertion(front_.empty());
    CGAL_assertion(!cand.empty());
    //CGAL_assertion(!all_in_front_);
    CGAL_assertion(step_ != 0);
    if (dprint) std::cout << "Growing front from " << ub_ << " with step " 
			  << step_ << "(" << recursive_count << ") ";

    //CGAL_assertion(ub_<end_split());
    if (cand.size() + front_.size() < max_front_size()) {
      if (dprint) std::cout << "Setting ub to end of " << end_time_ << std::endl;
      front_.splice(front_.end(), cand);
      return;
    }

    //CGAL_assertion(!too_big(ub_));

    /*if ( CGAL::compare(end_priority(), Priority(ub_)) == CGAL::SMALLER) {
      all_in_front_=true;
      //ub_=end_split();
      }*/

    CGAL_assertion_code(unsigned int num=)
      select(cand, front_/*, all_in_front_*/);
    CGAL_assertion(front_.size() >= num);
    /*if (all_in_front_) {
      make_inf(cand, cand.begin(), cand.end());
      } else*/
    if (front_.empty()) {
      if (recursive_count > 10) {
	// do something
	std::cerr << "Too many recursions " << std::endl;
	//all_in_front_=true;
	ub_=CGAL::to_interval(end_time_).second; //CGAL::to_interval(end_time_).second*2;// std::numeric_limits<double>::infinity();
	/*unsigned int num=*/ select(cand, front_/*, all_in_front_*/);
	select(back_, front_);
	make_inf(cand, cand.begin(), cand.end());
	make_inf(back_, back_.begin(), back_.end());
	//grow_front(cand, recursive_count+1);
      } else {
	if (dprint) {
	  std::cout << "undershot." << std::endl;
	  write(front_, std::cout);
	  std::cout << "-- " << ub_ << "--\n";
	  write(cand, std::cout);
	  std::cout << "--\n";
	  write(back_, std::cout);
	}
	ub_+= step_;
       	step_*=2.0;
	CGAL_assertion(step_!=0);
	cand.splice(cand.end(), back_);
	grow_front(cand, recursive_count+1);
      }
    } else {
      //      unsigned int ncand= cand.size();
      back_.splice(back_.begin(), cand);
      if (front_.size() >  max_front_size()) {
	if (recursive_count > 10) {
	  //std::cerr << "Gave up on isolating front. Let with size " << front_.size() << " ub=" << ub_ << "step=" << step_ <<  "\n";
	  return;
	} else {
	  // keep the bit length shortish
	  double frac= .75+.25*max_front_size()/static_cast<double>(front_.size()+1);
	  double ostep= step_;
	  CGAL_assertion(frac < 1.1);
	  CGAL_assertion(frac >= 0.0);
	  //double rfrac= std::ceil(frac*256.0)/256.0;
	  step_*=frac;
	  //else nstep = step_*.6;
	  //CGAL_assertion(nstep >0);
	  cand.swap(front_);
	  //ub_=lb_;
	  if (step_ == 0) {
	    CGAL_ERROR( "underflow in queue ");
	    CGAL_ERROR_WRITE(write(LOG_STREAM));
	    CGAL_assertion(cand.empty());
	    step_=.0000001;
	    return;
	  } else {
	    //CGAL_assertion(!all_in_front_);
	    
	    if (dprint) {
	      std::cout << "...overshot" << std::endl;
	      write(front_, std::cout);
	      std::cout << "-- " << ub_ << "(" <<step_ << ")" << "--\n";
	      write(cand, std::cout);
	      std::cout << "--\n";
	      write(back_, std::cout);
	    }
	    ub_=ub_-ostep+step_;  
	    grow_front(cand, recursive_count+1);
	  }
	}
      }
      else {
	if (dprint) std::cout << std::endl;
      }
    }
    CGAL_postcondition(cand.empty());
  }

  void grow_front() {
    //++growth__;
    //std::cout << "Growing front from " << ub_ << std::endl;
    //CGAL_assertion(is_valid());
    CGAL_precondition(!back_.empty());
    CGAL_precondition(front_.empty());
    CGAL_assertion_code(std::size_t sz= front_.size()+back_.size()+ inf_.size());
    Queue cand;
    cand.splice(cand.begin(), back_);
    ub_ += step_;
    grow_front(cand);
    set_front(front_.begin(), front_.end(), Item::FRONT);
    front_.sort();
    ub_= CGAL::to_interval(front_.back().time()).second;
    // hmmmm, now I should make a second pass to merge. Ick.

    CGAL_assertion(sz==front_.size()+back_.size() + inf_.size());
    CGAL_assertion(audit());
    //std::cout << "to " << ub_ << std::endl;
  }

  void shrink_front() {
    //++shrink__;
    //std::cout << "Shrinking front from " << ub_ << std::endl;
    typename Queue::iterator it=front_.begin();
    unsigned int mf= max_front_size();
    for (unsigned int i=0; i < mf; ++i) {
      ++it;
    }

    // use tii_ so that it does not subdivide

    double split =  CGAL::to_interval(it->time()).second;
    if (!INF && (CGAL::compare(end_priority(), it->time())==CGAL::SMALLER
		 || CGAL::compare(end_priority(), Priority(split))==CGAL::SMALLER)) {
      std::cout << "Past end in Queue " << end_priority() << ", "
		<< it->time() << ", " << Priority(split) << std::endl;
      //all_in_front_=true;
      ub_= CGAL::to_interval(end_time_).second;
      set_front(back_.begin(), back_.end(), Item::FRONT);
      front_.splice(front_.end(), back_);

      while (it != front_.end()) {
	if (CGAL::compare(it->time(), end_priority())==CGAL::LARGER) break;
      }
      set_front(it, front_.end(), Item::INF);
      std::vector<Item*> temp;
      for (typename Queue::iterator c= it; c != front_.end(); ++c) {
	temp.push_back(&*c);
#ifndef NDEBUG
	inf_.push_back(&*c);
#endif
	//std::cout << "Dropping inf event " << &*c << std::endl;
      }
      front_.erase(it, front_.end());

   
    
      for (unsigned int i=0; i< temp.size(); ++i) {
	unmake_event(temp[i]);
      }
    } else {
      while (CGAL::compare(it->time(), Priority(split)) != CGAL::LARGER
	     && it != front_.end()) ++it;
      
      if (it != front_.end()) {
	
	set_front(it, front_.end(), Item::BACK);
	back_.splice(back_.begin(), front_, it, front_.end());
	double oub=ub_;
	//all_in_front_=false;
	ub_ = split;
	//CGAL_assertion(!too_big(ub_));
	//CGAL_assertion(ub_ <= end_split());
	step_= std::abs(ub_-oub);
	if (step_<=0) {
	  /*if (dprint) std::cout << "fixing step since " << oub 
	    << " equals new bound" << std::endl;*/
	  CGAL_ERROR("Roundoff error in split " << split << " " << ub_ << " "
			     <<  oub);
	  step_=.00001;
	}
	//CGAL_assertion(!all_in_front_);
	/*CGAL_postcondition_code(if (step_<0) std::cerr << step_ << std::endl;);
	CGAL_postcondition_code(if (step_<0) std::cerr << ub_ << std::endl;);
	CGAL_postcondition_code(if (step_<0) std::cerr << oub << std::endl;);
	CGAL_postcondition_code(if (step_<0) std::cerr << front_.back().time() << std::endl;);
	CGAL_postcondition_code(if (step_==0) for (typename Queue::const_iterator it=front_.begin(); it != front_.end(); ++it) std::cout << *it << std::endl);
	CGAL_postcondition(step_>=0);*/
      }
    }
  }

  /*NT end_split() const
    {
    return end_split_;
    }*/

  /*const Priority& end_time() const
    {
    return end_time_;
    }*/

  unsigned int max_front_size() const
  {
    return TARGET;
    //return (std::max)(4U, static_cast<unsigned int>(std::sqrt(static_cast<double>(front_.size()+back_.size()))));
  }

  //typename FK::To_isolating_interval tii_;
  Queue front_, back_;
#ifndef NDEBUG
  std::vector<Key> inf_;
#endif
  double ub_;
  double step_;
  //bool all_in_front_;
  Priority end_time_;
  Key null_event_;
  //NT end_split_;
};

template <class D, unsigned int T, bool INF>
std::ostream &operator<<(std::ostream &out, const Two_list_pointer_event_queue<D, INF, T> &q)
{
  q.write(out);
  return out;
}


} } //namespace CGAL::Kinetic
#endif
