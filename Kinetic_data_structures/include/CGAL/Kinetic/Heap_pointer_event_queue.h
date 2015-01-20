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

#ifndef CGAL_KINETIC_HEAP_POINTER_EVENT_QUEUE_H
#define CGAL_KINETIC_HEAP_POINTER_EVENT_QUEUE_H
#include <CGAL/Kinetic/basic.h>
#include <iostream>
#include <vector>
#include <utility>
#include <functional>
#include <CGAL/assertions.h>
#include <iostream>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/internal/infinity_or_max.h>
#include <algorithm>

namespace CGAL { namespace Kinetic { namespace internal {

template <class Priority>
class Heap_pointer_event_queue_item_handle;

// The interface for an item stored in the ::Heap_pointer_event_queue
template <class Priority>
class Heap_pointer_event_queue_item: public Ref_counted<Heap_pointer_event_queue_item<Priority> >
{
  typedef Ref_counted<Heap_pointer_event_queue_item<Priority> > P;
public:

  /* struct Key: public typename P::Handle {
    Key(){}
    Key(Item_handle h): P::Handle(h){}
    bool is_valid() const {
      return this->get() != NULL;
    }
    };*/
  typedef Heap_pointer_event_queue_item_handle<Priority>  Key;

  Heap_pointer_event_queue_item():bin_(-1), time_(infinity_or_max<Priority>(Priority(0))){}
  Heap_pointer_event_queue_item(int bin, const Priority &t): bin_(bin), time_(t){}

  virtual std::ostream& write(std::ostream &out) const =0;
  const Priority& time() const {return time_;};
  virtual void process() =0;
  virtual void audit(Key) const=0;
  virtual void *kds() const =0;
  virtual CGAL::Comparison_result compare_concurrent(Key a, Key b) const =0;
  void set_bin(int bin) const { bin_=bin;}
  int bin() const {return bin_;};
  virtual ~Heap_pointer_event_queue_item(){}
  //virtual const void *event() const=0;
private:
  mutable int bin_;
  Priority time_;
};

CGAL_OUTPUT1(Heap_pointer_event_queue_item)


template <class Priority>
class Heap_pointer_event_queue_item_handle: public Ref_counted<Heap_pointer_event_queue_item<Priority> >::Handle
{
  typedef typename Ref_counted<Heap_pointer_event_queue_item<Priority> >::Handle P;
public:
  std::ostream &write(std::ostream &out) const {
    return P::operator*().write(out);
  }
  template <class T>
  Heap_pointer_event_queue_item_handle(T* t): P(t){}
  Heap_pointer_event_queue_item_handle(){}
  Heap_pointer_event_queue_item_handle(const P&p): P(p){}
};

CGAL_OUTPUT1(Heap_pointer_event_queue_item_handle)

// The how a dummy item is stored in the ::Heap_pointer_event_queue
/*
  One dummy item is used to represent events which will never happen.
*/
template <class Priority>
class Heap_pointer_event_queue_dummy_item: public Heap_pointer_event_queue_item<Priority>
{
  typedef Heap_pointer_event_queue_item<Priority> P;
public:
  Heap_pointer_event_queue_dummy_item(): P(-2, internal::infinity_or_max(Priority())){}
  virtual void process() {
  }
  virtual void *kds() const {return NULL;}
  virtual CGAL::Comparison_result compare_concurrent(typename P::Key a, 
						     typename P::Key b) const {
    if (a < b) return CGAL::SMALLER;
    else if (b < a) return CGAL::LARGER;
    else return CGAL::EQUAL;
    //return CGAL::compare(a,b);
  }
  virtual std::ostream& write(std::ostream &out) const
  {
    out << "Never.";
    return out;
  }
  virtual void audit(typename P::Key) const {
    std::cout << "Auditing a dummy event" << std::endl;
  }
  virtual ~Heap_pointer_event_queue_dummy_item(){}
};

// The how a real item is stored in the ::Heap_pointer_event_queue
/*
  This just stores an object of type Event and passes the virtual calls on to it.

  The object is reference counted so you don't have to worry about the
  queue deleting it or not.
*/
template <class Priority, class Event>
class Heap_pointer_event_queue_item_rep: public internal::Heap_pointer_event_queue_item<Priority>
{
  typedef internal::Heap_pointer_event_queue_item<Priority> P;
public:
  //typedef CGAL::Ref_counted_pointer<Heap_pointer_event_queue_item_rep<Priority, Event> > Pointer;
  //typedef CGAL::Ref_counted_pointer<const Heap_pointer_event_queue_item_rep<Priority, Event> > Const_pointer;
  Heap_pointer_event_queue_item_rep(): internal::Heap_pointer_event_queue_item<Priority>(){}
  Heap_pointer_event_queue_item_rep(const Priority &t, const Event &e,
				    unsigned int bin): internal::Heap_pointer_event_queue_item<Priority>(bin, t),
						       event_(e){}
  virtual void process() {
    event_.process();
  }
  virtual void *kds() const {return event_.kds();}
  virtual CGAL::Comparison_result compare_concurrent(typename P::Key a, 
						     typename P::Key b) const {
    return event_.compare_concurrent(a,b);
  }
  virtual void audit(typename P::Key k) const {
    event_.audit(k);
  }
  virtual std::ostream& write(std::ostream &out) const
  {
    out << "(";
    event_.write(out);
    out << ")";
    return out;
  }
  // Access the actual event
  /*
    There is no non-const access so you can't change the time.
  */
  const Event &event() const
  {
    return event_;
  }
  Event &event()
  {
    return event_;
  }
  void set_event(const Event &e) {
    event_=e;
  }
  virtual ~Heap_pointer_event_queue_item_rep(){}

protected:
  Event event_;
};

} } } //namespace CGAL::Kinetic::internal

namespace CGAL { namespace Kinetic {

template <class Priority> class Bin_pointer_event_queue;

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
template <class FK, bool INF=false>
class Heap_pointer_event_queue
{
public:
  typedef typename FK::Root Priority;
  friend class Bin_pointer_event_queue<Priority>;
  typedef Heap_pointer_event_queue<Priority> This;
  typedef internal::Heap_pointer_event_queue_item<Priority> Item;
  typedef typename Item::Handle Item_handle;
  typedef enum Child {FIRST=0, SECOND=1}
    Child;

  class Compare
  {
  public:
    Compare(){}
    bool operator()(Item_handle i0, Item_handle i1) const
    {
      CGAL::Comparison_result cr= CGAL::compare(i0->time(), i1->time());
      if (cr == CGAL::SMALLER) return true;
      else if (cr == CGAL::LARGER) return false; 
      else {
	if (i0->kds() < i1->kds()) return true;
	else if (i0->kds() > i1->kds()) return false;
	else {
	  cr= i0->compare_concurrent(Key(i0), Key(i1));
	  if (cr == CGAL::SMALLER) return true;
	  else if (cr == CGAL::LARGER) return false;
	  else {
	    //CGAL_error();
	    return false;
	  }
	}
	
      }
    }
  };


  //typedef Priority_t Priority;

  //! The key to identify something in the queue.
  /*!
    It uses a prettified version of the ref counted pointer.
  */
  typedef typename Item::Key Key;
  /*struct Key: public Item_handle {
    Key(){};
    Key(Item*i): Item_handle(i){}
    Key(Item_handle ih): Item_handle(ih){}
    static Key null() {
    return Key();
    }
    };*/
  //typedef Item_handle Key;

  //! Construct it with a suggested size of sz.
  Heap_pointer_event_queue(const Priority&, const Priority &end, FK, int sz=1000): end_(end) {
    queue_.reserve(sz);
    null_event_= Key(new internal::Heap_pointer_event_queue_dummy_item<Priority>());
  }

  Heap_pointer_event_queue(const Priority&, FK, int sz=1000) {
    queue_.reserve(sz);
    null_event_= Key(new internal::Heap_pointer_event_queue_dummy_item<Priority>());
  }

  const Priority& end_priority() const
  {
    return end_;
  }

  void set_interval(const Priority &, const Priority &e)
  {
    end_=e;
  }

  bool is_after_end(const Priority &t) const {
    if (INF) return false;
    else return  CGAL::compare(t,end_priority()) == CGAL::LARGER;
  }
  void audit_events() const {
    for (typename  std::vector<Item_handle>::const_iterator it= queue_.begin(); 
	 it != queue_.end(); ++it) {
      (* it)->audit(Key(*it));
    }
  }

  void audit_event(Key k) const {
    k->audit(k);
  }

  //! insert value_type into the queue and return a reference to it
  /*!
    The reference can be used to update or erase it.
  */
  template <class E>
  Key insert(const Priority &t, const E & e) {
    if (!is_after_end(t)) {
      Item_handle k= new internal::Heap_pointer_event_queue_item_rep<Priority, E>(t, e, static_cast<unsigned int>(queue_.size()));
      queue_.push_back(k);
      bubble_up(static_cast<unsigned int>(queue_.size()-1));
      //std::push_heap(queue_.begin(), queue_.end());
      CGAL_expensive_postcondition(is_valid());
      CGAL_postcondition(k->bin() != -1);
      return Key(k);
    } else {
      return null_event_;
    }
  }

  //! remove the event referenced by item from the queue
  void erase(const Key &item) {
    if (item == end_key()) return;
    CGAL_expensive_precondition(item != Key());
    CGAL_expensive_precondition(is_in_heap(item));
    int bin= item->bin();
    //if (bin ==-1) return;
    CGAL_precondition(static_cast<unsigned int> (bin) < queue_.size() && bin >= 0);

    // this is a bit more work than strictly necessary
    swap(static_cast<unsigned int>(queue_.size()-1), bin);
    queue_.pop_back();
    if (static_cast<unsigned int>(bin) < queue_.size()) {
      bubble(bin);
    }
    item->set_bin(-1);
    CGAL_expensive_postcondition(is_valid());
  }

  //! The pointer type to use to access an event of type E.

  /*
    template <class E>
    struct Event_pointer {
    typedef const Heap_pointer_event_queue_item_rep<Priority, E>* Pointer;
    struct Handle: public  Heap_pointer_event_queue_item_rep<Priority, E>::Const_pointer {
    Handle(Pointer p):
    Heap_pointer_event_queue_item_rep<Priority, E>
    ::Const_pointer(reinterpret_cast<const Heap_pointer_event_queue_item_rep<Priority, E>*>(p)){
    }
    Handle(){}

    };
    };*/

  //! Access a ref counted pointer corresponding to the event with item.
  /*!
    Note, there is a reinterpret cast going on in here, so the type better be correct.

    If you store the return value, please wrap it in an Event_pointer
    to make reference counting work.
  */
  /*template <class E>
    typename Event_pointer<E>::Pointer event(const Key &item, const E &) const {
    return reinterpret_cast<typename Event_pointer<E>::Pointer >( item.get());
    }*/

  template <class E>
  const E& get(const Key &item) const
  {
    CGAL_precondition(item && item != null_event_);
    return reinterpret_cast<internal::Heap_pointer_event_queue_item_rep<Priority, E>*>( item.get())->event();
  }


  template <class E>
  E& get(Key item)
  {
    CGAL_precondition(item && item != null_event_);
    typename internal::Heap_pointer_event_queue_item<Priority> *ptr= item.get();
    typename internal::Heap_pointer_event_queue_item_rep<Priority, E>* nptr
      = reinterpret_cast<internal::Heap_pointer_event_queue_item_rep<Priority, E>*>(ptr);
    return nptr->event();
  }

  //! Replace the event referenced by item with a new event referenced by ne
  /*!  They must have exactly the same times associated with
    them. This is checked when expensive checks are turned on.
  */
  template <class NE>
  Key set(const Key &item, const NE &ne) {
    CGAL_expensive_precondition(is_in_heap(item));
    CGAL_precondition(item && item != null_event_);
    //CGAL_expensive_precondition(item.get()->time() == ne.time());
    CGAL_expensive_precondition(item);
    int bin = item.get()->bin();
    CGAL_expensive_precondition(Key(queue_[bin])==item);
    queue_[bin]= new internal::Heap_pointer_event_queue_item_rep<Priority, NE>(item.get()->time(), ne, bin);
    CGAL_expensive_postcondition(is_valid());
    return queue_[bin];
  }

  //! Get the time of the next event to be processed.
  /*!  It is OK to call this if the queue is empty. In that case it
    will just return infinity.
  */
  Priority front_priority() const
  {
    CGAL_precondition(!queue_.empty());
    return queue_.front()->time();
  }

  //! Access the time of a particular event
  const Priority& priority(const Key &item) const
  {
    CGAL_precondition(bool(item));
    return item->time();
  }

  //! empty
  bool empty() const
  {
    return queue_.empty();
  }

  //! Clear everything
  void clear() {
    // ref counted pointers are nice.
    queue_.clear();
  }



  //! Remove the next event from the queue and process it.
  /*!
    Processing means that the process() method of the event object is called.
  */
  void process_front() {
    CGAL_precondition(!empty());
    //if (queue_.front()->time() < end_priority()) {
    //CGAL_precondition_code(Item_handle k= queue_.front());
    CGAL_LOG(Log::SOME, "Processing " << queue_.front() << std::endl);
    Item_handle ih= queue_.front();
    pop_front();
    //std::pop_heap(queue_.begin(), queue_.end());
    ih->process();
    //CGAL_expensive_postcondition(is_valid());
    /*}
      else {
      clear();
      }*/
  }

  //! debugging
  bool print() const
  {
    write(std::cout);
    return true;
  }
  /*
    void process_final_events() {
    bool do_i_use_this;
    std::vector<Item_handle> c;
    c.swap(final_events_);
    for (unsigned int i=0; i< c.size(); ++i){
    c[i]->set_processed(true);
    }
    }*/
  //! Write the events in order to out
  /*!  When warnings are turned on, it will print if two events occur
    at the same time.
  */
  bool write(std::ostream &out) const
  {
    //std::cout << "Not writing queue.\n";
    //return true;
#if 0
    std::vector<Item_handle> bins;

    for (unsigned int i=0; i< queue_.size(); ++i) {
      bins.push_back(queue_[i]);
    }

    std::sort(bins.begin(), bins.end(), compare_);

    typename std::vector<Item_handle>::const_iterator curi= bins.begin(), endi= bins.end();
#else
    typename std::vector<Item_handle>::const_iterator curi= queue_.begin(), endi= queue_.end();
#endif
    for (; curi != endi; ++curi) {
      Priority t=(*curi)->time();
      // HACK HACK because CGAL::to_double(t) won't compile with gcc 3.4
      std::pair<double,double> d= CGAL::to_interval(t);
      //CGAL::to_double(t);
      out << "<" << d.first << "..." << d.second << ": ";
      (*curi)->write(out);
      out << ">\n";
    }
    out << std::endl;
#if 0
    for (unsigned int i=1; i< bins.size(); ++i) {
      if (bins[i-1]->time()== bins[i]->time()) {
	out << "Warning, two concurrent events: " << bins[i] << " and "
	    << bins[i-1] << " at time " << bins[i]->time() << std::endl;
      }
    }
#endif

    /*if (!final_events.empty()) {
      out << "Finally: \n";
      for (unsigned int i=0; i< final_events_.size(); ++i){
	final_events_[i]->write(out) << std::endl;
      }
      }*/

    
    //out << std::endl;
    return false;
  }

  Key end_key() const
  {
    return null_event_;
  }

  bool contains(Key k) const {
    return is_in_heap(k);
  }


  Key front() const {
    return Key(queue_.front());
  }


protected:
  //! Stores the priorities and data and a refersence back to the _queue
  std::vector<Item_handle> queue_;
  Compare compare_;
  //std::vector<Item_handle> final_events_;
  Key null_event_;
  Priority end_;

  bool is_in_heap(Key k) const
  {
    for (unsigned int i=0; i< queue_.size(); ++i) {
      if (queue_[i]==k) return true;
    }
    return false;
  }

  //! Exchange the position of two bins
  void swap(unsigned int i, unsigned int j) {
    std::swap(queue_[i], queue_[j]);
    int bin= queue_[i]->bin();
    queue_[i]->set_bin(queue_[j]->bin());
    queue_[j]->set_bin(bin);
  }

  //! heap order parent
  int parent(unsigned int i) const
  {
    return (static_cast<int>(i)-1)/2;
  };

  //! heap order child
  unsigned int child(unsigned int i, Child which) const
  {
    return 2*i+1+which;
  }

  //! Compare the items of two indices
  bool less_items(unsigned int i, unsigned int j) const
  {
    if (static_cast<unsigned int>(j) >= queue_.size()) return true;
    else if (static_cast<unsigned int>(i) >= queue_.size()) return false;
    else {
      return compare_(queue_[i], queue_[j]);
    }
  }
  //! Check if item i is less than its c'th child
  bool less_than_child(unsigned int  i, unsigned int c) const
  {
    int ch= child(i,c);
    return less_items(i, ch);
  }
  //! check of item i is less than its parent
  bool less_than_parent(unsigned int i) const
  {
    if (i==0) return false;
    else return less_items(i, parent(i));
  }

  //! remove the front of the queue and fix the heap property
  void pop_front() {
    CGAL_expensive_precondition(!empty());
    Item_handle item= queue_.front();
    swap(0, static_cast<unsigned int>(queue_.size()-1));
    queue_.back()->set_bin(-1);
    queue_.pop_back();
    if (!queue_.empty()) bubble_down(0);
  }

  //! Make sure that item i is not less than its parents and fix
  void bubble_down(unsigned int i) {
    CGAL_expensive_precondition(i<queue_.size());
    while (i < queue_.size()) {
      unsigned int lc= child(i,FIRST), rc= child(i,SECOND);
      if (!less_items(rc, lc)) {        // make the sorting stable
	if (less_items(lc, i)) {
	  swap(i, lc);
	  i=lc;
	} else break;
      }
      else {
	if (less_items(rc, i)) {
	  swap(i, rc);
	  i=rc;
	} else break;
      }
    }
  }
  //! make sure the heap order is respected above item i
  void bubble_up(unsigned int i) {
    CGAL_expensive_precondition(i< queue_.size());
    while (less_than_parent(i)) {
      swap(i, parent(i));
      i=parent(i);
    }
  }
  //! both buttles
  void bubble(unsigned int i) {
    bubble_up(i);
    bubble_down(i);
  }

  //! debugging
  bool is_valid() const
  {
#ifndef CGAL_NO_ASSERTIONS
    for (unsigned int i=0; i< queue_.size(); ++i) {
      int back= queue_[i]->bin();
      CGAL_assertion(static_cast<unsigned int>(back)==i || write(std::cerr));
      CGAL_assertion(!less_than_parent(i) || write(std::cerr));
      //CGAL_assertion(less_than_child(i, FIRST));
      //CGAL_assertion(less_than_child(i, SECOND));
    }
#endif
    return true;
  }

};

template <class D, bool INF>
std::ostream &operator<<(std::ostream &out, const Heap_pointer_event_queue<D, INF> &q)
{
  q.write(out);
  return out;
}


} } //namespace CGAL::Kinetic
#endif
