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

#ifndef CGAL_KDS_QUEUE_H
#define CGAL_KDS_QUEUE_H
#include <CGAL/KDS/basic.h>
#include <iostream>
#include <vector>
#include <utility>
#include <functional>
#include <CGAL/assertions.h>
#include <iostream>
#include <CGAL/KDS/Ref_counted.h>
#include <CGAL/KDS/internal/infinity_or_max.h>
#include <algorithm>



CGAL_KDS_BEGIN_INTERNAL_NAMESPACE


// The interface for an item stored in the ::Heap_pointer_event_queue
template <class Priority>
class Heap_pointer_event_queue_item: public Ref_counted<Heap_pointer_event_queue_item<Priority> > {
public:
  Heap_pointer_event_queue_item():bin_(-1), time_(infinity_or_max<Priority>(Priority(0))){}
  Heap_pointer_event_queue_item(unsigned int bin, const Priority &t): bin_(bin), time_(t){}

  virtual void write(std::ostream &out) const =0;
  const Priority& time() const {return time_;};
  virtual void process(const Priority &t) =0;
  void set_bin(int bin) const { bin_=bin;}
  int bin() const {return bin_;};
  virtual ~Heap_pointer_event_queue_item(){}
  //virtual const void *event() const=0;
private:
  mutable int bin_;
  Priority time_;
};

template <class Priority>
inline std::ostream& operator<<(std::ostream &out, const Heap_pointer_event_queue_item<Priority> &i) {
  i.write(out);
  return i;
}


// The how a dummy item is stored in the ::Heap_pointer_event_queue
/*
  One dummy item is used to represent events which will never happen.
*/
template <class Priority>
class Heap_pointer_event_queue_dummy_item: public Heap_pointer_event_queue_item<Priority> {
public:
  Heap_pointer_event_queue_dummy_item(): Heap_pointer_event_queue_item<Priority>(-2, internal::infinity_or_max(Priority())){}
  virtual void process(const Priority &){
  }
  virtual void write(std::ostream &out) const {
    out << "Never.";
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
class Heap_pointer_event_queue_item_rep: public internal::Heap_pointer_event_queue_item<Priority> {
public:
  //typedef CGAL::Ref_counted_pointer<Heap_pointer_event_queue_item_rep<Priority, Event> > Pointer;
  //typedef CGAL::Ref_counted_pointer<const Heap_pointer_event_queue_item_rep<Priority, Event> > Const_pointer;
  Heap_pointer_event_queue_item_rep(): internal::Heap_pointer_event_queue_item<Priority>(){}
  Heap_pointer_event_queue_item_rep(const Priority &t, const Event &e, 
			       unsigned int bin): internal::Heap_pointer_event_queue_item<Priority>(bin, t), 
						  event_(e){}
  virtual void process(const Priority &t){
    event_.process(t);
  }
  virtual void write(std::ostream &out) const {
    out << event_;
  }
  // Access the actual event
  /*
    There is no non-const access so you can't change the time. 
  */
  const Event &event() const {
    return event_;
  }
  void set_event(const Event &e) {
    event_=e;
  }
  virtual ~Heap_pointer_event_queue_item_rep(){}

protected:
  Event event_;
};

CGAL_KDS_END_INTERNAL_NAMESPACE


CGAL_KDS_BEGIN_NAMESPACE;


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
template <class Priority_t>
class Heap_pointer_event_queue {
  friend class Bin_pointer_event_queue<Priority_t>;
  typedef Heap_pointer_event_queue<Priority_t> This;
  typedef internal::Heap_pointer_event_queue_item<Priority_t> Item;
  typedef typename Item::Pointer Item_handle;
  typedef enum Child {FIRST=0, SECOND=1} Child;

  class Compare {
  public:
    Compare(){}
    bool operator()(Item_handle i0, Item_handle i1) const {
      return i0->time() < i1->time();
    }
  };
public:

  typedef Priority_t Priority;

  //! The key to identify something in the queue.
  /*!
    It uses a prettified version of the ref counted pointer.
  */
  typedef Item_handle Key;
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
  Heap_pointer_event_queue(const Priority&, int sz=1000) {
    queue_.reserve(sz);
    null_event_= Key(new internal::Heap_pointer_event_queue_dummy_item<Priority>());
  }

  //! insert value_type into the queue and return a reference to it
  /*!
    The reference can be used to update or erase it.
  */
  template <class E>
  Key insert(const Priority &t, const E & e){
    Item_handle k= new internal::Heap_pointer_event_queue_item_rep<Priority, E>(t, e, queue_.size());
    queue_.push_back(k);
    bubble_up(queue_.size()-1);
    //std::push_heap(queue_.begin(), queue_.end());
    CGAL_expensive_postcondition(is_valid());
    CGAL_postcondition(k->bin() != -1);
    return Key(k);
  }

  //! remove the event referenced by item from the queue
  void erase(const Key &item){
    if (item == end_key()) return;
    CGAL_expensive_precondition(item);
    CGAL_expensive_precondition(is_in_heap(item));
    int bin= item->bin();
    //if (bin ==-1) return;
    CGAL_precondition(static_cast<unsigned int> (bin) < queue_.size() && bin >= 0);

    // this is a bit more work than strictly necessary
    swap(queue_.size()-1, bin);
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
  const E& get(const Key &item) const {
    CGAL_precondition(item && item != null_event_);
    return reinterpret_cast<internal::Heap_pointer_event_queue_item_rep<Priority, E>*>( item.get())->event();
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
  Priority front_priority() const {
    CGAL_precondition(!queue_.empty());
    return queue_.front()->time();
  }

  //! Access the time of a particular event
  Priority priority(const Key &item) const {
    CGAL_precondition(item);
    return item->time();
  }

  //! empty
  bool empty() const {
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
    CGAL_precondition_code(Item_handle k= queue_.front());
    CGAL_KDS_LOG(LOG_LOTS, "Processing " << queue_.front() << std::endl);
    Item_handle ih= queue_.front();
    pop_front();
    //std::pop_heap(queue_.begin(), queue_.end());
    ih->process(ih->time());
    CGAL_expensive_postcondition(is_valid());
  }

  //! debugging
  bool print() const {
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
  bool write(std::ostream &out) const {
    //std::cout << "Not writing queue.\n";
    //return true;
#if 1
    std::vector<Item_handle> bins;
    
    for (unsigned int i=0; i< queue_.size(); ++i){
      bins.push_back(queue_[i]);
    }

    std::sort(bins.begin(), bins.end(), compare_);
    
    typename std::vector<Item_handle>::const_iterator curi= bins.begin(), endi= bins.end();
#else
    std::sort(queue_.begin(), queue_.end(), compare_);
    
    typename std::vector<Item_handle>::const_iterator curi= queue_.begin(), endi= queue_.end();
#endif
    for (; curi != endi; ++curi){
      Priority t=(*curi)->time();
      // HACK HACK because CGAL::to_double(t) won't compile with gcc 3.4
      double d= CGAL::to_double(t); //CGAL::to_double(t);
      out << "<" << d << ": ";
      (*curi)->write(out);
      out << ">\n";
    }
    out << std::endl;
#if 0
    for (unsigned int i=1; i< bins.size(); ++i){
      if (bins[i-1]->time()== bins[i]->time()){
	out << "Warning, two concurrent events: " << bins[i] << " and " 
	    << bins[i-1] << " at time " << bins[i]->time() << std::endl;
      }
    }
#endif
    //out << std::endl;
    return false;
  }

  Key end_key() const {
    return null_event_;
  }
  

protected:
  //! Stores the priorities and data and a refersence back to the _queue
  std::vector<Item_handle> queue_;
  Compare compare_;
  std::vector<Item_handle> final_events_;
  Key null_event_;

  bool is_in_heap(Key k) const {
    for (unsigned int i=0; i< queue_.size(); ++i){
      if (queue_[i]==k) return true;
    }
    return false;
  }


  //! Exchange the position of two bins
  void swap(unsigned int i, unsigned int j){
    std::swap(queue_[i], queue_[j]);
    int bin= queue_[i]->bin();
    queue_[i]->set_bin(queue_[j]->bin());
    queue_[j]->set_bin(bin);
  }

  //! heap order parent
  int parent(unsigned int i) const {
    return (static_cast<int>(i)-1)/2;
  };

  //! heap order child
  unsigned int child(unsigned int i, Child which) const {
    return 2*i+1+which;
  }
    
  //! Compare the items of two indices
  bool less_items(unsigned int i, unsigned int j) const {
    if (static_cast<unsigned int>(j) >= queue_.size()) return true;
    else if (static_cast<unsigned int>(i) >= queue_.size()) return false;
    else {
      return compare_(queue_[i], queue_[j]);
    }
  }
  //! Check if item i is less than its c'th child
  bool less_than_child(unsigned int  i, unsigned int c) const {
    int ch= child(i,c);
    return less_items(i, ch);
  }
  //! check of item i is less than its parent
  bool less_than_parent(unsigned int i) const {
    if (i==0) return false;
    else return less_items(i, parent(i));
  }

  //! Make sure that item i is not less than its parents and fix
  void bubble_down(unsigned int i) {
    CGAL_expensive_precondition(i<queue_.size());
    while (i < queue_.size()){
      unsigned int lc= child(i,FIRST), rc= child(i,SECOND);
      if (!less_items(rc, lc)){ // make the sorting stable
	if (less_items(lc, i)){
	  swap(i, lc);
	  i=lc;
	} else break;
      } else {
	if (less_items(rc, i)) {
	  swap(i, rc);
	  i=rc;
	} else break;
      }
    }
  }
  //! make sure the heap order is respected above item i
  void bubble_up(unsigned int i){
    CGAL_expensive_precondition(i< queue_.size());
    while (less_than_parent(i)){
      swap(i, parent(i));
      i=parent(i);
    }
  }
  //! both buttles
  void bubble(unsigned int i){
    bubble_up(i);
    bubble_down(i);
  }
 //! remove the front of the queue and fix the heap property
  void pop_front() {
    CGAL_expensive_precondition(!empty());
    Item_handle item= queue_.front();
    swap(0, queue_.size()-1);
    queue_.back()->set_bin(-1);
    queue_.pop_back();   
    if (!queue_.empty()) bubble_down(0);
  }

  //! debugging
  bool is_valid() const {
    for (unsigned int i=0; i< queue_.size(); ++i){
      //int bin= i;
      int back= queue_[i]->bin();
      CGAL_assertion(static_cast<unsigned int>(back)==i || write(std::cerr));
      CGAL_assertion(!less_than_parent(i) || write(std::cerr));
      //CGAL_assertion(less_than_child(i, FIRST));
      //CGAL_assertion(less_than_child(i, SECOND));
    }
    return true;
  }

};

template <class D>
std::ostream &operator<<(std::ostream &out, const Heap_pointer_event_queue<D> &q){
  q.write(out);
  return out;
}

CGAL_KDS_END_NAMESPACE;


#endif
