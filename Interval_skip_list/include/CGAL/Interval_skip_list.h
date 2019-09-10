// Copyright (c) 2003 GeometryFactory
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Andreas Fabri


#ifndef CGAL_INTERVAL_SKIP_LIST_H
#define CGAL_INTERVAL_SKIP_LIST_H

#include <CGAL/license/Interval_skip_list.h>


#include <CGAL/basic.h>
#include <list>
#include <iostream>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <boost/random/variate_generator.hpp>


//#define CGAL_ISL_USE_CCC
#define CGAL_ISL_USE_LIST

#if defined(CGAL_ISL_USE_CCC) || ! defined(CGAL_ISL_USE_LIST)
#include <CGAL/Compact_container.h>
#endif

namespace CGAL {

  class Interval;
  template <class Interval_>
  class IntervalList;

  template <class Interval_>
  class Interval_skip_list;

  template <class Interval_>
  class IntervalListElt;

  template <class Interval_>
  class IntervalSLnode;

  const int MAX_FORWARD = 48; 	// Maximum number of forward pointers



  template <class Interval_>
  class IntervalSLnode  // interval skip list node
  {
    typedef Interval_ Interval;
    typedef typename Interval::Value Value;
    bool is_header;
    typedef Interval* Interval_handle;

    Value key;
    IntervalSLnode** forward;  // array of forward pointers
    IntervalList<Interval>**   markers;  // array of interval markers, 
                                         // one for each pointer
    IntervalList<Interval>* eqMarkers;   // markers for node itself
    int ownerCount;  // number of interval end points with value equal to key
    int topLevel;  // index of top level of forward pointers in this node.
    // Levels are numbered 0..topLevel.
  public:
    friend class Interval_skip_list<Interval>;

    IntervalSLnode(const Value& searchKey, int levels);  // constructor
    IntervalSLnode(int levels);  // constructor for the header

    IntervalSLnode* get_next();


    void print(std::ostream& os) const;

    const Value& 
    getValue()
    {
      return key;
    }
    
    // number of levels of this node
    int 
    level() const 
    {
      return(topLevel+1);
    }
  
    bool
    isHeader() const
    {
      return is_header;
    }

    void deleteMarks(IntervalList<Interval>* l);
    
    ~IntervalSLnode();  // destructor
  };


template <class Interval_>
class Interval_for_container : public Interval_
    {
      private:
      void * p;
    public:
      Interval_for_container(const Interval_& i)
	: Interval_(i), p(NULL)
      {}
      
      void *   for_compact_container() const { return p; }
      void * & for_compact_container()       { return p; }
    };


  template <class Interval_>
  class Interval_skip_list
  {
  private:
    typedef Interval_ Interval;
    typedef typename Interval::Value Value;
    boost::rand48 random;

#ifdef CGAL_ISL_USE_LIST
    std::list<Interval> container;

    typedef typename std::list<Interval>::iterator Interval_handle;
#else
    Compact_container<Interval_for_container<Interval> > container;
    typedef typename Compact_container<Interval_for_container<Interval> >::iterator 
      Interval_handle;
#endif

#ifdef CGAL_ISL_USE_CCC
    typedef typename Compact_container<IntervalListElt<Interval> >::iterator ILE_handle;
#else
    typedef IntervalListElt<Interval>* ILE_handle;
#endif

    int maxLevel;
    IntervalSLnode<Interval>* header;

    int randomLevel();  // choose a new node level at random

    // place markers for Interval I.  I must have been inserted in the list.
    // left is the left endpoint of I and right is the right endpoint if I.
    // *** needs to be fixed:
    void placeMarkers(IntervalSLnode<Interval>* left, 
		      IntervalSLnode<Interval>* right, 
		      const Interval_handle& I);


    // remove markers for Interval I
    void removeMarkers(const Interval_handle& I);  


    // adjust markers after insertion of x with update vector "update"
    void adjustMarkersOnInsert(IntervalSLnode<Interval>* x, 
			       IntervalSLnode<Interval>** update);


    // adjust markers to prepare for deletion of x, which has update vector
    // "update"
    void adjustMarkersOnDelete(IntervalSLnode<Interval>* x, 
			       IntervalSLnode<Interval>** update);


    // remove node x, which has updated vector update.
    void remove(IntervalSLnode<Interval>* x, 
		IntervalSLnode<Interval>** update);


    // remove markers for Interval I starting at left, the left endpoint
    // of I, and and stopping at the right endpoint of I.
    Interval_handle removeMarkers(IntervalSLnode<Interval>* left, 
				  const Interval& I);


    // Remove markers for interval m from the edges and nodes on the
    // level i path from l to r.
    void removeMarkFromLevel(const Interval& m, int i,
			     IntervalSLnode<Interval> *l, 
			     IntervalSLnode<Interval>* r);


    // Search for search key, and return a pointer to the 
    // intervalSLnode x found, as well as setting the update vector 
    // showing pointers into x. 
    IntervalSLnode<Interval>* search(const Value& searchKey, 
				     IntervalSLnode<Interval>** update);

  
    // insert a new single value 
    // into list, returning a pointer to its location. 
    IntervalSLnode<Interval>* insert(const Value& searchKey);


    // insert an interval into list 
    void insert(const Interval_handle& I);

  public:

    friend class IntervalSLnode<Interval>;

    Interval_skip_list();  
    
    template <class InputIterator>
    Interval_skip_list(InputIterator b, InputIterator e)
    {
      maxLevel = 0;
      header = new IntervalSLnode<Interval>(MAX_FORWARD);
      for (int i = 0; i< MAX_FORWARD; i++) {
	header->forward[i] = 0;
      }
      for(; b!= e; ++b){
	insert(*b);
      }
    }


    ~Interval_skip_list(); 


    void clear();

    int size() const 
    {
      return container.size();
    }

 
    // return node containing
    // Value if found, otherwise null
    IntervalSLnode<Interval>* search(const Value& searchKey); 



    template <class OutputIterator>
    OutputIterator 
    find_intervals(const Value& searchKey, OutputIterator out )
    {
      IntervalSLnode<Interval>* x = header;
      for(int i=maxLevel; 
	  i >= 0 && (x->isHeader() || (x->key != searchKey)); i--) {
	while (x->forward[i] != 0 && (searchKey >= x->forward[i]->key)) {
	  x = x->forward[i];
	}
	// Pick up markers on edge as you drop down a level, unless you are at 
	// the searchKey node already, in which case you pick up the
	// eqMarkers just prior to exiting loop.
	if(!x->isHeader() && (x->key != searchKey)) {
	  out = x->markers[i]->copy(out);  
	} else if (!x->isHeader()) { // we're at searchKey
	  out = x->eqMarkers->copy(out);
	}
      }
      return out;
    }
    
    bool
    is_contained(const Value& searchKey) const
    {
      IntervalSLnode<Interval>* x = header;
      for(int i=maxLevel; 
	  i >= 0 && (x->isHeader() || (x->key != searchKey)); i--) {
	while (x->forward[i] != 0 && (searchKey >= x->forward[i]->key)) {
	  x = x->forward[i];
	}
	// Pick up markers on edge as you drop down a level, unless you are at 
	// the searchKey node already, in which case you pick up the
	// eqMarkers just prior to exiting loop.
	if(!x->isHeader() && (x->key != searchKey)) {
	  return true;  
	} else if (!x->isHeader()) { // we're at searchKey
	  return true;
	}
      }
      return false;
    }
    

    void insert(const Interval& I);

    template <class InputIterator>
    int insert(InputIterator b, InputIterator e)
    {
      int i = 0;
      for(; b!= e; ++b){
	insert(*b);
	++i;
      }
      return i;
    }


    bool remove(const Interval& I);  // delete an interval from list
    void print(std::ostream& os) const;
    void printOrdered(std::ostream& os) const;

#ifdef CGAL_ISL_USE_LIST
    typedef typename std::list<Interval>::const_iterator const_iterator;
#else
    typedef typename 
    Compact_container<Interval_for_container<Interval> >::const_iterator
                                                             const_iterator;
#endif
    const_iterator begin() const
    {
      return container.begin();
    }
    const_iterator end() const
    {
      return container.end();
    }
    
  };



  template <class Interval_>
  class IntervalList
  {
    typedef Interval_ Interval;
    typedef typename Interval::Value Value;
    //typedef Interval* Interval_handle;
#ifdef CGAL_ISL_USE_LIST

    typedef typename std::list<Interval>::iterator Interval_handle;
#else
    typedef typename Compact_container<Interval_for_container<Interval> >::iterator 
      Interval_handle;
#endif

#ifdef CGAL_ISL_USE_CCC
    typedef typename Compact_container<IntervalListElt<Interval> >::iterator ILE_handle;
#else
    typedef IntervalListElt<Interval>* ILE_handle;
#endif


    ILE_handle header;

#ifdef CGAL_ISL_USE_CCC
    static Compact_container<IntervalListElt<Interval> > compact_container;
#endif
    typedef std::allocator<IntervalListElt<Interval> > Alloc;
    Alloc alloc;

  public:

    friend class IntervalListElt<Interval>;

    IntervalList();

    void insert(const Interval_handle& I);

    bool remove(const Interval& I, Interval_handle& res);

    void remove(const Interval& I);

    void removeAll(IntervalList* l);


    ILE_handle create_list_element(const Interval_handle& I)
    {
#ifdef CGAL_ISL_USE_CCC
      IntervalListElt<Interval> e(I);
      ILE_handle it = compact_container.insert(e);
      return it;
#else
      IntervalListElt<Interval> *elt_ptr = alloc.allocate(1);
#ifdef CGAL_CXX11
      std::allocator_traits<Alloc>::construct(alloc,elt_ptr, I);
#else      
      alloc.construct(elt_ptr, I);
#endif
      return elt_ptr;
      //return new IntervalListElt<Interval>(I);
#endif
    }

    void erase_list_element(ILE_handle I)
    {      
#ifdef CGAL_ISL_USE_CCC
      compact_container.erase(I);
#else

#ifdef CGAL_CXX11
      std::allocator_traits<Alloc>::destroy(alloc,I);
#else
      alloc.destroy(I);
#endif
      alloc.deallocate(I,1);
      //delete I;
#endif
    }

    ILE_handle get_first();

    ILE_handle get_next(ILE_handle element);

    void copy(IntervalList* from); // add contents of "from" to self
 
 
    template <class OutputIterator>
    OutputIterator
    copy(OutputIterator out) const
    {
      ILE_handle e = header;
      while(e!= NULL) { 
	out = *(e->I);
	++out;
	e = e->next;
      }
      return out;
    }
    
    bool contains(const Interval_handle& I) const;

    void clear();  // delete elements of self to make self an empty list.

    void print(std::ostream& os) const;

    ~IntervalList();
  };

#ifdef CGAL_ISL_USE_CCC
  template <class Interval_>
  Compact_container<IntervalListElt<Interval_> > 
     IntervalList<Interval_>::compact_container;
#endif
 



  template <class Interval_>
  class IntervalListElt
  {
    typedef Interval_ Interval;
#ifdef CGAL_ISL_USE_LIST

    typedef typename std::list<Interval>::iterator Interval_handle;
#else
    typedef typename Compact_container<Interval_for_container<Interval> >::iterator 
      Interval_handle;
#endif

#ifdef CGAL_ISL_USE_CCC
    typedef typename Compact_container<IntervalListElt<Interval> >::iterator ILE_handle;
#else
    typedef IntervalListElt<Interval>* ILE_handle;
#endif

    Interval_handle I;
    ILE_handle next;
#ifdef CGAL_ISL_USE_CCC
    void* p;
#endif
  public:
#ifdef CGAL_ISL_USE_CCC
    void *   for_compact_container() const { return p; }
    void * & for_compact_container()       { return p; }
#endif

    bool operator==(const IntervalListElt& e)
    {
      return ( ((*I) == (*(e.I))) && (next == e.next));
    }
    
    friend class IntervalList<Interval>;


    IntervalListElt();

    IntervalListElt(const Interval_handle& anInterval);

    ~IntervalListElt();

    void 
    set_next(ILE_handle nextElt)
    {
      next = nextElt;
    }

    ILE_handle get_next()
    {
      return next;
    }

    Interval_handle getInterval()
    {
      return I;
    }

    void print(std::ostream& os) const;
  };


  template <class Interval>
  IntervalSLnode<Interval>::IntervalSLnode(const Value& searchKey, int levels)
    : is_header(false)
  {
    // levels is actually one less than the real number of levels
    key = searchKey;
    topLevel = levels;
    forward = new IntervalSLnode*[levels+1];
    markers = new IntervalList<Interval>*[levels+1];
    eqMarkers = new IntervalList<Interval>();
    ownerCount = 0;
    for(int i=0; i<=levels; i++) {
      forward[i] = 0;
      // initialize an empty interval list
      markers[i] = new IntervalList<Interval>(); 
    }
  }


  template <class Interval>
  IntervalSLnode<Interval>::IntervalSLnode(int levels)
    : is_header(true)
  {
    // levels is actually one less than the real number of levels
    topLevel = levels;
    forward = new IntervalSLnode*[levels+1];
    markers = new IntervalList<Interval>*[levels+1];
    eqMarkers = new IntervalList<Interval>();
    ownerCount = 0;
    for(int i=0; i<=levels; i++) {
      forward[i] = 0;
      // initialize an empty interval list
      markers[i] = new IntervalList<Interval>(); 
    }
  }

  template <class Interval>
  Interval_skip_list<Interval>::Interval_skip_list()
  {
    maxLevel = 0;
    header = new IntervalSLnode<Interval>(MAX_FORWARD);
    for (int i = 0; i< MAX_FORWARD; i++) {
      header->forward[i] = 0;
    }
  }

  template <class Interval>
  Interval_skip_list<Interval>::~Interval_skip_list()
  {
    while(header != 0){
      IntervalSLnode<Interval>* next = header->get_next();
      delete header;
      header = next;
    }
  }

  template <class Interval>
  void
  Interval_skip_list<Interval>::clear()
  {
    while(header != 0){
      IntervalSLnode<Interval>* next = header->get_next();
      delete header;
      header = next;
    }
    header = new IntervalSLnode<Interval>(MAX_FORWARD);
    for (int i = 0; i< MAX_FORWARD; i++) {
      header->forward[i] = 0;
    }
    maxLevel = 0;

  }

  template <class Interval>
  IntervalSLnode<Interval>* IntervalSLnode<Interval>::get_next()
  {
    return(forward[0]);
  }
  template <class Interval>
  void Interval_skip_list<Interval>::print(std::ostream& os) const
  {
    os << "\nAn Interval_skip_list:  \n";
    os << "|container| == " << container.size() << std::endl;
    IntervalSLnode<Interval>* n = header->get_next();

    while( n != 0 ) {
      n->print(os);
      n = n->get_next();
    }
  }

  template <class Interval>
  std::ostream& operator<<(std::ostream& os, 
			   const Interval_skip_list<Interval>& isl)
  {
    isl.print(os);
    return os;
  }


  template <class Interval>
  void Interval_skip_list<Interval>::printOrdered(std::ostream& os) const
  {
    IntervalSLnode<Interval>* n = header->get_next();
    os << "values in list:  ";
    while( n != 0 ) {
      os << n->key << " ";
      n = n->get_next();
    }
    os << std::endl;
  }

template <class Interval>
  void IntervalList<Interval>::copy(IntervalList* from)
  {
    ILE_handle e = from->header;
    while(e!=NULL) { 
      insert(e->I);
      e = e->next;
    }
  }


template <class Interval>
  void IntervalList<Interval>::clear()
  {
    ILE_handle x = header;
    ILE_handle y; 
    while(x!= NULL) { // was 0
      y = x;
      x = x->next;
      erase_list_element(y);
    }
    header=0;
  }

  template <class Interval>
  IntervalSLnode<Interval>* 
  Interval_skip_list<Interval>::insert(const Value& searchKey)
  {
    // array for maintaining update pointers 
    IntervalSLnode<Interval>* update[MAX_FORWARD]; 
    IntervalSLnode<Interval>* x;
    int i;

    // Find location of searchKey, building update vector indicating
    // pointers to change on insertion.
    x = search(searchKey,update);
    if(x==0 || (x->key != searchKey)) {
      // put a new node in the list for this searchKey
      int newLevel = randomLevel();
      if (newLevel > maxLevel){
	for(i=maxLevel+1; i<=newLevel; i++){
	  update[i] = header;
	  header->markers[i]->clear();
	}
	maxLevel = newLevel;
      }
      x = new IntervalSLnode<Interval>(searchKey, newLevel);

      // add x to the list
      for(i=0; i<=newLevel; i++) {
	x->forward[i] = update[i]->forward[i];
	update[i]->forward[i] = x;
      }

      // adjust markers to maintain marker invariant
      this->adjustMarkersOnInsert(x,update);
    }
    // else, the searchKey is in the list already, and x points to it.
    return(x);
  }



  // Adjust markers on this IS-list to maintain marker invariant now that
  // node x has just been inserted, with update vector `update.'

  template <class Interval>
  void 
  Interval_skip_list<Interval>::adjustMarkersOnInsert
                                    (IntervalSLnode<Interval>* x,
				     IntervalSLnode<Interval>** update)
  {
    // Phase 1:  place markers on edges leading out of x as needed.

    // Starting at bottom level, place markers on outgoing level i edge of x.
    // If a marker has to be promoted from level i to i+1 of higher, place it
    // in the promoted set at each step.

    IntervalList<Interval> promoted;  
    // list of intervals that identify markers being
    // promoted, initially empty.

    IntervalList<Interval> newPromoted; 
    // temporary set to hold newly promoted markers.
    
    IntervalList<Interval> removePromoted;  
    // holding place for elements to be removed  from promoted list.

    IntervalList<Interval> tempMarkList;  // temporary mark list
    ILE_handle m;
    int i;

    for(i=0; (i<= x->level() - 2) && x->forward[i+1]!=0; i++) {
      IntervalList<Interval>* markList = update[i]->markers[i];
      for(m = markList->get_first(); m != NULL ; m = markList->get_next(m)) {
	if(m->getInterval()->contains_interval(x->key,x->forward[i+1]->key)) { 
	  // promote m
	  
	  // remove m from level i path from x->forward[i] to x->forward[i+1]
	  removeMarkFromLevel(*m->getInterval(),
			      i,
			      x->forward[i],
			      x->forward[i+1]);
	  // add m to newPromoted
	  newPromoted.insert(m->getInterval());
	} else {
	  // place m on the level i edge out of x
	  x->markers[i]->insert(m->getInterval());
	  // do *not* place m on x->forward[i]; it must already be there. 
	}
      }
      
      for(m = promoted.get_first(); m != NULL; m = promoted.get_next(m)) {
	if(!m->getInterval()->contains_interval(x->key, x->forward[i+1]->key)){
	  // Then m does not need to be promoted higher.
	  // Place m on the level i edge out of x and remove m from promoted.
	  x->markers[i]->insert(m->getInterval());
	  // mark x->forward[i] if needed
	  if(m->getInterval()->contains(x->forward[i]->key))
	    x->forward[i]->eqMarkers->insert(m->getInterval());
	  removePromoted.insert(m->getInterval());
	} else { 
	  // continue to promote m
	  // Remove m from the level i path from x->forward[i]
	  // to x->forward[i+1].
	  removeMarkFromLevel(*(m->getInterval()),
			      i,
			      x->forward[i],
			      x->forward[i+1]);
	}
      }
      promoted.removeAll(&removePromoted);
      removePromoted.clear();
      promoted.copy(&newPromoted);
      newPromoted.clear();
    }
    // Combine the promoted set and updated[i]->markers[i]
    // and install them as the set of markers on the top edge out of x
    // that is non-null.  
    
    x->markers[i]->copy(&promoted);
    x->markers[i]->copy(update[i]->markers[i]);
    for(m=promoted.get_first(); m!=NULL; m=promoted.get_next(m))
      if(m->getInterval()->contains(x->forward[i]->key))
        x->forward[i]->eqMarkers->insert(m->getInterval());
    
    // Phase 2:  place markers on edges leading into x as needed.
    
    // Markers on edges leading into x may need to be promoted as high as
    // the top edge coming into x, but never higher.
    
    promoted.clear();
    
    for (i=0; (i <= x->level() - 2) && !update[i+1]->isHeader(); i++) {
      tempMarkList.copy(update[i]->markers[i]);
      for(m = tempMarkList.get_first(); 
	  m != NULL; 
	  m = tempMarkList.get_next(m)){
	if(m->getInterval()->contains_interval(update[i+1]->key,x->key)) {
	  // m needs to be promoted
	  // add m to newPromoted
	  newPromoted.insert(m->getInterval());
	  
	  // Remove m from the path of level i edges between updated[i+1]
	  // and x (it will be on all those edges or else the invariant
	  // would have previously been violated.
	  removeMarkFromLevel(*(m->getInterval()),i,update[i+1],x);
	}
      }
      tempMarkList.clear();  // reclaim storage
      
      for(m = promoted.get_first(); m != NULL; m = promoted.get_next(m)) {
	if (!update[i]->isHeader() && 
	    m->getInterval()->contains_interval(update[i]->key,x->key) &&
	    !update[i+1]->isHeader() &&
	    ! m->getInterval()->contains_interval(update[i+1]->key,x->key) ) {
	  // Place m on the level i edge between update[i] and x, and
	  // remove m from promoted.
	  update[i]->markers[i]->insert(m->getInterval());
	  // mark update[i] if needed
	  if(m->getInterval()->contains(update[i]->key))
	    update[i]->eqMarkers->insert(m->getInterval());
	  removePromoted.insert(m->getInterval());
	} else {
	  // Strip m from the level i path from update[i+1] to x.
	  removeMarkFromLevel(*(m->getInterval()),i,update[i+1],x);
	}
	
      }
      // remove non-promoted marks from promoted
      promoted.removeAll(&removePromoted);
      removePromoted.clear();  // reclaim storage
      
      // add newPromoted to promoted and make newPromoted empty
      promoted.copy(&newPromoted);
      newPromoted.clear();     
    }
    
    /* Assertion:  i=x->level()-1 OR update[i+1] is the header.
       
       If i=x->level()-1 then either x has only one level, or the top-level
       pointer into x must not be from the header, since otherwise we would
       have stopped on the previous iteration.  If x has 1 level, then
       promoted is empty.  If x has 2 or more levels, and i!=x->level()-1,
       then the edge on the next level up (level i+1) is from the header.  In
       any of these cases, all markers in the promoted set should be
       deposited on the current level i edge into x.  An edge out of the
       header should never be marked.  Note that in the case where x has only
       1 level, we try to copy the contents of the promoted set onto the
       marker set of the edge out of the header into x at level i=0, but of
       course, the promoted set will be empty in this case, so no markers
       will be placed on the edge.  */

    update[i]->markers[i]->copy(&promoted);
    for(m=promoted.get_first(); m!=NULL; m=promoted.get_next(m))
      if(m->getInterval()->contains(update[i]->key))
	update[i]->eqMarkers->insert(m->getInterval());

    // Place markers on x for all intervals the cross x.
    // (Since x is a new node, every marker comming into x must also leave x).
    for(i=0; i<x->level(); i++)
      x->eqMarkers->copy(x->markers[i]);
    
    promoted.clear(); // reclaim storage
    
  } // end adjustMarkersOnInsert

  template <class Interval>
  void
  Interval_skip_list<Interval>::adjustMarkersOnDelete
                                        (IntervalSLnode<Interval>* x,
					 IntervalSLnode<Interval>** update)
  {
    // x is node being deleted.  It is still in the list.
    // update is the update vector for x.
    IntervalList<Interval> demoted;
    IntervalList<Interval> newDemoted;
    IntervalList<Interval> tempRemoved;
    ILE_handle m;
    int i;
    IntervalSLnode<Interval> *y;

    // Phase 1:  lower markers on edges to the left of x as needed.

    for(i=x->level()-1; i>=0; i--){
      // find marks on edge into x at level i to be demoted
      for(m=update[i]->markers[i]->get_first(); m!=NULL; 
	  m=update[i]->markers[i]->get_next(m)){
	if(x->forward[i]==0 ||
	   ! m->getInterval()->contains_interval(update[i]->key,
						 x->forward[i]->key)){
	  newDemoted.insert(m->getInterval());
	}
      }
      // Remove newly demoted marks from edge.
      update[i]->markers[i]->removeAll(&newDemoted);
      // NOTE:  update[i]->eqMarkers is left unchanged because any markers
      // there before demotion must be there afterwards.

      // Place previously demoted marks on this level as needed.
      for(m=demoted.get_first(); m!=NULL; m=demoted.get_next(m)){
	// Place mark on level i from update[i+1] to update[i], not including 
	// update[i+1] itself, since it already has a mark if it needs one.
	for(y=update[i+1]; y!=0 && y!=update[i]; y=y->forward[i]) {
	  if (y!=update[i+1] && m->getInterval()->contains(y->key)) 
	    y->eqMarkers->insert(m->getInterval());
	  y->markers[i]->insert(m->getInterval());
	}
	if(y!=0 && y!=update[i+1] && m->getInterval()->contains(y->key)) 
	  y->eqMarkers->insert(m->getInterval());

	// if this is the lowest level m needs to be placed on,
	// then place m on the level i edge out of update[i]
	// and remove m from the demoted set.
	if(x->forward[i]!=0 &&
	   m->getInterval()->contains_interval(update[i]->key,
					       x->forward[i]->key))
	  {
	    update[i]->markers[i]->insert(m->getInterval());
	    tempRemoved.insert(m->getInterval());
	  }
      }
      demoted.removeAll(&tempRemoved);
      tempRemoved.clear();
      demoted.copy(&newDemoted);
      newDemoted.clear();
    }

    // Phase 2:  lower markers on edges to the right of D as needed
  
    demoted.clear();
    // newDemoted is already empty

    for(i=x->level()-1; i>=0; i--){
      for(m=x->markers[i]->get_first(); m!=NULL ; m=x->markers[i]->get_next(m)){
	if(x->forward[i]!=0 && 
	   (update[i]->isHeader() ||
	    !m->getInterval()->contains_interval(update[i]->key,
						 x->forward[i]->key)))
	  {
	    newDemoted.insert(m->getInterval());
	  }
      }

      for(m=demoted.get_first(); m!= NULL; m=demoted.get_next(m)){
	// Place mark on level i from x->forward[i] to x->forward[i+1].
	// Don't place a mark directly on x->forward[i+1] since it is already
	// marked.
	for(y=x->forward[i];y!=x->forward[i+1];y=y->forward[i]){
	  y->eqMarkers->insert(m->getInterval());
	  y->markers[i]->insert(m->getInterval());
	}

	if(x->forward[i]!=0 && !update[i]->isHeader() &&
	   m->getInterval()->contains_interval(update[i]->key,
					       x->forward[i]->key))
	  {
	    tempRemoved.insert(m->getInterval());
	  }
      }
      demoted.removeAll(&tempRemoved);
      demoted.copy(&newDemoted);
      newDemoted.clear();
    }
  }  // end adjustMarkersOnDelete

  template <class Interval>
  IntervalSLnode<Interval>::~IntervalSLnode()
  {
    for(int i = 0; i<=topLevel; i++)
      delete markers[i];
    delete [] forward;
    delete [] markers;
    delete eqMarkers;
  }

  template <class Interval>
  bool Interval_skip_list<Interval>::remove(const Interval& I)
  {
    // arrays for maintaining update pointers 
    IntervalSLnode<Interval>* update[MAX_FORWARD]; 

    IntervalSLnode<Interval>* left = search(I.inf(),update);
    if(left==0 || left->ownerCount <= 0) {
      return false;
    }

    Interval_handle ih = removeMarkers(left,I);
    container.erase(ih);
    left->ownerCount--;
    if(left->ownerCount == 0) remove(left,update);

    // Note:  we search for right after removing left since some
    // of left's forward pointers may point to right.  We don't
    // want any pointers of update vector pointing to a node that is gone.

    IntervalSLnode<Interval>* right = search(I.sup(),update);
    if(right==0 || right->ownerCount <= 0) {
      return false;
    }
    right->ownerCount--;
    if(right->ownerCount == 0) remove(right,update);
    return true;
  }

  template <class Interval>
  void 
  Interval_skip_list<Interval>::remove(IntervalSLnode<Interval>* x, 
				     IntervalSLnode<Interval>** update)
  {
    // Remove interval skip list node x.  The markers that the interval
    // x belongs to have already been removed.

    adjustMarkersOnDelete(x,update);

    // now splice out x.
    for(int i=0; i<=x->level()-1; i++)
      update[i]->forward[i] = x->forward[i];

    // and finally deallocate it
    delete x;
  }


  template <class Interval>
  IntervalSLnode<Interval>* 
  Interval_skip_list<Interval>::search(const Value& searchKey)
  {
    IntervalSLnode<Interval>* x = header;
    for(int i=maxLevel; i >= 0; i--) {
      while (x->forward[i] != 0 && (x->forward[i]->key < searchKey)) {
	x = x->forward[i];
      }
    }
    x = x->forward[0];
    if(x != NULL && (x->key == searchKey))
      return(x);
    else
      return(NULL);
  }

  template <class Interval>
  IntervalSLnode<Interval>* 
  Interval_skip_list<Interval>::search(const Value& searchKey, 
				     IntervalSLnode<Interval>** update)
  {
    IntervalSLnode<Interval>* x = header;
    // Find location of searchKey, building update vector indicating
    // pointers to change on insertion.
    for(int i=maxLevel; i >= 0; i--) {
      while (x->forward[i] != 0 && (x->forward[i]->key < searchKey)) {
	x = x->forward[i];
      }
      update[i] = x;
    }
    x = x->forward[0];
    return(x);
  }



  template <class Interval>
  void Interval_skip_list<Interval>::insert(const Interval_handle& I)
    // insert an interval into list
  {
    // insert end points of interval
    IntervalSLnode<Interval>* left = this->insert(I->inf());
    IntervalSLnode<Interval>* right = this->insert(I->sup());
    left->ownerCount++;
    right->ownerCount++;

    // place markers on interval
    this->placeMarkers(left,right,I);
  }

  template <class Interval>
  void
  Interval_skip_list<Interval>::insert(const Interval& I)
  {
#ifdef CGAL_ISL_USE_LIST
    container.push_front(I);
    Interval_handle ih = container.begin();
#else
    Interval_for_container<Interval> ifc(I);
    Interval_handle ih = container.insert(ifc);
#endif
    insert(ih);
  }


  template <class Interval>
  void 
  Interval_skip_list<Interval>::placeMarkers(IntervalSLnode<Interval>* left, 
					   IntervalSLnode<Interval>* right, 
					   const Interval_handle& I)
  {
    // Place markers for the interval I.  left is the left endpoint
    // of I and right is the right endpoint of I, so it isn't necessary
    // to search to find the endpoints.

    IntervalSLnode<Interval>* x = left;
    if (I->contains(x->key)) x->eqMarkers->insert(I);
    int i = 0;  // start at level 0 and go up
    while(x->forward[i]!=0 && I->contains_interval(x->key,x->forward[i]->key)){
      // find level to put mark on
      while(i!=x->level()-1 
            && x->forward[i+1] != 0
            && I->contains_interval(x->key,x->forward[i+1]->key))
	i++;
      // Mark current level i edge since it is the highest edge out of
      // x that contains I, except in the case where current level i edge
      // is null, in which case it should never be marked.
      if (x->forward[i] != 0) { 
	x->markers[i]->insert(I);  
	x = x->forward[i];
	// Add I to eqMarkers set on node unless currently at right endpoint
	// of I and I doesn't contain right endpoint.
	if (I->contains(x->key)) x->eqMarkers->insert(I);
      }
    }

    // mark non-ascending path
    while(x->key != right->key) {
      // find level to put mark on
      while(i!=0 && (x->forward[i] == 0 || 
		     !I->contains_interval(x->key,x->forward[i]->key)))
	i--;
      // At this point, we can assert that i=0 or x->forward[i]!=0 and 
      // I contains 
      // (x->key,x->forward[i]->key).  In addition, x is between left and 
      // right so i=0 implies I contains (x->key,x->forward[i]->key).
      // Hence, the interval must be marked.  Note that it is impossible
      // for us to be at the end of the list because x->key is not equal
      // to right->key.
      x->markers[i]->insert(I);
      x = x->forward[i];
      if (I->contains(x->key)) x->eqMarkers->insert(I);     
    }
  }  // end placeMarkers

  template <class Interval>
  typename Interval_skip_list<Interval>::Interval_handle 
  Interval_skip_list<Interval>::removeMarkers(IntervalSLnode<Interval>* left, 
					    const Interval& I)
  {
    // Remove markers for interval I, which has left as it's left
    // endpoint,  following a staircase pattern.

    //    Interval_handle res=0, tmp=0; // af: assignment not possible with std::list
    Interval_handle res, tmp;
    // remove marks from ascending path
    IntervalSLnode<Interval>* x = left;
    if (I.contains(x->key)) {
      if(x->eqMarkers->remove(I, tmp)){
	res = tmp;
      }
    }
    int i = 0;  // start at level 0 and go up
    while(x->forward[i]!=0 && I.contains_interval(x->key,x->forward[i]->key)) {
      // find level to take mark from
      while(i!=x->level()-1 
            && x->forward[i+1] != 0
            && I.contains_interval(x->key,x->forward[i+1]->key))
	i++;
      // Remove mark from current level i edge since it is the highest edge out
      // of x that contains I, except in the case where current level i edge
      // is null, in which case there are no markers on it.
      if (x->forward[i] != 0) { 
	if(x->markers[i]->remove(I, tmp)){
	  res = tmp;
	}
	x = x->forward[i];
	// remove I from eqMarkers set on node unless currently at right 
	// endpoint of I and I doesn't contain right endpoint.
	if (I.contains(x->key)){
	  if(x->eqMarkers->remove(I, tmp)){
	    res = tmp;
	  }
	}
      }
    }

    // remove marks from non-ascending path
    while(x->key != I.sup()) {
      // find level to remove mark from
      while(i!=0 && (x->forward[i] == 0 || 
		     ! I.contains_interval(x->key,x->forward[i]->key)))
	i--;
      // At this point, we can assert that i=0 or x->forward[i]!=0 and 
      // I contains 
      // (x->key,x->forward[i]->key).  In addition, x is between left and 
      // right so i=0 implies I contains (x->key,x->forward[i]->key).
      // Hence, the interval is marked and the mark must be removed.  
      // Note that it is impossible for us to be at the end of the list 
      // because x->key is not equal to right->key.
      if(x->markers[i]->remove(I, tmp)){
	res = tmp;
      }
      x = x->forward[i];
      if (I.contains(x->key)){
	if(x->eqMarkers->remove(I, tmp)){
	  res = tmp;
	}     
      }
    }
    CGAL_assertion(*res == I);
    return res;
  }

  template <class Interval>
  void 
  Interval_skip_list<Interval>::removeMarkFromLevel(const Interval& m, int i,
						  IntervalSLnode<Interval> *l, 
						  IntervalSLnode<Interval>* r)
  {
    IntervalSLnode<Interval> *x;
    for(x=l; x!=0 && x!=r; x=x->forward[i]) {
      x->markers[i]->remove(m);
      x->eqMarkers->remove(m);
    }
    if(x!=0) x->eqMarkers->remove(m);
  }


  template <class Interval>
  int
  Interval_skip_list<Interval>::randomLevel()
  {
    boost::geometric_distribution<> proba(0.5);
    boost::variate_generator<boost::rand48&, boost::geometric_distribution<> > die(random, proba);
    
    return (std::min)(die(), (int)maxLevel)+1;
}


  template <class Interval>
  void IntervalSLnode<Interval>::print(std::ostream& os) const
  {
    int i;
    os << "IntervalSLnode key:  ";
    if (! is_header) {
     os << key;
    }else {
      os << "HEADER";
    }
    os << "\n";
    os << "number of levels: " << level() << std::endl;
    os << "owning intervals:  ";
    os << "ownerCount = " << ownerCount << std::endl;
    os <<  std::endl;
    os << "forward pointers:\n";
    for(i=0; i<=topLevel; i++)
      {
	os << "forward[" << i << "] = ";
	if(forward[i] != NULL) {
	  os << forward[i]->getValue();
	} else {
	  os << "NULL";
	}
	os << std::endl;
      }
    os << "markers:\n";
    for(i=0; i<=topLevel; i++)
      {
	os << "markers[" << i << "] = ";
	if(markers[i] != NULL) {
	  markers[i]->print(os);
	} else {
	  os << "NULL";
	}
	os << "\n";
      }
    os << "EQ markers:  ";
    eqMarkers->print(os);
    os << std::endl << std::endl;
  }


  template <class Interval>
  void IntervalList<Interval>::insert(const Interval_handle& I)
  {
    ILE_handle temp = create_list_element(I);
    temp->next = header;
    header = temp;
  }


  template <class Interval>
  inline
  bool
  IntervalList<Interval>::remove(const Interval& I, Interval_handle& res)
  {
    ILE_handle x, last;
    x = header; last = NULL;
    while(x != NULL && *(x->getInterval()) != I) {
      last = x;
      x = x->next;
    } 
    if(x==NULL) {
      return false;
    } else if (last==NULL) {
      header = x->next;
      res = x->getInterval();
      erase_list_element(x);
    } else {
      last->next = x->next;
      res = x->getInterval();
      erase_list_element(x);
    }
    return true;
  }


  template <class Interval>
  void
  IntervalList<Interval>::remove(const Interval& I)
  {
    ILE_handle x, last;
    x = header; last = NULL;
    while(x != NULL && *(x->getInterval()) != I) {
      last = x;
      x = x->next;
    }
    if(x==NULL) {
      return ;
    } else if (last==NULL) {
      header = x->next;
      erase_list_element(x);
    } else {
      last->next = x->next;
      erase_list_element(x);
    }
  }

  template <class Interval>
  void IntervalList<Interval>::removeAll(IntervalList<Interval> *l)
  {
    ILE_handle x;
    for(x=l->get_first(); x!=NULL; x=l->get_next(x))
      this->remove(*(x->getInterval()));
  }

  // We need the default constructor for the compact_container
  template <class Interval>
  inline 
  IntervalListElt<Interval>::IntervalListElt()
    : next(NULL)
  {}


  template <class Interval>
  inline 
  IntervalListElt<Interval>::IntervalListElt(const Interval_handle& anInterval)
    : I(anInterval), next(NULL)
  {}

  template <class Interval>
  inline 
  IntervalListElt<Interval>::~IntervalListElt()
  {}


  template <class Interval>
  inline
  typename IntervalList<Interval>::ILE_handle
  IntervalList<Interval>::get_first()
  {
    return header;
  }

  template <class Interval>
  inline
  
  typename IntervalList<Interval>::ILE_handle
  IntervalList<Interval>::get_next(ILE_handle element)
  {
    return element->next;
  }

  template <class Interval>
  void IntervalList<Interval>::print(std::ostream& os) const
  {
    ILE_handle e = header;
    while(e != NULL) {
      e->print(os);
      e = e->get_next();
    }
  }

  template <class Interval>
  void IntervalListElt<Interval>::print(std::ostream& /*os*/) const
  {
    /*
    if(I == 0) {
      os << "NULL";
    } else {
      os << *I;
    }
    */
  }

  template <class Interval>
  inline 
  bool IntervalList<Interval>::contains(const Interval_handle& I) const
  {
    ILE_handle x = header;
    while(x!=0 && I != x->I)
      x = x->next;
    if (x==NULL)
      return false;
    else
      return true;
  }



  template <class Interval>
  inline IntervalList<Interval>::IntervalList()
    :  header(NULL)
  {}


  template <class Interval>
  inline IntervalList<Interval>::~IntervalList()
  {
    this->clear();
  }


} // namespace CGAL

#endif // CGAL_INTERVAL_SKIP_LIST_H
