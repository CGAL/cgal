// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/DS_Container.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_DS_CONTAINER_H
#define CGAL_DS_CONTAINER_H

#include <CGAL/basic.h>
#include <vector>
#include <memory>

// What we want is a kind of STL like container/allocator with the following
// properties, by decreasing priority :
// - a fast get_new_element() and release_element().
// - a compact representation, which should be N*sizeof(Elt) + o(N) (modulo
//   alignment requirements).
// - a bidirectional iterator over the (non-free = used) elements.
// - size() to know how many (non-free) elements we have.
// - find() if an element is stored in the container or not.  This doesn't
//   really need to be fast, but if we can... (used by is_cell()...).
// - would be nice to have a temporary_free_list (still active elements, but
//   which are going to be freed soon).  Probably it prevents compactness.
// - all this is expected especially when there are not so much free objects
//   compared to the allocated elements.
// - eventually something to copy this data structure, providing a way to
//   update the pointers (give access to a hash_map, at least a function that
//   converts an old pointer to the new one ?)

// Internal representation :
// - Elements are allocated by arrays of DS_Container_allocation_size, and
//   pointers to these arrays are stored in a vector<>.
// - Representation of freeness :
//   - The iterator needs to be able to test if a given element is free or not.
//     (running over the free list each time could kill performance, but maybe
//     we can only do that after a first check that there's a NULL
//     somewhere else... or maybe not NULL, but a pointer that is definitely
//     not possible, like 1, or 1 + &this_container ?).  Maybe we should store
//     union {
//       ds_cell;
//       struct { void *list_ptr; int magic_key;}
//       } Elt;
//     in which case we _must_ use the STL allocators...
//     I wonder if it's still valid C++, otherwise which destructor is run.. ?
//   - We need to have fast access to a free element.

// So, 2 possibilities in practice :
// - Elt must be able to store a list pointer, and answer if it's free or not.
//   Given some constraints in practive, this can be achieved by :
// - squatting a particular zone of the Elt which is not NULL when used to be
//   able to tell if it's free.  And squating another (or the same) to store
//   the list pointer.

// So either Elt has an additional pointer field, or it has some way to store
// the necessary information.  Both things should be provided by an
// object passed to DS_Container giving access to the list pointer, and if
// NULL, then it means the Elt is free.

// TODO :
// - Add a const_iterator
// - .merge(DS_container &) that merges two containers into one.
//   (useful for divide and conquer and co ?)
// - The magic_key is not clean, we should preferably require
//   some member functions from Elt to store the necessary information
//   (a bit like In_place_list does).

CGAL_BEGIN_NAMESPACE

class Free_elt {
    static const unsigned magic0;

    unsigned key;
    Free_elt * ptr;
public:

    Free_elt * next() const {
	return ptr;
    }

    void set_next(Free_elt *p) {
	ptr = p;
    }

    void mark_free() {
	key = magic0;
    }

    void unmark_free() {
	key = 0;
    }

    bool seems_free() const {
	return key == magic0;
    }
};

const unsigned Free_elt::magic0 = 0xc9a1c9a1;

const int DS_Container_allocation_size = 1024;
// const int DS_Container_allocation_size = 8192;

template < class DSC > class DS_Container_iterator;

template < class Elt, class Alloc = CGAL_ALLOCATOR(Elt) > //, size_t=1024> ???
class DS_Container
{
    typedef DS_Container<Elt, Alloc> Self;
public:

    typedef Elt                         value_type;
    typedef Elt&                        reference;
    typedef const Elt &                 const_reference;
    typedef DS_Container_iterator<Self> iterator;
    typedef std::size_t                 size_type;

    friend class DS_Container_iterator<Self>;

    DS_Container()
    {
	CGAL_assertion(sizeof(Free_elt) <= sizeof(Elt));
	init_free_list();
    }

    // The copy constructor and assignment operators preserve
    // the iterator order.
    DS_Container(const DS_Container &c)
    {
	std::copy(c.begin(), c.end(), std::back_inserter(*this));
    }

    DS_Container & operator=(const DS_Container &c)
    {
	if (&c != this) {
	    clear();
	    std::copy(c.begin(), c.end(), std::back_inserter(*this));
	}
	return *this;
    }

    void push_back(const Elt &e)
    {
	*get_new_element() = e;
    }

    // Compared to push_back(), we avoid a copy, and get the pointer...
    Elt * get_new_element()
    {
        if (is_free_list_empty()) {
          array_vect.push_back(alloc.allocate(DS_Container_allocation_size));
          for (int i=DS_Container_allocation_size-1; i>=0; --i)
              put_on_free_list((Free_elt *) &array_vect.back()[i]);
  	}

	Free_elt * ret = free_list.next();
	free_list.set_next(ret->next());
	ret->unmark_free();
	// alloc.construct((Elt *) ret, Elt()); // Creates/copies temporary :(
	new ((void *) ret) Elt();
	return (Elt *) ret;
    }

    void release_element(Elt *x) // erase() is the std naming ?
    {
	CGAL_assertion(!is_free(x));
	alloc.destroy(x);
	put_on_free_list((Free_elt *)x);
    }

    bool is_element(Elt *x) const
    {
	// Compare address with the beginning of each array.
	for (typename Array_vect::const_iterator it = array_vect.begin();
		it != array_vect.end(); ++it)
	    if (x >= *it && (x - *it) < DS_Container_allocation_size)
		return !is_free(x);
	return false;
    }

    iterator begin() const
    {
	return iterator(*this);
    }

    iterator end() const
    {
	return iterator(*this, 0);
    }

    ~DS_Container()
    {
	clear();
    }

    void swap(Self &c)
    {
	std::swap(alloc, c.alloc);
	std::swap(array_vect, c.array_vect);
        std::swap(free_list, c.free_list);
    }

    void clear()
    {
	for (iterator it = begin(); it != end(); ++it)
	    alloc.destroy(&*it);
	for (typename Array_vect::iterator ait = array_vect.begin();
		ait != array_vect.end(); ++ait)
	    alloc.deallocate(*ait, DS_Container_allocation_size);
	array_vect.clear();
        init_free_list();
    }

    size_type size() const
    {
	size_type number_of_free_elts = 0;
	for (Free_elt *p=free_list.next(); p!=NULL; p=p->next())
	    ++number_of_free_elts;
	return array_vect.size()*DS_Container_allocation_size
	     - number_of_free_elts;
    }

private:

    void put_on_free_list(Free_elt *x)
    {
	x->mark_free();
	x->set_next(free_list.next());
	free_list.set_next(x);
    }

    bool is_in_free_list(const Free_elt *x) const
    {
	for (Free_elt *p=free_list.next(); p!=NULL; p=p->next())
	    if (p == x)
		return true;
	return false;
    }

    bool is_free(const Elt *x) const
    {
	// seems_free() may answer true when it's not, for unlucky people.
	// Therefore we can add this _really_ expensive assertion.
	// CGAL_expensive_assertion(!((const Free_elt *) x)->seems_free()
		                 // || is_in_free_list((const Free_elt *) x));
	return ((const Free_elt *) x)->seems_free();
    }

    bool is_free_list_empty() const
    {
	return free_list.next() == NULL;
    }

    void init_free_list()
    {
	free_list.unmark_free();
	free_list.set_next(NULL);
    }

    typedef std::vector<Elt *> Array_vect;

    const Array_vect & array() const
    {
	return array_vect;
    }

    Alloc alloc;
    Array_vect array_vect;
    Free_elt free_list;
};

template < class DSC >
class DS_Container_iterator
{
    typedef DS_Container_iterator<DSC>        Self;
    typedef typename DSC::value_type          Elt;

public:
    typedef typename DSC::value_type          value_type;
    typedef value_type *                      pointer;
    typedef value_type &                      reference;
    typedef std::size_t                       size_type;
    typedef std::ptrdiff_t                    difference_type;
    typedef std::bidirectional_iterator_tag   iterator_category;

    DS_Container_iterator()
    {}

    DS_Container_iterator(const DSC &d)
	: cont(&d), vect(0), index(0)
    {
	if (cont->array().size()>0 && is_free())
	    operator++();
    }

    DS_Container_iterator(const DSC &d, int)
	: cont(&d), vect(cont->array().size()), index(0)
    {}

    bool operator==(const Self & it) const
    {
	return (it.vect == vect) && (index == it.index) && (it.cont == cont);
    }

    bool operator!=(const Self & it) const
    {
	return !(*this == it);
    }

    Self & operator++()
    {
	CGAL_assertion(vect <= cont->array().size());

	if (vect == cont->array().size()) {
	    CGAL_assertion(index == 0);
	    return *this;
	}

	do {
	    if (index < DS_Container_allocation_size-1)
		++index;
	    else {
	        index = 0;
	        ++vect;
		if (vect == cont->array().size())
		    break;
	    }
	} while (is_free()); // Skip the free Elts.
	return *this;
    }

    Self & operator--()
    {
	do {
	    if (index > 0)
		--index;
	    else {
		if (vect == 0)
		    break;
	        index = DS_Container_allocation_size-1;
	        --vect;
	    }
	} while (is_free()); // Skip the free Elts.
	return *this;
    }

    Self operator++(int)
    {
        Self tmp(*this);
        ++(*this);
        return tmp;
    }

    Self operator--(int)
    {
        Self tmp(*this);
        --(*this);
        return tmp;
    }

    Elt & operator*() const
    {
	CGAL_assertion(!is_free());
        return cont->array()[vect][index];
    }

    Elt * operator->() const
    {
	CGAL_assertion(!is_free());
        return &(cont->array()[vect][index]);
    }

private:

    bool is_free() const
    {
	return cont->is_free(&(cont->array()[vect][index]));
    }

    const DSC * cont;                         // The container
    typename DSC::Array_vect::size_type vect; // index of the current array.
    unsigned int index;                       // index in the current array.
};

CGAL_END_NAMESPACE

#endif // CGAL_DS_CONTAINER_H
