// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Union_find.h
// package       : Union_find
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>,
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
// maintainer    : Lutz Kettner <kettner@mpi-sb.mpg.de>
// coordinator   : Lutz Kettner <kettner@mpi-sb.mpg.de>
//
// implementation: Union-find implemented with path compression
// ============================================================================


#ifndef CGAL_UNION_FIND_H
#define CGAL_UNION_FIND_H

#include <CGAL/basic.h>
#include <CGAL/memory.h>


CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class PTR_, class V_, class R_, class P_>
class UF_forward_iterator {
    PTR_   m_p;
public:
    // should be private and Union_find<...> a friend.
    PTR_   ptr() const { return m_p; }  

    typedef UF_forward_iterator<PTR_, V_, R_, P_> Self;
    typedef V_                                    value_type;
    typedef R_                                    reference;
    typedef P_                                    pointer;
    typedef std::forward_iterator_tag             iterator_category;

    UF_forward_iterator() : m_p(0) {}
    UF_forward_iterator(PTR_ p) : m_p(p) {}

    // Allows construction of const_iterator from iterator
    template <class PPTR, class VV, class RR, class PP>
    UF_forward_iterator(const UF_forward_iterator<PPTR, VV, RR, PP>& i)
        : m_p(i.ptr()) {}

    bool      operator==( const Self& x) const { return m_p == x.m_p; }
    bool      operator!=( const Self& x) const { return !(*this == x); }
    reference operator*()  const { return m_p->value; }
    pointer   operator->() const { return &(m_p->value); }
    Self&     operator++() {
                  CGAL_assertion(m_p);
                  m_p = m_p->next;
                  return *this;
    }
    Self      operator++(int) {
                  Self tmp = *this;
                  ++*this;
                  return tmp;
    }
};

} // CGALi


// Union-Find with path compression.
// --------------------------------------------------------------------
// An instance of the data type Union_find is a partition of values of
// type |T| into disjoint sets. The type |A| has to be a model of the
// allocator concept as defined in the C++ standard.

// Union_find is implemented with union by rank and path compression.
// The running time for $m$ set operations on $n$ elements is 
// $O(n \alpha(m,n))$ where $\alpha(m,n)$ is the extremely slow
// growing inverse of Ackermann's function.

template <typename T, typename A = CGAL_ALLOCATOR(T) > 
class Union_find {
    struct Union_find_struct {
        typedef Union_find_struct* pointer;
        // friend class Union_find<T,A>;
        mutable pointer      up;
        pointer              next;
        std::size_t          size;
        T                    value;
        Union_find_struct( pointer n, const T& x)
            : up(0), next(n), size(1), value(x) {}
    };

public:
    typedef Union_find<T,A>                                  Self;
    typedef Union_find_struct*                               pointer;
    typedef const Union_find_struct*                         const_pointer;

    typedef T                                                value_type;
    typedef CGALi::UF_forward_iterator< pointer, T, T&, T*>  iterator;
    typedef iterator                                         handle;
    typedef CGALi::UF_forward_iterator< const_pointer, T, const T&, const T*>
                                                             const_iterator;
    typedef const_iterator                                   const_handle;

#ifdef _MSC_VER
    typedef CGAL_ALLOCATOR(Union_find_struct)                allocator;
#else
    typedef typename A::template rebind<Union_find_struct>   Rebind;
    typedef typename Rebind::other                           allocator;
#endif

private:
    pointer      first;
    std::size_t  sets;
    std::size_t  values;
    allocator    alloc;

    // Private decl. and no definition, thus no copy constructor
    // and no assignment for this class.
    Union_find(const Self&);
    Self& operator=(const Self&);

    pointer find( pointer p) const {
        CGAL_assertion(p);
        pointer r = p;
        while (r->up) 
            r = r->up; // now r is the root;
        while (p->up) {
            pointer u = p->up;
            p->up = r; // path compression: assign root r as new parent
            p = u;     // this would need the 'mutable' for the up field
        }              // if we would have a const_pointer, see the cast
        return r;      // in the fct. below. We keep the mutable as reminder.
    }
    const_pointer find( const_pointer p) const {
        return find( const_cast<pointer>(p));
    }
    bool is_valid(const_handle v) const { return v != const_handle(0); }

public:
    Union_find() : first(0), sets(0), values(0) {}
    ~Union_find() { clear(); }

    allocator   get_allocator() const { return alloc; }

    std::size_t number_of_sets() const { return sets; }
    // returns the number of disjoint sets

    std::size_t size() const { return values; }
    // returns the number of values

    std::size_t bytes() const {
    // returns the memory consumed
        return values * sizeof(Union_find_struct) + sizeof( Self); 
    }

    std::size_t size( const_handle p) const { return find(p).ptr()->size; }
    // returns the size of the set containing $p$

    void clear();
    // reinitializes to an empty partition

    handle make_set(const T& x);
    // creates a new singleton set containing |x| and returns a handle to it

    handle push_back(const T& x) { return make_set(x); }
    // same as |make_set(x)|

    template <class Forward_iterator>
    void insert( Forward_iterator first, Forward_iterator beyond) {
    // insert the range of values referenced by |[first,beyond)|.
    // Precond: value type of |Forward_iterator| is |T|.
        while (first != beyond)
            push_back(*first++);
    }

    handle       find( handle p)       const { return find(p.ptr()); }
    const_handle find( const_handle p) const { return find(p.ptr()); }
    // returns a canonical handle of the set that contains |p|,
    // i.e., |P.same_set(p,q)| iff  |P.find(p)| and |P.find(q)| return
    // the same handle. Precond: |p| is a handle in the union find structure.

    void unify_sets(handle p, handle q);
    // unites the sets of partition containing $p$ and $q$.
    // Precond: $p$ and $q$ are in the partition.

    bool same_set( const_handle p, const_handle q) const {
    // returns true iff $p$ and $q$ belong to the same set.
    // Precond: $p$ and $q$ are in the partition.
        return find(p) == find(q); 
    }

    iterator begin() { return iterator(first); }
    iterator end()   { return iterator(0); }

    const_iterator begin() const { return const_iterator(first); }
    const_iterator end()   const { return const_iterator(0); }
};

template <typename T, typename A>
typename Union_find<T,A>::handle Union_find<T,A>::make_set(const T& x) {
    pointer tmp = first;
    first = alloc.allocate(1);
    alloc.construct( first, Union_find_struct(tmp,x));
    ++sets;
    ++values;
    return handle( first);
}

template <typename T, typename A>
void Union_find<T,A>::clear() {
    while (first) { 
        pointer tmp = first->next;
        alloc.destroy(first);
        alloc.deallocate(first,1);
        first = tmp;
    }
    sets   = 0;
    values = 0;
}

template <typename T, typename A>
void Union_find<T,A>::unify_sets( handle p, handle q) {
    CGAL_assertion( is_valid(p) && is_valid(q));
    pointer pit = find( p.ptr());
    pointer qit = find( q.ptr());
    if (pit == qit)
        return;
    std::size_t sp = pit->size;
    std::size_t sq = qit->size;
    if (sp > sq) 
        std::swap(pit,qit); // now sp <= sq
    pit->up = qit;  // linking roots
    qit->size += pit->size; // updating size
    ++sets;
}

CGAL_END_NAMESPACE

#endif // CGAL_UNION_FIND_H

