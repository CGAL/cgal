// Copyright (c) 1997-2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>,
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>


#ifndef CGAL_UNION_FIND_H
#define CGAL_UNION_FIND_H

#include <CGAL/basic.h>
#include <CGAL/memory.h>
#include <cstddef>

namespace CGAL {

namespace internal {

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
    bool      operator<(const Self& other) const{ return m_p<other.m_p; }
    bool      operator>(const Self& other) const{ return m_p>other.m_p; }
    bool      operator<=(const Self& other) const{ return m_p<=other.m_p; }
    bool      operator>=(const Self& other) const{ return m_p>=other.m_p; }

    reference operator*()  const { return m_p->value; }
    pointer   operator->() const { return &(m_p->value); }
    Self&     operator++() {
                  CGAL_assertion(m_p != 0);
                  m_p = m_p->next;
                  return *this;
    }
    Self      operator++(int) {
                  Self tmp = *this;
                  ++*this;
                  return tmp;
    }
};

} // internal


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
    typedef T&                                               reference; 
    typedef const T&                                         const_reference; 

    typedef internal::UF_forward_iterator< pointer, T, T&, T*>  iterator;
    typedef iterator                                         handle;
    typedef internal::UF_forward_iterator< const_pointer, T, const T&, const T*>
                                                             const_iterator;
    typedef const_iterator                                   const_handle;

#ifdef _MSC_VER
    typedef CGAL_ALLOCATOR(Union_find_struct)                allocator;
#else
    typedef typename A::template rebind<Union_find_struct>   Rebind;
    typedef typename Rebind::other                           allocator;
#endif

private:
    pointer      m_first;
    std::size_t  sets;
    std::size_t  values;
    allocator    alloc;

    // Private decl. and no definition, thus no copy constructor
    // and no assignment for this class.
    Union_find(const Self&);
    Self& operator=(const Self&);

    pointer find( pointer p) const {
        CGAL_assertion(p != 0);
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
    const_pointer find( const_pointer p ) const {
        return find( const_cast<pointer>(p));
    }
    bool is_valid(const_handle v) const { return v != const_handle(0); }

public:
    Union_find() : m_first(0), sets(0), values(0) {}
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

    const_handle find( const_handle p ) const {
    // returns a canonical handle of the set that contains |p|,
    // i.e., |P.same_set(p,q)| iff  |P.find(p)| and |P.find(q)| return
    // the same handle. Precond: |p| is a handle in the union find structure.
        return find(p.ptr());
    }

    void unify_sets(handle p, handle q);
    // unites the sets of partition containing $p$ and $q$.
    // Precond: $p$ and $q$ are in the partition.

    bool same_set( const_handle p, const_handle q) const {
    // returns true iff $p$ and $q$ belong to the same set.
    // Precond: $p$ and $q$ are in the partition.
        return find(p) == find(q); 
    }

    iterator begin() { return iterator(m_first); }
    iterator end()   { return iterator(0); }

    const_iterator begin() const { return const_iterator(m_first); }
    const_iterator end()   const { return const_iterator(0); }
};

template <typename T, typename A>
typename Union_find<T,A>::handle Union_find<T,A>::make_set(const T& x) {
    pointer tmp = m_first;
    m_first = alloc.allocate(1);
#ifdef CGAL_CXX11
    std::allocator_traits<allocator>::construct(alloc, m_first, tmp, x);
#else
    alloc.construct( m_first, Union_find_struct(tmp,x));
#endif
    ++sets;
    ++values;
    return handle( m_first);
}

template <typename T, typename A>
void Union_find<T,A>::clear() {
    while (m_first) { 
        pointer tmp = m_first->next;
#ifdef CGAL_CXX11
        std::allocator_traits<allocator>::destroy(alloc, m_first);
#else
        alloc.destroy(m_first);
#endif
        alloc.deallocate(m_first,1);
        m_first = tmp;
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
    --sets;
}

} //namespace CGAL

#endif // CGAL_UNION_FIND_H
