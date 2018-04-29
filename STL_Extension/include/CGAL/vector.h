// Copyright (c) 1997, 1998, 1999, 2000  
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
// Author(s)     : Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>

#ifndef CGAL_VECTOR_H
#define CGAL_VECTOR_H 1

#include <CGAL/basic.h>
#include <CGAL/memory.h>
#include <iterator>
#include <algorithm>
#include <memory>
#include <cstddef>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_pointer.hpp>

namespace CGAL {

namespace internal {

// We give the vector container class a class based iterator implementation.
// It ensures that iterator_traits work on compilers not supporting
// partial specializations and it guarantees that default initialization
// initializes the internal pointer to 0. Allows explicit construction 
// from a pointer.

template < class T, class Ref, class Ptr>
class vector_iterator {
private:
    Ptr ptr;
public:
    typedef vector_iterator< T, Ref, Ptr>      Self;
    typedef T                                  value_type;
    typedef Ref                                reference;
    typedef Ptr                                pointer;
    typedef std::ptrdiff_t                     difference_type;
    typedef std::random_access_iterator_tag    iterator_category;

    // CREATION
    // --------
    vector_iterator() : ptr(0) {}                // explicitly set to 0
    explicit vector_iterator( Ptr p) : ptr(p) {} // construction from pointer

    // Allows construction of const_iterator from iterator
    template < class A, class B, class C>
    vector_iterator( const vector_iterator<A,B,C>& i) : ptr(i.operator->()) {}

    // OPERATIONS Forward Category
    // ---------------------------

    bool      operator==( const Self& i) const { return ( ptr == i.ptr); }
    bool      operator!=( const Self& i) const { return !(*this == i); }
    reference operator*()                const { return *ptr; }
    pointer   operator->()               const { return ptr; }
    Self& operator++() {
        ++ptr;
        return *this;
    }
    Self  operator++(int) {
        Self tmp = *this;
        ++*this;
        return tmp;
    }

    // OPERATIONS Bidirectional Category
    // ---------------------------------
    Self& operator--() {
        --ptr;
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

    // OPERATIONS Random Access Category
    // ---------------------------------
    Self& operator+=( difference_type n) {
        ptr += n;
        return *this;
    }
    Self  operator+( difference_type n) const {
        Self tmp = *this;
        return tmp += n;
    }
    Self& operator-=( difference_type n) { return operator+=( -n); }
    Self  operator-( difference_type n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    difference_type  operator-( const Self& i) const { return ptr - i.ptr; }
    reference  operator[]( difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    bool operator< ( const Self& i) const { return ( ptr < i.ptr); }
    bool operator> ( const Self& i) const { return i < *this;    }
    bool operator<=( const Self& i) const { return !(i < *this); }
    bool operator>=( const Self& i) const { return !(*this < i); }

    vector_iterator<  T,
                      typename boost::remove_const<
                        typename boost::remove_reference<Ref>::type
                      >::type&,
                      typename boost::remove_const<
                          typename boost::remove_pointer<Ptr>::type
                      >::type* >
    remove_const() const
    {
      typedef typename boost::remove_const<
                typename boost::remove_pointer<Ptr>::type
              >::type* Ptr_no_c;
      return  vector_iterator< T,
                     typename boost::remove_const<typename boost::remove_reference<Ref>::type>::type&,
                     Ptr_no_c>
              ( const_cast<Ptr_no_c>(ptr) );
    }
};

template < class T, class Ref, class Ptr> inline
vector_iterator<T,Ref,Ptr> 
operator+( std::ptrdiff_t n, vector_iterator<T,Ref,Ptr> i) {
    return i += n;
}


template < class T, class Alloc = CGAL_ALLOCATOR(T)>
class vector {
public:
    typedef Alloc           Allocator;
    typedef Alloc           allocator_type; // STL compliant

    // Note: the standard requires the following types to be equivalent
    // to T, T*, const T*, T&, const T&, size_t, and ptrdiff_t, respectively.
    // So we don't pass these types to the iterators explicitly.
#ifdef CGAL_CXX11
  typedef typename std::allocator_traits<Allocator>::value_type            value_type;
  typedef typename std::allocator_traits<Allocator>::pointer               pointer;
  typedef typename std::allocator_traits<Allocator>::const_pointer         const_pointer;
  typedef typename std::allocator_traits<Allocator>::size_type             size_type;
  typedef typename std::allocator_traits<Allocator>::difference_type       difference_type;
#else
    typedef typename Allocator::value_type           value_type;
    typedef typename Allocator::pointer              pointer;
    typedef typename Allocator::const_pointer        const_pointer;
    typedef typename Allocator::size_type            size_type;
    typedef typename Allocator::difference_type      difference_type;
#endif

    typedef value_type&                              reference;
    typedef const value_type&                        const_reference;
    typedef std::random_access_iterator_tag          iterator_category;
    typedef vector_iterator< T, reference, pointer>  iterator;
    typedef vector_iterator< T, const_reference, const_pointer>
                                                     const_iterator;
    typedef vector< T, Alloc>                        Self;

    typedef std::reverse_iterator<iterator>          reverse_iterator;
    typedef std::reverse_iterator<const_iterator>    const_reverse_iterator;

protected:

    Allocator alloc;

    iterator start_;
    iterator finish;
    iterator end_of_storage;

    // ALLOCATION AND CONSTRUCTION HELPERS
    void construct( iterator i, const T& x) {
#ifdef CGAL_CXX11
      std::allocator_traits<Allocator>::construct(alloc,&*i, x);
#else
      alloc.construct(&*i, x);
#endif
    }
  
    void destroy( iterator i) {
#ifdef CGAL_CXX11
      std::allocator_traits<Allocator>::destroy(alloc,&*i);
#else
      alloc.destroy( &*i);
#endif
    }
  
    void destroy( iterator first, iterator last) {
        // destroy in reverse order than construction
        while ( last != first) {
            --last;
            destroy( last);
        }
    }
    void deallocate() {
        if ( start_ != iterator() )
            alloc.deallocate( &*start_, end_of_storage - start_ );
    }

protected:
    // pointer versions of begin()/end() to call the various
    // standard algorithms with the (possibly) more efficient pointers.
    pointer         pbegin()         { return &*start_; }
    const_pointer   pbegin()   const { return &*start_; }
    pointer         pend()           { return &*finish; }
    const_pointer   pend()     const { return &*finish; }

public:
    // ACCESS
    // ------
    iterator        begin()          { return start_; }
    const_iterator  begin()    const { return start_; }
    iterator        end()            { return finish; }
    const_iterator  end()      const { return finish; }
    size_type       size()     const { return size_type(end() - begin()); }
    size_type       max_size() const { return size_type(-1) / sizeof(T); }
    size_type       capacity() const {
                        return size_type(end_of_storage - start_);
    }
    bool            empty()    const { return begin() == end(); }

    reference       front()          { return *begin(); }
    const_reference front()    const { return *begin(); }
    reference       back()           { return *(end() - 1); }
    const_reference back()     const { return *(end() - 1); }
    reference       operator[] ( size_type n)       { return *(begin() + n); }
    const_reference operator[] ( size_type n) const { return *(begin() + n); }
    reference       at( size_type n)                { return *(begin() + n); }
    const_reference at( size_type n)          const { return *(begin() + n); }

    Allocator       get_allocator() const { return alloc; }

    reverse_iterator       rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(end());
    }
    reverse_iterator       rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const {
        return const_reverse_iterator(begin());
    }


    // COMPARISON
    // ----------
    bool      operator==( const Self& y) const {
        return size() == y.size() && std::equal( pbegin(), pend(), y.pbegin());
    }
    bool      operator!=( const Self& y) const { return !(*this == y); }
    bool operator< ( const Self& y) const {
        return std::lexicographical_compare( pbegin(),   pend(),
                                             y.pbegin(), y.pend());
    }
    bool operator> ( const Self& y) const { return y < *this;    }
    bool operator<=( const Self& y) const { return !(y < *this); }
    bool operator>=( const Self& y) const { return !(*this < y); }

    // CREATION
    // --------
    explicit vector()
        : start_(0), finish(0), end_of_storage(0) {}
    explicit vector( const Alloc& a)
        : start_(0), finish(0), end_of_storage(0) { alloc = a; }
    explicit vector( size_type n, const T& val) { fill_initialize(n, val); }
    explicit vector( size_type n) { fill_initialize(n, T()); }

    vector( const Self& x) {
        start_ = allocate_and_copy( x.end() - x.begin(), x.begin(), x.end());
        finish = start_ + (x.end() - x.begin());
        end_of_storage = finish;
    }

    template <class InputIterator>
    vector( InputIterator first, InputIterator last, const Alloc& a = Alloc())
        : start_(0), finish(0), end_of_storage(0)
    {
        alloc = a;
        typedef std::iterator_traits<InputIterator> Traits;
        typedef typename Traits::iterator_category  iterator_category;
        range_initialize( first, last, iterator_category());
    }

    ~vector() { 
        destroy( start_, finish);
        deallocate();
    }

    vector<T, Alloc>& operator=(const Self& x) {
        if (&x != this) {
            if ( x.size() > capacity()) {
                iterator tmp = allocate_and_copy( x.end() - x.begin(),
                                                  x.begin(),
                                                  x.end());
                destroy( start_, finish);
                deallocate();
                start_ = tmp;
                end_of_storage = start_ + (x.end() - x.begin());
            } else if (size() >= x.size()) {
                iterator i = std::copy( x.begin(), x.end(), begin());
                destroy( i, finish);
            } else {
                std::copy( x.begin(), x.begin() + size(), begin());
                std::uninitialized_copy(x.pbegin() + size(), x.pend(), pend());
            }
            finish = start_ + x.size();
        }
        return *this;
    }

    void swap( Self& x) {
        std::swap( start_, x.start_);
        std::swap( finish, x.finish);
        std::swap( end_of_storage, x.end_of_storage);
    }

    void reserve( size_type n) {
        if ( capacity() < n) {
            const size_type old_size = size();
            iterator tmp = allocate_and_copy( n, start_, finish);
            destroy(start_, finish);
            deallocate();
            start_ = tmp;
            finish = tmp + old_size;
            end_of_storage = start_ + n;
        }
    }

    // INSERTION
    // ---------
    void push_back( const T& x) {
        if ( finish != end_of_storage) {
            construct( finish, x);
            ++finish;
        } else {
            insert_aux( end(), x);
        }
    }

    iterator insert( iterator position, const T& x) {
        size_type n = position - begin();
        if (finish != end_of_storage && position == end()) {
            construct( finish, x);
            ++finish;
        } else {
            insert_aux( position, x);
        }
        return begin() + n;
    }
    iterator insert(iterator position) { return insert( position, T()); }

    template <class InputIterator>
    void insert( iterator position, InputIterator first, InputIterator last) {
        typedef std::iterator_traits<InputIterator> Traits;
        typedef typename Traits::iterator_category  iterator_category;
        range_insert( position, first, last, iterator_category());
    }
    void insert( iterator pos, size_type n, const T& x);

    // REMOVAL
    // -------
    void pop_back() {
        --finish;
        destroy( finish);
    }
    iterator erase( iterator position) {
        if (position + 1 != end())
            std::copy( position + 1, finish, position);
        --finish;
        destroy(finish);
        return position;
    }
    iterator erase( iterator first, iterator last) {
        iterator i = std::copy( last, finish, first);
        destroy( i, finish);
        finish = finish - (last - first);
        return first;
    }
    void clear() { erase( begin(), end()); }

    // ASSIGNMENT
    // ----------
    template <class InputIterator>
    void assign( InputIterator first, InputIterator last) {
        clear();
        insert( begin(), first, last);
    }
    void assign( size_type n, const T& u) {
        clear();
        insert( begin(), n, u);
    }

    void resize( size_type new_size, const T& x) {
        if (new_size < size()) 
            erase( begin() + new_size, end());
        else
            insert( end(), new_size - size(), x);
    }
    void resize( size_type new_size) { resize( new_size, T()); }

protected:
    // INTERNAL
    // --------
    void insert_aux( iterator position, const T& x);

    void fill_initialize( size_type n, const T& value) {
        start_ = allocate_and_fill(n, value);
        finish = start_ + n;
        end_of_storage = finish;
    }

    iterator allocate_and_fill( size_type n, const T& x) {
        iterator result = iterator( alloc.allocate(n));
        try {
            std::uninitialized_fill_n( &*result, n, x);
            return result;
        }
        catch(...) { 
            alloc.deallocate( &*result, n);
            throw;
        }
    }

    template <class ForwardIterator>
    iterator allocate_and_copy( size_type n,
                                ForwardIterator first,
                                ForwardIterator last) {
        iterator result = iterator( alloc.allocate(n));
        try {
            std::uninitialized_copy( first, last, &*result);
            return result;
        }
        catch(...) { 
            alloc.deallocate( &*result, n);
            throw;
        }
    }

    template <class InputIterator>
    void range_initialize(InputIterator first,
                          InputIterator last,
                          std::input_iterator_tag) {
        for ( ; first != last; ++first)
            push_back(*first);
    }

    // This function is only called by the constructor.  We have to worry
    //  about resource leaks, but not about maintaining invariants.
    template <class ForwardIterator>
    void range_initialize( ForwardIterator first,
                           ForwardIterator last,
                           std::forward_iterator_tag) {
        size_type n = std::distance( first, last);
        start_ = allocate_and_copy( n, first, last);
        finish = start_ + n;
        end_of_storage = finish;
    }

    template <class InputIterator>
    void range_insert( iterator pos,
                       InputIterator first,
                       InputIterator last,
                       std::input_iterator_tag) {
        for ( ; first != last; ++first) {
            pos = insert( pos, *first);
            ++pos;
        }
    }

    template <class ForwardIterator>
    void range_insert( iterator position,
                       ForwardIterator first,
                       ForwardIterator last,
                       std::forward_iterator_tag) {
        if (first != last) {
            size_type n = std::distance(first, last);
            if ( size_type(end_of_storage - finish) >= n) {
                const size_type elems_after = finish - position;
                iterator old_finish = finish;
                if (elems_after > n) {
                    std::uninitialized_copy( pend() - n, pend(), pend());
                    finish += n;
                    std::copy_backward( position, old_finish - n, old_finish);
                    std::copy( first, last, position);
                } else {
                    ForwardIterator mid = first;
                    std::advance( mid, elems_after);
                    std::uninitialized_copy( mid, last, pend());
                    finish += n - elems_after;
                    std::uninitialized_copy( position, old_finish, pend());
                    finish += elems_after;
                    std::copy( first, mid, position);
                }
            } else {
                const size_type old_size = size();
                const size_type len = old_size + (std::max)( old_size, n);
                iterator new_start = iterator( alloc.allocate(len));
                iterator new_finish = new_start;
                try {
                    new_finish = iterator( 
                        std::uninitialized_copy(start_, position,&*new_start));
                    new_finish = iterator(
                        std::uninitialized_copy( first, last, &*new_finish));
                    new_finish = iterator( 
                        std::uninitialized_copy(position,finish,&*new_finish));
                }
                catch(...) {
                    destroy( new_start, new_finish);
                    alloc.deallocate( &*new_start, len);
                    throw;
                }
                destroy( start_, finish);
                deallocate();
                start_ = new_start;
                finish = new_finish;
                end_of_storage = new_start + len;
            }
        }
    }
}; // class vector



template <class T, class Alloc>
inline void swap( vector<T, Alloc>& x, vector<T, Alloc>& y) {
    x.swap(y);
}

template <class T, class Alloc>
void vector<T, Alloc>::insert_aux( iterator position, const T& x) {
    if ( finish != end_of_storage) {
        construct( finish, *(finish - 1));
        ++finish;
        T x_copy = x;
        std::copy_backward( position, finish - 2, finish - 1);
        *position = x_copy;
    } else {
        const size_type old_size = size();
        const size_type len = old_size != 0 ? 2 * old_size : 1;
        iterator new_start = iterator( alloc.allocate(len));
        iterator new_finish = new_start;
        try {
            new_finish = iterator(
                std::uninitialized_copy(start_, position, &*new_start));
            construct( new_finish, x);
            ++new_finish;
            new_finish = iterator(
                std::uninitialized_copy(position,finish,&*new_finish));
        }
        catch(...) {
            destroy( new_start, new_finish); 
            alloc.deallocate( &*new_start, len);
            throw;
        }
        destroy( begin(), end());
        deallocate();
        start_ = new_start;
        finish = new_finish;
        end_of_storage = new_start + len;
    }
}


template <class T, class Alloc>
void vector<T, Alloc>::insert( iterator position, size_type n, const T& x) {
    if (n != 0) {
        if ( size_type(end_of_storage - finish) >= n) {
            T x_copy = x;
            const size_type elems_after = finish - position;
            iterator old_finish = finish;
            if (elems_after > n) {
                std::uninitialized_copy( pend() - n, pend(), pend());
                finish += n;
                std::copy_backward( position, old_finish - n, old_finish);
                std::fill( position, position + n, x_copy);
            } else {
                std::uninitialized_fill_n( pend(), n - elems_after, x_copy);
                finish += n - elems_after;
                std::uninitialized_copy( position, old_finish, pend());
                finish += elems_after;
                std::fill(position, old_finish, x_copy);
            }
        } else {
            const size_type old_size = size();        
            const size_type len = old_size + (std::max)(old_size, n);
            iterator new_start = iterator( alloc.allocate(len));
            iterator new_finish = new_start;
            try {
                new_finish = iterator(
                    std::uninitialized_copy( start_, position, &*new_start));
                std::uninitialized_fill_n( &*new_finish, n, x);
                new_finish += n;
                new_finish = iterator( 
                    std::uninitialized_copy( position, finish, &*new_finish));
            }
            catch(...) {
                destroy( new_start, new_finish);
                alloc.deallocate( &*new_start, len);
                throw;
            }
            destroy( start_, finish);
            deallocate();
            start_ = new_start;
            finish = new_finish;
            end_of_storage = new_start + len;
        }
    }
}

} // namespace internal

} //namespace CGAL

#endif // CGAL_VECTOR_H //
