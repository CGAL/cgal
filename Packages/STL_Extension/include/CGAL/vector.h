// ============================================================================
//
// Copyright (c) 1997, 1998, 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/vector.h
// package       : $CGAL_Package: STL_Extension $
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
// maintainer    : Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
// coordinator   : ETH Zurich
//
// A vector container class: internal replacement that works with declared
// but not yet defined types as they occur with cyclic type dependencies
// with templates.
// ============================================================================

#ifndef CGAL_VECTOR_H
#define CGAL_VECTOR_H 1

#include <CGAL/basic.h>
#include <CGAL/memory.h>
#include <iterator>
#include <algorithm>
#include <memory>


CGAL_BEGIN_NAMESPACE

namespace CGALi {


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
    vector_iterator( const vector_iterator<A,B,C>& i) : ptr( &*i) {}

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
    typedef typename Allocator::value_type           value_type;
    typedef typename Allocator::pointer              pointer;
    typedef typename Allocator::const_pointer        const_pointer;
    typedef typename Allocator::reference            reference;
    typedef typename Allocator::const_reference      const_reference;
    typedef typename Allocator::size_type            size_type;
    typedef typename Allocator::difference_type      difference_type;
    typedef std::random_access_iterator_tag          iterator_category;
    typedef vector_iterator< T, reference, pointer>  iterator;
    typedef vector_iterator< T, const_reference, const_pointer>
                                                     const_iterator;
    typedef vector< T, Alloc>                        Self;

#if defined(__sun) && defined(__SUNPRO_CC)
    typedef std::reverse_iterator< iterator,
                                   typename iterator::iterator_category,
                                   typename iterator::value_type,
                                   typename iterator::reference,
                                   typename iterator::pointer,
                                   typename iterator::difference_type
                                   > reverse_iterator;
    typedef std::reverse_iterator< const_iterator,
                                   typename const_iterator::iterator_category,
                                   typename const_iterator::value_type,
                                   typename const_iterator::reference,
                                   typename const_iterator::pointer,
                                   typename const_iterator::difference_type
                                   > const_reverse_iterator;
#else
    typedef std::reverse_iterator< iterator >       reverse_iterator;
    typedef std::reverse_iterator< const_iterator > const_reverse_iterator;
#endif // defined(__sun) && defined(__SUNPRO_CC)

protected:
#ifndef _MSC_VER
    // Somehow the static initialization does not work correctly for MSVC
    // ---> strange linker errors
    static
#endif // _MSC_VER
    Allocator alloc;

    iterator start;
    iterator finish;
    iterator end_of_storage;

    // ALLOCATION AND CONSTRUCTION HELPERS
    void construct( iterator i, const T& x) { alloc.construct( &*i, x);}
    void destroy( iterator i) { alloc.destroy( &*i); }
    void destroy( iterator first, iterator last) {
        // destroy in reverse order than construction
        while ( last != first) {
            --last;
            destroy( last);
        }
    }
    void deallocate() {
        if ( &*start)
            alloc.deallocate( &*start, end_of_storage - start );
    }

public:
    // ACCESS
    // ------
    iterator        begin()          { return start; }
    const_iterator  begin()    const { return start; }
    iterator        end()            { return finish; }
    const_iterator  end()      const { return finish; }
    size_type       size()     const { return size_type(end() - begin()); }
    size_type       max_size() const { return size_type(-1) / sizeof(T); }
    size_type       capacity() const {
                        return size_type(end_of_storage - start);
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
        return size() == y.size() && std::equal( begin(), end(), y.begin());
    }
    bool      operator!=( const Self& y) const { return !(*this == y); }
    bool operator< ( const Self& y) const {
        return std::lexicographical_compare( begin(),   end(),
                                             y.begin(), y.end());
    }
    bool operator> ( const Self& y) const { return y < *this;    }
    bool operator<=( const Self& y) const { return !(y < *this); }
    bool operator>=( const Self& y) const { return !(*this < y); }

    // CREATION
    // --------
    explicit vector()
        : start(0), finish(0), end_of_storage(0) {}
    explicit vector( const Alloc& a)
        : start(0), finish(0), end_of_storage(0) { alloc = a; }
    explicit vector( size_type n, const T& val) { fill_initialize(n, val); }
    explicit vector( size_type n) { fill_initialize(n, T()); }

    vector( const Self& x) {
        start = allocate_and_copy( x.end() - x.begin(), x.begin(), x.end());
        finish = start + (x.end() - x.begin());
        end_of_storage = finish;
    }

    template <class InputIterator>
    vector( InputIterator first, InputIterator last, const Alloc& a = Alloc())
        : start(0), finish(0), end_of_storage(0)
    {
        alloc = a;
        typedef std::iterator_traits<InputIterator> Traits;
        typedef typename Traits::iterator_category  iterator_category;
        range_initialize( first, last, iterator_category());
    }

    ~vector() { 
        destroy( start, finish);
        deallocate();
    }

    vector<T, Alloc>& operator=(const Self& x) {
        if (&x != this) {
            if ( x.size() > capacity()) {
                iterator tmp = allocate_and_copy( x.end() - x.begin(),
                                                  x.begin(),
                                                  x.end());
                destroy( start, finish);
                deallocate();
                start = tmp;
                end_of_storage = start + (x.end() - x.begin());
            } else if (size() >= x.size()) {
                iterator i = std::copy( x.begin(), x.end(), begin());
                destroy( i, finish);
            } else {
                std::copy( x.begin(), x.begin() + size(), start);
                std::uninitialized_copy( x.begin() + size(), x.end(), finish);
            }
            finish = start + x.size();
        }
        return *this;
    }

    void swap( Self& x) {
        std::swap( start, x.start);
        std::swap( finish, x.finish);
        std::swap( end_of_storage, x.end_of_storage);
    }

    void reserve( size_type n) {
        if ( capacity() < n) {
            const size_type old_size = size();
            iterator tmp = allocate_and_copy( n, start, finish);
            destroy(start, finish);
            deallocate();
            start = tmp;
            finish = tmp + old_size;
            end_of_storage = start + n;
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
        start = allocate_and_fill(n, value);
        finish = start + n;
        end_of_storage = finish;
    }

    iterator allocate_and_fill( size_type n, const T& x) {
        iterator result = iterator( alloc.allocate(n));
        try {
            std::uninitialized_fill_n( result, n, x);
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
            std::uninitialized_copy( first, last, result);
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
        start = allocate_and_copy( n, first, last);
        finish = start + n;
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
                    std::uninitialized_copy( finish - n, finish, finish);
                    finish += n;
                    std::copy_backward( position, old_finish - n, old_finish);
                    std::copy( first, last, position);
                } else {
                    ForwardIterator mid = first;
                    std::advance( mid, elems_after);
                    std::uninitialized_copy( mid, last, finish);
                    finish += n - elems_after;
                    std::uninitialized_copy( position, old_finish, finish);
                    finish += elems_after;
                    std::copy( first, mid, position);
                }
            } else {
                const size_type old_size = size();
                const size_type len = old_size + std::max( old_size, n);
                iterator new_start = iterator( alloc.allocate(len));
                iterator new_finish = new_start;
                try {
                    new_finish = std::uninitialized_copy( start, position,
                                                          new_start);
                    new_finish = std::uninitialized_copy( first, last,
                                                          new_finish);
                    new_finish = std::uninitialized_copy( position, finish,
                                                          new_finish);
                }
                catch(...) {
                    destroy( new_start, new_finish);
                    alloc.deallocate( &*new_start, len);
                    throw;
                }
                destroy( start, finish);
                deallocate();
                start = new_start;
                finish = new_finish;
                end_of_storage = new_start + len;
            }
        }
    }
}; // class vector

#ifndef _MSC_VER
// init static member allocator object
template <class T, class Alloc>
Alloc vector< T, Alloc>::alloc = Alloc();
#endif // _MSC_VER


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
            new_finish = std::uninitialized_copy( start, position, new_start);
            construct( new_finish, x);
            ++new_finish;
            new_finish = std::uninitialized_copy(position, finish, new_finish);
        }
        catch(...) {
            destroy( new_start, new_finish); 
            alloc.deallocate( &*new_start, len);
            throw;
        }
        destroy( begin(), end());
        deallocate();
        start = new_start;
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
                std::uninitialized_copy( finish - n, finish, finish);
                finish += n;
                std::copy_backward( position, old_finish - n, old_finish);
                std::fill( position, position + n, x_copy);
            } else {
                std::uninitialized_fill_n( finish, n - elems_after, x_copy);
                finish += n - elems_after;
                std::uninitialized_copy( position, old_finish, finish);
                finish += elems_after;
                std::fill(position, old_finish, x_copy);
            }
        } else {
            const size_type old_size = size();        
            const size_type len = old_size + std::max(old_size, n);
            iterator new_start = iterator( alloc.allocate(len));
            iterator new_finish = new_start;
            try {
                new_finish = std::uninitialized_copy( start, position,
                                                      new_start);
                new_finish = std::uninitialized_fill_n( new_finish, n, x);
                new_finish = std::uninitialized_copy( position, finish,
                                                      new_finish);
            }
            catch(...) {
                destroy( new_start, new_finish);
                alloc.deallocate( &*new_start, len);
                throw;
            }
            destroy( start, finish);
            deallocate();
            start = new_start;
            finish = new_finish;
            end_of_storage = new_start + len;
        }
    }
}

} // namespace CGALi

CGAL_END_NAMESPACE


#endif // CGAL_VECTOR_H //
// EOF //
