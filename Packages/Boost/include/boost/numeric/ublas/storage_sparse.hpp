//
//  Copyright (c) 2000-2002
//  Joerg Walter, Mathias Koch
//
//  Permission to use, copy, modify, distribute and sell this software
//  and its documentation for any purpose is hereby granted without fee,
//  provided that the above copyright notice appear in all copies and
//  that both that copyright notice and this permission notice appear
//  in supporting documentation.  The authors make no representations
//  about the suitability of this software for any purpose.
//  It is provided "as is" without express or implied warranty.
//
//  The authors gratefully acknowledge the support of
//  GeNeSys mbH & Co. KG in producing this work.
//

#ifndef BOOST_UBLAS_STORAGE_SPARSE_H
#define BOOST_UBLAS_STORAGE_SPARSE_H

#include <algorithm>
#include <map>
#include <set>

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/iterator.hpp>
#include <boost/numeric/ublas/storage.hpp>

namespace boost { namespace numeric { namespace ublas {

    namespace detail {

        template<class I, class T, class C>
        BOOST_UBLAS_INLINE
        I lower_bound (const I &begin, const I &end, const T &t, C compare) {
            // t <= *begin <=> ! (*begin < t)
            if (begin == end || ! compare (*begin, t))
                return begin;
            if (compare (*(end - 1), t))
                return end;
            return std::lower_bound (begin, end, t, compare);
        }
        template<class I, class T, class C>
        BOOST_UBLAS_INLINE
        I upper_bound (const I &begin, const I &end, const T &t, C compare) {
            if (begin == end || compare (t, *begin))
                return begin;
            // (*end - 1) <= t <=> ! (t < *end)
            if (! compare (t, *(end - 1)))
                return end;
            return std::upper_bound (begin, end, t, compare);
        }

        template<class P>
        struct less_pair {
            BOOST_UBLAS_INLINE
            bool operator () (const P &p1, const P &p2) {
                return p1.first < p2.first;
            }
        };
        template<class T>
        struct less_triple {
            BOOST_UBLAS_INLINE
            bool operator () (const T &t1, const T &t2) {
                return t1.first.first < t2.first.first ||
                       (t1.first.first == t2.first.first && t1.first.second < t2.first.second);
            }
        };

    }

#ifdef BOOST_UBLAS_STRICT_STORAGE_SPARSE
    template<class D>
    struct sparse_storage_element_traits {
        typedef typename D::index_type index_type;
        typedef typename D::data_const_reference data_const_reference;
        typedef typename D::data_reference data_reference;
    };
    template<>
    struct sparse_storage_element_traits<float> {
        typedef std::size_t index_type;
        typedef void data_const_reference;
        typedef void data_reference;
    };
    template<>
    struct sparse_storage_element_traits<double> {
        typedef std::size_t index_type;
        typedef void data_const_reference;
        typedef void data_reference;
    };
#ifdef BOOST_UBLAS_USE_LONG_DOUBLE
    template<>
    struct sparse_storage_element_traits<long double> {
        typedef std::size_t index_type;
        typedef void data_const_reference;
        typedef void data_reference;
    };
#endif
    template<>
    struct sparse_storage_element_traits<std::complex<float> > {
        typedef std::size_t index_type;
        typedef void data_const_reference;
        typedef void data_reference;
    };
    template<>
    struct sparse_storage_element_traits<std::complex<double> > {
        typedef std::size_t index_type;
        typedef void data_const_reference;
        typedef void data_reference;
    };
#ifdef BOOST_UBLAS_USE_LONG_DOUBLE
    template<>
    struct sparse_storage_element_traits<std::complex<long double> > {
        typedef std::size_t index_type;
        typedef void data_const_reference;
        typedef void data_reference;
    };
#endif

    template<class A>
    class sparse_storage_element:
       public container_reference<A> {
    public:
        typedef A array_type;
        typedef typename A::index_type index_type;
        typedef typename A::data_value_type data_value_type;
        // typedef const data_value_type &data_const_reference;
        typedef typename type_traits<data_value_type>::const_reference data_const_reference;
        typedef data_value_type &data_reference;
        typedef std::pair<index_type, data_value_type> value_type;
        typedef std::pair<index_type, data_value_type> *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        sparse_storage_element (array_type &a, pointer it):
            container_reference<array_type> (a), it_ (it), i_ (it->first), d_ (it->second), dirty_ (false) {}
        BOOST_UBLAS_INLINE
        sparse_storage_element (array_type &a, index_type i):
            container_reference<array_type> (a), it_ (), i_ (i), d_ (), dirty_ (false) {
            pointer it = (*this) ().find (i_);
            if (it == (*this) ().end ())
                it = (*this) ().insert ((*this) ().end (), value_type (i_, d_));
            d_ = it->second;
        }
        BOOST_UBLAS_INLINE
        ~sparse_storage_element () {
            if (dirty_) {
                if (! it_)
                    it_ = (*this) ().find (i_);
                BOOST_UBLAS_CHECK (it_ != (*this) ().end (), internal_logic ());
                it_->second = d_;
            }
        }

        // Element access
        BOOST_UBLAS_INLINE
        typename sparse_storage_element_traits<data_value_type>::data_const_reference
        operator [] (typename sparse_storage_element_traits<data_value_type>::index_type i) const {
            return d_ [i];
        }
#ifdef BOOST_UBLAS_DEPRECATED
        BOOST_UBLAS_INLINE
        typename sparse_storage_element_traits<data_value_type>::data_reference
        operator [] (typename sparse_storage_element_traits<data_value_type>::index_type i) {
            dirty_ = true;
            return d_ [i];
        }
#endif

        // Assignment
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator = (const D &d) {
            d_ = d;
            dirty_ = true;
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator = (const sparse_storage_element &p) {
            d_ = p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator += (const D &d) {
            d_ += d;
            dirty_ = true;
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator += (const sparse_storage_element &p) {
            d_ += p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator -= (const D &d) {
            d_ -= d;
            dirty_ = true;
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator -= (const sparse_storage_element &p) {
            d_ -= p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator *= (const D &d) {
            d_ *= d;
            dirty_ = true;
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator *= (const sparse_storage_element &p) {
            d_ *= p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator /= (const D &d) {
            d_ /= d;
            dirty_ = true;
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator /= (const sparse_storage_element &p) {
            d_ /= p.d_;
            dirty_ = true;
            return *this;
        }

        // Conversion
        BOOST_UBLAS_INLINE
        operator data_const_reference () const {
            return d_;
        }
#ifdef BOOST_UBLAS_DEPRECATED
        BOOST_UBLAS_INLINE
        operator data_reference () {
            dirty_ = true;
            return d_;
        }
#endif

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (sparse_storage_element p) {
            if (this != &p) {
                dirty_ = true;
                p.dirty_ = true;
                std::swap (d_, p.d_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (sparse_storage_element p1, sparse_storage_element p2) {
            p1.swap (p2);
        }
#endif

    private:
        pointer it_;
        index_type i_;
        data_value_type d_;
        bool dirty_;
    };
#endif

    // Map array
    template<class I, class T>
    class map_array {
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef I index_type;
        typedef T data_value_type;
        // typedef const T &data_const_reference;
        typedef typename type_traits<T>::const_reference data_const_reference;
#ifndef BOOST_UBLAS_STRICT_STORAGE_SPARSE
        typedef T &data_reference;
#else
        typedef sparse_storage_element<map_array> data_reference;
#endif
        typedef std::pair<I, T> value_type;
        typedef const std::pair<I, T> &const_reference;
        typedef std::pair<I, T> &reference;
        typedef const std::pair<I, T> *const_pointer;
        typedef std::pair<I, T> *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        map_array ():
            capacity_ (0), data_ (new value_type [0]), size_ (0) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            std::fill (data_, data_ + size_, value_type ());
        }
        BOOST_UBLAS_INLINE
        map_array (no_init):
            capacity_ (0), data_ (new value_type [0]), size_ (0) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
        }
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        map_array (size_type size):
            capacity_ (size), data_ (new value_type [size]), size_ (0) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            std::fill (data_, data_ + size_, value_type ());
        }
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        map_array (size_type size, no_init):
            capacity_ (size), data_ (new value_type [size]), size_ (0) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
        }
        BOOST_UBLAS_INLINE
        map_array (const map_array &a):
            capacity_ (a.size_), data_ (new value_type [a.size_]), size_ (a.size_) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            *this = a;
        }
        BOOST_UBLAS_INLINE
        ~map_array () {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            delete [] data_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size) {
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
            if (size > capacity_) {
                pointer data = new value_type [size << 1];
                // Assuming std compliant allocator as requested during review.
                // if (! data)
                //     throw std::bad_alloc ();
                // if (! data_)
                //     throw std::bad_alloc ();
                std::copy (data_, data_ + std::min (size, size_), data);
                std::fill (data + std::min (size, size_), data + size, value_type ());
                delete [] data_;
                capacity_ = size << 1;
                data_ = data;
            }
            size_ = size;
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
        }

        // Reserving
        BOOST_UBLAS_INLINE
        void reserve (size_type capacity) {
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
            if (capacity > capacity_) {
                pointer data = new value_type [capacity];
                // Assuming std compliant allocator as requested during review.
                // if (! data)
                //     throw std::bad_alloc ();
                // if (! data_)
                //     throw std::bad_alloc ();
                std::copy (data_, data_ + size_, data);
                delete [] data_;
                capacity_ = capacity;
                data_ = data;
            }
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
        }

        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        data_reference operator [] (index_type i) {
#ifndef BOOST_UBLAS_STRICT_STORAGE_SPARSE
            pointer it = find (i);
            if (it == end ())
                it = insert (end (), value_type (i, data_value_type ()));
            BOOST_UBLAS_CHECK (it != end (), internal_logic ());
            return it->second;
#else
            return data_reference (*this, i);
#endif
        }

        // Assignment
        BOOST_UBLAS_INLINE
        map_array &operator = (const map_array &a) {
            // Too unusual semantic.
            // Thanks to Michael Stevens for spotting this.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                resize (a.size_);
                std::copy (a.data_, a.data_ + a.size_, data_);
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        map_array &assign_temporary (map_array &a) {
            swap (a);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (map_array &a) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                std::swap (capacity_, a.capacity_);
                std::swap (data_, a.data_);
                std::swap (size_, a.size_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (map_array &a1, map_array &a2) {
            a1.swap (a2);
        }
#endif

        // Element insertion and deletion
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer push_back (pointer it, const value_type &p) {
            if (size () == 0 || (it = end () - 1)->first < p.first) {
                resize (size () + 1);
                *(it = end () - 1) = p;
                return it;
            }
            // Raising exceptions abstracted as requested during review.
            // throw external_logic ();
            external_logic ().raise ();
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer insert (pointer it, const value_type &p) {
            it = detail::lower_bound (begin (), end (), p, detail::less_pair<value_type> ());
            difference_type n = it - begin ();
            BOOST_UBLAS_CHECK (size () == 0 || size () == size_type (n) || it->first != p.first, external_logic ());
            resize (size () + 1);
            it = begin () + n;
            std::copy_backward (it, end () - 1, end ());
            *it = p;
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void insert (pointer it, pointer it1, pointer it2) {
#ifdef BOOST_UBLAS_BOUNDS_CHECK
            while (it1 != it2) {
                insert (it, *it1);
                ++ it1;
            }
#else
            difference_type n = it - begin ();
            resize (size () + it2 - it1);
            it = begin () + n;
            std::copy (it1, it2, it);
            std::sort (begin (), end (), detail::less_pair<value_type> ());
#endif
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void erase (pointer it) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            // Fixed by George Katsirelos.
            // (*it).second = data_value_type ();
            std::copy (it + 1, end (), it);
            resize (size () - 1);
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void erase (pointer it1, pointer it2) {
            BOOST_UBLAS_CHECK (begin () <= it1 && it1 < it2 && it2 <= end (), bad_index ());
            // Fixed by George Katsirelos.
            // while (it1 != it2) {
            //     BOOST_UBLAS_CHECK (begin () <= it1 && it1 < end (), bad_index ());
            //     (*it1).second = data_value_type ();
            //     ++ it1;
            // }
            std::copy (it2, end (), it1);
            resize (size () - (it2 - it1));
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void clear () {
            resize (0);
        }

        // Element lookup
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_pointer find (index_type i) const {
#ifdef BOOST_UBLAS_DEPRECATED
            std::pair<const_pointer, const_pointer> pit;
            pit = std::equal_range (begin (), end (), value_type (i, data_value_type ()), detail::less_pair<value_type> ());
            if (pit.first->first == i)
                return pit.first;
            else if (pit.second->first == i)
                return pit.second;
            else
                return end ();
#else
            const_pointer it (detail::lower_bound (begin (), end (), value_type (i, data_value_type ()), detail::less_pair<value_type> ()));
            if (it == end () || it->first != i)
                it = end ();
            return it;
#endif
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer find (index_type i) {
#ifdef BOOST_UBLAS_DEPRECATED
            std::pair<pointer, pointer> pit;
            pit = std::equal_range (begin (), end (), value_type (i, data_value_type ()), detail::less_pair<value_type> ());
            if (pit.first->first == i)
                return pit.first;
            else if (pit.second->first == i)
                return pit.second;
            else
                return end ();
#else
            pointer it (detail::lower_bound (begin (), end (), value_type (i, data_value_type ()), detail::less_pair<value_type> ()));
            if (it == end () || it->first != i)
                it = end ();
            return it;
#endif
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_pointer lower_bound (index_type i) const {
            return detail::lower_bound (begin (), end (), value_type (i, data_value_type ()), detail::less_pair<value_type> ());
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer lower_bound (index_type i) {
            return detail::lower_bound (begin (), end (), value_type (i, data_value_type ()), detail::less_pair<value_type> ());
        }
#ifdef BOOST_UBLAS_DEPRECATED
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_pointer upper_bound (index_type i) const {
            return detail::upper_bound (begin (), end (), value_type (i, data_value_type ()), detail::less_pair<value_type> ());
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer upper_bound (index_type i) {
            return detail::upper_bound (begin (), end (), value_type (i, data_value_type ()), detail::less_pair<value_type> ());
        }
#endif

        // Iterators simply are pointers.

        typedef const_pointer const_iterator;

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return data_ + size_;
        }

        typedef pointer iterator;

        BOOST_UBLAS_INLINE
        iterator begin () {
            return data_;
        }
        BOOST_UBLAS_INLINE
        iterator end () {
            return data_ + size_;
        }

        // Reverse iterators

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef std::reverse_iterator<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator rbegin () const {
            return const_reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator rend () const {
            return const_reverse_iterator (begin ());
        }

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef std::reverse_iterator<iterator, value_type, reference> reverse_iterator;
#else
        typedef std::reverse_iterator<iterator> reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        reverse_iterator rbegin () {
            return reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator rend () {
            return reverse_iterator (begin ());
        }

    private:
        size_type capacity_;
        pointer data_;
        size_type size_;
    };

    namespace detail {
        using namespace boost::numeric::ublas;

#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
        template<class A>
        struct map_traits {};

        template<class I, class T>
        struct map_traits<map_array<I, T> > {
            typedef typename map_array<I, T>::data_reference reference;
        };

        template<class I, class T>
        struct map_traits<std::map<I, T> > {
            typedef typename std::map<I, T>::mapped_type &reference;

        };
#endif

        // Some helpers for map_array

        template<class I, class T>
        BOOST_UBLAS_INLINE
        void reserve (map_array<I, T> &a, typename map_array<I, T>::size_type capacity) {
            a.reserve (capacity);
        }

#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
        template<class I, class T>
        BOOST_UBLAS_INLINE
        typename map_array<I, T>::data_reference make_reference (map_array<I, T> &a, typename map_array<I, T>::iterator it) {
            return reference (a, it);
        }
#endif

        // Some helpers for std::map

        template<class I, class T>
        BOOST_UBLAS_INLINE
        void reserve (std::map<I, T> &/* a */, typename std::map<I, T>::size_type /* capacity */) {}

#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
        template<class I, class T>
        BOOST_UBLAS_INLINE
        typename std::map<I, T>::mapped_type &make_reference (std::map<I, T> &/* a */, typename std::map<I, T>::iterator it) {
            return (*it).second;
        }
#endif

    }

    // This specialization is missing in Dinkumware's STL?!
    template<class I, class T, class F>
    BOOST_UBLAS_INLINE
    void swap (std::map<I, T, F> &a1, std::map<I, T, F> &a2) {
        // Too unusual semantic.
        // BOOST_UBLAS_CHECK (&a1 != &a2, external_logic ());
        if (&a1 != &a2)
            a1.swap (a2);
    }

    // Set array
    template<class I>
    class set_array {
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef I index_type;
        typedef I value_type;
        typedef const I &const_reference;
        typedef I &reference;
        typedef const I *const_pointer;
        typedef I *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        set_array ():
            capacity_ (0), data_ (new value_type [0]), size_ (0) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            std::fill (data_, data_ + size_, value_type ());
        }
        BOOST_UBLAS_INLINE
        set_array (no_init):
            capacity_ (0), data_ (new value_type [0]), size_ (0) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
        }
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        set_array (size_type size):
            capacity_ (size), data_ (new value_type [size]), size_ (0) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            std::fill (data_, data_ + size_, value_type ());
        }
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        set_array (size_type size, no_init):
            capacity_ (size), data_ (new value_type [size]), size_ (0) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
        }
        BOOST_UBLAS_INLINE
        set_array (const set_array &a):
            capacity_ (a.size_), data_ (new value_type [a.size_]), size_ (a.size_) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            *this = a;
        }
        BOOST_UBLAS_INLINE
        ~set_array () {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            delete [] data_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size) {
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
            if (size > capacity_) {
                pointer data = new value_type [size << 1];
                // Assuming std compliant allocator as requested during review.
                // if (! data)
                //     throw std::bad_alloc ();
                // if (! data_)
                //     throw std::bad_alloc ();
                std::copy (data_, data_ + std::min (size, size_), data);
                std::fill (data + std::min (size, size_), data + size, value_type ());
                delete [] data_;
                capacity_ = size << 1;
                data_ = data;
            }
            size_ = size;
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
        }

        // Reserving
        BOOST_UBLAS_INLINE
        void reserve (size_type capacity) {
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
            if (capacity > capacity_) {
                pointer data = new value_type [capacity];
                // Assuming std compliant allocator as requested during review.
                // if (! data)
                //     throw std::bad_alloc ();
                // if (! data_)
                //     throw std::bad_alloc ();
                std::copy (data_, data_ + size_, data);
                delete [] data_;
                capacity_ = capacity;
                data_ = data;
            }
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
        }

        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator [] (index_type i) {
            pointer it = find (i);
            if (it == end ())
                it = insert (end (), i);
            BOOST_UBLAS_CHECK (it != end (), internal_logic ());
            return *it;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        set_array &operator = (const set_array &a) {
            // Too unusual semantic.
            // Thanks to Michael Stevens for spotting this.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                resize (a.size_);
                std::copy (a.data_, a.data_ + a.size_, data_);
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        set_array &assign_temporary (set_array &a) {
            swap (a);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (set_array &a) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                std::swap (capacity_, a.capacity_);
                std::swap (data_, a.data_);
                std::swap (size_, a.size_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (set_array &a1, set_array &a2) {
            a1.swap (a2);
        }
#endif

        // Element insertion and deletion
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer push_back (pointer it, const value_type &p) {
            if (size () == 0 || (*(it = end () - 1)) < p) {
                resize (size () + 1);
                *(it = end () - 1) = p;
                return it;
            }
            // Raising exceptions abstracted as requested during review.
            // throw external_logic ();
            external_logic ().raise ();
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer insert (pointer it, const value_type &p) {
            it = detail::lower_bound (begin (), end (), p, std::less<value_type> ());
            difference_type n = it - begin ();
            BOOST_UBLAS_CHECK (size () == 0 || size () == size_type (n) || *it != p, external_logic ());
            resize (size () + 1);
            it = begin () + n;
            std::copy_backward (it, end () - 1, end ());
            *it = p;
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void insert (pointer it, pointer it1, pointer it2) {
#ifdef BOOST_UBLAS_BOUNDS_CHECK
            while (it1 != it2) {
                insert (it, *it1);
                ++ it1;
            }
#else
            difference_type n = it - begin ();
            resize (size () + it2 - it1);
            it = begin () + n;
            std::copy (it1, it2, it);
            std::sort (begin (), end (), std::less<value_type> ());
#endif
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void erase (pointer it) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            std::copy (it + 1, end (), it);
            resize (size () - 1);
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void erase (pointer it1, pointer it2) {
            BOOST_UBLAS_CHECK (begin () <= it1 && it1 < it2 && it2 <= end (), bad_index ());
            std::copy (it2, end (), it1);
            resize (size () - (it2 - it1));
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void clear () {
            resize (0);
        }

        // Element lookup
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_pointer find (index_type i) const {
#ifdef BOOST_UBLAS_DEPRECATED
            std::pair<const_pointer, const_pointer> pit;
            pit = std::equal_range (begin (), end (), i, std::less<value_type> ());
            if (*pit.first == i)
                return pit.first;
            else if (*pit.second == i)
                return pit.second;
            else
                return end ();
#else
            const_pointer it (detail::lower_bound (begin (), end (), i, std::less<value_type> ()));
            if (it == end () || *it != i)
                it = end ();
            return it;
#endif
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer find (index_type i) {
#ifdef BOOST_UBLAS_DEPRECATED
            std::pair<pointer, pointer> pit;
            pit = std::equal_range (begin (), end (), i, std::less<value_type> ());
            if (*pit.first == i)
                return pit.first;
            else if (*pit.second == i)
                return pit.second;
            else
                return end ();
#else
            pointer it (detail::lower_bound (begin (), end (), i, std::less<value_type> ()));
            if (it == end () || *it != i)
                it = end ();
            return it;
#endif
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_pointer lower_bound (index_type i) const {
            return detail::lower_bound (begin (), end (), i, std::less<value_type> ());
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer lower_bound (index_type i) {
            return detail::lower_bound (begin (), end (), i, std::less<value_type> ());
        }
#ifdef BOOST_UBLAS_DEPRECATED
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_pointer upper_bound (index_type i) const {
            return detail::upper_bound (begin (), end (), i, std::less<value_type> ());
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer upper_bound (index_type i) {
            return detail::upper_bound (begin (), end (), i, std::less<value_type> ());
        }
#endif

        // Iterators simply are pointers.

        typedef const_pointer const_iterator;

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return data_ + size_;
        }

        typedef pointer iterator;

        BOOST_UBLAS_INLINE
        iterator begin () {
            return data_;
        }
        BOOST_UBLAS_INLINE
        iterator end () {
            return data_ + size_;
        }

        // Reverse iterators

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef std::reverse_iterator<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator rbegin () const {
            return const_reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator rend () const {
            return const_reverse_iterator (begin ());
        }

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef std::reverse_iterator<iterator, value_type, reference> reverse_iterator;
#else
        typedef std::reverse_iterator<iterator> reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        reverse_iterator rbegin () {
            return reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator rend () {
            return reverse_iterator (begin ());
        }

    private:
        size_type capacity_;
        pointer data_;
        size_type size_;
    };

    // This specialization is missing in Dinkumware's STL?!
    template<class I, class F>
    BOOST_UBLAS_INLINE
    void swap (std::set<I, F> &a1, std::set<I, F> &a2) {
        // Too unusual semantic.
        // BOOST_UBLAS_CHECK (&a1 != &a2, external_logic ());
        if (&a1 != &a2)
            a1.swap (a2);
    }

}}}

#endif


