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

#ifndef BOOST_UBLAS_STORAGE_H
#define BOOST_UBLAS_STORAGE_H

#include <algorithm>
#include <valarray>
#include <vector>

// #define BOOST_UBLAS_SIMPLE_ARRAY_ADAPTOR
#ifndef BOOST_UBLAS_SIMPLE_ARRAY_ADAPTOR
#include <boost/shared_array.hpp>
#endif

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/iterator.hpp>
#include <boost/numeric/ublas/traits.hpp>

namespace boost { namespace numeric { namespace ublas {

#ifndef BOOST_UBLAS_USE_FAST_SAME
// FIXME: for performance reasons we better use macros
//    template<class T>
//    BOOST_UBLAS_INLINE
//    const T &same_impl (const T &size1, const T &size2) {
//        BOOST_UBLAS_CHECK (size1 == size2, bad_argument ());
//        return std::min (size1, size2);
//    }
// #define BOOST_UBLAS_SAME(size1, size2) same_impl ((size1), (size2))
    template<class T>
    BOOST_UBLAS_INLINE
    // Kresimir Fresl and Dan Muller reported problems with COMO.
    // We better change the signature instead of libcomo ;-)
    // const T &same_impl_ex (const T &size1, const T &size2, const char *file, int line) {
    T same_impl_ex (const T &size1, const T &size2, const char *file, int line) {
        BOOST_UBLAS_CHECK_EX (size1 == size2, file, line, bad_argument ());
        return std::min (size1, size2);
    }
#define BOOST_UBLAS_SAME(size1, size2) same_impl_ex ((size1), (size2), __FILE__, __LINE__)
#else
// FIXME: for performance reasons we better use macros
//    template<class T>
//    BOOST_UBLAS_INLINE
//    const T &same_impl (const T &size1, const T &size2) {
//        return size1;
//    }
// #define BOOST_UBLAS_SAME(size1, size2) same_impl ((size1), (size2))
#define BOOST_UBLAS_SAME(size1, size2) (size1)
#endif

    struct no_init {};

    // Unbounded array
    template<class T>
    class unbounded_array {
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef const T *const_pointer;
        typedef T *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        unbounded_array ():
            size_ (0), data_ (new value_type [0]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            std::fill (data_, data_ + size_, value_type ());
        }
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        unbounded_array (no_init):
            size_ (0), data_ (new value_type [0]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
        }
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        unbounded_array (size_type size):
            size_ (size), data_ (new value_type [size]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            std::fill (data_, data_ + size_, value_type ());
        }
        BOOST_UBLAS_INLINE
        unbounded_array (size_type size, no_init):
            size_ (size), data_ (new value_type [size]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
        }
        BOOST_UBLAS_INLINE
        unbounded_array (const unbounded_array &a):
            size_ (a.size_), data_ (new value_type [a.size_]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            *this = a;
        }
        BOOST_UBLAS_INLINE
        ~unbounded_array () {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            delete [] data_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            if (size != size_) {
                pointer data = new value_type [size];
                // Assuming std compliant allocator as requested during review.
                // if (! data)
                //     throw std::bad_alloc ();
                // if (! data_)
                //     throw std::bad_alloc ();
                if (preserve) {
                    std::copy (data_, data_ + std::min (size, size_), data);
                    std::fill (data + std::min (size, size_), data + size, value_type ());
                }
                delete [] data_;
                size_ = size;
                data_ = data;
            }
        }

        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }
        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }

        // Assignment
        BOOST_UBLAS_INLINE
        unbounded_array &operator = (const unbounded_array &a) {
            // Too unusual semantic.
            // Thanks to Michael Stevens for spotting this.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == a.size_, bad_size ());
                resize (a.size_, false);
                std::copy (a.data_, a.data_ + a.size_, data_);
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        unbounded_array &assign_temporary (unbounded_array &a) {
            swap (a);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (unbounded_array &a) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == a.size_, bad_size ());
                std::swap (size_, a.size_);
                std::swap (data_, a.data_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (unbounded_array &a1, unbounded_array &a2) {
            a1.swap (a2);
        }
#endif

        // Element insertion and deletion
        BOOST_UBLAS_INLINE
        pointer insert (pointer it, const value_type &t) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            BOOST_UBLAS_CHECK (*it == value_type (), external_logic ());
            *it = t;
            return it;
        }
        BOOST_UBLAS_INLINE
        void insert (pointer it, pointer it1, pointer it2) {
            while (it1 != it2) {
                BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
                BOOST_UBLAS_CHECK (*it == value_type (), external_logic ());
                *it = *it1;
                ++ it, ++ it1;
            }
        }
        BOOST_UBLAS_INLINE
        void erase (pointer it) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            *it = value_type ();
        }
        BOOST_UBLAS_INLINE
        void erase (pointer it1, pointer it2) {
            while (it1 != it2) {
                BOOST_UBLAS_CHECK (begin () <= it1 && it1 < end (), bad_index ());
                *it1 = value_type ();
                ++ it1;
            }
        }
        BOOST_UBLAS_INLINE
        void clear () {
            erase (begin (), end ());
        }

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
        size_type size_;
        pointer data_;
    };

    // Bounded array
    template<class T, std::size_t N>
    class bounded_array {
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef const T *const_pointer;
        typedef T *pointer;

        // Construction and destruction
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        bounded_array ():
            // Kresimir Fresl suggested to change the default back to the template argument.
            // size_ (0) /* , data_ () */ {
            size_ (N) /* , data_ () */ {
            std::fill (data_, data_ + size_, value_type ());
        }
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        bounded_array (no_init):
            // Kresimir Fresl suggested to change the default back to the template argument.
            // size_ (0) /* , data_ () */ {
            size_ (N) /* , data_ () */ {}
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        bounded_array (size_type size):
            size_ (size) /* , data_ () */ {
            if (size_ > N)
                // Raising exceptions abstracted as requested during review.
                // throw std::bad_alloc ();
                bad_size ().raise ();
            std::fill (data_, data_ + size_, value_type ());
        }
        BOOST_UBLAS_INLINE
        bounded_array (size_type size, no_init):
            size_ (size) /* , data_ () */ {
            if (size_ > N)
                // Raising exceptions abstracted as requested during review.
                // throw std::bad_alloc ();
                bad_size ().raise ();
        }
        BOOST_UBLAS_INLINE
        bounded_array (const bounded_array &a):
            size_ (a.size_) /* , data_ () */ {
            if (size_ > N)
                // Raising exceptions abstracted as requested during review.
                // throw std::bad_alloc ();
                bad_size ().raise ();
            *this = a;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            if (size > N)
                // Raising exceptions abstracted as requested during review.
                // throw std::bad_alloc ();
                bad_size ().raise ();
            if (preserve)
                std::fill (data_ + std::min (size, size_), data_ + size, value_type ());
            size_ = size;
        }

        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }
        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }

        // Assignment
        BOOST_UBLAS_INLINE
        bounded_array &operator = (const bounded_array &a) {
            // Too unusual semantic.
            // Thanks to Michael Stevens for spotting this.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == a.size_, bad_size ());
                resize (a.size_, false);
                std::copy (a.data_, a.data_ + a.size_, data_);
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        bounded_array &assign_temporary (bounded_array &a) { 
            *this = a;
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (bounded_array &a) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == a.size_, bad_size ());
                std::swap (size_, a.size_);
                std::swap_ranges (data_, data_ + std::max (size_, a.size_), a.data_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (bounded_array &a1, bounded_array &a2) {
            a1.swap (a2);
        }
#endif

        // Element insertion and deletion
        BOOST_UBLAS_INLINE
        pointer insert (pointer it, const value_type &t) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            BOOST_UBLAS_CHECK (*it == value_type (), external_logic ());
            *it = t;
            return it;
        }
        BOOST_UBLAS_INLINE
        void insert (pointer it, pointer it1, pointer it2) {
            while (it1 != it2) {
                BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
                BOOST_UBLAS_CHECK (*it == value_type (), external_logic ());
                *it = *it1;
                ++ it, ++ it1;
            }
        }
        BOOST_UBLAS_INLINE
        void erase (pointer it) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            *it = value_type ();
        }
        BOOST_UBLAS_INLINE
        void erase (pointer it1, pointer it2) {
            while (it1 != it2) {
                BOOST_UBLAS_CHECK (begin () <= it1 && it1 < end (), bad_index ());
                *it1 = value_type ();
                ++ it1;
            }
        }
        BOOST_UBLAS_INLINE
        void clear () {
            erase (begin (), end ());
        }

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
        size_type size_;
        BOOST_UBLAS_ALIGN_16 value_type data_ [N];
    };

#ifdef BOOST_UBLAS_SIMPLE_ARRAY_ADAPTOR

    // Array adaptor
    template<class T>
    class array_adaptor {
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef const T *const_pointer;
        typedef T *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        array_adaptor ():
            size_ (0), own_ (true), data_ (new value_type [0]) {
            std::fill (data_, data_ + size_, value_type ());
        }
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        array_adaptor (no_init):
            size_ (0), own_ (true), data_ (new value_type [0]) {}
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        array_adaptor (size_type size):
            size_ (size), own_ (true), data_ (new value_type [size]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            std::fill (data_, data_ + size_, value_type ());
        }
        BOOST_UBLAS_INLINE
        array_adaptor (size_type size, no_init):
            size_ (size), own_ (true), data_ (new value_type [size]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
        }
        BOOST_UBLAS_INLINE
        array_adaptor (size_type size, pointer data):
            size_ (size), own_ (false), data_ (data) {}
#ifdef BOOST_UBLAS_DEEP_COPY
        BOOST_UBLAS_INLINE
        array_adaptor (const array_adaptor &a):
            size_ (a.size_), own_ (true), data_ (new value_type [a.size_]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            *this = a;
        }
#else
        BOOST_UBLAS_INLINE
        array_adaptor (const array_adaptor &a):
            size_ (a.size_), own_ (a.own_), data_ (a.data_) {
            if (own_)
                // Raising exceptions abstracted as requested during review.
                // throw std::bad_alloc ();
                external_logic ().raise ();
        }
#endif
        BOOST_UBLAS_INLINE
        ~array_adaptor () {
            if (own_) {
                // Assuming std compliant allocator as requested during review.
                // if (! data_)
                //     throw std::bad_alloc ();
                delete [] data_;
            }
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            if (size != size_) {
                pointer data = new value_type [size];
                // Assuming std compliant allocator as requested during review.
                // if (! data)
                //     throw std::bad_alloc ();
                // if (! data_)
                //     throw std::bad_alloc ();
                if (preserve) {
                    std::copy (data_, data_ + std::min (size, size_), data);
                    std::fill (data + std::min (size, size_), data + size, value_type ());
                }
                if (own_)
                    delete [] data_;
                size_ = size;
                own_ = true;
                data_ = data;
            }
        }
        BOOST_UBLAS_INLINE
        void resize (size_type size, pointer data, bool preserve = true) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_)
            //     throw std::bad_alloc ();
            if (preserve) {
                std::copy (data_, data_ + std::min (size, size_), data);
                std::fill (data + std::min (size, size_), data + size, value_type ());
            }
            if (own_)
                delete [] data_;
            size_ = size;
            own_ = false;
            data_ = data;
        }

        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }
        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }

        // Assignment
        BOOST_UBLAS_INLINE
        array_adaptor &operator = (const array_adaptor &a) {
            // Too unusual semantic.
            // Thanks to Michael Stevens for spotting this.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == a.size_, bad_size ());
                resize (a.size_, false);
                std::copy (a.data_, a.data_ + a.size_, data_);
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        array_adaptor &assign_temporary (array_adaptor &a) {
            if (own_ && a.own_)
                swap (a);
            else
                *this = a;
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (array_adaptor &a) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == a.size_, bad_size ());
                std::swap (size_, a.size_);
                std::swap (own_, a.own_);
                std::swap (data_, a.data_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (array_adaptor &a1, array_adaptor &a2) {
            a1.swap (a2);
        }
#endif

        // Element insertion and deletion
        BOOST_UBLAS_INLINE
        pointer insert (pointer it, const value_type &t) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            BOOST_UBLAS_CHECK (*it == value_type (), external_logic ());
            *it = t;
            return it;
        }
        BOOST_UBLAS_INLINE
        void insert (pointer it, pointer it1, pointer it2) {
            while (it1 != it2) {
                BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
                BOOST_UBLAS_CHECK (*it == value_type (), external_logic ());
                *it = *it1;
                ++ it, ++ it1;
            }
        }
        BOOST_UBLAS_INLINE
        void erase (pointer it) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            *it = value_type ();
        }
        BOOST_UBLAS_INLINE
        void erase (pointer it1, pointer it2) {
            while (it1 != it2) {
                BOOST_UBLAS_CHECK (begin () <= it1 && it1 < end (), bad_index ());
                *it1 = value_type ();
                ++ it1;
            }
        }
        BOOST_UBLAS_INLINE
        void clear () {
            erase (begin (), end ());
        }

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
        size_type size_;
        bool own_;
        pointer data_;
    };

#else

    template<class T>
    struct leaker {
        typedef void result_type;
        typedef T *argument_type;

        BOOST_UBLAS_INLINE
        result_type operator () (argument_type x) {}
    };

    // Array adaptor
    template<class T>
    class array_adaptor {
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef const T *const_pointer;
        typedef T *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        array_adaptor ():
            size_ (0), own_ (true), data_ (new value_type [0]) {
            std::fill (data_.get (), data_.get () + size_, value_type ());
        }
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        array_adaptor (no_init):
            size_ (0), own_ (true), data_ (new value_type [0]) {}
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        array_adaptor (size_type size):
            size_ (size), own_ (true), data_ (new value_type [size]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_.get ())
            //     throw std::bad_alloc ();
            std::fill (data_.get (), data_.get () + size_, value_type ());
        }
        BOOST_UBLAS_INLINE
        array_adaptor (size_type size, no_init):
            size_ (size), own_ (true), data_ (new value_type [size]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_.get ())
            //     throw std::bad_alloc ();
        }
        BOOST_UBLAS_INLINE
        array_adaptor (size_type size, pointer data):
            size_ (size), own_ (false), data_ (data, leaker<value_type> ()) {}
#ifdef BOOST_UBLAS_DEEP_COPY
        BOOST_UBLAS_INLINE
        array_adaptor (const array_adaptor &a):
            size_ (a.size_), own_ (true), data_ (new value_type [a.size_]) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_.get ())
            //     throw std::bad_alloc ();
            *this = a;
        }
#else
        BOOST_UBLAS_INLINE
        array_adaptor (const array_adaptor &a):
            size_ (a.size_), own_ (a.own_), data_ (a.data_) {}
#endif
        BOOST_UBLAS_INLINE
        ~array_adaptor () {
            // Assuming std compliant allocator as requested during review.
            // if (! data_.get ())
            //     throw std::bad_alloc ();
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            if (size != size_) {
                shared_array<value_type> data (new value_type [size]);
                // Assuming std compliant allocator as requested during review.
                // if (! data.get ())
                //     throw std::bad_alloc ();
                // if (! data_.get ())
                //     throw std::bad_alloc ();
                if (preserve) {
                    std::copy (data_.get (), data_.get () + std::min (size, size_), data.get ());
                    std::fill (data.get () + std::min (size, size_), data.get () + size, value_type ());
                }
                size_ = size;
                data_ = data;
            }
        }
        BOOST_UBLAS_INLINE
        void resize (size_type size, pointer data, bool preserve = true) {
            // Assuming std compliant allocator as requested during review.
            // if (! data_.get ())
            //     throw std::bad_alloc ();
            if (preserve) {
                std::copy (data_.get (), data_.get () + std::min (size, size_), data);
                std::fill (data + std::min (size, size_), data + size, value_type ());
            }
            size_ = size;
            data_ = data;
        }

        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }
        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }

        // Assignment
        BOOST_UBLAS_INLINE
        array_adaptor &operator = (const array_adaptor &a) {
            // Too unusual semantic.
            // Thanks to Michael Stevens for spotting this.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == a.size_, bad_size ());
                resize (a.size_, false);
                std::copy (a.data_.get (), a.data_.get () + a.size_, data_.get ());
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        array_adaptor &assign_temporary (array_adaptor &a) {
            if (own_ && a.own_)
                swap (a);
            else
                *this = a;
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (array_adaptor &a) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &a, external_logic ());
            if (this != &a) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == a.size_, bad_size ());
                std::swap (size_, a.size_);
                std::swap (own_, a.own_);
                std::swap (data_, a.data_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (array_adaptor &a1, array_adaptor &a2) {
            a1.swap (a2);
        }
#endif

        // Element insertion and deletion
        BOOST_UBLAS_INLINE
        pointer insert (pointer it, const value_type &t) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            BOOST_UBLAS_CHECK (*it == value_type (), external_logic ());
            *it = t;
            return it;
        }
        BOOST_UBLAS_INLINE
        void insert (pointer it, pointer it1, pointer it2) {
            while (it1 != it2) {
                BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
                BOOST_UBLAS_CHECK (*it == value_type (), external_logic ());
                *it = *it1;
                ++ it, ++ it1;
            }
        }
        BOOST_UBLAS_INLINE
        void erase (pointer it) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            *it = value_type ();
        }
        BOOST_UBLAS_INLINE
        void erase (pointer it1, pointer it2) {
            while (it1 != it2) {
                BOOST_UBLAS_CHECK (begin () <= it1 && it1 < end (), bad_index ());
                *it1 = value_type ();
                ++ it1;
            }
        }
        BOOST_UBLAS_INLINE
        void clear () {
            erase (begin (), end ());
        }

        // Iterators simply are pointers.

        typedef const_pointer const_iterator;

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return data_.get ();
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return data_.get () + size_;
        }

        typedef pointer iterator;

        BOOST_UBLAS_INLINE
        iterator begin () {
            return data_.get ();
        }
        BOOST_UBLAS_INLINE
        iterator end () {
            return data_.get () + size_;
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
        size_type size_;
        bool own_;
        shared_array<value_type> data_;
    };

#endif

    namespace detail {
        using namespace boost::numeric::ublas;

        // Some helpers for unbounded_array

        template<class T>
        BOOST_UBLAS_INLINE
        void resize (unbounded_array<T> &a, typename unbounded_array<T>::size_type size, bool preserve) {
            a.resize (size, preserve);
        }

        // Some helpers for bounded_array

        template<class T, std::size_t N>
        BOOST_UBLAS_INLINE
        void resize (bounded_array<T, N> &a, typename bounded_array<T, N>::size_type size, bool preserve) {
            a.resize (size, preserve);
        }

        // Some helpers for array_adaptor

        template<class T>
        BOOST_UBLAS_INLINE
        void resize (array_adaptor<T> &a, typename array_adaptor<T>::size_type size, bool preserve) {
            a.resize (size, preserve);
        }

        // Some helpers for std::vector

        template<class T>
        BOOST_UBLAS_INLINE
        void resize (std::vector<T> &a, typename std::vector<T>::size_type size, bool /* preserve */) {
            a.resize (size);
        }

    }

    // Range class
    class range {
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef difference_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const difference_type *const_pointer;
        typedef difference_type *pointer;
        typedef size_type const_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        range ():
            start_ (), size_ () {}
        BOOST_UBLAS_INLINE
        range (size_type start, size_type stop):
            start_ (start), size_ (stop - start) {
            BOOST_UBLAS_CHECK (start <= stop, bad_size ());
        }

        BOOST_UBLAS_INLINE
        size_type start () const {
            return start_;
        }
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return start_ + i;
        }

        // Composition
        BOOST_UBLAS_INLINE
        range compose (const range &r) const {
            BOOST_UBLAS_CHECK (r.start_ + r.size_ <= size_, bad_size ());
            return range (start_ + r.start_, start_ + r.start_ + r.size_);
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const range &r) const {
            return start_ == r.start_ && size_ == r.size_;
        }
        BOOST_UBLAS_INLINE
        bool operator != (const range &r) const {
            return ! (*this == r);
        }

        // Iterator simply is a index.

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator<range, std::random_access_iterator_tag> const_iterator;
#else
        class const_iterator:
            public container_const_reference<range>,
            public random_access_iterator_base<std::random_access_iterator_tag,
                                               const_iterator, value_type> {
        public:
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef range::difference_type difference_type;
            typedef range::value_type value_type;
            typedef range::const_reference reference;
            typedef range::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<range> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const range &r, const const_iterator_type &it):
                container_const_reference<range> (r), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator &operator ++ () {
                ++ it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -- () {
                -- it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator += (difference_type n) {
                it_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -= (difference_type n) {
                it_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator &it) const {
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK ((*this) ().start () <= it_, bad_index ());
                BOOST_UBLAS_CHECK (it_ < (*this) ().start () + (*this) ().size (), bad_index ());
                return it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_ - (*this) ().start ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                // Comeau recommends...
                this->assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) () == it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) () == it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator_type it_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return const_iterator (*this, start_);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return const_iterator (*this, start_ + size_);
        }

        // Reverse iterator

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

        BOOST_UBLAS_INLINE
        range preprocess (size_type size) const {
            if (*this != all ())
                return *this;
            return range (0, size);
        }
        static
        BOOST_UBLAS_INLINE
        range all () {
            return range (0, size_type (-1));
        }

    private:
        size_type start_;
        size_type size_;
    };

    // Slice class
    class slice {
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef difference_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const difference_type *const_pointer;
        typedef difference_type *pointer;
        typedef difference_type const_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        slice ():
            start_ (), stride_ (), size_ () {}
        BOOST_UBLAS_INLINE
        slice (size_type start, difference_type stride, size_type size):
            start_ (start), stride_ (stride), size_ (size) {}

        BOOST_UBLAS_INLINE
        size_type start () const {
            return start_;
        }
        BOOST_UBLAS_INLINE
        difference_type stride () const {
            return stride_;
        }
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return start_ + i * stride_;
        }

        // Composition
        BOOST_UBLAS_INLINE
        slice compose (const range &r) const {
            BOOST_UBLAS_CHECK (r.start () + r.size () <= size_, bad_size ());
            return slice (start_ + stride_ * r.start (), stride_, r.size ());
        }
        BOOST_UBLAS_INLINE
        slice compose (const slice &s) const {
            BOOST_UBLAS_CHECK (s.start_ + s.stride_ * (s.size_ - (s.size_ > 0)) <= size_, bad_size ());
            return slice (start_ + stride_ * s.start_, stride_ * s.stride_, s.size_);
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const slice &s) const {
            return start_ == s.start_ && stride_ == s.stride_ && size_ == s.size_; 
        }
        BOOST_UBLAS_INLINE
        bool operator != (const slice &s) const {
            return ! (*this == s);
        }

        // Iterator simply is a index.

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator<slice, std::random_access_iterator_tag> const_iterator;
#else
        class const_iterator:
            public container_const_reference<slice>,
            public random_access_iterator_base<std::random_access_iterator_tag,
                                               const_iterator, value_type> {
        public:
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef slice::difference_type difference_type;
            typedef slice::value_type value_type;
            typedef slice::const_reference reference;
            typedef slice::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<slice> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const slice &s, const const_iterator_type &it):
                container_const_reference<slice> (s), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator &operator ++ () {
                it_ += (*this) ().stride ();
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -- () {
                it_ -= (*this) ().stride ();
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator += (difference_type n) {
                it_ += n * (*this) ().stride ();
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -= (difference_type n) {
                it_ -= n * (*this) ().stride ();
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator &it) const {
                // Stationary slice must remain at start OR stride must be non zero
                // Thanks to Michael Stevens for this extension.
                BOOST_UBLAS_CHECK ((it_ == difference_type ((*this) ().start ()) && it.it_ == difference_type ((*this) ().start ())) ||
                                   (*this) ().stride () != 0, divide_by_zero ());
                return (*this) ().stride () != 0 ? (it_ - it.it_) / (*this) ().stride () : 0;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                // Stationary slice must remain at start OR index must remain within size
                // Thanks to Michael Stevens for this extension.
                BOOST_UBLAS_CHECK (((*this) ().stride () == 0 && it_ == difference_type ((*this) ().start ())) ||
                                   index () < (*this) ().size (), bad_index ());
                return it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                // Stationary slice must remain at start OR stride must be non zero
                // Thanks to Michael Stevens for this extension.
                BOOST_UBLAS_CHECK (it_ == difference_type ((*this) ().start ()) ||
                                   (*this) ().stride () != 0, divide_by_zero ());
                return (*this) ().stride () != 0 ? (it_ - difference_type ((*this) ().start ())) / (*this) ().stride () : 0;
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                // Comeau recommends...
                this->assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) () == it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) () == it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator_type it_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return const_iterator (*this, start_);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return const_iterator (*this, start_ + stride_ * size_);
        }

        // Reverse iterator

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

        BOOST_UBLAS_INLINE
        slice preprocess (size_type size) const {
            if (*this != all ())
                return *this;
            return slice (0, 1, size);
        }
        static
        BOOST_UBLAS_INLINE
        slice all () {
            return slice (0, 1, size_type (-1));
        }

    private:
        size_type start_;
        difference_type stride_;
        size_type size_;
    };

    // Indirect array class
    template<class A>
    class indirect_array {
    public:
        typedef A array_type;
        typedef const A const_array_type;
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef typename A::value_type value_type;
        typedef typename A::const_reference const_reference;
        typedef typename A::reference reference;
        typedef typename A::const_pointer const_pointer;
        typedef typename A::pointer pointer;
        typedef difference_type const_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        indirect_array ():
            size_ (), data_ () {}
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        indirect_array (size_type size):
            size_ (size), data_ (size) {}
        BOOST_UBLAS_INLINE
        indirect_array (size_type size, const array_type &data):
            size_ (size), data_ (data) {}
        BOOST_UBLAS_INLINE
        indirect_array (pointer start, pointer stop):
            size_ (stop - start), data_ (stop - start) {
            std::copy (start, stop, data_.begin ());
        }

        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }
        BOOST_UBLAS_INLINE
        const_array_type data () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        array_type data () {
            return data_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }

        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            return (*this) (i);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) {
            return (*this) (i);
        }

        // Composition
        BOOST_UBLAS_INLINE
        indirect_array compose (const range &r) const {
            BOOST_UBLAS_CHECK (r.start () + r.size () <= size_, bad_size ());
            array_type data (r.size ());
            for (size_type i = 0; i < r.size (); ++ i)
                data [i] = data_ [r.start () + i];
            return indirect_array (r.size (), data);
        }
        BOOST_UBLAS_INLINE
        indirect_array compose (const slice &s) const {
            BOOST_UBLAS_CHECK (s.start () + s.stride () * (s.size () - (s.size () > 0)) <= size (), bad_size ());
            array_type data (s.size ());
            for (size_type i = 0; i < s.size (); ++ i)
                data [i] = data_ [s.start () + s.stride () * i];
            return indirect_array (s.size (), data);
        }
        BOOST_UBLAS_INLINE
        indirect_array compose (const indirect_array &ia) const {
            array_type data (ia.size_);
            for (size_type i = 0; i < ia.size_; ++ i) {
                BOOST_UBLAS_CHECK (ia.data_ [i] <= size_, bad_size ());
                data [i] = data_ [ia.data_ [i]];
            }
            return indirect_array (ia.size_, data);
        }

        // Comparison
        template<class OA>
        BOOST_UBLAS_INLINE
        bool operator == (const indirect_array<OA> &ia) const {
            if (size_ != ia.size_)
                return false;
            for (size_type i = 0; i < BOOST_UBLAS_SAME (size_, ia.size_); ++ i)
                if (data_ [i] != ia.data_ [i])
                    return false;
            return true;
        }
        template<class OA>
        BOOST_UBLAS_INLINE
        bool operator != (const indirect_array<OA> &ia) const {
            return ! (*this == ia);
        }

        // Iterator simply is a index.

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator<indirect_array, std::random_access_iterator_tag> const_iterator;
#else
        class const_iterator:
            public container_const_reference<indirect_array>,
            public random_access_iterator_base<std::random_access_iterator_tag,
                                               const_iterator, value_type> {
        public:
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename indirect_array::difference_type difference_type;
            typedef typename indirect_array::value_type value_type;
            typedef typename indirect_array::const_reference reference;
            typedef typename indirect_array::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<indirect_array> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const indirect_array &ia, const const_iterator_type &it):
                container_const_reference<indirect_array> (ia), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator &operator ++ () {
                ++ it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -- () {
                -- it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator += (difference_type n) {
                it_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -= (difference_type n) {
                it_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator &it) const {
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return (*this) () (it_);
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_;
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                // Comeau recommends...
                this->assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) () == it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) () == it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator_type it_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return const_iterator (*this, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return const_iterator (*this, size_);
        }

        // Reverse iterator

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

        BOOST_UBLAS_INLINE
        indirect_array preprocess (size_type size) const {
            if (this != &all_)
                return *this;
            indirect_array ia (size);
            for (size_type i = 0; i < size; ++ i)
               ia (i) = i;
            return ia;
        }
        static
        BOOST_UBLAS_INLINE
        const indirect_array &all () {
            return all_;
        }

    private:
        size_type size_;
        array_type data_;
        static indirect_array all_;
    };

    template<class A>
    indirect_array<A> indirect_array<A>::all_;

    // Gunter Winkler contributed the classes index_pair, index_pair_array,
    // index_triple and index_triple_array to enable inplace sort of parallel arrays.

    template <class V>
    class index_pair :
        private boost::noncopyable,
        public container_reference<V> {
    public:
        typedef index_pair<V> self_type;
        typedef typename V::size_type size_type;

        BOOST_UBLAS_INLINE
        index_pair(V& v, size_type i) :
            boost::noncopyable (), container_reference<V>(v), i_(i),
            v1_(v.data1_[i]), v2_(v.data2_[i]), dirty_(false) {}
        BOOST_UBLAS_INLINE
        index_pair(const self_type& rhs) :
            boost::noncopyable (), container_reference<V>(rhs()), i_(rhs.i_),
            v1_(rhs.v1_), v2_(rhs.v2_), dirty_(false) {}
        BOOST_UBLAS_INLINE
        ~index_pair() {
            if (dirty_) {
                (*this)().data1_[i_] = v1_;
                (*this)().data2_[i_] = v2_;
            }
        }

        BOOST_UBLAS_INLINE
        self_type& operator=(const self_type& rhs) {
            v1_ = rhs.v1_;
            v2_ = rhs.v2_;
            dirty_ = true;
            return *this;
        }

        BOOST_UBLAS_INLINE
        void swap(self_type rhs) {
            self_type tmp(rhs);
            rhs = *this;
            *this = tmp;
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap(self_type lhs, self_type rhs) {
            lhs.swap(rhs);
        }
#endif

        BOOST_UBLAS_INLINE
        bool compare(const self_type& rhs) const {
            return (v1_ < rhs.v1_);
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend bool operator < (const self_type& lhs, const self_type& rhs) {
            return lhs.compare(rhs);
        }
#endif

    private:
        size_type i_;
        typename V::value1_type v1_;
        typename V::value2_type v2_;
        bool dirty_;
    };

#ifdef BOOST_UBLAS_NO_MEMBER_FRIENDS
    template<class V>
    BOOST_UBLAS_INLINE
    void swap(index_pair<V> lhs, index_pair<V> rhs) {
        lhs.swap(rhs);
    }
#endif

#if defined (BOOST_UBLAS_NO_MEMBER_FRIENDS) || defined (BOOST_MSVC)
    template<class V>
    BOOST_UBLAS_INLINE
    bool operator < (const index_pair<V>& lhs, const index_pair<V>& rhs) {
        return lhs.compare(rhs);
    }
#endif

    template <class V1, class V2>
    class index_pair_array:
        private boost::noncopyable {
    public:
        typedef index_pair_array<V1, V2> self_type;
        typedef typename V1::value_type value1_type;
        typedef typename V2::value_type value2_type;

        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef index_pair<self_type> value_type;
        typedef value_type reference;
        // typedef const value_type& const_reference;
        typedef const value_type const_reference;

        BOOST_UBLAS_INLINE
        index_pair_array(size_type size, V1& data1, V2& data2) :
              size_(size),data1_(data1),data2_(data2) {}

        BOOST_UBLAS_INLINE
        size_type size() const {
            return size_;
        }

        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            return value_type((*this), i);
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) {
            return value_type((*this), i);
        }

        typedef indexed_iterator<self_type, std::random_access_iterator_tag> iterator;
        typedef indexed_const_iterator<self_type, std::random_access_iterator_tag> const_iterator;

        BOOST_UBLAS_INLINE
        iterator begin() {
            return iterator( (*this), 0);
        }
        BOOST_UBLAS_INLINE
        iterator end() {
            return iterator( (*this), size());
        }

        BOOST_UBLAS_INLINE
        const_iterator begin() const {
            return const_iterator( (*this), 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end() const {
            return const_iterator( (*this), size());
        }

        // unnecessary function:
        BOOST_UBLAS_INLINE
        bool compare(size_type i1, size_type i2) const {
            return data1_[i1] < data1_[i2];
        }

        // gives a large speedup
        BOOST_UBLAS_INLINE
        friend void iter_swap(const iterator& lhs, const iterator& rhs) {
            const size_type i1 = lhs.index();
            const size_type i2 = rhs.index();
            std::swap(lhs().data1_[i1], rhs().data1_[i2]);
            std::swap(lhs().data2_[i1], rhs().data2_[i2]);
        }

    private:
        size_type size_;
        V1& data1_;
        V2& data2_;

        // friend class value_type;
        friend class index_pair<self_type>;
    };

    template <class M>
    class index_triple :
        private boost::noncopyable,
        public container_reference<M> {
    public:
        typedef index_triple<M> self_type;
        typedef typename M::size_type size_type;

        BOOST_UBLAS_INLINE
        index_triple(M& m, size_type i) :
            boost::noncopyable (), container_reference<M>(m), i_(i),
            v1_(m.data1_[i]), v2_(m.data2_[i]), v3_(m.data3_[i]), dirty_(false) {}
        BOOST_UBLAS_INLINE
        index_triple(const self_type& rhs) :
            boost::noncopyable (), container_reference<M>(rhs()), i_(rhs.i_),
            v1_(rhs.v1_), v2_(rhs.v2_), v3_(rhs.v3_), dirty_(false) {}
        BOOST_UBLAS_INLINE
        ~index_triple() {
            if (dirty_) {
                (*this)().data1_[i_] = v1_;
                (*this)().data2_[i_] = v2_;
                (*this)().data3_[i_] = v3_;
            }
        }

        BOOST_UBLAS_INLINE
        self_type& operator=(const self_type& rhs) {
            v1_ = rhs.v1_;
            v2_ = rhs.v2_;
            v3_ = rhs.v3_;
            dirty_ = true;
            return *this;
        }

        BOOST_UBLAS_INLINE
        void swap(self_type rhs) {
            self_type tmp(rhs);
            rhs = *this;
            *this = tmp;
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap(self_type lhs, self_type rhs) {
            lhs.swap(rhs);
        }
#endif

        BOOST_UBLAS_INLINE
        bool compare(const self_type& rhs) const {
            return ((v1_ < rhs.v1_) ||
                    (v1_ == rhs.v1_ && v2_ < rhs.v2_));
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend bool operator < (const self_type& lhs, const self_type& rhs) {
            return lhs.compare(rhs);
        }
#endif

    private:
        size_type i_;
        typename M::value1_type v1_;
        typename M::value2_type v2_;
        typename M::value3_type v3_;
        bool dirty_;
    };

#ifdef BOOST_UBLAS_NO_MEMBER_FRIENDS
    template<class M>
    BOOST_UBLAS_INLINE
    void swap(index_triple<M> lhs, index_triple<M> rhs) {
        lhs.swap(rhs);
    }
#endif

#if defined (BOOST_UBLAS_NO_MEMBER_FRIENDS) || defined (BOOST_MSVC)
    template<class M>
    BOOST_UBLAS_INLINE
    bool operator < (const index_triple<M>& lhs, const index_triple<M>& rhs) {
        return lhs.compare(rhs);
    }
#endif

    template <class V1, class V2, class V3>
    class index_triple_array:
        private boost::noncopyable {
    public:
        typedef index_triple_array<V1, V2, V3> self_type;
        typedef typename V1::value_type value1_type;
        typedef typename V2::value_type value2_type;
        typedef typename V3::value_type value3_type;

        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef index_triple<self_type> value_type;
        typedef value_type reference;
        // typedef const value_type& const_reference;
        typedef const value_type const_reference;

        BOOST_UBLAS_INLINE
        index_triple_array(size_type size, V1& data1, V2& data2, V3& data3) :
              size_(size),data1_(data1),data2_(data2),data3_(data3) {}

        BOOST_UBLAS_INLINE
        size_type size() const {
            return size_;
        }

        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            return value_type((*this), i);
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) {
            return value_type((*this), i);
        }

        typedef indexed_iterator<self_type, std::random_access_iterator_tag> iterator;
        typedef indexed_const_iterator<self_type, std::random_access_iterator_tag> const_iterator;

        BOOST_UBLAS_INLINE
        iterator begin() {
            return iterator( (*this), 0);
        }
        BOOST_UBLAS_INLINE
        iterator end() {
            return iterator( (*this), size());
        }

        BOOST_UBLAS_INLINE
        const_iterator begin() const {
            return const_iterator( (*this), 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end() const {
            return const_iterator( (*this), size());
        }

        // unnecessary function:
        BOOST_UBLAS_INLINE
        bool compare(size_type i1, size_type i2) const {
            return ((data1_[i1] < data1_[i2]) ||
                    (data1_[i1] == data1_[i2] && data2_[i1] < data2_[i2]));
        }

        // gives a large speedup
        BOOST_UBLAS_INLINE
        friend void iter_swap(const iterator& lhs, const iterator& rhs) {
            const size_type i1 = lhs.index();
            const size_type i2 = rhs.index();
            std::swap(lhs().data1_[i1], rhs().data1_[i2]);
            std::swap(lhs().data2_[i1], rhs().data2_[i2]);
            std::swap(lhs().data3_[i1], rhs().data3_[i2]);
        }

    private:
        size_type size_;
        V1& data1_;
        V2& data2_;
        V3& data3_;

        // friend class value_type;
        friend class index_triple<self_type>;
    };

}}}

#endif



