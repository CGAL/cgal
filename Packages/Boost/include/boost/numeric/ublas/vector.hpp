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

#ifndef BOOST_UBLAS_VECTOR_H
#define BOOST_UBLAS_VECTOR_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_assign.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

// Iterators based on ideas of Jeremy Siek

namespace boost { namespace numeric { namespace ublas {

    // Array based vector class
    template<class T, class A>
    class vector:
        public vector_expression<vector<T, A> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<vector<T, A> >::operator ();
#endif
        typedef concrete_tag simd_category;
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef const T *const_pointer;
        typedef T *pointer;
        typedef A array_type;
        typedef const A const_array_type;
        typedef const vector<T, A> const_self_type;
        typedef vector<T, A> self_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const_self_type> const_closure_type;
#else
        typedef const vector_reference<const_self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef typename A::const_iterator const_iterator_type;
        typedef typename A::iterator iterator_type;
        typedef dense_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        vector ():
            vector_expression<self_type> (),
            size_ (0), data_ (0) {}
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        vector (size_type size):
            vector_expression<self_type> (),
            size_ (size), data_ (0) {
            resize (size);
        }
        BOOST_UBLAS_INLINE
        vector (size_type size, const array_type &data):
            vector_expression<self_type> (),
            size_ (size), data_ (data) {}
        BOOST_UBLAS_INLINE
        vector (const vector &v):
            vector_expression<self_type> (),
            size_ (v.size_), data_ (v.data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        vector (const vector_expression<AE> &ae):
            vector_expression<self_type> (),
            size_ (ae ().size ()), data_ (0) {
#ifndef BOOST_UBLAS_TYPE_CHECK
            resize (ae ().size (), false);
#else
            resize (ae ().size (), true);
#endif
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
        }

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }
        BOOST_UBLAS_INLINE
        const_array_type &data () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        array_type &data () {
            return data_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            size_ = size;
            detail::resize (data (), size, preserve);
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            return data () [i];
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) {
            return data () [i];
        }

        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            return (*this) (i);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) {
            return (*this) (i);
        }

        // Assignment
        BOOST_UBLAS_INLINE
        vector &operator = (const vector &v) {
            // Precondition for container relaxed as requested during review.
            // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
            size_ = v.size_;
            data () = v.data ();
            return *this;
        }
        BOOST_UBLAS_INLINE
        vector &assign_temporary (vector &v) {
            swap (v);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector &operator = (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (ae));
#else
            // return assign (self_type (ae));
            self_type temporary (ae);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector &reset (const vector_expression<AE> &ae) {
            self_type temporary (ae);
            resize (temporary.size (), false);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector &assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector &operator += (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (*this + ae));
#else
            // return assign (self_type (*this + ae));
            self_type temporary (*this + ae);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector &plus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector &operator -= (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (*this - ae));
#else
            // return assign (self_type (*this - ae));
            self_type temporary (*this - ae);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector &minus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        vector &operator *= (const AT &at) {
            vector_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        vector &operator /= (const AT &at) {
            vector_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                std::swap (size_, v.size_);
                data ().swap (v.data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (vector &v1, vector &v2) {
            v1.swap (v2);
        }
#endif

        // Element insertion and erasure
        // These functions should work with std::vector.
        // Thanks to Kresimir Fresl for spotting this.
        BOOST_UBLAS_INLINE
        void insert (size_type i, const_reference t) {
            // FIXME: only works for EqualityComparable value types.
            // BOOST_UBLAS_CHECK (data () [i] == value_type (), bad_index ());
            // data ().insert (data ().begin () + i, t);
            data () [i] = t;
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i) {
            // data ().erase (data ().begin () + i);
            data () [i] = value_type ();
        }
        BOOST_UBLAS_INLINE
        void clear () {
            // data ().clear ();
            std::fill (data ().begin (), data ().end (), value_type ());
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_iterator<self_type, dense_random_access_iterator_tag> iterator;
        typedef indexed_const_iterator<self_type, dense_random_access_iterator_tag> const_iterator;
#else
        class const_iterator;
        class iterator;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator (*this, data ().begin () + i);
#else
            return const_iterator (*this, i);
#endif
        }
        BOOST_UBLAS_INLINE
        iterator find (size_type i) {
#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return iterator (*this, data ().begin () + i);
#else
            return iterator (*this, i);
#endif
        }

        // Iterators simply are pointers.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator:
            public container_const_reference<vector>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               const_iterator, value_type> {
        public:
            typedef dense_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename vector::difference_type difference_type;
            typedef typename vector::value_type value_type;
            typedef typename vector::const_reference reference;
            typedef typename vector::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &v, const const_iterator_type &it):
                container_const_reference<self_type> (v), it_ (it) {}
            BOOST_UBLAS_INLINE
#ifndef BOOST_UBLAS_QUALIFIED_TYPENAME
            const_iterator (const iterator &it):
#else
            const_iterator (const typename self_type::iterator &it):
#endif
                container_const_reference<self_type> (it ()), it_ (it.it_) {}

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
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return *it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_ - (*this) ().begin ().it_;
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                container_const_reference<self_type>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator_type it_;

            friend class iterator;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size_);
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator:
            public container_reference<vector>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               iterator, value_type> {
        public:
            typedef dense_random_access_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename vector::difference_type difference_type;
            typedef typename vector::value_type value_type;
            typedef typename vector::reference reference;
            typedef typename vector::pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator ():
                container_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator (self_type &v, const iterator_type &it):
                container_reference<self_type> (v), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator &operator ++ () {
                ++ it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator &operator -- () {
                -- it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator &operator += (difference_type n) {
                it_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator &operator -= (difference_type n) {
                it_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return *it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_ - (*this) ().begin ().it_;
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator &operator = (const iterator &it) {
                container_reference<self_type>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            iterator_type it_;

            friend class const_iterator;
        };
#endif

        BOOST_UBLAS_INLINE
        iterator begin () {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        iterator end () {
            return find (size_);
        }

        // Reverse iterator

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
        typedef reverse_iterator_base<const_iterator> const_reverse_iterator;
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
        typedef reverse_iterator_base<iterator, value_type, reference> reverse_iterator;
#else
        typedef reverse_iterator_base<iterator> reverse_iterator;
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
        array_type data_;
    };

    // Bounded vector class
    template<class T, std::size_t N>
    class bounded_vector:
        public vector<T, bounded_array<T, N> > {
    public:
#ifndef BOOST_UBLAS_NO_DERIVED_HELPERS
        BOOST_UBLAS_USING vector<T, bounded_array<T, N> >::operator =;
#endif
        BOOST_STATIC_CONSTANT (std::size_t,  max_size = N);
        typedef vector<T, bounded_array<T, N> > vector_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        bounded_vector ():
            vector_type (N) {}
        BOOST_UBLAS_INLINE
        bounded_vector (std::size_t size):
            vector_type (size) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        bounded_vector (const vector_expression<AE> &ae):
            vector_type (ae) {}
        BOOST_UBLAS_INLINE
        ~bounded_vector () {}

        // Assignment
        BOOST_UBLAS_INLINE
        bounded_vector &operator = (const bounded_vector &v) {
            vector_type::operator = (v);
            return *this;
        }
    };

    // Unit vector class
    template<class T>
    class unit_vector:
        public vector_expression<unit_vector<T> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<unit_vector<T> >::operator ();
#endif
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef const T *const_pointer;
        typedef T *pointer;
        typedef const unit_vector<T> const_self_type;
        typedef unit_vector<T> self_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const_self_type> const_closure_type;
#else
        typedef const vector_reference<const_self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef size_type const_iterator_type;
        typedef packed_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        unit_vector ():
            vector_expression<self_type> (),
            size_ (0), index_ (0) {}
        BOOST_UBLAS_INLINE
        unit_vector (size_type size, size_type index):
            vector_expression<self_type> (),
            size_ (size), index_ (index) {}
        BOOST_UBLAS_INLINE
        unit_vector (const unit_vector &v):
            vector_expression<self_type> (),
            size_ (v.size_), index_ (v.index_) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_; 
        }
        BOOST_UBLAS_INLINE
        size_type index () const { 
            return index_; 
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            size_ = size;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            return i == index_ ? one_ : zero_; 
        }

        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            return (*this) (i); 
        }

        // Assignment
        BOOST_UBLAS_INLINE
        unit_vector &operator = (const unit_vector &v) {
            // Precondition for container relaxed as requested during review.
            // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
            size_ = v.size_;
            index_ = v.index_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        unit_vector &assign_temporary (unit_vector &v) { 
            swap (v);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (unit_vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                std::swap (size_, v.size_);
                std::swap (index_, v.index_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (unit_vector &v1, unit_vector &v2) {
            v1.swap (v2);
        }
#endif

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator<self_type, packed_random_access_iterator_tag> iterator;
        typedef indexed_const_iterator<self_type, packed_random_access_iterator_tag> const_iterator;
#else
        class const_iterator;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
            i = std::max (i, index_);
            i = std::min (i, index_ + 1);
            return const_iterator (*this, i);
        }

        // Iterators simply are pointers.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator:
            public container_const_reference<unit_vector>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename unit_vector::difference_type difference_type;
            typedef typename unit_vector::value_type value_type;
            typedef typename unit_vector::const_reference reference;
            typedef typename unit_vector::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<unit_vector> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const unit_vector &v, const const_iterator_type &it):
                container_const_reference<unit_vector> (v), it_ (it) {}

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
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return (*this) () (index ());
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_;
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                container_const_reference<unit_vector>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator_type it_;
        };

        typedef const_iterator iterator;
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size_);
        }

        // Reverse iterator

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
        typedef reverse_iterator_base<const_iterator> const_reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator rbegin () const {
            return const_reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator rend () const {
            return const_reverse_iterator (begin ());
        }

    private:
        size_type size_;
        size_type index_;
        static value_type zero_;
        static value_type one_;
    };

    template<class T>
    typename unit_vector<T>::value_type unit_vector<T>::zero_ =
        unit_vector<T>::value_type ();
    template<class T>
    typename unit_vector<T>::value_type unit_vector<T>::one_ =
        unit_vector<T>::value_type (1);

    // Zero vector class
    template<class T>
    class zero_vector:
        public vector_expression<zero_vector<T> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<zero_vector<T> >::operator ();
#endif
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef const T *const_pointer;
        typedef T *pointer;
        typedef const zero_vector<T> const_self_type;
        typedef zero_vector<T> self_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const_self_type> const_closure_type;
#else
        typedef const vector_reference<const_self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef size_type const_iterator_type;
        typedef sparse_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        zero_vector ():
            vector_expression<self_type> (),
            size_ (0) {}
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        zero_vector (size_type size):
            vector_expression<self_type> (),
            size_ (size) {}
        BOOST_UBLAS_INLINE
        zero_vector (const zero_vector &v):
            vector_expression<self_type> (),
            size_ (v.size_) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            size_ = size;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type /* i */) const {
            return zero_;
        }

        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            return (*this) (i);
        }

        // Assignment
        BOOST_UBLAS_INLINE
        zero_vector &operator = (const zero_vector &v) {
            // Precondition for container relaxed as requested during review.
            // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
            size_ = v.size_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        zero_vector &assign_temporary (zero_vector &v) {
            swap (v);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (zero_vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                std::swap (size_, v.size_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (zero_vector &v1, zero_vector &v2) {
            v1.swap (v2);
        }
#endif

        class const_iterator;

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
            return const_iterator (*this, i);
        }

        // Iterators simply are pointers.

        class const_iterator:
            public container_const_reference<zero_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename zero_vector::difference_type difference_type;
            typedef typename zero_vector::value_type value_type;
            typedef typename zero_vector::const_reference reference;
            typedef typename zero_vector::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &v, const const_iterator_type &it):
                container_const_reference<self_type> (v), it_ (it) {}

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

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return (*this) () (index ());
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_;
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                container_const_reference<self_type>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }

        private:
            const_iterator_type it_;
        };

        typedef const_iterator iterator;

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size_);
        }

        // Reverse iterator

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
        typedef reverse_iterator_base<const_iterator> const_reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator rbegin () const {
            return const_reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator rend () const {
            return const_reverse_iterator (begin ());
        }

    private:
        size_type size_;
        static value_type zero_;
    };

    template<class T>
    typename zero_vector<T>::value_type zero_vector<T>::zero_ =
        zero_vector<T>::value_type ();

    // Scalar vector class
    template<class T>
    class scalar_vector:
        public vector_expression<scalar_vector<T> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<scalar_vector<T> >::operator ();
#endif
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef const T *const_pointer;
        typedef T *pointer;
        typedef const scalar_vector<T> const_self_type;
        typedef scalar_vector<T> self_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const_self_type> const_closure_type;
#else
        typedef const vector_reference<const_self_type> const_closure_type;
#endif
        typedef size_type const_iterator_type;
        typedef dense_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        scalar_vector ():
            vector_expression<self_type> (),
            size_ (0), value_ () {}
        BOOST_UBLAS_INLINE
        scalar_vector (size_type size, const value_type &value):
            vector_expression<self_type> (),
            size_ (size), value_ (value) {}
        BOOST_UBLAS_INLINE
        scalar_vector (const scalar_vector &v):
            vector_expression<self_type> (),
            size_ (v.size_), value_ (v.value_) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const { 
            return size_; 
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            size_ = size;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            return value_; 
        }

        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const { 
            return value_; 
        }

        // Assignment
        BOOST_UBLAS_INLINE
        scalar_vector &operator = (const scalar_vector &v) {
            // Precondition for container relaxed as requested during review.
            // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
            size_ = v.size_;
            value_ = v.value_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        scalar_vector &assign_temporary (scalar_vector &v) { 
            swap (v);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (scalar_vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                std::swap (size_, v.size_);
                std::swap (value_, v.value_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (scalar_vector &v1, scalar_vector &v2) {
            v1.swap (v2);
        }
#endif

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator<self_type, dense_random_access_iterator_tag> iterator;
        typedef indexed_const_iterator<self_type, dense_random_access_iterator_tag> const_iterator;
#else
        class const_iterator;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
            return const_iterator (*this, i);
        }

        // Iterators simply are pointers.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator:
            public container_const_reference<scalar_vector>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               const_iterator, value_type> {
        public:
            typedef dense_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename scalar_vector::difference_type difference_type;
            typedef typename scalar_vector::value_type value_type;
            typedef typename scalar_vector::const_reference reference;
            typedef typename scalar_vector::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<scalar_vector> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const scalar_vector &v, const const_iterator_type &it):
                container_const_reference<scalar_vector> (v), it_ (it) {}

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
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return (*this) () (index ());
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_;
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                container_const_reference<scalar_vector>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator_type it_;
        };

        typedef const_iterator iterator;
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size_);
        }

        // Reverse iterator

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
        typedef reverse_iterator_base<const_iterator> const_reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator rbegin () const {
            return const_reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator rend () const {
            return const_reverse_iterator (begin ());
        }

    private:
        size_type size_;
        value_type value_;
    };

    // Array based vector class 
    template<class T, std::size_t N>
    class c_vector:
        public vector_expression<c_vector<T, N> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<c_vector<T, N> >::operator ();
#endif
        typedef concrete_tag simd_category;
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef const T *const_pointer;
        typedef T *pointer;
        typedef const c_vector<T, N> const_self_type;
        typedef c_vector<T, N> self_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const_self_type> const_closure_type;
#else
        typedef const vector_reference<const_self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef const T *const_iterator_type;
        typedef T *iterator_type;
        typedef dense_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        c_vector ():
            size_ (N) /* , data_ () */ {}
        BOOST_UBLAS_EXPLICIT BOOST_UBLAS_INLINE
        c_vector (size_type size):
            size_ (size) /* , data_ () */ {
            if (size_ > N) 
                // Raising exceptions abstracted as requested during review.
                // throw std::bad_alloc ();
                bad_size ().raise ();
        }
        BOOST_UBLAS_INLINE
        c_vector (const c_vector &v): 
            size_ (v.size_) /* , data_ () */ {
            if (size_ > N) 
                // Raising exceptions abstracted as requested during review.
                // throw std::bad_alloc ();
                bad_size ().raise ();
            *this = v;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        c_vector (const vector_expression<AE> &ae):
            size_ (ae ().size ()) /* , data_ () */ {
            if (size_ > N)
                // Raising exceptions abstracted as requested during review.
                // throw std::bad_alloc ();
                bad_size ().raise ();
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
        }

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }
        BOOST_UBLAS_INLINE
        const_pointer data () const {
            return reinterpret_cast<const_pointer> (data_);
        }
        BOOST_UBLAS_INLINE
        pointer data () {
            return reinterpret_cast<pointer> (data_);
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            if (size > N)
                // Raising exceptions abstracted as requested during review.
                // throw std::bad_alloc ();
                bad_size ().raise ();
            // The content of the array is intentionally not copied.
            size_ = size;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_,  bad_index ());
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

        // Assignment
        BOOST_UBLAS_INLINE
        c_vector &operator = (const c_vector &v) {
            // Precondition for container relaxed as requested during review.
            // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
            size_ = v.size_;
            std::copy (v.data_, v.data_ + v.size_, data_);
            return *this;
        }
        BOOST_UBLAS_INLINE
        c_vector &assign_temporary (c_vector &v) {
            swap (v);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        c_vector &operator = (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (ae));
#else
            // return assign (self_type (ae));
            self_type temporary (ae);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        c_vector &reset (const vector_expression<AE> &ae) {
            self_type temporary (ae);
            resize (temporary.size (), false);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        c_vector &assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        c_vector &operator += (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (*this + ae));
#else
            // return assign (self_type (*this + ae));
            self_type temporary (*this + ae);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        c_vector &plus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        c_vector &operator -= (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (*this - ae));
#else
            // return assign (self_type (*this - ae));
            self_type temporary (*this - ae);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        c_vector &minus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        c_vector &operator *= (const AT &at) {
            vector_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        c_vector &operator /= (const AT &at) {
            vector_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (c_vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                std::swap (size_, v.size_);
                std::swap_ranges (data_, data_ + size_, v.data_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (c_vector &v1, c_vector &v2) {
            v1.swap (v2);
        }
#endif

        // Element insertion and erasure
        BOOST_UBLAS_INLINE
        void insert (size_type i, const_reference t) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            BOOST_UBLAS_CHECK (data_ [i] == value_type (), bad_index ());
            data_ [i] = t;
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            data_ [i] = value_type ();
        }
        BOOST_UBLAS_INLINE
        void clear () {
            std::fill (data_, data_ + size_, value_type ());
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_iterator<self_type, dense_random_access_iterator_tag> iterator;
        typedef indexed_const_iterator<self_type, dense_random_access_iterator_tag> const_iterator;
#else
        class const_iterator;
        class iterator;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator (*this, &data_ [i]);
#else
            return const_iterator (*this, i);
#endif
        }
        BOOST_UBLAS_INLINE
        iterator find (size_type i) {
#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return iterator (*this, &data_ [i]);
#else
            return iterator (*this, i);
#endif
        }

        // Iterators simply are pointers.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator:
            public container_const_reference<c_vector>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               const_iterator, value_type> {
        public:
            typedef dense_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename c_vector::difference_type difference_type;
            typedef typename c_vector::value_type value_type;
            typedef typename c_vector::const_reference reference;
            typedef typename c_vector::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &v, const const_iterator_type &it):
                container_const_reference<self_type> (v), it_ (it) {}
            BOOST_UBLAS_INLINE
#ifndef BOOST_UBLAS_QUALIFIED_TYPENAME
            const_iterator (const iterator &it):
#else
            const_iterator (const typename self_type::iterator &it):
#endif
                container_const_reference<self_type> (it ()), it_ (it.it_) {}

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
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return *it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                const self_type &v = (*this) ();
                return it_ - v.begin ().it_;
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                container_const_reference<self_type>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator_type it_;

            friend class iterator;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size_);
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator:
            public container_reference<c_vector>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               iterator, value_type> {
        public:
            typedef dense_random_access_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename c_vector::difference_type difference_type;
            typedef typename c_vector::value_type value_type;
            typedef typename c_vector::reference reference;
            typedef typename c_vector::pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator ():
                container_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator (self_type &v, const iterator_type &it):
                container_reference<self_type> (v), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator &operator ++ () {
                ++ it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator &operator -- () {
                -- it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator &operator += (difference_type n) {
                it_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator &operator -= (difference_type n) {
                it_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return *it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                self_type &v = (*this) ();
                return it_ - v.begin ().it_;
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator &operator = (const iterator &it) {
                container_reference<self_type>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            iterator_type it_;

            friend class const_iterator;
        };
#endif

        BOOST_UBLAS_INLINE
        iterator begin () {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        iterator end () {
            return find (size_);
        }

        // Reverse iterator

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
        typedef reverse_iterator_base<const_iterator> const_reverse_iterator;
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
        typedef reverse_iterator_base<iterator, value_type, reference> reverse_iterator;
#else
        typedef reverse_iterator_base<iterator> reverse_iterator;
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
        value_type data_ [N];
    };

}}}

#endif



















