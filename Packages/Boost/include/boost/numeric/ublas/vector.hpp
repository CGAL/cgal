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

#include <vector>

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
        typedef typename A::size_type size_type;
        typedef typename A::difference_type difference_type;
        typedef T value_type;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef A array_type;
    private:
        typedef T *pointer;
        typedef vector<T, A> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const self_type> const_closure_type;
#else
        typedef const vector_reference<const self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef self_type vector_temporary_type;
        typedef dense_tag storage_category;
        typedef concrete_tag simd_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        vector ():
            vector_expression<self_type> (),
            data_ () {}
        explicit BOOST_UBLAS_INLINE
        vector (size_type size):
            vector_expression<self_type> (),
            data_ (size) {
        }
        BOOST_UBLAS_INLINE
        vector (size_type size, const array_type &data):
            vector_expression<self_type> (),
            data_ (data) {}
        BOOST_UBLAS_INLINE
        vector (const vector &v):
            vector_expression<self_type> (),
            data_ (v.data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        vector (const vector_expression<AE> &ae):
            vector_expression<self_type> (),
            data_ (ae ().size ()) {
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
        }

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const {
            return data_.size ();
        }
        BOOST_UBLAS_INLINE
        const array_type &data () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        array_type &data () {
            return data_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            if (preserve)
                data ().resize (size, BOOST_UBLAS_TYPENAME A::value_type (0));
            else
                data ().resize (size);
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
            data () = v.data ();
            return *this;
        }
#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
        template<class A2>          // Generic vector assignment without temporary
        BOOST_UBLAS_INLINE
        vector &operator = (const vector<T, A2> &v) {
            resize (v.size ());
            assign (v);
            return *this;
        }
#endif
        BOOST_UBLAS_INLINE
        vector &assign_temporary (vector &v) {
            swap (v);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector &operator = (const vector_expression<AE> &ae) {
            // return assign (self_type (ae));
            self_type temporary (ae);
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
            // return assign (self_type (*this + ae));
            self_type temporary (*this + ae);
            return assign_temporary (temporary);
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
            // return assign (self_type (*this - ae));
            self_type temporary (*this - ae);
            return assign_temporary (temporary);
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
            if (this != &v) {
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
            // BOOST_UBLAS_CHECK (data () [i] == value_type (0), bad_index ());
            // Previously: data ().insert (data ().begin () + i, t);
            data () [i] = t;
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i) {
            // Previously: data ().erase (data ().begin () + i);
            data () [i] = value_type (0);
        }
        BOOST_UBLAS_INLINE
        void clear () {
            // Previously: data ().clear ();
            std::fill (data ().begin (), data ().end (), value_type (0));
        }

        // Iterator types
    private:
        // Use the storage array iterator
        typedef typename A::const_iterator const_iterator_type;
        typedef typename A::iterator iterator_type;

    public:
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


#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator:
            public container_const_reference<vector>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               const_iterator, value_type, difference_type> {
        public:
            typedef dense_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename vector::difference_type difference_type;
            typedef typename vector::value_type value_type;
            typedef typename vector::const_reference reference;
            typedef const typename vector::pointer pointer;
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
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_, bad_index ());
                return *it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_, bad_index ());
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
            return find (data_.size ());
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator:
            public container_reference<vector>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               iterator, value_type, difference_type> {
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
                BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_ , bad_index ());
                return *it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_ , bad_index ());
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
            return find (data_.size ());
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
        array_type data_;
    };

    // Bounded vector class
    template<class T, std::size_t N>
    class bounded_vector:
        public vector<T, bounded_array<T, N> > {

        typedef vector<T, bounded_array<T, N> > vector_type;
    public:
        typedef typename vector_type::size_type size_type;
        BOOST_STATIC_CONSTANT (size_type,  max_size = N);

        // Construction and destruction
        BOOST_UBLAS_INLINE
        bounded_vector ():
            vector_type (N) {}
        BOOST_UBLAS_INLINE
        bounded_vector (size_type size):
            vector_type (size) {}
        BOOST_UBLAS_INLINE
        bounded_vector (const bounded_vector &v):
            vector_type (v) {}
#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
        template<class A2>              // Allow vector<T,bounded_array<N> construction
        BOOST_UBLAS_INLINE
        bounded_vector (const vector<T, A2> &v):
            vector_type (v) {}
#endif
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
        template<class A2>              // Generic vector assignment
        BOOST_UBLAS_INLINE
        bounded_vector &operator = (const vector<T, A2> &v) {
            vector_type::operator = (v);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        bounded_vector &operator = (const vector_expression<AE> &ae) {
            vector_type::operator = (ae);
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
        typedef const T &const_reference;
        typedef T &reference;
    private:
        typedef const T *const_pointer;
        typedef unit_vector<T> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const self_type> const_closure_type;
#else
        typedef const vector_reference<const self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef packed_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        unit_vector ():
            vector_expression<self_type> (),
            size_ (0), index_ (0) {}
        BOOST_UBLAS_INLINE
        explicit unit_vector (size_type size, size_type index = 0):
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
            if (i == index_)
                return one_;
            else
                return zero_;
        }

        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            return (*this) (i);
        }

        // Assignment
        BOOST_UBLAS_INLINE
        unit_vector &operator = (const unit_vector &v) {
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
            if (this != &v) {
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

        // Iterator types
    private:
        // Use an index
        typedef size_type const_iterator_type;
    public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator<self_type, packed_random_access_iterator_tag> iterator;
        typedef indexed_const_iterator<self_type, packed_random_access_iterator_tag> const_iterator;
#else
        class const_iterator;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
            i = (std::max) (i, index_);
            i = (std::min) (i, index_ + 1);
            return const_iterator (*this, i);
        }

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
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (it_ < (*this) ().size (), bad_index ());
                return (*this) () (index ());
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                BOOST_UBLAS_CHECK (it_ < (*this) ().size (), bad_index ());
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
        typedef const value_type const_value_type;
        static const_value_type zero_;
        static const_value_type one_;
    };

    template<class T>
    typename unit_vector<T>::const_value_type unit_vector<T>::zero_ 
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME unit_vector<T>::const_value_type
#endif
        (0);
    template<class T>
    typename unit_vector<T>::const_value_type unit_vector<T>::one_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME unit_vector<T>::const_value_type
#endif
        (1);

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
        typedef const T &const_reference;
        typedef T &reference;
    private:
        typedef const T *const_pointer;
        typedef zero_vector<T> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const self_type> const_closure_type;
#else
        typedef const vector_reference<const self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef sparse_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        zero_vector ():
            vector_expression<self_type> (),
            size_ (0) {}
        explicit BOOST_UBLAS_INLINE
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
            if (this != &v) {
                std::swap (size_, v.size_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (zero_vector &v1, zero_vector &v2) {
            v1.swap (v2);
        }
#endif

        // Iterator types
    private:
        // Use an index
        typedef size_type const_iterator_type;
    public:
        class const_iterator;

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
            return const_iterator (*this, i);
        }


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
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (it_ < (*this) ().size (), bad_index ());
                return (*this) () (index ());
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                BOOST_UBLAS_CHECK (it_ < (*this) ().size (), bad_index ());
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
        typedef const value_type const_value_type;
        static const_value_type zero_;
    };

    template<class T>
    typename zero_vector<T>::const_value_type zero_vector<T>::zero_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME zero_vector<T>::const_value_type
#endif
        (0);

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
        typedef const T &const_reference;
        typedef T &reference;
    private:
        typedef const T *const_pointer;
        typedef scalar_vector<T> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const self_type> const_closure_type;
#else
        typedef const vector_reference<const self_type> const_closure_type;
#endif
        typedef dense_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        scalar_vector ():
            vector_expression<self_type> (),
            size_ (0), value_ () {}
        BOOST_UBLAS_INLINE
        explicit scalar_vector (size_type size, const value_type &value = value_type(1)):
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
            if (this != &v) {
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

        // Iterator types
    private:
        // Use an index
        typedef size_type const_iterator_type;

    public:
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
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (it_ < (*this) ().size (), bad_index ());
                return (*this) () (index ());
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                BOOST_UBLAS_CHECK (it_ < (*this) ().size (), bad_index ());
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
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        typedef const T &const_reference;
        typedef T &reference;
        typedef value_type array_type[N];
        typedef T *pointer;
        typedef const T *const_pointer;
    private:
        typedef c_vector<T, N> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const self_type> const_closure_type;
#else
        typedef const vector_reference<const self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef self_type vector_temporary_type;
        typedef dense_tag storage_category;
        typedef concrete_tag simd_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        c_vector ():
            size_ (N) /* , data_ () */ {}
        explicit BOOST_UBLAS_INLINE
        c_vector (size_type size):
            size_ (size) /* , data_ () */ {
            if (size_ > N)
                bad_size ().raise ();
            std::fill (data_, data_ + size_, value_type (0));
        }
        BOOST_UBLAS_INLINE
        c_vector (const c_vector &v):
            size_ (v.size_) /* , data_ () */ {
            if (size_ > N)
                bad_size ().raise ();
            *this = v;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        c_vector (const vector_expression<AE> &ae):
            size_ (ae ().size ()) /* , data_ () */ {
            if (size_ > N)
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
            return data_;
        }
        BOOST_UBLAS_INLINE
        pointer data () {
            return data_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            if (size > N)
                bad_size ().raise ();
            std::fill (data_ + (std::min) (size, size_), data_ + size, value_type (0));
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
            // return assign (self_type (ae));
            self_type temporary (ae);
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
            // return assign (self_type (*this + ae));
            self_type temporary (*this + ae);
            return assign_temporary (temporary);
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
            // return assign (self_type (*this - ae));
            self_type temporary (*this - ae);
            return assign_temporary (temporary);
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
            BOOST_UBLAS_CHECK (data_ [i] == value_type (0), bad_index ());
            data_ [i] = t;
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            data_ [i] = value_type (0);
        }
        BOOST_UBLAS_INLINE
        void clear () {
            std::fill (data_, data_ + size_, value_type (0));
        }

        // Iterator types
    private:
        // Use pointers for iterator
        typedef const_pointer const_iterator_type;
        typedef pointer iterator_type;

    public:
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
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_, bad_index ());
                return *it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_, bad_index ());
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
                BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_, bad_index ());
                return *it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_, bad_index ());
                // EDG won't allow const self_type it doesn't allow friend access to it_
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
        array_type data_;
    };

}}}

#endif
