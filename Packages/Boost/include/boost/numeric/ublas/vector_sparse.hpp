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

#ifndef BOOST_UBLAS_VECTOR_SPARSE_H
#define BOOST_UBLAS_VECTOR_SPARSE_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/storage_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

// Iterators based on ideas of Jeremy Siek

namespace boost { namespace numeric { namespace ublas {

#ifdef BOOST_UBLAS_STRICT_VECTOR_SPARSE

    template<class V>
    class sparse_vector_element:
       public container_reference<V> {
    public:
        typedef V vector_type;
        typedef typename V::size_type size_type;
        typedef typename V::value_type value_type;
        // typedef const value_type &const_reference;
        typedef typename type_traits<value_type>::const_reference const_reference;
        typedef value_type &reference;
        typedef value_type *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        sparse_vector_element (const value_type &d):
            container_reference<vector_type> (), it_ (), i_ (), d_ (d), dirty_ (false) {
            external_logic ().raise ();
        }
        BOOST_UBLAS_INLINE
        sparse_vector_element (vector_type &v, pointer it, size_type i):
            container_reference<vector_type> (v), it_ (it), i_ (i), d_ (*it), dirty_ (false) {}
        BOOST_UBLAS_INLINE
        sparse_vector_element (vector_type &v, size_type i):
            container_reference<vector_type> (v), it_ (), i_ (i), d_ (), dirty_ (false) {
            pointer it = (*this) ().find_element (i_);
            if (it)
                d_ = *it;
        }
        BOOST_UBLAS_INLINE
        sparse_vector_element (const sparse_vector_element &p):
            container_reference<vector_type> (p), it_ (p.it_), i_ (p.i_), d_ (p.d_), dirty_ (p.dirty_) {}
        BOOST_UBLAS_INLINE
        ~sparse_vector_element () {
            if (dirty_) {
                if (! it_) {
                    it_ = (*this) ().find_element (i_);
                    if (! it_)
                        (*this) ().insert (i_, d_);
                    else
                        *it_ = d_;
                } else
                    *it_ = d_;
            }
        }

        // Assignment
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator = (const D &d) {
            d_ = d;
            dirty_ = true;
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator = (const sparse_vector_element &p) {
            d_ = p.d_;
            dirty_ = true;
            return *this;
        }
        template<class OV>
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator = (const sparse_vector_element<OV> &p) {
            d_ = p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator += (const D &d) {
            d_ += d;
            dirty_ = true;
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator += (const sparse_vector_element &p) {
            d_ += p.d_;
            dirty_ = true;
            return *this;
        }
        template<class OV>
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator += (const sparse_vector_element<OV> &p) {
            d_ += p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator -= (const D &d) {
            d_ -= d;
            dirty_ = true;
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator -= (const sparse_vector_element &p) {
            d_ -= p.d_;
            dirty_ = true;
            return *this;
        }
        template<class OV>
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator -= (const sparse_vector_element<OV> &p) {
            d_ -= p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator *= (const D &d) {
            d_ *= d;
            dirty_ = true;
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator *= (const sparse_vector_element &p) {
            d_ *= p.d_;
            dirty_ = true;
            return *this;
        }
        template<class OV>
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator *= (const sparse_vector_element<OV> &p) {
            d_ *= p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator /= (const D &d) {
            d_ /= d;
            dirty_ = true;
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator /= (const sparse_vector_element &p) {
            d_ /= p.d_;
            dirty_ = true;
            return *this;
        }
        template<class OV>
        BOOST_UBLAS_INLINE
        sparse_vector_element &operator /= (const sparse_vector_element<OV> &p) {
            d_ /= p.d_;
            dirty_ = true;
            return *this;
        }

        // Conversion
        BOOST_UBLAS_INLINE
        operator const_reference () const {
            return d_;
        }
#ifdef BOOST_UBLAS_DEPRECATED
        BOOST_UBLAS_INLINE
        operator reference () {
            dirty_ = true;
            return d_;
        }
#endif

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (sparse_vector_element p) {
            if (this != &p) {
                dirty_ = true;
                p.dirty_ = true;
                std::swap (d_, p.d_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (sparse_vector_element p1, sparse_vector_element p2) {
            p1.swap (p2);
        }
#endif

    private:
        pointer it_;
        size_type i_;
        value_type d_;
        bool dirty_;
    };

    template<class V>
    struct type_traits<sparse_vector_element<V> > {
        typedef typename V::value_type element_type;
        typedef type_traits<sparse_vector_element<V> > self_type;
        typedef typename type_traits<element_type>::value_type value_type;
        typedef typename type_traits<element_type>::const_reference const_reference;
        typedef sparse_vector_element<V> reference;
        typedef typename type_traits<element_type>::real_type real_type;
        typedef typename type_traits<element_type>::precision_type precision_type;

        BOOST_STATIC_CONSTANT (std::size_t, plus_complexity = type_traits<element_type>::plus_complexity);
        BOOST_STATIC_CONSTANT (std::size_t, multiplies_complexity = type_traits<element_type>::multiplies_complexity);

        static
        BOOST_UBLAS_INLINE
        real_type real (const_reference t) {
            return type_traits<element_type>::real (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type imag (const_reference t) {
            return type_traits<element_type>::imag (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type conj (const_reference t) {
            return type_traits<element_type>::conj (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type abs (const_reference t) {
            return type_traits<element_type>::abs (t);
        }
        static
        BOOST_UBLAS_INLINE
        value_type sqrt (const_reference t) {
            return type_traits<element_type>::sqrt (t);
        }

        static
        BOOST_UBLAS_INLINE
        real_type norm_1 (const_reference t) {
            return type_traits<element_type>::norm_1 (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_2 (const_reference t) {
            return type_traits<element_type>::norm_2 (t);
        }
        static
        BOOST_UBLAS_INLINE
        real_type norm_inf (const_reference t) {
            return type_traits<element_type>::norm_inf (t);
        }

        static
        BOOST_UBLAS_INLINE
        bool equals (const_reference t1, const_reference t2) {
            return type_traits<element_type>::equals (t1, t2);
        }
    };

    template<class V1, class T2>
    struct promote_traits<sparse_vector_element<V1>, T2> {
        typedef typename promote_traits<typename sparse_vector_element<V1>::value_type, T2>::promote_type promote_type;
    };
    template<class T1, class V2>
    struct promote_traits<T1, sparse_vector_element<V2> > {
        typedef typename promote_traits<T1, typename sparse_vector_element<V2>::value_type>::promote_type promote_type;
    };
    template<class V1, class V2>
    struct promote_traits<sparse_vector_element<V1>, sparse_vector_element<V2> > {
        typedef typename promote_traits<typename sparse_vector_element<V1>::value_type,
                                        typename sparse_vector_element<V2>::value_type>::promote_type promote_type;
    };

#endif

    // Array based sparse vector class
    template<class T, class A>
    class sparse_vector:
        public vector_expression<sparse_vector<T, A> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<sparse_vector<T, A> >::operator ();
#endif
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
#if ! defined (BOOST_UBLAS_STRICT_STORAGE_SPARSE) && ! defined (BOOST_UBLAS_STRICT_VECTOR_SPARSE)
        typedef T &reference;
#elif defined (BOOST_UBLAS_STRICT_VECTOR_SPARSE)
        typedef sparse_vector_element<sparse_vector<T, A> > reference;
#elif defined (BOOST_UBLAS_STRICT_STORAGE_SPARSE)
        typedef typename map_traits<A>::reference reference;
#endif
        typedef const T *const_pointer;
        typedef T *pointer;
        typedef A array_type;
        typedef const A const_array_type;
        typedef const sparse_vector<T, A> const_self_type;
        typedef sparse_vector<T, A> self_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const_self_type> const_closure_type;
#else
        typedef const vector_reference<const_self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef typename A::const_iterator const_iterator_type;
        typedef typename A::iterator iterator_type;
        typedef sparse_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        sparse_vector ():
            vector_expression<self_type> (),
            size_ (0), non_zeros_ (0), data_ () {}
        BOOST_UBLAS_INLINE
        sparse_vector (size_type size, size_type non_zeros = 0):
            vector_expression<self_type> (),
            size_ (size), non_zeros_ (non_zeros), data_ () {
            reserve (non_zeros_);
        }
        BOOST_UBLAS_INLINE
        sparse_vector (const sparse_vector &v):
            vector_expression<self_type> (),
            size_ (v.size_), non_zeros_ (v.non_zeros_), data_ (v.data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector (const vector_expression<AE> &ae, size_type non_zeros = 0):
            vector_expression<self_type> (),
            size_ (ae ().size ()), non_zeros_ (non_zeros), data_ () {
            reserve (non_zeros_);
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
        }

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }
        BOOST_UBLAS_INLINE
        size_type non_zeros () const {
            return data_.size ();
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
        void resize (size_type size, size_type non_zeros = 0) {
            size_ = size;
            non_zeros_ = std::max (non_zeros, size_type (1));
            non_zeros_ = std::min (non_zeros_, size_);
            detail::reserve (data (), non_zeros_);
            data ().clear ();
        }

        // Reserving
        BOOST_UBLAS_INLINE
        void reserve (size_type non_zeros = 0) {
            non_zeros_ = std::max (non_zeros, size_type (1));
            non_zeros_ = std::min (non_zeros_, size_);
            detail::reserve (data (), non_zeros_);
        }

        // Proxy support
#ifdef BOOST_UBLAS_STRICT_VECTOR_SPARSE
        pointer find_element (size_type i) {
            iterator_type it (data ().find (i));
            if (it == data ().end () || (*it).first != i)
                return 0;
            return &(*it).second;
        }
#endif

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            const_iterator_type it (data ().find (i));
            if (it == data ().end () || (*it).first != i)
                return zero_;
            return (*it).second;
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
#ifndef BOOST_UBLAS_STRICT_VECTOR_SPARSE
            return data () [i];
#else
            return reference (*this, i);
#endif
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
        sparse_vector &operator = (const sparse_vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                size_ = v.size_;
                non_zeros_ = v.non_zeros_;
                data () = v.data ();
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_vector &assign_temporary (sparse_vector &v) {
            swap (v);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector &operator = (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (ae, non_zeros_));
#else
            // return assign (self_type (ae, non_zeros_));
            self_type temporary (ae, non_zeros_);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector &reset (const vector_expression<AE> &ae) {
            self_type temporary (ae, non_zeros_);
            resize (temporary.size (), non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector &assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector &operator += (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (*this + ae, non_zeros_));
#else
            // return assign (self_type (*this + ae, non_zeros_));
            self_type temporary (*this + ae, non_zeros_);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector &plus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector &operator -= (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (*this - ae, non_zeros_));
#else
            // return assign (self_type (*this - ae, non_zeros_));
            self_type temporary (*this - ae, non_zeros_);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector &minus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        sparse_vector &operator *= (const AT &at) {
            vector_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        sparse_vector &operator /= (const AT &at) {
            vector_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (sparse_vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                // BOOST_UBLAS_CHECK (non_zeros_ == v.non_zeros_, bad_size ());
                std::swap (size_, v.size_);
                std::swap (non_zeros_, v.non_zeros_);
                data ().swap (v.data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (sparse_vector &v1, sparse_vector &v2) {
            v1.swap (v2);
        }
#endif

        // Element insertion and erasure
        BOOST_UBLAS_INLINE
        void insert (size_type i, const_reference t) {
            BOOST_UBLAS_CHECK (data ().find (i) == data ().end (), bad_index ());
            data ().insert (data ().end (), std::pair<size_type, value_type> (i, t));
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i) {
            // FIXME: shouldn't we use const_iterator_type here?
            iterator_type it = data ().find (i);
            if (it == data ().end ())
                return;
            data ().erase (it);
        }
        BOOST_UBLAS_INLINE
        void clear () {
            data ().clear ();
        }

        class const_iterator;
        class iterator;

        // Element lookup
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
            return const_iterator (*this, data ().lower_bound (i));
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator find (size_type i) {
            return iterator (*this, data ().lower_bound (i));
        }

        // Iterators simply are pointers.

        class const_iterator:
            public container_const_reference<sparse_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename sparse_vector::difference_type difference_type;
            typedef typename sparse_vector::value_type value_type;
            typedef typename sparse_vector::const_reference reference;
            typedef typename sparse_vector::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &v, const const_iterator_type &it):
                container_const_reference<self_type> (v), it_ (it) {}
#ifndef BOOST_UBLAS_QUALIFIED_TYPENAME
            BOOST_UBLAS_INLINE
            const_iterator (const iterator &it):
                container_const_reference<self_type> (it ()), it_ (it.it_) {}
#else
            BOOST_UBLAS_INLINE
            const_iterator (const typename self_type::iterator &it):
                container_const_reference<self_type> (it ()), it_ (it.it_) {}
#endif
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
                return (*it_).second;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return (*it_).first;
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
                // FIXME: won't work with vector of vectors.
                // BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return &(*this) () == &it () && it_ == it.it_;
            }

        private:
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size_);
        }

        class iterator:
            public container_reference<sparse_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename sparse_vector::difference_type difference_type;
            typedef typename sparse_vector::value_type value_type;
            typedef typename sparse_vector::reference reference;
            typedef typename sparse_vector::pointer pointer;
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

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
#if ! defined (BOOST_UBLAS_STRICT_STORAGE_SPARSE) && ! defined (BOOST_UBLAS_STRICT_VECTOR_SPARSE)
                return (*it_).second;
#elif defined (BOOST_UBLAS_STRICT_VECTOR_SPARSE)
                return reference ((*this) (), &(*it_).second, index ());
#elif defined (BOOST_UBLAS_STRICT_STORAGE_SPARSE)
                return detail::make_reference ((*this) ().data (), it_);
#endif
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return (*it_).first;
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
                // FIXME: won't work with vector of vectors.
                // BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return &(*this) () == &it () && it_ == it.it_;
            }

        private:
            iterator_type it_;

            friend class const_iterator;
        };

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
        size_type non_zeros_;
        array_type data_;
        static value_type zero_;
    };

    template<class T, class A>
    typename sparse_vector<T, A>::value_type sparse_vector<T, A>::zero_ =
        sparse_vector<T, A>::value_type ();

    // Array based sparse vector class
    // Thanks to Kresimir Fresl for extending this to cover different index bases.
    template<class T, std::size_t IB, class IA, class TA>
    class compressed_vector:
        public vector_expression<compressed_vector<T, IB, IA, TA> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<compressed_vector<T, IB, IA, TA> >::operator ();
#endif
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
#ifndef BOOST_UBLAS_STRICT_VECTOR_SPARSE
        typedef T &reference;
#else
        typedef sparse_vector_element<compressed_vector<T, IB, IA, TA> > reference;
#endif
        typedef const T *const_pointer;
        typedef T *pointer;
        typedef IA index_array_type;
        typedef TA value_array_type;
        typedef const compressed_vector<T, IB, IA, TA> const_self_type;
        typedef compressed_vector<T, IB, IA, TA> self_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const_self_type> const_closure_type;
#else
        typedef const vector_reference<const_self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef typename IA::const_iterator const_iterator_type;
        typedef typename IA::iterator iterator_type;
        typedef sparse_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        compressed_vector ():
            vector_expression<self_type> (),
            size_ (0), non_zeros_ (0), filled_ (0),
            index_data_ (), value_data_ () {}
        BOOST_UBLAS_INLINE
        compressed_vector (size_type size, size_type non_zeros = 0):
            vector_expression<self_type> (),
            size_ (size), non_zeros_ (non_zeros), filled_ (0),
            index_data_ (non_zeros), value_data_ (non_zeros) {
            reserve (non_zeros_);
        }
        BOOST_UBLAS_INLINE
        compressed_vector (const compressed_vector &v):
            vector_expression<self_type> (),
            size_ (v.size_), non_zeros_ (v.non_zeros_), filled_ (v.filled_),
            index_data_ (v.index_data_), value_data_ (v.value_data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_vector (const vector_expression<AE> &ae, size_type non_zeros = 0):
            vector_expression<self_type> (),
            size_ (ae ().size ()), non_zeros_ (non_zeros), filled_ (0),
            index_data_ (non_zeros), value_data_ (non_zeros) {
            reserve (non_zeros_, false);
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
        }

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }
        BOOST_UBLAS_INLINE
        size_type non_zeros () const {
            return non_zeros_;
        }
        BOOST_UBLAS_INLINE
        size_type filled () const {
            return filled_;
        }
        BOOST_UBLAS_INLINE
        static size_type index_base () {
            return index_base_;
        }
        BOOST_UBLAS_INLINE
        const index_array_type &index_data () const {
            return index_data_;
        }
        BOOST_UBLAS_INLINE
        index_array_type &index_data () {
            return index_data_;
        }
        BOOST_UBLAS_INLINE
        const value_array_type &value_data () const {
            return value_data_;
        }
        BOOST_UBLAS_INLINE
        value_array_type &value_data () {
            return value_data_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, size_type non_zeros = 0, bool preserve = true) {
            size_ = size;
            non_zeros_ = std::max (non_zeros, size_type (1));
            non_zeros_ = std::min (non_zeros_, size_);
            filled_ = 0;
            detail::resize (index_data (), non_zeros_, preserve);
            detail::resize (value_data (), non_zeros_, preserve);
        }

        // Reserving
        BOOST_UBLAS_INLINE
        void reserve (size_type non_zeros = 0, bool preserve = true) {
            non_zeros_ = std::max (non_zeros, size_type (1));
            non_zeros_ = std::min (non_zeros_, size_);
            detail::resize (index_data (), non_zeros_, preserve);
            detail::resize (value_data (), non_zeros_, preserve);
        }

        // Proxy support
#ifdef BOOST_UBLAS_STRICT_VECTOR_SPARSE
        pointer find_element (size_type i) {
            iterator_type it (detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
            if (it == index_data ().begin () + filled_ || *it != k_based (i))
                return 0;
            return &value_data () [it - index_data ().begin ()];
        }
#endif

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            const_iterator_type it (detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
            if (it == index_data ().begin () + filled_ || *it != k_based (i))
                return zero_;
            return value_data () [it - index_data ().begin ()];
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
#ifndef BOOST_UBLAS_STRICT_VECTOR_SPARSE
            iterator_type it (detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
            if (it == index_data ().begin () + filled_ || *it != k_based (i)) {
                insert (i, value_type ());
                it = detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ());
            }
            return value_data () [it - index_data ().begin ()];
#else
            return reference (*this, i);
#endif
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
        compressed_vector &operator = (const compressed_vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                size_ = v.size_;
                non_zeros_ = v.non_zeros_;
                filled_ = v.filled_;
                detail::resize (index_data (), non_zeros_, false);
                detail::resize (value_data (), non_zeros_, false);
                index_data () = v.index_data ();
                value_data () = v.value_data ();
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        compressed_vector &assign_temporary (compressed_vector &v) { 
            swap (v);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_vector &operator = (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (ae, non_zeros_));
#else
            // return assign (self_type (ae, non_zeros_));
            self_type temporary (ae, non_zeros_);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_vector &reset (const vector_expression<AE> &ae) {
            self_type temporary (ae, non_zeros_);
            resize (temporary.size (), non_zeros_, false);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_vector &assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_vector &operator += (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (*this + ae, non_zeros_));
#else
            // return assign (self_type (*this + ae, non_zeros_));
            self_type temporary (*this + ae, non_zeros_);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_vector &plus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_vector &operator -= (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (*this - ae, non_zeros_));
#else
            // return assign (self_type (*this - ae, non_zeros_));
            self_type temporary (*this - ae, non_zeros_);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_vector &minus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        compressed_vector &operator *= (const AT &at) {
            vector_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        compressed_vector &operator /= (const AT &at) {
            vector_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (compressed_vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                // BOOST_UBLAS_CHECK (non_zeros_ == v.non_zeros_, bad_size ());
                std::swap (size_, v.size_);
                std::swap (non_zeros_, v.non_zeros_);
                std::swap (filled_, v.filled_);
                index_data ().swap (v.index_data ());
                value_data ().swap (v.value_data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (compressed_vector &v1, compressed_vector &v2) {
            v1.swap (v2);
        }
#endif

        // Element insertion and erasure
        BOOST_UBLAS_INLINE
        void push_back (size_type i, const_reference t) {
            if (filled_ >= non_zeros_)
                reserve (2 * non_zeros_);
            if (filled_ == 0 || index_data () [filled_ - 1] < k_based (i)) {
                ++ filled_;
                index_data () [filled_ - 1] = k_based (i);
                value_data () [filled_ - 1] = t;
                return;
            }
            // Raising exceptions abstracted as requested during review.
            // throw external_logic ();
            external_logic ().raise ();
        }
        BOOST_UBLAS_INLINE
        void insert (size_type i, const_reference t) {
            if (filled_ >= non_zeros_)
                reserve (2 * non_zeros_);
            iterator_type it (detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
            difference_type n = it - index_data ().begin ();
            BOOST_UBLAS_CHECK (filled_ == 0 || filled_ == size_type (n) || *it != k_based (i), external_logic ());
            ++ filled_;
            it = index_data ().begin () + n;
            std::copy_backward (it, index_data ().begin () + filled_ - 1, index_data ().begin () + filled_);
            *it = k_based (i);
            typename value_array_type::iterator itt (value_data ().begin () + n);
            std::copy_backward (itt, value_data ().begin () + filled_ - 1, value_data ().begin () + filled_);
            *itt = t;
        }
        BOOST_UBLAS_INLINE
        void pop_back () {
            BOOST_UBLAS_CHECK (filled_ > 0, external_logic ());
            -- filled_;
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i) {
            iterator_type it (detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
            difference_type n = it - index_data ().begin ();
            if (filled_ > size_type (n) && *it == k_based (i)) {
                std::copy (it + 1, index_data ().begin () + filled_, it);
                typename value_array_type::iterator itt (value_data ().begin () + n);
                std::copy (itt + 1, value_data ().begin () + filled_, itt);
                -- filled_;
            }
        }
        BOOST_UBLAS_INLINE
        void clear () {
            filled_ = 0;
        }

        class const_iterator;
        class iterator;

        // Element lookup
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
            return const_iterator (*this, detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator find (size_type i) {
            return iterator (*this, detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
        }

        // Iterators simply are pointers.

        class const_iterator:
            public container_const_reference<compressed_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename compressed_vector::difference_type difference_type;
            typedef typename compressed_vector::value_type value_type;
            typedef typename compressed_vector::const_reference reference;
            typedef typename compressed_vector::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &v, const const_iterator_type &it):
                container_const_reference<self_type> (v), it_ (it) {}
#ifndef BOOST_UBLAS_QUALIFIED_TYPENAME
            BOOST_UBLAS_INLINE
            const_iterator (const iterator &it):
                container_const_reference<self_type> (it ()), it_ (it.it_) {}
#else
            BOOST_UBLAS_INLINE
            const_iterator (const typename self_type::iterator &it):
                container_const_reference<self_type> (it ()), it_ (it.it_) {}
#endif
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
                return (*this) ().value_data () [it_ - (*this) ().index_data ().begin ()];
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return (*this) ().zero_based (*it_);
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
                // FIXME: won't work with vector of vectors.
                // BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return &(*this) () == &it () && it_ == it.it_;
            }

        private:
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size_);
        }

        class iterator:
            public container_reference<compressed_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename compressed_vector::difference_type difference_type;
            typedef typename compressed_vector::value_type value_type;
            typedef typename compressed_vector::reference reference;
            typedef typename compressed_vector::pointer pointer;
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

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
#if ! defined (BOOST_UBLAS_STRICT_VECTOR_SPARSE)
                return (*this) ().value_data () [it_ - (*this) ().index_data ().begin ()];
#else
                return reference ((*this) (), &(*this) ().value_data () [it_ - (*this) ().index_data ().begin ()], index ());
#endif
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return (*this) ().zero_based (*it_);
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
                // FIXME: won't work with vector of vectors.
                // BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return &(*this) () == &it () && it_ == it.it_;
            }

        private:
            iterator_type it_;

            friend class const_iterator;
        };

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
        BOOST_STATIC_CONSTANT (size_type, index_base_ = IB);
        size_type size_;
        size_type non_zeros_;
        size_type filled_;
        index_array_type index_data_;
        value_array_type value_data_;
        static value_type zero_;

        BOOST_UBLAS_INLINE
        static size_type zero_based (size_type k_based_index) {
            return k_based_index - index_base_;
        }
        BOOST_UBLAS_INLINE
        static size_type k_based (size_type zero_based_index) {
            return zero_based_index + index_base_;
        }

        friend class iterator;
        friend class const_iterator;
    };

    template<class T, std::size_t IB, class IA, class TA>
    typename compressed_vector<T, IB, IA, TA>::value_type compressed_vector<T, IB, IA, TA>::zero_ =
        compressed_vector<T, IB, IA, TA>::value_type ();

    // Array based sparse vector class
    // Thanks to Kresimir Fresl for extending this to cover different index bases.
    template<class T, std::size_t IB, class IA, class TA>
    class coordinate_vector:
        public vector_expression<coordinate_vector<T, IB, IA, TA> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<coordinate_vector<T, IB, IA, TA> >::operator ();
#endif
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
#ifndef BOOST_UBLAS_STRICT_VECTOR_SPARSE
        typedef T &reference;
#else
        typedef sparse_vector_element<coordinate_vector<T, IB, IA, TA> > reference;
#endif
        typedef const T *const_pointer;
        typedef T *pointer;
        typedef IA index_array_type;
        typedef TA value_array_type;
        typedef const coordinate_vector<T, IB, IA, TA> const_self_type;
        typedef coordinate_vector<T, IB, IA, TA> self_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const vector_const_reference<const_self_type> const_closure_type;
#else
        typedef const vector_reference<const_self_type> const_closure_type;
#endif
        typedef vector_reference<self_type> closure_type;
        typedef typename IA::const_iterator const_iterator_type;
        typedef typename IA::iterator iterator_type;
        typedef sparse_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        coordinate_vector ():
            vector_expression<self_type> (),
            size_ (0), non_zeros_ (0), filled_ (0),
            sorted_ (true), index_data_ (), value_data_ () {}
        BOOST_UBLAS_INLINE
        coordinate_vector (size_type size, size_type non_zeros = 0):
            vector_expression<self_type> (),
            size_ (size), non_zeros_ (non_zeros), filled_ (0),
            sorted_ (true), index_data_ (non_zeros), value_data_ (non_zeros) {
            reserve (non_zeros_);
        }
        BOOST_UBLAS_INLINE
        coordinate_vector (const coordinate_vector &v):
            vector_expression<self_type> (),
            size_ (v.size_), non_zeros_ (v.non_zeros_), filled_ (v.filled_),
            sorted_ (v.sorted_), index_data_ (v.index_data_), value_data_ (v.value_data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_vector (const vector_expression<AE> &ae, size_type non_zeros = 0):
            vector_expression<self_type> (),
            size_ (ae ().size ()), non_zeros_ (non_zeros), filled_ (0),
            sorted_ (true), index_data_ (non_zeros), value_data_ (non_zeros) {
            reserve (non_zeros_, false);
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
        }

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }
        BOOST_UBLAS_INLINE
        size_type non_zeros () const {
            return non_zeros_;
        }
        BOOST_UBLAS_INLINE
        size_type filled () const {
            return filled_;
        }
        BOOST_UBLAS_INLINE
        static size_type index_base () {
            return index_base_;
        }
        BOOST_UBLAS_INLINE
        const index_array_type &index_data () const {
            return index_data_;
        }
        BOOST_UBLAS_INLINE
        index_array_type &index_data () {
            return index_data_;
        }
        BOOST_UBLAS_INLINE
        const value_array_type &value_data () const {
            return value_data_;
        }
        BOOST_UBLAS_INLINE
        value_array_type &value_data () {
            return value_data_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, size_type non_zeros = 0, bool preserve = true) {
            size_ = size;
            non_zeros_ = std::max (non_zeros, size_type (1));
            // FIX: coordinate_vector may contain duplicate elements.
            // non_zeros_ = std::min (non_zeros_, size_);
            filled_ = 0;
            detail::resize (index_data (), non_zeros_, preserve);
            detail::resize (value_data (), non_zeros_, preserve);
        }

        // Reserving
        BOOST_UBLAS_INLINE
        void reserve (size_type non_zeros = 0, bool preserve = true) {
            non_zeros_ = std::max (non_zeros, size_type (1));
            // FIX: coordinate_vector may contain duplicate elements.
            // non_zeros_ = std::min (non_zeros_, size_);
            detail::resize (index_data (), non_zeros_, preserve);
            detail::resize (value_data (), non_zeros_, preserve);
        }

        // Proxy support
#ifdef BOOST_UBLAS_STRICT_VECTOR_SPARSE
        pointer find_element (size_type i) {
            sort ();
            iterator_type it (detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
            if (it == index_data ().begin () + filled_ || *it != k_based (i))
                return 0;
            return &value_data () [it - index_data ().begin ()];
        }
#endif

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            sort ();
            const_iterator_type it (detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
            if (it == index_data ().begin () + filled_ || *it != k_based (i))
                return zero_;
            return value_data () [it - index_data ().begin ()];
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
#ifndef BOOST_UBLAS_STRICT_VECTOR_SPARSE
            sort ();
            iterator_type it (detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
            if (it == index_data ().begin () + filled_ || *it != k_based (i)) {
                insert (i, value_type ());
                sort ();
                it = detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ());
            }
            return value_data () [it - index_data ().begin ()];
#else
            return reference (*this, i);
#endif
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
        coordinate_vector &operator = (const coordinate_vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                size_ = v.size_;
                non_zeros_ = v.non_zeros_;
                filled_ = v.filled_;
                sorted_ = v.sorted_;
                detail::resize (index_data (), non_zeros_, false);
                detail::resize (value_data (), non_zeros_, false);
                index_data () = v.index_data ();
                value_data () = v.value_data ();
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        coordinate_vector &assign_temporary (coordinate_vector &v) {
            swap (v);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_vector &operator = (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (ae, non_zeros_));
#else
            // return assign (self_type (ae, non_zeros_));
            self_type temporary (ae, non_zeros_);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_vector &reset (const vector_expression<AE> &ae) {
            self_type temporary (ae, non_zeros_);
            resize (temporary.size (), non_zeros_, false);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_vector &assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_vector &operator += (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (*this + ae, non_zeros_));
#else
            // return assign (self_type (*this + ae, non_zeros_));
            self_type temporary (*this + ae, non_zeros_);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_vector &plus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_vector &operator -= (const vector_expression<AE> &ae) {
#ifdef BOOST_UBLAS_MUTABLE_TEMPORARY
            return assign_temporary (self_type (*this - ae, non_zeros_));
#else
            // return assign (self_type (*this - ae, non_zeros_));
            self_type temporary (*this - ae, non_zeros_);
            return assign_temporary (temporary);
#endif
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_vector &minus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        coordinate_vector &operator *= (const AT &at) {
            vector_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        coordinate_vector &operator /= (const AT &at) {
            vector_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (coordinate_vector &v) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &v, external_logic ());
            if (this != &v) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == v.size_, bad_size ());
                // BOOST_UBLAS_CHECK (non_zeros_ == v.non_zeros_, bad_size ());
                std::swap (size_, v.size_);
                std::swap (non_zeros_, v.non_zeros_);
                std::swap (filled_, v.filled_);
                std::swap (sorted_, v.sorted_);
                index_data ().swap (v.index_data ());
                value_data ().swap (v.value_data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (coordinate_vector &v1, coordinate_vector &v2) {
            v1.swap (v2);
        }
#endif

        // Sorting
        BOOST_UBLAS_INLINE
        void sort () const {
            if (! sorted_ && filled_ > 0) {
                index_pair_array<index_array_type, value_array_type>
                    ipa (filled_, index_data_, value_data_);
                std::sort (ipa.begin (), ipa.end ());
                // FIX: check for duplicates
                size_type filled = 1;
                for (size_type i = 1; i < filled_; ++ i) {
                    if (index_data_ [filled - 1] != index_data_ [i]) {
                        ++ filled;
                        if (filled - 1 != i) {
                            index_data_ [filled - 1] = index_data_ [i];
                            value_data_ [filled - 1] = value_data_ [i];
                        }
                    } else {
                        value_data_ [filled - 1] += value_data_ [i];
                    }
                }
                filled_ = filled;
                sorted_ = true;
            }
        }

        // Element insertion and erasure
        BOOST_UBLAS_INLINE
        void push_back (size_type i, const_reference t) {
            if (filled_ >= non_zeros_)
                reserve (2 * non_zeros_);
            if (filled_ == 0 || index_data () [filled_ - 1] < k_based (i)) {
                ++ filled_;
                index_data () [filled_ - 1] = k_based (i);
                value_data () [filled_ - 1] = t;
                return;
            }
            // Raising exceptions abstracted as requested during review.
            // throw external_logic ();
            external_logic ().raise ();
        }
        BOOST_UBLAS_INLINE
        void insert (size_type i, const_reference t) {
            if (filled_ >= non_zeros_)
                reserve (2 * non_zeros_);
            ++ filled_;
            index_data () [filled_ - 1] = k_based (i);
            value_data () [filled_ - 1] = t;
            sorted_ = false;
        }
        BOOST_UBLAS_INLINE
        void pop_back () {
            BOOST_UBLAS_CHECK (filled_ > 0, external_logic ());
            -- filled_;
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i) {
            sort ();
            iterator_type it (detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
            difference_type n = it - index_data ().begin ();
            if (filled_ > size_type (n) && *it == k_based (i)) {
                std::copy (it + 1, index_data ().begin () + filled_, it);
                typename value_array_type::iterator itt (value_data ().begin () + n);
                std::copy (itt + 1, value_data ().begin () + filled_, itt);
                -- filled_;
            }
        }
        BOOST_UBLAS_INLINE
        void clear () {
            filled_ = 0;
        }

        class const_iterator;
        class iterator;

        // Element lookup
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
            sort ();
            return const_iterator (*this, detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator find (size_type i) {
            sort ();
            return iterator (*this, detail::lower_bound (index_data ().begin (), index_data ().begin () + filled_, k_based (i), std::less<size_type> ()));
        }

        // Iterators simply are pointers.

        class const_iterator:
            public container_const_reference<coordinate_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename coordinate_vector::difference_type difference_type;
            typedef typename coordinate_vector::value_type value_type;
            typedef typename coordinate_vector::const_reference reference;
            typedef typename coordinate_vector::const_pointer pointer;
#endif

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &v, const const_iterator_type &it):
                container_const_reference<self_type> (v), it_ (it) {}
#ifndef BOOST_UBLAS_QUALIFIED_TYPENAME
            BOOST_UBLAS_INLINE
            const_iterator (const iterator &it):
                container_const_reference<self_type> (it ()), it_ (it.it_) {}
#else
            BOOST_UBLAS_INLINE
            const_iterator (const typename self_type::iterator &it):
                container_const_reference<self_type> (it ()), it_ (it.it_) {}
#endif
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
                return (*this) ().value_data () [it_ - (*this) ().index_data ().begin ()];
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return (*this) ().zero_based (*it_);
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
                // FIXME: won't work with vector of vectors.
                // BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return &(*this) () == &it () && it_ == it.it_;
            }

        private:
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size_);
        }

        class iterator:
            public container_reference<coordinate_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename coordinate_vector::difference_type difference_type;
            typedef typename coordinate_vector::value_type value_type;
            typedef typename coordinate_vector::reference reference;
            typedef typename coordinate_vector::pointer pointer;
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

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
#if ! defined (BOOST_UBLAS_STRICT_VECTOR_SPARSE)
                return (*this) ().value_data () [it_ - (*this) ().index_data ().begin ()];
#else
                return reference ((*this) (), &(*this) ().value_data () [it_ - (*this) ().index_data ().begin ()], index ());
#endif
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return (*this) ().zero_based (*it_);
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
                // FIXME: won't work with vector of vectors.
                // BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return &(*this) () == &it () && it_ == it.it_;
            }

        private:
            iterator_type it_;

            friend class const_iterator;
        };

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
        BOOST_STATIC_CONSTANT (size_type, index_base_ = IB);
        size_type size_;
        size_type non_zeros_;
        mutable size_type filled_;
        mutable bool sorted_;
        mutable index_array_type index_data_;
        mutable value_array_type value_data_;
        static value_type zero_;

        BOOST_UBLAS_INLINE
        static size_type zero_based (size_type k_based_index) {
            return k_based_index - index_base_;
        }
        BOOST_UBLAS_INLINE
        static size_type k_based (size_type zero_based_index) {
            return zero_based_index + index_base_;
        }

        friend class iterator;
        friend class const_iterator;
    };

    template<class T, std::size_t IB, class IA, class TA>
    typename coordinate_vector<T, IB, IA, TA>::value_type coordinate_vector<T, IB, IA, TA>::zero_ =
        coordinate_vector<T, IB, IA, TA>::value_type ();

}}}

#endif




