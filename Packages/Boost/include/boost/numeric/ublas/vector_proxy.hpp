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

#ifndef BOOST_UBLAS_VECTOR_PROXY_H
#define BOOST_UBLAS_VECTOR_PROXY_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>

// Iterators based on ideas of Jeremy Siek

namespace boost { namespace numeric { namespace ublas {

    // Vector based range class
    template<class V>
    class vector_range:
        public vector_expression<vector_range<V> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<vector_range<V> >::operator ();
#endif
        typedef const V const_vector_type;
        typedef V vector_type;
        typedef typename V::simd_category simd_category;
        typedef typename V::size_type size_type;
        typedef typename V::difference_type difference_type;
        typedef typename V::value_type value_type;
#ifndef BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
        typedef typename V::const_reference const_reference;
        typedef typename V::reference reference;
#else
        typedef typename V::const_reference const_reference;
        typedef typename boost::mpl::if_<boost::is_const<V>,
                                          typename V::const_reference,
                                          typename V::reference>::type reference;
#endif
#ifndef BOOST_UBLAS_CT_PROXY_CLOSURE_TYPEDEFS
        typedef typename V::closure_type vector_closure_type;
#else
        typedef typename boost::mpl::if_<boost::is_const<V>,
                                          typename V::const_closure_type,
                                          typename V::closure_type>::type vector_closure_type;
#endif
    private:
        typedef vector_range<vector_type> self_type;
    public:
        typedef basic_range<size_type, difference_type> range_type;
        typedef const self_type const_closure_type;
        typedef self_type closure_type;
        typedef typename V::vector_temporary_type vector_temporary_type;
        typedef typename storage_restrict_traits<typename V::storage_category,
                                                 dense_proxy_tag>::storage_category storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        vector_range ():
            data_ (nil_), r_ () {}
        BOOST_UBLAS_INLINE
        vector_range (vector_type &data, const range_type &r):
            data_ (data), r_ (r.preprocess (data.size ())) {
            // Early checking of preconditions here.
            // BOOST_UBLAS_CHECK (r_.start () <= data_.size () &&
            //                   r_.start () + r_.size () <= data_.size (), bad_index ());
        }
        BOOST_UBLAS_INLINE
        vector_range (const vector_closure_type &data, const range_type &r, bool):
            data_ (data), r_ (r.preprocess (data.size ())) {
            // Early checking of preconditions here.
            // BOOST_UBLAS_CHECK (r_.start () <= data_.size () &&
            //                    r_.start () + r_.size () <= data_.size (), bad_index ());
        }

        // Accessors
        BOOST_UBLAS_INLINE
        size_type start () const {
            return r_.start ();
        }
        BOOST_UBLAS_INLINE
        size_type size () const {
            return r_.size ();
        }
        BOOST_UBLAS_INLINE
        const vector_closure_type &data () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        vector_closure_type &data () {
            return data_;
        }

        // Element access
#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            return data_ (r_ (i));
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) {
            return data_ (r_ (i));
        }

        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            return (*this) (i);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) {
            return (*this) (i);
        }
#else
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) const {
            return data_ (r_ (i));
        }

        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) const {
            return (*this) (i);
        }
#endif

        // ISSUE can this be done in free project function?
        // Although a const function can create a non-const proxy to a non-const object
        // Critical is that vector_type and data_ (vector_closure_type) are const correct
        BOOST_UBLAS_INLINE
        vector_range<vector_type> project (const range_type &r) const {
            return vector_range<vector_type> (data_, r_.compose (r.preprocess (data_.size ())), false);
        }

        // Assignment
        BOOST_UBLAS_INLINE
        vector_range &operator = (const vector_range &vr) {
            // FIXME: the ranges could be differently sized.
            // std::copy (vr.begin (), vr.end (), begin ());
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (vr));
            return *this;
        }
        BOOST_UBLAS_INLINE
        vector_range &assign_temporary (vector_range &vr) {
            // FIXME: this is suboptimal.
            // return *this = vr;
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, value_type> (), *this, vr);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_range &operator = (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_range &assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_range &operator += (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (*this + ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_range &plus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_plus_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_range &operator -= (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (*this - ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_range &minus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_minus_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        vector_range &operator *= (const AT &at) {
            vector_assign_scalar (scalar_multiplies_assign<BOOST_UBLAS_TYPENAME iterator::reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        vector_range &operator /= (const AT &at) {
            vector_assign_scalar (scalar_divides_assign<BOOST_UBLAS_TYPENAME iterator::reference, AT> (), *this, at);
            return *this;
        }

        // Closure comparison
        BOOST_UBLAS_INLINE
        bool same_closure (const vector_range &vr) const {
            return (*this).data_.same_closure (vr.data_);
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const vector_range &vr) const {
            return (*this).data_ == vr.data_ && r_ == vr.r_;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (vector_range vr) {
            if (this != &vr) {
                BOOST_UBLAS_CHECK (size () == vr.size (), bad_size ());
                // Sparse ranges may be nonconformant now.
                // std::swap_ranges (begin (), end (), vr.begin ());
                vector_swap (scalar_swap<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME iterator::reference> (), *this, vr);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (vector_range vr1, vector_range vr2) {
            vr1.swap (vr2);
        }
#endif

        // Iterator types
    private:
#ifndef BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
        // Reuse the vector iterator
        typedef typename V::const_iterator const_iterator_type;
        typedef typename V::iterator iterator_type;
#else
        typedef typename V::const_iterator const_iterator_type;
        typedef typename boost::mpl::if_<boost::is_const<V>,
                                          typename V::const_iterator,
                                          typename V::iterator>::type iterator_type;
#endif

    public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_iterator<vector_range<vector_type>,
                                 BOOST_UBLAS_TYPENAME iterator_type::iterator_category> iterator;
        typedef indexed_const_iterator<vector_range<vector_type>,
                                       BOOST_UBLAS_TYPENAME const_iterator_type::iterator_category> const_iterator;
#else
        class const_iterator;
        class iterator;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
            const_iterator_type it (data_.find (start () + i));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator (*this, it.index ());
#else
            return const_iterator (*this, it);
#endif
        }
        BOOST_UBLAS_INLINE
        iterator find (size_type i) {
            iterator_type it (data_.find (start () + i));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return iterator (*this, it.index ());
#else
            return iterator (*this, it);
#endif
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator:
            public container_const_reference<vector_range>,
#ifndef BOOST_UBLAS_NO_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename const_iterator_type::iterator_category>::template
                        iterator_base<const_iterator, value_type>::type {
#else
            public random_access_iterator_base<typename V::const_iterator::iterator_category,
                                               const_iterator, value_type> {
#endif
        public:
            typedef typename const_iterator_type::iterator_category iterator_category;
            typedef typename const_iterator_type::difference_type difference_type;
            typedef typename const_iterator_type::value_type value_type;
            typedef typename const_iterator_type::reference reference;
            typedef typename const_iterator_type::pointer pointer;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &vr, const const_iterator_type &it):
                container_const_reference<self_type> (vr), it_ (it) {}
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return *it_;
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_.index () - (*this) ().start ();
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator_type it_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size ());
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator:
            public container_reference<vector_range>,
#ifndef BOOST_UBLAS_NO_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename iterator_type::iterator_category>::template
                        iterator_base<iterator, value_type>::type {
#else
            public random_access_iterator_base<typename V::iterator::iterator_category,
                                               iterator, value_type> {
#endif
        public:
            typedef typename iterator_type::iterator_category iterator_category;
            typedef typename iterator_type::difference_type difference_type;
            typedef typename iterator_type::value_type value_type;
            typedef typename iterator_type::reference reference;
            typedef typename iterator_type::pointer pointer;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator ():
                container_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator (self_type &vr, const iterator_type &it):
                container_reference<self_type> (vr), it_ (it) {}

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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
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
                return it_.index () - (*this) ().start ();
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
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
            return find (size ());
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
        vector_closure_type data_;
        range_type r_;
        static vector_type nil_;
    };

    template<class V>
    typename vector_range<V>::vector_type vector_range<V>::nil_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME vector_range<V>::vector_type ()
#endif
    ;

    // Projections
    template<class V>
    BOOST_UBLAS_INLINE
    vector_range<V> project (V &data, const typename vector_range<V>::range_type &r) {
        return vector_range<V> (data, r);
    }
#ifdef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
    template<class V>
    BOOST_UBLAS_INLINE
    const vector_range<const V> project_const (const V &data, const typename vector_range<V>::range_type &r) {
        return vector_range<const V> (data, r);
    }
#else
    template<class V>
    BOOST_UBLAS_INLINE
    const vector_range<const V> project (const V &data, const typename vector_range<V>::range_type &r) {
        // ISSUE was: return vector_range<V> (const_cast<V &> (data), r);
        return vector_range<const V> (data, r);
    }
    template<class V>
    BOOST_UBLAS_INLINE
    vector_range<V> project (vector_range<V> &data, const typename vector_range<V>::range_type &r) {
        return data.project (r);
    }
    template<class V>
    BOOST_UBLAS_INLINE
    const vector_range<V> project (const vector_range<V> &data, const typename vector_range<V>::range_type &r) {
        return data.project (r);
    }
#endif

    // Vector based slice class
    template<class V>
    class vector_slice:
        public vector_expression<vector_slice<V> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<vector_slice<V> >::operator ();
#endif
        typedef const V const_vector_type;
        typedef V vector_type;
        typedef typename V::simd_category simd_category;
        typedef typename V::size_type size_type;
        typedef typename V::difference_type difference_type;
        typedef typename V::value_type value_type;
#ifndef BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
        typedef typename V::const_reference const_reference;
        typedef typename V::reference reference;
#else
        typedef typename V::const_reference const_reference;
        typedef typename boost::mpl::if_<boost::is_const<V>,
                                          typename V::const_reference,
                                          typename V::reference>::type reference;
#endif
#ifndef BOOST_UBLAS_CT_PROXY_CLOSURE_TYPEDEFS
        typedef typename V::closure_type vector_closure_type;
#else
        typedef typename boost::mpl::if_<boost::is_const<V>,
                                          typename V::const_closure_type,
                                          typename V::closure_type>::type vector_closure_type;
#endif
    private:
        typedef vector_slice<vector_type> self_type;
    public:
        typedef basic_range<size_type, difference_type> range_type;
        typedef basic_slice<size_type, difference_type> slice_type;
        typedef const self_type const_closure_type;
        typedef self_type closure_type;
        typedef typename V::vector_temporary_type vector_temporary_type;
        typedef typename storage_restrict_traits<typename V::storage_category,
                                                 dense_proxy_tag>::storage_category storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        vector_slice ():
            data_ (nil_), s_ () {}
        BOOST_UBLAS_INLINE
        vector_slice (vector_type &data, const slice_type &s):
            data_ (data), s_ (s.preprocess (data.size ())) {
            // Early checking of preconditions here.
            // BOOST_UBLAS_CHECK (s_.start () <= data_.size () &&
            //                    s_.start () + s_.stride () * (s_.size () - (s_.size () > 0)) <= data_.size (), bad_index ());
        }
        BOOST_UBLAS_INLINE
        vector_slice (const vector_closure_type &data, const slice_type &s, int):
            data_ (data), s_ (s.preprocess (data.size ())) {
            // Early checking of preconditions here.
            // BOOST_UBLAS_CHECK (s_.start () <= data_.size () &&
            //                    s_.start () + s_.stride () * (s_.size () - (s_.size () > 0)) <= data_.size (), bad_index ());
        }

        // Accessors
        BOOST_UBLAS_INLINE
        size_type start () const {
            return s_.start ();
        }
        BOOST_UBLAS_INLINE
        difference_type stride () const {
            return s_.stride ();
        }
        BOOST_UBLAS_INLINE
        size_type size () const {
            return s_.size ();
        }
        BOOST_UBLAS_INLINE
        const vector_closure_type &data () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        vector_closure_type &data () {
            return data_;
        }

        // Element access
#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            return data_ (s_ (i));
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) {
            return data_ (s_ (i));
        }

        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            return (*this) (i);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) {
            return (*this) (i);
        }
#else
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) const {
            return data_ (s_ (i));
        }

        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) const {
            return (*this) (i);
        }
#endif

        // ISSUE can this be done in free project function?
        // Although a const function can create a non-const proxy to a non-const object
        // Critical is that vector_type and data_ (vector_closure_type) are const correct
        BOOST_UBLAS_INLINE
        vector_slice<vector_type> project (const range_type &r) const {
            return vector_slice<vector_type>  (data_, s_.compose (r.preprocess (data_.size ())), false);
        }
        BOOST_UBLAS_INLINE
        vector_slice<vector_type> project (const slice_type &s) const {
            return vector_slice<vector_type>  (data_, s_.compose (s.preprocess (data_.size ())), false);
        }

        // Assignment
        BOOST_UBLAS_INLINE
        vector_slice &operator = (const vector_slice &vs) {
            // FIXME: the slices could be differently sized.
            // std::copy (vs.begin (), vs.end (), begin ());
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (vs));
            return *this;
        }
        BOOST_UBLAS_INLINE
        vector_slice &assign_temporary (vector_slice &vs) {
            // FIXME: this is suboptimal.
            // return *this = vs;
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, value_type> (), *this, vs);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_slice &operator = (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_slice &assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_slice &operator += (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (*this + ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_slice &plus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_plus_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_slice &operator -= (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (*this - ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_slice &minus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_minus_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        vector_slice &operator *= (const AT &at) {
            vector_assign_scalar (scalar_multiplies_assign<BOOST_UBLAS_TYPENAME iterator::reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        vector_slice &operator /= (const AT &at) {
            vector_assign_scalar (scalar_divides_assign<BOOST_UBLAS_TYPENAME iterator::reference, AT> (), *this, at);
            return *this;
        }

        // Closure comparison
        BOOST_UBLAS_INLINE
        bool same_closure (const vector_slice &vr) const {
            return (*this).data_.same_closure (vr.data_);
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const vector_slice &vs) const {
            return (*this).data_ == vs.data_ && s_ == vs.s_;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (vector_slice vs) {
            if (this != &vs) {
                BOOST_UBLAS_CHECK (size () == vs.size (), bad_size ());
                // Sparse ranges may be nonconformant now.
                // std::swap_ranges (begin (), end (), vs.begin ());
                vector_swap (scalar_swap<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME iterator::reference> (), *this, vs);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (vector_slice vs1, vector_slice vs2) {
            vs1.swap (vs2);
        }
#endif

        // Iterator types
    private:
        // Use slice as an index - FIXME this fails for packed assignment
        typedef BOOST_UBLAS_TYPENAME slice_type::const_iterator const_iterator_type;
        typedef BOOST_UBLAS_TYPENAME slice_type::const_iterator iterator_type;

    public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_iterator<vector_slice<vector_type>,
                                 BOOST_UBLAS_TYPENAME vector_type::iterator::iterator_category> iterator;
        typedef indexed_const_iterator<vector_slice<vector_type>,
                                       BOOST_UBLAS_TYPENAME vector_type::const_iterator::iterator_category> const_iterator;
#else
        class const_iterator;
        class iterator;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator (*this, i);
#else
            return const_iterator (*this, s_.begin () + i);
#endif
        }
        BOOST_UBLAS_INLINE
        iterator find (size_type i) {
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return iterator (*this, i);
#else
            return iterator (*this, s_.begin () + i);
#endif
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator:
            public container_const_reference<vector_slice>,
#ifndef BOOST_UBLAS_NO_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename V::const_iterator::iterator_category>::template
                        iterator_base<const_iterator, value_type>::type {
#else
            public random_access_iterator_base<typename V::const_iterator::iterator_category,
                                               const_iterator, value_type> {
#endif
        public:
            typedef typename V::const_iterator::iterator_category iterator_category;
            typedef typename V::const_iterator::difference_type difference_type;
            typedef typename V::const_iterator::value_type value_type;
            typedef typename V::const_reference reference;    //FIXME due to indexing access
            typedef typename V::const_iterator::pointer pointer;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &vs, const const_iterator_type &it):
                container_const_reference<self_type> (vs), it_ (it) {}
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                // FIXME replace find with at_element
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return (*this) ().data_ (*it_);
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_.index ();
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator_type it_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size ());
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator:
            public container_reference<vector_slice>,
#ifndef BOOST_UBLAS_NO_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename V::iterator::iterator_category>::template
                        iterator_base<iterator, value_type>::type {
#else
            public random_access_iterator_base<typename V::iterator::iterator_category,
                                               iterator, value_type> {
#endif
        public:
            typedef typename V::iterator::iterator_category iterator_category;
            typedef typename V::iterator::difference_type difference_type;
            typedef typename V::iterator::value_type value_type;
            typedef typename V::reference reference;    //FIXME due to indexing access
            typedef typename V::iterator::pointer pointer;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator ():
                container_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator (self_type &vs, const iterator_type &it):
                container_reference<self_type> (vs), it_ (it) {}

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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                // FIXME replace find with at_element
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return (*this) ().data_ (*it_);
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_.index ();
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
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
            return find (size ());
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

#ifndef BOOST_UBLAS_NESTED_CLASS_DR45
    private:
#endif
        vector_closure_type data_;
        slice_type s_;
        static vector_type nil_;
    };

    template<class V>
    typename vector_slice<V>::vector_type vector_slice<V>::nil_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME vector_slice<V>::vector_type ()
#endif
    ;
    
    // Projections
    template<class V>
    BOOST_UBLAS_INLINE
    vector_slice<V> project (V &data, const typename vector_slice<V>::slice_type &s) {
        return vector_slice<V> (data, s);
    }
#ifdef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
    template<class V>
    BOOST_UBLAS_INLINE
    const vector_slice<const V> project_const (const V &data, const typename vector_slice<V>::slice_type &s) {
        return vector_slice<const V> (data, s);
    }
#else
    template<class V>
    BOOST_UBLAS_INLINE
    const vector_slice<const V> project (const V &data, const typename vector_slice<V>::slice_type &s) {
        // ISSUE was: return vector_slice<V> (const_cast<V &> (data), s);
        return vector_slice<const V> (data, s);
    }
    template<class V>
    BOOST_UBLAS_INLINE
    vector_slice<V> project (vector_slice<V> &data, const typename vector_slice<V>::slice_type &s) {
        return data.project (s);
    }
    template<class V>
    BOOST_UBLAS_INLINE
    const vector_slice<V> project (const vector_slice<V> &data, const typename vector_slice<V>::slice_type &s) {
        return data.project (s);
    }
#ifndef BOOST_UBLAS_MSVC71_FUNCTION_TEMPLATE_ORDERING
    template<class V>
    BOOST_UBLAS_INLINE
    vector_slice<V> project (vector_slice<V> &data, const typename vector_slice<V>::range_type &r) {
        return data.project (r);
    }
    template<class V>
    BOOST_UBLAS_INLINE
    const vector_slice<V> project (const vector_slice<V> &data, const typename vector_slice<V>::range_type &r) {
        return data.project (r);
    }
#endif
#endif

    // Vector based indirection class
    // Contributed by Toon Knapen.
    // Extended and optimized by Kresimir Fresl.
    template<class V, class IA>
    class vector_indirect:
        public vector_expression<vector_indirect<V, IA> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<vector_indirect<V, IA> >::operator ();
#endif
        typedef const V const_vector_type;
        typedef V vector_type;
        typedef const IA const_indirect_array_type;
        typedef IA indirect_array_type;
        typedef typename V::size_type size_type;
        typedef typename V::difference_type difference_type;
        typedef typename V::value_type value_type;
#ifndef BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
        typedef typename V::const_reference const_reference;
        typedef typename V::reference reference;
#else
        typedef typename V::const_reference const_reference;
        typedef typename boost::mpl::if_<boost::is_const<V>,
                                          typename V::const_reference,
                                          typename V::reference>::type reference;
#endif
#ifndef BOOST_UBLAS_CT_PROXY_CLOSURE_TYPEDEFS
        typedef typename V::closure_type vector_closure_type;
#else
        typedef typename boost::mpl::if_<boost::is_const<V>,
                                          typename V::const_closure_type,
                                          typename V::closure_type>::type vector_closure_type;
#endif
    private:
        typedef vector_indirect<vector_type, indirect_array_type> self_type;
    public:
        typedef basic_range<size_type, difference_type> range_type;
        typedef basic_slice<size_type, difference_type> slice_type;
        typedef const self_type const_closure_type;
        typedef self_type closure_type;
        typedef typename V::vector_temporary_type vector_temporary_type;
        typedef typename storage_restrict_traits<typename V::storage_category,
                                                 dense_proxy_tag>::storage_category storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        vector_indirect ():
            data_ (nil_), ia_ () {}
        BOOST_UBLAS_INLINE
        vector_indirect (vector_type &data, size_type size):
            data_ (data), ia_ (size) {}
        BOOST_UBLAS_INLINE
        vector_indirect (vector_type &data, const indirect_array_type &ia):
            data_ (data), ia_ (ia.preprocess (data.size ())) {}
        BOOST_UBLAS_INLINE
        vector_indirect (const vector_closure_type &data, const indirect_array_type &ia, int):
            data_ (data), ia_ (ia.preprocess (data.size ())) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const {
            return ia_.size ();
        }
        BOOST_UBLAS_INLINE
        const vector_closure_type &data () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        vector_closure_type &data () {
            return data_;
        }
        BOOST_UBLAS_INLINE
        const_indirect_array_type &indirect () const {
            return ia_;
        }
        BOOST_UBLAS_INLINE
        indirect_array_type &indirect () {
            return ia_;
        }

        // Element access
#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            return data_ (ia_ (i));
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) {
            return data_ (ia_ (i));
        }

        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            return (*this) (i);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) {
            return (*this) (i);
        }
#else
        BOOST_UBLAS_INLINE
        reference operator () (size_type i) const {
            return data_ (ia_ (i));
        }

        BOOST_UBLAS_INLINE
        reference operator [] (size_type i) const {
            return (*this) (i);
        }
#endif

        // ISSUE can this be done in free project function?
        // Although a const function can create a non-const proxy to a non-const object
        // Critical is that vector_type and data_ (vector_closure_type) are const correct
        BOOST_UBLAS_INLINE
        vector_indirect<vector_type, indirect_array_type> project (const range_type &r) const {
            return vector_indirect<vector_type, indirect_array_type> (data_, ia_.compose (r.preprocess (data_.size ())), 0);
        }
        BOOST_UBLAS_INLINE
        vector_indirect<vector_type, indirect_array_type> project (const slice_type &s) const {
            return vector_indirect<vector_type, indirect_array_type> (data_, ia_.compose (s.preprocess (data_.size ())), 0);
        }
        BOOST_UBLAS_INLINE
        vector_indirect<vector_type, indirect_array_type> project (const indirect_array_type &ia) const {
            return vector_indirect<vector_type, indirect_array_type> (data_, ia_.compose (ia.preprocess (data_.size ())), 0);
        }

        // Assignment
        BOOST_UBLAS_INLINE
        vector_indirect &operator = (const vector_indirect &vi) {
            // FIXME: the indirect_arrays could be differently sized.
            // std::copy (vi.begin (), vi.end (), begin ());
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (vi));
            return *this;
        }
        BOOST_UBLAS_INLINE
        vector_indirect &assign_temporary (vector_indirect &vi) {
            // FIXME: this is suboptimal.
            // return *this = vi;
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, value_type> (), *this, vi);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_indirect &operator = (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_indirect &assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_indirect &operator += (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (*this + ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_indirect &plus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_plus_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_indirect &operator -= (const vector_expression<AE> &ae) {
            vector_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME vector_temporary_type::value_type> (), *this, vector_temporary_type (*this - ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        vector_indirect &minus_assign (const vector_expression<AE> &ae) {
            vector_assign (scalar_minus_assign<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        vector_indirect &operator *= (const AT &at) {
            vector_assign_scalar (scalar_multiplies_assign<BOOST_UBLAS_TYPENAME iterator::reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        vector_indirect &operator /= (const AT &at) {
            vector_assign_scalar (scalar_divides_assign<BOOST_UBLAS_TYPENAME iterator::reference, AT> (), *this, at);
            return *this;
        }

        // Closure comparison
        BOOST_UBLAS_INLINE
        bool same_closure (const vector_indirect &vr) const {
return true;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const vector_indirect &vi) const {
            return (*this).data_ == vi.data_ && ia_ == vi.ia_;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (vector_indirect vi) {
            if (this != &vi) {
                BOOST_UBLAS_CHECK (size () == vi.size (), bad_size ());
                // Sparse ranges may be nonconformant now.
                // std::swap_ranges (begin (), end (), vi.begin ());
                vector_swap (scalar_swap<BOOST_UBLAS_TYPENAME iterator::reference, BOOST_UBLAS_TYPENAME iterator::reference> (), *this, vi);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (vector_indirect vi1, vector_indirect vi2) {
            vi1.swap (vi2);
        }
#endif

        // Iterator types
    private:
        // Use indirect array as an index - FIXME this fails for packed assignment
        typedef typename IA::const_iterator const_iterator_type;
        typedef typename IA::const_iterator iterator_type;

    public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_iterator<vector_indirect<vector_type, indirect_array_type>,
                                 BOOST_UBLAS_TYPENAME vector_type::iterator::iterator_category> iterator;
        typedef indexed_const_iterator<vector_indirect<vector_type, indirect_array_type>,
                                       BOOST_UBLAS_TYPENAME vector_type::const_iterator::iterator_category> const_iterator;
#else
        class const_iterator;
        class iterator;
#endif
        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator (*this, i);
#else
            return const_iterator (*this, ia_.begin () + i);
#endif
        }
        BOOST_UBLAS_INLINE
        iterator find (size_type i) {
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return iterator (*this, i);
#else
            return iterator (*this, ia_.begin () + i);
#endif
        }

        // Iterators simply are indices.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator:
            public container_const_reference<vector_indirect>,
#ifndef BOOST_UBLAS_NO_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename V::const_iterator::iterator_category>::template
                        iterator_base<const_iterator, value_type>::type {
#else
            public random_access_iterator_base<typename V::const_iterator::iterator_category,
                                               const_iterator, value_type> {
#endif
        public:
            typedef typename V::const_iterator::iterator_category iterator_category;
            typedef typename V::const_iterator::difference_type difference_type;
            typedef typename V::const_iterator::value_type value_type;
            typedef typename V::const_reference reference;    //FIXME due to indexing access
            typedef typename V::const_iterator::pointer pointer;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &vi, const const_iterator_type &it):
                container_const_reference<self_type> (vi), it_ (it) {}
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                // FIXME replace find with at_element
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return (*this) ().data_ (*it_);
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_.index ();
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator_type it_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return find (0);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return find (size ());
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator:
            public container_reference<vector_indirect>,
#ifndef BOOST_UBLAS_NO_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename V::iterator::iterator_category>::template
                        iterator_base<iterator, value_type>::type {
#else
            public random_access_iterator_base<typename V::iterator::iterator_category,
                                               iterator, value_type> {
#endif
        public:
            typedef typename V::iterator::iterator_category iterator_category;
            typedef typename V::iterator::difference_type difference_type;
            typedef typename V::iterator::value_type value_type;
            typedef typename V::reference reference;    //FIXME due to indexing access
            typedef typename V::iterator::pointer pointer;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator ():
                container_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator (self_type &vi, const iterator_type &it):
                container_reference<self_type> (vi), it_ (it) {}

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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                // FIXME replace find with at_element
                BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
                return (*this) ().data_ (*it_);
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it_.index ();
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
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
            return find (size ());
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

#ifndef BOOST_UBLAS_NESTED_CLASS_DR45
    private:
#endif
        vector_closure_type data_;
        indirect_array_type ia_;
        static vector_type nil_;
    };

    template<class V, class IA>
    typename vector_indirect<V, IA>::vector_type vector_indirect<V, IA>::nil_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME vector_indirect<V, IA>::vector_type ()
#endif
    ;

    // Projections

    // These signatures are too general for MSVC
    // template<class V, class IA>
    // BOOST_UBLAS_INLINE
    // vector_indirect<V, IA> project (V &data, const IA &ia) {
    //     return vector_indirect<V, IA> (data, ia);
    // }
#ifdef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
    // template<class V, class IA>
    // BOOST_UBLAS_INLINE
    // const vector_indirect<const V, IA> project_const (const V &data, const IA &ia) {
    //     return vector_indirect<const V, IA> (data, ia);
    // }
#else
    // template<class V, class IA>
    // BOOST_UBLAS_INLINE
    // const vector_indirect<const V, IA> project (const V &data, const IA &ia) {
    //     // ISSUE was: return vector_indirect<V, IA> (const_cast<V &> (data), ia)
    //     return vector_indirect<const V, IA> (data, ia);
    // }
    // template<class V, class IA>
    // BOOST_UBLAS_INLINE
    // vector_indirect<V, IA> project (vector_indirect<V, IA> &data, const IA &ia) {
    //     return data.project (ia);
    // }
    // template<class V, class IA>
    // BOOST_UBLAS_INLINE
    // const vector_indirect<V, IA> project (const vector_indirect<V, IA> &data, const IA &ia) {
    //     return data.project (ia);
    // }
#endif

    template<class V, class A>
    BOOST_UBLAS_INLINE
    vector_indirect<V, indirect_array<A> > project (V &data, const indirect_array<A> &ia) {
        return vector_indirect<V, indirect_array<A> > (data, ia);
    }
#ifdef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
    template<class V, class A>
    BOOST_UBLAS_INLINE
    const vector_indirect<const V, indirect_array<A> > project_const (const V &data, const indirect_array<A> &ia) {
        return vector_indirect<const V, indirect_array<A> > (data, ia);
    }
#else
    template<class V, class A>
    BOOST_UBLAS_INLINE
    const vector_indirect<const V, indirect_array<A> > project (const V &data, const indirect_array<A> &ia) {
        // ISSUE was: return vector_indirect<V, indirect_array<A> > (const_cast<V &> (data), ia)
        return vector_indirect<const V, indirect_array<A> > (data, ia);
    }
    template<class V, class IA>
    BOOST_UBLAS_INLINE
    vector_indirect<V, IA> project (vector_indirect<V, IA> &data, const typename vector_indirect<V, IA>::range_type &r) {
        return data.project (r);
    }
    template<class V, class IA>
    BOOST_UBLAS_INLINE
    const vector_indirect<V, IA> project (const vector_indirect<V, IA> &data, const typename vector_indirect<V, IA>::range_type &r) {
        return data.project (r);
    }
    template<class V, class IA>
    BOOST_UBLAS_INLINE
    vector_indirect<V, IA> project (vector_indirect<V, IA> &data, const typename vector_indirect<V, IA>::slice_type &s) {
        return data.project (s);
    }
    template<class V, class IA>
    BOOST_UBLAS_INLINE
    const vector_indirect<V, IA> project (const vector_indirect<V, IA> &data, const typename vector_indirect<V, IA>::slice_type &s) {
        return data.project (s);
    }
    template<class V, class A>
    BOOST_UBLAS_INLINE
    vector_indirect<V, indirect_array<A> > project (vector_indirect<V, indirect_array<A> > &data, const indirect_array<A> &ia) {
        return data.project (ia);
    }
    template<class V, class A>
    BOOST_UBLAS_INLINE
    const vector_indirect<V, indirect_array<A> > project (const vector_indirect<V, indirect_array<A> > &data, const indirect_array<A> &ia) {
        return data.project (ia);
    }
#endif

}}}

#endif
