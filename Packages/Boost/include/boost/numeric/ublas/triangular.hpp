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

#ifndef BOOST_UBLAS_TRIANGULAR_H
#define BOOST_UBLAS_TRIANGULAR_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Iterators based on ideas of Jeremy Siek

namespace boost { namespace numeric { namespace ublas {

    // Array based triangular matrix class
    template<class T, class F1, class F2, class A>
    class triangular_matrix:
        public matrix_expression<triangular_matrix<T, F1, F2, A> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<triangular_matrix<T, F1, F2, A> >::operator ();
#endif
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        // typedef const T &const_reference;
        typedef typename type_traits<T>::const_reference const_reference;
        typedef T &reference;
        typedef const T *const_pointer;
        typedef T *pointer;
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef A array_type;
        typedef const A const_array_type;
        typedef const triangular_matrix<T, F1, F2, A> const_self_type;
        typedef triangular_matrix<T, F1, F2, A> self_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const matrix_const_reference<const_self_type> const_closure_type;
#else
        typedef const matrix_reference<const_self_type> const_closure_type;
#endif
        typedef matrix_reference<self_type> closure_type;
        typedef packed_proxy_tag storage_category;
        typedef typename F1::packed_category packed_category;
        typedef typename F2::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        triangular_matrix ():
            matrix_expression<self_type> (),
            size1_ (0), size2_ (0), data_ (0) {}
        BOOST_UBLAS_INLINE
        triangular_matrix (size_type size1, size_type size2):
            matrix_expression<self_type> (),
            size1_ (size1), size2_ (size2), data_ (0) {
            resize (size1, size2);
        }
        BOOST_UBLAS_INLINE
        triangular_matrix (size_type size1, size_type size2, const array_type &data):
            matrix_expression<self_type> (),
            size1_ (size1), size2_ (size2), data_ (data) {}
        BOOST_UBLAS_INLINE
        triangular_matrix (const triangular_matrix &m):
            matrix_expression<self_type> (),
            size1_ (m.size1_), size2_ (m.size2_), data_ (m.data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_matrix (const matrix_expression<AE> &ae):
            matrix_expression<self_type> (),
            size1_ (ae ().size1 ()), size2_ (ae ().size2 ()), data_ (0) {
#ifndef BOOST_UBLAS_TYPE_CHECK
            resize (ae ().size1 (), ae ().size2 (), false);
#else
            resize (ae ().size1 (), ae ().size2 (), true);
#endif
            matrix_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
        }

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const {
            return size1_;
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const {
            return size2_;
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
        void resize (size_type size1, size_type size2, bool preserve = true) {
            size1_ = size1;
            size2_ = size2;
            detail::resize (data (), functor1_type::packed_size (size1, size2), preserve);
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size1_, bad_index ());
            BOOST_UBLAS_CHECK (j < size2_, bad_index ());
            if (functor1_type::other (i, j))
                return data () [functor1_type::element (functor2_type (), i, size1_, j, size2_)];
            else if (functor1_type::one (i, j))
                return one_;
            else
                return zero_;
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size1_, bad_index ());
            BOOST_UBLAS_CHECK (j < size2_, bad_index ());
            if (functor1_type::other (i, j))
                return data () [functor1_type::element (functor2_type (), i, size1_, j, size2_)];
            else if (functor1_type::one (i, j)) {
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
                // Raising exceptions abstracted as requested during review.
                // throw external_logic ();
                external_logic ().raise ();
#endif
                return one_;
            } else {
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
                // Raising exceptions abstracted as requested during review.
                // throw external_logic ();
                external_logic ().raise ();
#endif
                return zero_;
            }
        }

        // Assignment
        BOOST_UBLAS_INLINE
        triangular_matrix &operator = (const triangular_matrix &m) {
            // Precondition for container relaxed as requested during review.
            // BOOST_UBLAS_CHECK (size1_ == m.size1_, bad_size ());
            // BOOST_UBLAS_CHECK (size2_ == m.size2_, bad_size ());
            size1_ = m.size1_;
            size2_ = m.size2_;
            data () = m.data ();
            return *this;
        }
        BOOST_UBLAS_INLINE
        triangular_matrix &assign_temporary (triangular_matrix &m) {
            swap (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_matrix &operator = (const matrix_expression<AE> &ae) {
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
        triangular_matrix &reset (const matrix_expression<AE> &ae) {
            self_type temporary (ae);
            resize (temporary.size1 (), temporary.size2 (), false);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_matrix &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_matrix& operator += (const matrix_expression<AE> &ae) {
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
        triangular_matrix &plus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_matrix& operator -= (const matrix_expression<AE> &ae) {
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
        triangular_matrix &minus_assign (const matrix_expression<AE> &ae) { 
            matrix_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae); 
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        triangular_matrix& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        triangular_matrix& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (triangular_matrix &m) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &m, external_logic ());
            if (this != &m) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size1_ == m.size1_, bad_size ());
                // BOOST_UBLAS_CHECK (size2_ == m.size2_, bad_size ());
                std::swap (size1_, m.size1_);
                std::swap (size2_, m.size2_);
                data ().swap (m.data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (triangular_matrix &m1, triangular_matrix &m2) {
            m1.swap (m2);
        }
#endif

        // Element insertion and erasure
        // These functions should work with std::vector.
        // Thanks to Kresimir Fresl for spotting this.
        BOOST_UBLAS_INLINE
        void insert (size_type i, size_type j, const_reference t) {
            BOOST_UBLAS_CHECK (i < size1_, bad_index ());
            BOOST_UBLAS_CHECK (j < size2_, bad_index ());
// FIXME: is this ugly check still needed?!
// #ifndef BOOST_UBLAS_USE_ET
//             if (t == value_type ())
//                 return;
// #endif
            if (functor1_type::other (i, j)) {
                size_type k = functor1_type::element (functor2_type (), i, size1_, j, size2_);
                BOOST_UBLAS_CHECK (type_traits<value_type>::equals (data () [k], value_type ()), bad_index ());
                // data ().insert (data ().begin () + k, t);
                data () [k] = t;
            } else {
                // Raising exceptions abstracted as requested during review.
                // throw external_logic ();
                external_logic ().raise ();
            }
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size1_, bad_index ());
            BOOST_UBLAS_CHECK (j < size2_, bad_index ());
            if (functor1_type::other (i, j)) {
                size_type k = functor1_type::element (functor2_type (), i, size1_, j, size2_);
                // data ().erase (data ().begin () + k);
                data () [k] = value_type ();
            }
        }
        BOOST_UBLAS_INLINE
        void clear () {
            // data ().clear ();
            std::fill (data ().begin (), data ().end (), value_type ());
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_iterator1<self_type, packed_random_access_iterator_tag> iterator1;
        typedef indexed_iterator2<self_type, packed_random_access_iterator_tag> iterator2;
        typedef indexed_const_iterator1<self_type, packed_random_access_iterator_tag> const_iterator1;
        typedef indexed_const_iterator2<self_type, packed_random_access_iterator_tag> const_iterator2;
#else
        class const_iterator1;
        class iterator1;
        class const_iterator2;
        class iterator2;
#endif
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
        typedef reverse_iterator_base1<iterator1, value_type, reference> reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
        typedef reverse_iterator_base2<iterator2, value_type, reference> reverse_iterator2;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
        typedef reverse_iterator_base1<iterator1> reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
        typedef reverse_iterator_base2<iterator2> reverse_iterator2;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j) const {
            if (rank == 1)
                i = functor1_type::restrict1 (i, j);
            return const_iterator1 (*this, i, j);
        }
        BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j) {
            if (rank == 1)
                i = functor1_type::mutable_restrict1 (i, j);
            return iterator1 (*this, i, j);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            if (rank == 1)
                j = functor1_type::restrict2 (i, j);
            return const_iterator2 (*this, i, j);
        }
        BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j) {
            if (rank == 1)
                j = functor1_type::mutable_restrict2 (i, j);
            return iterator2 (*this, i, j);
        }

        // Iterators simply are indices.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<triangular_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator1, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename triangular_matrix::difference_type difference_type;
            typedef typename triangular_matrix::value_type value_type;
            typedef typename triangular_matrix::const_reference reference;
            typedef typename triangular_matrix::const_pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &m, size_type it1, size_type it2):
                container_const_reference<self_type> (m), it1_ (it1), it2_ (it2) {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const iterator1 &it):
                container_const_reference<self_type> (it ()), it1_ (it.it1_), it2_ (it.it2_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                ++ it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                -- it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator += (difference_type n) {
                it1_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -= (difference_type n) {
                it1_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ - it.it1_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return (*this) () (it1_, it2_);
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 begin () const {
                return (*this) ().find2 (1, it1_, 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 end () const {
                return (*this) ().find2 (1, it1_, (*this) ().size2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_reverse_iterator2 rbegin () const {
                return const_reverse_iterator2 (end ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_reverse_iterator2 rend () const {
                return const_reverse_iterator2 (begin ());
            }
#endif

            // Indices
            BOOST_UBLAS_INLINE
            size_type index1 () const {
                return it1_;
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it2_;
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
                it2_ = it.it2_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ == it.it1_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ < it.it1_;
            }

        private:
            size_type it1_;
            size_type it2_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator1 begin1 () const {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator1 end1 () const {
            return find1 (0, size1_, 0);
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator1:
            public container_reference<triangular_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator1, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename triangular_matrix::difference_type difference_type;
            typedef typename triangular_matrix::value_type value_type;
            typedef typename triangular_matrix::reference reference;
            typedef typename triangular_matrix::pointer pointer;
#endif
            typedef iterator2 dual_iterator_type;
            typedef reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator1 ():
                container_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            iterator1 (self_type &m, size_type it1, size_type it2):
                container_reference<self_type> (m), it1_ (it1), it2_ (it2) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator1 &operator ++ () {
                ++ it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator1 &operator -- () {
                -- it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator1 &operator += (difference_type n) {
                it1_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator1 &operator -= (difference_type n) {
                it1_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ - it.it1_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return (*this) () (it1_, it2_);
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 begin () const {
                return (*this) ().find2 (1, it1_, 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 end () const {
                return (*this) ().find2 (1, it1_, (*this) ().size2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            reverse_iterator2 rbegin () const {
                return reverse_iterator2 (end ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            reverse_iterator2 rend () const {
                return reverse_iterator2 (begin ());
            }
#endif

            // Indices
            BOOST_UBLAS_INLINE
            size_type index1 () const {
                return it1_;
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it2_;
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            iterator1 &operator = (const iterator1 &it) {
                container_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
                it2_ = it.it2_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ == it.it1_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ < it.it1_;
            }

        private:
            size_type it1_;
            size_type it2_;

            friend class const_iterator1;
        };
#endif

        BOOST_UBLAS_INLINE
        iterator1 begin1 () {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        iterator1 end1 () {
            return find1 (0, size1_, 0);
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator2:
            public container_const_reference<triangular_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator2, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename triangular_matrix::difference_type difference_type;
            typedef typename triangular_matrix::value_type value_type;
            typedef typename triangular_matrix::const_reference reference;
            typedef typename triangular_matrix::const_pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &m, size_type it1, size_type it2):
                container_const_reference<self_type> (m), it1_ (it1), it2_ (it2) {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const iterator2 &it):
                container_const_reference<self_type> (it ()), it1_ (it.it1_), it2_ (it.it2_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                ++ it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                -- it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator += (difference_type n) {
                it2_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -= (difference_type n) {
                it2_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ - it.it2_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return (*this) () (it1_, it2_);
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 begin () const {
                return (*this) ().find1 (1, 0, it2_);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 end () const {
                return (*this) ().find1 (1, (*this) ().size1 (), it2_);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_reverse_iterator1 rbegin () const {
                return const_reverse_iterator1 (end ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_reverse_iterator1 rend () const {
                return const_reverse_iterator1 (begin ());
            }
#endif

            // Indices
            BOOST_UBLAS_INLINE
            size_type index1 () const {
                return it1_;
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it2_;
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
                it2_ = it.it2_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ == it.it2_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ < it.it2_;
            }

        private:
            size_type it1_;
            size_type it2_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator2 begin2 () const {
            return find2 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 end2 () const {
            return find2 (0, 0, size2_);
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator2:
            public container_reference<triangular_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator2, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename triangular_matrix::difference_type difference_type;
            typedef typename triangular_matrix::value_type value_type;
            typedef typename triangular_matrix::reference reference;
            typedef typename triangular_matrix::pointer pointer;
#endif
            typedef iterator1 dual_iterator_type;
            typedef reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator2 ():
                container_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            iterator2 (self_type &m, size_type it1, size_type it2):
                container_reference<self_type> (m), it1_ (it1), it2_ (it2) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator2 &operator ++ () {
                ++ it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator2 &operator -- () {
                -- it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator2 &operator += (difference_type n) {
                it2_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator2 &operator -= (difference_type n) {
                it2_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ - it.it2_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return (*this) () (it1_, it2_);
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 begin () const {
                return (*this) ().find1 (1, 0, it2_);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 end () const {
                return (*this) ().find1 (1, (*this) ().size1 (), it2_);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            reverse_iterator1 rbegin () const {
                return reverse_iterator1 (end ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            reverse_iterator1 rend () const {
                return reverse_iterator1 (begin ());
            }
#endif

            // Indices
            BOOST_UBLAS_INLINE
            size_type index1 () const {
                return it1_;
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it2_;
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            iterator2 &operator = (const iterator2 &it) {
                container_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
                it2_ = it.it2_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ == it.it2_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ < it.it2_;
            }

        private:
            size_type it1_;
            size_type it2_;

            friend class const_iterator2;
        };
#endif

        BOOST_UBLAS_INLINE
        iterator2 begin2 () {
            return find2 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        iterator2 end2 () {
            return find2 (0, 0, size2_);
        }

        // Reverse iterators

        BOOST_UBLAS_INLINE
        const_reverse_iterator1 rbegin1 () const {
            return const_reverse_iterator1 (end1 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator1 rend1 () const {
            return const_reverse_iterator1 (begin1 ());
        }

        BOOST_UBLAS_INLINE
        reverse_iterator1 rbegin1 () {
            return reverse_iterator1 (end1 ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator1 rend1 () {
            return reverse_iterator1 (begin1 ());
        }

        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

        BOOST_UBLAS_INLINE
        reverse_iterator2 rbegin2 () {
            return reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator2 rend2 () {
            return reverse_iterator2 (begin2 ());
        }

    private:
        size_type size1_;
        size_type size2_;
        array_type data_;
        static value_type zero_;
        static value_type one_;
    };

    template<class T, class F1, class F2, class A>
    typename triangular_matrix<T, F1, F2, A>::value_type triangular_matrix<T, F1, F2, A>::zero_ =
        triangular_matrix<T, F1, F2, A>::value_type ();
    template<class T, class F1, class F2, class A>
    typename triangular_matrix<T, F1, F2, A>::value_type triangular_matrix<T, F1, F2, A>::one_ =
        triangular_matrix<T, F1, F2, A>::value_type (1);

    // Triangular matrix adaptor class
    template<class M, class F>
    class triangular_adaptor:
        public matrix_expression<triangular_adaptor<M, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<triangular_adaptor<M, F> >::operator ();
#endif
        typedef const M const_matrix_type;
        typedef M matrix_type;
        typedef F functor_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
#ifndef BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
        typedef typename M::const_reference const_reference;
        typedef typename M::reference reference;
        typedef typename M::const_pointer const_pointer;
        typedef typename M::pointer pointer;
#else
        typedef typename M::const_reference const_reference;
        typedef typename boost::mpl::if_c<boost::is_const<M>::value,
                                          typename M::const_reference,
                                          typename M::reference>::type reference;
        typedef typename M::const_pointer const_pointer;
        typedef typename boost::mpl::if_c<boost::is_const<M>::value,
                                          typename M::const_pointer,
                                          typename M::pointer>::type pointer;
#endif
#ifndef BOOST_UBLAS_CT_PROXY_CLOSURE_TYPEDEFS
        typedef typename M::closure_type matrix_closure_type;
#else
        typedef typename boost::mpl::if_c<boost::is_const<M>::value,
                                          typename M::const_closure_type,
                                          typename M::closure_type>::type matrix_closure_type;
#endif
        typedef const triangular_adaptor<M, F> const_self_type;
        typedef triangular_adaptor<M, F> self_type;
        typedef const_self_type const_closure_type;
        typedef self_type closure_type;
#ifndef BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
        typedef typename M::const_iterator1 const_iterator1_type;
        typedef typename M::iterator1 iterator1_type;
        typedef typename M::const_iterator2 const_iterator2_type;
        typedef typename M::iterator2 iterator2_type;
#else
        typedef typename M::const_iterator1 const_iterator1_type;
        typedef typename boost::mpl::if_c<boost::is_const<M>::value,
                                          typename M::const_iterator1,
                                          typename M::iterator1>::type iterator1_type;
        typedef typename M::const_iterator2 const_iterator2_type;
        typedef typename boost::mpl::if_c<boost::is_const<M>::value,
                                          typename M::const_iterator2,
                                          typename M::iterator2>::type iterator2_type;
#endif
        typedef typename storage_restrict_traits<typename M::storage_category,
                                                 packed_proxy_tag>::storage_category storage_category;
        typedef typename F::packed_category packed_category;
        typedef typename M::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        triangular_adaptor ():
            matrix_expression<self_type> (),
            data_ (nil_) {}
        BOOST_UBLAS_INLINE
        triangular_adaptor (matrix_type &data):
            matrix_expression<self_type> (),
            data_ (data) {}
        BOOST_UBLAS_INLINE
        triangular_adaptor (const triangular_adaptor &m):
            matrix_expression<self_type> (),
            data_ (m.data_) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const {
            return data_.size1 ();
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const {
            return data_.size2 ();
        }
        BOOST_UBLAS_INLINE
        const matrix_closure_type &data () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        matrix_closure_type &data () {
            return data_;
        }

#ifdef BOOST_UBLAS_DEPRECATED
        // Resetting
        BOOST_UBLAS_INLINE
        void reset (matrix_type &data) {
            // References are not retargetable.
            // Thanks to Michael Stevens for spotting this.
            // data_ = data;
            data_.reset (data);
        }
#endif

        // Element access
#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
            if (functor_type::other (i, j))
                return data () (i, j);
            else if (functor_type::one (i, j))
                return one_;
            else
                return zero_;
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
            if (functor_type::other (i, j))
                return data () (i, j);
            else if (functor_type::one (i, j)) {
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
                // Raising exceptions abstracted as requested during review.
                // throw external_logic ();
                external_logic ().raise ();
#endif
                return one_;
            } else {
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
                // Raising exceptions abstracted as requested during review.
                // throw external_logic ();
                external_logic ().raise ();
#endif
                return zero_;
            }
        }
#else
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
            if (functor_type::other (i, j))
                return data () (i, j);
            else if (functor_type::one (i, j)) {
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
                // Raising exceptions abstracted as requested during review.
                // throw external_logic ();
                external_logic ().raise ();
#endif
                return one_;
            } else {
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
                // Raising exceptions abstracted as requested during review.
                // throw external_logic ();
                external_logic ().raise ();
#endif
                return zero_;
            }
        }
#endif

        // Assignment
        BOOST_UBLAS_INLINE
        triangular_adaptor &operator = (const triangular_adaptor &m) {
            matrix_assign (scalar_assign<reference, value_type> (), *this, m);
            return *this;
        }
        BOOST_UBLAS_INLINE
        triangular_adaptor &assign_temporary (triangular_adaptor &m) {
            *this = m;
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_adaptor &operator = (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, value_type> (), *this, matrix<value_type> (ae)); 
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_adaptor &assign (const matrix_expression<AE> &ae) { 
            matrix_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae); 
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_adaptor& operator += (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, value_type> (), *this, matrix<value_type> (*this + ae)); 
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_adaptor &plus_assign (const matrix_expression<AE> &ae) { 
            matrix_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae); 
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_adaptor& operator -= (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, value_type> (), *this, matrix<value_type> (*this - ae)); 
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        triangular_adaptor &minus_assign (const matrix_expression<AE> &ae) { 
            matrix_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae); 
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        triangular_adaptor& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        triangular_adaptor& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Comparison
        bool operator == (const triangular_adaptor &ta) const {
            return (*this).data () == ta.data ();
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (triangular_adaptor &m) {
            // Too unusual semantic.
            // BOOST_UBLAS_CHECK (this != &m, external_logic ());
            if (this != &m)
                matrix_swap (scalar_swap<reference, reference> (), *this, m); 
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (triangular_adaptor &m1, triangular_adaptor &m2) {
            m1.swap (m2);
        }
#endif

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_iterator1<self_type, packed_random_access_iterator_tag> iterator1;
        typedef indexed_iterator2<self_type, packed_random_access_iterator_tag> iterator2;
        typedef indexed_const_iterator1<self_type, packed_random_access_iterator_tag> const_iterator1;
        typedef indexed_const_iterator2<self_type, packed_random_access_iterator_tag> const_iterator2;
#else
        class const_iterator1;
        class iterator1;
        class const_iterator2;
        class iterator2;
#endif
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
        typedef reverse_iterator_base1<iterator1, value_type, reference> reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
        typedef reverse_iterator_base2<iterator2, value_type, reference> reverse_iterator2;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
        typedef reverse_iterator_base1<iterator1> reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
        typedef reverse_iterator_base2<iterator2> reverse_iterator2;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j) const {
            if (rank == 1)
                i = functor_type::restrict1 (i, j);
            return const_iterator1 (*this, data ().find1 (rank, i, j));
        }
        BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j) {
            if (rank == 1)
                i = functor_type::mutable_restrict1 (i, j);
            return iterator1 (*this, data ().find1 (rank, i, j));
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            if (rank == 1)
                j = functor_type::restrict2 (i, j);
            return const_iterator2 (*this, data ().find2 (rank, i, j));
        }
        BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j) {
            if (rank == 1)
                j = functor_type::mutable_restrict2 (i, j);
            return iterator2 (*this, data ().find2 (rank, i, j));
        }

        // Iterators simply are indices.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<triangular_adaptor>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator1, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename const_iterator1_type::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename const_iterator1_type::difference_type difference_type;
            typedef typename const_iterator1_type::value_type value_type;
            typedef typename const_iterator1_type::reference reference;
            typedef typename const_iterator1_type::pointer pointer;
#else
            typedef typename iterator_restrict_traits<typename M::const_iterator1::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
            typedef const_reference reference;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), it1_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &m, const const_iterator1_type &it1):
                container_const_reference<self_type> (m), it1_ (it1) {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const iterator1 &it):
                container_const_reference<self_type> (it ()), it1_ (it.it1_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                ++ it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                -- it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator += (difference_type n) {
                it1_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -= (difference_type n) {
                it1_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it1_ - it.it1_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                size_type i = index1 ();
                size_type j = index2 ();
                BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
                if (functor_type::other (i, j))
                    return *it1_;
                else
                    return (*this) () (i, j);
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 begin () const {
                return (*this) ().find2 (1, index1 (), 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 end () const {
                return (*this) ().find2 (1, index1 (), (*this) ().size2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_reverse_iterator2 rbegin () const {
                return const_reverse_iterator2 (end ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_reverse_iterator2 rend () const {
                return const_reverse_iterator2 (begin ());
            }
#endif

            // Indices
            BOOST_UBLAS_INLINE
            size_type index1 () const {
                return it1_.index1 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it1_.index2 ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it1_ == it.it1_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it1_ < it.it1_;
            }

        private:
            const_iterator1_type it1_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator1 begin1 () const {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator1 end1 () const {
            return find1 (0, size1 (), 0);
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator1:
            public container_reference<triangular_adaptor>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator1, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename iterator1_type::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename iterator1_type::difference_type difference_type;
            typedef typename iterator1_type::value_type value_type;
            typedef typename iterator1_type::reference reference;
            typedef typename iterator1_type::pointer pointer;
#else
            typedef typename iterator_restrict_traits<typename M::iterator1::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
#endif
            typedef iterator2 dual_iterator_type;
            typedef reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator1 ():
                container_reference<self_type> (), it1_ () {}
            BOOST_UBLAS_INLINE
            iterator1 (self_type &m, const iterator1_type &it1):
                container_reference<self_type> (m), it1_ (it1) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator1 &operator ++ () {
                ++ it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator1 &operator -- () {
                -- it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator1 &operator += (difference_type n) {
                it1_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator1 &operator -= (difference_type n) {
                it1_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it1_ - it.it1_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                size_type i = index1 ();
                size_type j = index2 ();
                BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
                if (functor_type::other (i, j))
                    return *it1_;
                else
                    return (*this) () (i, j);
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 begin () const {
                return (*this) ().find2 (1, index1 (), 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 end () const {
                return (*this) ().find2 (1, index1 (), (*this) ().size2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            reverse_iterator2 rbegin () const {
                return reverse_iterator2 (end ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            reverse_iterator2 rend () const {
                return reverse_iterator2 (begin ());
            }
#endif

            // Indices
            BOOST_UBLAS_INLINE
            size_type index1 () const {
                return it1_.index1 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it1_.index2 ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator1 &operator = (const iterator1 &it) {
                container_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it1_ == it.it1_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it1_ < it.it1_;
            }

        private:
            iterator1_type it1_;

            friend class const_iterator1;
        };
#endif

        BOOST_UBLAS_INLINE
        iterator1 begin1 () {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        iterator1 end1 () {
            return find1 (0, size1 (), 0);
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator2:
            public container_const_reference<triangular_adaptor>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator2, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename const_iterator2_type::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename const_iterator2_type::difference_type difference_type;
            typedef typename const_iterator2_type::value_type value_type;
            typedef typename const_iterator2_type::reference reference;
            typedef typename const_iterator2_type::pointer pointer;
#else
            typedef typename iterator_restrict_traits<typename M::const_iterator2::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
            typedef const_reference reference;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &m, const const_iterator2_type &it2):
                container_const_reference<self_type> (m), it2_ (it2) {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const iterator2 &it):
                container_const_reference<self_type> (it ()), it2_ (it.it2_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                ++ it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                -- it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator += (difference_type n) {
                it2_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -= (difference_type n) {
                it2_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it2_ - it.it2_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                size_type i = index1 ();
                size_type j = index2 ();
                BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
                if (functor_type::other (i, j))
                    return *it2_;
                else
                    return (*this) () (i, j);
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 begin () const {
                return (*this) ().find1 (1, 0, index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 end () const {
                return (*this) ().find1 (1, (*this) ().size1 (), index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_reverse_iterator1 rbegin () const {
                return const_reverse_iterator1 (end ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_reverse_iterator1 rend () const {
                return const_reverse_iterator1 (begin ());
            }
#endif

            // Indices
            BOOST_UBLAS_INLINE
            size_type index1 () const {
                return it2_.index1 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it2_.index2 ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                it2_ = it.it2_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it2_ == it.it2_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it2_ < it.it2_;
            }

        private:
            const_iterator2_type it2_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator2 begin2 () const {
            return find2 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 end2 () const {
            return find2 (0, 0, size2 ());
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator2:
            public container_reference<triangular_adaptor>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator2, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename iterator2_type::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename iterator2_type::difference_type difference_type;
            typedef typename iterator2_type::value_type value_type;
            typedef typename iterator2_type::reference reference;
            typedef typename iterator2_type::pointer pointer;
#else
            typedef typename iterator_restrict_traits<typename M::iterator2::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
#endif
            typedef iterator1 dual_iterator_type;
            typedef reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator2 ():
                container_reference<self_type> (), it2_ () {}
            BOOST_UBLAS_INLINE
            iterator2 (self_type &m, const iterator2_type &it2):
                container_reference<self_type> (m), it2_ (it2) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator2 &operator ++ () {
                ++ it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator2 &operator -- () {
                -- it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator2 &operator += (difference_type n) {
                it2_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator2 &operator -= (difference_type n) {
                it2_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it2_ - it.it2_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                size_type i = index1 ();
                size_type j = index2 ();
                BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
                if (functor_type::other (i, j))
                    return *it2_;
                else
                    return (*this) () (i, j);
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 begin () const {
                return (*this) ().find1 (1, 0, index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 end () const {
                return (*this) ().find1 (1, (*this) ().size1 (), index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            reverse_iterator1 rbegin () const {
                return reverse_iterator1 (end ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            reverse_iterator1 rend () const {
                return reverse_iterator1 (begin ());
            }
#endif

            // Indices
            BOOST_UBLAS_INLINE
            size_type index1 () const {
                return it2_.index1 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it2_.index2 ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator2 &operator = (const iterator2 &it) {
                container_reference<self_type>::assign (&it ());
                it2_ = it.it2_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it2_ == it.it2_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it2_ < it.it2_;
            }

        private:
            iterator2_type it2_;

            friend class const_iterator2;
        };
#endif

        BOOST_UBLAS_INLINE
        iterator2 begin2 () {
            return find2 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        iterator2 end2 () {
            return find2 (0, 0, size2 ());
        }

        // Reverse iterators

        BOOST_UBLAS_INLINE
        const_reverse_iterator1 rbegin1 () const {
            return const_reverse_iterator1 (end1 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator1 rend1 () const {
            return const_reverse_iterator1 (begin1 ());
        }

        BOOST_UBLAS_INLINE
        reverse_iterator1 rbegin1 () {
            return reverse_iterator1 (end1 ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator1 rend1 () {
            return reverse_iterator1 (begin1 ());
        }

        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

        BOOST_UBLAS_INLINE
        reverse_iterator2 rbegin2 () {
            return reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator2 rend2 () {
            return reverse_iterator2 (begin2 ());
        }

    private:
        matrix_closure_type data_;
        static matrix_type nil_;
        static value_type zero_;
        static value_type one_;
    };

    template<class M, class F>
    typename triangular_adaptor<M, F>::matrix_type triangular_adaptor<M, F>::nil_;
    template<class M, class F>
    typename triangular_adaptor<M, F>::value_type triangular_adaptor<M, F>::zero_ =
        triangular_adaptor<M, F>::value_type ();
    template<class M, class F>
    typename triangular_adaptor<M, F>::value_type triangular_adaptor<M, F>::one_ =
        triangular_adaptor<M, F>::value_type (1);

    template<class E1, class E2>
    struct matrix_vector_solve_traits {
        typedef typename promote_traits<typename E1::value_type, typename E2::value_type>::promote_type promote_type;
        typedef vector<promote_type> result_type;
    };

    // Operations:
    //  n * (n - 1) / 2 + n = n * (n + 1) / 2 multiplications,
    //  n * (n - 1) / 2 additions

    // Dense (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        lower_tag, column_major_tag, dense_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size (), bad_size ());
        size_type size = e2 ().size ();
        for (size_type n = 0; n < size; ++ n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e2 () (n) /= e1 () (n, n);
            if (t != value_type ()) {
                for (size_type m = n + 1; m < size; ++ m)
                    e2 () (m) -= e1 () (m, n) * t;
            }
        }
    }
    // Packed (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        lower_tag, column_major_tag, packed_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size (), bad_size ());
        size_type size = e2 ().size ();
        for (size_type n = 0; n < size; ++ n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e2 () (n) /= e1 () (n, n);
            if (t != value_type ()) {
                typename E1::const_iterator1 it1e1 (e1 ().find1 (1, n + 1, n));
                typename E1::const_iterator1 it1e1_end (e1 ().find1 (1, e1 ().size1 (), n));
                difference_type m (it1e1_end - it1e1);
                while (-- m >= 0)
                    e2 () (it1e1.index1 ()) -= *it1e1 * t, ++ it1e1;
            }
        }
    }
    // Sparse (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        lower_tag, column_major_tag, unknown_storage_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size (), bad_size ());
        size_type size = e2 ().size ();
        for (size_type n = 0; n < size; ++ n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e2 () (n) /= e1 () (n, n);
            if (t != value_type ()) {
                typename E1::const_iterator1 it1e1 (e1 ().find1 (1, n + 1, n));
                typename E1::const_iterator1 it1e1_end (e1 ().find1 (1, e1 ().size1 (), n));
                while (it1e1 != it1e1_end)
                    e2 () (it1e1.index1 ()) -= *it1e1 * t, ++ it1e1;
            }
        }
    }
    // Redirectors :-)
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        lower_tag, column_major_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category storage_category;
        inplace_solve (e1, e2,
                       lower_tag (), column_major_tag (), storage_category ());
    }
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        lower_tag, row_major_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category storage_category;
        inplace_solve (e2, trans (e1),
                       upper_tag (), row_major_tag (), storage_category ());
    }
    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        lower_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::orientation_category orientation_category;
        inplace_solve (e1, e2,
                       lower_tag (), orientation_category ());
    }
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        unit_lower_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::orientation_category orientation_category;
        inplace_solve (triangular_adaptor<const E1, unit_lower> (e1 ()), e2,
                       unit_lower_tag (), orientation_category ());
    }

    // Dense (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        upper_tag, column_major_tag, dense_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size (), bad_size ());
        size_type size = e2 ().size ();
        for (difference_type n = size - 1; n >= 0; -- n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e2 () (n) /= e1 () (n, n);
            if (t != value_type ()) {
                for (difference_type m = n - 1; m >= 0; -- m)
                    e2 () (m) -= e1 () (m, n) * t;
            }
        }
    }
    // Packed (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        upper_tag, column_major_tag, packed_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size (), bad_size ());
        size_type size = e2 ().size ();
        for (difference_type n = size - 1; n >= 0; -- n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e2 () (n) /= e1 () (n, n);
            if (t != value_type ()) {
                typename E1::const_reverse_iterator1 it1e1 (e1 ().find1 (1, n, n));
                typename E1::const_reverse_iterator1 it1e1_rend (e1 ().find1 (1, 0, n));
                difference_type m (it1e1_rend - it1e1);
                while (-- m >= 0)
                    e2 () (it1e1.index1 ()) -= *it1e1 * t, ++ it1e1;
            }
        }
    }
    // Sparse (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        upper_tag, column_major_tag, unknown_storage_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size (), bad_size ());
        size_type size = e2 ().size ();
        for (difference_type n = size - 1; n >= 0; -- n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e2 () (n) /= e1 () (n, n);
            if (t != value_type ()) {
                typename E1::const_reverse_iterator1 it1e1 (e1 ().find1 (1, n, n));
                typename E1::const_reverse_iterator1 it1e1_rend (e1 ().find1 (1, 0, n));
                while (it1e1 != it1e1_rend)
                    e2 () (it1e1.index1 ()) -= *it1e1 * t, ++ it1e1;
            }
        }
    }
    // Redirectors :-)
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        upper_tag, column_major_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category storage_category;
        inplace_solve (e1, e2,
                       upper_tag (), column_major_tag (), storage_category ());
    }
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        upper_tag, row_major_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category storage_category;
        inplace_solve (e2, trans (e1),
                       lower_tag (), row_major_tag (), storage_category ());
    }
    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        upper_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::orientation_category orientation_category;
        inplace_solve (e1, e2,
                       upper_tag (), orientation_category ());
    }
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, vector_expression<E2> &e2,
                        unit_upper_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::orientation_category orientation_category;
        inplace_solve (triangular_adaptor<const E1, unit_upper> (e1 ()), e2,
                       unit_upper_tag (), orientation_category ());
    }

#ifdef BOOST_UBLAS_DEPRECATED
    template<class E1, class E2, class C>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1,
                        vector_expression<E2> &e2,
                        C) {
        inplace_solve (e1, e2, C ());
    }
#endif
    template<class E1, class E2, class C>
    BOOST_UBLAS_INLINE
    typename matrix_vector_solve_traits<E1, E2>::result_type
    solve (const matrix_expression<E1> &e1,
           const vector_expression<E2> &e2,
           C) {
        typename matrix_vector_solve_traits<E1, E2>::result_type r (e2);
        inplace_solve (e1, r, C ());
        return r;
    }

    // Dense (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        lower_tag, row_major_tag, dense_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E1::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E1::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size () == e2 ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (e2 ().size1 () == e2 ().size2 (), bad_size ());
        size_type size = e1 ().size ();
        for (difference_type n = size - 1; n >= 0; -- n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e2 () (n, n) != value_type (), singular ());
#else
            if (e2 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e1 () (n) /= e2 () (n, n);
            if (t != value_type ()) {
                for (difference_type m = n - 1; m >= 0; -- m)
                    e1 () (m) -= t * e2 () (n, m);
            }
        }
    }
    // Packed (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        lower_tag, row_major_tag, packed_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E1::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E1::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size () == e2 ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (e2 ().size1 () == e2 ().size2 (), bad_size ());
        size_type size = e1 ().size ();
        for (difference_type n = size - 1; n >= 0; -- n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e2 () (n, n) != value_type (), singular ());
#else
            if (e2 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e1 () (n) /= e2 () (n, n);
            if (t != value_type ()) {
                typename E2::const_reverse_iterator2 it2e2 (e2 ().find2 (1, n, n));
                typename E2::const_reverse_iterator2 it2e2_rend (e2 ().find2 (1, n, 0));
                difference_type m (it2e2_rend - it2e2);
                while (-- m >= 0)
                    e1 () (it2e2.index2 ()) -= *it2e2 * t, ++ it2e2;
            }
        }
    }
    // Sparse (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        lower_tag, row_major_tag, unknown_storage_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E1::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E1::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size () == e2 ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (e2 ().size1 () == e2 ().size2 (), bad_size ());
        size_type size = e1 ().size ();
        for (difference_type n = size - 1; n >= 0; -- n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e2 () (n, n) != value_type (), singular ());
#else
            if (e2 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e1 () (n) /= e2 () (n, n);
            if (t != value_type ()) {
                typename E2::const_reverse_iterator2 it2e2 (e2 ().find2 (1, n, n));
                typename E2::const_reverse_iterator2 it2e2_rend (e2 ().find2 (1, n, 0));
                while (it2e2 != it2e2_rend)
                    e1 () (it2e2.index2 ()) -= *it2e2 * t, ++ it2e2;
            }
        }
    }
    // Redirectors :-)
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        lower_tag, row_major_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category storage_category;
        inplace_solve (e1, e2,
                       lower_tag (), row_major_tag (), storage_category ());
    }
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        lower_tag, column_major_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category storage_category;
        inplace_solve (trans (e2), e1,
                       upper_tag (), row_major_tag (), storage_category ());
    }
    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        lower_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::orientation_category orientation_category;
        inplace_solve (e1, e2,
                       lower_tag (), orientation_category ());
    }
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        unit_lower_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::orientation_category orientation_category;
        inplace_solve (e1, triangular_adaptor<const E2, unit_lower> (e2 ()),
                       unit_lower_tag (), orientation_category ());
    }

    // Dense (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        upper_tag, row_major_tag, dense_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E1::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E1::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size () == e2 ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (e2 ().size1 () == e2 ().size2 (), bad_size ());
        size_type size = e1 ().size ();
        for (size_type n = 0; n < size; ++ n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e2 () (n, n) != value_type (), singular ());
#else
            if (e2 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e1 () (n) /= e2 () (n, n);
            if (t != value_type ()) {
                for (size_type m = n + 1; m < size; ++ m)
                    e1 () (m) -= t * e2 () (n, m);
            }
        }
    }
    // Packed (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        upper_tag, row_major_tag, packed_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E1::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E1::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size () == e2 ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (e2 ().size1 () == e2 ().size2 (), bad_size ());
        size_type size = e1 ().size ();
        for (size_type n = 0; n < size; ++ n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e2 () (n, n) != value_type (), singular ());
#else
            if (e2 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e1 () (n) /= e2 () (n, n);
            if (t != value_type ()) {
                typename E2::const_iterator2 it2e2 (e2 ().find2 (1, n, n + 1));
                typename E2::const_iterator2 it2e2_end (e2 ().find2 (1, n, e2 ().size2 ()));
                difference_type m (it2e2_end - it2e2);
                while (-- m >= 0)
                    e1 () (it2e2.index2 ()) -= *it2e2 * t, ++ it2e2;
            }
        }
    }
    // Sparse (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        upper_tag, row_major_tag, unknown_storage_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E1::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E1::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size () == e2 ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (e2 ().size1 () == e2 ().size2 (), bad_size ());
        size_type size = e1 ().size ();
        for (size_type n = 0; n < size; ++ n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e2 () (n, n) != value_type (), singular ());
#else
            if (e2 () (n, n) == value_type ())
                singular ().raise ();
#endif
            value_type t = e1 () (n) /= e2 () (n, n);
            if (t != value_type ()) {
                typename E2::const_iterator2 it2e2 (e2 ().find2 (1, n, n + 1));
                typename E2::const_iterator2 it2e2_end (e2 ().find2 (1, n, e2 ().size2 ()));
                while (it2e2 != it2e2_end)
                    e1 () (it2e2.index2 ()) -= *it2e2 * t, ++ it2e2;
            }
        }
    }
    // Redirectors :-)
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        upper_tag, row_major_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category storage_category;
        inplace_solve (e1, e2,
                       upper_tag (), row_major_tag (), storage_category ());
    }
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        upper_tag, column_major_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category storage_category;
        inplace_solve (trans (e2), e1,
                       lower_tag (), row_major_tag (), storage_category ());
    }
    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        upper_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::orientation_category orientation_category;
        inplace_solve (e1, e2,
                       upper_tag (), orientation_category ());
    }
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1, const matrix_expression<E2> &e2,
                        unit_upper_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::orientation_category orientation_category;
        inplace_solve (e1, triangular_adaptor<const E2, unit_upper> (e2 ()),
                       unit_upper_tag (), orientation_category ());
    }

#ifdef BOOST_UBLAS_DEPRECATED
    template<class E1, class E2, class C>
    BOOST_UBLAS_INLINE
    void inplace_solve (vector_expression<E1> &e1,
                        const matrix_expression<E2> &e2,
                        C) {
        inplace_solve (e1 (), e2, C ());
    }
#endif
    template<class E1, class E2, class C>
    BOOST_UBLAS_INLINE
    typename matrix_vector_solve_traits<E1, E2>::result_type
    solve (const vector_expression<E1> &e1,
           const matrix_expression<E2> &e2,
           C) {
        typename matrix_vector_solve_traits<E1, E2>::result_type r (e1);
        inplace_solve (r, e2, C ());
        return r;
    }

    template<class E1, class E2>
    struct matrix_matrix_solve_traits {
        typedef typename promote_traits<typename E1::value_type, typename E2::value_type>::promote_type promote_type;
        typedef matrix<promote_type> result_type;
    };

    // Operations:
    //  k * n * (n - 1) / 2 + k * n = k * n * (n + 1) / 2 multiplications,
    //  k * n * (n - 1) / 2 additions

    // Dense (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, matrix_expression<E2> &e2,
                        lower_tag, dense_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size1 (), bad_size ());
        size_type size1 = e2 ().size1 ();
        size_type size2 = e2 ().size2 ();
        for (size_type n = 0; n < size1; ++ n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            for (size_type l = 0; l < size2; ++ l) {
                value_type t = e2 () (n, l) /= e1 () (n, n);
                if (t != value_type ()) {
                    for (size_type m = n + 1; m < size1; ++ m)
                        e2 () (m, l) -= e1 () (m, n) * t;
                }
            }
        }
    }
    // Packed (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, matrix_expression<E2> &e2,
                        lower_tag, packed_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size1 (), bad_size ());
        size_type size1 = e2 ().size1 ();
        size_type size2 = e2 ().size2 ();
        for (size_type n = 0; n < size1; ++ n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            for (size_type l = 0; l < size2; ++ l) {
                value_type t = e2 () (n, l) /= e1 () (n, n);
                if (t != value_type ()) {
                    typename E1::const_iterator1 it1e1 (e1 ().find1 (1, n + 1, n));
                    typename E1::const_iterator1 it1e1_end (e1 ().find1 (1, e1 ().size1 (), n));
                    difference_type m (it1e1_end - it1e1);
                    while (-- m >= 0)
                        e2 () (it1e1.index1 (), l) -= *it1e1 * t, ++ it1e1;
                }
            }
        }
    }
    // Sparse (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, matrix_expression<E2> &e2,
                        lower_tag, unknown_storage_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size1 (), bad_size ());
        size_type size1 = e2 ().size1 ();
        size_type size2 = e2 ().size2 ();
        for (size_type n = 0; n < size1; ++ n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            for (size_type l = 0; l < size2; ++ l) {
                value_type t = e2 () (n, l) /= e1 () (n, n);
                if (t != value_type ()) {
                    typename E1::const_iterator1 it1e1 (e1 ().find1 (1, n + 1, n));
                    typename E1::const_iterator1 it1e1_end (e1 ().find1 (1, e1 ().size1 (), n));
                    while (it1e1 != it1e1_end)
                        e2 () (it1e1.index1 (), l) -= *it1e1 * t, ++ it1e1;
                }
            }
        }
    }
    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, matrix_expression<E2> &e2,
                        lower_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category dispatch_category;
        inplace_solve (e1, e2,
                       lower_tag (), dispatch_category ());
    }
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, matrix_expression<E2> &e2,
                        unit_lower_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category dispatch_category;
        inplace_solve (triangular_adaptor<const E1, unit_lower> (e1 ()), e2,
                       unit_lower_tag (), dispatch_category ());
    }

    // Dense (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, matrix_expression<E2> &e2,
                        upper_tag, dense_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size1 (), bad_size ());
        size_type size1 = e2 ().size1 ();
        size_type size2 = e2 ().size2 ();
        for (difference_type n = size1 - 1; n >= 0; -- n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            for (difference_type l = size2 - 1; l >= 0; -- l) {
                value_type t = e2 () (n, l) /= e1 () (n, n);
                if (t != value_type ()) {
                    for (difference_type m = n - 1; m >= 0; -- m)
                        e2 () (m, l) -= e1 () (m, n) * t;
                }
            }
        }
    }
    // Packed (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, matrix_expression<E2> &e2,
                        upper_tag, packed_proxy_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size1 (), bad_size ());
        size_type size1 = e2 ().size1 ();
        size_type size2 = e2 ().size2 ();
        for (difference_type n = size1 - 1; n >= 0; -- n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            for (difference_type l = size2 - 1; l >= 0; -- l) {
                value_type t = e2 () (n, l) /= e1 () (n, n);
                if (t != value_type ()) {
                    typename E1::const_reverse_iterator1 it1e1 (e1 ().find1 (1, n, n));
                    typename E1::const_reverse_iterator1 it1e1_rend (e1 ().find1 (1, 0, n));
                    difference_type m (it1e1_rend - it1e1);
                    while (-- m >= 0)
                        e2 () (it1e1.index1 (), l) -= *it1e1 * t, ++ it1e1;
                }
            }
        }
    }
    // Sparse (proxy) case
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, matrix_expression<E2> &e2,
                        upper_tag, unknown_storage_tag) {
        typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
        typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
        typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;

        BOOST_UBLAS_CHECK (e1 ().size1 () == e1 ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size1 (), bad_size ());
        size_type size1 = e2 ().size1 ();
        size_type size2 = e2 ().size2 ();
        for (difference_type n = size1 - 1; n >= 0; -- n) {
#ifndef BOOST_UBLAS_SINGULAR_CHECK
            BOOST_UBLAS_CHECK (e1 () (n, n) != value_type (), singular ());
#else
            if (e1 () (n, n) == value_type ())
                singular ().raise ();
#endif
            for (difference_type l = size2 - 1; l >= 0; -- l) {
                value_type t = e2 () (n, l) /= e1 () (n, n);
                if (t != value_type ()) {
                    typename E1::const_reverse_iterator1 it1e1 (e1 ().find1 (1, n, n));
                    typename E1::const_reverse_iterator1 it1e1_rend (e1 ().find1 (1, 0, n));
                    while (it1e1 != it1e1_rend)
                        e2 () (it1e1.index1 (), l) -= *it1e1 * t, ++ it1e1;
                }
            }
        }
    }
    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, matrix_expression<E2> &e2,
                        upper_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category dispatch_category;
        inplace_solve (e1, e2,
                       upper_tag (), dispatch_category ());
    }
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1, matrix_expression<E2> &e2,
                        unit_upper_tag) {
        typedef BOOST_UBLAS_TYPENAME E1::storage_category dispatch_category;
        inplace_solve (triangular_adaptor<const E1, unit_upper> (e1 ()), e2,
                       unit_upper_tag (), dispatch_category ());
    }

#ifdef BOOST_UBLAS_DEPRECATED
    template<class E1, class E2, class C>
    BOOST_UBLAS_INLINE
    void inplace_solve (const matrix_expression<E1> &e1,
                        matrix_expression<E2> &e2,
                        C) {
        inplace_solve (e1, e2 (), C ());
    }
#endif
    template<class E1, class E2, class C>
    BOOST_UBLAS_INLINE
    typename matrix_matrix_solve_traits<E1, E2>::result_type
    solve (const matrix_expression<E1> &e1,
           const matrix_expression<E2> &e2,
           C) {
        typename matrix_matrix_solve_traits<E1, E2>::result_type r (e2);
        inplace_solve (e1, r, C ());
        return r;
    }

}}}

#endif





