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

#ifndef BOOST_UBLAS_SYMMETRIC_H
#define BOOST_UBLAS_SYMMETRIC_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Iterators based on ideas of Jeremy Siek
// Symmetric matrices are square. Thanks to Peter Schmitteckert for spotting this.

namespace boost { namespace numeric { namespace ublas {

    template<class M>
    bool is_symmetric (const M &m) {
        typedef typename M::size_type size_type;

        if (m.size1 () != m.size2 ())
            return false;
        size_type size = BOOST_UBLAS_SAME (m.size1 (), m.size2 ());
        for (size_type i = 0; i < size; ++ i) {
            for (size_type j = i; j < size; ++ j) {
                if (m (i, j) != m (j, i))
                    return false;
            }
        }
        return true;
    }

    // Array based symmetric matrix class
    template<class T, class F1, class F2, class A>
    class symmetric_matrix:
        public matrix_expression<symmetric_matrix<T, F1, F2, A> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<symmetric_matrix<T, F1, F2, A> >::operator ();
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
        typedef const symmetric_matrix<T, F1, F2, A> const_self_type;
        typedef symmetric_matrix<T, F1, F2, A> self_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const matrix_const_reference<const_self_type> const_closure_type;
#else
        typedef const matrix_reference<const_self_type> const_closure_type;
#endif
        typedef matrix_reference<self_type> closure_type;
        typedef packed_tag storage_category;
        typedef typename F1::packed_category packed_category;
        typedef typename F2::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        symmetric_matrix ():
            matrix_expression<self_type> (),
            size_ (0), data_ (0) {}
        BOOST_UBLAS_INLINE
        symmetric_matrix (size_type size):
            matrix_expression<self_type> (),
            size_ (BOOST_UBLAS_SAME (size, size)), data_ (0) {
            resize (size);
        }
        BOOST_UBLAS_INLINE
        symmetric_matrix (size_type size1, size_type size2):
            matrix_expression<self_type> (),
            size_ (BOOST_UBLAS_SAME (size1, size2)), data_ (0) {
            resize (size1, size2);
        }
        BOOST_UBLAS_INLINE
        symmetric_matrix (size_type size, const array_type &data):
            matrix_expression<self_type> (),
            size_ (size), data_ (data) {}
        BOOST_UBLAS_INLINE
        symmetric_matrix (const symmetric_matrix &m):
            matrix_expression<self_type> (),
            size_ (m.size_), data_ (m.data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_matrix (const matrix_expression<AE> &ae):
            matrix_expression<self_type> (),
            size_ (BOOST_UBLAS_SAME (ae ().size1 (), ae ().size2 ())), data_ (0) {
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
            return size_;
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const {
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
            size_ = BOOST_UBLAS_SAME (size, size);
            detail::resize (data (), functor1_type::packed_size (size, size), preserve);
        }
        BOOST_UBLAS_INLINE
        void resize (size_type size1, size_type size2, bool preserve = true) {
            size_ = BOOST_UBLAS_SAME (size1, size2);
            detail::resize (data (), functor1_type::packed_size (size1, size2), preserve);
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            BOOST_UBLAS_CHECK (j < size_, bad_index ());
            if (functor1_type::other (i, j))
                return data () [functor1_type::element (functor2_type (), i, size_, j, size_)];
            else
                return data () [functor1_type::element (functor2_type (), j, size_, i, size_)];
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            BOOST_UBLAS_CHECK (j < size_, bad_index ());
            if (functor1_type::other (i, j))
                return data () [functor1_type::element (functor2_type (), i, size_, j, size_)];
            else
                return data () [functor1_type::element (functor2_type (), j, size_, i, size_)];
        }

        // Assignment
        BOOST_UBLAS_INLINE
        symmetric_matrix &operator = (const symmetric_matrix &m) {
            // Precondition for container relaxed as requested during review.
            // BOOST_UBLAS_CHECK (size_ == m.size_, bad_size ());
            size_ = m.size_;
            data () = m.data ();
            return *this;
        }
        BOOST_UBLAS_INLINE
        symmetric_matrix &assign_temporary (symmetric_matrix &m) {
            swap (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_matrix &operator = (const matrix_expression<AE> &ae) {
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
        symmetric_matrix &reset (const matrix_expression<AE> &ae) {
            self_type temporary (ae);
            resize (temporary.size1 (), temporary.size2 (), false);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_matrix &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_matrix& operator += (const matrix_expression<AE> &ae) {
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
        symmetric_matrix &plus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_matrix& operator -= (const matrix_expression<AE> &ae) {
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
        symmetric_matrix &minus_assign (const matrix_expression<AE> &ae) { 
            matrix_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae); 
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        symmetric_matrix& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        symmetric_matrix& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (symmetric_matrix &m) {
            // Too unusual semantic
            // BOOST_UBLAS_CHECK (this != &m, external_logic ());
            if (this != &m) {
                // Precondition for container relaxed as requested during review.
                // BOOST_UBLAS_CHECK (size_ == m.size_, bad_size ());
                std::swap (size_, m.size_);
                data ().swap (m.data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (symmetric_matrix &m1, symmetric_matrix &m2) {
            m1.swap (m2);
        }
#endif

        // Element insertion and erasure
        // These functions should work with std::vector.
        // Thanks to Kresimir Fresl for spotting this.
        BOOST_UBLAS_INLINE
        void insert (size_type i, size_type j, const_reference t) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            BOOST_UBLAS_CHECK (j < size_, bad_index ());
// FIXME: is this ugly check still needed?!
// #ifndef BOOST_UBLAS_USE_ET
//             if (t == value_type ())
//                 return;
// #endif
            if (functor1_type::other (i, j)) {
                size_type k = functor1_type::element (functor2_type (), i, size_, j, size_);
                BOOST_UBLAS_CHECK (type_traits<value_type>::equals (data () [k], value_type ()) ||
                                   type_traits<value_type>::equals (data () [k], t), bad_index ());
                // data ().insert (data ().begin () + k, t);
                data () [k] = t;
            } else {
                size_type k = functor1_type::element (functor2_type (), j, size_, i, size_);
                BOOST_UBLAS_CHECK (type_traits<value_type>::equals (data () [k], value_type ()) ||
                                   type_traits<value_type>::equals (data () [k], t), bad_index ());
                // data ().insert (data ().begin () + k, t);
                data () [k] = t;
            }
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            BOOST_UBLAS_CHECK (j < size_, bad_index ());
            if (functor1_type::other (i, j)) {
                size_type k = functor1_type::element (functor2_type (), i, size_, j, size_);
                // data ().erase (data ().begin () + k));
                data () [functor1_type::element (functor2_type (), i, size_, j, size_)] = value_type ();
            } else {
                size_type k = functor1_type::element (functor2_type (), j, size_, i, size_);
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
        typedef indexed_const_iterator1<self_type, dense_random_access_iterator_tag> const_iterator1;
        typedef indexed_const_iterator2<self_type, dense_random_access_iterator_tag> const_iterator2;
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
        const_iterator1 find1 (int /* rank */, size_type i, size_type j) const {
            return const_iterator1 (*this, i, j);
        }
        BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j) {
            if (rank == 1)
                i = functor1_type::restrict1 (i, j);
            return iterator1 (*this, i, j);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int /* rank */, size_type i, size_type j) const {
            return const_iterator2 (*this, i, j);
        }
        BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j) {
            if (rank == 1)
                j = functor1_type::restrict2 (i, j);
            return iterator2 (*this, i, j);
        }

        // Iterators simply are indices.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<symmetric_matrix>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               const_iterator1, value_type> {
        public:
            typedef dense_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename symmetric_matrix::difference_type difference_type;
            typedef typename symmetric_matrix::value_type value_type;
            typedef typename symmetric_matrix::const_reference reference;
            typedef typename symmetric_matrix::const_pointer pointer;
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
            return find1 (0, size_, 0);
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator1:
            public container_reference<symmetric_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator1, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename symmetric_matrix::difference_type difference_type;
            typedef typename symmetric_matrix::value_type value_type;
            typedef typename symmetric_matrix::reference reference;
            typedef typename symmetric_matrix::pointer pointer;
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
            return find1 (0, size_, 0);
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator2:
            public container_const_reference<symmetric_matrix>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               const_iterator2, value_type> {
        public:
            typedef dense_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename symmetric_matrix::difference_type difference_type;
            typedef typename symmetric_matrix::value_type value_type;
            typedef typename symmetric_matrix::const_reference reference;
            typedef typename symmetric_matrix::const_pointer pointer;
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
            return find2 (0, 0, size_);
        }

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class iterator2:
            public container_reference<symmetric_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator2, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename symmetric_matrix::difference_type difference_type;
            typedef typename symmetric_matrix::value_type value_type;
            typedef typename symmetric_matrix::reference reference;
            typedef typename symmetric_matrix::pointer pointer;
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
            return find2 (0, 0, size_);
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
        size_type size_;
        array_type data_;
    };

    // Symmetric matrix adaptor class
    template<class M, class F>
    class symmetric_adaptor:
        public matrix_expression<symmetric_adaptor<M, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<symmetric_adaptor<M, F> >::operator ();
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
        typedef const symmetric_adaptor<M, F> const_self_type;
        typedef symmetric_adaptor<M, F> self_type;
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
        symmetric_adaptor ():
            matrix_expression<self_type> (),
            data_ (nil_) {
            BOOST_UBLAS_CHECK (data_.size1 () == data_.size2 (), bad_size ());
        }
        BOOST_UBLAS_INLINE
        symmetric_adaptor (matrix_type &data):
            matrix_expression<self_type> (),
            data_ (data) {
            BOOST_UBLAS_CHECK (data_.size1 () == data_.size2 (), bad_size ());
        }
        BOOST_UBLAS_INLINE
        symmetric_adaptor (const symmetric_adaptor &m):
            matrix_expression<self_type> (),
            data_ (m.data_) {
            BOOST_UBLAS_CHECK (data_.size1 () == data_.size2 (), bad_size ());
        }

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
            BOOST_UBLAS_CHECK (data.size1 () == data.size2 (), bad_size ());
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
            else
                return data () (j, i);
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
            if (functor_type::other (i, j))
                return data () (i, j);
            else
                return data () (j, i);
        }
#else
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
            if (functor_type::other (i, j))
                return data () (i, j);
            else
                return data () (j, i);
        }
#endif

        // Assignment
        BOOST_UBLAS_INLINE
        symmetric_adaptor &operator = (const symmetric_adaptor &m) {
            matrix_assign (scalar_assign<reference, value_type> (), *this, m, functor_type ());
            return *this;
        }
        BOOST_UBLAS_INLINE
        symmetric_adaptor &assign_temporary (symmetric_adaptor &m) {
            *this = m;
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_adaptor &operator = (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, value_type> (), *this, matrix<value_type> (ae), functor_type ());
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_adaptor &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae, functor_type ());
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_adaptor& operator += (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, value_type> (), *this, matrix<value_type> (*this + ae), functor_type ()); 
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_adaptor &plus_assign (const matrix_expression<AE> &ae) { 
            matrix_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae, functor_type ());
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_adaptor& operator -= (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, value_type> (), *this, matrix<value_type> (*this - ae), functor_type ()); 
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        symmetric_adaptor &minus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae, functor_type ()); 
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        symmetric_adaptor& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        symmetric_adaptor& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Comparison
        bool operator == (const symmetric_adaptor &sa) const {
            return (*this).data () == sa.data ();
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (symmetric_adaptor &m) {
            // Too unusual semantic
            // BOOST_UBLAS_CHECK (this != &m, external_logic ());
            if (this != &m)
                matrix_swap (scalar_swap<reference, reference> (), *this, m, functor_type ());
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (symmetric_adaptor &m1, symmetric_adaptor &m2) {
            m1.swap (m2);
        }
#endif

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_iterator1<self_type, packed_random_access_iterator_tag> iterator1;
        typedef indexed_iterator2<self_type, packed_random_access_iterator_tag> iterator2;
        typedef indexed_const_iterator1<self_type, dense_random_access_iterator_tag> const_iterator1;
        typedef indexed_const_iterator2<self_type, dense_random_access_iterator_tag> const_iterator2;
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
            if (functor_type::other (i, j)) {
                if (functor_type::other (size1 (), j)) {
                    return const_iterator1 (*this, 0, 0,
                                            data ().find1 (rank, i, j), data ().find1 (rank, size1 (), j),
                                            data ().find2 (rank, size2 (), size1 ()), data ().find2 (rank, size2 (), size1 ()));
                } else {
                    return const_iterator1 (*this, 0, 1,
                                            data ().find1 (rank, i, j), data ().find1 (rank, j, j),
                                            data ().find2 (rank, j, j), data ().find2 (rank, j, size1 ()));
                }
            } else {
                if (functor_type::other (size1 (), j)) {
                    return const_iterator1 (*this, 1, 0,
                                            data ().find1 (rank, j, j), data ().find1 (rank, size1 (), j),
                                            data ().find2 (rank, j, i), data ().find2 (rank, j, j));
                } else {
                    return const_iterator1 (*this, 1, 1,
                                            data ().find1 (rank, size1 (), size2 ()), data ().find1 (rank, size1 (), size2 ()),
                                            data ().find2 (rank, j, i), data ().find2 (rank, j, size1 ()));
                }
            }
        }
        BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j) {
            if (rank == 1)
                i = functor_type::restrict1 (i, j);
            return iterator1 (*this, data ().find1 (rank, i, j));
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            if (functor_type::other (i, j)) {
                if (functor_type::other (i, size2 ())) {
                    return const_iterator2 (*this, 1, 1,
                                            data ().find1 (rank, size2 (), size1 ()), data ().find1 (rank, size2 (), size1 ()),
                                            data ().find2 (rank, i, j), data ().find2 (rank, i, size2 ()));
                } else {
                    return const_iterator2 (*this, 1, 0,
                                            data ().find1 (rank, i, i), data ().find1 (rank, size2 (), i),
                                            data ().find2 (rank, i, j), data ().find2 (rank, i, i));
                }
            } else {
                if (functor_type::other (i, size2 ())) {
                    return const_iterator2 (*this, 0, 1,
                                            data ().find1 (rank, j, i), data ().find1 (rank, i, i),
                                            data ().find2 (rank, i, i), data ().find2 (rank, i, size2 ()));
                } else {
                    return const_iterator2 (*this, 0, 0,
                                            data ().find1 (rank, j, i), data ().find1 (rank, size2 (), i),
                                            data ().find2 (rank, size1 (), size2 ()), data ().find2 (rank, size2 (), size2 ()));
                }
            }
        }
        BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j) {
            if (rank == 1)
                j = functor_type::restrict2 (i, j);
            return iterator2 (*this, data ().find2 (rank, i, j));
        }

        // Iterators simply are indices.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<symmetric_adaptor>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               const_iterator1, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename const_iterator1_type::iterator_category,
                                                      dense_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename const_iterator1_type::difference_type difference_type;
            typedef typename const_iterator1_type::value_type value_type;
            typedef typename const_iterator1_type::reference reference;
            typedef typename const_iterator1_type::pointer pointer;
#else
            typedef typename iterator_restrict_traits<typename M::const_iterator1::iterator_category,
                                                      dense_random_access_iterator_tag>::iterator_category iterator_category;
            typedef const_reference reference;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (),
                begin_ (-1), end_ (-1), current_ (-1),
                it1_begin_ (), it1_end_ (), it1_ (),
                it2_begin_ (), it2_end_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &m, int begin, int end,
                             const const_iterator1_type &it1_begin, const const_iterator1_type &it1_end,
                             const const_iterator2_type &it2_begin, const const_iterator2_type &it2_end):
                container_const_reference<self_type> (m),
                begin_ (begin), end_ (end), current_ (begin),
                it1_begin_ (it1_begin), it1_end_ (it1_end), it1_ (it1_begin_),
                it2_begin_ (it2_begin), it2_end_ (it2_end), it2_ (it2_begin_) {
                if (current_ == 0 && it1_ == it1_end_)
                    current_ = 1;
                if (current_ == 1 && it2_ == it2_end_)
                    current_ = 0;
                if ((current_ == 0 && it1_ == it1_end_) ||
                    (current_ == 1 && it2_ == it2_end_))
                    current_ = end_;
                BOOST_UBLAS_CHECK (current_ == end_ ||
                                   (current_ == 0 && it1_ != it1_end_) ||
                                   (current_ == 1 && it2_ != it2_end_), internal_logic ());
            }
            BOOST_UBLAS_INLINE
            const_iterator1 (const iterator1 &it):
                container_const_reference<self_type> (it ()),
                begin_ (it.begin_), end_ (it.end_), current_ (it.current_),
                it1_begin_ (it.it1_begin_), it1_end_ (it.it1_end_), it1_ (it.it1_),
                it2_begin_ (it.it2_begin_), it2_end_ (it.it2_end_), it2_ (it.it2_) {
                BOOST_UBLAS_CHECK (current_ == end_ ||
                                   (current_ == 0 && it1_ != it1_end_) ||
                                   (current_ == 1 && it2_ != it2_end_), internal_logic ());
            }

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    BOOST_UBLAS_CHECK (it1_ != it1_end_, internal_logic ());
                    ++ it1_;
                    if (it1_ == it1_end_ && end_ == 1) {
                        it2_ = it2_begin_;
                        current_ = 1;
                    }
                } else /* if (current_ == 1) */ {
                    BOOST_UBLAS_CHECK (it2_ != it2_end_, internal_logic ());
                    ++ it2_;
                    if (it2_ == it2_end_ && end_ == 0) {
                        it1_ = it1_begin_;
                        current_ = 0;
                    }
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    if (it1_ == it1_begin_ && begin_ == 1) {
                        it2_ = it2_end_;
                        BOOST_UBLAS_CHECK (it2_ != it2_begin_, internal_logic ());
                        -- it2_;
                        current_ = 1;
                    } else {
                        -- it1_;
                    }
                } else /* if (current_ == 1) */ {
                    if (it2_ == it2_begin_ && begin_ == 0) {
                        it1_ = it1_end_;
                        BOOST_UBLAS_CHECK (it1_ != it1_begin_, internal_logic ());
                        -- it1_;
                        current_ = 0;
                    } else {
                        -- it2_;
                    }
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator += (difference_type n) {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    size_type d = std::min (n, it1_end_ - it1_);
                    it1_ += d;
                    n -= d;
                    if (n > 0 || (end_ == 1 && it1_ == it1_end_)) {
                        BOOST_UBLAS_CHECK (end_ == 1, external_logic ());
                        d = std::min (n, it2_end_ - it2_begin_);
                        it2_ = it2_begin_ + d;
                        n -= d;
                        current_ = 1;
                    }
                } else /* if (current_ == 1) */ {
                    size_type d = std::min (n, it2_end_ - it2_);
                    it2_ += d;
                    n -= d;
                    if (n > 0 || (end_ == 0 && it2_ == it2_end_)) {
                        BOOST_UBLAS_CHECK (end_ == 0, external_logic ());
                        d = std::min (n, it1_end_ - it1_begin_);
                        it1_ = it1_begin_ + d;
                        n -= d;
                        current_ = 0;
                    }
                }
                BOOST_UBLAS_CHECK (n == 0, external_logic ());
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -= (difference_type n) {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    size_type d = std::min (n, it1_ - it1_begin_);
                    it1_ -= d;
                    n -= d;
                    if (n > 0) {
                        BOOST_UBLAS_CHECK (end_ == 1, external_logic ());
                        d = std::min (n, it2_end_ - it2_begin_);
                        it2_ = it2_end_ - d;
                        n -= d;
                        current_ = 1;
                    }
                } else /* if (current_ == 1) */ {
                    size_type d = std::min (n, it2_ - it2_begin_);
                    it2_ -= d;
                    n -= d;
                    if (n > 0) {
                        BOOST_UBLAS_CHECK (end_ == 0, external_logic ());
                        d = std::min (n, it1_end_ - it1_begin_);
                        it1_ = it1_end_ - d;
                        n -= d;
                        current_ = 0;
                    }
                }
                BOOST_UBLAS_CHECK (n == 0, external_logic ());
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                BOOST_UBLAS_CHECK (it.current_ == 0 || it.current_ == 1, internal_logic ());
                BOOST_UBLAS_CHECK (/* begin_ == it.begin_ && */ end_ == it.end_, internal_logic ());
                if (current_ == 0 && it.current_ == 0) {
                    return it1_ - it.it1_;
                } else if (current_ == 0 && it.current_ == 1) {
                    if (end_ == 1 && it.end_ == 1) {
                        return (it1_ - it.it1_end_) + (it.it2_begin_ - it.it2_);
                    } else /* if (end_ == 0 && it.end_ == 0) */ {
                        return (it1_ - it.it1_begin_) + (it.it2_end_ - it.it2_);
                    }

                } else if (current_ == 1 && it.current_ == 0) {
                    if (end_ == 1 && it.end_ == 1) {
                        return (it2_ - it.it2_begin_) + (it.it1_end_ - it.it1_);
                    } else /* if (end_ == 0 && it.end_ == 0) */ {
                        return (it2_ - it.it2_end_) + (it.it1_begin_ - it.it1_);
                    }
                } else /* if (current_ == 1 && it.current_ == 1) */ {
                    return it2_ - it.it2_;
                }
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    BOOST_UBLAS_CHECK (it1_ != it1_end_, internal_logic ());
                    return *it1_;
                } else /* if (current_ == 1) */ {
                    BOOST_UBLAS_CHECK (it2_ != it2_end_, internal_logic ());
                    return *it2_;
                }
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
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    BOOST_UBLAS_CHECK (it1_ != it1_end_, internal_logic ());
                    return it1_.index1 ();
                } else /* if (current_ == 1) */ {
                    BOOST_UBLAS_CHECK (it2_ != it2_end_, internal_logic ());
                    return it2_.index2 ();
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    BOOST_UBLAS_CHECK (it1_ != it1_end_, internal_logic ());
                    return it1_.index2 ();
                } else /* if (current_ == 1) */ {
                    BOOST_UBLAS_CHECK (it2_ != it2_end_, internal_logic ());
                    return it2_.index1 ();
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                begin_ = it.begin_;
                end_ = it.end_;
                current_ = it.current_;
                it1_begin_ = it.it1_begin_;
                it1_end_ = it.it1_end_;
                it1_ = it.it1_;
                it2_begin_ = it.it2_begin_;
                it2_end_ = it.it2_end_;
                it2_ = it.it2_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                BOOST_UBLAS_CHECK (it.current_ == 0 || it.current_ == 1, internal_logic ());
                BOOST_UBLAS_CHECK (/* begin_ == it.begin_ && */ end_ == it.end_, internal_logic ());
                return (current_ == 0 && it.current_ == 0 && it1_ == it.it1_) ||
                       (current_ == 1 && it.current_ == 1 && it2_ == it.it2_);
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it - *this > 0;
            }

        private:
            int begin_;
            int end_;
            int current_;
            const_iterator1_type it1_begin_;
            const_iterator1_type it1_end_;
            const_iterator1_type it1_;
            const_iterator2_type it2_begin_;
            const_iterator2_type it2_end_;
            const_iterator2_type it2_;
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
            public container_reference<symmetric_adaptor>,
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
                return *it1_;
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
            public container_const_reference<symmetric_adaptor>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               const_iterator2, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename const_iterator2_type::iterator_category,
                                                      dense_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename const_iterator2_type::difference_type difference_type;
            typedef typename const_iterator2_type::value_type value_type;
            typedef typename const_iterator2_type::reference reference;
            typedef typename const_iterator2_type::pointer pointer;
#else
            typedef typename iterator_restrict_traits<typename M::const_iterator2::iterator_category,
                                                      dense_random_access_iterator_tag>::iterator_category iterator_category;
            typedef const_reference reference;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (),
                begin_ (-1), end_ (-1), current_ (-1),
                it1_begin_ (), it1_end_ (), it1_ (),
                it2_begin_ (), it2_end_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &m, int begin, int end,
                             const const_iterator1_type &it1_begin, const const_iterator1_type &it1_end,
                             const const_iterator2_type &it2_begin, const const_iterator2_type &it2_end):
                container_const_reference<self_type> (m),
                begin_ (begin), end_ (end), current_ (begin),
                it1_begin_ (it1_begin), it1_end_ (it1_end), it1_ (it1_begin_),
                it2_begin_ (it2_begin), it2_end_ (it2_end), it2_ (it2_begin_) {
                if (current_ == 0 && it1_ == it1_end_)
                    current_ = 1;
                if (current_ == 1 && it2_ == it2_end_)
                    current_ = 0;
                if ((current_ == 0 && it1_ == it1_end_) ||
                    (current_ == 1 && it2_ == it2_end_))
                    current_ = end_;
                BOOST_UBLAS_CHECK (current_ == end_ ||
                                   (current_ == 0 && it1_ != it1_end_) ||
                                   (current_ == 1 && it2_ != it2_end_), internal_logic ());
            }
            BOOST_UBLAS_INLINE
            const_iterator2 (const iterator2 &it):
                container_const_reference<self_type> (it ()),
                begin_ (it.begin_), end_ (it.end_), current_ (it.current_),
                it1_begin_ (it.it1_begin_), it1_end_ (it.it1_end_), it1_ (it.it1_),
                it2_begin_ (it.it2_begin_), it2_end_ (it.it2_end_), it2_ (it.it2_) {
                BOOST_UBLAS_CHECK (current_ == end_ ||
                                   (current_ == 0 && it1_ != it1_end_) ||
                                   (current_ == 1 && it2_ != it2_end_), internal_logic ());
            }

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    BOOST_UBLAS_CHECK (it1_ != it1_end_, internal_logic ());
                    ++ it1_;
                    if (it1_ == it1_end_ && end_ == 1) {
                        it2_ = it2_begin_;
                        current_ = 1;
                    }
                } else /* if (current_ == 1) */ {
                    BOOST_UBLAS_CHECK (it2_ != it2_end_, internal_logic ());
                    ++ it2_;
                    if (it2_ == it2_end_ && end_ == 0) {
                        it1_ = it1_begin_;
                        current_ = 0;
                    }
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    if (it1_ == it1_begin_ && begin_ == 1) {
                        it2_ = it2_end_;
                        BOOST_UBLAS_CHECK (it2_ != it2_begin_, internal_logic ());
                        -- it2_;
                        current_ = 1;
                    } else {
                        -- it1_;
                    }
                } else /* if (current_ == 1) */ {
                    if (it2_ == it2_begin_ && begin_ == 0) {
                        it1_ = it1_end_;
                        BOOST_UBLAS_CHECK (it1_ != it1_begin_, internal_logic ());
                        -- it1_;
                        current_ = 0;
                    } else {
                        -- it2_;
                    }
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator += (difference_type n) {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    size_type d = std::min (n, it1_end_ - it1_);
                    it1_ += d;
                    n -= d;
                    if (n > 0 || (end_ == 1 && it1_ == it1_end_)) {
                        BOOST_UBLAS_CHECK (end_ == 1, external_logic ());
                        d = std::min (n, it2_end_ - it2_begin_);
                        it2_ = it2_begin_ + d;
                        n -= d;
                        current_ = 1;
                    }
                } else /* if (current_ == 1) */ {
                    size_type d = std::min (n, it2_end_ - it2_);
                    it2_ += d;
                    n -= d;
                    if (n > 0 || (end_ == 0 && it2_ == it2_end_)) {
                        BOOST_UBLAS_CHECK (end_ == 0, external_logic ());
                        d = std::min (n, it1_end_ - it1_begin_);
                        it1_ = it1_begin_ + d;
                        n -= d;
                        current_ = 0;
                    }
                }
                BOOST_UBLAS_CHECK (n == 0, external_logic ());
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -= (difference_type n) {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    size_type d = std::min (n, it1_ - it1_begin_);
                    it1_ -= d;
                    n -= d;
                    if (n > 0) {
                        BOOST_UBLAS_CHECK (end_ == 1, external_logic ());
                        d = std::min (n, it2_end_ - it2_begin_);
                        it2_ = it2_end_ - d;
                        n -= d;
                        current_ = 1;
                    }
                } else /* if (current_ == 1) */ {
                    size_type d = std::min (n, it2_ - it2_begin_);
                    it2_ -= d;
                    n -= d;
                    if (n > 0) {
                        BOOST_UBLAS_CHECK (end_ == 0, external_logic ());
                        d = std::min (n, it1_end_ - it1_begin_);
                        it1_ = it1_end_ - d;
                        n -= d;
                        current_ = 0;
                    }
                }
                BOOST_UBLAS_CHECK (n == 0, external_logic ());
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                BOOST_UBLAS_CHECK (it.current_ == 0 || it.current_ == 1, internal_logic ());
                BOOST_UBLAS_CHECK (/* begin_ == it.begin_ && */ end_ == it.end_, internal_logic ());
                if (current_ == 0 && it.current_ == 0) {
                    return it1_ - it.it1_;
                } else if (current_ == 0 && it.current_ == 1) {
                    if (end_ == 1 && it.end_ == 1) {
                        return (it1_ - it.it1_end_) + (it.it2_begin_ - it.it2_);
                    } else /* if (end_ == 0 && it.end_ == 0) */ {
                        return (it1_ - it.it1_begin_) + (it.it2_end_ - it.it2_);
                    }

                } else if (current_ == 1 && it.current_ == 0) {
                    if (end_ == 1 && it.end_ == 1) {
                        return (it2_ - it.it2_begin_) + (it.it1_end_ - it.it1_);
                    } else /* if (end_ == 0 && it.end_ == 0) */ {
                        return (it2_ - it.it2_end_) + (it.it1_begin_ - it.it1_);
                    }
                } else /* if (current_ == 1 && it.current_ == 1) */ {
                    return it2_ - it.it2_;
                }
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    BOOST_UBLAS_CHECK (it1_ != it1_end_, internal_logic ());
                    return *it1_;
                } else /* if (current_ == 1) */ {
                    BOOST_UBLAS_CHECK (it2_ != it2_end_, internal_logic ());
                    return *it2_;
                }
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
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    BOOST_UBLAS_CHECK (it1_ != it1_end_, internal_logic ());
                    return it1_.index2 ();
                } else /* if (current_ == 1) */ {
                    BOOST_UBLAS_CHECK (it2_ != it2_end_, internal_logic ());
                    return it2_.index1 ();
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    BOOST_UBLAS_CHECK (it1_ != it1_end_, internal_logic ());
                    return it1_.index1 ();
                } else /* if (current_ == 1) */ {
                    BOOST_UBLAS_CHECK (it2_ != it2_end_, internal_logic ());
                    return it2_.index2 ();
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                begin_ = it.begin_;
                end_ = it.end_;
                current_ = it.current_;
                it1_begin_ = it.it1_begin_;
                it1_end_ = it.it1_end_;
                it1_ = it.it1_;
                it2_begin_ = it.it2_begin_;
                it2_end_ = it.it2_end_;
                it2_ = it.it2_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                BOOST_UBLAS_CHECK (it.current_ == 0 || it.current_ == 1, internal_logic ());
                BOOST_UBLAS_CHECK (/* begin_ == it.begin_ && */ end_ == it.end_, internal_logic ());
                return (current_ == 0 && it.current_ == 0 && it1_ == it.it1_) ||
                       (current_ == 1 && it.current_ == 1 && it2_ == it.it2_);
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it - *this > 0;
            }

        private:
            int begin_;
            int end_;
            int current_;
            const_iterator1_type it1_begin_;
            const_iterator1_type it1_end_;
            const_iterator1_type it1_;
            const_iterator2_type it2_begin_;
            const_iterator2_type it2_end_;
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
            public container_reference<symmetric_adaptor>,
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
                return *it2_;
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
    };

    template<class M, class F>
    typename symmetric_adaptor<M, F>::matrix_type symmetric_adaptor<M, F>::nil_;

}}}

#endif




