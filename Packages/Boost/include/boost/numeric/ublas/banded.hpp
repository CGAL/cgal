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

#ifndef BOOST_UBLAS_BANDED_H
#define BOOST_UBLAS_BANDED_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Iterators based on ideas of Jeremy Siek

namespace boost { namespace numeric { namespace ublas {

    // Array based banded matrix class
    template<class T, class F, class A>
    class banded_matrix:
        public matrix_expression<banded_matrix<T, F, A> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<banded_matrix<T, F, A> >::operator ();
#endif
        typedef typename A::size_type size_type;
        typedef typename A::difference_type difference_type;
        typedef T value_type;
        typedef const T &const_reference;
        typedef T &reference;
        typedef A array_type;
    private:
        typedef T *pointer;
        typedef F functor_type;
        typedef banded_matrix<T, F, A> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const matrix_const_reference<const self_type> const_closure_type;
#else
        typedef const matrix_reference<const self_type> const_closure_type;
#endif
        typedef matrix_reference<self_type> closure_type;
        typedef vector<T, A> vector_temporary_type;
        typedef matrix<T, F, A> matrix_temporary_type;  // general sub-matrix
        typedef packed_tag storage_category;
        typedef typename F::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        banded_matrix ():
            matrix_expression<self_type> (),
            size1_ (0), size2_ (0),
            lower_ (0), upper_ (0), data_ (0) {}
        BOOST_UBLAS_INLINE
        banded_matrix (size_type size1, size_type size2, size_type lower = 0, size_type upper = 0):
            matrix_expression<self_type> (),
            size1_ (size1), size2_ (size2),
            lower_ (lower), upper_ (upper), data_ ((std::max) (size1, size2) * (lower + 1 + upper)) {
        }
        BOOST_UBLAS_INLINE
        banded_matrix (size_type size1, size_type size2, size_type lower, size_type upper, const array_type &data):
            matrix_expression<self_type> (),
            size1_ (size1), size2_ (size2),
            lower_ (lower), upper_ (upper), data_ (data) {}
        BOOST_UBLAS_INLINE
        banded_matrix (const banded_matrix &m):
            matrix_expression<self_type> (),
            size1_ (m.size1_), size2_ (m.size2_),
            lower_ (m.lower_), upper_ (m.upper_), data_ (m.data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_matrix (const matrix_expression<AE> &ae, size_type lower = 0, size_type upper = 0):
            matrix_expression<self_type> (),
            size1_ (ae ().size1 ()), size2_ (ae ().size2 ()),
            lower_ (lower), upper_ (upper),
            data_ (std::max (size1_, size2_) * (lower_ + 1 + upper_)) {
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
        size_type lower () const {
            return lower_;
        }
        BOOST_UBLAS_INLINE
        size_type upper () const {
            return upper_;
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
        void resize (size_type size1, size_type size2, size_type lower = 0, size_type upper = 0, bool preserve = true) {
            size1_ = size1;
            size2_ = size2;
            lower_ = lower;
            upper_ = upper;
            if (preserve) {
                self_type temporary (size1, size2, lower, upper);
                // FIXME use matrix_resize_preserve on conformant compilers
                // detail::matrix_resize_preserve<functor_type> (*this, temporary, size_, size_);
                assign_temporary (temporary);
            }
            else
                data ().resize ((std::max) (size1, size2) * (lower + 1 + upper));
        }

        BOOST_UBLAS_INLINE
        void resize_packed_preserve (size_type size1, size_type size2, size_type lower = 0, size_type upper = 0) {
            size1_ = size1;
            size2_ = size2;
            lower_ = lower;
            upper_ = upper;
            data ().resize ((std::max) (size1, size2) * (lower + 1 + upper), value_type (0));
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size1_, bad_index ());
            BOOST_UBLAS_CHECK (j < size2_, bad_index ());
#ifdef BOOST_UBLAS_OWN_BANDED
            size_type k = (std::max) (i, j);
            size_type l = lower_ + j - i;
            if (k < (std::max) (size1_, size2_) &&
                l < lower_ + 1 + upper_)
                return data () [functor_type::element (k, (std::max) (size1_, size2_),
                                                       l, lower_ + 1 + upper_)];
#else
            size_type k = j;
            size_type l = upper_ + i - j;
            if (k < size2_ &&
                l < lower_ + 1 + upper_)
                return data () [functor_type::element (k, size2_,
                                                       l, lower_ + 1 + upper_)];
#endif
            return zero_;
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size1_, bad_index ());
            BOOST_UBLAS_CHECK (j < size2_, bad_index ());
#ifdef BOOST_UBLAS_OWN_BANDED
            size_type k = (std::max) (i, j);
            size_type l = lower_ + j - i;
            if (k < (std::max) (size1_, size2_) &&
                l < lower_ + 1 + upper_)
                return data () [functor_type::element (k, (std::max) (size1_, size2_),
                                                       l, lower_ + 1 + upper_)];
#else
            size_type k = j;
            size_type l = upper_ + i - j;
            if (k < size2_ &&
                l < lower_ + 1 + upper_)
                return data () [functor_type::element (k, size2_,
                                                       l, lower_ + 1 + upper_)];
#endif
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
            bad_index ().raise ();
#endif
            return const_cast<reference>(zero_);
        }

        // Assignment
        BOOST_UBLAS_INLINE
        banded_matrix &operator = (const banded_matrix &m) {
            size1_ = m.size1_;
            size2_ = m.size2_;
            lower_ = m.lower_;
            upper_ = m.upper_;
            data () = m.data ();
            return *this;
        }
        BOOST_UBLAS_INLINE
        banded_matrix &assign_temporary (banded_matrix &m) {
            swap (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_matrix &operator = (const matrix_expression<AE> &ae) {
            // return assign (self_type (ae, lower_, upper_));
            self_type temporary (ae, lower_, upper_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_matrix &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_matrix& operator += (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this + ae, lower_, upper_));
            self_type temporary (*this + ae, lower_, upper_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_matrix &plus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_matrix& operator -= (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this - ae, lower_, upper_));
            self_type temporary (*this - ae, lower_, upper_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_matrix &minus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        banded_matrix& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        banded_matrix& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (banded_matrix &m) {
            if (this != &m) {
                std::swap (size1_, m.size1_);
                std::swap (size2_, m.size2_);
                std::swap (lower_, m.lower_);
                std::swap (upper_, m.upper_);
                data ().swap (m.data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (banded_matrix &m1, banded_matrix &m2) {
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
#ifdef BOOST_UBLAS_OWN_BANDED
            size_type k = (std::max) (i, j);
            size_type l = lower_ + j - i;
            BOOST_UBLAS_CHECK (type_traits<value_type>::equals (data () [functor_type::element (k, std::max (size1_, size2_),
                                                                                                l, lower_ + 1 + upper_)], value_type (0)), bad_index ());
            // data ().insert (data ().begin () + functor_type::element (k, (std::max) (size1_, size2_),
            //                                                           l, lower_ + 1 + upper_), t);
            data () [functor_type::element (k, (std::max) (size1_, size2_),
                                            l, lower_ + 1 + upper_)] = t;
#else
            size_type k = j;
            size_type l = upper_ + i - j;
            BOOST_UBLAS_CHECK (type_traits<value_type>::equals (data () [functor_type::element (k, size2_,
                                                                                                l, lower_ + 1 + upper_)], value_type (0)), bad_index ());
            // data ().insert (data ().begin () + functor_type::element (k, size2_,
            //                                                           l, lower_ + 1 + upper_), t);
            data () [functor_type::element (k, size2_,
                                            l, lower_ + 1 + upper_)] = t;
#endif
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size1_, bad_index ());
            BOOST_UBLAS_CHECK (j < size2_, bad_index ());
#ifdef BOOST_UBLAS_OWN_BANDED
            size_type k = (std::max) (i, j);
            size_type l = lower_ + j - i;
            // data ().erase (data ().begin () + functor_type::element (k, (std::max) (size1_, size2_),
            //                                                         l, lower_ + 1 + upper_));
            data () [functor_type::element (k, (std::max) (size1_, size2_),
                                            l, lower_ + 1 + upper_)] = value_type (0);
#else
            size_type k = j;
            size_type l = upper_ + i - j;
            // data ().erase (data ().begin () + functor_type::element (k, size2_,
            //                                                          l, lower_ + 1 + upper_));
            data () [functor_type::element (k, size2_,
                                            l, lower_ + 1 + upper_)] = value_type (0);
#endif
        }
        BOOST_UBLAS_INLINE
        void clear () {
            // data ().clear ();
            std::fill (data ().begin (), data ().end (), value_type (0));
        }

        // Iterator types
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
            if (rank == 1) {
                size_type lower_i = (std::max) (difference_type (j - upper_), difference_type (0));
                i = (std::max) (i, lower_i);
                size_type upper_i = (std::min) (j + 1 + lower_, size1_);
                i = (std::min) (i, upper_i);
            }
            return const_iterator1 (*this, i, j);
        }
        BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j) {
            if (rank == 1) {
                size_type lower_i = (std::max) (difference_type (j - upper_), difference_type (0));
                i = (std::max) (i, lower_i);
                size_type upper_i = (std::min) (j + 1 + lower_, size1_);
                i = (std::min) (i, upper_i);
            }
            return iterator1 (*this, i, j);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            if (rank == 1) {
                size_type lower_j = (std::max) (difference_type (i - lower_), difference_type (0));
                j = (std::max) (j, lower_j);
                size_type upper_j = (std::min) (i + 1 + upper_, size2_);
                j = (std::min) (j, upper_j);
            }
            return const_iterator2 (*this, i, j);
        }
        BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j) {
            if (rank == 1) {
                size_type lower_j = (std::max) (difference_type (i - lower_), difference_type (0));
                j = (std::max) (j, lower_j);
                size_type upper_j = (std::min) (i + 1 + upper_, size2_);
                j = (std::min) (j, upper_j);
            }
            return iterator2 (*this, i, j);
        }

        // Iterators simply are indices.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<banded_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator1, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename banded_matrix::value_type value_type;
            typedef typename banded_matrix::difference_type difference_type;
            typedef typename banded_matrix::const_reference reference;
            typedef const typename banded_matrix::pointer pointer;
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
            const_reference operator * () const {
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
            public container_reference<banded_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator1, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename banded_matrix::value_type value_type;
            typedef typename banded_matrix::difference_type difference_type;
            typedef typename banded_matrix::reference reference;
            typedef typename banded_matrix::pointer pointer;
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
            public container_const_reference<banded_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator2, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename banded_matrix::value_type value_type;
            typedef typename banded_matrix::difference_type difference_type;
            typedef typename banded_matrix::const_reference reference;
            typedef const typename banded_matrix::pointer pointer;
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
            const_reference operator * () const {
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
            public container_reference<banded_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator2, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename banded_matrix::value_type value_type;
            typedef typename banded_matrix::difference_type difference_type;
            typedef typename banded_matrix::reference reference;
            typedef typename banded_matrix::pointer pointer;
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
        size_type lower_;
        size_type upper_;
        array_type data_;
        typedef const value_type const_value_type;
        static const_value_type zero_;
    };

    template<class T, class F, class A>
    typename banded_matrix<T, F, A>::const_value_type banded_matrix<T, F, A>::zero_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME banded_matrix<T>::const_value_type
#endif
        (0);

    // Diagonal matrix class
    template<class T, class F, class A>
    class diagonal_matrix:
        public banded_matrix<T, F, A> {
    public:
        typedef typename A::size_type size_type;
        typedef banded_matrix<T, F, A> matrix_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        diagonal_matrix ():
            matrix_type () {}
        BOOST_UBLAS_INLINE
        diagonal_matrix (size_type size):
            matrix_type (size, size) {}
        BOOST_UBLAS_INLINE
        diagonal_matrix (size_type size1, size_type size2):
            matrix_type (size1, size2) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        diagonal_matrix (const matrix_expression<AE> &ae):
            matrix_type (ae) {}
        BOOST_UBLAS_INLINE
        ~diagonal_matrix () {}

        // Assignment
        BOOST_UBLAS_INLINE
        diagonal_matrix &operator = (const diagonal_matrix &m) {
            matrix_type::operator = (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        diagonal_matrix &operator = (const matrix_expression<AE> &ae) {
            matrix_type::operator = (ae);
            return *this;
        }
    };

    // Banded matrix adaptor class
    template<class M>
    class banded_adaptor:
        public matrix_expression<banded_adaptor<M> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<banded_adaptor<M> >::operator ();
#endif
        typedef const M const_matrix_type;
        typedef M matrix_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
#ifndef BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
        typedef typename M::const_reference const_reference;
        typedef typename M::reference reference;
#else
        typedef typename M::const_reference const_reference;
        typedef typename boost::mpl::if_<boost::is_const<M>,
                                          typename M::const_reference,
                                          typename M::reference>::type reference;
#endif
#ifndef BOOST_UBLAS_CT_PROXY_CLOSURE_TYPEDEFS
        typedef typename M::closure_type matrix_closure_type;
#else
        typedef typename boost::mpl::if_<boost::is_const<M>,
                                          typename M::const_closure_type,
                                          typename M::closure_type>::type matrix_closure_type;
#endif
    private:
        typedef banded_adaptor<M> self_type;
    public:
        typedef const self_type const_closure_type;
        typedef self_type closure_type;
        typedef typename M::vector_temporary_type vector_temporary_type;
        typedef typename M::matrix_temporary_type matrix_temporary_type;
        typedef typename storage_restrict_traits<typename M::storage_category,
                                                 packed_proxy_tag>::storage_category storage_category;
        typedef typename M::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        banded_adaptor ():
            matrix_expression<self_type> (),
            data_ (nil_), lower_ (0), upper_ (0) {}
        BOOST_UBLAS_INLINE
        banded_adaptor (matrix_type &data, size_type lower = 0, size_type upper = 0):
            matrix_expression<self_type> (),
            data_ (data), lower_ (lower), upper_ (upper) {}
        BOOST_UBLAS_INLINE
        banded_adaptor (const banded_adaptor &m):
            matrix_expression<self_type> (),
            data_ (m.data_), lower_ (m.lower_), upper_ (m.upper_) {}

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
        size_type lower () const {
            return lower_;
        }
        BOOST_UBLAS_INLINE
        size_type upper () const {
            return upper_;
        }
        BOOST_UBLAS_INLINE
        const matrix_closure_type &data () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        matrix_closure_type &data () {
            return data_;
        }

        // Element access
#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
#ifdef BOOST_UBLAS_OWN_BANDED
            size_type k = (std::max) (i, j);
            size_type l = lower_ + j - i;
            if (k < (std::max) (size1 (), size2 ()) &&
                l < lower_ + 1 + upper_)
                return data () (i, j);
#else
            size_type k = j;
            size_type l = upper_ + i - j;
            if (k < size2 () &&
                l < lower_ + 1 + upper_)
                return data () (i, j);
#endif
            return zero_;
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
#ifdef BOOST_UBLAS_OWN_BANDED
            size_type k = (std::max) (i, j);
            size_type l = lower_ + j - i;
            if (k < (std::max) (size1 (), size2 ()) &&
                l < lower_ + 1 + upper_)
                return data () (i, j);
#else
            size_type k = j;
            size_type l = upper_ + i - j;
            if (k < size2 () &&
                l < lower_ + 1 + upper_)
                return data () (i, j);
#endif
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
            bad_index ().raise ();
#endif
            return const_cast<reference>(zero_);
        }
#else
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
#ifdef BOOST_UBLAS_OWN_BANDED
            size_type k = (std::max) (i, j);
            size_type l = lower_ + j - i;
            if (k < (std::max) (size1 (), size2 ()) &&
                l < lower_ + 1 + upper_)
                return data () (i, j);
#else
            size_type k = j;
            size_type l = upper_ + i - j;
            if (k < size2 () &&
                l < lower_ + 1 + upper_)
                return data () (i, j);
#endif
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
            bad_index ().raise ();
#endif
            return const_cast<reference>(zero_);
        }
#endif

        // Assignment
        BOOST_UBLAS_INLINE
        banded_adaptor &operator = (const banded_adaptor &m) {
            matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, value_type> (), *this, m);
            return *this;
        }
        BOOST_UBLAS_INLINE
        banded_adaptor &assign_temporary (banded_adaptor &m) {
            *this = m;
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_adaptor &operator = (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, value_type> (), *this, matrix<value_type> (ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_adaptor &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_adaptor& operator += (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, value_type> (), *this, matrix<value_type> (*this + ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_adaptor &plus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_plus_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_adaptor& operator -= (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, value_type> (), *this, matrix<value_type> (*this - ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        banded_adaptor &minus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_minus_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        banded_adaptor& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        banded_adaptor& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, AT> (), *this, at);
            return *this;
        }

        // Closure comparison
        BOOST_UBLAS_INLINE
        bool same_closure (const banded_adaptor &ba) const {
            return (*this).data ().same_closure (ba.data ());
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (banded_adaptor &m) {
            if (this != &m) {
                BOOST_UBLAS_CHECK (lower_ == m.lower_, bad_size ());
                BOOST_UBLAS_CHECK (upper_ == m.upper_, bad_size ());
                matrix_swap (scalar_swap<BOOST_UBLAS_TYPENAME iterator1_type::reference, BOOST_UBLAS_TYPENAME iterator1_type::reference> (), *this, m);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (banded_adaptor &m1, banded_adaptor &m2) {
            m1.swap (m2);
        }
#endif

        // Iterator types
    private:
        // Use the matrix iterator
#ifndef BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
        typedef typename M::const_iterator1 const_iterator1_type;
        typedef typename M::iterator1 iterator1_type;
        typedef typename M::const_iterator2 const_iterator2_type;
        typedef typename M::iterator2 iterator2_type;
#else
        typedef typename M::const_iterator1 const_iterator1_type;
        typedef typename boost::mpl::if_<boost::is_const<M>,
                                          typename M::const_iterator1,
                                          typename M::iterator1>::type iterator1_type;
        typedef typename M::const_iterator2 const_iterator2_type;
        typedef typename boost::mpl::if_<boost::is_const<M>,
                                          typename M::const_iterator2,
                                          typename M::iterator2>::type iterator2_type;
#endif

    public:
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
            if (rank == 1) {
                size_type lower_i = (std::max) (difference_type (j - upper_), difference_type (0));
                i = (std::max) (i, lower_i);
                size_type upper_i = (std::min) (j + 1 + lower_, size1 ());
                i = (std::min) (i, upper_i);
            }
            return const_iterator1 (*this, data ().find1 (rank, i, j));
        }
        BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j) {
            if (rank == 1) {
                size_type lower_i = (std::max) (difference_type (j - upper_), difference_type (0));
                i = (std::max) (i, lower_i);
                size_type upper_i = (std::min) (j + 1 + lower_, size1 ());
                i = (std::min) (i, upper_i);
            }
            return iterator1 (*this, data ().find1 (rank, i, j));
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            if (rank == 1) {
                size_type lower_j = (std::max) (difference_type (i - lower_), difference_type (0));
                j = (std::max) (j, lower_j);
                size_type upper_j = (std::min) (i + 1 + upper_, size2 ());
                j = (std::min) (j, upper_j);
            }
            return const_iterator2 (*this, data ().find2 (rank, i, j));
        }
        BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j) {
            if (rank == 1) {
                size_type lower_j = (std::max) (difference_type (i - lower_), difference_type (0));
                j = (std::max) (j, lower_j);
                size_type upper_j = (std::min) (i + 1 + upper_, size2 ());
                j = (std::min) (j, upper_j);
            }
            return iterator2 (*this, data ().find2 (rank, i, j));
        }

        // Iterators simply are indices.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<banded_adaptor>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator1, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename const_iterator1_type::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename const_iterator1_type::value_type value_type;
            typedef typename const_iterator1_type::difference_type difference_type;
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
            const_reference operator * () const {
                size_type i = index1 ();
                size_type j = index2 ();
                BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
#ifdef BOOST_UBLAS_OWN_BANDED
                size_type k = (std::max) (i, j);
                size_type l = (*this) ().lower () + j - i;
                if (k < (std::max) ((*this) ().size1 (), (*this) ().size2 ()) &&
                    l < (*this) ().lower () + 1 + (*this) ().upper ())
                    return *it1_;
#else
                size_type k = j;
                size_type l = (*this) ().upper () + i - j;
                if (k < (*this) ().size2 () &&
                    l < (*this) ().lower () + 1 + (*this) ().upper ())
                    return *it1_;
#endif
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
            public container_reference<banded_adaptor>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator1, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename iterator1_type::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename iterator1_type::value_type value_type;
            typedef typename iterator1_type::difference_type difference_type;
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
#ifdef BOOST_UBLAS_OWN_BANDED
                size_type k = (std::max) (i, j);
                size_type l = (*this) ().lower () + j - i;
                if (k < (std::max) ((*this) ().size1 (), (*this) ().size2 ()) &&
                    l < (*this) ().lower () + 1 + (*this) ().upper ())
                    return *it1_;
#else
                size_type k = j;
                size_type l = (*this) ().upper () + i - j;
                if (k < (*this) ().size2 () &&
                    l < (*this) ().lower () + 1 + (*this) ().upper ())
                    return *it1_;
#endif
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
            public container_const_reference<banded_adaptor>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator2, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename const_iterator2_type::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename const_iterator2_type::value_type value_type;
            typedef typename const_iterator2_type::difference_type difference_type;
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
            const_reference operator * () const {
                size_type i = index1 ();
                size_type j = index2 ();
                BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
#ifdef BOOST_UBLAS_OWN_BANDED
                size_type k = (std::max) (i, j);
                size_type l = (*this) ().lower () + j - i;
                if (k < (std::max) ((*this) ().size1 (), (*this) ().size2 ()) &&
                    l < (*this) ().lower () + 1 + (*this) ().upper ())
                    return *it2_;
#else
                size_type k = j;
                size_type l = (*this) ().upper () + i - j;
                if (k < (*this) ().size2 () &&
                    l < (*this) ().lower () + 1 + (*this) ().upper ())
                    return *it2_;
#endif
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
            public container_reference<banded_adaptor>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator2, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename iterator2_type::iterator_category,
                                                      packed_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename iterator2_type::value_type value_type;
            typedef typename iterator2_type::difference_type difference_type;
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
#ifdef BOOST_UBLAS_OWN_BANDED
                size_type k = (std::max) (i, j);
                size_type l = (*this) ().lower () + j - i;
                if (k < (std::max) ((*this) ().size1 (), (*this) ().size2 ()) &&
                    l < (*this) ().lower () + 1 + (*this) ().upper ())
                    return *it2_;
#else
                size_type k = j;
                size_type l = (*this) ().upper () + i - j;
                if (k < (*this) ().size2 () &&
                    l < (*this) ().lower () + 1 + (*this) ().upper ())
                    return *it2_;
#endif
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
        size_type lower_;
        size_type upper_;
        static matrix_type nil_;
        typedef const value_type const_value_type;
        static const_value_type zero_;
    };

    template<class M>
    typename banded_adaptor<M>::matrix_type banded_adaptor<M>::nil_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME banded_adaptor<M>::matrix_type ()
#endif
    ;
    template<class M>
    typename banded_adaptor<M>::const_value_type banded_adaptor<M>::zero_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME banded_adaptor<M>::const_value_type
#endif
        (0);

    // Diagonal matrix adaptor class
    template<class M>
    class diagonal_adaptor:
        public banded_adaptor<M> {
    public:
        typedef M matrix_type;
        typedef banded_adaptor<M> adaptor_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        diagonal_adaptor ():
            adaptor_type () {}
        BOOST_UBLAS_INLINE
        diagonal_adaptor (matrix_type &data):
            adaptor_type (data) {}
        BOOST_UBLAS_INLINE
        ~diagonal_adaptor () {}

        // Assignment
        BOOST_UBLAS_INLINE
        diagonal_adaptor &operator = (const diagonal_adaptor &m) {
            adaptor_type::operator = (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        diagonal_adaptor &operator = (const matrix_expression<AE> &ae) {
            adaptor_type::operator = (ae);
            return *this;
        }
    };

}}}

#endif
