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

#ifndef BOOST_UBLAS_MATRIX_SPARSE_H
#define BOOST_UBLAS_MATRIX_SPARSE_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/storage_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Iterators based on ideas of Jeremy Siek

namespace boost { namespace numeric { namespace ublas {

#ifdef BOOST_UBLAS_STRICT_MATRIX_SPARSE

    template<class M>
    class sparse_matrix_element:
       public container_reference<M> {
    public:
        typedef M matrix_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;
        typedef const value_type &const_reference;
        typedef value_type *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        sparse_matrix_element (matrix_type &m, size_type i, size_type j):
            container_reference<matrix_type> (m), i_ (i), j_ (j), d_ (), dirty_ (false) {
            pointer it = (*this) ().find_element (i_, j_);
            if (it)
                d_ = *it;
        }
        BOOST_UBLAS_INLINE
        sparse_matrix_element (const sparse_matrix_element &p):
            container_reference<matrix_type> (p), i_ (p.i_), d_ (p.d_), dirty_ (p.dirty_) {}
        BOOST_UBLAS_INLINE
        ~sparse_matrix_element () {
            if (dirty_) {
                pointer it = (*this) ().find_element (i_, j_);
                if (! it)
                    (*this) ().insert (i_, j_, d_);
                else
                    *it = d_;
            }
        }

        // Assignment
        BOOST_UBLAS_INLINE
        sparse_matrix_element &operator = (const sparse_matrix_element &p) {
            // Overide the implict copy assignment
            d_ = p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_matrix_element &operator = (const D &d) {
            d_ = d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_matrix_element &operator += (const D &d) {
            d_ += d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_matrix_element &operator -= (const D &d) {
            d_ -= d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_matrix_element &operator *= (const D &d) {
            d_ *= d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_matrix_element &operator /= (const D &d) {
            d_ /= d;
            dirty_ = true;
            return *this;
        }

        // Comparison
        template<class D>
        BOOST_UBLAS_INLINE
        bool operator == (const D &d) const {
            return d_ == d;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        bool operator != (const D &d) const {
            return d_ != d;
        }


        // Conversion
        BOOST_UBLAS_INLINE
        operator const_reference () const {
            return d_;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (sparse_matrix_element p) {
            if (this != &p) {
                dirty_ = true;
                p.dirty_ = true;
                std::swap (d_, p.d_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (sparse_matrix_element p1, sparse_matrix_element p2) {
            p1.swap (p2);
        }
#endif

    private:
        size_type i_;
        size_type j_;
        value_type d_;
        bool dirty_;
    };

    template<class M>
    struct type_traits<sparse_matrix_element<M> > {
        typedef typename M::value_type element_type;
        typedef type_traits<sparse_matrix_element<M> > self_type;
        typedef typename type_traits<element_type>::value_type value_type;
        typedef typename type_traits<element_type>::const_reference const_reference;
        typedef sparse_matrix_element<M> reference;
        typedef typename type_traits<element_type>::real_type real_type;
        typedef typename type_traits<element_type>::precision_type precision_type;

        BOOST_STATIC_CONSTANT (unsigned, plus_complexity = type_traits<element_type>::plus_complexity);
        BOOST_STATIC_CONSTANT (unsigned, multiplies_complexity = type_traits<element_type>::multiplies_complexity);

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

    template<class M1, class T2>
    struct promote_traits<sparse_matrix_element<M1>, T2> {
        typedef typename promote_traits<typename sparse_matrix_element<M1>::value_type, T2>::promote_type promote_type;
    };
    template<class T1, class M2>
    struct promote_traits<T1, sparse_matrix_element<M2> > {
        typedef typename promote_traits<T1, typename sparse_matrix_element<M2>::value_type>::promote_type promote_type;
    };
    template<class M1, class M2>
    struct promote_traits<sparse_matrix_element<M1>, sparse_matrix_element<M2> > {
        typedef typename promote_traits<typename sparse_matrix_element<M1>::value_type,
                                        typename sparse_matrix_element<M2>::value_type>::promote_type promote_type;
    };

#endif


    // Array based sparse matrix class
    template<class T, class F, class A>
    class sparse_matrix:
        public matrix_expression<sparse_matrix<T, F, A> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<sparse_matrix<T, F, A> >::operator ();
#endif
        typedef typename A::size_type size_type;
        typedef typename A::difference_type difference_type;
        typedef T value_type;
        typedef A array_type;
        typedef const T &const_reference;
#ifndef BOOST_UBLAS_STRICT_MATRIX_SPARSE
        typedef BOOST_UBLAS_TYPENAME detail::map_traits<A, T>::reference reference;
#else
        typedef sparse_matrix_element<sparse_matrix<T, F, A> > reference;
#endif
    private:
        typedef T &true_reference;
        typedef T *pointer;
        typedef F functor_type;
        typedef sparse_matrix<T, F, A> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const matrix_const_reference<const self_type> const_closure_type;
#else
        typedef const matrix_reference<const self_type> const_closure_type;
#endif
        typedef matrix_reference<self_type> closure_type;
        typedef sparse_vector<T, A> vector_temporary_type;
        typedef self_type matrix_temporary_type;
        typedef sparse_tag storage_category;
        typedef typename F::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        sparse_matrix ():
            matrix_expression<self_type> (),
            size1_ (0), size2_ (0), data_ () {}
        BOOST_UBLAS_INLINE
        sparse_matrix (size_type size1, size_type size2, size_type non_zeros = 0):
            matrix_expression<self_type> (),
            size1_ (size1), size2_ (size2), data_ () {
            detail::map_reserve (data (), restrict_nz (non_zeros));
        }
        BOOST_UBLAS_INLINE
        sparse_matrix (const sparse_matrix &m):
            matrix_expression<self_type> (),
            size1_ (m.size1_), size2_ (m.size2_), data_ (m.data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_matrix (const matrix_expression<AE> &ae, size_type non_zeros = 0):
            matrix_expression<self_type> (),
            size1_ (ae ().size1 ()), size2_ (ae ().size2 ()), data_ () {
            detail::map_reserve (data (), restrict_nz (non_zeros));
            matrix_assign (scalar_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
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
        size_type non_zeros () const {
            return detail::map_capacity (data ());
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
    private:
        BOOST_UBLAS_INLINE
        size_type restrict_nz (size_type non_zeros) const {
            // Guarding against overflow - Thanks to Alexei Novakov for the hint.
            // non_zeros_ = (std::min) (non_zeros, size1_ * size2_);
            if (size1_ > 0 && non_zeros / size1_ >= size2_)
                non_zeros = size1_ * size2_;
            return non_zeros;
        }
    public:
        BOOST_UBLAS_INLINE
        void resize (size_type size1, size_type size2, bool preserve = true) {
            // FIXME preserve unimplemented
            BOOST_UBLAS_CHECK (!preserve, internal_logic ());
            size1_ = size1;
            size2_ = size2;
            data ().clear ();
        }

        // Reserving
        BOOST_UBLAS_INLINE
        void reserve (size_type non_zeros, bool preserve = true) {
            detail::map_reserve (data (), restrict_nz (non_zeros));
        }

        // Proxy support
#ifdef BOOST_UBLAS_STRICT_MATRIX_SPARSE
        BOOST_UBLAS_INLINE
        pointer find_element (size_type i, size_type j) {
            iterator_type it (data ().find (functor_type::element (i, size1_, j, size2_)));
            if (it == data ().end () || (*it).first != functor_type::element (i, size1_, j, size2_))
                return 0;
            return &(*it).second;
        }
#endif

        // Element access
        BOOST_UBLAS_INLINE
        const_reference at_element (size_type i, size_type j) const {
            const_iterator_type it (data ().find (functor_type::element (i, size1_, j, size2_)));
            if (it == data ().end () || (*it).first != functor_type::element (i, size1_, j, size2_))
                return zero_;
            return (*it).second;
        }
        BOOST_UBLAS_INLINE
        true_reference at_element (size_type i, size_type j) {
            return data () [functor_type::element (i, size1_, j, size2_)];
        }
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return at_element (i, j);
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
#ifndef BOOST_UBLAS_STRICT_MATRIX_SPARSE
            return at_element (i, j);
#else
            return reference (*this, i, j);
#endif
        }

        // Assignment
        BOOST_UBLAS_INLINE
        sparse_matrix &operator = (const sparse_matrix &m) {
            if (this != &m) {
                size1_ = m.size1_;
                size2_ = m.size2_;
                data () = m.data ();
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_matrix &assign_temporary (sparse_matrix &m) {
            swap (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_matrix &operator = (const matrix_expression<AE> &ae) {
            self_type temporary (ae, detail::map_capacity (data ()));
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_matrix &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_matrix& operator += (const matrix_expression<AE> &ae) {
            self_type temporary (*this + ae, detail::map_capacity (data ()));
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_matrix &plus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_plus_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_matrix& operator -= (const matrix_expression<AE> &ae) {
            self_type temporary (*this - ae, detail::map_capacity (data ()));
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_matrix &minus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_minus_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        sparse_matrix& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<true_reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        sparse_matrix& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<true_reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (sparse_matrix &m) {
            if (this != &m) {
                std::swap (size1_, m.size1_);
                std::swap (size2_, m.size2_);
                data ().swap (m.data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (sparse_matrix &m1, sparse_matrix &m2) {
            m1.swap (m2);
        }
#endif

        // Element insertion and erasure
        BOOST_UBLAS_INLINE
        void insert (size_type i, size_type j, const_reference t) {
            BOOST_UBLAS_CHECK (data ().find (functor_type::element (i, size1_, j, size2_)) == data ().end (), bad_index ());
            data ().insert (data ().end (), BOOST_UBLAS_TYPENAME array_type::value_type (functor_type::element (i, size1_, j, size2_), t));
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i, size_type j) {
            // FIXME: shouldn't we use const_iterator_type here?
            iterator_type it = data ().find (functor_type::element (i, size1_, j, size2_));
            if (it == data ().end ())
                return;
            data ().erase (it);
        }
        BOOST_UBLAS_INLINE
        void clear () {
            data ().clear ();
        }

        // Iterator types
    private:
        // Use storage iterator
        typedef typename A::const_iterator const_iterator_type;
        typedef typename A::iterator iterator_type;

    public:
        class const_iterator1;
        class iterator1;
        class const_iterator2;
        class iterator2;
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
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j, int direction = 1) const {
            const_iterator_type it (data ().lower_bound (functor_type::address (i, size1_, j, size2_)));
            const_iterator_type it_end (data ().end ());
            size_type index1 = size_type (-1);
            size_type index2 = size_type (-1);
            while (rank == 1 && it != it_end) {
                index1 = functor_type::index1 ((*it).first, size1_, size2_);
                index2 = functor_type::index2 ((*it).first, size1_, size2_);
                if (direction > 0) {
                    if ((index1 >= i && index2 == j) || (i >= size1_))
                        break;
                    ++ i;
                } else /* if (direction < 0) */ {
                    if ((index1 <= i && index2 == j) || (i == 0))
                        break;
                    -- i;
                }
                it = data ().lower_bound (functor_type::address (i, size1_, j, size2_));
            }
            if (rank == 1 && index2 != j) {
                if (direction > 0)
                    i = size1_;
                else /* if (direction < 0) */
                    i = 0;
                rank = 0;
            }
            return const_iterator1 (*this, rank, i, j, it);
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j, int direction = 1) {
            iterator_type it (data ().lower_bound (functor_type::address (i, size1_, j, size2_)));
            iterator_type it_end (data ().end ());
            size_type index1 = size_type (-1);
            size_type index2 = size_type (-1);
            while (rank == 1 && it != it_end) {
                index1 = functor_type::index1 ((*it).first, size1_, size2_);
                index2 = functor_type::index2 ((*it).first, size1_, size2_);
                if (direction > 0) {
                    if ((index1 >= i && index2 == j) || (i >= size1_))
                        break;
                    ++ i;
                } else /* if (direction < 0) */ {
                    if ((index1 <= i && index2 == j) || (i == 0))
                        break;
                    -- i;
                }
                it = data ().lower_bound (functor_type::address (i, size1_, j, size2_));
            }
            if (rank == 1 && index2 != j) {
                if (direction > 0)
                    i = size1_;
                else /* if (direction < 0) */
                    i = 0;
                rank = 0;
            }
            return iterator1 (*this, rank, i, j, it);
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j, int direction = 1) const {
            const_iterator_type it (data ().lower_bound (functor_type::address (i, size1_, j, size2_)));
            const_iterator_type it_end (data ().end ());
            size_type index1 = size_type (-1);
            size_type index2 = size_type (-1);
            while (rank == 1 && it != it_end) {
                index1 = functor_type::index1 ((*it).first, size1_, size2_);
                index2 = functor_type::index2 ((*it).first, size1_, size2_);
                if (direction > 0) {
                    if ((index2 >= j && index1 == i) || (j >= size2_))
                        break;
                    ++ j;
                } else /* if (direction < 0) */ {
                    if ((index2 <= j && index1 == i) || (j == 0))
                        break;
                    -- j;
                }
                it = data ().lower_bound (functor_type::address (i, size1_, j, size2_));
            }
            if (rank == 1 && index1 != i) {
                if (direction > 0)
                    j = size2_;
                else /* if (direction < 0) */
                    j = 0;
                rank = 0;
            }
            return const_iterator2 (*this, rank, i, j, it);
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j, int direction = 1) {
            iterator_type it (data ().lower_bound (functor_type::address (i, size1_, j, size2_)));
            iterator_type it_end (data ().end ());
            size_type index1 = size_type (-1);
            size_type index2 = size_type (-1);
            while (rank == 1 && it != it_end) {
                index1 = functor_type::index1 ((*it).first, size1_, size2_);
                index2 = functor_type::index2 ((*it).first, size1_, size2_);
                if (direction > 0) {
                    if ((index2 >= j && index1 == i) || (j >= size2_))
                        break;
                    ++ j;
                } else /* if (direction < 0) */ {
                    if ((index2 <= j && index1 == i) || (j == 0))
                        break;
                    -- j;
                }
                it = data ().lower_bound (functor_type::address (i, size1_, j, size2_));
            }
            if (rank == 1 && index1 != i) {
                if (direction > 0)
                    j = size2_;
                else /* if (direction < 0) */
                    j = 0;
                rank = 0;
            }
            return iterator2 (*this, rank, i, j, it);
        }


        class const_iterator1:
            public container_const_reference<sparse_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator1, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename sparse_matrix::value_type value_type;
            typedef typename sparse_matrix::difference_type difference_type;
            typedef typename sparse_matrix::const_reference reference;
            typedef const typename sparse_matrix::pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), rank_ (), i_ (), j_ (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &m, int rank, size_type i, size_type j, const const_iterator_type &it):
                container_const_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), it_ (it) {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const iterator1 &it):
                container_const_reference<self_type> (it ()), rank_ (it.rank_), i_ (it.i_), j_ (it.j_), it_ (it.it_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                if (rank_ == 1 && functor_type::fast1 ())
                    ++ it_;
                else
                    *this = (*this) ().find1 (rank_, index1 () + 1, j_, 1);
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                if (rank_ == 1 && functor_type::fast1 ())
                    -- it_;
                else
                    *this = (*this) ().find1 (rank_, index1 () - 1, j_, -1);
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*it_).second;
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 begin () const {
                const self_type &m = (*this) ();
                return m.find2 (1, index1 (), 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 end () const {
                const self_type &m = (*this) ();
                return m.find2 (1, index1 (), m.size2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    const self_type &m = (*this) ();
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*it_).first, m.size1 (), m.size2 ()) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*it_).first, m.size1 (), m.size2 ());
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    const self_type &m = (*this) ();
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*it_).first, m.size1 (), m.size2 ()) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*it_).first, m.size1 (), m.size2 ());
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator1 begin1 () const {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator1 end1 () const {
            return find1 (0, size1_, 0);
        }

        class iterator1:
            public container_reference<sparse_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator1, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename sparse_matrix::value_type value_type;
            typedef typename sparse_matrix::difference_type difference_type;
            typedef typename sparse_matrix::true_reference reference;
            typedef typename sparse_matrix::pointer pointer;
#endif
            typedef iterator2 dual_iterator_type;
            typedef reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator1 ():
                container_reference<self_type> (), rank_ (), i_ (), j_ (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator1 (self_type &m, int rank, size_type i, size_type j, const iterator_type &it):
                container_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator1 &operator ++ () {
                if (rank_ == 1 && functor_type::fast1 ())
                    ++ it_;
                else
                    *this = (*this) ().find1 (rank_, index1 () + 1, j_, 1);
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator1 &operator -- () {
                if (rank_ == 1 && functor_type::fast1 ())
                    -- it_;
                else
                    *this = (*this) ().find1 (rank_, index1 () - 1, j_, -1);
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*it_).second;
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 begin () const {
                self_type &m = (*this) ();
                return m.find2 (1, index1 (), 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 end () const {
                self_type &m = (*this) ();
                return m.find2 (1, index1 (), m.size2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    const self_type &m = (*this) ();
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*it_).first, m.size1 (), m.size2 ()) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*it_).first, m.size1 (), m.size2 ());
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    const self_type &m = (*this) ();
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*it_).first, m.size1 (), m.size2 ()) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*it_).first, m.size1 (), m.size2 ());
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator1 &operator = (const iterator1 &it) {
                container_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            iterator_type it_;

            friend class const_iterator1;
        };

        BOOST_UBLAS_INLINE
        iterator1 begin1 () {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        iterator1 end1 () {
            return find1 (0, size1_, 0);
        }

        class const_iterator2:
            public container_const_reference<sparse_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator2, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename sparse_matrix::value_type value_type;
            typedef typename sparse_matrix::difference_type difference_type;
            typedef typename sparse_matrix::const_reference reference;
            typedef const typename sparse_matrix::pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), rank_ (), i_ (), j_ (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &m, int rank, size_type i, size_type j, const const_iterator_type &it):
                container_const_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), it_ (it) {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const iterator2 &it):
                container_const_reference<self_type> (it ()), rank_ (it.rank_), i_ (it.i_), j_ (it.j_), it_ (it.it_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                if (rank_ == 1 && functor_type::fast2 ())
                    ++ it_;
                else
                    *this = (*this) ().find2 (rank_, i_, index2 () + 1, 1);
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                if (rank_ == 1 && functor_type::fast2 ())
                    -- it_;
                else
                    *this = (*this) ().find2 (rank_, i_, index2 () - 1, -1);
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*it_).second;
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 begin () const {
                const self_type &m = (*this) ();
                return m.find1 (1, 0, index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 end () const {
                const self_type &m = (*this) ();
                return m.find1 (1, m.size1 (), index2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    const self_type &m = (*this) ();
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*it_).first, m.size1 (), m.size2 ()) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*it_).first, m.size1 (), m.size2 ());
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    const self_type &m = (*this) ();
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*it_).first, m.size1 (), m.size2 ()) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*it_).first, m.size1 (), m.size2 ());
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator2 begin2 () const {
            return find2 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 end2 () const {
            return find2 (0, 0, size2_);
        }

        class iterator2:
            public container_reference<sparse_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator2, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename sparse_matrix::value_type value_type;
            typedef typename sparse_matrix::difference_type difference_type;
            typedef typename sparse_matrix::true_reference reference;
            typedef typename sparse_matrix::pointer pointer;
#endif
            typedef iterator1 dual_iterator_type;
            typedef reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator2 ():
                container_reference<self_type> (), rank_ (), i_ (), j_ (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator2 (self_type &m, int rank, size_type i, size_type j, const iterator_type &it):
                container_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator2 &operator ++ () {
                if (rank_ == 1 && functor_type::fast2 ())
                    ++ it_;
                else
                    *this = (*this) ().find2 (rank_, i_, index2 () + 1, 1);
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator2 &operator -- () {
                if (rank_ == 1 && functor_type::fast2 ())
                    -- it_;
                else
                    *this = (*this) ().find2 (rank_, i_, index2 () - 1, -1);
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*it_).second;
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 begin () const {
                self_type &m = (*this) ();
                return m.find1 (1, 0, index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 end () const {
                self_type &m = (*this) ();
                return m.find1 (1, m.size1 (), index2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    const self_type &m = (*this) ();
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*it_).first, m.size1 (), m.size2 ()) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*it_).first, m.size1 (), m.size2 ());
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    const self_type &m = (*this) ();
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*it_).first, m.size1 (), m.size2 ()) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*it_).first, m.size1 (), m.size2 ());
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator2 &operator = (const iterator2 &it) {
                container_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            iterator_type it_;

            friend class const_iterator2;
        };

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
        static const value_type zero_;
    };

    template<class T, class F, class A>
    const typename sparse_matrix<T, F, A>::value_type sparse_matrix<T, F, A>::zero_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME sparse_matrix<T, F, A>::value_type
#endif
        (0);

    // Array based sparse matrix class
    template<class T, class F, class A>
    class sparse_vector_of_sparse_vector:
        public matrix_expression<sparse_vector_of_sparse_vector<T, F, A> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<sparse_vector_of_sparse_vector<T, F, A> >::operator ();
#endif
        typedef typename A::size_type size_type;
        typedef typename A::difference_type difference_type;
        typedef T value_type;
        typedef const T &const_reference;
#ifndef BOOST_UBLAS_STRICT_MATRIX_SPARSE
        typedef typename detail::map_traits<typename A::data_value_type, T>::reference reference;
#else
        typedef sparse_matrix_element<sparse_vector_of_sparse_vector<T, F, A> > reference;
#endif
    private:
        typedef T &true_reference;
        typedef T *pointer;
        typedef A array_type;
        typedef const A const_array_type;
        typedef F functor_type;
        typedef sparse_vector_of_sparse_vector<T, F, A> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const matrix_const_reference<const self_type> const_closure_type;
#else
        typedef const matrix_reference<const self_type> const_closure_type;
#endif
        typedef matrix_reference<self_type> closure_type;
        typedef sparse_vector<T, typename A::value_type> vector_temporary_type;
        typedef self_type matrix_temporary_type;
        typedef typename A::value_type::second_type vector_data_value_type;
        typedef sparse_tag storage_category;
        typedef typename F::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector ():
            matrix_expression<self_type> (),
            size1_ (0), size2_ (0), non_zeros_ (0), data_ () {
            data_ [functor_type::size1 (size1_, size2_)] = vector_data_value_type ();
        }
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector (size_type size1, size_type size2, size_type non_zeros = 0):
            matrix_expression<self_type> (),
            size1_ (size1), size2_ (size2), non_zeros_ (non_zeros), data_ () {
            data_ [functor_type::size1 (size1_, size2_)] = vector_data_value_type ();
        }
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector (const sparse_vector_of_sparse_vector &m):
            matrix_expression<self_type> (),
            size1_ (m.size1_), size2_ (m.size2_), non_zeros_ (m.non_zeros_), data_ (m.data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector (const matrix_expression<AE> &ae, size_type non_zeros = 0):
            matrix_expression<self_type> (),
            size1_ (ae ().size1 ()), size2_ (ae ().size2 ()), non_zeros_ (non_zeros), data_ () {
            data_ [functor_type::size1 (size1_, size2_)] = vector_data_value_type ();
            matrix_assign (scalar_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
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
        size_type non_zeros () const {
            size_type non_zeros = 0;
            for (vector_const_iterator_type itv = data_ ().begin (); itv != data_ ().end (); ++ itv)
                non_zeros += detail::map_capacity (*itv);
            return non_zeros;
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
            // FIXME preserve unimplemented
            BOOST_UBLAS_CHECK (!preserve, internal_logic ());
            size1_ = size1;
            size2_ = size2;
            data ().clear ();
            data () [functor_type::size1 (size1_, size2_)] = vector_data_value_type ();
        }

        // Proxy support
#ifdef BOOST_UBLAS_STRICT_MATRIX_SPARSE
        BOOST_UBLAS_INLINE
        pointer find_element (size_type i, size_type j) {
            vector_iterator_type itv (data ().find (functor_type::element1 (i, size1_, j, size2_)));
            if (itv == data ().end () || (*itv).first != functor_type::element1 (i, size1_, j, size2_))
                return 0;
            iterator_type it ((*itv).second.find (functor_type::element2 (i, size1_, j, size2_)));
            if (it == (*itv).second.end () || (*it).first != functor_type::element2 (i, size1_, j, size2_))
                return 0;
            return &(*it).second;
        }
#endif

        // Element access
        BOOST_UBLAS_INLINE
        const_reference at_element (size_type i, size_type j) const {
            vector_const_iterator_type itv (data ().find (functor_type::element1 (i, size1_, j, size2_)));
            if (itv == data ().end () || (*itv).first != functor_type::element1 (i, size1_, j, size2_))
                return zero_;
            const_iterator_type it ((*itv).second.find (functor_type::element2 (i, size1_, j, size2_)));
            if (it == (*itv).second.end () || (*it).first != functor_type::element2 (i, size1_, j, size2_))
                return zero_;
            return (*it).second;
        }
        BOOST_UBLAS_INLINE
        true_reference at_element (size_type i, size_type j) {
            return data () [functor_type::element1 (i, size1_, j, size2_)] [functor_type::element2 (i, size1_, j, size2_)];
        }
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return at_element (i, j);
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
#ifndef BOOST_UBLAS_STRICT_MATRIX_SPARSE
            return at_element (i, j);
#else
            return reference (*this, i, j);
#endif
        }

        // Assignment
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector &operator = (const sparse_vector_of_sparse_vector &m) {
            if (this != &m) {
                size1_ = m.size1_;
                size2_ = m.size2_;
                non_zeros_ = m.non_zeros_;
                data () = m.data ();
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector &assign_temporary (sparse_vector_of_sparse_vector &m) {
            swap (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector &operator = (const matrix_expression<AE> &ae) {
            // return assign (self_type (ae, non_zeros_));
            self_type temporary (ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector& operator += (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this + ae, non_zeros_));
            self_type temporary (*this + ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector &plus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_plus_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector& operator -= (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this - ae, non_zeros_));
            self_type temporary (*this - ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector &minus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_minus_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<true_reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        sparse_vector_of_sparse_vector& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<true_reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (sparse_vector_of_sparse_vector &m) {
            if (this != &m) {
                std::swap (size1_, m.size1_);
                std::swap (size2_, m.size2_);
                std::swap (non_zeros_, m.non_zeros_);
                data ().swap (m.data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (sparse_vector_of_sparse_vector &m1, sparse_vector_of_sparse_vector &m2) {
            m1.swap (m2);
        }
#endif

        // Element insertion and erasure
        BOOST_UBLAS_INLINE
        void insert (size_type i, size_type j, const_reference t) {
            vector_iterator_type itv (data ().find (functor_type::element1 (i, size1_, j, size2_)));
            if (itv == data ().end ())
                itv = data ().insert (data ().end (), BOOST_UBLAS_TYPENAME array_type::value_type (functor_type::element1 (i, size1_, j, size2_), vector_data_value_type ()));
            BOOST_UBLAS_CHECK ((*itv).second.find (functor_type::element2 (i, size1_, j, size2_)) == (*itv).second.end (), bad_index ());
            (*itv).second.insert ((*itv).second.end (), BOOST_UBLAS_TYPENAME array_type::value_type::second_type::value_type (functor_type::element2 (i, size1_, j, size2_), t));
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i, size_type j) {
            vector_iterator_type itv (data ().find (functor_type::element1 (i, size1_, j, size2_)));
            if (itv == data ().end ())
                return;
            iterator_type it ((*itv).second.find (functor_type::element2 (i, size1_, j, size2_)));
            if (it == (*itv).second.end ())
                return;
            (*itv).second.erase (it);
        }
        BOOST_UBLAS_INLINE
        void clear () {
            data ().clear ();
            data_ [functor_type::size1 (size1_, size2_)] = vector_data_value_type ();
        }

        // Iterator types
    private:
        // Use storage iterators
        typedef typename A::const_iterator vector_const_iterator_type;
        typedef typename A::iterator vector_iterator_type;
        typedef typename A::value_type::second_type::const_iterator const_iterator_type;
        typedef typename A::value_type::second_type::iterator iterator_type;

    public:
        class const_iterator1;
        class iterator1;
        class const_iterator2;
        class iterator2;
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
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j, int direction = 1) const {
            BOOST_UBLAS_CHECK (data ().begin () != data ().end (), internal_logic ());
            for (;;) {
                vector_const_iterator_type itv (data ().lower_bound (functor_type::address1 (i, size1_, j, size2_)));
                vector_const_iterator_type itv_end (data ().end ());
                if (itv == itv_end)
                    return const_iterator1 (*this, rank, i, j, itv_end, (*(-- itv)).second.end ());

                const_iterator_type it ((*itv).second.lower_bound (functor_type::address2 (i, size1_, j, size2_)));
                const_iterator_type it_end ((*itv).second.end ());
                if (rank == 0)
                    return const_iterator1 (*this, rank, i, j, itv, it);
                if (it != it_end && (*it).first == functor_type::address2 (i, size1_, j, size2_))
                    return const_iterator1 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast1 ()) {
                        if (it == it_end)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        i = (*it).first;
                    } else {
                        if (i >= size1_)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        ++ i;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast1 ()) {
                        if (it == (*itv).second.begin ())
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        -- it;
                        i = (*it).first;
                    } else {
                        if (i == 0)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        -- i;
                    }
                }
            }
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j, int direction = 1) {
            BOOST_UBLAS_CHECK (data ().begin () != data ().end (), internal_logic ());
            for (;;) {
                vector_iterator_type itv (data ().lower_bound (functor_type::address1 (i, size1_, j, size2_)));
                vector_iterator_type itv_end (data ().end ());
                if (itv == itv_end)
                    return iterator1 (*this, rank, i, j, itv_end, (*(-- itv)).second.end ());

                iterator_type it ((*itv).second.lower_bound (functor_type::address2 (i, size1_, j, size2_)));
                iterator_type it_end ((*itv).second.end ());
                if (rank == 0)
                    return iterator1 (*this, rank, i, j, itv, it);
                if (it != it_end && (*it).first == functor_type::address2 (i, size1_, j, size2_))
                    return iterator1 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast1 ()) {
                        if (it == it_end)
                            return iterator1 (*this, rank, i, j, itv, it);
                        i = (*it).first;
                    } else {
                        if (i >= size1_)
                            return iterator1 (*this, rank, i, j, itv, it);
                        ++ i;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast1 ()) {
                        if (it == (*itv).second.begin ())
                            return iterator1 (*this, rank, i, j, itv, it);
                        -- it;
                        i = (*it).first;
                    } else {
                        if (i == 0)
                            return iterator1 (*this, rank, i, j, itv, it);
                        -- i;
                    }
                }
            }
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j, int direction = 1) const {
            BOOST_UBLAS_CHECK (data ().begin () != data ().end (), internal_logic ());
            for (;;) {
                vector_const_iterator_type itv (data ().lower_bound (functor_type::address1 (i, size1_, j, size2_)));
                vector_const_iterator_type itv_end (data ().end ());
                if (itv == itv_end)
                    return const_iterator2 (*this, rank, i, j, itv_end, (*(-- itv)).second.end ());

                const_iterator_type it ((*itv).second.lower_bound (functor_type::address2 (i, size1_, j, size2_)));
                const_iterator_type it_end ((*itv).second.end ());
                if (rank == 0)
                    return const_iterator2 (*this, rank, i, j, itv, it);
                if (it != it_end && (*it).first == functor_type::address2 (i, size1_, j, size2_))
                    return const_iterator2 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast2 ()) {
                        if (it == it_end)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        j = (*it).first;
                    } else {
                        if (j >= size2_)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        ++ j;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast2 ()) {
                        if (it == (*itv).second.begin ())
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        -- it;
                        j = (*it).first;
                    } else {
                        if (j == 0)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        -- j;
                    }
                }
            }
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j, int direction = 1) {
            BOOST_UBLAS_CHECK (data ().begin () != data ().end (), internal_logic ());
            for (;;) {
                vector_iterator_type itv (data ().lower_bound (functor_type::address1 (i, size1_, j, size2_)));
                vector_iterator_type itv_end (data ().end ());
                if (itv == itv_end)
                    return iterator2 (*this, rank, i, j, itv_end, (*(-- itv)).second.end ());

                iterator_type it ((*itv).second.lower_bound (functor_type::address2 (i, size1_, j, size2_)));
                iterator_type it_end ((*itv).second.end ());
                if (rank == 0)
                    return iterator2 (*this, rank, i, j, itv, it);
                if (it != it_end && (*it).first == functor_type::address2 (i, size1_, j, size2_))
                    return iterator2 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast2 ()) {
                        if (it == it_end)
                            return iterator2 (*this, rank, i, j, itv, it);
                        j = (*it).first;
                    } else {
                        if (j >= size2_)
                            return iterator2 (*this, rank, i, j, itv, it);
                        ++ j;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast2 ()) {
                        if (it == (*itv).second.begin ())
                            return iterator2 (*this, rank, i, j, itv, it);
                        -- it;
                        j = (*it).first;
                    } else {
                        if (j == 0)
                            return iterator2 (*this, rank, i, j, itv, it);
                        -- j;
                    }
                }
            }
        }


        class const_iterator1:
            public container_const_reference<sparse_vector_of_sparse_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator1, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename sparse_vector_of_sparse_vector::value_type value_type;
            typedef typename sparse_vector_of_sparse_vector::difference_type difference_type;
            typedef typename sparse_vector_of_sparse_vector::const_reference reference;
            typedef const typename sparse_vector_of_sparse_vector::pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &m, int rank, size_type i, size_type j, const vector_const_iterator_type &itv, const const_iterator_type &it):
                container_const_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const iterator1 &it):
                container_const_reference<self_type> (it ()), rank_ (it.rank_), i_ (it.i_), j_ (it.j_), itv_ (it.itv_), it_ (it.it_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                if (rank_ == 1 && functor_type::fast1 ())
                    ++ it_;
                else {
                    const self_type &m = (*this) ();
                    i_ = index1 () + 1;
                    if (rank_ == 1 && ++ itv_ == m.end1 ().itv_)
                        *this = m.find1 (rank_, i_, j_, 1);
                    else if (rank_ == 1) {
                        it_ = (*itv_).second.begin ();
                        if (it_ == (*itv_).second.end () || index2 () != j_)
                            *this = m.find1 (rank_, i_, j_, 1);
                    }
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                if (rank_ == 1 && functor_type::fast1 ())
                    -- it_;
                else {
                    const self_type &m = (*this) ();
                    i_ = index1 () - 1;
                    if (rank_ == 1 && -- itv_ == m.end1 ().itv_)
                        *this = m.find1 (rank_, i_, j_, -1);
                    else if (rank_ == 1) {
                        it_ = (*itv_).second.begin ();
                        if (it_ == (*itv_).second.end () || index2 () != j_)
                            *this = m.find1 (rank_, i_, j_, -1);
                    }
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*it_).second;
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 begin () const {
                const self_type &m = (*this) ();
                return m.find2 (1, index1 (), 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 end () const {
                const self_type &m = (*this) ();
                return m.find2 (1, index1 (), m.size2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*itv_).first, (*it_).first) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*itv_).first, (*it_).first);
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*itv_).first, (*it_).first) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*itv_).first, (*it_).first);
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_const_iterator_type itv_;
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator1 begin1 () const {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator1 end1 () const {
            return find1 (0, size1_, 0);
        }

        class iterator1:
            public container_reference<sparse_vector_of_sparse_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator1, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename sparse_vector_of_sparse_vector::value_type value_type;
            typedef typename sparse_vector_of_sparse_vector::difference_type difference_type;
            typedef typename sparse_vector_of_sparse_vector::true_reference reference;
            typedef typename sparse_vector_of_sparse_vector::pointer pointer;
#endif
            typedef iterator2 dual_iterator_type;
            typedef reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator1 ():
                container_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator1 (self_type &m, int rank, size_type i, size_type j, const vector_iterator_type &itv, const iterator_type &it):
                container_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator1 &operator ++ () {
                if (rank_ == 1 && functor_type::fast1 ())
                    ++ it_;
                else {
                    self_type &m = (*this) ();
                    i_ = index1 () + 1;
                    if (rank_ == 1 && ++ itv_ == m.end1 ().itv_)
                        *this = m.find1 (rank_, i_, j_, 1);
                    else if (rank_ == 1) {
                        it_ = (*itv_).second.begin ();
                        if (it_ == (*itv_).second.end () || index2 () != j_)
                            *this = m.find1 (rank_, i_, j_, 1);
                    }
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator1 &operator -- () {
                if (rank_ == 1 && functor_type::fast1 ())
                    -- it_;
                else {
                    self_type &m = (*this) ();
                    i_ = index1 () - 1;
                    if (rank_ == 1 && -- itv_ == m.end1 ().itv_)
                        *this = m.find1 (rank_, i_, j_, -1);
                    else if (rank_ == 1) {
                        it_ = (*itv_).second.begin ();
                        if (it_ == (*itv_).second.end () || index2 () != j_)
                            *this = m.find1 (rank_, i_, j_, -1);
                    }
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*it_).second;
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 begin () const {
                self_type &m = (*this) ();
                return m.find2 (1, index1 (), 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 end () const {
                self_type &m = (*this) ();
                return m.find2 (1, index1 (), m.size2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*itv_).first, (*it_).first) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*itv_).first, (*it_).first);
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*itv_).first, (*it_).first) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*itv_).first, (*it_).first);
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator1 &operator = (const iterator1 &it) {
                container_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_iterator_type itv_;
            iterator_type it_;

            friend class const_iterator1;
        };

        BOOST_UBLAS_INLINE
        iterator1 begin1 () {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        iterator1 end1 () {
            return find1 (0, size1_, 0);
        }

        class const_iterator2:
            public container_const_reference<sparse_vector_of_sparse_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator2, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename sparse_vector_of_sparse_vector::value_type value_type;
            typedef typename sparse_vector_of_sparse_vector::difference_type difference_type;
            typedef typename sparse_vector_of_sparse_vector::const_reference reference;
            typedef const typename sparse_vector_of_sparse_vector::pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &m, int rank, size_type i, size_type j, const vector_const_iterator_type &itv, const const_iterator_type &it):
                container_const_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const iterator2 &it):
                container_const_reference<self_type> (it ()), rank_ (it.rank_), i_ (it.i_), j_ (it.j_), itv_ (it.itv_), it_ (it.it_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                if (rank_ == 1 && functor_type::fast2 ())
                    ++ it_;
                else {
                    const self_type &m = (*this) ();
                    j_ = index2 () + 1;
                    if (rank_ == 1 && ++ itv_ == m.end2 ().itv_)
                        *this = m.find2 (rank_, i_, j_, 1);
                    else if (rank_ == 1) {
                        it_ = (*itv_).second.begin ();
                        if (it_ == (*itv_).second.end () || index1 () != i_)
                            *this = m.find2 (rank_, i_, j_, 1);
                    }
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                if (rank_ == 1 && functor_type::fast2 ())
                    -- it_;
                else {
                    const self_type &m = (*this) ();
                    j_ = index2 () - 1;
                    if (rank_ == 1 && -- itv_ == m.end2 ().itv_)
                        *this = m.find2 (rank_, i_, j_, -1);
                    else if (rank_ == 1) {
                        it_ = (*itv_).second.begin ();
                        if (it_ == (*itv_).second.end () || index1 () != i_)
                            *this = m.find2 (rank_, i_, j_, -1);
                    }
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*it_).second;
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 begin () const {
                const self_type &m = (*this) ();
                return m.find1 (1, 0, index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 end () const {
                const self_type &m = (*this) ();
                return m.find1 (1, m.size1 (), index2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*itv_).first, (*it_).first) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*itv_).first, (*it_).first);
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*itv_).first, (*it_).first) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*itv_).first, (*it_).first);
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_const_iterator_type itv_;
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator2 begin2 () const {
            return find2 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 end2 () const {
            return find2 (0, 0, size2_);
        }

        class iterator2:
            public container_reference<sparse_vector_of_sparse_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator2, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename sparse_vector_of_sparse_vector::value_type value_type;
            typedef typename sparse_vector_of_sparse_vector::difference_type difference_type;
            typedef typename sparse_vector_of_sparse_vector::true_reference reference;
            typedef typename sparse_vector_of_sparse_vector::pointer pointer;
#endif
            typedef iterator1 dual_iterator_type;
            typedef reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator2 ():
                container_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator2 (self_type &m, int rank, size_type i, size_type j, const vector_iterator_type &itv, const iterator_type &it):
                container_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator2 &operator ++ () {
                if (rank_ == 1 && functor_type::fast2 ())
                    ++ it_;
                else {
                    self_type &m = (*this) ();
                    j_ = index2 () + 1;
                    if (rank_ == 1 && ++ itv_ == m.end2 ().itv_)
                        *this = m.find2 (rank_, i_, j_, 1);
                    else if (rank_ == 1) {
                        it_ = (*itv_).second.begin ();
                        if (it_ == (*itv_).second.end () || index1 () != i_)
                            *this = m.find2 (rank_, i_, j_, 1);
                    }
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator2 &operator -- () {
                if (rank_ == 1 && functor_type::fast2 ())
                    -- it_;
                else {
                    self_type &m = (*this) ();
                    j_ = index2 () - 1;
                    if (rank_ == 1 && -- itv_ == m.end2 ().itv_)
                        *this = m.find2 (rank_, i_, j_, -1);
                    else if (rank_ == 1) {
                        it_ = (*itv_).second.begin ();
                        if (it_ == (*itv_).second.end () || index1 () != i_)
                            *this = m.find2 (rank_, i_, j_, -1);
                    }
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*it_).second;
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 begin () const {
                self_type &m = (*this) ();
                return m.find1 (1, 0, index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 end () const {
                self_type &m = (*this) ();
                return m.find1 (1, m.size1 (), index2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*itv_).first, (*it_).first) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*itv_).first, (*it_).first);
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*itv_).first, (*it_).first) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*itv_).first, (*it_).first);
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator2 &operator = (const iterator2 &it) {
                container_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_iterator_type itv_;
            iterator_type it_;

            friend class const_iterator2;
        };

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
        size_type non_zeros_;
        array_type data_;
        static const value_type zero_;
    };

    template<class T, class F, class A>
    const typename sparse_vector_of_sparse_vector<T, F, A>::value_type sparse_vector_of_sparse_vector<T, F, A>::zero_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME sparse_vector_of_sparse_vector<T, F, A>::value_type
#endif
        (0);


    // Array based sparse matrix class
    // Thanks to Kresimir Fresl for extending this to cover different index bases.
    template<class T, class F, std::size_t IB, class IA, class TA>
    class compressed_matrix:
        public matrix_expression<compressed_matrix<T, F, IB, IA, TA> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<compressed_matrix<T, F, IB, IA, TA> >::operator ();
#endif
        // ISSUE require type consistency check for IA TA and IA::value_type
        typedef typename IA::size_type size_type;
        typedef typename IA::difference_type difference_type;
        typedef T value_type;
        typedef const T &const_reference;
#ifndef BOOST_UBLAS_STRICT_MATRIX_SPARSE
        typedef T &reference;
#else
        typedef sparse_matrix_element<compressed_matrix<T, F, IB, IA, TA> > reference;
#endif
        typedef IA index_array_type;
        typedef TA value_array_type;
    private:
        typedef T &true_reference;
        typedef T *pointer;
        typedef F functor_type;
        typedef compressed_matrix<T, F, IB, IA, TA> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const matrix_const_reference<const self_type> const_closure_type;
#else
        typedef const matrix_reference<const self_type> const_closure_type;
#endif
        typedef matrix_reference<self_type> closure_type;
        typedef compressed_vector<T, IB, IA, TA> vector_temporary_type;
        typedef self_type matrix_temporary_type;
        typedef sparse_tag storage_category;
        typedef typename F::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        compressed_matrix ():
            matrix_expression<self_type> (),
            size1_ (0), size2_ (0), non_zeros_ (restrict_nz (0)),
            filled1_ (1), filled2_ (0),
            index1_data_ (functor_type::size1 (size1_, size2_) + 1),
            index2_data_ (non_zeros_), value_data_ (non_zeros_) {
            index1_data_ [filled1_ - 1] = k_based (filled2_);
        }
        BOOST_UBLAS_INLINE
        compressed_matrix (size_type size1, size_type size2, size_type non_zeros = 0):
            matrix_expression<self_type> (),
            size1_ (size1), size2_ (size2), non_zeros_ (restrict_nz (non_zeros)),
            filled1_ (1), filled2_ (0),
            index1_data_ (functor_type::size1 (size1_, size2_) + 1),
            index2_data_ (non_zeros_), value_data_ (non_zeros_) {
            index1_data_ [filled1_ - 1] = k_based (filled2_);
        }
        BOOST_UBLAS_INLINE
        compressed_matrix (const compressed_matrix &m):
            matrix_expression<self_type> (),
            size1_ (m.size1_), size2_ (m.size2_), non_zeros_ (m.non_zeros_),
            filled1_ (m.filled1_), filled2_ (m.filled2_),
            index1_data_ (m.index1_data_),
            index2_data_ (m.index2_data_), value_data_ (m.value_data_) {
            BOOST_UBLAS_CHECK (index1_data_ [filled1_ - 1] == k_based (filled2_), internal_logic ());
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_matrix (const matrix_expression<AE> &ae, size_type non_zeros = 0):
            matrix_expression<self_type> (),
            size1_ (ae ().size1 ()), size2_ (ae ().size2 ()), non_zeros_ (restrict_nz (non_zeros)),
            filled1_ (1), filled2_ (0),
            index1_data_ (functor_type::size1 (ae ().size1 (), ae ().size2 ()) + 1),
            index2_data_ (non_zeros_), value_data_ (non_zeros_) {
            index1_data_ [filled1_ - 1] = k_based (filled2_);
            matrix_assign (scalar_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
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
        size_type non_zeros () const {
            return non_zeros_;
        }
        BOOST_UBLAS_INLINE
        size_type filled () const {
            return filled2_;
        }
        BOOST_UBLAS_INLINE
        size_type &filled1 () {
            return filled1_;
        }
        BOOST_UBLAS_INLINE
        size_type &filled2 () {
            return filled2_;
        }
        BOOST_UBLAS_INLINE
        static size_type index_base () {
            return IB;
        }
        BOOST_UBLAS_INLINE
        const index_array_type &index1_data () const {
            return index1_data_;
        }
        BOOST_UBLAS_INLINE
        index_array_type &index1_data () {
            return index1_data_;
        }
        BOOST_UBLAS_INLINE
        const index_array_type &index2_data () const {
            return index2_data_;
        }
        BOOST_UBLAS_INLINE
        index_array_type &index2_data () {
            return index2_data_;
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
    private:
        BOOST_UBLAS_INLINE
        size_type restrict_nz (size_type non_zeros) const {
            non_zeros = (std::max) (non_zeros, (std::min) (size1_, size2_));
            // Guarding against overflow.
            // Thanks to Alexei Novakov for the hint.
            // non_zeros_ = (std::min) (non_zeros, size1_ * size2_);
            if (size1_ > 0 && non_zeros / size1_ >= size2_)
                non_zeros = size1_ * size2_;
            return non_zeros;
        }
    public:
        BOOST_UBLAS_INLINE
        void resize (size_type size1, size_type size2, bool preserve = true) {
            // FIXME preserve unimplemented
            BOOST_UBLAS_CHECK (!preserve, internal_logic ());
            size1_ = size1;
            size2_ = size2;
            non_zeros_ = restrict_nz (non_zeros_);
            filled1_ = 1;
            filled2_ = 0;
            index1_data ().resize (functor_type::size1 (size1_, size2_) + 1);
            index2_data ().resize (non_zeros_);
            value_data ().resize (non_zeros_);
            index1_data_ [filled1_ - 1] = k_based (filled2_);
        }

        // Reserving
        BOOST_UBLAS_INLINE
        void reserve (size_type non_zeros, bool preserve = true) {
            non_zeros_ = restrict_nz (non_zeros);
            if (preserve) {
                index2_data ().resize (non_zeros_, size_type ());
                value_data ().resize (non_zeros_, value_type ());
                filled1_ = (std::min) (non_zeros_ + 1, filled1_);
                filled2_ = (std::min) (non_zeros_, filled2_);
            }
            else {
                index2_data ().resize (non_zeros_);
                value_data ().resize (non_zeros_);
                filled1_ = 1;
                filled2_ = 0;
            }
            BOOST_UBLAS_CHECK (index1_data () [filled1_ - 1] == k_based (filled2_), internal_logic ());
       }

        // Proxy support
#ifdef BOOST_UBLAS_STRICT_MATRIX_SPARSE
        BOOST_UBLAS_INLINE
        pointer find_element (size_type i, size_type j) {
            size_type element1 (functor_type::element1 (i, size1_, j, size2_));
            size_type element2 (functor_type::element2 (i, size1_, j, size2_));
            if (filled1_ <= element1 + 1)
                return 0;
            vector_const_iterator_type itv (index1_data ().begin () + element1);
            const_iterator_type it_begin (index2_data ().begin () + zero_based (*itv));
            const_iterator_type it_end (index2_data ().begin () + zero_based (*(itv + 1)));
            const_iterator_type it (detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ()));
            if (it == it_end || *it != k_based (element2))
                return 0;
            return &value_data () [it - index2_data ().begin ()];
        }
#endif

        // Element access
        BOOST_UBLAS_INLINE
        const_reference at_element (size_type i, size_type j) const {
            size_type element1 (functor_type::element1 (i, size1_, j, size2_));
            size_type element2 (functor_type::element2 (i, size1_, j, size2_));
            if (filled1_ <= element1 + 1)
                return zero_;
            vector_const_iterator_type itv (index1_data ().begin () + element1);
            const_iterator_type it_begin (index2_data ().begin () + zero_based (*itv));
            const_iterator_type it_end (index2_data ().begin () + zero_based (*(itv + 1)));
            const_iterator_type it (detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ()));
            if (it == it_end || *it != k_based (element2))
                return zero_;
            return value_data () [it - index2_data ().begin ()];
        }
        BOOST_UBLAS_INLINE
        true_reference at_element (size_type i, size_type j) {
            size_type element1 (functor_type::element1 (i, size1_, j, size2_));
            size_type element2 (functor_type::element2 (i, size1_, j, size2_));
            if (filled1_ <= element1 + 1)
                insert (i, j, value_type (0));
            vector_const_iterator_type itv (index1_data ().begin () + element1);
            const_iterator_type it_begin (index2_data ().begin () + zero_based (*itv));
            const_iterator_type it_end (index2_data ().begin () + zero_based (*(itv + 1)));
            const_iterator_type it (detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ()));
            if (it == it_end || *it != k_based (element2)) {
                insert (i, j, value_type (0));
                it_begin = index2_data ().begin () + zero_based (*itv);
                it_end = index2_data ().begin () + zero_based (*(itv + 1));
                it = detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ());
            }
            return value_data () [it - index2_data ().begin ()];
        }
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return at_element (i, j);
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
#ifndef BOOST_UBLAS_STRICT_MATRIX_SPARSE
            return at_element (i, j);
#else
            return reference (*this, i, j);
#endif
        }

        // Assignment
        BOOST_UBLAS_INLINE
        compressed_matrix &operator = (const compressed_matrix &m) {
            if (this != &m) {
                size1_ = m.size1_;
                size2_ = m.size2_;
                non_zeros_ = m.non_zeros_;
                filled1_ = m.filled1_;
                filled2_ = m.filled2_;
                index1_data () = m.index1_data ();
                index2_data () = m.index2_data ();
                value_data () = m.value_data ();
                BOOST_UBLAS_CHECK (functor_type::size1 (size1_, size2_) + 1 == index1_data ().size (), internal_logic ());
                BOOST_UBLAS_CHECK (non_zeros_ == index2_data ().size (), internal_logic ());
                BOOST_UBLAS_CHECK (non_zeros_ == value_data ().size (), internal_logic ());
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        compressed_matrix &assign_temporary (compressed_matrix &m) {
            swap (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_matrix &operator = (const matrix_expression<AE> &ae) {
            // return assign (self_type (ae, non_zeros_));
            self_type temporary (ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_matrix &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_matrix& operator += (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this + ae, non_zeros_));
            self_type temporary (*this + ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_matrix &plus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_plus_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_matrix& operator -= (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this - ae, non_zeros_));
            self_type temporary (*this - ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        compressed_matrix &minus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_minus_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        compressed_matrix& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<true_reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        compressed_matrix& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<true_reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (compressed_matrix &m) {
            if (this != &m) {
                std::swap (size1_, m.size1_);
                std::swap (size2_, m.size2_);
                std::swap (non_zeros_, m.non_zeros_);
                std::swap (filled1_, m.filled1_);
                std::swap (filled2_, m.filled2_);
                index1_data ().swap (m.index1_data ());
                index2_data ().swap (m.index2_data ());
                value_data ().swap (m.value_data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (compressed_matrix &m1, compressed_matrix &m2) {
            m1.swap (m2);
        }
#endif

        // Element insertion and erasure
        BOOST_UBLAS_INLINE
        void push_back (size_type i, size_type j, const_reference t) {
            BOOST_UBLAS_CHECK (index1_data () [filled1_ - 1] == k_based (filled2_), internal_logic ());
            if (filled2_ >= non_zeros_)
                reserve (2 * non_zeros_, true);
            size_type element1 = functor_type::element1 (i, size1_, j, size2_);
            size_type element2 = functor_type::element2 (i, size1_, j, size2_);
            while (filled1_ < element1 + 2) {
                index1_data () [filled1_] = k_based (filled2_);
                ++ filled1_;
            }
            if (filled1_ == element1 + 2 &&
                (filled2_ == zero_based (index1_data () [filled1_ - 2]) ||
                 index2_data () [filled2_ - 1] < k_based (element2))) {
                ++ filled2_;
                index1_data () [filled1_ - 1] = k_based (filled2_);
                index2_data () [filled2_ - 1] = k_based (element2);
                value_data () [filled2_ - 1] = t;
                BOOST_UBLAS_CHECK (index1_data () [filled1_ - 1] == k_based (filled2_), internal_logic ());
                return;
            }
            external_logic ().raise ();
        }
        BOOST_UBLAS_INLINE
        void insert (size_type i, size_type j, const_reference t) {
            BOOST_UBLAS_CHECK (index1_data () [filled1_ - 1] == k_based (filled2_), internal_logic ());
            if (filled2_ >= non_zeros_)
                reserve (2 * non_zeros_, true);
            size_type element1 = functor_type::element1 (i, size1_, j, size2_);
            size_type element2 = functor_type::element2 (i, size1_, j, size2_);
            while (filled1_ < element1 + 2) {
                index1_data () [filled1_] = k_based (filled2_);
                ++ filled1_;
            }
            vector_iterator_type itv (index1_data ().begin () + element1);
            iterator_type it_begin (index2_data ().begin () + zero_based (*itv));
            iterator_type it_end (index2_data ().begin () + zero_based (*(itv + 1)));
            iterator_type it (detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ()));
            difference_type n = it - index2_data ().begin ();
            BOOST_UBLAS_CHECK (it == it_end || *it != k_based (element2), external_logic ());
            ++ filled2_;
            it = index2_data ().begin () + n;
            std::copy_backward (it, index2_data ().begin () + filled2_ - 1, index2_data ().begin () + filled2_);
            *it = k_based (element2);
            typename value_array_type::iterator itt (value_data ().begin () + n);
            std::copy_backward (itt, value_data ().begin () + filled2_ - 1, value_data ().begin () + filled2_);
            *itt = t;
            while (element1 + 1 < filled1_) {
                ++ index1_data () [element1 + 1];
                ++ element1;
            }
            BOOST_UBLAS_CHECK (index1_data () [filled1_ - 1] == k_based (filled2_), internal_logic ());
        }
        BOOST_UBLAS_INLINE
        void pop_back () {
            BOOST_UBLAS_CHECK (filled1_ > 0 && filled2_ > 0, external_logic ());
            BOOST_UBLAS_CHECK (index1_data () [filled1_ - 1] == k_based (filled2_), internal_logic ());
            -- filled2_;
            while (index1_data () [filled1_ - 2] > k_based (filled2_)) {
                index1_data () [filled1_ - 1] = 0;
                -- filled1_;
            }
            -- index1_data () [filled1_ - 1];
            BOOST_UBLAS_CHECK (index1_data () [filled1_ - 1] == k_based (filled2_), internal_logic ());
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (index1_data () [filled1_ - 1] == k_based (filled2_), internal_logic ());
            size_type element1 = functor_type::element1 (i, size1_, j, size2_);
            size_type element2 = functor_type::element2 (i, size1_, j, size2_);
            if (element1 + 1 > filled1_)
                return;
            vector_iterator_type itv (index1_data ().begin () + element1);
            iterator_type it_begin (index2_data ().begin () + zero_based (*itv));
            iterator_type it_end (index2_data ().begin () + zero_based (*(itv + 1)));
            iterator_type it (detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ()));
            if (it != it_end && *it == k_based (element2)) {
                difference_type n = it - index2_data ().begin ();
                std::copy (it + 1, index2_data ().begin () + filled2_, it);
                typename value_array_type::iterator itt (value_data ().begin () + n);
                std::copy (itt + 1, value_data ().begin () + filled2_, itt);
                -- filled2_;
                while (index1_data () [filled1_ - 2] > k_based (filled2_)) {
                    index1_data () [filled1_ - 1] = 0;
                    -- filled1_;
                }
                while (element1 + 1 < filled1_) {
                    -- index1_data () [element1 + 1];
                    ++ element1;
                }
            }
            BOOST_UBLAS_CHECK (index1_data () [filled1_ - 1] == k_based (filled2_), internal_logic ());
        }
        BOOST_UBLAS_INLINE
        void clear () {
            filled1_ = 1;
            filled2_ = 0;
            index1_data_ [filled1_ - 1] = k_based (filled2_);
        }

        // Iterator types
    private:
        // Use index array iterator
        typedef typename IA::const_iterator vector_const_iterator_type;
        typedef typename IA::iterator vector_iterator_type;
        typedef typename IA::const_iterator const_iterator_type;
        typedef typename IA::iterator iterator_type;

    public:
        class const_iterator1;
        class iterator1;
        class const_iterator2;
        class iterator2;
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
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j, int direction = 1) const {
            for (;;) {
                size_type address1 (functor_type::address1 (i, size1_, j, size2_));
                size_type address2 (functor_type::address2 (i, size1_, j, size2_));
                vector_const_iterator_type itv (index1_data ().begin () + (std::min) (filled1_ - 1, address1));
                if (filled1_ <= address1 + 1)
                    return const_iterator1 (*this, rank, i, j, itv, index2_data ().begin () + filled2_);

                const_iterator_type it_begin (index2_data ().begin () + zero_based (*itv));
                const_iterator_type it_end (index2_data ().begin () + zero_based (*(itv + 1)));

                const_iterator_type it (detail::lower_bound (it_begin, it_end, k_based (address2), std::less<size_type> ()));
                if (rank == 0)
                    return const_iterator1 (*this, rank, i, j, itv, it);
                if (it != it_end && zero_based (*it) == address2)
                    return const_iterator1 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast1 ()) {
                        if (it == it_end)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        i = zero_based (*it);
                    } else {
                        if (i >= size1_)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        ++ i;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast1 ()) {
                        if (it == index2_data ().begin () + zero_based (*itv))
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        i = zero_based (*(it - 1));
                    } else {
                        if (i == 0)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        -- i;
                    }
                }
            }
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j, int direction = 1) {
            for (;;) {
                size_type address1 (functor_type::address1 (i, size1_, j, size2_));
                size_type address2 (functor_type::address2 (i, size1_, j, size2_));
                vector_iterator_type itv (index1_data ().begin () + (std::min) (filled1_ - 1, address1));
                if (filled1_ <= address1 + 1)
                    return iterator1 (*this, rank, i, j, itv, index2_data ().begin () + filled2_);

                iterator_type it_begin (index2_data ().begin () + zero_based (*itv));
                iterator_type it_end (index2_data ().begin () + zero_based (*(itv + 1)));

                iterator_type it (detail::lower_bound (it_begin, it_end, k_based (address2), std::less<size_type> ()));
                if (rank == 0)
                    return iterator1 (*this, rank, i, j, itv, it);
                if (it != it_end && zero_based (*it) == address2)
                    return iterator1 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast1 ()) {
                        if (it == it_end)
                            return iterator1 (*this, rank, i, j, itv, it);
                        i = zero_based (*it);
                    } else {
                        if (i >= size1_)
                            return iterator1 (*this, rank, i, j, itv, it);
                        ++ i;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast1 ()) {
                        if (it == index2_data ().begin () + zero_based (*itv))
                            return iterator1 (*this, rank, i, j, itv, it);
                        i = zero_based (*(it - 1));
                    } else {
                        if (i == 0)
                            return iterator1 (*this, rank, i, j, itv, it);
                        -- i;
                    }
                }
            }
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j, int direction = 1) const {
            for (;;) {
                size_type address1 (functor_type::address1 (i, size1_, j, size2_));
                size_type address2 (functor_type::address2 (i, size1_, j, size2_));
                vector_const_iterator_type itv (index1_data ().begin () + (std::min) (filled1_ - 1, address1));
                if (filled1_ <= address1 + 1)
                    return const_iterator2 (*this, rank, i, j, itv, index2_data ().begin () + filled2_);

                const_iterator_type it_begin (index2_data ().begin () + zero_based (*itv));
                const_iterator_type it_end (index2_data ().begin () + zero_based (*(itv + 1)));

                const_iterator_type it (detail::lower_bound (it_begin, it_end, k_based (address2), std::less<size_type> ()));
                if (rank == 0)
                    return const_iterator2 (*this, rank, i, j, itv, it);
                if (it != it_end && zero_based (*it) == address2)
                    return const_iterator2 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast2 ()) {
                        if (it == it_end)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        j = zero_based (*it);
                    } else {
                        if (j >= size2_)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        ++ j;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast2 ()) {
                        if (it == index2_data ().begin () + zero_based (*itv))
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        j = zero_based (*(it - 1));
                    } else {
                        if (j == 0)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        -- j;
                    }
                }
            }
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j, int direction = 1) {
            for (;;) {
                size_type address1 (functor_type::address1 (i, size1_, j, size2_));
                size_type address2 (functor_type::address2 (i, size1_, j, size2_));
                vector_iterator_type itv (index1_data ().begin () + (std::min) (filled1_ - 1, address1));
                if (filled1_ <= address1 + 1)
                    return iterator2 (*this, rank, i, j, itv, index2_data ().begin () + filled2_);

                iterator_type it_begin (index2_data ().begin () + zero_based (*itv));
                iterator_type it_end (index2_data ().begin () + zero_based (*(itv + 1)));

                iterator_type it (detail::lower_bound (it_begin, it_end, k_based (address2), std::less<size_type> ()));
                if (rank == 0)
                    return iterator2 (*this, rank, i, j, itv, it);
                if (it != it_end && zero_based (*it) == address2)
                    return iterator2 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast2 ()) {
                        if (it == it_end)
                            return iterator2 (*this, rank, i, j, itv, it);
                        j = zero_based (*it);
                    } else {
                        if (j >= size2_)
                            return iterator2 (*this, rank, i, j, itv, it);
                        ++ j;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast2 ()) {
                        if (it == index2_data ().begin () + zero_based (*itv))
                            return iterator2 (*this, rank, i, j, itv, it);
                        j = zero_based (*(it - 1));
                    } else {
                        if (j == 0)
                            return iterator2 (*this, rank, i, j, itv, it);
                        -- j;
                    }
                }
            }
        }


        class const_iterator1:
            public container_const_reference<compressed_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator1, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename compressed_matrix::value_type value_type;
            typedef typename compressed_matrix::difference_type difference_type;
            typedef typename compressed_matrix::const_reference reference;
            typedef const typename compressed_matrix::pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &m, int rank, size_type i, size_type j, const vector_const_iterator_type &itv, const const_iterator_type &it):
                container_const_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const iterator1 &it):
                container_const_reference<self_type> (it ()), rank_ (it.rank_), i_ (it.i_), j_ (it.j_), itv_ (it.itv_), it_ (it.it_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                if (rank_ == 1 && functor_type::fast1 ())
                    ++ it_;
                else {
                    i_ = index1 () + 1;
                    if (rank_ == 1)
                        *this = (*this) ().find1 (rank_, i_, j_, 1);
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                if (rank_ == 1 && functor_type::fast1 ())
                    -- it_;
                else {
                    i_ = index1 () - 1;
                    if (rank_ == 1)
                        *this = (*this) ().find1 (rank_, i_, j_, -1);
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*this) ().value_data () [it_ - (*this) ().index2_data ().begin ()];
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 begin () const {
                const self_type &m = (*this) ();
                return m.find2 (1, index1 (), 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 end () const {
                const self_type &m = (*this) ();
                return m.find2 (1, index1 (), m.size2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_)) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_));
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_)) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_));
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_const_iterator_type itv_;
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator1 begin1 () const {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator1 end1 () const {
            return find1 (0, size1_, 0);
        }

        class iterator1:
            public container_reference<compressed_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator1, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename compressed_matrix::value_type value_type;
            typedef typename compressed_matrix::difference_type difference_type;
            typedef typename compressed_matrix::true_reference reference;
            typedef typename compressed_matrix::pointer pointer;
#endif
            typedef iterator2 dual_iterator_type;
            typedef reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator1 ():
                container_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator1 (self_type &m, int rank, size_type i, size_type j, const vector_iterator_type &itv, const iterator_type &it):
                container_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator1 &operator ++ () {
                if (rank_ == 1 && functor_type::fast1 ())
                    ++ it_;
                else {
                    i_ = index1 () + 1;
                    if (rank_ == 1)
                        *this = (*this) ().find1 (rank_, i_, j_, 1);
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator1 &operator -- () {
                if (rank_ == 1 && functor_type::fast1 ())
                    -- it_;
                else {
                    i_ = index1 () - 1;
                    if (rank_ == 1)
                        *this = (*this) ().find1 (rank_, i_, j_, -1);
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*this) ().value_data () [it_ - (*this) ().index2_data ().begin ()];
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 begin () const {
                self_type &m = (*this) ();
                return m.find2 (1, index1 (), 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 end () const {
                self_type &m = (*this) ();
                return m.find2 (1, index1 (), m.size2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_)) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_));
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_)) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_));
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator1 &operator = (const iterator1 &it) {
                container_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_iterator_type itv_;
            iterator_type it_;

            friend class const_iterator1;
        };

        BOOST_UBLAS_INLINE
        iterator1 begin1 () {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        iterator1 end1 () {
            return find1 (0, size1_, 0);
        }

        class const_iterator2:
            public container_const_reference<compressed_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator2, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename compressed_matrix::value_type value_type;
            typedef typename compressed_matrix::difference_type difference_type;
            typedef typename compressed_matrix::const_reference reference;
            typedef const typename compressed_matrix::pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &m, int rank, size_type i, size_type j, const vector_const_iterator_type itv, const const_iterator_type &it):
                container_const_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const iterator2 &it):
                container_const_reference<self_type> (it ()), rank_ (it.rank_), i_ (it.i_), j_ (it.j_), itv_ (it.itv_), it_ (it.it_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                if (rank_ == 1 && functor_type::fast2 ())
                    ++ it_;
                else {
                    j_ = index2 () + 1;
                    if (rank_ == 1)
                        *this = (*this) ().find2 (rank_, i_, j_, 1);
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                if (rank_ == 1 && functor_type::fast2 ())
                    -- it_;
                else {
                    j_ = index2 () - 1;
                    if (rank_ == 1)
                        *this = (*this) ().find2 (rank_, i_, j_, -1);
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*this) ().value_data () [it_ - (*this) ().index2_data ().begin ()];
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 begin () const {
                const self_type &m = (*this) ();
                return m.find1 (1, 0, index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 end () const {
                const self_type &m = (*this) ();
                return m.find1 (1, m.size1 (), index2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_)) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_));
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_)) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_));
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_const_iterator_type itv_;
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator2 begin2 () const {
            return find2 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 end2 () const {
            return find2 (0, 0, size2_);
        }

        class iterator2:
            public container_reference<compressed_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator2, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename compressed_matrix::value_type value_type;
            typedef typename compressed_matrix::difference_type difference_type;
            typedef typename compressed_matrix::true_reference reference;
            typedef typename compressed_matrix::pointer pointer;
#endif
            typedef iterator1 dual_iterator_type;
            typedef reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator2 ():
                container_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator2 (self_type &m, int rank, size_type i, size_type j, const vector_iterator_type &itv, const iterator_type &it):
                container_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator2 &operator ++ () {
                if (rank_ == 1 && functor_type::fast2 ())
                    ++ it_;
                else {
                    j_ = index2 () + 1;
                    if (rank_ == 1)
                        *this = (*this) ().find2 (rank_, i_, j_, 1);
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator2 &operator -- () {
                if (rank_ == 1 && functor_type::fast2 ())
                    -- it_;
                else {
                    j_ = index2 ();
                    if (rank_ == 1)
                        *this = (*this) ().find2 (rank_, i_, j_, -1);
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*this) ().value_data () [it_ - (*this) ().index2_data ().begin ()];
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 begin () const {
                self_type &m = (*this) ();
                return m.find1 (1, 0, index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 end () const {
                self_type &m = (*this) ();
                return m.find1 (1, m.size1 (), index2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_)) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_));
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_)) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 (itv_ - (*this) ().index1_data ().begin (), (*this) ().zero_based (*it_));
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator2 &operator = (const iterator2 &it) {
                container_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_iterator_type itv_;
            iterator_type it_;

            friend class const_iterator2;
        };

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
        size_type non_zeros_;
        size_type filled1_;
        size_type filled2_;
        index_array_type index1_data_;
        index_array_type index2_data_;
        value_array_type value_data_;
        static const value_type zero_;

        BOOST_UBLAS_INLINE
        static size_type zero_based (size_type k_based_index) {
            return k_based_index - IB;
        }
        BOOST_UBLAS_INLINE
        static size_type k_based (size_type zero_based_index) {
            return zero_based_index + IB;
        }

        friend class iterator1;
        friend class iterator2;
        friend class const_iterator1;
        friend class const_iterator2;
    };

    template<class T, class F, std::size_t IB, class IA, class TA>
    const typename compressed_matrix<T, F, IB, IA, TA>::value_type compressed_matrix<T, F, IB, IA, TA>::zero_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME compressed_matrix<T, F, IB, IA, TA>::value_type
#endif
        (0);


    // Array based sparse matrix class
    // Thanks to Kresimir Fresl for extending this to cover different index bases.
    template<class T, class F, std::size_t IB, class IA, class TA>
    class coordinate_matrix:
        public matrix_expression<coordinate_matrix<T, F, IB, IA, TA> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<coordinate_matrix<T, F, IB, IA, TA> >::operator ();
#endif
        // ISSUE require type consistency check for IA TA and IA::value_type
        typedef typename IA::size_type size_type;
        typedef typename IA::difference_type difference_type;
        typedef T value_type;
        typedef const T &const_reference;
#ifndef BOOST_UBLAS_STRICT_MATRIX_SPARSE
        typedef T &reference;
#else
        typedef sparse_matrix_element<coordinate_matrix<T, F, IB, IA, TA> > reference;
#endif
        typedef IA index_array_type;
        typedef TA value_array_type;
    private:
        typedef T &true_reference;
        typedef T *pointer;
        typedef F functor_type;
        typedef coordinate_matrix<T, F, IB, IA, TA> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const matrix_const_reference<const self_type> const_closure_type;
#else
        typedef const matrix_reference<const self_type> const_closure_type;
#endif
        typedef matrix_reference<self_type> closure_type;
        typedef coordinate_vector<T, IB, IA, TA> vector_temporary_type;
        typedef self_type matrix_temporary_type;
        typedef sparse_tag storage_category;
        typedef typename F::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        coordinate_matrix ():
            matrix_expression<self_type> (),
            size1_ (0), size2_ (0), non_zeros_ (restrict_nz (0)),
            filled_ (0),
            sorted_ (true), index1_data_ (non_zeros_),
            index2_data_ (non_zeros_), value_data_ (non_zeros_) {}
        BOOST_UBLAS_INLINE
        coordinate_matrix (size_type size1, size_type size2, size_type non_zeros = 0):
            matrix_expression<self_type> (),
            size1_ (size1), size2_ (size2), non_zeros_ (restrict_nz (non_zeros)),
            filled_ (0),
            sorted_ (true), index1_data_ (non_zeros_),
            index2_data_ (non_zeros_), value_data_ (non_zeros_) {
        }
        BOOST_UBLAS_INLINE
        coordinate_matrix (const coordinate_matrix &m):
            matrix_expression<self_type> (),
            size1_ (m.size1_), size2_ (m.size2_), non_zeros_ (m.non_zeros_),
            filled_ (m.filled_),
            sorted_ (m.sorted_), index1_data_ (m.index1_data_),
            index2_data_ (m.index2_data_), value_data_ (m.value_data_) {
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_matrix (const matrix_expression<AE> &ae, size_type non_zeros = 0):
            matrix_expression<self_type> (),
            size1_ (ae ().size1 ()), size2_ (ae ().size2 ()), non_zeros_ (restrict_nz (non_zeros)),
            filled_ (0),
            sorted_ (true), index1_data_ (non_zeros_),
            index2_data_ (non_zeros_), value_data_ (non_zeros_) {
            matrix_assign (scalar_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
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
        size_type non_zeros () const {
            return non_zeros_;
        }
        BOOST_UBLAS_INLINE
        size_type filled () const {
            return filled_;
        }
        BOOST_UBLAS_INLINE
        static size_type index_base () {
            return IB;
        }
        BOOST_UBLAS_INLINE
        const index_array_type &index1_data () const {
            return index1_data_;
        }
        BOOST_UBLAS_INLINE
        index_array_type &index1_data () {
            return index1_data_;
        }
        BOOST_UBLAS_INLINE
        const index_array_type &index2_data () const {
            return index2_data_;
        }
        BOOST_UBLAS_INLINE
        index_array_type &index2_data () {
            return index2_data_;
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
    private:
        BOOST_UBLAS_INLINE
        size_type restrict_nz (size_type non_zeros) const {
            // minimum non_zeros
            non_zeros = (std::max) (non_zeros, (std::min) (size1_, size2_));
            // ISSUE no maximum as coordinate may contain inserted duplicates
            return non_zeros;
        }
    public:
        BOOST_UBLAS_INLINE
        void resize (size_type size1, size_type size2, bool preserve = true) {
            // FIXME preserve unimplemented
            BOOST_UBLAS_CHECK (!preserve, internal_logic ());
            size1_ = size1;
            size2_ = size2;
            non_zeros_ = restrict_nz (non_zeros_);
            index1_data ().resize (non_zeros_);
            index2_data ().resize (non_zeros_);
            value_data ().resize (non_zeros_);
            filled_ = 0;
            BOOST_UBLAS_CHECK (filled_ <= non_zeros_, internal_logic ());
        }

        // Reserving
        BOOST_UBLAS_INLINE
        void reserve (size_type non_zeros, bool preserve = true) {
            sort ();    // remove duplicate elements
            non_zeros_ = restrict_nz (non_zeros);
            if (preserve) {
                index1_data ().resize (non_zeros_, size_type ());
                index2_data ().resize (non_zeros_, size_type ());
                value_data ().resize (non_zeros_, value_type ());
                filled_ = (std::min) (non_zeros_, filled_);
            }
            else {
                index1_data ().resize (non_zeros_);
                index2_data ().resize (non_zeros_);
                value_data ().resize (non_zeros_);
                filled_ = 0;
            }
            BOOST_UBLAS_CHECK (filled_ <= non_zeros_, internal_logic ());
        }

        // Proxy support
#ifdef BOOST_UBLAS_STRICT_MATRIX_SPARSE
        BOOST_UBLAS_INLINE
        pointer find_element (size_type i, size_type j) {
            sort ();
            size_type element1 (functor_type::element1 (i, size1_, j, size2_));
            size_type element2 (functor_type::element2 (i, size1_, j, size2_));
            vector_const_iterator_type itv_begin (detail::lower_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ()));
            vector_const_iterator_type itv_end (detail::upper_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ()));
            if (itv_begin == itv_end)
                return 0;
            const_iterator_type it_begin (index2_data ().begin () + (itv_begin - index1_data ().begin ()));
            const_iterator_type it_end (index2_data ().begin () + (itv_end - index1_data ().begin ()));
            const_iterator_type it (detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ()));
            if (it == it_end || *it != k_based (element2))
                return 0;
            return &value_data () [it - index2_data ().begin ()];
        }
#endif

        // Element access
        BOOST_UBLAS_INLINE
        const_reference at_element (size_type i, size_type j) const {
            sort ();
            size_type element1 (functor_type::element1 (i, size1_, j, size2_));
            size_type element2 (functor_type::element2 (i, size1_, j, size2_));
            vector_const_iterator_type itv_begin (detail::lower_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ()));
            vector_const_iterator_type itv_end (detail::upper_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ()));
            if (itv_begin == itv_end)
                return zero_;
            const_iterator_type it_begin (index2_data ().begin () + (itv_begin - index1_data ().begin ()));
            const_iterator_type it_end (index2_data ().begin () + (itv_end - index1_data ().begin ()));
            const_iterator_type it (detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ()));
            if (it == it_end || *it != k_based (element2))
                return zero_;
            return value_data () [it - index2_data ().begin ()];
        }
        BOOST_UBLAS_INLINE
        true_reference at_element (size_type i, size_type j) {
            sort ();
            size_type element1 (functor_type::element1 (i, size1_, j, size2_));
            size_type element2 (functor_type::element2 (i, size1_, j, size2_));
            vector_const_iterator_type itv_begin (detail::lower_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ()));
            vector_const_iterator_type itv_end (detail::upper_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ()));
            if (itv_begin == itv_end) {
                insert (i, j, value_type (0));
                sort ();
                itv_begin = detail::lower_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ());
                itv_end = detail::upper_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ());
            }
            const_iterator_type it_begin (index2_data ().begin () + (itv_begin - index1_data ().begin ()));
            const_iterator_type it_end (index2_data ().begin () + (itv_end - index1_data ().begin ()));
            const_iterator_type it (detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ()));
            if (it == it_end || *it != k_based (element2)) {
                insert (i, j, value_type (0));
                sort ();
                itv_begin = detail::lower_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ());
                itv_end = detail::upper_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ());
                it_begin = index2_data ().begin () + (itv_begin - index1_data ().begin ());
                it_end = index2_data ().begin () + (itv_end - index1_data ().begin ());
                it = detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ());
            }
            return value_data () [it - index2_data ().begin ()];
        }
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return at_element (i, j);
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
#ifndef BOOST_UBLAS_STRICT_MATRIX_SPARSE
            return at_element (i, j);
#else
            return reference (*this, i, j);
#endif
        }

        // Assignment
        BOOST_UBLAS_INLINE
        coordinate_matrix &operator = (const coordinate_matrix &m) {
            if (this != &m) {
                size1_ = m.size1_;
                size2_ = m.size2_;
                non_zeros_ = m.non_zeros_;
                filled_ = m.filled_;
                sorted_ = m.sorted_;
                index1_data () = m.index1_data ();
                index2_data () = m.index2_data ();
                value_data () = m.value_data ();
                BOOST_UBLAS_CHECK (non_zeros_ == index1_data ().size (), internal_logic ());
                BOOST_UBLAS_CHECK (non_zeros_ == index2_data ().size (), internal_logic ());
                BOOST_UBLAS_CHECK (non_zeros_ == value_data ().size (), internal_logic ());
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        coordinate_matrix &assign_temporary (coordinate_matrix &m) {
            swap (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_matrix &operator = (const matrix_expression<AE> &ae) {
            // return assign (self_type (ae, non_zeros_));
            self_type temporary (ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_matrix &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_matrix& operator += (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this + ae, non_zeros_));
            self_type temporary (*this + ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_matrix &plus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_plus_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_matrix& operator -= (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this - ae, non_zeros_));
            self_type temporary (*this - ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        coordinate_matrix &minus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_minus_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        coordinate_matrix& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<true_reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        coordinate_matrix& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<true_reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (coordinate_matrix &m) {
            if (this != &m) {
                std::swap (size1_, m.size1_);
                std::swap (size2_, m.size2_);
                std::swap (non_zeros_, m.non_zeros_);
                std::swap (filled_, m.filled_);
                std::swap (sorted_, m.sorted_);
                index1_data ().swap (m.index1_data ());
                index2_data ().swap (m.index2_data ());
                value_data ().swap (m.value_data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (coordinate_matrix &m1, coordinate_matrix &m2) {
            m1.swap (m2);
        }
#endif

        // Sorting and remove duplicates
        BOOST_UBLAS_INLINE
        void sort () const {
            if (! sorted_ && filled_ > 0) {
                index_triple_array<index_array_type, index_array_type, value_array_type>
                    ita (filled_, index1_data_, index2_data_, value_data_);
                std::sort (ita.begin (), ita.end ());
                // ISSUE: unusual semantics - sum values of duplicates
                size_type filled = 1;
                for (size_type i = 1; i < filled_; ++ i) {
                    if (index1_data_ [filled - 1] != index1_data_ [i] ||
                        index2_data_ [filled - 1] != index2_data_ [i]) {
                        ++ filled;
                        if (filled - 1 != i) {
                            index1_data_ [filled - 1] = index1_data_ [i];
                            index2_data_ [filled - 1] = index2_data_ [i];
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
        void push_back (size_type i, size_type j, const_reference t) {
            size_type element1 = functor_type::element1 (i, size1_, j, size2_);
            size_type element2 = functor_type::element2 (i, size1_, j, size2_);
            if (filled_ == 0 ||
                index1_data () [filled_ - 1] < k_based (element1) ||
                (index1_data () [filled_ - 1] == k_based (element1) && index2_data () [filled_ - 1] < k_based (element2))) {
                if (filled_ >= non_zeros_)
                    reserve (2 * filled_, true);
                BOOST_UBLAS_CHECK (filled_ < non_zeros_, internal_logic ());
                index1_data () [filled_] = k_based (element1);
                index2_data () [filled_] = k_based (element2);
                value_data () [filled_] = t;
                ++ filled_;
                return;
            }
            external_logic ().raise ();
        }
        BOOST_UBLAS_INLINE
        void insert (size_type i, size_type j, const_reference t) {
            if (filled_ >= non_zeros_)
                reserve (2 * filled_, true);
            BOOST_UBLAS_CHECK (filled_ < non_zeros_, internal_logic ());
            size_type element1 = functor_type::element1 (i, size1_, j, size2_);
            size_type element2 = functor_type::element2 (i, size1_, j, size2_);
            index1_data () [filled_] = k_based (element1);
            index2_data () [filled_] = k_based (element2);
            value_data () [filled_] = t;
            ++ filled_;
            sorted_ = false;
        }
        BOOST_UBLAS_INLINE
        void pop_back () {
            BOOST_UBLAS_CHECK (filled_ > 0, external_logic ());
            -- filled_;
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i, size_type j) {
            size_type element1 = functor_type::element1 (i, size1_, j, size2_);
            size_type element2 = functor_type::element2 (i, size1_, j, size2_);
            sort ();
            vector_iterator_type itv_begin (detail::lower_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ()));
            vector_iterator_type itv_end (detail::upper_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (element1), std::less<size_type> ()));
            iterator_type it_begin (index2_data ().begin () + (itv_begin - index1_data ().begin ()));
            iterator_type it_end (index2_data ().begin () + (itv_end - index1_data ().begin ()));
            iterator_type it (detail::lower_bound (it_begin, it_end, k_based (element2), std::less<size_type> ()));
            if (it != it_end && *it == k_based (element2)) {
                difference_type n = it - index2_data ().begin ();
                vector_iterator_type itv (index1_data ().begin () + n);
                std::copy (itv + 1, index1_data ().begin () + filled_, itv);
                std::copy (it + 1, index2_data ().begin () + filled_, it);
                typename value_array_type::iterator itt (value_data ().begin () + n);
                std::copy (itt + 1, value_data ().begin () + filled_, itt);
                -- filled_;
            }
        }
        BOOST_UBLAS_INLINE
        void clear () {
            filled_ = 0;
        }

        // Iterator types
    private:
        // Use index array iterator
        typedef typename IA::const_iterator vector_const_iterator_type;
        typedef typename IA::iterator vector_iterator_type;
        typedef typename IA::const_iterator const_iterator_type;
        typedef typename IA::iterator iterator_type;

    public:
        class const_iterator1;
        class iterator1;
        class const_iterator2;
        class iterator2;
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
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j, int direction = 1) const {
            sort ();
            for (;;) {
                size_type address1 (functor_type::address1 (i, size1_, j, size2_));
                size_type address2 (functor_type::address2 (i, size1_, j, size2_));
                vector_const_iterator_type itv_begin (detail::lower_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (address1), std::less<size_type> ()));
                vector_const_iterator_type itv_end (detail::upper_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (address1), std::less<size_type> ()));

                const_iterator_type it_begin (index2_data ().begin () + (itv_begin - index1_data ().begin ()));
                const_iterator_type it_end (index2_data ().begin () + (itv_end - index1_data ().begin ()));

                const_iterator_type it (detail::lower_bound (it_begin, it_end, k_based (address2), std::less<size_type> ()));
                vector_const_iterator_type itv (index1_data ().begin () + (it - index2_data ().begin ()));
                if (rank == 0)
                    return const_iterator1 (*this, rank, i, j, itv, it);
                if (it != it_end && zero_based (*it) == address2)
                    return const_iterator1 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast1 ()) {
                        if (it == it_end)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        i = zero_based (*it);
                    } else {
                        if (i >= size1_)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        ++ i;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast1 ()) {
                        if (it == index2_data ().begin () + zero_based (*itv))
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        i = zero_based (*(it - 1));
                    } else {
                        if (i == 0)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        -- i;
                    }
                }
            }
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j, int direction = 1) {
            sort ();
            for (;;) {
                size_type address1 (functor_type::address1 (i, size1_, j, size2_));
                size_type address2 (functor_type::address2 (i, size1_, j, size2_));
                vector_iterator_type itv_begin (detail::lower_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (address1), std::less<size_type> ()));
                vector_iterator_type itv_end (detail::upper_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (address1), std::less<size_type> ()));

                iterator_type it_begin (index2_data ().begin () + (itv_begin - index1_data ().begin ()));
                iterator_type it_end (index2_data ().begin () + (itv_end - index1_data ().begin ()));

                iterator_type it (detail::lower_bound (it_begin, it_end, k_based (address2), std::less<size_type> ()));
                vector_iterator_type itv (index1_data ().begin () + (it - index2_data ().begin ()));
                if (rank == 0)
                    return iterator1 (*this, rank, i, j, itv, it);
                if (it != it_end && zero_based (*it) == address2)
                    return iterator1 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast1 ()) {
                        if (it == it_end)
                            return iterator1 (*this, rank, i, j, itv, it);
                        i = zero_based (*it);
                    } else {
                        if (i >= size1_)
                            return iterator1 (*this, rank, i, j, itv, it);
                        ++ i;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast1 ()) {
                        if (it == index2_data ().begin () + zero_based (*itv))
                            return iterator1 (*this, rank, i, j, itv, it);
                        i = zero_based (*(it - 1));
                    } else {
                        if (i == 0)
                            return iterator1 (*this, rank, i, j, itv, it);
                        -- i;
                    }
                }
            }
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j, int direction = 1) const {
            sort ();
            for (;;) {
                size_type address1 (functor_type::address1 (i, size1_, j, size2_));
                size_type address2 (functor_type::address2 (i, size1_, j, size2_));
                vector_const_iterator_type itv_begin (detail::lower_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (address1), std::less<size_type> ()));
                vector_const_iterator_type itv_end (detail::upper_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (address1), std::less<size_type> ()));

                const_iterator_type it_begin (index2_data ().begin () + (itv_begin - index1_data ().begin ()));
                const_iterator_type it_end (index2_data ().begin () + (itv_end - index1_data ().begin ()));

                const_iterator_type it (detail::lower_bound (it_begin, it_end, k_based (address2), std::less<size_type> ()));
                vector_const_iterator_type itv (index1_data ().begin () + (it - index2_data ().begin ()));
                if (rank == 0)
                    return const_iterator2 (*this, rank, i, j, itv, it);
                if (it != it_end && zero_based (*it) == address2)
                    return const_iterator2 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast2 ()) {
                        if (it == it_end)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        j = zero_based (*it);
                    } else {
                        if (j >= size2_)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        ++ j;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast2 ()) {
                        if (it == index2_data ().begin () + zero_based (*itv))
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        j = zero_based (*(it - 1));
                    } else {
                        if (j == 0)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        -- j;
                    }
                }
            }
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j, int direction = 1) {
            sort ();
            for (;;) {
                size_type address1 (functor_type::address1 (i, size1_, j, size2_));
                size_type address2 (functor_type::address2 (i, size1_, j, size2_));
                vector_iterator_type itv_begin (detail::lower_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (address1), std::less<size_type> ()));
                vector_iterator_type itv_end (detail::upper_bound (index1_data ().begin (), index1_data ().begin () + filled_, k_based (address1), std::less<size_type> ()));

                iterator_type it_begin (index2_data ().begin () + (itv_begin - index1_data ().begin ()));
                iterator_type it_end (index2_data ().begin () + (itv_end - index1_data ().begin ()));

                iterator_type it (detail::lower_bound (it_begin, it_end, k_based (address2), std::less<size_type> ()));
                vector_iterator_type itv (index1_data ().begin () + (it - index2_data ().begin ()));
                if (rank == 0)
                    return iterator2 (*this, rank, i, j, itv, it);
                if (it != it_end && zero_based (*it) == address2)
                    return iterator2 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast2 ()) {
                        if (it == it_end)
                            return iterator2 (*this, rank, i, j, itv, it);
                        j = zero_based (*it);
                    } else {
                        if (j >= size2_)
                            return iterator2 (*this, rank, i, j, itv, it);
                        ++ j;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast2 ()) {
                        if (it == index2_data ().begin () + zero_based (*itv))
                            return iterator2 (*this, rank, i, j, itv, it);
                        j = zero_based (*(it - 1));
                    } else {
                        if (j == 0)
                            return iterator2 (*this, rank, i, j, itv, it);
                        -- j;
                    }
                }
            }
        }


        class const_iterator1:
            public container_const_reference<coordinate_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator1, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename coordinate_matrix::value_type value_type;
            typedef typename coordinate_matrix::difference_type difference_type;
            typedef typename coordinate_matrix::const_reference reference;
            typedef const typename coordinate_matrix::pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &m, int rank, size_type i, size_type j, const vector_const_iterator_type &itv, const const_iterator_type &it):
                container_const_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const iterator1 &it):
                container_const_reference<self_type> (it ()), rank_ (it.rank_), i_ (it.i_), j_ (it.j_), itv_ (it.itv_), it_ (it.it_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                if (rank_ == 1 && functor_type::fast1 ())
                    ++ it_;
                else {
                    i_ = index1 () + 1;
                    if (rank_ == 1)
                        *this = (*this) ().find1 (rank_, i_, j_, 1);
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                if (rank_ == 1 && functor_type::fast1 ())
                    -- it_;
                else {
                    i_ = index1 () - 1;
                    if (rank_ == 1)
                        *this = (*this) ().find1 (rank_, i_, j_, -1);
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*this) ().value_data () [it_ - (*this) ().index2_data ().begin ()];
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 begin () const {
                const self_type &m = (*this) ();
                return m.find2 (1, index1 (), 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator2 end () const {
                const self_type &m = (*this) ();
                return m.find2 (1, index1 (), m.size2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_)) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_));
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_)) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_));
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_const_iterator_type itv_;
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator1 begin1 () const {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator1 end1 () const {
            return find1 (0, size1_, 0);
        }

        class iterator1:
            public container_reference<coordinate_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator1, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename coordinate_matrix::value_type value_type;
            typedef typename coordinate_matrix::difference_type difference_type;
            typedef typename coordinate_matrix::true_reference reference;
            typedef typename coordinate_matrix::pointer pointer;
#endif
            typedef iterator2 dual_iterator_type;
            typedef reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator1 ():
                container_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator1 (self_type &m, int rank, size_type i, size_type j, const vector_iterator_type &itv, const iterator_type &it):
                container_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator1 &operator ++ () {
                if (rank_ == 1 && functor_type::fast1 ())
                    ++ it_;
                else {
                    i_ = index1 () + 1;
                    if (rank_ == 1)
                        *this = (*this) ().find1 (rank_, i_, j_, 1);
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator1 &operator -- () {
                if (rank_ == 1 && functor_type::fast1 ())
                    -- it_;
                else {
                    i_ = index1 () - 1;
                    if (rank_ == 1)
                        *this = (*this) ().find1 (rank_, i_, j_, -1);
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*this) ().value_data () [it_ - (*this) ().index2_data ().begin ()];
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 begin () const {
                self_type &m = (*this) ();
                return m.find2 (1, index1 (), 0);
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator2 end () const {
                self_type &m = (*this) ();
                return m.find2 (1, index1 (), m.size2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_)) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_));
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_)) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_));
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator1 &operator = (const iterator1 &it) {
                container_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_iterator_type itv_;
            iterator_type it_;

            friend class const_iterator1;
        };

        BOOST_UBLAS_INLINE
        iterator1 begin1 () {
            return find1 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        iterator1 end1 () {
            return find1 (0, size1_, 0);
        }

        class const_iterator2:
            public container_const_reference<coordinate_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator2, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename coordinate_matrix::value_type value_type;
            typedef typename coordinate_matrix::difference_type difference_type;
            typedef typename coordinate_matrix::const_reference reference;
            typedef const typename coordinate_matrix::pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &m, int rank, size_type i, size_type j, const vector_const_iterator_type itv, const const_iterator_type &it):
                container_const_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const iterator2 &it):
                container_const_reference<self_type> (it ()), rank_ (it.rank_), i_ (it.i_), j_ (it.j_), itv_ (it.itv_), it_ (it.it_) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                if (rank_ == 1 && functor_type::fast2 ())
                    ++ it_;
                else {
                    j_ = index2 () + 1;
                    if (rank_ == 1)
                        *this = (*this) ().find2 (rank_, i_, j_, 1);
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                if (rank_ == 1 && functor_type::fast2 ())
                    -- it_;
                else {
                    j_ = index2 () - 1;
                    if (rank_ == 1)
                        *this = (*this) ().find2 (rank_, i_, j_, -1);
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*this) ().value_data () [it_ - (*this) ().index2_data ().begin ()];
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 begin () const {
                const self_type &m = (*this) ();
                return m.find1 (1, 0, index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            const_iterator1 end () const {
                const self_type &m = (*this) ();
                return m.find1 (1, m.size1 (), index2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_)) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_));
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_)) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_));
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_const_iterator_type itv_;
            const_iterator_type it_;
        };

        BOOST_UBLAS_INLINE
        const_iterator2 begin2 () const {
            return find2 (0, 0, 0);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 end2 () const {
            return find2 (0, 0, size2_);
        }

        class iterator2:
            public container_reference<coordinate_matrix>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator2, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename coordinate_matrix::value_type value_type;
            typedef typename coordinate_matrix::difference_type difference_type;
            typedef typename coordinate_matrix::true_reference reference;
            typedef typename coordinate_matrix::pointer pointer;
#endif
            typedef iterator1 dual_iterator_type;
            typedef reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            iterator2 ():
                container_reference<self_type> (), rank_ (), i_ (), j_ (), itv_ (), it_ () {}
            BOOST_UBLAS_INLINE
            iterator2 (self_type &m, int rank, size_type i, size_type j, const vector_iterator_type &itv, const iterator_type &it):
                container_reference<self_type> (m), rank_ (rank), i_ (i), j_ (j), itv_ (itv), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            iterator2 &operator ++ () {
                if (rank_ == 1 && functor_type::fast2 ())
                    ++ it_;
                else {
                    j_ = index2 () + 1;
                    if (rank_ == 1)
                        *this = (*this) ().find2 (rank_, i_, j_, 1);
                }
                return *this;
            }
            BOOST_UBLAS_INLINE
            iterator2 &operator -- () {
                if (rank_ == 1 && functor_type::fast2 ())
                    -- it_;
                else {
                    j_ = index2 ();
                    if (rank_ == 1)
                        *this = (*this) ().find2 (rank_, i_, j_, -1);
                }
                return *this;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
                BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
                if (rank_ == 1) {
                    return (*this) ().value_data () [it_ - (*this) ().index2_data ().begin ()];
                } else {
                    return (*this) ().at_element (i_, j_);
                }
            }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 begin () const {
                self_type &m = (*this) ();
                return m.find1 (1, 0, index2 ());
            }
            BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
            typename self_type::
#endif
            iterator1 end () const {
                self_type &m = (*this) ();
                return m.find1 (1, m.size1 (), index2 ());
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
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index1 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_)) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_));
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_)) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 ((*this) ().zero_based (*itv_), (*this) ().zero_based (*it_));
                } else {
                    return j_;
                }
            }

            // Assignment
            BOOST_UBLAS_INLINE
            iterator2 &operator = (const iterator2 &it) {
                container_reference<self_type>::assign (&it ());
                rank_ = it.rank_;
                i_ = it.i_;
                j_ = it.j_;
                itv_ = it.itv_;
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // BOOST_UBLAS_CHECK (rank_ == it.rank_, internal_logic ());
                if (rank_ == 1 || it.rank_ == 1) {
                    return it_ == it.it_;
                } else {
                    return i_ == it.i_ && j_ == it.j_;
                }
            }

        private:
            int rank_;
            size_type i_;
            size_type j_;
            vector_iterator_type itv_;
            iterator_type it_;

            friend class const_iterator2;
        };

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
        size_type non_zeros_;
        mutable size_type filled_;
        mutable bool sorted_;
        mutable index_array_type index1_data_;
        mutable index_array_type index2_data_;
        mutable value_array_type value_data_;
        static const value_type zero_;

        BOOST_UBLAS_INLINE
        static size_type zero_based (size_type k_based_index) {
            return k_based_index - IB;
        }
        BOOST_UBLAS_INLINE
        static size_type k_based (size_type zero_based_index) {
            return zero_based_index + IB;
        }

        friend class iterator1;
        friend class iterator2;
        friend class const_iterator1;
        friend class const_iterator2;
    };

    template<class T, class F, std::size_t IB, class IA, class TA>
    const typename coordinate_matrix<T, F, IB, IA, TA>::value_type coordinate_matrix<T, F, IB, IA, TA>::zero_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME coordinate_matrix<T, F, IB, IA, TA>::value_type
#endif
        (0);

}}}

#endif
