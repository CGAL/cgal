//
//  Copyright (c) 2003
//  Gunter Winkler, Joerg Walter
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
#ifndef BOOST_UBLAS_ENABLE_EXPERIMENTAL
#error class generalized_vector_of_vector is experiment and currently does not work
#endif

#ifndef BOOST_UBLAS_VECTOR_OF_VECTOR_H
#define BOOST_UBLAS_VECTOR_OF_VECTOR_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/storage_sparse.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

// Iterators based on ideas of Jeremy Siek

namespace boost { namespace numeric { namespace ublas {

    // Array based sparse matrix class
    template<class T, class F, class A>
    class generalized_vector_of_vector:
        public matrix_expression<generalized_vector_of_vector<T, F, A> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<generalized_vector_of_vector<T, F, A> >::operator ();
#endif
        typedef typename A::size_type size_type;
        typedef typename A::difference_type difference_type;
        typedef T value_type;
        typedef const T &const_reference;
#ifndef BOOST_UBLAS_STRICT_VECTOR_SPARSE
        typedef T &reference;
#else
        typedef sparse_vector_element<typename A::value_type> reference;
#endif
        typedef A array_type;
    private:
        typedef T *pointer;
        typedef F functor_type;
        typedef generalized_vector_of_vector<T, F, A> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const matrix_const_reference<const self_type> const_closure_type;
#else
        typedef const matrix_reference<const self_type> const_closure_type;
#endif
        typedef matrix_reference<self_type> closure_type;
        typedef typename A::value_type vector_data_value_type;
        typedef sparse_tag storage_category;
        typedef typename F::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector ():
            size1_ (0), size2_ (0), non_zeros_ (0), data_ (1) {
            for (size_type i = 0; i < functor_type::size1 (size1_, size2_); ++ i)
                static_cast<vector_data_value_type &> (data_ [i]).resize (functor_type::size2 (size1_, size2_));
            data_ [functor_type::size1 (size1_, size2_)] = vector_data_value_type ();
        }
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector (size_type size1, size_type size2, size_type non_zeros = 0):
            size1_ (size1), size2_ (size2), non_zeros_ (non_zeros), data_ (functor_type::size1 (size1_, size2_) + 1) {
            for (size_type i = 0; i < functor_type::size1 (size1_, size2_); ++ i)
                static_cast<vector_data_value_type &> (data_ [i]).resize (functor_type::size2 (size1_, size2_));
            data_ [functor_type::size1 (size1_, size2_)] = vector_data_value_type ();
        }
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector (const generalized_vector_of_vector &m):
            size1_ (m.size1_), size2_ (m.size2_), non_zeros_ (m.non_zeros_), data_ (m.data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector (const matrix_expression<AE> &ae, size_type non_zeros = 0):
            size1_ (ae ().size1 ()), size2_ (ae ().size2 ()), non_zeros_ (non_zeros), data_ (functor_type::size1 (size1_, size2_) + 1) {
            for (size_type i = 0; i < functor_type::size1 (size1_, size2_); ++ i)
                static_cast<vector_data_value_type &> (data_ [i]).resize (functor_type::size2 (size1_, size2_));
            data_ [functor_type::size1 (size1_, size2_)] = vector_data_value_type ();
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
        size_type non_zeros () const {
            size_type non_zeros = 0;
            for (vector_const_iterator_type itv = data_ ().begin (); itv != data_ ().end (); ++ itv)
                non_zeros += (*itv).size ();
            return non_zeros;
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
        void resize (size_type size1, size_type size2, size_type non_zeros = 0) {
            size1_ = size1;
            size2_ = size2;
            non_zeros_ = non_zeros;
            data ().resize (functor_type::size1 (size1_, size2_) + 1);
            for (size_type i = 0; i < functor_type::size1 (size1_, size2_); ++ i)
                static_cast<vector_data_value_type &> (data_ [i]).resize (functor_type::size2 (size1_, size2_));
            data () [functor_type::size1 (size1_, size2_)] = vector_data_value_type ();
        }

        // Proxy support
#ifdef BOOST_UBLAS_STRICT_VECTOR_SPARSE
        pointer find_element (size_type i, size_type j) {
            vector_iterator_type itv (data ().find (functor_type::element1 (i, size1_, j, size2_)));
            if (itv == data ().end () || itv.index () != functor_type::element1 (i, size1_, j, size2_))
                return 0;
            iterator_type it (static_cast<vector_data_value_type &> (*itv).find (functor_type::element2 (i, size1_, j, size2_)));
            if (it == static_cast<vector_data_value_type &> (*itv).end () || it.index () != functor_type::element2 (i, size1_, j, size2_))
                return 0;
            return &static_cast<value_type &> (*it);
        }
#endif

        // Element access
        BOOST_UBLAS_INLINE
        const_reference at_element (size_type i, size_type j) const {
            vector_const_iterator_type itv (data ().find (functor_type::element1 (i, size1_, j, size2_)));
            if (itv == data ().end () || itv.index () != functor_type::element1 (i, size1_, j, size2_))
                return zero_;
            const_iterator_type it (static_cast<const vector_data_value_type &> (*itv).find (functor_type::element2 (i, size1_, j, size2_)));
            if (it == static_cast<const vector_data_value_type &> (*itv).end () || it.index () != functor_type::element2 (i, size1_, j, size2_))
                return zero_;
            return static_cast<const value_type &> (*it);
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
            return reference (this->data () [functor_type::element1 (i, size1_, j, size2_)], functor_type::element2 (i, size1_, j, size2_));
#endif
        }

        // Assignment
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector &operator = (const generalized_vector_of_vector &m) {
            if (this != &m) {
                size1_ = m.size1_;
                size2_ = m.size2_;
                non_zeros_ = m.non_zeros_;
                data () = m.data ();
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector &assign_temporary (generalized_vector_of_vector &m) {
            swap (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector &operator = (const matrix_expression<AE> &ae) {
            // return assign (self_type (ae, non_zeros_));
            self_type temporary (ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae); 
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector& operator += (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this + ae, non_zeros_));
            self_type temporary (*this + ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector &plus_assign (const matrix_expression<AE> &ae) { 
            matrix_assign (scalar_plus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector& operator -= (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this - ae, non_zeros_));
            self_type temporary (*this - ae, non_zeros_);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector &minus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_minus_assign<reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector& operator *= (const AT &at) {
            matrix_assign_scalar (scalar_multiplies_assign<reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        generalized_vector_of_vector& operator /= (const AT &at) {
            matrix_assign_scalar (scalar_divides_assign<reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (generalized_vector_of_vector &m) {
            if (this != &m) {
                std::swap (size1_, m.size1_);
                std::swap (size2_, m.size2_);
                std::swap (non_zeros_, m.non_zeros_);
                data ().swap (m.data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (generalized_vector_of_vector &m1, generalized_vector_of_vector &m2) {
            m1.swap (m2);
        }
#endif

        // Sorting
        void sort () {
            vector_iterator_type itv (data ().begin ());
            vector_iterator_type itv_end (data ().end ());
            while (itv != itv_end) {
                (*itv).sort ();
                ++ itv;
            }
        }

        // Element insertion and erasure
        BOOST_UBLAS_INLINE
        void insert (size_type i, size_type j, const_reference t) {
            vector_iterator_type itv (data ().find (functor_type::element1 (i, size1_, j, size2_)));
            if (itv == data ().end ()) {
                data ().insert (functor_type::element1 (i, size1_, j, size2_), vector_data_value_type (functor_type::size2 (size1_, size2_)));
                itv = data ().find (functor_type::element1 (i, size1_, j, size2_));
            }
            // FIXME: should be allowed for coordinate_vector.
            // BOOST_UBLAS_CHECK (static_cast<vector_data_value_type &> (*itv).find (functor_type::element2 (i, size1_, j, size2_)) == static_cast<vector_data_value_type &> (*itv).end (), bad_index ());
            static_cast<vector_data_value_type &> (*itv).insert (functor_type::element2 (i, size1_, j, size2_), t);
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i, size_type j) {
            vector_iterator_type itv (data ().find (functor_type::element1 (i, size1_, j, size2_)));
            if (itv == data ().end ())
                return;
            static_cast<vector_data_value_type &> (*itv).erase (functor_type::element2 (i, size1_, j, size2_));
        }
        BOOST_UBLAS_INLINE
        void clear () {
            data ().resize (functor_type::size1 (size1_, size2_) + 1);
            for (size_type i = 0; i < functor_type::size1 (size1_, size2_); ++ i)
                static_cast<vector_data_value_type &> (data_ [i]).resize (functor_type::size2 (size1_, size2_));
            data_ [functor_type::size1 (size1_, size2_)] = vector_data_value_type ();
        }

        // Iterator types
    private:
        // Use vector iterator
        typedef typename A::const_iterator vector_const_iterator_type;
        typedef typename A::iterator vector_iterator_type;
        typedef typename A::value_type::const_iterator const_iterator_type;
        typedef typename A::value_type::iterator iterator_type;

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
                vector_const_iterator_type itv (data ().find (functor_type::address1 (i, size1_, j, size2_)));
                vector_const_iterator_type itv_end (data ().end ());
                if (itv == itv_end)
                    return const_iterator1 (*this, rank, i, j, itv_end, static_cast<const vector_data_value_type &> (*(-- itv)).end ());

                const_iterator_type it (static_cast<const vector_data_value_type &> (*itv).find (functor_type::address2 (i, size1_, j, size2_)));
                const_iterator_type it_end (static_cast<const vector_data_value_type &> (*itv).end ());
                if (rank == 0)
                    return const_iterator1 (*this, rank, i, j, itv, it);
                if (it != it_end && *it == functor_type::address2 (i, size1_, j, size2_))
                    return const_iterator1 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast1 ()) {
                        if (it == it_end)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        i = *it;
                    } else {
                        if (i >= size1_)
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        ++ i;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast1 ()) {
                        if (it == static_cast<const vector_data_value_type &> (*itv).begin ())
                            return const_iterator1 (*this, rank, i, j, itv, it);
                        i = *(it - 1);
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
                vector_iterator_type itv (data ().find (functor_type::address1 (i, size1_, j, size2_)));
                vector_iterator_type itv_end (data ().end ());
                if (itv == itv_end)
                    return iterator1 (*this, rank, i, j, itv_end, static_cast<vector_data_value_type &> (*(-- itv)).end ());

                iterator_type it (static_cast<vector_data_value_type &> (*itv).find (functor_type::address2 (i, size1_, j, size2_)));
                iterator_type it_end (static_cast<vector_data_value_type &> (*itv).end ());
                if (rank == 0)
                    return iterator1 (*this, rank, i, j, itv, it);
                if (it != it_end && *it == functor_type::address2 (i, size1_, j, size2_))
                    return iterator1 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast1 ()) {
                        if (it == it_end)
                            return iterator1 (*this, rank, i, j, itv, it);
                        i = *it;
                    } else {
                        if (i >= size1_)
                            return iterator1 (*this, rank, i, j, itv, it);
                        ++ i;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast1 ()) {
                        if (it == static_cast<const vector_data_value_type &> (*itv).begin ())
                            return iterator1 (*this, rank, i, j, itv, it);
                        i = *(it - 1);
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
                vector_const_iterator_type itv (data ().find (functor_type::address1 (i, size1_, j, size2_)));
                vector_const_iterator_type itv_end (data ().end ());
                if (itv == itv_end)
                    return const_iterator2 (*this, rank, i, j, itv_end, static_cast<const vector_data_value_type &> (*(-- itv)).end ());

                const_iterator_type it (static_cast<const vector_data_value_type &> (*itv).find (functor_type::address2 (i, size1_, j, size2_)));
                const_iterator_type it_end (static_cast<const vector_data_value_type &> (*itv).end ());
                if (rank == 0)
                    return const_iterator2 (*this, rank, i, j, itv, it);
                if (it != it_end && *it == functor_type::address2 (i, size1_, j, size2_))
                    return const_iterator2 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast2 ()) {
                        if (it == it_end)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        j = *it;
                    } else {
                        if (j >= size2_)
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        ++ j;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast2 ()) {
                        if (it == static_cast<const vector_data_value_type &> (*itv).begin ())
                            return const_iterator2 (*this, rank, i, j, itv, it);
                        j = *(it - 1);
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
                vector_iterator_type itv (data ().find (functor_type::address1 (i, size1_, j, size2_)));
                vector_iterator_type itv_end (data ().end ());
                if (itv == itv_end)
                    return iterator2 (*this, rank, i, j, itv_end, static_cast<vector_data_value_type &> (*(-- itv)).end ());

                iterator_type it (static_cast<vector_data_value_type &> (*itv).find (functor_type::address2 (i, size1_, j, size2_)));
                iterator_type it_end (static_cast<vector_data_value_type &> (*itv).end ());
                if (rank == 0)
                    return iterator2 (*this, rank, i, j, itv, it);
                if (it != it_end && *it == functor_type::address2 (i, size1_, j, size2_))
                    return iterator2 (*this, rank, i, j, itv, it);
                if (direction > 0) {
                    if (functor_type::fast2 ()) {
                        if (it == it_end)
                            return iterator2 (*this, rank, i, j, itv, it);
                        j = *it;
                    } else {
                        if (j >= size2_)
                            return iterator2 (*this, rank, i, j, itv, it);
                        ++ j;
                    }
                } else /* if (direction < 0)  */ {
                    if (functor_type::fast2 ()) {
                        if (it == static_cast<const vector_data_value_type &> (*itv).begin ())
                            return iterator2 (*this, rank, i, j, itv, it);
                        j = *(it - 1);
                    } else {
                        if (j == 0)
                            return iterator2 (*this, rank, i, j, itv, it);
                        -- j;
                    }
                }
            }
        }


        class const_iterator1:
            public container_const_reference<generalized_vector_of_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator1, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename generalized_vector_of_vector::difference_type difference_type;
            typedef typename generalized_vector_of_vector::value_type value_type;
            typedef typename generalized_vector_of_vector::const_reference reference;
            typedef const typename generalized_vector_of_vector::pointer pointer;
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
                        it_ = static_cast<const vector_data_value_type &> (*itv_).begin ();
                        if (it_ == static_cast<const vector_data_value_type &> (*itv_).end () || index2 () != j_)
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
                        it_ = static_cast<const vector_data_value_type &> (*itv_).begin ();
                        if (it_ == static_cast<const vector_data_value_type &> (*itv_).end () || index2 () != j_)
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
                    return static_cast<const value_type &> (*it_);
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
                    BOOST_UBLAS_CHECK (functor_type::index1 (itv_.index (), it_.index ()) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 (itv_.index (), it_.index ());
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 (itv_.index (), it_.index ()) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 (itv_.index (), it_.index ());
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
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
            public container_reference<generalized_vector_of_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator1, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename generalized_vector_of_vector::difference_type difference_type;
            typedef typename generalized_vector_of_vector::value_type value_type;
            typedef typename generalized_vector_of_vector::true_reference reference;
            typedef typename generalized_vector_of_vector::pointer pointer;
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
                        it_ = static_cast<vector_data_value_type &> (*itv_).begin ();
                        if (it_ == static_cast<vector_data_value_type &> (*itv_).end () || index2 () != j_)
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
                        it_ = static_cast<vector_data_value_type &> (*itv_).begin ();
                        if (it_ == static_cast<vector_data_value_type &> (*itv_).end () || index2 () != j_)
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
                    return static_cast<value_type &> (*it_);
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
                    BOOST_UBLAS_CHECK (functor_type::index1 (itv_.index (), it_.index ()) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 (itv_.index (), it_.index ());
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find1 (0, (*this) ().size1 (), j_), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 (itv_.index (), it_.index ()) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 (itv_.index (), it_.index ());
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
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
            public container_const_reference<generalized_vector_of_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               const_iterator2, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename generalized_vector_of_vector::difference_type difference_type;
            typedef typename generalized_vector_of_vector::value_type value_type;
            typedef typename generalized_vector_of_vector::const_reference reference;
            typedef const typename generalized_vector_of_vector::pointer pointer;
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
                        it_ = static_cast<const vector_data_value_type &> (*itv_).begin ();
                        if (it_ == static_cast<const vector_data_value_type &> (*itv_).end () || index1 () != i_)
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
                        it_ = static_cast<const vector_data_value_type &> (*itv_).begin ();
                        if (it_ == static_cast<const vector_data_value_type &> (*itv_).end () || index1 () != i_)
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
                    return static_cast<const value_type &> (*it_);
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
                    BOOST_UBLAS_CHECK (functor_type::index1 (itv_.index (), it_.index ()) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 (itv_.index (), it_.index ());
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 (itv_.index (), it_.index ()) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 (itv_.index (), it_.index ());
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
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
            public container_reference<generalized_vector_of_vector>,
            public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                               iterator2, value_type> {
        public:
            typedef sparse_bidirectional_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename generalized_vector_of_vector::difference_type difference_type;
            typedef typename generalized_vector_of_vector::value_type value_type;
            typedef typename generalized_vector_of_vector::true_reference reference;
            typedef typename generalized_vector_of_vector::pointer pointer;
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
                        it_ = static_cast<vector_data_value_type &> (*itv_).begin ();
                        if (it_ == static_cast<vector_data_value_type &> (*itv_).end () || index1 () != i_)
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
                        it_ = static_cast<vector_data_value_type &> (*itv_).begin ();
                        if (it_ == static_cast<vector_data_value_type &> (*itv_).end () || index1 () != i_)
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
                    return static_cast<value_type &> (*it_);
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
                    BOOST_UBLAS_CHECK (functor_type::index1 (itv_.index (), it_.index ()) < (*this) ().size1 (), bad_index ());
                    return functor_type::index1 (itv_.index (), it_.index ());
                } else {
                    return i_;
                }
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                BOOST_UBLAS_CHECK (*this != (*this) ().find2 (0, i_, (*this) ().size2 ()), bad_index ());
                if (rank_ == 1) {
                    BOOST_UBLAS_CHECK (functor_type::index2 (itv_.index (), it_.index ()) < (*this) ().size2 (), bad_index ());
                    return functor_type::index2 (itv_.index (), it_.index ());
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
                BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
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
    const typename generalized_vector_of_vector<T, F, A>::value_type generalized_vector_of_vector<T, F, A>::zero_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME generalized_vector_of_vector<T, F, A>::value_type
#endif
        (0);

}}}

#endif
