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

#ifndef BOOST_UBLAS_MATRIX_EXPRESSION_H
#define BOOST_UBLAS_MATRIX_EXPRESSION_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/functional.hpp>
#include <boost/numeric/ublas/noalias.hpp>

// Expression templates based on ideas of Todd Veldhuizen and Geoffrey Furnish
// Iterators based on ideas of Jeremy Siek

namespace boost { namespace numeric { namespace ublas {

    // Base class for the Barton Nackman trick
    template<class E>
    struct matrix_expression:
        private boost::nonassignable {
        BOOST_STATIC_CONSTANT (int, complexity = 0);
        typedef E expression_type;
        typedef matrix_tag type_category;
        typedef abstract_tag simd_category;
        // FIXME: Why doesn't this work?
        // typedef typename E::size_type size_type;
        typedef std::size_t size_type;
        typedef noalias_proxy<E> noalias_proxy_type;
        typedef const matrix_row<const E> const_matrix_row_type;
        typedef matrix_row<E> matrix_row_type;
        typedef const matrix_column<const E> const_matrix_column_type;
        typedef matrix_column<E> matrix_column_type;
        typedef const matrix_range<const E> const_matrix_range_type;
        typedef matrix_range<E> matrix_range_type;
        typedef const matrix_slice<const E> const_matrix_slice_type;
        typedef matrix_slice<E> matrix_slice_type;
        typedef const matrix_indirect<const E> const_matrix_indirect_type;
        typedef matrix_indirect<E> matrix_indirect_type;

        // This class could define an common interface for all
        // statically derived expression type classes.
        // Due to a compiler deficiency - one can not reference class typedefs of E
        // on MSVC 6.0 (error C2027) - we only implement the casts.

        BOOST_UBLAS_INLINE
        const expression_type &operator () () const {
            return *static_cast<const expression_type *> (this);
        }
        BOOST_UBLAS_INLINE
        expression_type &operator () () {
            return *static_cast<expression_type *> (this);
        }

        BOOST_UBLAS_INLINE
        noalias_proxy_type noalias () {
            return noalias_proxy_type (operator () ());
        }
        BOOST_UBLAS_INLINE
        const_matrix_row_type operator [] (size_type i) const {
            return const_matrix_row_type (operator () (), i);
        }
        BOOST_UBLAS_INLINE
        matrix_row_type operator [] (size_type i) {
            return matrix_row_type (operator () (), i);
        }
        BOOST_UBLAS_INLINE
        const_matrix_row_type row (size_type i) const {
            return const_matrix_row_type (operator () (), i);
        }
        BOOST_UBLAS_INLINE
        matrix_row_type row (size_type i) {
            return matrix_row_type (operator () (), i);
        }
        BOOST_UBLAS_INLINE
        const_matrix_column_type column (size_type j) const {
            return const_matrix_column_type (operator () (), j);
        }
        BOOST_UBLAS_INLINE
        matrix_column_type column (size_type j) {
            return matrix_column_type (operator () (), j);
        }
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_INLINE
        const_matrix_range_type operator () (const range &r1, const range &r2) const {
            return const_matrix_range_type (operator () (), r1, r2);
        }
        BOOST_UBLAS_INLINE
        matrix_range_type operator () (const range &r1, const range &r2) {
            return matrix_range_type (operator () (), r1, r2);
        }
        BOOST_UBLAS_INLINE
        const_matrix_slice_type operator () (const slice &s1, const slice &s2) const {
            return const_matrix_slice_type (operator () (), s1, s2);
        }
        BOOST_UBLAS_INLINE
        matrix_slice_type operator () (const slice &s1, const slice &s2) {
            return matrix_slice_type (operator () (), s1, s2);
        }
        template<class A>
        BOOST_UBLAS_INLINE
        const_matrix_indirect_type operator () (const indirect_array<A> &ia1, const indirect_array<A> &ia2) const {
            return const_matrix_indirect_type (operator () (), ia1, ia2);
        }
        template<class A>
        BOOST_UBLAS_INLINE
        matrix_indirect_type operator () (const indirect_array<A> &ia1, const indirect_array<A> &ia2) {
            return matrix_indirect_type (operator () (), ia1, ia2);
        }
#else
        BOOST_UBLAS_INLINE
        const_matrix_range_type project (const range &r1, const range &r2) const {
            return const_matrix_range_type (operator () (), r1, r2);
        }
        BOOST_UBLAS_INLINE
        matrix_range_type project (const range &r1, const range &r2) {
            return matrix_range_type (operator () (), r1, r2);
        }
        BOOST_UBLAS_INLINE
        const_matrix_slice_type project (const slice &s1, const slice &s2) const {
            return const_matrix_slice_type (operator () (), s1, s2);
        }
        BOOST_UBLAS_INLINE
        matrix_slice_type project (const slice &s1, const slice &s2) {
            return matrix_slice_type (operator () (), s1, s2);
        }
        template<class A>
        BOOST_UBLAS_INLINE
        const_matrix_indirect_type project (const indirect_array<A> &ia1, const indirect_array<A> &ia2) const {
            return const_matrix_indirect_type (operator () (), ia1, ia2);
        }
        template<class A>
        BOOST_UBLAS_INLINE
        matrix_indirect_type project (const indirect_array<A> &ia1, const indirect_array<A> &ia2) {
            return matrix_indirect_type (operator () (), ia1, ia2);
        }
#endif
    };

#ifdef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
    struct iterator1_tag {};
    struct iterator2_tag {};

    template<class I>
    BOOST_UBLAS_INLINE
    typename I::dual_iterator_type begin (const I &it, iterator1_tag) {
        return it ().find2 (1, it.index1 (), 0);
    }
    template<class I>
    BOOST_UBLAS_INLINE
    typename I::dual_iterator_type end (const I &it, iterator1_tag) {
        return it ().find2 (1, it.index1 (), it ().size2 ());
    }
    template<class I>
    BOOST_UBLAS_INLINE
    typename I::dual_reverse_iterator_type rbegin (const I &it, iterator1_tag) {
        return typename I::dual_reverse_iterator_type (end (it, iterator1_tag ()));
    }
    template<class I>
    BOOST_UBLAS_INLINE
    typename I::dual_reverse_iterator_type rend (const I &it, iterator1_tag) {
        return typename I::dual_reverse_iterator_type (begin (it, iterator1_tag ()));
    }

    template<class I>
    BOOST_UBLAS_INLINE
    typename I::dual_iterator_type begin (const I &it, iterator2_tag) {
        return it ().find1 (1, 0, it.index2 ());
    }
    template<class I>
    BOOST_UBLAS_INLINE
    typename I::dual_iterator_type end (const I &it, iterator2_tag) {
        return it ().find1 (1, it ().size1 (), it.index2 ());
    }
    template<class I>
    BOOST_UBLAS_INLINE
    typename I::dual_reverse_iterator_type rbegin (const I &it, iterator2_tag) {
        return typename I::dual_reverse_iterator_type (end (it, iterator2_tag ()));
    }
    template<class I>
    BOOST_UBLAS_INLINE
    typename I::dual_reverse_iterator_type rend (const I &it, iterator2_tag) {
        return typename I::dual_reverse_iterator_type (begin (it, iterator2_tag ()));
    }
#endif

#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
    template<class E>
    class matrix_const_reference:
        public matrix_expression<matrix_const_reference<E> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<matrix_const_reference<E> >::operator ();
#endif
        typedef E expression_type;
        typedef typename E::simd_category simd_category;
        typedef typename E::size_type size_type;
        typedef typename E::difference_type difference_type;
        typedef typename E::value_type value_type;
        typedef typename E::const_reference const_reference;
        typedef const_reference reference;
        typedef typename E::const_pointer const_pointer;
        typedef const_pointer pointer;
        typedef typename E::orientation_category orientation_category;
        typedef typename E::const_iterator1 const_iterator1_type;
        typedef typename E::const_iterator2 const_iterator2_type;
        typedef unknown_storage_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_const_reference ():
            e_ (nil_) {}
        BOOST_UBLAS_INLINE
        matrix_const_reference (const expression_type &e):
            e_ (e) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const {
            return e_.size1 ();
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const {
            return e_.size2 ();
        }
        BOOST_UBLAS_INLINE
        const expression_type &expression () const {
            return e_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return expression () (i, j);
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const matrix_const_reference &mr) const {
            return &(*this).expression () == &mr.expression ();
        }

        typedef const_iterator1_type const_iterator1;
        typedef const_iterator1 iterator1;
        typedef const_iterator2_type const_iterator2;
        typedef const_iterator2 iterator2;

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j) const {
            return expression ().find1 (rank, i, j);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            return expression ().find2 (rank, i, j);
        }

        // Iterators are the iterators of the referenced expression.

        BOOST_UBLAS_INLINE
        const_iterator1 begin1 () const {
            return expression ().begin1 ();
        }
        BOOST_UBLAS_INLINE
        const_iterator1 end1 () const {
            return expression ().end1 ();
        }

        BOOST_UBLAS_INLINE
        const_iterator2 begin2 () const {
            return expression ().begin2 ();
        }
        BOOST_UBLAS_INLINE
        const_iterator2 end2 () const {
            return expression ().end2 ();
        }

        // Reverse iterators

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator1 rbegin1 () const {
            return const_reverse_iterator1 (end1 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator1 rend1 () const {
            return const_reverse_iterator1 (begin1 ());
        }

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> reverse_iterator1;
#else
        typedef reverse_iterator_base1<const_iterator1> reverse_iterator1;
#endif

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
#else
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> reverse_iterator2;
#else
        typedef reverse_iterator_base2<const_iterator2> reverse_iterator2;
#endif

    private:
        const expression_type &e_;
        static expression_type nil_;
    };

    template<class E>
    typename matrix_const_reference<E>::expression_type matrix_const_reference<E>::nil_;
#endif

    template<class E>
    class matrix_reference:
        public matrix_expression<matrix_reference<E> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<matrix_reference<E> >::operator ();
#endif
        typedef E expression_type;
        typedef typename E::simd_category simd_category;
        typedef typename E::size_type size_type;
        typedef typename E::difference_type difference_type;
        typedef typename E::value_type value_type;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef typename E::const_reference const_reference;
        typedef typename E::reference reference;
        typedef typename E::const_pointer const_pointer;
        typedef typename E::pointer pointer;
#else
        typedef typename E::const_reference const_reference;
        typedef typename boost::mpl::if_c<boost::is_const<E>::value,
                                          typename E::const_reference,
                                          typename E::reference>::type reference;
        typedef typename E::const_pointer const_pointer;
        typedef typename boost::mpl::if_c<boost::is_const<E>::value,
                                          typename E::const_pointer,
                                          typename E::pointer>::type pointer;
#endif
        typedef const matrix_reference<E> const_self_type;
        typedef matrix_reference<E> self_type;
        typedef const_self_type const_closure_type;
        typedef const_closure_type closure_type;
        typedef typename E::orientation_category orientation_category;
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef typename E::const_iterator1 const_iterator1_type;
        typedef typename E::iterator1 iterator1_type;
        typedef typename E::const_iterator2 const_iterator2_type;
        typedef typename E::iterator2 iterator2_type;
#else
        typedef typename E::const_iterator1 const_iterator1_type;
        typedef typename boost::mpl::if_c<boost::is_const<E>::value,
                                          typename E::const_iterator1,
                                          typename E::iterator1>::type iterator1_type;
        typedef typename E::const_iterator2 const_iterator2_type;
        typedef typename boost::mpl::if_c<boost::is_const<E>::value,
                                          typename E::const_iterator2,
                                          typename E::iterator2>::type iterator2_type;
#endif
        typedef unknown_storage_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_reference ():
              e_ (nil_) {}
        BOOST_UBLAS_INLINE
        matrix_reference (expression_type &e):
              e_ (e) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const {
            return e_.size1 ();
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const {
            return e_.size2 ();
        }
        BOOST_UBLAS_INLINE
        const expression_type &expression () const {
            return e_;
        }
        BOOST_UBLAS_INLINE
        expression_type &expression () {
            return e_;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const matrix_reference &mr) const {
            return &(*this).expression () == &mr.expression ();
        }

        // Resizing
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
        BOOST_UBLAS_INLINE
        void resize (size_type size1, size_type size2) {
            expression ().resize (size1, size2);
        }
#else
        BOOST_UBLAS_INLINE
        void resize (size_type size1, size_type size2) const {
            expression ().resize (size1, size2);
        }
#endif

        // Element access
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return expression () (i, j);
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
            return expression () (i, j);
        }
#else
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) const {
            return e_ (i, j);
        }
#endif

        // Assignment
        template<class AE>
        BOOST_UBLAS_INLINE
        matrix_reference &operator = (const matrix_expression<AE> &ae) {
            expression ().operator = (ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        matrix_reference &reset (const matrix_expression<AE> &ae) {
            expression ().reset (ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        matrix_reference &assign (const matrix_expression<AE> &ae) {
            expression ().assign (ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        matrix_reference &operator += (const matrix_expression<AE> &ae) {
            expression ().operator += (ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        matrix_reference &plus_assign (const matrix_expression<AE> &ae) {
            expression ().plus_assign (ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        matrix_reference &operator -= (const matrix_expression<AE> &ae) {
            expression ().operator -= (ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        matrix_reference &minus_assign (const matrix_expression<AE> &ae) {
            expression ().minus_assign (ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        matrix_reference &operator *= (const AT &at) {
            expression ().operator *= (at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        matrix_reference &operator /= (const AT &at) {
            expression ().operator /= (at);
            return *this;
        }

        typedef const_iterator1_type const_iterator1;
        typedef iterator1_type iterator1;
        typedef const_iterator2_type const_iterator2;
        typedef iterator2_type iterator2;

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j) const {
            return expression ().find1 (rank, i, j);
        }
        BOOST_UBLAS_INLINE
        iterator1 find1 (int rank, size_type i, size_type j) {
            return expression ().find1 (rank, i, j);
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            return expression ().find2 (rank, i, j);
        }
        BOOST_UBLAS_INLINE
        iterator2 find2 (int rank, size_type i, size_type j) {
            return expression ().find2 (rank, i, j);
        }

        // Iterators are the iterators of the referenced expression.

        BOOST_UBLAS_INLINE
        const_iterator1 begin1 () const {
            return expression ().begin1 ();
        }
        BOOST_UBLAS_INLINE
        const_iterator1 end1 () const {
            return expression ().end1 ();
        }

        BOOST_UBLAS_INLINE
        iterator1 begin1 () {
            return expression ().begin1 ();
        }
        BOOST_UBLAS_INLINE
        iterator1 end1 () {
            return expression ().end1 ();
        }

        BOOST_UBLAS_INLINE
        const_iterator2 begin2 () const {
            return expression ().begin2 ();
        }
        BOOST_UBLAS_INLINE
        const_iterator2 end2 () const {
            return expression ().end2 ();
        }

        BOOST_UBLAS_INLINE
        iterator2 begin2 () {
            return expression ().begin2 ();
        }
        BOOST_UBLAS_INLINE
        iterator2 end2 () {
            return expression ().end2 ();
        }

        // Reverse iterators

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator1 rbegin1 () const {
            return const_reverse_iterator1 (end1 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator1 rend1 () const {
            return const_reverse_iterator1 (begin1 ());
        }

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<iterator1, value_type, reference> reverse_iterator1;
#else
        typedef reverse_iterator_base1<iterator1> reverse_iterator1;
#endif

        BOOST_UBLAS_INLINE
        reverse_iterator1 rbegin1 () {
            return reverse_iterator1 (end1 ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator1 rend1 () {
            return reverse_iterator1 (begin1 ());
        }

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
#else
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base2<iterator2, value_type, reference> reverse_iterator2;
#else
        typedef reverse_iterator_base2<iterator2> reverse_iterator2;
#endif

        BOOST_UBLAS_INLINE
        reverse_iterator2 rbegin2 () {
            return reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator2 rend2 () {
            return reverse_iterator2 (begin2 ());
        }

    private:
        expression_type &e_;
        static expression_type nil_;
    };

    template<class E>
    typename matrix_reference<E>::expression_type matrix_reference<E>::nil_;

    template<class E1, class E2, class F>
    class vector_matrix_binary:
        public matrix_expression<vector_matrix_binary<E1, E2, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<vector_matrix_binary<E1, E2, F> >::operator ();
#endif
        typedef E1 expression1_type;
        typedef E2 expression2_type;
        typedef F functor_type;
        typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
        typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
        typedef typename F::result_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const value_type *const_pointer;
        typedef const_pointer pointer;
        typedef typename E1::const_closure_type expression1_closure_type;
        typedef typename E2::const_closure_type expression2_closure_type;
        typedef const vector_matrix_binary<E1, E2, F> const_self_type;
        typedef vector_matrix_binary<E1, E2, F> self_type;
        typedef const_self_type const_closure_type;
        typedef const_closure_type closure_type;
        typedef unknown_orientation_tag orientation_category;
        typedef typename E1::const_iterator const_iterator1_type;
        typedef typename E2::const_iterator const_iterator2_type;
        typedef unknown_storage_tag storage_category;

        // Construction and destruction 
        BOOST_UBLAS_INLINE
        vector_matrix_binary ():
            e1_ (), e2_ () {}
        BOOST_UBLAS_INLINE
        vector_matrix_binary (const expression1_type &e1, const expression2_type &e2): 
            e1_ (e1), e2_ (e2) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const {
            return e1_.size ();
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const { 
            return e2_.size ();
        }
        BOOST_UBLAS_INLINE
        const expression1_closure_type &expression1 () const {
            return e1_;
        }
        BOOST_UBLAS_INLINE
        const expression2_closure_type &expression2 () const {
            return e2_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return functor_type () (e1_ (i), e2_ (j));
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef typename iterator_restrict_traits<typename const_iterator1_type::iterator_category,
                                                  typename const_iterator2_type::iterator_category>::iterator_category iterator_category;
        typedef indexed_const_iterator1<const_closure_type, iterator_category> const_iterator1;
        typedef const_iterator1 iterator1;
        typedef indexed_const_iterator2<const_closure_type, iterator_category> const_iterator2;
        typedef const_iterator2 iterator2;
#else
        class const_iterator1;
        typedef const_iterator1 iterator1;
        class const_iterator2;
        typedef const_iterator2 iterator2;
#endif
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j) const {
            const_iterator1_type it1 (e1_.find (i));
            const_iterator1_type it1_end (e1_.find (size1 ()));
            const_iterator2_type it2 (e2_.find (j));
            const_iterator2_type it2_end (e2_.find (size2 ()));
            if (it2 == it2_end || (rank == 1 && (it2.index () != j || *it2 == value_type ()))) {
                it1 = it1_end;
                it2 = it2_end;
            }
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator1 (*this, it1.index (), it2.index ());
#else
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            return const_iterator1 (*this, it1, it2, it2 != it2_end ? *it2 : value_type ());
#else
            return const_iterator1 (*this, it1, it2);
#endif
#endif
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            const_iterator2_type it2 (e2_.find (j));
            const_iterator2_type it2_end (e2_.find (size2 ()));
            const_iterator1_type it1 (e1_.find (i));
            const_iterator1_type it1_end (e1_.find (size1 ()));
            if (it1 == it1_end || (rank == 1 && (it1.index () != i || *it1 == value_type ()))) {
                it2 = it2_end;
                it1 = it1_end;
            }
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator2 (*this, it1.index (), it2.index ());
#else
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            return const_iterator2 (*this, it1, it2, it1 != it1_end ? *it1 : value_type ());
#else
            return const_iterator2 (*this, it1, it2);
#endif
#endif
        }

        // Iterators enhance the iterators of the referenced expressions
        // with the binary functor.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<vector_matrix_binary>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
                                                                          typename E2::const_iterator::iterator_category>::iterator_category>::template
                iterator_base<const_iterator1, value_type>::type {
#else
            public random_access_iterator_base<typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
                                                                                 typename E2::const_iterator::iterator_category>::iterator_category,
                                               const_iterator1, value_type> {
#endif
        public:
            typedef typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
                                                      typename E2::const_iterator::iterator_category>::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename vector_matrix_binary::difference_type difference_type;
            typedef typename vector_matrix_binary::value_type value_type;
            typedef typename vector_matrix_binary::const_reference reference;
            typedef typename vector_matrix_binary::const_pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), it1_ (), it2_ (), t2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &vmb, const const_iterator1_type &it1, const const_iterator2_type &it2, value_type t2):
                container_const_reference<self_type> (vmb), it1_ (it1), it2_ (it2), t2_ (t2) {}
#else
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &vmb, const const_iterator1_type &it1, const const_iterator2_type &it2):
                container_const_reference<self_type> (vmb), it1_ (it1), it2_ (it2) {}
#endif

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
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (*it1_, t2_);
#else
                return functor_type () (*it1_, *it2_);
#endif
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
                return it1_.index ();
            }
            BOOST_UBLAS_INLINE
            size_type  index2 () const {
                return it2_.index ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
                it2_ = it.it2_;
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                t2_ = it.t2_;
#endif
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
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            const_iterator1_type it1_;
            // Mutable due to assignment
            /* const */ const_iterator2_type it2_;
            value_type t2_;
#else
            const_iterator1_type it1_;
            const_iterator2_type it2_;
#endif
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
        class const_iterator2:
            public container_const_reference<vector_matrix_binary>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
                                                                          typename E2::const_iterator::iterator_category>::iterator_category>::template
                iterator_base<const_iterator2, value_type>::type {
#else
            public random_access_iterator_base<typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
                                                                                 typename E2::const_iterator::iterator_category>::iterator_category,
                                               const_iterator2, value_type> {
#endif
        public:
            typedef typename iterator_restrict_traits<typename E1::const_iterator::iterator_category, 
                                                      typename E2::const_iterator::iterator_category>::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename vector_matrix_binary::difference_type difference_type;
            typedef typename vector_matrix_binary::value_type value_type;
            typedef typename vector_matrix_binary::const_reference reference;
            typedef typename vector_matrix_binary::const_pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), it1_ (), it2_ (), t1_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &vmb, const const_iterator1_type &it1, const const_iterator2_type &it2, value_type t1):
                container_const_reference<self_type> (vmb), it1_ (it1), it2_ (it2), t1_ (t1) {}
#else
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &vmb, const const_iterator1_type &it1, const const_iterator2_type &it2):
                container_const_reference<self_type> (vmb), it1_ (it1), it2_ (it2) {}
#endif

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
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (t1_, *it2_);
#else
                return functor_type () (*it1_, *it2_);
#endif
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
                return it1_.index ();
            }
            BOOST_UBLAS_INLINE
            size_type  index2 () const {
                return it2_.index ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
                it2_ = it.it2_;
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                t1_ = it.t1_;
#endif
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
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            // Mutable due to assignment
            /* const */ const_iterator1_type it1_;
            const_iterator2_type it2_;
            value_type t1_;
#else
            const_iterator1_type it1_;
            const_iterator2_type it2_;
#endif
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
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

    private:
        expression1_closure_type e1_;
        expression2_closure_type e2_;
    };

    template<class E1, class E2, class F>
    struct vector_matrix_binary_traits {
        typedef vector_matrix_binary<E1, E2, F> expression_type;
#ifdef BOOST_UBLAS_USE_ET
        typedef expression_type result_type; 
#else
        typedef matrix<typename F::result_type> result_type;
#endif
    };

    // (outer_prod (v1, v2)) [i] [j] = v1 [i] * v2 [j]
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename vector_matrix_binary_traits<E1, E2, scalar_multiplies<typename E1::value_type, typename E2::value_type> >::result_type
    outer_prod (const vector_expression<E1> &e1,
                const vector_expression<E2> &e2) {
        BOOST_STATIC_ASSERT (E1::complexity == 0 && E2::complexity == 0);
        typedef BOOST_UBLAS_TYPENAME vector_matrix_binary_traits<E1, E2, scalar_multiplies<BOOST_UBLAS_TYPENAME E1::value_type, BOOST_UBLAS_TYPENAME E2::value_type> >::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    template<class E, class F>
    class matrix_unary1:
        public matrix_expression<matrix_unary1<E, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<matrix_unary1<E, F> >::operator ();
#endif
        typedef E expression_type;
        typedef F functor_type;
        typedef typename E::size_type size_type;
        typedef typename E::difference_type difference_type;
        typedef typename F::result_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const value_type *const_pointer;
        typedef const_pointer pointer;
        typedef typename E::const_closure_type expression_closure_type;
        typedef const matrix_unary1<E, F> const_self_type;
        typedef matrix_unary1<E, F> self_type;
        typedef const_self_type const_closure_type;
        typedef const_closure_type closure_type;
        typedef typename E::orientation_category orientation_category;
        typedef typename E::const_iterator1 const_iterator1_type;
        typedef typename E::const_iterator2 const_iterator2_type;
        typedef unknown_storage_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_unary1 ():
            e_ () {}
        BOOST_UBLAS_INLINE
        matrix_unary1 (const expression_type &e):
            e_ (e) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const {
            return e_.size1 ();
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const {
            return e_.size2 ();
        }
        BOOST_UBLAS_INLINE
        const expression_closure_type &expression () const {
            return e_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return functor_type () (e_ (i, j));
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator1<const_closure_type, typename const_iterator1_type::iterator_category> const_iterator1;
        typedef const_iterator1 iterator1;
        typedef indexed_const_iterator2<const_closure_type, typename const_iterator2_type::iterator_category> const_iterator2;
        typedef const_iterator2 iterator2;
#else
        class const_iterator1;
        typedef const_iterator1 iterator1;
        class const_iterator2;
        typedef const_iterator2 iterator2;
#endif
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j) const {
            const_iterator1_type it1 (e_.find1 (rank, i, j));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator1 (*this, it1.index1 (), it1.index2 ());
#else
            return const_iterator1 (*this, it1);
#endif
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            const_iterator2_type it2 (e_.find2 (rank, i, j));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator2 (*this, it2.index1 (), it2.index2 ());
#else
            return const_iterator2 (*this, it2);
#endif
        }

        // Iterators enhance the iterators of the referenced expression
        // with the unary functor.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<matrix_unary1>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename E::const_iterator1::iterator_category>::template
                iterator_base<const_iterator1, value_type>::type {
#else
            public random_access_iterator_base<typename E::const_iterator1::iterator_category,
                                               const_iterator1, value_type> {
#endif
        public:
            typedef typename E::const_iterator1::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_unary1::difference_type difference_type;
            typedef typename matrix_unary1::value_type value_type;
            typedef typename matrix_unary1::const_reference reference;
            typedef typename matrix_unary1::const_pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &mu, const const_iterator1_type &it):
                container_const_reference<self_type> (mu), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                ++ it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                -- it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator += (difference_type n) {
                it_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -= (difference_type n) {
                it_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return functor_type () (*it_);
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
                return it_.index1 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it_.index2 ();
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator1_type it_;
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
        class const_iterator2:
            public container_const_reference<matrix_unary1>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename E::const_iterator2::iterator_category>::template
                iterator_base<const_iterator2, value_type>::type {
#else
            public random_access_iterator_base<typename E::const_iterator2::iterator_category,
                                               const_iterator2, value_type> {
#endif                                               
        public:
            typedef typename E::const_iterator2::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_unary1::difference_type difference_type;
            typedef typename matrix_unary1::value_type value_type;
            typedef typename matrix_unary1::const_reference reference;
            typedef typename matrix_unary1::const_pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &mu, const const_iterator2_type &it):
                container_const_reference<self_type> (mu), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                ++ it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                -- it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator += (difference_type n) {
                it_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -= (difference_type n) {
                it_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return functor_type () (*it_);
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
                return it_.index1 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it_.index2 ();
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator2_type it_;
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
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

    private:
        expression_closure_type e_;
    };

    template<class E, class F>
    struct matrix_unary1_traits {
        typedef matrix_unary1<E, F> expression_type;
#ifdef BOOST_UBLAS_USE_ET
        typedef expression_type result_type; 
#else
        typedef matrix<typename F::result_type> result_type;
#endif
    };

    // (- m) [i] [j] = - m [i] [j]
    template<class E>
    BOOST_UBLAS_INLINE
    typename matrix_unary1_traits<E, scalar_negate<typename E::value_type> >::result_type
    operator - (const matrix_expression<E> &e) {
        typedef BOOST_UBLAS_TYPENAME matrix_unary1_traits<E, scalar_negate<BOOST_UBLAS_TYPENAME E::value_type> >::expression_type expression_type;
        return expression_type (e ());
    }

    // (conj m) [i] [j] = conj (m [i] [j])
    template<class E> 
    BOOST_UBLAS_INLINE
    typename matrix_unary1_traits<E, scalar_conj<typename E::value_type> >::result_type
    conj (const matrix_expression<E> &e) {
        typedef BOOST_UBLAS_TYPENAME matrix_unary1_traits<E, scalar_conj<BOOST_UBLAS_TYPENAME E::value_type> >::expression_type expression_type;
        return expression_type (e ());
    }

    // (real m) [i] [j] = real (m [i] [j])
    template<class E> 
    BOOST_UBLAS_INLINE
    typename matrix_unary1_traits<E, scalar_real<typename E::value_type> >::result_type
    real (const matrix_expression<E> &e) {
        typedef BOOST_UBLAS_TYPENAME matrix_unary1_traits<E, scalar_real<BOOST_UBLAS_TYPENAME E::value_type> >::expression_type expression_type;
        return expression_type (e ());
    }

    // (imag m) [i] [j] = imag (m [i] [j])
    template<class E> 
    BOOST_UBLAS_INLINE
    typename matrix_unary1_traits<E, scalar_imag<typename E::value_type> >::result_type
    imag (const matrix_expression<E> &e) {
        typedef BOOST_UBLAS_TYPENAME matrix_unary1_traits<E, scalar_imag<BOOST_UBLAS_TYPENAME E::value_type> >::expression_type expression_type;
        return expression_type (e ());
    }

    template<class E, class F>
    class matrix_unary2:
        public matrix_expression<matrix_unary2<E, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<matrix_unary2<E, F> >::operator ();
#endif
        typedef E expression_type;
        typedef F functor_type;
        typedef typename E::size_type size_type;
        typedef typename E::difference_type difference_type;
        typedef typename F::result_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const value_type *const_pointer;
        typedef const_pointer pointer;
        typedef typename E::const_closure_type expression_closure_type;
        typedef const matrix_unary2<E, F> const_self_type;
        typedef matrix_unary2<E, F> self_type;
        typedef const_self_type const_closure_type;
        typedef const_closure_type closure_type;
        typedef typename E::orientation_category orientation_category;
        typedef typename E::const_iterator1 const_iterator2_type;
        typedef typename E::const_iterator2 const_iterator1_type;
        typedef unknown_storage_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_unary2 (): 
            e_ () {}
        BOOST_UBLAS_INLINE
        matrix_unary2 (const expression_type &e): 
            e_ (e) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const {
            return e_.size2 ();
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const {
            return e_.size1 ();
        }
        BOOST_UBLAS_INLINE
        const expression_closure_type &expression () const {
            return e_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return functor_type () (e_ (j, i)); 
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator1<const_closure_type, typename const_iterator1_type::iterator_category> const_iterator1;
        typedef const_iterator1 iterator1;
        typedef indexed_const_iterator2<const_closure_type, typename const_iterator2_type::iterator_category> const_iterator2;
        typedef const_iterator2 iterator2;
#else
        class const_iterator1;
        typedef const_iterator1 iterator1;
        class const_iterator2;
        typedef const_iterator2 iterator2;
#endif
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j) const {
            const_iterator1_type it1 (e_.find2 (rank, j, i));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator1 (*this, it1.index2 (), it1.index1 ());
#else
            return const_iterator1 (*this, it1);
#endif
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            const_iterator2_type it2 (e_.find1 (rank, j, i));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator2 (*this, it2.index2 (), it2.index1 ());
#else
            return const_iterator2 (*this, it2);
#endif
        }

        // Iterators enhance the iterators of the referenced expression
        // with the unary functor.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<matrix_unary2>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename E::const_iterator2::iterator_category>::template
                iterator_base<const_iterator1, value_type>::type {
#else
            public random_access_iterator_base<typename E::const_iterator2::iterator_category,
                                               const_iterator1, value_type> {
#endif
        public:
            typedef typename E::const_iterator2::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_unary2::difference_type difference_type;
            typedef typename matrix_unary2::value_type value_type;
            typedef typename matrix_unary2::const_reference reference;
            typedef typename matrix_unary2::const_pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &mu, const const_iterator1_type &it):
                container_const_reference<self_type> (mu), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                ++ it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                -- it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator += (difference_type n) {
                it_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -= (difference_type n) {
                it_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return functor_type () (*it_);
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
                return it_.index2 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it_.index1 ();
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator1_type it_;
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
        class const_iterator2:
            public container_const_reference<matrix_unary2>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename E::const_iterator1::iterator_category>::template
                iterator_base<const_iterator2, value_type>::type {
#else
            public random_access_iterator_base<typename E::const_iterator1::iterator_category,
                                               const_iterator2, value_type> {
#endif
        public:
            typedef typename E::const_iterator1::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_unary2::difference_type difference_type;
            typedef typename matrix_unary2::value_type value_type;
            typedef typename matrix_unary2::const_reference reference;
            typedef typename matrix_unary2::const_pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &mu, const const_iterator2_type &it):
                container_const_reference<self_type> (mu), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                ++ it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                -- it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator += (difference_type n) {
                it_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -= (difference_type n) {
                it_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return functor_type () (*it_);
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
                return it_.index2 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it_.index1 ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_iterator2_type it_;
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
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

    private:
        expression_closure_type e_;
    };

    template<class E, class F>
    struct matrix_unary2_traits {
        typedef matrix_unary2<E, F> expression_type;
#ifdef BOOST_UBLAS_USE_ET
        typedef expression_type result_type; 
#else
        typedef matrix<typename F::result_type> result_type;
#endif
    };

    // (trans m) [i] [j] = m [j] [i]
    template<class E>
    BOOST_UBLAS_INLINE
    typename matrix_unary2_traits<E, scalar_identity<typename E::value_type> >::result_type
    trans (const matrix_expression<E> &e) {
        typedef BOOST_UBLAS_TYPENAME matrix_unary2_traits<E, scalar_identity<BOOST_UBLAS_TYPENAME E::value_type> >::expression_type expression_type;
        return expression_type (e ());
    }

    // (herm m) [i] [j] = conj (m [j] [i])
    template<class E> 
    BOOST_UBLAS_INLINE
    typename matrix_unary2_traits<E, scalar_conj<typename E::value_type> >::result_type
    herm (const matrix_expression<E> &e) {
        typedef BOOST_UBLAS_TYPENAME matrix_unary2_traits<E, scalar_conj<BOOST_UBLAS_TYPENAME E::value_type> >::expression_type expression_type;
        return expression_type (e ());
    }

    template<class E1, class E2, class F>
    class matrix_binary:
        public matrix_expression<matrix_binary<E1, E2, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<matrix_binary<E1, E2, F> >::operator ();
#endif
        typedef E1 expression1_type;
        typedef E2 expression2_type;
        typedef F functor_type;
        typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
        typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
        typedef typename F::result_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const value_type *const_pointer;
        typedef const_pointer pointer;
        typedef typename E1::const_closure_type expression1_closure_type;
        typedef typename E2::const_closure_type expression2_closure_type;
        typedef const matrix_binary<E1, E2, F> const_self_type;
        typedef matrix_binary<E1, E2, F> self_type;
        typedef const_self_type const_closure_type;
        typedef const_closure_type closure_type;
        typedef unknown_orientation_tag orientation_category;
        typedef typename E1::const_iterator1 const_iterator11_type;
        typedef typename E1::const_iterator2 const_iterator12_type;
        typedef typename E2::const_iterator1 const_iterator21_type;
        typedef typename E2::const_iterator2 const_iterator22_type;
        typedef unknown_storage_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_binary (): 
            e1_ (), e2_ () {}
        BOOST_UBLAS_INLINE
        matrix_binary (const E1 &e1, const E2 &e2): 
            e1_ (e1), e2_ (e2) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const { 
            return BOOST_UBLAS_SAME (e1_.size1 (), e2_.size1 ());
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const { 
            return BOOST_UBLAS_SAME (e1_.size2 (), e2_.size2 ());
        }
        BOOST_UBLAS_INLINE
        const expression1_closure_type &expression1 () const {
            return e1_;
        }
        BOOST_UBLAS_INLINE
        const expression2_closure_type &expression2 () const {
            return e2_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const { 
            return functor_type () (e1_ (i, j), e2_ (i, j)); 
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef typename iterator_restrict_traits<typename const_iterator11_type::iterator_category,
                                                  typename const_iterator21_type::iterator_category>::iterator_category iterator_category1;
        typedef indexed_const_iterator1<const_closure_type, iterator_category1> const_iterator1;
        typedef const_iterator1 iterator1;
        typedef typename iterator_restrict_traits<typename const_iterator12_type::iterator_category,
                                                  typename const_iterator22_type::iterator_category>::iterator_category iterator_category2;
        typedef indexed_const_iterator2<const_closure_type, iterator_category2> const_iterator2;
        typedef const_iterator2 iterator2;
#else
        class const_iterator1;
        typedef const_iterator1 iterator1;
        class const_iterator2;
        typedef const_iterator2 iterator2;
#endif
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j) const {
            const_iterator11_type it11 (e1_.find1 (rank, i, j));
            const_iterator11_type it11_end (e1_.find1 (rank, size1 (), j));
            const_iterator21_type it21 (e2_.find1 (rank, i, j));
            const_iterator21_type it21_end (e2_.find1 (rank, size1 (), j));
            BOOST_UBLAS_CHECK (rank == 0 || it11 == it11_end || it11.index2 () == j, internal_logic ())
            BOOST_UBLAS_CHECK (rank == 0 || it21 == it21_end || it21.index2 () == j, internal_logic ())
            i = std::min (it11 != it11_end ? it11.index1 () : size1 (),
                          it21 != it21_end ? it21.index1 () : size1 ());
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator1 (*this, i, j);
#else
            return const_iterator1 (*this, i, j, it11, it11_end, it21, it21_end);
#endif
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            const_iterator12_type it12 (e1_.find2 (rank, i, j));
            const_iterator12_type it12_end (e1_.find2 (rank, i, size2 ()));
            const_iterator22_type it22 (e2_.find2 (rank, i, j));
            const_iterator22_type it22_end (e2_.find2 (rank, i, size2 ()));
            BOOST_UBLAS_CHECK (rank == 0 || it12 == it12_end || it12.index1 () == i, internal_logic ())
            BOOST_UBLAS_CHECK (rank == 0 || it22 == it22_end || it22.index1 () == i, internal_logic ())
            j = std::min (it12 != it12_end ? it12.index2 () : size2 (),
                          it22 != it22_end ? it22.index2 () : size2 ());
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator2 (*this, i, j);
#else
            return const_iterator2 (*this, i, j, it12, it12_end, it22, it22_end);
#endif
        }

        // Iterators enhance the iterators of the referenced expression
        // with the binary functor.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<matrix_binary>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                                          typename E2::const_iterator1::iterator_category>::iterator_category>::template
                iterator_base<const_iterator1, value_type>::type {
#else
            public random_access_iterator_base<typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                                                 typename E2::const_iterator1::iterator_category>::iterator_category,
                                               const_iterator1, value_type> {
#endif
        public:
            typedef typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                      typename E2::const_iterator1::iterator_category>::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_binary::difference_type difference_type;
            typedef typename matrix_binary::value_type value_type;
            typedef typename matrix_binary::const_reference reference;
            typedef typename matrix_binary::const_pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), i_ (), j_ (), it1_ (), it1_end_ (), it2_ (), it2_end_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &mb, size_type i, size_type j,
                             const const_iterator11_type &it1, const const_iterator11_type &it1_end,
                             const const_iterator21_type &it2, const const_iterator21_type &it2_end):
                container_const_reference<self_type> (mb), i_ (i), j_ (j), it1_ (it1), it1_end_ (it1_end), it2_ (it2), it2_end_ (it2_end) {}

            // Dense specializations
            BOOST_UBLAS_INLINE
            void increment (dense_random_access_iterator_tag) {
                ++ i_, ++ it1_, ++ it2_;
            }
            BOOST_UBLAS_INLINE
            void decrement (dense_random_access_iterator_tag) {
                -- i_, -- it1_, -- it2_;
            }
            BOOST_UBLAS_INLINE
            value_type dereference (dense_random_access_iterator_tag) const {
                return functor_type () (*it1_, *it2_);
            }

            // Packed specializations
            BOOST_UBLAS_INLINE
            void increment (packed_random_access_iterator_tag) {
                if (it1_ != it1_end_)
                    if (it1_.index1 () <= i_)
                        ++ it1_;
                if (it2_ != it2_end_)
                    if (it2_.index1 () <= i_)
                        ++ it2_;
                ++ i_;
            }
            BOOST_UBLAS_INLINE
            void decrement (packed_random_access_iterator_tag) {
                if (it1_ != it1_end_)
                    if (i_ <= it1_.index1 ())
                        -- it1_;
                if (it2_ != it2_end_)
                    if (i_ <= it2_.index1 ())
                        -- it2_;
                -- i_;
            }
            BOOST_UBLAS_INLINE
            value_type dereference (packed_random_access_iterator_tag) const {
                value_type t1 = value_type ();
                if (it1_ != it1_end_) {
                    BOOST_UBLAS_CHECK (it1_.index2 () == j_, internal_logic ());
                    if (it1_.index1 () == i_)
                        t1 = *it1_;
                }
                value_type t2 = value_type ();
                if (it2_ != it2_end_) {
                    BOOST_UBLAS_CHECK (it2_.index2 () == j_, internal_logic ());
                    if (it2_.index1 () == i_)
                        t2 = *it2_;
                }
                return functor_type () (t1, t2);
            }

            // Sparse specializations
            BOOST_UBLAS_INLINE
            void increment (sparse_bidirectional_iterator_tag) {
                size_type index1 = (*this) ().size1 ();
                if (it1_ != it1_end_) {
                    if (it1_.index1 () <= i_)
                        ++ it1_;
                    if (it1_ != it1_end_)
                        index1 = it1_.index1 ();
                }
                size_type index2 = (*this) ().size1 ();
                if (it2_ != it2_end_)
                    if (it2_.index1 () <= i_)
                        ++ it2_;
                    if (it2_ != it2_end_) {
                        index2 = it2_.index1 ();
                }
                i_ = std::min (index1, index2);
            }
            BOOST_UBLAS_INLINE
            void decrement (sparse_bidirectional_iterator_tag) {
                size_type index1 = (*this) ().size1 ();
                if (it1_ != it1_end_) {
                    if (i_ <= it1_.index1 ())
                        -- it1_;
                    if (it1_ != it1_end_)
                        index1 = it1_.index1 ();
                }
                size_type index2 = (*this) ().size1 ();
                if (it2_ != it2_end_) {
                    if (i_ <= it2_.index1 ())
                        -- it2_;
                    if (it2_ != it2_end_)
                        index2 = it2_.index1 ();
                }
                i_ = std::max (index1, index2);
            }
            BOOST_UBLAS_INLINE
            value_type dereference (sparse_bidirectional_iterator_tag) const {
                value_type t1 = value_type ();
                if (it1_ != it1_end_) {
                    BOOST_UBLAS_CHECK (it1_.index2 () == j_, internal_logic ());
                    if (it1_.index1 () == i_)
                        t1 = *it1_;
                }
                value_type t2 = value_type ();
                if (it2_ != it2_end_) {
                    BOOST_UBLAS_CHECK (it2_.index2 () == j_, internal_logic ());
                    if (it2_.index1 () == i_)
                        t2 = *it2_;
                }
                return functor_type () (t1, t2);
            }

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                increment (iterator_category ());
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                decrement (iterator_category ());
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator += (difference_type n) {
                i_ += n, it1_ += n, it2_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -= (difference_type n) {
                i_ -= n, it1_ -= n, it2_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (index2 () == it.index2 (), external_logic ());
                return index1 () - it.index1 ();
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return dereference (iterator_category ());
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
                return i_;
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                // if (it1_ != it1_end_ && it2_ != it2_end_)
                //    return BOOST_UBLAS_SAME (it1_.index2 (), it2_.index2 ());
                // else
                    return j_;
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                i_ = it.i_;
                j_ = it.j_;
                it1_ = it.it1_;
                it1_end_ = it.it1_end_;
                it2_ = it.it2_;
                it2_end_ = it.it2_end_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (index2 () == it.index2 (), external_logic ());
                return index1 () == it.index1 ();
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (index2 () == it.index2 (), external_logic ());
                return index1 () < it.index1 ();
            }

        private:
            size_type i_;
            size_type j_;
            const_iterator11_type it1_;
            const_iterator11_type it1_end_;
            const_iterator21_type it2_;
            const_iterator21_type it2_end_;
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
        class const_iterator2:
            public container_const_reference<matrix_binary>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename iterator_restrict_traits<typename E1::const_iterator2::iterator_category,
                                                                          typename E2::const_iterator2::iterator_category>::iterator_category>::template
                iterator_base<const_iterator2, value_type>::type {
#else
            public random_access_iterator_base<typename iterator_restrict_traits<typename E1::const_iterator2::iterator_category,
                                                                                 typename E2::const_iterator2::iterator_category>::iterator_category,
                                               const_iterator2, value_type> {
#endif
        public:
            typedef typename iterator_restrict_traits<typename E1::const_iterator2::iterator_category,
                                                      typename E2::const_iterator2::iterator_category>::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_binary::difference_type difference_type;
            typedef typename matrix_binary::value_type value_type;
            typedef typename matrix_binary::const_reference reference;
            typedef typename matrix_binary::const_pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), i_ (), j_ (), it1_ (), it1_end_ (), it2_ (), it2_end_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &mb, size_type i, size_type j,
                             const const_iterator12_type &it1, const const_iterator12_type &it1_end,
                             const const_iterator22_type &it2, const const_iterator22_type &it2_end):
                container_const_reference<self_type> (mb), i_ (i), j_ (j), it1_ (it1), it1_end_ (it1_end), it2_ (it2), it2_end_ (it2_end) {}

            // Dense access specializations
            BOOST_UBLAS_INLINE
            void increment (dense_random_access_iterator_tag) {
                ++ j_, ++ it1_, ++ it2_;
            }
            BOOST_UBLAS_INLINE
            void decrement (dense_random_access_iterator_tag) {
                -- j_, -- it1_, -- it2_;
            }
            BOOST_UBLAS_INLINE
            value_type dereference (dense_random_access_iterator_tag) const {
                return functor_type () (*it1_, *it2_);
            }

            // Packed specializations
            BOOST_UBLAS_INLINE
            void increment (packed_random_access_iterator_tag) {
                if (it1_ != it1_end_)
                    if (it1_.index2 () <= j_)
                        ++ it1_;
                if (it2_ != it2_end_)
                    if (it2_.index2 () <= j_)
                        ++ it2_;
                ++ j_;
            }
            BOOST_UBLAS_INLINE
            void decrement (packed_random_access_iterator_tag) {
                if (it1_ != it1_end_)
                    if (j_ <= it1_.index2 ())
                        -- it1_;
                if (it2_ != it2_end_)
                    if (j_ <= it2_.index2 ())
                        -- it2_;
                -- j_;
            }
            BOOST_UBLAS_INLINE
            value_type dereference (packed_random_access_iterator_tag) const {
                value_type t1 = value_type ();
                if (it1_ != it1_end_) {
                    BOOST_UBLAS_CHECK (it1_.index1 () == i_, internal_logic ());
                    if (it1_.index2 () == j_)
                        t1 = *it1_;
                }
                value_type t2 = value_type ();
                if (it2_ != it2_end_) {
                    BOOST_UBLAS_CHECK (it2_.index1 () == i_, internal_logic ());
                    if (it2_.index2 () == j_)
                        t2 = *it2_;
                }
                return functor_type () (t1, t2);
            }

            // Sparse specializations
            BOOST_UBLAS_INLINE
            void increment (sparse_bidirectional_iterator_tag) {
                size_type index1 = (*this) ().size2 ();
                if (it1_ != it1_end_) {
                    if (it1_.index2 () <= j_)
                        ++ it1_;
                    if (it1_ != it1_end_)
                        index1 = it1_.index2 ();
                }
                size_type index2 = (*this) ().size2 ();
                if (it2_ != it2_end_) {
                    if (it2_.index2 () <= j_)
                        ++ it2_;
                    if (it2_ != it2_end_)
                        index2 = it2_.index2 ();
                }
                j_ = std::min (index1, index2);
            }
            BOOST_UBLAS_INLINE
            void decrement (sparse_bidirectional_iterator_tag) {
                size_type index1 = (*this) ().size2 ();
                if (it1_ != it1_end_) {
                    if (j_ <= it1_.index2 ())
                        -- it1_;
                    if (it1_ != it1_end_)
                        index1 = it1_.index2 ();
                }
                size_type index2 = (*this) ().size2 ();
                if (it2_ != it2_end_) {
                    if (j_ <= it2_.index2 ())
                        -- it2_;
                    if (it2_ != it2_end_)
                        index2 = it2_.index2 ();
                }
                j_ = std::max (index1, index2);
            }
            BOOST_UBLAS_INLINE
            value_type dereference (sparse_bidirectional_iterator_tag) const {
                value_type t1 = value_type ();
                if (it1_ != it1_end_) {
                    BOOST_UBLAS_CHECK (it1_.index1 () == i_, internal_logic ());
                    if (it1_.index2 () == j_)
                        t1 = *it1_;
                }
                value_type t2 = value_type ();
                if (it2_ != it2_end_) {
                    BOOST_UBLAS_CHECK (it2_.index1 () == i_, internal_logic ());
                    if (it2_.index2 () == j_)
                        t2 = *it2_;
                }
                return functor_type () (t1, t2);
            }

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                increment (iterator_category ());
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                decrement (iterator_category ());
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator += (difference_type n) {
                j_ += n, it1_ += n, it2_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -= (difference_type n) {
                j_ -= n, it1_ -= n, it2_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (index1 () == it.index1 (), external_logic ());
                return index2 () - it.index2 ();
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return dereference (iterator_category ());
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
                // if (it1_ != it1_end_ && it2_ != it2_end_)
                //    return BOOST_UBLAS_SAME (it1_.index1 (), it2_.index1 ());
                // else
                    return i_;
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return j_;
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                i_ = it.i_;
                j_ = it.j_;
                it1_ = it.it1_;
                it1_end_ = it.it1_end_;
                it2_ = it.it2_;
                it2_end_ = it.it2_end_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (index1 () == it.index1 (), external_logic ());
                return index2 () == it.index2 ();
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                BOOST_UBLAS_CHECK (index1 () == it.index1 (), external_logic ());
                return index2 () < it.index2 ();
            }

        private:
            size_type i_;
            size_type j_;
            const_iterator12_type it1_;
            const_iterator12_type it1_end_;
            const_iterator22_type it2_;
            const_iterator22_type it2_end_;
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
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

    private:
        expression1_closure_type e1_;
        expression2_closure_type e2_;
    };

    template<class E1, class E2, class F>
    struct matrix_binary_traits {
        typedef matrix_binary<E1, E2, F> expression_type;
#ifdef BOOST_UBLAS_USE_ET
        typedef expression_type result_type; 
#else
        typedef matrix<typename F::result_type> result_type;
#endif
    };

    // (m1 + m2) [i] [j] = m1 [i] [j] + m2 [i] [j]
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_binary_traits<E1, E2, scalar_plus<typename E1::value_type,
                                                      typename E2::value_type> >::result_type
    operator + (const matrix_expression<E1> &e1,
                const matrix_expression<E2> &e2) {
        typedef BOOST_UBLAS_TYPENAME matrix_binary_traits<E1, E2, scalar_plus<BOOST_UBLAS_TYPENAME E1::value_type,
                                                                              BOOST_UBLAS_TYPENAME E2::value_type> >::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    // (m1 - m2) [i] [j] = m1 [i] [j] - m2 [i] [j]
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_binary_traits<E1, E2, scalar_minus<typename E1::value_type,
                                                       typename E2::value_type> >::result_type
    operator - (const matrix_expression<E1> &e1,
                const matrix_expression<E2> &e2) {
        typedef BOOST_UBLAS_TYPENAME matrix_binary_traits<E1, E2, scalar_minus<BOOST_UBLAS_TYPENAME E1::value_type,
                                                                               BOOST_UBLAS_TYPENAME E2::value_type> >::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    // (m1 * m2) [i] [j] = m1 [i] [j] * m2 [i] [j]
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_binary_traits<E1, E2, scalar_multiplies<typename E1::value_type,
                                                            typename E2::value_type> >::result_type
    element_prod (const matrix_expression<E1> &e1,
                  const matrix_expression<E2> &e2) {
        typedef BOOST_UBLAS_TYPENAME matrix_binary_traits<E1, E2, scalar_multiplies<BOOST_UBLAS_TYPENAME E1::value_type,
                                                                                    BOOST_UBLAS_TYPENAME E2::value_type> >::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    // (m1 / m2) [i] [j] = m1 [i] [j] / m2 [i] [j]
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_binary_traits<E1, E2, scalar_divides<typename E1::value_type,
                                                       typename E2::value_type> >::result_type
    element_div (const matrix_expression<E1> &e1,
                 const matrix_expression<E2> &e2) {
        typedef BOOST_UBLAS_TYPENAME matrix_binary_traits<E1, E2, scalar_divides<BOOST_UBLAS_TYPENAME E1::value_type,
                                                                                 BOOST_UBLAS_TYPENAME E2::value_type> >::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    template<class E1, class E2, class F>
    class matrix_binary_scalar1:
        public matrix_expression<matrix_binary_scalar1<E1, E2, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<matrix_binary_scalar1<E1, E2, F> >::operator ();
#endif
        typedef E1 expression1_type;
        typedef E2 expression2_type;
        typedef F functor_type;
        typedef typename E2::size_type size_type;
        typedef typename E2::difference_type difference_type;
        typedef typename F::result_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const value_type *const_pointer;
        typedef const_pointer pointer;
        typedef E1 expression1_closure_type;
        typedef typename E2::const_closure_type expression2_closure_type;
        typedef const matrix_binary_scalar1<E1, E2, F> const_self_type;
        typedef matrix_binary_scalar1<E1, E2, F> self_type;
        typedef const_self_type const_closure_type;
        typedef const_closure_type closure_type;
        typedef typename E2::orientation_category orientation_category;
        typedef typename E1::value_type const_iterator1_type;
        typedef typename E2::const_iterator1 const_iterator21_type;
        typedef typename E2::const_iterator2 const_iterator22_type;
        typedef unknown_storage_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_binary_scalar1 ():
            e1_ (), e2_ () {}
        BOOST_UBLAS_INLINE
        matrix_binary_scalar1 (const expression1_type &e1, const expression2_type &e2):
            e1_ (e1), e2_ (e2) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const {
            return e2_.size1 ();
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const {
            return e2_.size2 ();
        }
        BOOST_UBLAS_INLINE
        const expression1_closure_type &expression1 () const {
            return e1_;
        }
        BOOST_UBLAS_INLINE
        const expression2_closure_type &expression2 () const {
            return e2_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return functor_type () (e1_, e2_ (i, j));
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator1<const_closure_type, typename const_iterator21_type::iterator_category> const_iterator1;
        typedef const_iterator1 iterator1;
        typedef indexed_const_iterator2<const_closure_type, typename const_iterator22_type::iterator_category> const_iterator2;
        typedef const_iterator2 iterator2;
#else
        class const_iterator1;
        typedef const_iterator1 iterator1;
        class const_iterator2;
        typedef const_iterator2 iterator2;
#endif
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j) const {
            const_iterator21_type it21 (e2_.find1 (rank, i, j));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator1 (*this, it21.index1 (), it21.index2 ());
#else
            return const_iterator1 (*this, e1_, it21);
#endif
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            const_iterator22_type it22 (e2_.find2 (rank, i, j));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator2 (*this, it22.index1 (), it22.index2 ());
#else
            return const_iterator2 (*this, e1_, it22);
#endif
        }

        // Iterators enhance the iterators of the referenced expression
        // with the binary functor.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<matrix_binary_scalar1>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename E2::const_iterator1::iterator_category>::template
                iterator_base<const_iterator1, value_type>::type {
#else
            public random_access_iterator_base<typename E2::const_iterator1::iterator_category,
                                               const_iterator1, value_type> {
#endif                                               
        public:
            typedef typename E2::const_iterator1::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_binary_scalar1::difference_type difference_type;
            typedef typename matrix_binary_scalar1::value_type value_type;
            typedef typename matrix_binary_scalar1::const_reference reference;
            typedef typename matrix_binary_scalar1::const_pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &mbs, const const_iterator1_type &it1, const const_iterator21_type &it2):
                container_const_reference<self_type> (mbs), it1_ (it1), it2_ (it2) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                ++ it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                -- it2_ ;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator += (difference_type n) {
                it2_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -= (difference_type n) {
                it2_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ - it.it2_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return functor_type () (it1_, *it2_);
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
                return it2_.index1 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it2_.index2 ();
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
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ == it.it2_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ < it.it2_;
            }

        private:
            const_iterator1_type it1_;
            const_iterator21_type it2_;
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
        class const_iterator2:
            public container_const_reference<matrix_binary_scalar1>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename E2::const_iterator2::iterator_category>::template
                iterator_base<const_iterator2, value_type>::type {
#else
            public random_access_iterator_base<typename E2::const_iterator2::iterator_category,
                                               const_iterator2, value_type> {
#endif                                               
        public:
            typedef typename E2::const_iterator2::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_binary_scalar1::difference_type difference_type;
            typedef typename matrix_binary_scalar1::value_type value_type;
            typedef typename matrix_binary_scalar1::const_reference reference;
            typedef typename matrix_binary_scalar1::const_pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &mbs, const const_iterator1_type &it1, const const_iterator22_type &it2):
                container_const_reference<self_type> (mbs), it1_ (it1), it2_ (it2) {}

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
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ - it.it2_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return functor_type () (it1_, *it2_);
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
                it1_ = it.it1_;
                it2_ = it.it2_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ == it.it2_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
                return it2_ < it.it2_;
            }

        private:
            const_iterator1_type it1_;
            const_iterator22_type it2_;
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
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

    private:
        expression1_closure_type e1_;
        expression2_closure_type e2_;
    };

    template<class T1, class E2, class F>
    struct matrix_binary_scalar1_traits {
        typedef matrix_binary_scalar1<scalar_const_reference<T1>, E2, F> expression_type;
#ifdef BOOST_UBLAS_USE_ET
        typedef expression_type result_type;
#else
        typedef matrix<typename F::result_type> result_type;
#endif
    };

    // (t * m) [i] [j] = t * m [i] [j]
    template<class T1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_binary_scalar1_traits<T1, E2, scalar_multiplies<T1, typename E2::value_type> >::result_type
    operator * (const T1 &e1,
                const matrix_expression<E2> &e2) {
        typedef BOOST_UBLAS_TYPENAME matrix_binary_scalar1_traits<T1, E2, scalar_multiplies<T1, BOOST_UBLAS_TYPENAME E2::value_type> >::expression_type expression_type;
        return expression_type (e1, e2 ());
    }

    template<class E1, class E2, class F>
    class matrix_binary_scalar2:
        public matrix_expression<matrix_binary_scalar2<E1, E2, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<matrix_binary_scalar2<E1, E2, F> >::operator ();
#endif
        typedef E1 expression1_type;
        typedef E2 expression2_type;
        typedef F functor_type;
        typedef typename E1::size_type size_type;
        typedef typename E1::difference_type difference_type;
        typedef typename F::result_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const value_type *const_pointer;
        typedef const_pointer pointer;
        typedef typename E1::const_closure_type expression1_closure_type;
        typedef E2 expression2_closure_type;
        typedef const matrix_binary_scalar2<E1, E2, F> const_self_type;
        typedef matrix_binary_scalar2<E1, E2, F> self_type;
        typedef const_self_type const_closure_type;
        typedef const_closure_type closure_type;
        typedef typename E1::orientation_category orientation_category;
        typedef typename E1::const_iterator1 const_iterator11_type;
        typedef typename E1::const_iterator2 const_iterator12_type;
        typedef typename E2::value_type const_iterator2_type;
        typedef unknown_storage_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_binary_scalar2 (): 
            e1_ (), e2_ () {}
        BOOST_UBLAS_INLINE
        matrix_binary_scalar2 (const expression1_type &e1, const expression2_type &e2): 
            e1_ (e1), e2_ (e2) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const { 
            return e1_.size1 ();
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const { 
            return e1_.size2 ();
        }
        BOOST_UBLAS_INLINE
        const expression1_closure_type &expression1 () const {
            return e1_;
        }
        BOOST_UBLAS_INLINE
        const expression2_closure_type &expression2 () const {
            return e2_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const { 
            return functor_type () (e1_ (i, j), e2_); 
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator1<const_closure_type, typename const_iterator11_type::iterator_category> const_iterator1;
        typedef const_iterator1 iterator1;
        typedef indexed_const_iterator2<const_closure_type, typename const_iterator12_type::iterator_category> const_iterator2;
        typedef const_iterator2 iterator2;
#else
        class const_iterator1;
        typedef const_iterator1 iterator1;
        class const_iterator2;
        typedef const_iterator2 iterator2;
#endif
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int rank, size_type i, size_type j) const {
            const_iterator11_type it11 (e1_.find1 (rank, i, j));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator1 (*this, it11.index1 (), it11.index2 ());
#else
            return const_iterator1 (*this, it11, e2_);
#endif
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int rank, size_type i, size_type j) const {
            const_iterator12_type it12 (e1_.find2 (rank, i, j));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator2 (*this, it12.index1 (), it12.index2 ());
#else
            return const_iterator2 (*this, it12, e2_);
#endif
        }

        // Iterators enhance the iterators of the referenced expression
        // with the binary functor.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<matrix_binary_scalar2>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename E1::const_iterator1::iterator_category>::template
                iterator_base<const_iterator1, value_type>::type {
#else
            public random_access_iterator_base<typename E1::const_iterator1::iterator_category,
                                               const_iterator1, value_type> {
#endif                                               
        public:
            typedef typename E1::const_iterator1::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_binary_scalar2::difference_type difference_type;
            typedef typename matrix_binary_scalar2::value_type value_type;
            typedef typename matrix_binary_scalar2::const_reference reference;
            typedef typename matrix_binary_scalar2::const_pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &mbs, const const_iterator11_type &it1, const const_iterator2_type &it2):
                container_const_reference<self_type> (mbs), it1_ (it1), it2_ (it2) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator1 &operator ++ () {
                ++ it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator1 &operator -- () {
                -- it1_ ;
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
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ - it.it1_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return functor_type () (*it1_, it2_);
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
                it2_ = it.it2_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ == it.it1_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator1 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ < it.it1_;
            }

        private:
            const_iterator11_type it1_;
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
        class const_iterator2:
            public container_const_reference<matrix_binary_scalar2>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename E1::const_iterator2::iterator_category>::template
                iterator_base<const_iterator2, value_type>::type {
#else
            public random_access_iterator_base<typename E1::const_iterator2::iterator_category,
                                               const_iterator2, value_type> {
#endif                                               
        public:
            typedef typename E1::const_iterator2::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_binary_scalar2::difference_type difference_type;
            typedef typename matrix_binary_scalar2::value_type value_type;
            typedef typename matrix_binary_scalar2::const_reference reference;
            typedef typename matrix_binary_scalar2::const_pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &mbs, const const_iterator12_type &it1, const const_iterator2_type &it2):
                container_const_reference<self_type> (mbs), it1_ (it1), it2_ (it2) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator2 &operator ++ () {
                ++ it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -- () {
                -- it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator += (difference_type n) {
                it1_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator2 &operator -= (difference_type n) {
                it1_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ - it.it1_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return functor_type () (*it1_, it2_);
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
                return it1_.index1 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it1_.index2 ();
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
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ == it.it1_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator2 &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                // FIXME: we shouldn't compare floats
                // BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
                return it1_ < it.it1_;
            }

        private:
            const_iterator12_type it1_;
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
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

    private:
        expression1_closure_type e1_;
        expression2_closure_type e2_;
    };

    template<class E1, class T2, class F>
    struct matrix_binary_scalar2_traits {
        typedef matrix_binary_scalar2<E1, scalar_const_reference<T2>, F> expression_type;
#ifdef BOOST_UBLAS_USE_ET
        typedef expression_type result_type; 
#else
        typedef matrix<typename F::result_type> result_type;
#endif
    };

    // (m * t) [i] [j] = m [i] [j] * t
    template<class E1, class T2>
    BOOST_UBLAS_INLINE
    typename matrix_binary_scalar2_traits<E1, T2, scalar_multiplies<typename E1::value_type, T2> >::result_type
    operator * (const matrix_expression<E1> &e1,
                const T2 &e2) {
        typedef BOOST_UBLAS_TYPENAME matrix_binary_scalar2_traits<E1, T2, scalar_multiplies<BOOST_UBLAS_TYPENAME E1::value_type, T2> >::expression_type expression_type;
        return expression_type (e1 (), e2);
    }

    // (m / t) [i] [j] = m [i] [j] / t
    template<class E1, class T2>
    BOOST_UBLAS_INLINE
    typename matrix_binary_scalar2_traits<E1, T2, scalar_divides<typename E1::value_type, T2> >::result_type
    operator / (const matrix_expression<E1> &e1,
                const T2 &e2) {
        typedef BOOST_UBLAS_TYPENAME matrix_binary_scalar2_traits<E1, T2, scalar_divides<BOOST_UBLAS_TYPENAME E1::value_type, T2> >::expression_type expression_type;
        return expression_type (e1 (), e2);
    }

    template<class E1, class E2, class F>
    class matrix_vector_binary1:
        public vector_expression<matrix_vector_binary1<E1, E2, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<matrix_vector_binary1<E1, E2, F> >::operator ();
#endif
        BOOST_STATIC_CONSTANT (int, complexity = 1);
        typedef E1 expression1_type;
        typedef E2 expression2_type;
        typedef F functor_type;
        typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
        typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
        typedef typename F::result_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const value_type *const_pointer;
        typedef const_pointer pointer;
        typedef typename E1::const_closure_type expression1_closure_type;
        typedef typename E2::const_closure_type expression2_closure_type;
        typedef const matrix_vector_binary1<E1, E2, F> const_self_type;
        typedef matrix_vector_binary1<E1, E2, F> self_type;
        typedef const_self_type const_closure_type;
        typedef const_closure_type closure_type;
        typedef typename E1::const_iterator1 const_iterator1_type;
        typedef typename E2::const_iterator const_iterator2_type;
        typedef unknown_storage_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_vector_binary1 (): 
            e1_ (), e2_ () {}
        BOOST_UBLAS_INLINE
        matrix_vector_binary1 (const expression1_type &e1, const expression2_type &e2): 
            e1_ (e1), e2_ (e2) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const { 
            return e1_.size1 (); 
        }
        BOOST_UBLAS_INLINE
        const expression1_closure_type &expression1 () const {
            return e1_;
        }
        BOOST_UBLAS_INLINE
        const expression2_closure_type &expression2 () const {
            return e2_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            return functor_type () (e1_, e2_, i); 
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator<const_closure_type, typename const_iterator1_type::iterator_category> const_iterator;
        typedef const_iterator iterator;
#else
        class const_iterator;
        typedef const_iterator iterator;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type i) const {
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            const_iterator1_type it1 (e1_.find1 (0, i, 0));
            return const_iterator (*this, it1.index1 ());
#else
            return const_iterator (*this, e1_.find1 (0, i, 0));
#endif
        }

        // Iterator simply is a pointer.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator:
            public container_const_reference<matrix_vector_binary1>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                                          typename E2::const_iterator::iterator_category>::iterator_category>::template
                iterator_base<const_iterator, value_type>::type {
#else
            public random_access_iterator_base<typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                                                 typename E2::const_iterator::iterator_category>::iterator_category,
                                               const_iterator, value_type> {
#endif                                               
        public:
            typedef typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category, 
                                                      typename E2::const_iterator::iterator_category>::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_vector_binary1::difference_type difference_type;
            typedef typename matrix_vector_binary1::value_type value_type;
            typedef typename matrix_vector_binary1::const_reference reference;
            typedef typename matrix_vector_binary1::const_pointer pointer;
#endif

            // Construction and destruction
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it1_ (), e2_begin_ (), e2_end_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &mvb, const const_iterator1_type &it1):
                container_const_reference<self_type> (mvb), it1_ (it1), e2_begin_ (mvb.expression2 ().begin ()), e2_end_ (mvb.expression2 ().end ()) {}
#else
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it1_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &mvb, const const_iterator1_type &it1):
                container_const_reference<self_type> (mvb), it1_ (it1) {}
#endif

            // Dense random access specialization
            BOOST_UBLAS_INLINE
            value_type dereference (dense_random_access_iterator_tag) const {
                const self_type &mvb = (*this) ();
#ifdef BOOST_UBLAS_USE_INDEXING
                return mvb (index ());
#elif BOOST_UBLAS_USE_ITERATING
                difference_type size = BOOST_UBLAS_SAME (mvb.expression1 ().size2 (), mvb.expression2 ().size ());
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (size, it1_.begin (), e2_begin_);
#else
                return functor_type () (size, it1_.begin (), mvb.expression2 ().begin ());
#endif
#else
                difference_type size = BOOST_UBLAS_SAME (mvb.expression1 ().size2 (), mvb.expression2 ().size ());
                if (size >= BOOST_UBLAS_ITERATOR_THRESHOLD)
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                    return functor_type () (size, it1_.begin (), e2_begin_);
#else
                    return functor_type () (size, it1_.begin (), mvb.expression2 ().begin ());
#endif
                else
                    return mvb (index ());
#endif
            }

            // Packed bidirectional specialization
            BOOST_UBLAS_INLINE
            value_type dereference (packed_random_access_iterator_tag) const {
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (it1_.begin (), it1_.end (), e2_begin_, e2_end_);
#else
                const self_type &mvb = (*this) ();
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                return functor_type () (it1_.begin (), it1_.end (),
                                        mvb.expression2 ().begin (), mvb.expression2 ().end ());
#else
                return functor_type () (boost::numeric::ublas::begin (it1_, iterator1_tag ()),
                                        boost::numeric::ublas::end (it1_, iterator1_tag ()),
                                        mvb.expression2 ().begin (), mvb.expression2 ().end ());
#endif
#endif
            }

            // Sparse bidirectional specialization
            BOOST_UBLAS_INLINE
            value_type dereference (sparse_bidirectional_iterator_tag) const {
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (it1_.begin (), it1_.end (), e2_begin_, e2_end_, sparse_bidirectional_iterator_tag ());
#else
                const self_type &mvb = (*this) ();
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                return functor_type () (it1_.begin (), it1_.end (),
                                        mvb.expression2 ().begin (), mvb.expression2 ().end (), sparse_bidirectional_iterator_tag ());
#else
                return functor_type () (boost::numeric::ublas::begin (it1_, iterator1_tag ()),
                                        boost::numeric::ublas::end (it1_, iterator1_tag ()),
                                        mvb.expression2 ().begin (), mvb.expression2 ().end (), sparse_bidirectional_iterator_tag ());
#endif
#endif
            }

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator &operator ++ () {
                ++ it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -- () {
                -- it1_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator += (difference_type n) {
                it1_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -= (difference_type n) {
                it1_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it1_ - it.it1_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return dereference (iterator_category ());
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it1_.index1 ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                container_const_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                e2_begin_ = it.e2_begin_;
                e2_end_ = it.e2_end_;
#endif
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it1_ == it.it1_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it1_ < it.it1_;
            }

        private:
            const_iterator1_type it1_;
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            // Mutable due to assignment
            /* const */ const_iterator2_type e2_begin_;
            /* const */ const_iterator2_type e2_end_;
#endif
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
        expression1_closure_type e1_;
        expression2_closure_type e2_;
    };

    template<class T1, class E1, class T2, class E2>
    struct matrix_vector_binary1_traits {
        typedef unknown_storage_tag storage_category;
        typedef row_major_tag orientation_category;
        typedef typename promote_traits<T1, T2>::promote_type promote_type;
        typedef matrix_vector_binary1<E1, E2, matrix_vector_prod1<T1, T2, promote_type> > expression_type;
#ifdef BOOST_UBLAS_USE_ET
        typedef expression_type result_type;
#else
        typedef vector<promote_type> result_type;
#endif
    };

    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_vector_binary1_traits<typename E1::value_type, E1,
                                          typename E2::value_type, E2>::result_type
    prod (const matrix_expression<E1> &e1,
          const vector_expression<E2> &e2,
          unknown_storage_tag,
          row_major_tag) {
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary1_traits<BOOST_UBLAS_TYPENAME E1::value_type, E1,
                                                                  BOOST_UBLAS_TYPENAME E2::value_type, E2>::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_vector_binary1_traits<typename E1::value_type, E1,
                                          typename E2::value_type, E2>::result_type
    prod (const matrix_expression<E1> &e1,
          const vector_expression<E2> &e2) {
        BOOST_STATIC_ASSERT (E2::complexity == 0);
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary1_traits<BOOST_UBLAS_TYPENAME E1::value_type, E1,
                                                                  BOOST_UBLAS_TYPENAME E2::value_type, E2>::storage_category storage_category;
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary1_traits<BOOST_UBLAS_TYPENAME E1::value_type, E1,
                                                                  BOOST_UBLAS_TYPENAME E2::value_type, E2>::orientation_category orientation_category;
        return prod (e1, e2, storage_category (), orientation_category ());
    }

    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_vector_binary1_traits<typename type_traits<typename E1::value_type>::precision_type, E1,
                                          typename type_traits<typename E2::value_type>::precision_type, E2>::result_type
    prec_prod (const matrix_expression<E1> &e1,
               const vector_expression<E2> &e2,
               unknown_storage_tag,
               row_major_tag) {
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary1_traits<BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E1::value_type>::precision_type, E1,
                                                                  BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E2::value_type>::precision_type, E2>::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_vector_binary1_traits<typename type_traits<typename E1::value_type>::precision_type, E1,
                                          typename type_traits<typename E2::value_type>::precision_type, E2>::result_type
    prec_prod (const matrix_expression<E1> &e1,
               const vector_expression<E2> &e2) {
        BOOST_STATIC_ASSERT (E2::complexity == 0);
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary1_traits<BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E1::value_type>::precision_type, E1,
                                                                  BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E2::value_type>::precision_type, E2>::storage_category storage_category;
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary1_traits<BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E1::value_type>::precision_type, E1,
                                                                  BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E2::value_type>::precision_type, E2>::orientation_category orientation_category;
        return prec_prod (e1, e2, storage_category (), orientation_category ());
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    prod (const matrix_expression<E1> &e1,
          const vector_expression<E2> &e2,
          V &v) {
        return v.assign (prod (e1, e2));
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    prec_prod (const matrix_expression<E1> &e1,
               const vector_expression<E2> &e2,
               V &v) {
        return v.assign (prec_prod (e1, e2));
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V
    prod (const matrix_expression<E1> &e1,
          const vector_expression<E2> &e2) {
        return V (prod (e1, e2));
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V
    prec_prod (const matrix_expression<E1> &e1,
               const vector_expression<E2> &e2) {
        return V (prec_prod (e1, e2));
    }

    template<class E1, class E2, class F>
    class matrix_vector_binary2:
        public vector_expression<matrix_vector_binary2<E1, E2, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING vector_expression<matrix_vector_binary2<E1, E2, F> >::operator ();
#endif
        BOOST_STATIC_CONSTANT (int, complexity = 1);
        typedef E1 expression1_type;
        typedef E2 expression2_type;
        typedef F functor_type;
        typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
        typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
        typedef typename F::result_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const value_type *const_pointer;
        typedef const_pointer pointer;
        typedef typename E1::const_closure_type expression1_closure_type;
        typedef typename E2::const_closure_type expression2_closure_type;
        typedef const matrix_vector_binary2<E1, E2, F> const_self_type;
        typedef matrix_vector_binary2<E1, E2, F> self_type;
        typedef const_self_type const_closure_type;
        typedef const_closure_type closure_type;
        typedef typename E1::const_iterator const_iterator1_type;
        typedef typename E2::const_iterator2 const_iterator2_type;
        typedef unknown_storage_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_vector_binary2 (): 
            e1_ (), e2_ () {}
        BOOST_UBLAS_INLINE
        matrix_vector_binary2 (const expression1_type &e1, const expression2_type &e2): 
            e1_ (e1), e2_ (e2) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size () const { 
            return e2_.size2 (); 
        }
        BOOST_UBLAS_INLINE
        const expression1_closure_type &expression1 () const {
            return e1_;
        }
        BOOST_UBLAS_INLINE
        const expression2_closure_type &expression2 () const {
            return e2_;
        }


        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type j) const { 
            return functor_type () (e1_, e2_, j); 
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator<const_closure_type, typename const_iterator2_type::iterator_category> const_iterator;
        typedef const_iterator iterator;
#else
        class const_iterator;
        typedef const_iterator iterator;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator find (size_type j) const {
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            const_iterator2_type it2 (e2_.find2 (0, 0, j));
            return const_iterator (*this, it2.index2 ());
#else
            return const_iterator (*this, e2_.find2 (0, 0, j));
#endif
        }

        // Iterator simply is a pointer.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator:
            public container_const_reference<matrix_vector_binary2>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
                                                                          typename E2::const_iterator2::iterator_category>::iterator_category>::template
                iterator_base<const_iterator, value_type>::type {
#else
            public random_access_iterator_base<typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
                                                                                 typename E2::const_iterator2::iterator_category>::iterator_category,
                                               const_iterator, value_type> {
#endif
        public:
            typedef typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
                                                      typename E2::const_iterator2::iterator_category>::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_vector_binary2::difference_type difference_type;
            typedef typename matrix_vector_binary2::value_type value_type;
            typedef typename matrix_vector_binary2::const_reference reference;
            typedef typename matrix_vector_binary2::const_pointer pointer;
#endif

            // Construction and destruction
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it2_ (), e1_begin_ (), e1_end_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &mvb, const const_iterator2_type &it2):
                container_const_reference<self_type> (mvb), it2_ (it2), e1_begin_ (mvb.expression1 ().begin ()), e1_end_ (mvb.expression1 ().end ()) {}
#else
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<self_type> (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const self_type &mvb, const const_iterator2_type &it2):
                container_const_reference<self_type> (mvb), it2_ (it2) {}
#endif

            // Dense random access specialization
            BOOST_UBLAS_INLINE
            value_type dereference (dense_random_access_iterator_tag) const {
                const self_type &mvb = (*this) ();
#ifdef BOOST_UBLAS_USE_INDEXING
                return mvb (index ());
#elif BOOST_UBLAS_USE_ITERATING
                difference_type size = BOOST_UBLAS_SAME (mvb.expression2 ().size1 (), mvb.expression1 ().size ());
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (size, e1_begin_, it2_.begin ());
#else
                return functor_type () (size, mvb.expression1 ().begin (), it2_.begin ());
#endif
#else
                difference_type size = BOOST_UBLAS_SAME (mvb.expression2 ().size1 (), mvb.expression1 ().size ());
                if (size >= BOOST_UBLAS_ITERATOR_THRESHOLD)
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                    return functor_type () (size, e1_begin_, it2_.begin ());
#else
                    return functor_type () (size, mvb.expression1 ().begin (), it2_.begin ());
#endif
                else
                    return mvb (index ());
#endif
            }

            // Packed bidirectional specialization
            BOOST_UBLAS_INLINE
            value_type dereference (packed_random_access_iterator_tag) const {
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (e1_begin_, e1_end_, it2_.begin (), it2_.end ());
#else
                const self_type &mvb = (*this) ();
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                return functor_type () (mvb.expression1 ().begin (), mvb.expression1 ().end (),
                                        it2_.begin (), it2_.end ());
#else
                return functor_type () (mvb.expression1 ().begin (), mvb.expression1 ().end (),
                                        boost::numeric::ublas::begin (it2_, iterator2_tag ()),
                                        boost::numeric::ublas::end (it2_, iterator2_tag ()));
#endif
#endif
            }

            // Sparse bidirectional specialization
            BOOST_UBLAS_INLINE
            value_type dereference (sparse_bidirectional_iterator_tag) const {
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (e1_begin_, e1_end_, it2_.begin (), it2_.end (), sparse_bidirectional_iterator_tag ());
#else
                const self_type &mvb = (*this) ();
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                return functor_type () (mvb.expression1 ().begin (), mvb.expression1 ().end (),
                                        it2_.begin (), it2_.end (), sparse_bidirectional_iterator_tag ());
#else
                return functor_type () (mvb.expression1 ().begin (), mvb.expression1 ().end (),
                                        boost::numeric::ublas::begin (it2_, iterator2_tag ()),
                                        boost::numeric::ublas::end (it2_, iterator2_tag ()), sparse_bidirectional_iterator_tag ());
#endif
#endif
            }

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator &operator ++ () {
                ++ it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -- () {
                -- it2_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator += (difference_type n) {
                it2_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -= (difference_type n) {
                it2_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it2_ - it.it2_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            reference operator * () const {
                return dereference (iterator_category ());
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                return it2_.index2 ();
            }

            // Assignment 
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                container_const_reference<self_type>::assign (&it ());
                it2_ = it.it2_;
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                e1_begin_ = it.e1_begin_;
                e1_end_ = it.e1_end_;
#endif
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it2_ == it.it2_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
                return it2_ < it.it2_;
            }

        private:
            const_iterator2_type it2_;
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            // Mutable due to assignment 
            /* const */ const_iterator1_type e1_begin_;
            /* const */ const_iterator1_type e1_end_;
#endif
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
        expression1_closure_type e1_;
        expression2_closure_type e2_;
    };

    template<class T1, class E1, class T2, class E2>
    struct matrix_vector_binary2_traits {
        typedef unknown_storage_tag storage_category;
        typedef column_major_tag orientation_category;
        typedef typename promote_traits<T1, T2>::promote_type promote_type;
        typedef matrix_vector_binary2<E1, E2, matrix_vector_prod2<T1, T2, promote_type> > expression_type;
#ifdef BOOST_UBLAS_USE_ET
        typedef expression_type result_type;
#else
        typedef vector<promote_type> result_type;
#endif
    };

    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_vector_binary2_traits<typename E1::value_type, E1,
                                          typename E2::value_type, E2>::result_type
    prod (const vector_expression<E1> &e1,
          const matrix_expression<E2> &e2,
          unknown_storage_tag,
          column_major_tag) {
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary2_traits<BOOST_UBLAS_TYPENAME E1::value_type, E1,
                                                                  BOOST_UBLAS_TYPENAME E2::value_type, E2>::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_vector_binary2_traits<typename E1::value_type, E1,
                                          typename E2::value_type, E2>::result_type
    prod (const vector_expression<E1> &e1,
          const matrix_expression<E2> &e2) {
        BOOST_STATIC_ASSERT (E1::complexity == 0);
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary2_traits<BOOST_UBLAS_TYPENAME E1::value_type, E1,
                                                                  BOOST_UBLAS_TYPENAME E2::value_type, E2>::storage_category storage_category;
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary2_traits<BOOST_UBLAS_TYPENAME E1::value_type, E1,
                                                                  BOOST_UBLAS_TYPENAME E2::value_type, E2>::orientation_category orientation_category;
        return prod (e1, e2, storage_category (), orientation_category ());
    }

    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_vector_binary2_traits<typename type_traits<typename E1::value_type>::precision_type, E1,
                                          typename type_traits<typename E2::value_type>::precision_type, E2>::result_type
    prec_prod (const vector_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               unknown_storage_tag,
               column_major_tag) {
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary2_traits<BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E1::value_type>::precision_type, E1,
                                                                  BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E2::value_type>::precision_type, E2>::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_vector_binary2_traits<typename type_traits<typename E1::value_type>::precision_type, E1,
                                          typename type_traits<typename E2::value_type>::precision_type, E2>::result_type
    prec_prod (const vector_expression<E1> &e1,
               const matrix_expression<E2> &e2) {
        BOOST_STATIC_ASSERT (E1::complexity == 0);
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary2_traits<BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E1::value_type>::precision_type, E1,
                                                                  BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E2::value_type>::precision_type, E2>::storage_category storage_category;
        typedef BOOST_UBLAS_TYPENAME matrix_vector_binary2_traits<BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E1::value_type>::precision_type, E1,
                                                                  BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E2::value_type>::precision_type, E2>::orientation_category orientation_category;
        return prec_prod (e1, e2, storage_category (), orientation_category ());
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    prod (const vector_expression<E1> &e1,
          const matrix_expression<E2> &e2,
          V &v) {
        return v.assign (prod (e1, e2));
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    prec_prod (const vector_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               V &v) {
        return v.assign (prec_prod (e1, e2));
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V
    prod (const vector_expression<E1> &e1,
          const matrix_expression<E2> &e2) {
        return V (prod (e1, e2));
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V
    prec_prod (const vector_expression<E1> &e1,
               const matrix_expression<E2> &e2) {
        return V (prec_prod (e1, e2));
    }

    template<class E1, class E2, class F>
    class matrix_matrix_binary:
        public matrix_expression<matrix_matrix_binary<E1, E2, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<matrix_matrix_binary<E1, E2, F> >::operator ();
#endif
        BOOST_STATIC_CONSTANT (int, complexity = 1);
        typedef E1 expression1_type;
        typedef E2 expression2_type;
        typedef F functor_type;
        typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
        typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
        typedef typename F::result_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const value_type *const_pointer;
        typedef const_pointer pointer;
        typedef typename E1::const_closure_type expression1_closure_type;
        typedef typename E2::const_closure_type expression2_closure_type;
        typedef const matrix_matrix_binary<E1, E2, F> const_self_type;
        typedef matrix_matrix_binary<E1, E2, F> self_type;
        typedef const_self_type const_closure_type;
        typedef const_closure_type closure_type;
        typedef unknown_orientation_tag orientation_category;
        typedef typename E1::const_iterator1 const_iterator11_type;
        typedef typename E1::const_iterator2 const_iterator12_type;
        typedef typename E2::const_iterator1 const_iterator21_type;
        typedef typename E2::const_iterator2 const_iterator22_type;
        typedef unknown_storage_tag storage_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_matrix_binary ():
            e1_ (), e2_ () {}
        BOOST_UBLAS_INLINE
        matrix_matrix_binary (const expression1_type &e1, const expression2_type &e2):
            e1_ (e1), e2_ (e2) {}

        // Accessors
        BOOST_UBLAS_INLINE
        size_type size1 () const {
            return e1_.size1 ();
        }
        BOOST_UBLAS_INLINE
        size_type size2 () const {
            return e2_.size2 ();
        }
        BOOST_UBLAS_INLINE
        const expression1_closure_type &expression1 () const {
            return e1_;
        }
        BOOST_UBLAS_INLINE
        const expression2_closure_type &expression2 () const {
            return e2_;
        }


        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            return functor_type () (e1_, e2_, i, j);
        }

#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef typename iterator_restrict_traits<typename const_iterator11_type::iterator_category,
                                                  typename const_iterator22_type::iterator_category>::iterator_category iterator_category;
        typedef indexed_const_iterator1<const_closure_type, iterator_category> const_iterator1;
        typedef const_iterator1 iterator1;
        typedef indexed_const_iterator2<const_closure_type, iterator_category> const_iterator2;
        typedef const_iterator2 iterator2;
#else
        class const_iterator1;
        typedef const_iterator1 iterator1;
        class const_iterator2;
        typedef const_iterator2 iterator2;
#endif
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<const_iterator1, value_type, const_reference> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2, value_type, const_reference> const_reverse_iterator2;
#else
        typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
        typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
#endif

        // Element lookup
        BOOST_UBLAS_INLINE
        const_iterator1 find1 (int /* rank */, size_type i, size_type j) const {
            // FIXME: sparse matrix tests fail!
            // const_iterator11_type it11 (e1_.find1 (rank, i, 0));
            const_iterator11_type it11 (e1_.find1 (0, i, 0));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator1 (*this, it11.index1 (), j);
#else
            // FIXME: sparse matrix tests fail!
            // const_iterator22_type it22 (e2_.find2 (rank, 0, j));
            const_iterator22_type it22 (e2_.find2 (0, 0, j));
            return const_iterator1 (*this, it11, it22);
#endif
        }
        BOOST_UBLAS_INLINE
        const_iterator2 find2 (int /* rank */, size_type i, size_type j) const {
            // FIXME: sparse matrix tests fail!
            // const_iterator22_type it22 (e2_.find2 (rank, 0, j));
            const_iterator22_type it22 (e2_.find2 (0, 0, j));
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
            return const_iterator2 (*this, i, it22.index2 ());
#else
            // FIXME: sparse matrix tests fail!
            // const_iterator11_type it11 (e1_.find1 (rank, i, 0));
            const_iterator11_type it11 (e1_.find1 (0, i, 0));
            return const_iterator2 (*this, it11, it22);
#endif
        }

        // Iterators simply are pointers.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
        class const_iterator1:
            public container_const_reference<matrix_matrix_binary>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                                          typename E2::const_iterator2::iterator_category>::iterator_category>::template
                iterator_base<const_iterator1, value_type>::type {
#else
            public random_access_iterator_base<typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                                                 typename E2::const_iterator2::iterator_category>::iterator_category,
                                               const_iterator1, value_type> {
#endif
        public:
            typedef typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                      typename E2::const_iterator2::iterator_category>::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_matrix_binary::difference_type difference_type;
            typedef typename matrix_matrix_binary::value_type value_type;
            typedef typename matrix_matrix_binary::const_reference reference;
            typedef typename matrix_matrix_binary::const_pointer pointer;
#endif
            typedef const_iterator2 dual_iterator_type;
            typedef const_reverse_iterator2 dual_reverse_iterator_type;

            // Construction and destruction
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), it1_ (), it2_ (), it2_begin_ (), it2_end_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &mmb, const const_iterator11_type &it1, const const_iterator22_type &it2):
                container_const_reference<self_type> (mmb), it1_ (it1), it2_ (it2), it2_begin_ (it2.begin ()), it2_end_ (it2.end ()) {}
#else
            BOOST_UBLAS_INLINE
            const_iterator1 ():
                container_const_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator1 (const self_type &mmb, const const_iterator11_type &it1, const const_iterator22_type &it2):
                container_const_reference<self_type> (mmb), it1_ (it1), it2_ (it2) {}
#endif

            // Random access specialization
            BOOST_UBLAS_INLINE
            value_type dereference (dense_random_access_iterator_tag) const {
                const self_type &mmb = (*this) ();
#ifdef BOOST_UBLAS_USE_INDEXING
                return mmb (index1 (), index2 ());
#elif BOOST_UBLAS_USE_ITERATING
                difference_type size = BOOST_UBLAS_SAME (mmb.expression1 ().size2 (), mmb.expression2 ().size1 ());
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (size, it1_.begin (), it2_begin_);
#else
                return functor_type () (size, it1_.begin (), it2_.begin ());
#endif
#else
                difference_type size = BOOST_UBLAS_SAME (mmb.expression1 ().size2 (), mmb.expression2 ().size1 ());
                if (size >= BOOST_UBLAS_ITERATOR_THRESHOLD)
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                    return functor_type () (size, it1_.begin (), it2_begin_);
#else
                    return functor_type () (size, it1_.begin (), it2_.begin ());
#endif
                else
                    return mmb (index1 (), index2 ());
#endif
            }

            // Packed bidirectional specialization
            BOOST_UBLAS_INLINE
            value_type dereference (packed_random_access_iterator_tag) const {
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (it1_.begin (), it1_.end (),
                                        it2_begin_, it2_end_, packed_random_access_iterator_tag ());
#else
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                return functor_type () (it1_.begin (), it1_.end (),
                                        it2_.begin (), it2_.end (), packed_random_access_iterator_tag ());
#else
                return functor_type () (boost::numeric::ublas::begin (it1_, iterator1_tag ()),
                                        boost::numeric::ublas::end (it1_, iterator1_tag ()),
                                        boost::numeric::ublas::begin (it2_, iterator2_tag ()),
                                        boost::numeric::ublas::end (it2_, iterator2_tag ()), packed_random_access_iterator_tag ());
#endif
#endif
            }

            // Sparse bidirectional specialization
            BOOST_UBLAS_INLINE
            value_type dereference (sparse_bidirectional_iterator_tag) const {
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (it1_.begin (), it1_.end (),
                                        it2_begin_, it2_end_, sparse_bidirectional_iterator_tag ());
#else
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                return functor_type () (it1_.begin (), it1_.end (),
                                        it2_.begin (), it2_.end (), sparse_bidirectional_iterator_tag ());
#else
                return functor_type () (boost::numeric::ublas::begin (it1_, iterator1_tag ()),
                                        boost::numeric::ublas::end (it1_, iterator1_tag ()),
                                        boost::numeric::ublas::begin (it2_, iterator2_tag ()),
                                        boost::numeric::ublas::end (it2_, iterator2_tag ()), sparse_bidirectional_iterator_tag ());
#endif
#endif
            }

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
                return dereference (iterator_category ());
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
                return it2_.index2 ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator1 &operator = (const const_iterator1 &it) {
                container_const_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
                it2_ = it.it2_;
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                it2_begin_ = it.it2_begin_;
                it2_end_ = it.it2_end_;
#endif
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
            const_iterator11_type it1_;
            // Mutable due to assignment
            /* const */ const_iterator22_type it2_;
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            /* const */ const_iterator21_type it2_begin_;
            /* const */ const_iterator21_type it2_end_;
#endif
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
        class const_iterator2:
            public container_const_reference<matrix_matrix_binary>,
#ifdef BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
            public iterator_base_traits<typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                                          typename E2::const_iterator2::iterator_category>::iterator_category>::template
                iterator_base<const_iterator2, value_type>::type {
#else
            public random_access_iterator_base<typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                                                 typename E2::const_iterator2::iterator_category>::iterator_category,
                                               const_iterator2, value_type> {
#endif
        public:
            typedef typename iterator_restrict_traits<typename E1::const_iterator1::iterator_category,
                                                      typename E2::const_iterator2::iterator_category>::iterator_category iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename matrix_matrix_binary::difference_type difference_type;
            typedef typename matrix_matrix_binary::value_type value_type;
            typedef typename matrix_matrix_binary::const_reference reference;
            typedef typename matrix_matrix_binary::const_pointer pointer;
#endif
            typedef const_iterator1 dual_iterator_type;
            typedef const_reverse_iterator1 dual_reverse_iterator_type;

            // Construction and destruction
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), it1_ (), it2_ (), it1_begin_ (), it1_end_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &mmb, const const_iterator11_type &it1, const const_iterator22_type &it2):
                container_const_reference<self_type> (mmb), it1_ (it1), it2_ (it2), it1_begin_ (it1.begin ()), it1_end_ (it1.end ()) {}
#else
            BOOST_UBLAS_INLINE
            const_iterator2 ():
                container_const_reference<self_type> (), it1_ (), it2_ () {}
            BOOST_UBLAS_INLINE
            const_iterator2 (const self_type &mmb, const const_iterator11_type &it1, const const_iterator22_type &it2):
                container_const_reference<self_type> (mmb), it1_ (it1), it2_ (it2) {}
#endif

            // Random access specialization
            BOOST_UBLAS_INLINE
            value_type dereference (dense_random_access_iterator_tag) const {
                const self_type &mmb = (*this) ();
#ifdef BOOST_UBLAS_USE_INDEXING
                return mmb (index1 (), index2 ());
#elif BOOST_UBLAS_USE_ITERATING
                difference_type size = BOOST_UBLAS_SAME (mmb.expression1 ().size2 (), mmb.expression2 ().size1 ());
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (size, it1_begin_, it2_.begin ());
#else
                return functor_type () (size, it1_.begin (), it2_.begin ());
#endif
#else
                difference_type size = BOOST_UBLAS_SAME (mmb.expression1 ().size2 (), mmb.expression2 ().size1 ());
                if (size >= BOOST_UBLAS_ITERATOR_THRESHOLD)
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                    return functor_type () (size, it1_begin_, it2_.begin ());
#else
                    return functor_type () (size, it1_.begin (), it2_.begin ());
#endif
                else
                    return mmb (index1 (), index2 ());
#endif
            }

            // Packed bidirectional specialization
            BOOST_UBLAS_INLINE
            value_type dereference (packed_random_access_iterator_tag) const {
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (it1_begin_, it1_end_,
                                        it2_.begin (), it2_.end (), packed_random_access_iterator_tag ());
#else
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                return functor_type () (it1_.begin (), it1_.end (),
                                        it2_.begin (), it2_.end (), packed_random_access_iterator_tag ());
#else
                return functor_type () (boost::numeric::ublas::begin (it1_, iterator1_tag ()),
                                        boost::numeric::ublas::end (it1_, iterator1_tag ()),
                                        boost::numeric::ublas::begin (it2_, iterator2_tag ()),
                                        boost::numeric::ublas::end (it2_, iterator2_tag ()), packed_random_access_iterator_tag ());
#endif
#endif
            }

            // Sparse bidirectional specialization
            BOOST_UBLAS_INLINE
            value_type dereference (sparse_bidirectional_iterator_tag) const {
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                return functor_type () (it1_begin_, it1_end_,
                                        it2_.begin (), it2_.end (), sparse_bidirectional_iterator_tag ());
#else
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                return functor_type () (it1_.begin (), it1_.end (),
                                        it2_.begin (), it2_.end (), sparse_bidirectional_iterator_tag ());
#else
                return functor_type () (boost::numeric::ublas::begin (it1_, iterator1_tag ()),
                                        boost::numeric::ublas::end (it1_, iterator1_tag ()),
                                        boost::numeric::ublas::begin (it2_, iterator2_tag ()),
                                        boost::numeric::ublas::end (it2_, iterator2_tag ()), sparse_bidirectional_iterator_tag ());
#endif
#endif
            }

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
                return dereference (iterator_category ());
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
                return it1_.index1 ();
            }
            BOOST_UBLAS_INLINE
            size_type index2 () const {
                return it2_.index2 ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator2 &operator = (const const_iterator2 &it) {
                container_const_reference<self_type>::assign (&it ());
                it1_ = it.it1_;
                it2_ = it.it2_;
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
                it1_begin_ = it.it1_begin_;
                it1_end_ = it.it1_end_;
#endif
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
            // Mutable due to assignment
            /* const */ const_iterator11_type it1_;
            const_iterator22_type it2_;
#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
            /* const */ const_iterator12_type it1_begin_;
            /* const */ const_iterator12_type it1_end_;
#endif
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
        const_reverse_iterator2 rbegin2 () const {
            return const_reverse_iterator2 (end2 ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator2 rend2 () const {
            return const_reverse_iterator2 (begin2 ());
        }

    private:
        expression1_closure_type e1_;
        expression2_closure_type e2_;
    };

    template<class T1, class E1, class T2, class E2>
    struct matrix_matrix_binary_traits {
        typedef unknown_storage_tag storage_category;
        typedef unknown_orientation_tag orientation_category;
        typedef typename promote_traits<T1, T2>::promote_type promote_type;
        typedef matrix_matrix_binary<E1, E2, matrix_matrix_prod<T1, T2, promote_type> > expression_type;
#ifdef BOOST_UBLAS_USE_ET
        typedef expression_type result_type;
#else
        typedef matrix<promote_type> result_type;
#endif
    };

    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_matrix_binary_traits<typename E1::value_type, E1,
                                         typename E2::value_type, E2>::result_type
    prod (const matrix_expression<E1> &e1,
          const matrix_expression<E2> &e2,
          unknown_storage_tag,
          unknown_orientation_tag) {
        typedef BOOST_UBLAS_TYPENAME matrix_matrix_binary_traits<BOOST_UBLAS_TYPENAME E1::value_type, E1,
                                                                 BOOST_UBLAS_TYPENAME E2::value_type, E2>::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_matrix_binary_traits<typename E1::value_type, E1,
                                         typename E2::value_type, E2>::result_type
    prod (const matrix_expression<E1> &e1,
          const matrix_expression<E2> &e2) {
        BOOST_STATIC_ASSERT (E1::complexity == 0 && E2::complexity == 0);
        typedef BOOST_UBLAS_TYPENAME matrix_matrix_binary_traits<BOOST_UBLAS_TYPENAME E1::value_type, E1,
                                                                 BOOST_UBLAS_TYPENAME E2::value_type, E2>::storage_category storage_category;
        typedef BOOST_UBLAS_TYPENAME matrix_matrix_binary_traits<BOOST_UBLAS_TYPENAME E1::value_type, E1,
                                                                 BOOST_UBLAS_TYPENAME E2::value_type, E2>::orientation_category orientation_category;
        return prod (e1, e2, storage_category (), orientation_category ());
    }

    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_matrix_binary_traits<typename type_traits<typename E1::value_type>::precision_type, E1,
                                         typename type_traits<typename E2::value_type>::precision_type, E2>::result_type
    prec_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               unknown_storage_tag,
               unknown_orientation_tag) {
        typedef BOOST_UBLAS_TYPENAME matrix_matrix_binary_traits<BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E1::value_type>::precision_type, E1,
                                                                 BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E2::value_type>::precision_type, E2>::expression_type expression_type;
        return expression_type (e1 (), e2 ());
    }

    // Dispatcher
    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    typename matrix_matrix_binary_traits<typename type_traits<typename E1::value_type>::precision_type, E1,
                                         typename type_traits<typename E2::value_type>::precision_type, E2>::result_type
    prec_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2) {
        BOOST_STATIC_ASSERT (E1::complexity == 0 && E2::complexity == 0);
        typedef BOOST_UBLAS_TYPENAME matrix_matrix_binary_traits<BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E1::value_type>::precision_type, E1,
                                                                 BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E2::value_type>::precision_type, E2>::storage_category storage_category;
        typedef BOOST_UBLAS_TYPENAME matrix_matrix_binary_traits<BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E1::value_type>::precision_type, E1,
                                                                 BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME E2::value_type>::precision_type, E2>::orientation_category orientation_category;
        return prec_prod (e1, e2, storage_category (), orientation_category ());
    }

    template<class M, class E1, class E2>
    BOOST_UBLAS_INLINE
    M &
    prod (const matrix_expression<E1> &e1,
          const matrix_expression<E2> &e2,
          M &m) {
        return m.assign (prod (e1, e2));
    }

    template<class M, class E1, class E2>
    BOOST_UBLAS_INLINE
    M &
    prec_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               M &m) {
        return m.assign (prec_prod (e1, e2));
    }

    template<class M, class E1, class E2>
    BOOST_UBLAS_INLINE
    M
    prod (const matrix_expression<E1> &e1,
          const matrix_expression<E2> &e2) {
        return M (prod (e1, e2));
    }

    template<class M, class E1, class E2>
    BOOST_UBLAS_INLINE
    M
    prec_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2) {
        return M (prec_prod (e1, e2));
    }

    template<class E, class F>
    class matrix_scalar_unary:
        public scalar_expression<typename F::result_type> {
    public:
        typedef E expression_type;
        typedef F functor_type;
        typedef typename F::result_type value_type;
        typedef typename E::const_closure_type expression_closure_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        matrix_scalar_unary ():
            e_ () {}
        BOOST_UBLAS_INLINE
        matrix_scalar_unary (const expression_type &e):
            e_ (e) {}

        // Accessors
        BOOST_UBLAS_INLINE
        const expression_closure_type &expression () const {
            return e_;
        }

        BOOST_UBLAS_INLINE
        operator value_type () const {
            return functor_type () (e_);
        }

    private:
        expression_closure_type e_;
    };

    template<class E, class F>
    struct matrix_scalar_unary_traits {
        typedef matrix_scalar_unary<E, F> expression_type;
#ifdef BOOST_UBLAS_USE_ET
         typedef expression_type result_type;
#else
         typedef typename F::result_type result_type;
#endif
    };

    template<class E>
    BOOST_UBLAS_INLINE
    typename matrix_scalar_unary_traits<E, matrix_norm_1<typename E::value_type> >::result_type
    norm_1 (const matrix_expression<E> &e) {
        typedef BOOST_UBLAS_TYPENAME matrix_scalar_unary_traits<E, matrix_norm_1<BOOST_UBLAS_TYPENAME E::value_type> >::expression_type expression_type;
        return expression_type (e ());
    }

    template<class E>
    BOOST_UBLAS_INLINE
    typename matrix_scalar_unary_traits<E, matrix_norm_frobenius<typename E::value_type> >::result_type
    norm_frobenius (const matrix_expression<E> &e) {
        typedef BOOST_UBLAS_TYPENAME matrix_scalar_unary_traits<E, matrix_norm_frobenius<BOOST_UBLAS_TYPENAME E::value_type> >::expression_type expression_type;
        return expression_type (e ());
    }

    template<class E>
    BOOST_UBLAS_INLINE
    typename matrix_scalar_unary_traits<E, matrix_norm_inf<typename E::value_type> >::result_type
    norm_inf (const matrix_expression<E> &e) {
        typedef BOOST_UBLAS_TYPENAME matrix_scalar_unary_traits<E, matrix_norm_inf<BOOST_UBLAS_TYPENAME E::value_type> >::expression_type expression_type;
        return expression_type (e ());
    }

}}}

#endif

























