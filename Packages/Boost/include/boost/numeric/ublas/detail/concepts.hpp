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

#ifndef BOOST_UBLAS_CONCEPTS_H
#define BOOST_UBLAS_CONCEPTS_H

#include <boost/concept_check.hpp>

// #define INTERNAL_STORAGE
// #define INTERNAL_STORAGE_SPARSE
// #define INTERNAL_VECTOR
// #define INTERNAL_VECTOR_PROXY
// #define INTERNAL_VECTOR_SPARSE
// #define INTERNAL_MATRIX
// #define INTERNAL_MATRIX_PROXY
// #define INTERNAL_BANDED
// #define INTERNAL_TRIANGULAR
// #define INTERNAL_SYMMETRIC
// #define INTERNAL_HERMITIAN
// #define INTERNAL_MATRIX_SPARSE
// #define INTERNAL_VECTOR_EXPRESSION
// #define INTERNAL_MATRIX_EXPRESSION

// #define EXTERNAL

// Concept checks based on ideas of Jeremy Siek

namespace boost { namespace numeric { namespace ublas {


    template<class T>
    struct AssignableConcept {
        typedef T value_type;

        static void constraints (value_type t) {
            // Copy Constructor
            value_type c1 (t);
            value_type c2 = t;
            // Assignment
            value_type a = t;
            std::swap (c1, c2);
            ignore_unused_variable_warning (a);
        }
    };

    template<class T>
    struct EqualityComparableConcept {
        typedef T value_type;

        static void constraints (const value_type t) {
            bool b;
            // Equality
            b = t == t;
            // Inequality
            b = t != t;
            ignore_unused_variable_warning (b);
        }
    };

    template<class T>
    struct LessThanComparableConcept {
        typedef T value_type;

        static void constraints (const value_type t) {
            bool b;
            b = t < t;
            b = t <= t;
            b = t >= t;
            b = t > t;
            ignore_unused_variable_warning (b);
        }
    };

    template<class T>
    struct DefaultConstructibleConcept {
        typedef T value_type;

        static void constraints () {
            // Default Constructor
            static value_type c1 = value_type ();
            static value_type c2;
            ignore_unused_variable_warning (c1);
            ignore_unused_variable_warning (c2);
        }
    };

    template<class I, class T = typename std::iterator_traits<I>::value_type>
    struct BidirectionalIteratorConcept {
        typedef I iterator_type;
        
        typedef typename std::iterator_traits<I>::iterator_category iterator_category;
        typedef typename std::iterator_traits<I>::difference_type difference_type;
        typedef typename std::iterator_traits<I>::value_type value_type;
        typedef typename std::iterator_traits<I>::reference reference;
        typedef typename std::iterator_traits<I>::pointer pointer;

        static void constraints () {
            AssignableConcept<iterator_type>::constraints (iterator_type ());
            EqualityComparableConcept<iterator_type>::constraints (iterator_type ());
            DefaultConstructibleConcept<iterator_type>::constraints ();
            iterator_type it = iterator_type ();

            // Associated types - assume constructable
            iterator_category c;
            difference_type d (0);
            pointer p (0);
            // Dereference
            reference r (*it);
            value_type t (r);
            // Member access
            // FIXME it->m;
            // Preincrement
            ++ it;
            // Postincrement
            it ++;
            // Predecrement
            -- it;
            // Postdecrement
            it --;
            ignore_unused_variable_warning (t);
            ignore_unused_variable_warning (c);
            ignore_unused_variable_warning (d);
            ignore_unused_variable_warning (p);
            ignore_unused_variable_warning (r);
        }
    };

    template<class I, class T = typename std::iterator_traits<I>::value_type>
    struct MutableBidirectionalIteratorConcept {
        typedef I iterator_type;
        typedef T value_type;

        static void constraints () {
            BidirectionalIteratorConcept<iterator_type, value_type>::constraints ();
            iterator_type it = iterator_type ();
            value_type t = value_type ();
            // Dereference assignment
            *it = t;
            ignore_unused_variable_warning (t);
        }
    };

    template<class I, class D = typename std::iterator_traits<I>::difference_type, class T = typename std::iterator_traits<I>::value_type>
    struct RandomAccessIteratorConcept {
        typedef I iterator_type;
        typedef D difference_type;
        typedef T value_type;

        static void constraints () {
            LessThanComparableConcept<iterator_type>::constraints (iterator_type ());
            BidirectionalIteratorConcept<iterator_type, value_type>::constraints ();
            iterator_type it = iterator_type (), it1 = iterator_type (), it2 = iterator_type ();
            difference_type n (0);
            value_type t;
            // Forward motion
            it += n;
            // Iterator addition
            it = it + n;
            iterator_type itp (it + n);
            // Backward motion
            it -= n;
            // Iterator subtraction
            it = it - n;
            iterator_type itm (it - n);
            // Difference
            n = it1 - it2;
            // Element operator
#ifdef BOOST_UBLAS_ITERATOR_IS_INDEXABLE
            t = it [n];
#endif
            t = *(it + n);
            ignore_unused_variable_warning (itp);
            ignore_unused_variable_warning (itm);
            ignore_unused_variable_warning (t);
        }
    };

    template<class I, class D = typename std::iterator_traits<I>::difference_type, class T = typename std::iterator_traits<I>::value_type>
    struct MutableRandomAccessIteratorConcept {
        typedef I iterator_type;
        typedef D difference_type;
        typedef T value_type;

        static void constraints () {
            MutableBidirectionalIteratorConcept<iterator_type, value_type>::constraints ();
            RandomAccessIteratorConcept<iterator_type, difference_type, value_type>::constraints ();
            iterator_type it = iterator_type ();
            difference_type n (0);
            value_type t = value_type ();
            // Element assignment
#ifdef BOOST_UBLAS_ITERATOR_IS_INDEXABLE
            it [n] = t;
#endif
            *(it + n) = t;
        }
    };

    template<class I>
    struct Indexed1DIteratorConcept {
        typedef I iterator_type;

        static void constraints () {
            iterator_type it = iterator_type ();
            // Index
            it.index ();
        }
    };

    template<class I>
    struct IndexedBidirectional1DIteratorConcept {
        typedef I iterator_type;

        static void constraints () {
            BidirectionalIteratorConcept<iterator_type>::constraints ();
            Indexed1DIteratorConcept<iterator_type>::constraints ();
        }
    };

    template<class I>
    struct MutableIndexedBidirectional1DIteratorConcept {
        typedef I iterator_type;

        static void constraints () {
            MutableBidirectionalIteratorConcept<iterator_type>::constraints ();
            Indexed1DIteratorConcept<iterator_type>::constraints ();
        }
    };

    template<class I>
    struct IndexedRandomAccess1DIteratorConcept {
        typedef I iterator_type;

        static void constraints () {
            RandomAccessIteratorConcept<iterator_type>::constraints ();
            Indexed1DIteratorConcept<iterator_type>::constraints ();
        }
    };

    template<class I>
    struct MutableIndexedRandomAccess1DIteratorConcept {
        typedef I iterator_type;

        static void constraints () {
            MutableRandomAccessIteratorConcept<iterator_type>::constraints ();
            Indexed1DIteratorConcept<iterator_type>::constraints ();
        }
    };

    template<class I>
    struct Indexed2DIteratorConcept {
        typedef I iterator_type;
        typedef typename I::dual_iterator_type dual_iterator_type;
        typedef typename I::dual_reverse_iterator_type dual_reverse_iterator_type;

        static void constraints () {
            iterator_type it = iterator_type ();
            // Indices
            it.index1 ();
            it.index2 ();
            // Iterator begin/end
            dual_iterator_type it_begin (it.begin ());
            dual_iterator_type it_end (it.end ());
            // Reverse iterator begin/end
            dual_reverse_iterator_type it_rbegin (it.rbegin ());
            dual_reverse_iterator_type it_rend (it.rend ());
            ignore_unused_variable_warning (it_begin);
            ignore_unused_variable_warning (it_end);
            ignore_unused_variable_warning (it_rbegin);
            ignore_unused_variable_warning (it_rend);
        }
    };

    template<class I1, class I2>
    struct IndexedBidirectional2DIteratorConcept {
        typedef I1 iterator1_type;
        typedef I2 iterator2_type;

        static void constraints () {
            BidirectionalIteratorConcept<iterator1_type>::constraints ();
            BidirectionalIteratorConcept<iterator2_type>::constraints ();
            Indexed2DIteratorConcept<iterator1_type>::constraints ();
            Indexed2DIteratorConcept<iterator2_type>::constraints ();
        }
    };

    template<class I1, class I2>
    struct MutableIndexedBidirectional2DIteratorConcept {
        typedef I1 iterator1_type;
        typedef I2 iterator2_type;

        static void constraints () {
            MutableBidirectionalIteratorConcept<iterator1_type>::constraints ();
            MutableBidirectionalIteratorConcept<iterator2_type>::constraints ();
            Indexed2DIteratorConcept<iterator1_type>::constraints ();
            Indexed2DIteratorConcept<iterator2_type>::constraints ();
        }
    };

    template<class I1, class I2>
    struct IndexedRandomAccess2DIteratorConcept {
        typedef I1 iterator1_type;
        typedef I2 iterator2_type;

        static void constraints () {
            RandomAccessIteratorConcept<iterator1_type>::constraints ();
            RandomAccessIteratorConcept<iterator2_type>::constraints ();
            Indexed2DIteratorConcept<iterator1_type>::constraints ();
            Indexed2DIteratorConcept<iterator2_type>::constraints ();
        }
    };

    template<class I1, class I2>
    struct MutableIndexedRandomAccess2DIteratorConcept {
        typedef I1 iterator1_type;
        typedef I2 iterator2_type;

        static void constraints () {
            MutableRandomAccessIteratorConcept<iterator1_type>::constraints ();
            MutableRandomAccessIteratorConcept<iterator2_type>::constraints ();
            Indexed2DIteratorConcept<iterator1_type>::constraints ();
            Indexed2DIteratorConcept<iterator2_type>::constraints ();
        }
    };

    template<class C>
    struct ContainerConcept {
        typedef C container_type;
        typedef typename C::size_type size_type;
        typedef typename C::const_iterator const_iterator_type;

        static void constraints () {
            DefaultConstructibleConcept<container_type>::constraints ();
            container_type c = container_type ();
            size_type n (0);
            // Beginning of range
            const_iterator_type cit_begin (c.begin ());
            // End of range
            const_iterator_type cit_end (c.end ());
            // Size
            n = c.size ();
            ignore_unused_variable_warning (cit_end);
            ignore_unused_variable_warning (cit_begin);
            ignore_unused_variable_warning (n);
        }
    };

    template<class C>
    struct MutableContainerConcept {
        typedef C container_type;
        typedef typename C::iterator iterator_type;

        static void constraints () {
            AssignableConcept<container_type>::constraints (container_type ());
            ContainerConcept<container_type>::constraints ();
            container_type c = container_type (), c1 = container_type (), c2 = container_type ();
            // Beginning of range
            iterator_type it_begin (c.begin ());
            // End of range
            iterator_type it_end (c.end ());
            // Swap
            c1.swap (c2);
            ignore_unused_variable_warning (it_end);
            ignore_unused_variable_warning (it_begin);
        }
    };

    template<class C>
    struct ReversibleContainerConcept {
        typedef C container_type;
        typedef typename C::const_reverse_iterator const_reverse_iterator_type;

        static void constraints () {
            ContainerConcept<container_type>::constraints ();
            const container_type cc = container_type ();
            // Beginning of reverse range
            const_reverse_iterator_type crit_begin (cc.rbegin ());
            // End of reverse range
            const_reverse_iterator_type crit_end (cc.rend ());
            ignore_unused_variable_warning (crit_end);
            ignore_unused_variable_warning (crit_begin);
        }
    };

    template<class C>
    struct MutableReversibleContainerConcept {
        typedef C container_type;
        typedef typename C::reverse_iterator reverse_iterator_type;

        static void constraints () {
            MutableContainerConcept<container_type>::constraints ();
            ReversibleContainerConcept<container_type>::constraints ();
            container_type c = container_type ();
            // Beginning of reverse range
            reverse_iterator_type rit_begin (c.rbegin ());
            // End of reverse range
            reverse_iterator_type rit_end (c.rend ());
            ignore_unused_variable_warning (rit_end);
            ignore_unused_variable_warning (rit_begin);
        }
    };

    template<class C>
    struct RandomAccessContainerConcept {
        typedef C container_type;
        typedef typename C::size_type size_type;
        typedef typename C::value_type value_type;

        static void constraints () {
            ReversibleContainerConcept<container_type>::constraints ();
            container_type c = container_type ();
            size_type n (0);
            value_type t = value_type ();
            // Element access
            t = c [n];
            ignore_unused_variable_warning (t);
        }
    };

    template<class C>
    struct MutableRandomAccessContainerConcept {
        typedef C container_type;
        typedef typename C::size_type size_type;
        typedef typename C::value_type value_type;

        static void constraints () {
            MutableReversibleContainerConcept<container_type>::constraints ();
            RandomAccessContainerConcept<container_type>::constraints ();
            container_type c = container_type ();
            size_type n (0);
            value_type t = value_type ();
            // Element access
            c [n] = t;
        }
    };

    template<class C>
    struct StorageArrayConcept {
        typedef C container_type;
        typedef typename C::size_type size_type;
        typedef typename C::value_type value_type;

        static void constraints () {
            RandomAccessContainerConcept<container_type>::constraints ();
            size_type n (0);
            // Sizing constructor
            container_type c = container_type (n);
            // Initialised sizing constructor
            container_type (n, value_type (5));
            ignore_unused_variable_warning (c);
        }
    };

    template<class C>
    struct MutableStorageArrayConcept {
        typedef C container_type;
        typedef typename C::size_type size_type;
        typedef typename C::value_type value_type;
        typedef typename C::iterator iterator_type;

        static void constraints () {
            MutableRandomAccessContainerConcept<container_type>::constraints ();
            size_type n (0);
            // Sizing constructor
            container_type c = container_type (n);
            // Initialised sizing constructor
            c = container_type (n, value_type (3));
            // Resize
            c.resize (n, value_type (5));
            // Resize - none preserving
            c.resize (n);
        }
    };

    template<class C>
    struct StorageSparseConcept {
        typedef C container_type;
        typedef typename C::size_type size_type;

        static void constraints () {
            ReversibleContainerConcept<container_type>::constraints ();
        }
    };

    template<class C>
    struct MutableStorageSparseConcept {
        typedef C container_type;
        typedef typename C::size_type size_type;
        typedef typename C::value_type value_type;
        typedef typename C::iterator iterator_type;

        static void constraints () {
            MutableReversibleContainerConcept<container_type>::constraints ();
            container_type c = container_type ();
            value_type t = value_type ();
            iterator_type it = iterator_type (), it1 = iterator_type (), it2 = iterator_type ();
            // Insert
            c.insert (it, t);
            // Erase
            c.erase (it);
            // Range erase
            c.erase (it1, it2);
            // Clear
            c.clear ();
        }
    };

    template<class G>
    struct IndexSetConcept {
        typedef G generator_type;
        typedef typename G::size_type size_type;
        typedef typename G::value_type value_type;

        static void constraints () {
            DefaultConstructibleConcept<generator_type>::constraints ();
            ReversibleContainerConcept<generator_type>::constraints ();
            generator_type g = generator_type ();
            size_type n (0);
            value_type t;
            // Element access
            t = g (n);
            ignore_unused_variable_warning (t);
        }
    };

    template<class S>
    struct ScalarExpressionConcept {
        typedef S scalar_type;
        typedef typename S::value_type value_type;

        static void constraints () {
            DefaultConstructibleConcept<scalar_type>::constraints ();
            scalar_type s = scalar_type ();
            value_type t;
            // Conversion
            t = s;
            ignore_unused_variable_warning (t);
        }
    };

    template<class V>
    struct VectorExpressionConcept {
        typedef V vector_type;
        typedef typename V::type_category type_category;
        typedef typename V::size_type size_type;
        typedef typename V::value_type value_type;
        typedef typename V::const_iterator const_iterator_type;
        typedef typename V::const_reverse_iterator const_reverse_iterator_type;

        static void constraints () {
            DefaultConstructibleConcept<vector_type>::constraints ();
            vector_type v = vector_type ();
            size_type n (0), i (0);
            value_type t;
            // Find (internal?)
            const_iterator_type cit (v.find (i));
            // Beginning of range
            const_iterator_type cit_begin (v.begin ());
            // End of range
            const_iterator_type cit_end (v.end ());
            // Size
            n = v.size ();
            // Beginning of reverse range
            const vector_type cv = vector_type ();
            const_reverse_iterator_type crit_begin (cv.rbegin ());
            // End of reverse range
            const_reverse_iterator_type crit_end (cv.rend ());
            // Element access
            t = v (i);
            ignore_unused_variable_warning (n);
            ignore_unused_variable_warning (cit);
            ignore_unused_variable_warning (cit_begin);
            ignore_unused_variable_warning (cit_end);
            ignore_unused_variable_warning (crit_begin);
            ignore_unused_variable_warning (crit_end);
            ignore_unused_variable_warning (t);
        }
    };

    template<class V>
    struct MutableVectorExpressionConcept {
        typedef V vector_type;
        typedef typename V::size_type size_type;
        typedef typename V::value_type value_type;
        typedef typename V::iterator iterator_type;
        typedef typename V::reverse_iterator reverse_iterator_type;

        static void constraints () {
            AssignableConcept<vector_type>::constraints (vector_type ());
            VectorExpressionConcept<vector_type>::constraints ();
            vector_type v = vector_type (), v1 = vector_type (), v2 = vector_type ();
            size_type i (0);
            value_type t = value_type ();
            // Find (internal?)
            iterator_type it (v.find (i));
            // Beginning of range
            iterator_type it_begin (v.begin ());
            // End of range
            iterator_type it_end (v.end ());
            // Swap
            v1.swap (v2);
            // Beginning of reverse range
            reverse_iterator_type rit_begin (v.rbegin ());
            // End of reverse range
            reverse_iterator_type rit_end (v.rend ());
            // Assignments
            v2 = v1;
            v2.assign (v1);
            v2 += v1;
            v2.plus_assign (v1);
            v2 -= v1;
            v2.minus_assign (v1);
            v *= t;
            ignore_unused_variable_warning (it);
            ignore_unused_variable_warning (it_begin);
            ignore_unused_variable_warning (it_end);
            ignore_unused_variable_warning (rit_begin);
            ignore_unused_variable_warning (rit_end);
        }
    };

    template<class M>
    struct MatrixExpressionConcept {
        typedef M matrix_type;
        typedef typename M::type_category type_category;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;
        typedef typename M::const_iterator1 const_iterator1_type;
        typedef typename M::const_iterator2 const_iterator2_type;
        typedef typename M::const_reverse_iterator1 const_reverse_iterator1_type;
        typedef typename M::const_reverse_iterator2 const_reverse_iterator2_type;

        static void constraints () {
            DefaultConstructibleConcept<matrix_type>::constraints ();
            matrix_type m = matrix_type ();
            size_type n (0), i (0), j (0);
            value_type t;
            // Find (internal?)
            const_iterator1_type cit1 (m.find1 (0, i, j));
            const_iterator2_type cit2 (m.find2 (0, i, j));
            // Beginning of range
            const_iterator1_type cit1_begin (m.begin1 ());
            const_iterator2_type cit2_begin (m.begin2 ());
            // End of range
            const_iterator1_type cit1_end (m.end1 ());
            const_iterator2_type cit2_end (m.end2 ());
            // Size
            n = m.size1 ();
            n = m.size2 ();
            // Beginning of reverse range
            const matrix_type cm = matrix_type ();
            const_reverse_iterator1_type crit1_begin (cm.rbegin1 ());
            const_reverse_iterator2_type crit2_begin (cm.rbegin2 ());
            // End of reverse range
            const_reverse_iterator1_type crit1_end (cm.rend1 ());
            const_reverse_iterator2_type crit2_end (cm.rend2 ());
            // Element access
            t = m (i, j);
            ignore_unused_variable_warning (n);
            ignore_unused_variable_warning (cit1);
            ignore_unused_variable_warning (cit2);
            ignore_unused_variable_warning (cit1_begin);
            ignore_unused_variable_warning (cit2_begin);
            ignore_unused_variable_warning (cit1_end);
            ignore_unused_variable_warning (cit2_end);
            ignore_unused_variable_warning (crit1_begin);
            ignore_unused_variable_warning (crit2_begin);
            ignore_unused_variable_warning (crit1_end);
            ignore_unused_variable_warning (crit2_end);
            ignore_unused_variable_warning (t);
        }
    };

    template<class M>
    struct MutableMatrixExpressionConcept {
        typedef M matrix_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;
        typedef typename M::iterator1 iterator1_type;
        typedef typename M::iterator2 iterator2_type;
        typedef typename M::reverse_iterator1 reverse_iterator1_type;
        typedef typename M::reverse_iterator2 reverse_iterator2_type;

        static void constraints () {
            AssignableConcept<matrix_type>::constraints (matrix_type ());
            MatrixExpressionConcept<matrix_type>::constraints ();
            matrix_type m = matrix_type (), m1 = matrix_type (), m2 = matrix_type ();
            size_type i (0), j (0);
            value_type t = value_type ();
            // Find (internal?)
            iterator1_type it1 (m.find1 (0, i, j));
            iterator2_type it2 (m.find2 (0, i, j));
            // Beginning of range
            iterator1_type it1_begin (m.begin1 ());
            iterator2_type it2_begin (m.begin2 ());
            // End of range
            iterator1_type it1_end (m.end1 ());
            iterator2_type it2_end (m.end2 ());
            // Swap
            m1.swap (m2);
            // Beginning of reverse range
            reverse_iterator1_type rit1_begin (m.rbegin1 ());
            reverse_iterator2_type rit2_begin (m.rbegin2 ());
            // End of reverse range
            reverse_iterator1_type rit1_end (m.rend1 ());
            reverse_iterator2_type rit2_end (m.rend2 ());
            // Assignments
            m2 = m1;
            m2.assign (m1);
            m2 += m1;
            m2.plus_assign (m1);
            m2 -= m1;
            m2.minus_assign (m1);
            m *= t;
            ignore_unused_variable_warning (it1);
            ignore_unused_variable_warning (it2);
            ignore_unused_variable_warning (it1_begin);
            ignore_unused_variable_warning (it2_begin);
            ignore_unused_variable_warning (it1_end);
            ignore_unused_variable_warning (it2_end);
            ignore_unused_variable_warning (rit1_begin);
            ignore_unused_variable_warning (rit2_begin);
            ignore_unused_variable_warning (rit1_end);
            ignore_unused_variable_warning (rit2_end);
        }
    };

    template<class V>
    struct VectorConcept {
        typedef V vector_type;
        typedef typename V::size_type size_type;

        static void constraints () {
            VectorExpressionConcept<vector_type>::constraints ();
            size_type n (0);
            // Sizing constructor
            vector_type v (n);
        }
    };

    template<class V>
    struct MutableVectorConcept {
        typedef V vector_type;
        typedef typename V::size_type size_type;
        typedef typename V::value_type value_type;

        static void constraints () {
            VectorConcept<vector_type>::constraints ();
            MutableVectorExpressionConcept<vector_type>::constraints ();
            size_type n (0);
            value_type t = value_type ();
            size_type i (0);
            // Sizing constructor
            vector_type v (n);
            // Insert
            v.insert (i, t);
            // Erase
            v.erase (i);
            // Clear
            v.clear ();
            // Resize
            v.resize (n);
        }
    };

    template<class M>
    struct MatrixConcept {
        typedef M matrix_type;
        typedef typename M::size_type size_type;

        static void constraints () {
            MatrixExpressionConcept<matrix_type>::constraints ();
            size_type n (0);
            // Sizing constructor
            matrix_type m (n, n);
        }
    };

    template<class M>
    struct MutableMatrixConcept {
        typedef M matrix_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;

        static void constraints () {
            MatrixConcept<matrix_type>::constraints ();
            MutableMatrixExpressionConcept<matrix_type>::constraints ();
            size_type n (0);
            value_type t = value_type ();
            size_type i (0), j (0);
            // Sizing constructor
            matrix_type m (n, n);
            // Insert
            m.insert (i, j, t);
            // Erase
            m.erase (i, j);
            // Clear
            m.clear ();
            // Resize
            m.resize (n, n);
            m.resize (n, n, false);
        }
    };

    template<class T>
    T
    ZeroElement (T);
    template<>
    float
    ZeroElement (float) {
        return 0.f;
    }
    template<>
    double
    ZeroElement (double) {
        return 0.;
    }
    template<>
    vector<float>
    ZeroElement (vector<float>) {
        return zero_vector<float> ();
    }
    template<>
    vector<double>
    ZeroElement (vector<double>) {
        return zero_vector<double> ();
    }
    template<>
    matrix<float>
    ZeroElement (matrix<float>) {
        return zero_matrix<float> ();
    }
    template<>
    matrix<double>
    ZeroElement (matrix<double>) {
        return zero_matrix<double> ();
    }
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300)
    template<>
    std::complex<float>
    ZeroElement (std::complex<float>) {
        return std::complex<float> (0.f);
    }
    template<>
    std::complex<double>
    ZeroElement (std::complex<double>) {
        return std::complex<double> (0.);
    }
    template<>
    vector<std::complex<float> >
    ZeroElement (vector<std::complex<float> >) {
        return zero_vector<std::complex<float> > ();
    }
    template<>
    vector<std::complex<double> >
    ZeroElement (vector<std::complex<double> >) {
        return zero_vector<std::complex<double> > ();
    }
    template<>
    matrix<std::complex<float> >
    ZeroElement (matrix<std::complex<float> >) {
        return zero_matrix<std::complex<float> > ();
    }
    template<>
    matrix<std::complex<double> >
    ZeroElement (matrix<std::complex<double> >) {
        return zero_matrix<std::complex<double> > ();
    }
#endif

    template<class T>
    T
    OneElement (T);
    template<>
    float
    OneElement (float) {
        return 1.f;
    }
    template<>
    double
    OneElement (double) {
        return 1.;
    }
    template<>
    matrix<float>
    OneElement (matrix<float>) {
        return identity_matrix<float> ();
    }
    template<>
    matrix<double>
    OneElement (matrix<double>) {
        return identity_matrix<double> ();
    }
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300)
    template<>
    std::complex<float>
    OneElement (std::complex<float>) {
        return std::complex<float> (1.f);
    }
    template<>
    std::complex<double>
    OneElement (std::complex<double>) {
        return std::complex<double> (1.);
    }
    template<>
    matrix<std::complex<float> >
    OneElement (matrix<std::complex<float> >) {
        return identity_matrix<std::complex<float> > ();
    }
    template<>
    matrix<std::complex<double> >
    OneElement (matrix<std::complex<double> >) {
        return identity_matrix<std::complex<double> > ();
    }
#endif

    template<class E1, class E2>
    bool
    operator == (const vector_expression<E1> &e1, const vector_expression<E2> &e2) {
        typedef BOOST_UBLAS_TYPENAME promote_traits<BOOST_UBLAS_TYPENAME E1::value_type,
                                                    BOOST_UBLAS_TYPENAME E2::value_type>::promote_type value_type;
        typedef BOOST_UBLAS_TYPENAME type_traits<value_type>::real_type real_type;
        return norm_inf (e1 - e2) == real_type (0);
    }
    template<class E1, class E2>
    bool
    operator == (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
        typedef BOOST_UBLAS_TYPENAME promote_traits<BOOST_UBLAS_TYPENAME E1::value_type,
                                                    BOOST_UBLAS_TYPENAME E2::value_type>::promote_type value_type;
        typedef BOOST_UBLAS_TYPENAME type_traits<value_type>::real_type real_type;
        return norm_inf (e1 - e2) == real_type (0);
    }

    template<class T>
    struct AdditiveAbelianGroupConcept {
        typedef T value_type;

        static void constraints () {
            bool r;
            value_type a = value_type (), b = value_type (), c = value_type ();
            r = (a + b) + c == a + (b + c);
            r = ZeroElement (value_type ()) + a == a;
            r = a + ZeroElement (value_type ()) == a;
            r = a + (- a) == ZeroElement (value_type ());
            r = (- a) + a == ZeroElement (value_type ());
            r = a + b == b + a;
            ignore_unused_variable_warning (r);
        }
    };

    template<class T>
    struct MultiplicativeAbelianGroupConcept {
        typedef T value_type;

        static void constraints () {
            bool r;
            value_type a = value_type (), b = value_type (), c = value_type ();
            r = (a * b) * c == a * (b * c);
            r = OneElement (value_type ()) * a == a;
            r = a * OneElement (value_type ()) == a;
            r = a * (OneElement (value_type ()) / a) == a;
            r = (OneElement (value_type ()) / a) * a == a;
            r = a * b == b * a;
            ignore_unused_variable_warning (r);
        }
    };

    template<class T>
    struct RingWithIdentityConcept {
        typedef T value_type;

        static void constraints () {
            AdditiveAbelianGroupConcept<value_type>::constraints ();
            bool r;
            value_type a = value_type (), b = value_type (), c = value_type ();
            r = (a * b) * c == a * (b * c);
            r = (a + b) * c == a * c + b * c;
            r = OneElement (value_type ()) * a == a;
            r = a * OneElement (value_type ()) == a;
            ignore_unused_variable_warning (r);
        }
        static void constraints (int) {
            AdditiveAbelianGroupConcept<value_type>::constraints ();
            bool r;
            value_type a = value_type (), b = value_type (), c = value_type ();
            r = prod (T (prod (a, b)), c) == prod (a, T (prod (b, c)));
            r = prod (a + b, c) == prod (a, c) + prod (b, c);
            r = prod (OneElement (value_type ()), a) == a;
            r = prod (a, OneElement (value_type ())) == a;
            ignore_unused_variable_warning (r);
        }
    };

    template<class T>
    struct CommutativeRingWithIdentityConcept {
        typedef T value_type;

        static void constraints () {
            RingWithIdentityConcept<value_type>::constraints ();
            bool r;
            value_type a = value_type (), b = value_type ();
            r = a * b == b * a;
            ignore_unused_variable_warning (r);
        }
    };

    template<class T>
    struct FieldConcept {
        typedef T value_type;

        static void constraints () {
            CommutativeRingWithIdentityConcept<value_type>::constraints ();
            bool r;
            value_type a = value_type ();
            r = a == ZeroElement (value_type ()) || a * (OneElement (value_type ()) / a) == a;
            r = a == ZeroElement (value_type ()) || (OneElement (value_type ()) / a) * a == a;
            ignore_unused_variable_warning (r);
        }
    };

    template<class T, class V>
    struct VectorSpaceConcept {
        typedef T value_type;
        typedef V vector_type;

        static void constraints () {
            FieldConcept<value_type>::constraints ();
            AdditiveAbelianGroupConcept<vector_type>::constraints ();
            bool r;
            value_type alpha = value_type (), beta = value_type ();
            vector_type a = vector_type (), b = vector_type ();
            r = alpha * (a + b) == alpha * a + alpha * b;
            r = (alpha + beta) * a == alpha * a + beta * a;
            r = (alpha * beta) * a == alpha * (beta * a);
            r = OneElement (value_type ()) * a == a;
            ignore_unused_variable_warning (r);
        }
    };

    template<class T, class V, class M>
    struct LinearOperatorConcept {
        typedef T value_type;
        typedef V vector_type;
        typedef M matrix_type;

        static void constraints () {
            VectorSpaceConcept<value_type, vector_type>::constraints ();
            bool r;
            value_type alpha = value_type (), beta = value_type ();
            vector_type a = vector_type (), b = vector_type ();
            matrix_type A = matrix_type ();
            r = prod (A, alpha * a + beta * b) == alpha * prod (A, a) + beta * prod (A, b);
            ignore_unused_variable_warning (r);
        }
    };

    void concept_checks () {

        // Storage Array
#if defined (INTERNAL) || defined (INTERNAL_STORAGE)
        StorageArrayConcept<const std::vector<double> >::constraints ();
        MutableStorageArrayConcept<std::vector<double> >::constraints ();
        RandomAccessIteratorConcept<std::vector<double>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<std::vector<double>::iterator>::constraints ();

        StorageArrayConcept<const bounded_array<double, 1> >::constraints ();
        MutableStorageArrayConcept<bounded_array<double, 1> >::constraints ();
        RandomAccessIteratorConcept<bounded_array<double, 1>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<bounded_array<double, 1>::iterator>::constraints ();

        StorageArrayConcept<const unbounded_array<double> >::constraints ();
        MutableStorageArrayConcept<unbounded_array<double> >::constraints ();
        RandomAccessIteratorConcept<unbounded_array<double>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<unbounded_array<double>::iterator>::constraints ();

        StorageArrayConcept<const array_adaptor<double> >::constraints ();
        MutableStorageArrayConcept<array_adaptor<double> >::constraints ();
        RandomAccessIteratorConcept<array_adaptor<double>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<array_adaptor<double>::iterator>::constraints ();

        IndexSetConcept<range>::constraints ();
        RandomAccessIteratorConcept<range::const_iterator>::constraints ();

        IndexSetConcept<slice>::constraints ();
        RandomAccessIteratorConcept<slice::const_iterator>::constraints ();

        IndexSetConcept<indirect_array<> >::constraints ();
        RandomAccessIteratorConcept<indirect_array<>::const_iterator>::constraints ();
#endif

        // Storage Sparse
#if defined (INTERNAL) || defined (INTERNAL_STORAGE_SPARSE)
        StorageSparseConcept<const map_array<std::size_t, double> >::constraints ();
        MutableStorageSparseConcept<map_array<std::size_t, double> >::constraints ();
        RandomAccessIteratorConcept<map_array<std::size_t, double>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<map_array<std::size_t, double>::iterator>::constraints ();

        StorageSparseConcept<const std::map<std::size_t, double> >::constraints ();
        MutableStorageSparseConcept<std::map<std::size_t, double> >::constraints ();
        BidirectionalIteratorConcept<std::map<std::size_t, double>::const_iterator>::constraints ();
                // Not value_type mutable
        BidirectionalIteratorConcept<std::map<std::size_t, double>::iterator>::constraints ();

#ifdef BOOST_UBLAS_DEPRACATED
        StorageSparseConcept<const set_array<std::size_t> >::constraints ();
        MutableStorageSparseConcept<set_array<std::size_t> >::constraints ();
        RandomAccessIteratorConcept<set_array<std::size_t>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<set_array<std::size_t>::iterator>::constraints ();

        StorageSparseConcept<const std::set<std::size_t> >::constraints ();
        MutableStorageSparseConcept<std::set<std::size_t> >::constraints ();
        BidirectionalIteratorConcept<std::set<std::size_t>::const_iterator>::constraints ();
        MutableBidirectionalIteratorConcept<std::set<std::size_t>::iterator>::constraints ();
#endif

#endif

        // Vector
#if defined (INTERNAL) || defined (INTERNAL_VECTOR)
        VectorConcept<const vector<double> >::constraints ();
        MutableVectorConcept<vector<double> >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector<double>::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector<double>::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector<double>::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector<double>::reverse_iterator>::constraints ();

        VectorConcept<unit_vector<double> >::constraints ();
        IndexedRandomAccess1DIteratorConcept<unit_vector<double>::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<unit_vector<double>::const_reverse_iterator>::constraints ();

        VectorConcept<zero_vector<double> >::constraints ();
        IndexedBidirectional1DIteratorConcept<zero_vector<double>::const_iterator>::constraints ();
        IndexedBidirectional1DIteratorConcept<zero_vector<double>::const_reverse_iterator>::constraints ();

        VectorConcept<scalar_vector<double> >::constraints ();
        IndexedRandomAccess1DIteratorConcept<scalar_vector<double>::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<scalar_vector<double>::const_reverse_iterator>::constraints ();

        VectorConcept<const c_vector<double, 1> >::constraints ();
        MutableVectorConcept<c_vector<double, 1> >::constraints ();
        IndexedRandomAccess1DIteratorConcept<c_vector<double, 1>::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<c_vector<double, 1>::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<c_vector<double, 1>::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<c_vector<double, 1>::reverse_iterator>::constraints ();
#endif

        // Vector Proxies
#if defined (INTERNAL) || defined (INTERNAL_VECTOR_PROXY)
        VectorExpressionConcept<const vector_range<const vector<double> > >::constraints ();
        MutableVectorExpressionConcept<vector_range<vector<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_range<vector<double> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_range<vector<double> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_range<vector<double> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_range<vector<double> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const vector_slice<const vector<double> > >::constraints ();
        MutableVectorExpressionConcept<vector_slice<vector<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_slice<vector<double> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_slice<vector<double> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_slice<vector<double> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_slice<vector<double> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const vector_indirect<const vector<double> > >::constraints ();
        MutableVectorExpressionConcept<vector_indirect<vector<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_indirect<vector<double> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_indirect<vector<double> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_indirect<vector<double> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_indirect<vector<double> >::reverse_iterator>::constraints ();
#endif

        // Sparse Vector
#if defined (INTERNAL) || defined (INTERNAL_VECTOR_SPARSE)
        VectorConcept<const sparse_vector<double> >::constraints ();
        MutableVectorConcept<sparse_vector<double> >::constraints ();
        IndexedBidirectional1DIteratorConcept<sparse_vector<double>::const_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<sparse_vector<double>::iterator>::constraints ();
        IndexedBidirectional1DIteratorConcept<sparse_vector<double>::const_reverse_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<sparse_vector<double>::reverse_iterator>::constraints ();

        VectorConcept<const compressed_vector<double> >::constraints ();
        MutableVectorConcept<compressed_vector<double> >::constraints ();
        IndexedBidirectional1DIteratorConcept<compressed_vector<double>::const_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<compressed_vector<double>::iterator>::constraints ();
        IndexedBidirectional1DIteratorConcept<compressed_vector<double>::const_reverse_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<compressed_vector<double>::reverse_iterator>::constraints ();

        VectorConcept<const coordinate_vector<double> >::constraints ();
        MutableVectorConcept<coordinate_vector<double> >::constraints ();
        IndexedBidirectional1DIteratorConcept<coordinate_vector<double>::const_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<coordinate_vector<double>::iterator>::constraints ();
        IndexedBidirectional1DIteratorConcept<coordinate_vector<double>::const_reverse_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<coordinate_vector<double>::reverse_iterator>::constraints ();
#endif

        // Matrix
#if defined (INTERNAL) || defined (INTERNAL_MATRIX)
        MatrixConcept<const matrix<double> >::constraints ();
        MutableMatrixConcept<matrix<double> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix<double>::const_iterator1,
                                             matrix<double>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix<double>::iterator1,
                                                    matrix<double>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix<double>::const_reverse_iterator1,
                                             matrix<double>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix<double>::reverse_iterator1,
                                                    matrix<double>::reverse_iterator2>::constraints ();

        MatrixConcept<const vector_of_vector<double> >::constraints ();
        MutableMatrixConcept<vector_of_vector<double> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<vector_of_vector<double>::const_iterator1,
                                             vector_of_vector<double>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<vector_of_vector<double>::iterator1,
                                                    vector_of_vector<double>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<vector_of_vector<double>::const_reverse_iterator1,
                                             vector_of_vector<double>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<vector_of_vector<double>::reverse_iterator1,
                                                    vector_of_vector<double>::reverse_iterator2>::constraints ();

        MatrixConcept<identity_matrix<double> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<identity_matrix<double>::const_iterator1,
                                             identity_matrix<double>::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<identity_matrix<double>::const_reverse_iterator1,
                                             identity_matrix<double>::const_reverse_iterator2>::constraints ();

        MatrixConcept<zero_matrix<double> >::constraints ();
        IndexedBidirectional2DIteratorConcept<zero_matrix<double>::const_iterator1,
                                              zero_matrix<double>::const_iterator2>::constraints ();
        IndexedBidirectional2DIteratorConcept<zero_matrix<double>::const_reverse_iterator1,
                                              zero_matrix<double>::const_reverse_iterator2>::constraints ();

        MatrixConcept<scalar_matrix<double> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<scalar_matrix<double>::const_iterator1,
                                             scalar_matrix<double>::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<scalar_matrix<double>::const_reverse_iterator1,
                                             scalar_matrix<double>::const_reverse_iterator2>::constraints ();

        MatrixConcept<const c_matrix<double, 1, 1> >::constraints ();
        MutableMatrixConcept<c_matrix<double, 1, 1> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<c_matrix<double, 1, 1>::const_iterator1,
                                             c_matrix<double, 1, 1>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<c_matrix<double, 1, 1>::iterator1,
                                                    c_matrix<double, 1, 1>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<c_matrix<double, 1, 1>::const_reverse_iterator1,
                                             c_matrix<double, 1, 1>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<c_matrix<double, 1, 1>::reverse_iterator1,
                                                    c_matrix<double, 1, 1>::reverse_iterator2>::constraints ();
#endif

        // Matrix Proxies
#if defined (INTERNAL) || defined (INTERNAL_MATRIX_PROXY)
        VectorExpressionConcept<const matrix_row<const matrix<double> > >::constraints ();
        MutableVectorExpressionConcept<matrix_row<matrix<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_row<matrix<double> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_row<matrix<double> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_row<matrix<double> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_row<matrix<double> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const matrix_column<const matrix<double> > >::constraints ();
        MutableVectorExpressionConcept<matrix_column<matrix<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_column<matrix<double> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_column<matrix<double> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_column<matrix<double> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_column<matrix<double> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const matrix_vector_range<const matrix<double> > >::constraints ();
        MutableVectorExpressionConcept<matrix_vector_range<matrix<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_range<matrix<double> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_range<matrix<double> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_range<matrix<double> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_range<matrix<double> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const matrix_vector_slice<const matrix<double> > >::constraints ();
        MutableVectorExpressionConcept<matrix_vector_slice<matrix<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_slice<matrix<double> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_slice<matrix<double> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_slice<matrix<double> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_slice<matrix<double> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const matrix_vector_indirect<const matrix<double>, vector<unsigned> > >::constraints ();
        MutableVectorExpressionConcept<matrix_vector_indirect<matrix<double>, vector<unsigned> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_indirect<matrix<double>, vector<unsigned> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_indirect<matrix<double>, vector<unsigned> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_indirect<matrix<double>, vector<unsigned> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_indirect<matrix<double>, vector<unsigned> >::reverse_iterator>::constraints ();

        MatrixExpressionConcept<const matrix_range<const matrix<double> > >::constraints ();
        MutableMatrixExpressionConcept<matrix_range<matrix<double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_range<matrix<double> >::const_iterator1,
                                             matrix_range<matrix<double> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_range<matrix<double> >::iterator1,
                                                    matrix_range<matrix<double> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_range<matrix<double> >::const_reverse_iterator1,
                                             matrix_range<matrix<double> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_range<matrix<double> >::reverse_iterator1,
                                                    matrix_range<matrix<double> >::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const matrix_slice<const matrix<double> > >::constraints ();
        MutableMatrixExpressionConcept<matrix_slice<matrix<double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_slice<matrix<double> >::const_iterator1,
                                             matrix_slice<matrix<double> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_slice<matrix<double> >::iterator1,
                                                    matrix_slice<matrix<double> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_slice<matrix<double> >::const_reverse_iterator1,
                                             matrix_slice<matrix<double> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_slice<matrix<double> >::reverse_iterator1,
                                                    matrix_slice<matrix<double> >::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const matrix_indirect<const matrix<double> > >::constraints ();
        MutableMatrixExpressionConcept<matrix_indirect<matrix<double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_indirect<matrix<double> >::const_iterator1,
                                             matrix_indirect<matrix<double> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_indirect<matrix<double> >::iterator1,
                                                    matrix_indirect<matrix<double> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_indirect<matrix<double> >::const_reverse_iterator1,
                                             matrix_indirect<matrix<double> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_indirect<matrix<double> >::reverse_iterator1,
                                                    matrix_indirect<matrix<double> >::reverse_iterator2>::constraints ();
#endif

        // Banded Matrix
#if defined (INTERNAL) || defined (INTERNAL_BANDED)
        MatrixConcept<const banded_matrix<double> >::constraints ();
        MutableMatrixConcept<banded_matrix<double> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<banded_matrix<double>::const_iterator1,
                                             banded_matrix<double>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<banded_matrix<double>::iterator1,
                                                    banded_matrix<double>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<banded_matrix<double>::const_reverse_iterator1,
                                             banded_matrix<double>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<banded_matrix<double>::reverse_iterator1,
                                                    banded_matrix<double>::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const banded_adaptor<const matrix<double> > >::constraints ();
        MutableMatrixExpressionConcept<banded_adaptor<matrix<double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<banded_adaptor<matrix<double> >::const_iterator1,
                                             banded_adaptor<matrix<double> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<banded_adaptor<matrix<double> >::iterator1,
                                                    banded_adaptor<matrix<double> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<banded_adaptor<matrix<double> >::const_reverse_iterator1,
                                             banded_adaptor<matrix<double> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<banded_adaptor<matrix<double> >::reverse_iterator1,
                                                    banded_adaptor<matrix<double> >::reverse_iterator2>::constraints ();
#endif

        // Triangular Matrix
#if defined (INTERNAL) || defined (INTERNAL_TRIANGULAR)
        MatrixConcept<const triangular_matrix<double> >::constraints ();
        MutableMatrixConcept<triangular_matrix<double> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<triangular_matrix<double>::const_iterator1,
                                             triangular_matrix<double>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<triangular_matrix<double>::iterator1,
                                                    triangular_matrix<double>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<triangular_matrix<double>::const_reverse_iterator1,
                                             triangular_matrix<double>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<triangular_matrix<double>::reverse_iterator1,
                                                    triangular_matrix<double>::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const triangular_adaptor<const matrix<double> > >::constraints ();
        MutableMatrixExpressionConcept<triangular_adaptor<matrix<double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<triangular_adaptor<matrix<double> >::const_iterator1,
                                             triangular_adaptor<matrix<double> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<triangular_adaptor<matrix<double> >::iterator1,
                                                    triangular_adaptor<matrix<double> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<triangular_adaptor<matrix<double> >::const_reverse_iterator1,
                                             triangular_adaptor<matrix<double> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<triangular_adaptor<matrix<double> >::reverse_iterator1,
                                                    triangular_adaptor<matrix<double> >::reverse_iterator2>::constraints ();
#endif

        // Symmetric Matrix
#if defined (INTERNAL) || defined (INTERNAL_SYMMETRIC)
        MatrixConcept<const symmetric_matrix<double> >::constraints ();
        MutableMatrixConcept<symmetric_matrix<double> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<symmetric_matrix<double>::const_iterator1,
                                             symmetric_matrix<double>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<symmetric_matrix<double>::iterator1,
                                                    symmetric_matrix<double>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<symmetric_matrix<double>::const_reverse_iterator1,
                                             symmetric_matrix<double>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<symmetric_matrix<double>::reverse_iterator1,
                                                    symmetric_matrix<double>::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const symmetric_adaptor<const matrix<double> > >::constraints ();
#ifndef SKIP_BAD
        // const_iterator (iterator) constructor is bad
        MutableMatrixExpressionConcept<symmetric_adaptor<matrix<double> > >::constraints ();
#endif
        IndexedRandomAccess2DIteratorConcept<symmetric_adaptor<matrix<double> >::const_iterator1,
                                             symmetric_adaptor<matrix<double> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<symmetric_adaptor<matrix<double> >::iterator1,
                                                    symmetric_adaptor<matrix<double> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<symmetric_adaptor<matrix<double> >::const_reverse_iterator1,
                                             symmetric_adaptor<matrix<double> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<symmetric_adaptor<matrix<double> >::reverse_iterator1,
                                                    symmetric_adaptor<matrix<double> >::reverse_iterator2>::constraints ();
#endif

        // Hermitian Matrix
#if defined (INTERNAL) || defined (INTERNAL_HERMITIAN)
        MatrixConcept<const hermitian_matrix<double> >::constraints ();
        MutableMatrixConcept<hermitian_matrix<double> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<hermitian_matrix<double>::const_iterator1,
                                             hermitian_matrix<double>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<hermitian_matrix<double>::iterator1,
                                                    hermitian_matrix<double>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<hermitian_matrix<double>::const_reverse_iterator1,
                                             hermitian_matrix<double>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<hermitian_matrix<double>::reverse_iterator1,
                                                    hermitian_matrix<double>::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const hermitian_adaptor<const matrix<double> > >::constraints ();
#ifndef SKIP_BAD
        // const_iterator (iterator) constructor is bad
        MutableMatrixExpressionConcept<hermitian_adaptor<matrix<double> > >::constraints ();
#endif
        IndexedRandomAccess2DIteratorConcept<hermitian_adaptor<matrix<double> >::const_iterator1,
                                             hermitian_adaptor<matrix<double> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<hermitian_adaptor<matrix<double> >::iterator1,
                                                    hermitian_adaptor<matrix<double> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<hermitian_adaptor<matrix<double> >::const_reverse_iterator1,
                                             hermitian_adaptor<matrix<double> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<hermitian_adaptor<matrix<double> >::reverse_iterator1,
                                                    hermitian_adaptor<matrix<double> >::reverse_iterator2>::constraints ();
#endif

        // Sparse Matrix
#if defined (INTERNAL) || defined (INTERNAL_MATRIX_SPARSE)
        MatrixConcept<const sparse_matrix<double> >::constraints ();
        MutableMatrixConcept<sparse_matrix<double> >::constraints ();
        IndexedBidirectional2DIteratorConcept<sparse_matrix<double>::const_iterator1,
                                              sparse_matrix<double>::const_iterator2>::constraints ();
        MutableIndexedBidirectional2DIteratorConcept<sparse_matrix<double>::iterator1,
                                                     sparse_matrix<double>::iterator2>::constraints ();
        IndexedBidirectional2DIteratorConcept<sparse_matrix<double>::const_reverse_iterator1,
                                              sparse_matrix<double>::const_reverse_iterator2>::constraints ();
        MutableIndexedBidirectional2DIteratorConcept<sparse_matrix<double>::reverse_iterator1,
                                                     sparse_matrix<double>::reverse_iterator2>::constraints ();

        MatrixConcept<const sparse_vector_of_sparse_vector<double> >::constraints ();
        MutableMatrixConcept<sparse_vector_of_sparse_vector<double> >::constraints ();
        IndexedBidirectional2DIteratorConcept<sparse_vector_of_sparse_vector<double>::const_iterator1,
                                              sparse_vector_of_sparse_vector<double>::const_iterator2>::constraints ();
        MutableIndexedBidirectional2DIteratorConcept<sparse_vector_of_sparse_vector<double>::iterator1,
                                                     sparse_vector_of_sparse_vector<double>::iterator2>::constraints ();
        IndexedBidirectional2DIteratorConcept<sparse_vector_of_sparse_vector<double>::const_reverse_iterator1,
                                              sparse_vector_of_sparse_vector<double>::const_reverse_iterator2>::constraints ();
        MutableIndexedBidirectional2DIteratorConcept<sparse_vector_of_sparse_vector<double>::reverse_iterator1,
                                                     sparse_vector_of_sparse_vector<double>::reverse_iterator2>::constraints ();

        MatrixConcept<const compressed_matrix<double> >::constraints ();
        MutableMatrixConcept<compressed_matrix<double> >::constraints ();
        IndexedBidirectional2DIteratorConcept<compressed_matrix<double>::const_iterator1,
                                              compressed_matrix<double>::const_iterator2>::constraints ();
        MutableIndexedBidirectional2DIteratorConcept<compressed_matrix<double>::iterator1,
                                                     compressed_matrix<double>::iterator2>::constraints ();
        IndexedBidirectional2DIteratorConcept<compressed_matrix<double>::const_reverse_iterator1,
                                              compressed_matrix<double>::const_reverse_iterator2>::constraints ();
        MutableIndexedBidirectional2DIteratorConcept<compressed_matrix<double>::reverse_iterator1,
                                                     compressed_matrix<double>::reverse_iterator2>::constraints ();

        MatrixConcept<const coordinate_matrix<double> >::constraints ();
        MutableMatrixConcept<coordinate_matrix<double> >::constraints ();
        IndexedBidirectional2DIteratorConcept<coordinate_matrix<double>::const_iterator1,
                                              coordinate_matrix<double>::const_iterator2>::constraints ();
        MutableIndexedBidirectional2DIteratorConcept<coordinate_matrix<double>::iterator1,
                                                     coordinate_matrix<double>::iterator2>::constraints ();
        IndexedBidirectional2DIteratorConcept<coordinate_matrix<double>::const_reverse_iterator1,
                                              coordinate_matrix<double>::const_reverse_iterator2>::constraints ();
        MutableIndexedBidirectional2DIteratorConcept<coordinate_matrix<double>::reverse_iterator1,
                                                     coordinate_matrix<double>::reverse_iterator2>::constraints ();
#endif

        // Scalar Expressions
#if defined (INTERNAL) || defined (INTERNAL_VECTOR_EXPRESSION)
        ScalarExpressionConcept<scalar_value<double > >::constraints ();
        ScalarExpressionConcept<scalar_reference<double > >::constraints ();

        // Vector Expressions
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        VectorExpressionConcept<vector_const_reference<vector<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_const_reference<vector<double> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_const_reference<vector<double> >::const_reverse_iterator>::constraints ();
#endif

        VectorExpressionConcept<vector_reference<vector<double> > >::constraints ();
        MutableVectorExpressionConcept<vector_reference<vector<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_reference<vector<double> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_reference<vector<double> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_reference<vector<double> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_reference<vector<double> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_unary<vector<double>, scalar_identity<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_unary<vector<double>, scalar_identity<double>  >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_unary<vector<double>, scalar_identity<double>  >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_binary<vector<double>, vector<double>, scalar_plus<double, double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary<vector<double>, vector<double>, scalar_plus<double, double> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary<vector<double>, vector<double>, scalar_plus<double, double> >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_binary_scalar1<double, vector<double>, scalar_multiplies<double, double>, scalar_reference<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar1<double, vector<double>, scalar_multiplies<double, double>, scalar_reference<double> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar1<double, vector<double>, scalar_multiplies<double, double>, scalar_reference<double> >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_binary_scalar2<vector<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar2<vector<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<double> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar2<vector<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<double> >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_binary_scalar1<scalar_value<double>, vector<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar1<scalar_value<double>, vector<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar1<scalar_value<double>, vector<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_binary_scalar2<vector<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar2<vector<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar2<vector<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_reverse_iterator>::constraints ();

        ScalarExpressionConcept<vector_scalar_unary<vector<double>, vector_sum<double> > >::constraints ();
        ScalarExpressionConcept<vector_scalar_unary<vector<double>, vector_norm_1<double> > >::constraints ();
        ScalarExpressionConcept<vector_scalar_unary<vector<double>, vector_norm_2<double> > >::constraints ();
        ScalarExpressionConcept<vector_scalar_unary<vector<double>, vector_norm_inf<double> > >::constraints ();

        ScalarExpressionConcept<vector_scalar_binary<vector<double>, vector<double>, vector_inner_prod<double, double, double> > >::constraints ();
#endif

        // Matrix Expressions
#if defined (INTERNAL) || defined (INTERNAL_MATRIX_EXPRESSION)
        MatrixExpressionConcept<matrix_reference<matrix<double> > >::constraints ();
        MutableMatrixExpressionConcept<matrix_reference<matrix<double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_reference<matrix<double> >::const_iterator1,
                                             matrix_reference<matrix<double> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_reference<matrix<double> >::iterator1,
                                                    matrix_reference<matrix<double> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_reference<matrix<double> >::const_reverse_iterator1,
                                             matrix_reference<matrix<double> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_reference<matrix<double> >::reverse_iterator1,
                                                    matrix_reference<matrix<double> >::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<vector_matrix_binary<vector<double>, vector<double>, scalar_multiplies<double, double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<vector_matrix_binary<vector<double>, vector<double>, scalar_multiplies<double, double> >::const_iterator1,
                                             vector_matrix_binary<vector<double>, vector<double>, scalar_multiplies<double, double> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<vector_matrix_binary<vector<double>, vector<double>, scalar_multiplies<double, double> >::const_reverse_iterator1,
                                             vector_matrix_binary<vector<double>, vector<double>, scalar_multiplies<double, double> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_unary1<matrix<double>, scalar_identity<double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_unary1<matrix<double>, scalar_identity<double> >::const_iterator1,
                                             matrix_unary1<matrix<double>, scalar_identity<double> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_unary1<matrix<double>, scalar_identity<double> >::const_reverse_iterator1,
                                             matrix_unary1<matrix<double>, scalar_identity<double> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_unary2<matrix<double>, scalar_identity<double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_unary2<matrix<double>, scalar_identity<double> >::const_iterator1,
                                             matrix_unary2<matrix<double>, scalar_identity<double> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_unary2<matrix<double>, scalar_identity<double> >::const_reverse_iterator1,
                                             matrix_unary2<matrix<double>, scalar_identity<double> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_binary<matrix<double>, matrix<double>, scalar_plus<double, double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary<matrix<double>, matrix<double>, scalar_plus<double, double> >::const_iterator1,
                                             matrix_binary<matrix<double>, matrix<double>, scalar_plus<double, double> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary<matrix<double>, matrix<double>, scalar_plus<double, double> >::const_reverse_iterator1,
                                             matrix_binary<matrix<double>, matrix<double>, scalar_plus<double, double> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_binary_scalar1<double, matrix<double>, scalar_multiplies<double, double>, scalar_reference<double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar1<double, matrix<double>, scalar_multiplies<double, double>, scalar_reference<double> >::const_iterator1,
                                             matrix_binary_scalar1<double, matrix<double>, scalar_multiplies<double, double>, scalar_reference<double> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar1<double, matrix<double>, scalar_multiplies<double, double>, scalar_reference<double> >::const_reverse_iterator1,
                                             matrix_binary_scalar1<double, matrix<double>, scalar_multiplies<double, double>, scalar_reference<double> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_binary_scalar2<matrix<double>, double, scalar_multiplies<double, double>, scalar_reference<double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar2<matrix<double>, double, scalar_multiplies<double, double>, scalar_reference<double> >::const_iterator1,
                                             matrix_binary_scalar2<matrix<double>, double, scalar_multiplies<double, double>, scalar_reference<double> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar2<matrix<double>, double, scalar_multiplies<double, double>, scalar_reference<double> >::const_reverse_iterator1,
                                             matrix_binary_scalar2<matrix<double>, double, scalar_multiplies<double, double>, scalar_reference<double> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_binary_scalar1<scalar_value<double>, matrix<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar1<scalar_value<double>, matrix<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_iterator1,
                                             matrix_binary_scalar1<scalar_value<double>, matrix<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar1<scalar_value<double>, matrix<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_reverse_iterator1,
                                             matrix_binary_scalar1<scalar_value<double>, matrix<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_binary_scalar2<matrix<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar2<matrix<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_iterator1,
                                             matrix_binary_scalar2<matrix<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar2<matrix<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_reverse_iterator1,
                                             matrix_binary_scalar2<matrix<double>, scalar_value<double>, scalar_multiplies<double, double>, scalar_reference<scalar_value<double> > >::const_reverse_iterator2>::constraints ();

        VectorExpressionConcept<matrix_vector_binary1<matrix<double>, vector<double>, matrix_vector_prod1<double, double, double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_binary1<matrix<double>, vector<double>, matrix_vector_prod1<double, double, double> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_binary1<matrix<double>, vector<double>, matrix_vector_prod1<double, double, double> >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<matrix_vector_binary2<vector<double>, matrix<double>, matrix_vector_prod2<double, double, double> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_binary2<vector<double>, matrix<double>, matrix_vector_prod2<double, double, double> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_binary2<vector<double>, matrix<double>, matrix_vector_prod2<double, double, double> >::const_reverse_iterator>::constraints ();

        MatrixExpressionConcept<matrix_matrix_binary<matrix<double>, matrix<double>, matrix_matrix_prod<double, double, double> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_matrix_binary<matrix<double>, matrix<double>, matrix_matrix_prod<double, double, double> >::const_iterator1,
                                             matrix_matrix_binary<matrix<double>, matrix<double>, matrix_matrix_prod<double, double, double> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_matrix_binary<matrix<double>, matrix<double>, matrix_matrix_prod<double, double, double> >::const_reverse_iterator1,
                                             matrix_matrix_binary<matrix<double>, matrix<double>, matrix_matrix_prod<double, double, double> >::const_reverse_iterator2>::constraints ();

        ScalarExpressionConcept<matrix_scalar_unary<matrix<double>, matrix_norm_1<double> > >::constraints ();
        ScalarExpressionConcept<matrix_scalar_unary<matrix<double>, matrix_norm_frobenius<double> > >::constraints ();
        ScalarExpressionConcept<matrix_scalar_unary<matrix<double>, matrix_norm_inf<double> > >::constraints ();
#endif

#ifdef EXTERNAL
        AdditiveAbelianGroupConcept<float>::constraints ();
        CommutativeRingWithIdentityConcept<float>::constraints ();
        FieldConcept<float>::constraints ();
        VectorSpaceConcept<float, vector<float> >::constraints ();
        RingWithIdentityConcept<matrix<float> >::constraints (0);
        VectorSpaceConcept<float, matrix<float> >::constraints ();
        LinearOperatorConcept<float, vector<float>, matrix<float> >::constraints ();

        AdditiveAbelianGroupConcept<double>::constraints ();
        CommutativeRingWithIdentityConcept<double>::constraints ();
        FieldConcept<double>::constraints ();
        VectorSpaceConcept<double, vector<double> >::constraints ();
        RingWithIdentityConcept<matrix<double> >::constraints (0);
        VectorSpaceConcept<double, matrix<double> >::constraints ();
        LinearOperatorConcept<double, vector<double>, matrix<double> >::constraints ();

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300)
        AdditiveAbelianGroupConcept<std::complex<float> >::constraints ();
        CommutativeRingWithIdentityConcept<std::complex<float> >::constraints ();
        FieldConcept<std::complex<float> >::constraints ();
        VectorSpaceConcept<std::complex<float>, vector<std::complex<float> > >::constraints ();
        RingWithIdentityConcept<matrix<std::complex<float> > >::constraints (0);
        VectorSpaceConcept<std::complex<float>, matrix<std::complex<float> > >::constraints ();
        LinearOperatorConcept<std::complex<float>, vector<std::complex<float> >, matrix<std::complex<float> > >::constraints ();

        AdditiveAbelianGroupConcept<std::complex<double> >::constraints ();
        CommutativeRingWithIdentityConcept<std::complex<double> >::constraints ();
        FieldConcept<std::complex<double> >::constraints ();
        VectorSpaceConcept<std::complex<double>, vector<std::complex<double> > >::constraints ();
        RingWithIdentityConcept<matrix<std::complex<double> > >::constraints (0);
        VectorSpaceConcept<std::complex<double>, matrix<std::complex<double> > >::constraints ();
        LinearOperatorConcept<std::complex<double>, vector<std::complex<double> >, matrix<std::complex<double> > >::constraints ();
#endif
#endif
    }

}}}

#endif
