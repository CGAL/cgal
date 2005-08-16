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

#ifndef _BOOST_UBLAS_CONCEPTS_
#define _BOOST_UBLAS_CONCEPTS_

#include <boost/concept_check.hpp>

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
        typedef I1 subiterator1_type;
        typedef I2 subiterator2_type;

        static void constraints () {
            BidirectionalIteratorConcept<subiterator1_type>::constraints ();
            BidirectionalIteratorConcept<subiterator2_type>::constraints ();
            Indexed2DIteratorConcept<subiterator1_type>::constraints ();
            Indexed2DIteratorConcept<subiterator2_type>::constraints ();
        }
    };

    template<class I1, class I2>
    struct MutableIndexedBidirectional2DIteratorConcept {
        typedef I1 subiterator1_type;
        typedef I2 subiterator2_type;

        static void constraints () {
            MutableBidirectionalIteratorConcept<subiterator1_type>::constraints ();
            MutableBidirectionalIteratorConcept<subiterator2_type>::constraints ();
            Indexed2DIteratorConcept<subiterator1_type>::constraints ();
            Indexed2DIteratorConcept<subiterator2_type>::constraints ();
        }
    };

    template<class I1, class I2>
    struct IndexedRandomAccess2DIteratorConcept {
        typedef I1 subiterator1_type;
        typedef I2 subiterator2_type;

        static void constraints () {
            RandomAccessIteratorConcept<subiterator1_type>::constraints ();
            RandomAccessIteratorConcept<subiterator2_type>::constraints ();
            Indexed2DIteratorConcept<subiterator1_type>::constraints ();
            Indexed2DIteratorConcept<subiterator2_type>::constraints ();
        }
    };

    template<class I1, class I2>
    struct MutableIndexedRandomAccess2DIteratorConcept {
        typedef I1 subiterator1_type;
        typedef I2 subiterator2_type;

        static void constraints () {
            MutableRandomAccessIteratorConcept<subiterator1_type>::constraints ();
            MutableRandomAccessIteratorConcept<subiterator2_type>::constraints ();
            Indexed2DIteratorConcept<subiterator1_type>::constraints ();
            Indexed2DIteratorConcept<subiterator2_type>::constraints ();
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
            AssignableConcept<generator_type>::constraints (generator_type ());
            ReversibleContainerConcept<generator_type>::constraints ();
            generator_type g = generator_type ();
            size_type n (0);
            value_type t;
            // Element access
            t = g (n);
            ignore_unused_variable_warning (t);
        }
    };

    template<class SE>
    struct ScalarExpressionConcept {
        typedef SE scalar_expression_type;
        typedef typename SE::value_type value_type;

        static void constraints () {
                scalar_expression_type *sp;
            scalar_expression_type s = *sp;
            value_type t;
            // Conversion
            t = s;
            ignore_unused_variable_warning (t);
        }
    };

    template<class VE>
    struct VectorExpressionConcept {
        typedef VE vector_expression_type;
        typedef typename VE::type_category type_category;
        typedef typename VE::size_type size_type;
        typedef typename VE::value_type value_type;
        typedef typename VE::const_iterator const_iterator_type;
        typedef typename VE::const_reverse_iterator const_reverse_iterator_type;

        static void constraints () {
                vector_expression_type *vp;
                const vector_expression_type *cvp;
            vector_expression_type v = *vp;
            const vector_expression_type cv = *cvp;
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

    template<class VE>
    struct MutableVectorExpressionConcept {
        typedef VE vector_expression_type;
        typedef typename VE::size_type size_type;
        typedef typename VE::value_type value_type;
        typedef typename VE::iterator iterator_type;
        typedef typename VE::reverse_iterator reverse_iterator_type;

        static void constraints () {
                vector_expression_type *vp;
            AssignableConcept<vector_expression_type>::constraints (*vp);
            VectorExpressionConcept<vector_expression_type>::constraints ();
            vector_expression_type v = *vp, v1 = *vp, v2 = *vp;
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

    template<class ME>
    struct MatrixExpressionConcept {
        typedef ME matrix_expression_type;
        typedef typename ME::type_category type_category;
        typedef typename ME::size_type size_type;
        typedef typename ME::value_type value_type;
        typedef typename ME::const_iterator1 const_subiterator1_type;
        typedef typename ME::const_iterator2 const_subiterator2_type;
        typedef typename ME::const_reverse_iterator1 const_reverse_subiterator1_type;
        typedef typename ME::const_reverse_iterator2 const_reverse_subiterator2_type;

        static void constraints () {
                matrix_expression_type *mp;
                const matrix_expression_type *cmp;
            matrix_expression_type m = *mp;
            const matrix_expression_type cm = *cmp;
            size_type n (0), i (0), j (0);
            value_type t;
            // Find (internal?)
            const_subiterator1_type cit1 (m.find1 (0, i, j));
            const_subiterator2_type cit2 (m.find2 (0, i, j));
            // Beginning of range
            const_subiterator1_type cit1_begin (m.begin1 ());
            const_subiterator2_type cit2_begin (m.begin2 ());
            // End of range
            const_subiterator1_type cit1_end (m.end1 ());
            const_subiterator2_type cit2_end (m.end2 ());
            // Size
            n = m.size1 ();
            n = m.size2 ();
            // Beginning of reverse range
            const_reverse_subiterator1_type crit1_begin (cm.rbegin1 ());
            const_reverse_subiterator2_type crit2_begin (cm.rbegin2 ());
            // End of reverse range
            const_reverse_subiterator1_type crit1_end (cm.rend1 ());
            const_reverse_subiterator2_type crit2_end (cm.rend2 ());
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

    template<class ME>
    struct MutableMatrixExpressionConcept {
        typedef ME matrix_expression_type;
        typedef typename ME::size_type size_type;
        typedef typename ME::value_type value_type;
        typedef typename ME::iterator1 subiterator1_type;
        typedef typename ME::iterator2 subiterator2_type;
        typedef typename ME::reverse_iterator1 reverse_subiterator1_type;
        typedef typename ME::reverse_iterator2 reverse_subiterator2_type;

        static void constraints () {
                matrix_expression_type *mp;
            AssignableConcept<matrix_expression_type>::constraints (*mp);
            MatrixExpressionConcept<matrix_expression_type>::constraints ();
            matrix_expression_type m = *mp, m1 = *mp, m2 = *mp;
            size_type i (0), j (0);
            value_type t = value_type ();
            // Find (internal?)
            subiterator1_type it1 (m.find1 (0, i, j));
            subiterator2_type it2 (m.find2 (0, i, j));
            // Beginning of range
            subiterator1_type it1_begin (m.begin1 ());
            subiterator2_type it2_begin (m.begin2 ());
            // End of range
            subiterator1_type it1_end (m.end1 ());
            subiterator2_type it2_end (m.end2 ());
            // Swap
            m1.swap (m2);
            // Beginning of reverse range
            reverse_subiterator1_type rit1_begin (m.rbegin1 ());
            reverse_subiterator2_type rit2_begin (m.rbegin2 ());
            // End of reverse range
            reverse_subiterator1_type rit1_end (m.rend1 ());
            reverse_subiterator2_type rit2_end (m.rend2 ());
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
        typedef typename V::value_type value_type;
        typedef const value_type *const_pointer;

        static void constraints () {
            VectorExpressionConcept<vector_type>::constraints ();
            size_type n (0);
            size_type i (0);
            // Sizing constructor
            vector_type v (n);
            // Element support
            const_pointer p = v.find_element (i);

            ignore_unused_variable_warning (p);
        }
    };

    template<class V>
    struct MutableVectorConcept {
        typedef V vector_type;
        typedef typename V::size_type size_type;
        typedef typename V::value_type value_type;
        typedef value_type *pointer;

        static void constraints () {
            VectorConcept<vector_type>::constraints ();
            MutableVectorExpressionConcept<vector_type>::constraints ();
            size_type n (0);
            value_type t = value_type ();
            size_type i (0);
            // Sizing constructor
            vector_type v (n);
            // Element support
            pointer p = v.find_element (i);
            // Element assignment
            value_type r = v.insert_element (i, t);
            v.insert_element (i, t) = r;
            // Zeroing
            v.clear ();
            // Resize
            v.resize (n);

            ignore_unused_variable_warning (p);
            ignore_unused_variable_warning (r);
        }
    };

    template<class V>
    struct SparseVectorConcept {
        typedef V vector_type;
        typedef typename V::size_type size_type;

        static void constraints () {
            VectorConcept<vector_type>::constraints ();
        }
    };

    template<class V>
    struct MutableSparseVectorConcept {
        typedef V vector_type;
        typedef typename V::size_type size_type;
        typedef typename V::value_type value_type;

        static void constraints () {
            SparseVectorConcept<vector_type>::constraints ();
            MutableVectorConcept<vector_type>::constraints ();
            size_type n (0);
            size_type i (0);
            // Sizing constructor
            vector_type v (n);
            // Element erasure
            v.erase_element (i);
        }
    };

    template<class M>
    struct MatrixConcept {
        typedef M matrix_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;
        typedef const value_type *const_pointer;

        static void constraints () {
            MatrixExpressionConcept<matrix_type>::constraints ();
            size_type n (0);
            size_type i (0), j (0);
            // Sizing constructor
            matrix_type m (n, n);
            // Element support
#ifndef SKIP_BAD
            const_pointer p = m.find_element (i, j);
#else
            const_pointer p;
            ignore_unused_variable_warning (i);
            ignore_unused_variable_warning (j);
#endif
            ignore_unused_variable_warning (p);
        }
    };

    template<class M>
    struct MutableMatrixConcept {
        typedef M matrix_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;
        typedef value_type *pointer;

        static void constraints () {
            MatrixConcept<matrix_type>::constraints ();
            MutableMatrixExpressionConcept<matrix_type>::constraints ();
            size_type n (0);
            value_type t = value_type ();
            size_type i (0), j (0);
            // Sizing constructor
            matrix_type m (n, n);
            // Element support
#ifndef SKIP_BAD
            pointer p = m.find_element (i, j);
            ignore_unused_variable_warning (i);
            ignore_unused_variable_warning (j);
#else
            pointer p;
#endif
            // Element assigment
            value_type r = m.insert_element (i, j, t);
            m.insert_element (i, j, t) = r;
            // Zeroing
            m.clear ();
            // Resize
            m.resize (n, n);
            m.resize (n, n, false);

            ignore_unused_variable_warning (p);
            ignore_unused_variable_warning (r);
        }
    };

    template<class M>
    struct SparseMatrixConcept {
        typedef M matrix_type;
        typedef typename M::size_type size_type;

        static void constraints () {
            MatrixConcept<matrix_type>::constraints ();
        }
    };

    template<class M>
    struct MutableSparseMatrixConcept {
        typedef M matrix_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;

        static void constraints () {
            SparseMatrixConcept<matrix_type>::constraints ();
            MutableMatrixConcept<matrix_type>::constraints ();
            size_type n (0);
            size_type i (0), j (0);
            // Sizing constructor
            matrix_type m (n, n);
            // Elemnent erasure
            m.erase_element (i, j);
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

    template<class E1, class E2>
    bool
    operator == (const vector_expression<E1> &e1, const vector_expression<E2> &e2) {
        typedef typename promote_traits<typename E1::value_type,
                                                    typename E2::value_type>::promote_type value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        return norm_inf (e1 - e2) == real_type/*zero*/();
    }
    template<class E1, class E2>
    bool
    operator == (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
        typedef typename promote_traits<typename E1::value_type,
                                                    typename E2::value_type>::promote_type value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        return norm_inf (e1 - e2) == real_type/*zero*/();
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

        // Allow tests to be group to keep down compiler storage requirement
#ifdef INTERAL
#define INTERNAL_STORAGE
#define INTERNAL_VECTOR
#define INTERNAL_MATRIX
#define INTERNAL_SPECIAL
#define INTERNAL_SPARSE
#define INTERNAL_EXPRESSION
#endif

        // Element value type for tests
        typedef float T;

        // Storage Array
#if defined (INTERNAL_STORAGE) || defined (INTERNAL_STORAGE_DENSE)
        StorageArrayConcept<const std::vector<T> >::constraints ();
        MutableStorageArrayConcept<std::vector<T> >::constraints ();
        RandomAccessIteratorConcept<std::vector<T>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<std::vector<T>::iterator>::constraints ();

        StorageArrayConcept<const bounded_array<T, 1> >::constraints ();
        MutableStorageArrayConcept<bounded_array<T, 1> >::constraints ();
        RandomAccessIteratorConcept<bounded_array<T, 1>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<bounded_array<T, 1>::iterator>::constraints ();

        StorageArrayConcept<const unbounded_array<T> >::constraints ();
        MutableStorageArrayConcept<unbounded_array<T> >::constraints ();
        RandomAccessIteratorConcept<unbounded_array<T>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<unbounded_array<T>::iterator>::constraints ();

        StorageArrayConcept<const array_adaptor<T> >::constraints ();
        MutableStorageArrayConcept<array_adaptor<T> >::constraints ();
        RandomAccessIteratorConcept<array_adaptor<T>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<array_adaptor<T>::iterator>::constraints ();

        IndexSetConcept<range>::constraints ();
        RandomAccessIteratorConcept<range::const_iterator>::constraints ();

        IndexSetConcept<slice>::constraints ();
        RandomAccessIteratorConcept<slice::const_iterator>::constraints ();

        IndexSetConcept<indirect_array<> >::constraints ();
        RandomAccessIteratorConcept<indirect_array<>::const_iterator>::constraints ();
#endif

        // Storage Sparse
#if defined (INTERNAL_STORAGE) || defined (INTERNAL_STORAGE_SPARSE)
        StorageSparseConcept<const map_array<std::size_t, T> >::constraints ();
        MutableStorageSparseConcept<map_array<std::size_t, T> >::constraints ();
        RandomAccessIteratorConcept<map_array<std::size_t, T>::const_iterator>::constraints ();
        MutableRandomAccessIteratorConcept<map_array<std::size_t, T>::iterator>::constraints ();

        StorageSparseConcept<const std::map<std::size_t, T> >::constraints ();
        MutableStorageSparseConcept<std::map<std::size_t, T> >::constraints ();
        BidirectionalIteratorConcept<std::map<std::size_t, T>::const_iterator>::constraints ();
                // Not value_type mutable
        BidirectionalIteratorConcept<std::map<std::size_t, T>::iterator>::constraints ();
#endif

        // Vector
#if defined (INTERNAL_VECTOR) || defined (INTERNAL_VECTOR_DENSE)
        VectorConcept<const vector<T> >::constraints ();
        MutableVectorConcept<vector<T> >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector<T>::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector<T>::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector<T>::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector<T>::reverse_iterator>::constraints ();

        VectorConcept<zero_vector<T> >::constraints ();
        IndexedBidirectional1DIteratorConcept<zero_vector<T>::const_iterator>::constraints ();
        IndexedBidirectional1DIteratorConcept<zero_vector<T>::const_reverse_iterator>::constraints ();

        VectorConcept<unit_vector<T> >::constraints ();
        IndexedBidirectional1DIteratorConcept<unit_vector<T>::const_iterator>::constraints ();
        IndexedBidirectional1DIteratorConcept<unit_vector<T>::const_reverse_iterator>::constraints ();

        VectorConcept<scalar_vector<T> >::constraints ();
        IndexedRandomAccess1DIteratorConcept<scalar_vector<T>::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<scalar_vector<T>::const_reverse_iterator>::constraints ();

        VectorConcept<const c_vector<T, 1> >::constraints ();
        MutableVectorConcept<c_vector<T, 1> >::constraints ();
        IndexedRandomAccess1DIteratorConcept<c_vector<T, 1>::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<c_vector<T, 1>::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<c_vector<T, 1>::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<c_vector<T, 1>::reverse_iterator>::constraints ();
#endif

        // Vector Proxies
#if defined (INTERNAL_VECTOR) || defined (INTERNAL_VECTOR_PROXY)
        VectorExpressionConcept<const vector_range<const vector<T> > >::constraints ();
        MutableVectorExpressionConcept<vector_range<vector<T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_range<vector<T> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_range<vector<T> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_range<vector<T> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_range<vector<T> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const vector_slice<const vector<T> > >::constraints ();
        MutableVectorExpressionConcept<vector_slice<vector<T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_slice<vector<T> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_slice<vector<T> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_slice<vector<T> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_slice<vector<T> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const vector_indirect<const vector<T> > >::constraints ();
        MutableVectorExpressionConcept<vector_indirect<vector<T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_indirect<vector<T> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_indirect<vector<T> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_indirect<vector<T> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_indirect<vector<T> >::reverse_iterator>::constraints ();
#endif

        // Sparse Vector
#if defined (INTERNAL_SPARSE) || defined (INTERNAL_VECTOR_SPARSE)
        SparseVectorConcept<const mapped_vector<T> >::constraints ();
        MutableSparseVectorConcept<mapped_vector<T> >::constraints ();
        IndexedBidirectional1DIteratorConcept<mapped_vector<T>::const_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<mapped_vector<T>::iterator>::constraints ();
        IndexedBidirectional1DIteratorConcept<mapped_vector<T>::const_reverse_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<mapped_vector<T>::reverse_iterator>::constraints ();

        SparseVectorConcept<const compressed_vector<T> >::constraints ();
        MutableSparseVectorConcept<compressed_vector<T> >::constraints ();
        IndexedBidirectional1DIteratorConcept<compressed_vector<T>::const_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<compressed_vector<T>::iterator>::constraints ();
        IndexedBidirectional1DIteratorConcept<compressed_vector<T>::const_reverse_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<compressed_vector<T>::reverse_iterator>::constraints ();

        SparseVectorConcept<const coordinate_vector<T> >::constraints ();
        MutableSparseVectorConcept<coordinate_vector<T> >::constraints ();
        IndexedBidirectional1DIteratorConcept<coordinate_vector<T>::const_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<coordinate_vector<T>::iterator>::constraints ();
        IndexedBidirectional1DIteratorConcept<coordinate_vector<T>::const_reverse_iterator>::constraints ();
        MutableIndexedBidirectional1DIteratorConcept<coordinate_vector<T>::reverse_iterator>::constraints ();
#endif

        // Matrix
#if defined (INTERNAL_MATRIX) || defined (INTERNAL_MATRIX_DENSE)
        MatrixConcept<const matrix<T> >::constraints ();
        MutableMatrixConcept<matrix<T> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix<T>::const_iterator1,
                                             matrix<T>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix<T>::iterator1,
                                                    matrix<T>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix<T>::const_reverse_iterator1,
                                             matrix<T>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix<T>::reverse_iterator1,
                                                    matrix<T>::reverse_iterator2>::constraints ();

        MatrixConcept<const vector_of_vector<T> >::constraints ();
        MutableMatrixConcept<vector_of_vector<T> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<vector_of_vector<T>::const_iterator1,
                                             vector_of_vector<T>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<vector_of_vector<T>::iterator1,
                                                    vector_of_vector<T>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<vector_of_vector<T>::const_reverse_iterator1,
                                             vector_of_vector<T>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<vector_of_vector<T>::reverse_iterator1,
                                                    vector_of_vector<T>::reverse_iterator2>::constraints ();

        MatrixConcept<zero_matrix<T> >::constraints ();
        IndexedBidirectional2DIteratorConcept<zero_matrix<T>::const_iterator1,
                                              zero_matrix<T>::const_iterator2>::constraints ();
        IndexedBidirectional2DIteratorConcept<zero_matrix<T>::const_reverse_iterator1,
                                              zero_matrix<T>::const_reverse_iterator2>::constraints ();

        MatrixConcept<identity_matrix<T> >::constraints ();
        IndexedBidirectional2DIteratorConcept<identity_matrix<T>::const_iterator1,
                                             identity_matrix<T>::const_iterator2>::constraints ();
        IndexedBidirectional2DIteratorConcept<identity_matrix<T>::const_reverse_iterator1,
                                             identity_matrix<T>::const_reverse_iterator2>::constraints ();

        MatrixConcept<scalar_matrix<T> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<scalar_matrix<T>::const_iterator1,
                                             scalar_matrix<T>::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<scalar_matrix<T>::const_reverse_iterator1,
                                             scalar_matrix<T>::const_reverse_iterator2>::constraints ();

        MatrixConcept<const c_matrix<T, 1, 1> >::constraints ();
        MutableMatrixConcept<c_matrix<T, 1, 1> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<c_matrix<T, 1, 1>::const_iterator1,
                                             c_matrix<T, 1, 1>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<c_matrix<T, 1, 1>::iterator1,
                                                    c_matrix<T, 1, 1>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<c_matrix<T, 1, 1>::const_reverse_iterator1,
                                             c_matrix<T, 1, 1>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<c_matrix<T, 1, 1>::reverse_iterator1,
                                                    c_matrix<T, 1, 1>::reverse_iterator2>::constraints ();
#endif

        // Matrix Proxies
#if defined (INTERNAL_MATRIX) || defined (INTERNAL_MATRIX_PROXY)
        VectorExpressionConcept<const matrix_row<const matrix<T> > >::constraints ();
        MutableVectorExpressionConcept<matrix_row<matrix<T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_row<matrix<T> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_row<matrix<T> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_row<matrix<T> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_row<matrix<T> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const matrix_column<const matrix<T> > >::constraints ();
        MutableVectorExpressionConcept<matrix_column<matrix<T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_column<matrix<T> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_column<matrix<T> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_column<matrix<T> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_column<matrix<T> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const matrix_vector_range<const matrix<T> > >::constraints ();
        MutableVectorExpressionConcept<matrix_vector_range<matrix<T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_range<matrix<T> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_range<matrix<T> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_range<matrix<T> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_range<matrix<T> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const matrix_vector_slice<const matrix<T> > >::constraints ();
        MutableVectorExpressionConcept<matrix_vector_slice<matrix<T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_slice<matrix<T> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_slice<matrix<T> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_slice<matrix<T> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_slice<matrix<T> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<const matrix_vector_indirect<const matrix<T>, vector<unsigned> > >::constraints ();
        MutableVectorExpressionConcept<matrix_vector_indirect<matrix<T>, vector<unsigned> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_indirect<matrix<T>, vector<unsigned> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_indirect<matrix<T>, vector<unsigned> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_indirect<matrix<T>, vector<unsigned> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<matrix_vector_indirect<matrix<T>, vector<unsigned> >::reverse_iterator>::constraints ();

        MatrixExpressionConcept<const matrix_range<const matrix<T> > >::constraints ();
        MutableMatrixExpressionConcept<matrix_range<matrix<T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_range<matrix<T> >::const_iterator1,
                                             matrix_range<matrix<T> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_range<matrix<T> >::iterator1,
                                                    matrix_range<matrix<T> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_range<matrix<T> >::const_reverse_iterator1,
                                             matrix_range<matrix<T> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_range<matrix<T> >::reverse_iterator1,
                                                    matrix_range<matrix<T> >::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const matrix_slice<const matrix<T> > >::constraints ();
        MutableMatrixExpressionConcept<matrix_slice<matrix<T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_slice<matrix<T> >::const_iterator1,
                                             matrix_slice<matrix<T> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_slice<matrix<T> >::iterator1,
                                                    matrix_slice<matrix<T> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_slice<matrix<T> >::const_reverse_iterator1,
                                             matrix_slice<matrix<T> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_slice<matrix<T> >::reverse_iterator1,
                                                    matrix_slice<matrix<T> >::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const matrix_indirect<const matrix<T> > >::constraints ();
        MutableMatrixExpressionConcept<matrix_indirect<matrix<T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_indirect<matrix<T> >::const_iterator1,
                                             matrix_indirect<matrix<T> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_indirect<matrix<T> >::iterator1,
                                                    matrix_indirect<matrix<T> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_indirect<matrix<T> >::const_reverse_iterator1,
                                             matrix_indirect<matrix<T> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_indirect<matrix<T> >::reverse_iterator1,
                                                    matrix_indirect<matrix<T> >::reverse_iterator2>::constraints ();
#endif

        // Banded Matrix
#if defined (INTERNAL_SPECIAL) || defined (INTERNAL_BANDED)
        MatrixConcept<const banded_matrix<T> >::constraints ();
        MutableMatrixConcept<banded_matrix<T> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<banded_matrix<T>::const_iterator1,
                                             banded_matrix<T>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<banded_matrix<T>::iterator1,
                                                    banded_matrix<T>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<banded_matrix<T>::const_reverse_iterator1,
                                             banded_matrix<T>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<banded_matrix<T>::reverse_iterator1,
                                                    banded_matrix<T>::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const banded_adaptor<const matrix<T> > >::constraints ();
        MutableMatrixExpressionConcept<banded_adaptor<matrix<T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<banded_adaptor<matrix<T> >::const_iterator1,
                                             banded_adaptor<matrix<T> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<banded_adaptor<matrix<T> >::iterator1,
                                                    banded_adaptor<matrix<T> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<banded_adaptor<matrix<T> >::const_reverse_iterator1,
                                             banded_adaptor<matrix<T> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<banded_adaptor<matrix<T> >::reverse_iterator1,
                                                    banded_adaptor<matrix<T> >::reverse_iterator2>::constraints ();
#endif

        // Triangular Matrix
#if defined (INTERNAL_SPECIAL) || defined (INTERNAL_TRIANGULAR)
        MatrixConcept<const triangular_matrix<T> >::constraints ();
        MutableMatrixConcept<triangular_matrix<T> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<triangular_matrix<T>::const_iterator1,
                                             triangular_matrix<T>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<triangular_matrix<T>::iterator1,
                                                    triangular_matrix<T>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<triangular_matrix<T>::const_reverse_iterator1,
                                             triangular_matrix<T>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<triangular_matrix<T>::reverse_iterator1,
                                                    triangular_matrix<T>::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const triangular_adaptor<const matrix<T> > >::constraints ();
        MutableMatrixExpressionConcept<triangular_adaptor<matrix<T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<triangular_adaptor<matrix<T> >::const_iterator1,
                                             triangular_adaptor<matrix<T> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<triangular_adaptor<matrix<T> >::iterator1,
                                                    triangular_adaptor<matrix<T> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<triangular_adaptor<matrix<T> >::const_reverse_iterator1,
                                             triangular_adaptor<matrix<T> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<triangular_adaptor<matrix<T> >::reverse_iterator1,
                                                    triangular_adaptor<matrix<T> >::reverse_iterator2>::constraints ();
#endif

        // Symmetric Matrix
#if defined (INTERNA_SPECIAL) || defined (INTERNAL_SYMMETRIC)
        MatrixConcept<const symmetric_matrix<T> >::constraints ();
        MutableMatrixConcept<symmetric_matrix<T> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<symmetric_matrix<T>::const_iterator1,
                                             symmetric_matrix<T>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<symmetric_matrix<T>::iterator1,
                                                    symmetric_matrix<T>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<symmetric_matrix<T>::const_reverse_iterator1,
                                             symmetric_matrix<T>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<symmetric_matrix<T>::reverse_iterator1,
                                                    symmetric_matrix<T>::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const symmetric_adaptor<const matrix<T> > >::constraints ();
#ifndef SKIP_BAD
        // const_iterator (iterator) constructor is bad
        MutableMatrixExpressionConcept<symmetric_adaptor<matrix<T> > >::constraints ();
#endif
        IndexedRandomAccess2DIteratorConcept<symmetric_adaptor<matrix<T> >::const_iterator1,
                                             symmetric_adaptor<matrix<T> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<symmetric_adaptor<matrix<T> >::iterator1,
                                                    symmetric_adaptor<matrix<T> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<symmetric_adaptor<matrix<T> >::const_reverse_iterator1,
                                             symmetric_adaptor<matrix<T> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<symmetric_adaptor<matrix<T> >::reverse_iterator1,
                                                    symmetric_adaptor<matrix<T> >::reverse_iterator2>::constraints ();
#endif

        // Hermitian Matrix
#if defined (INTERNAL_SPECIAL) || defined (INTERNAL_HERMITIAN)
        MatrixConcept<const hermitian_matrix<T> >::constraints ();
        MutableMatrixConcept<hermitian_matrix<T> >::constraints ();
        IndexedRandomAccess2DIteratorConcept<hermitian_matrix<T>::const_iterator1,
                                             hermitian_matrix<T>::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<hermitian_matrix<T>::iterator1,
                                                    hermitian_matrix<T>::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<hermitian_matrix<T>::const_reverse_iterator1,
                                             hermitian_matrix<T>::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<hermitian_matrix<T>::reverse_iterator1,
                                                    hermitian_matrix<T>::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<const hermitian_adaptor<const matrix<T> > >::constraints ();
#ifndef SKIP_BAD
        // const_iterator (iterator) constructor is bad
        MutableMatrixExpressionConcept<hermitian_adaptor<matrix<T> > >::constraints ();
#endif
        IndexedRandomAccess2DIteratorConcept<hermitian_adaptor<matrix<T> >::const_iterator1,
                                             hermitian_adaptor<matrix<T> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<hermitian_adaptor<matrix<T> >::iterator1,
                                                    hermitian_adaptor<matrix<T> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<hermitian_adaptor<matrix<T> >::const_reverse_iterator1,
                                             hermitian_adaptor<matrix<T> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<hermitian_adaptor<matrix<T> >::reverse_iterator1,
                                                    hermitian_adaptor<matrix<T> >::reverse_iterator2>::constraints ();
#endif

        // Sparse Matrix
#if defined (INTERNAL_SPARSE) || defined (INTERNAL_MATRIX_SPARSE)
        {
            typedef mapped_matrix<T> container_model;
            SparseMatrixConcept<const container_model>::constraints ();
            MutableSparseMatrixConcept<container_model>::constraints ();
            IndexedBidirectional2DIteratorConcept<container_model::const_iterator1, container_model::const_iterator2>::constraints ();
            MutableIndexedBidirectional2DIteratorConcept<container_model::iterator1, container_model::iterator2>::constraints ();
            IndexedBidirectional2DIteratorConcept<container_model::const_reverse_iterator1, container_model::const_reverse_iterator2>::constraints ();
            MutableIndexedBidirectional2DIteratorConcept<container_model::reverse_iterator1, container_model::reverse_iterator2>::constraints ();
        }
        {
            typedef mapped_vector_of_mapped_vector<T> container_model;
            SparseMatrixConcept<const container_model>::constraints ();
            MutableSparseMatrixConcept<container_model>::constraints ();
            IndexedBidirectional2DIteratorConcept<container_model::const_iterator1, container_model::const_iterator2>::constraints ();
            MutableIndexedBidirectional2DIteratorConcept<container_model::iterator1, container_model::iterator2>::constraints ();
            IndexedBidirectional2DIteratorConcept<container_model::const_reverse_iterator1, container_model::const_reverse_iterator2>::constraints ();
            MutableIndexedBidirectional2DIteratorConcept<container_model::reverse_iterator1, container_model::reverse_iterator2>::constraints ();
        }
        {
            typedef compressed_matrix<T> container_model;
            SparseMatrixConcept<const container_model>::constraints ();
            MutableSparseMatrixConcept<container_model>::constraints ();
            IndexedBidirectional2DIteratorConcept<container_model::const_iterator1, container_model::const_iterator2>::constraints ();
            MutableIndexedBidirectional2DIteratorConcept<container_model::iterator1, container_model::iterator2>::constraints ();
            IndexedBidirectional2DIteratorConcept<container_model::const_reverse_iterator1, container_model::const_reverse_iterator2>::constraints ();
            MutableIndexedBidirectional2DIteratorConcept<container_model::reverse_iterator1, container_model::reverse_iterator2>::constraints ();
        }
        {
            typedef coordinate_matrix<T> container_model;
            SparseMatrixConcept<const container_model>::constraints ();
            MutableSparseMatrixConcept<container_model>::constraints ();
            IndexedBidirectional2DIteratorConcept<container_model::const_iterator1, container_model::const_iterator2>::constraints ();
            MutableIndexedBidirectional2DIteratorConcept<container_model::iterator1, container_model::iterator2>::constraints ();
            IndexedBidirectional2DIteratorConcept<container_model::const_reverse_iterator1, container_model::const_reverse_iterator2>::constraints ();
            MutableIndexedBidirectional2DIteratorConcept<container_model::reverse_iterator1, container_model::reverse_iterator2>::constraints ();
        }
        {
            typedef generalized_vector_of_vector<T, row_major, vector< coordinate_vector<T> > > container_model;
            SparseMatrixConcept<const container_model>::constraints ();
            MutableSparseMatrixConcept<container_model>::constraints ();
            IndexedBidirectional2DIteratorConcept<container_model::const_iterator1, container_model::const_iterator2>::constraints ();
            MutableIndexedBidirectional2DIteratorConcept<container_model::iterator1, container_model::iterator2>::constraints ();
            IndexedBidirectional2DIteratorConcept<container_model::const_reverse_iterator1, container_model::const_reverse_iterator2>::constraints ();
            MutableIndexedBidirectional2DIteratorConcept<container_model::reverse_iterator1, container_model::reverse_iterator2>::constraints ();
        }

#endif

        // Scalar Expressions
#if defined (INTERNAL_EXPRESSION) || defined (INTERNAL_VECTOR_EXPRESSION)
        ScalarExpressionConcept<scalar_value<T > >::constraints ();
        ScalarExpressionConcept<scalar_reference<T > >::constraints ();

        // Vector Expressions
        VectorExpressionConcept<vector_reference<vector<T> > >::constraints ();
        MutableVectorExpressionConcept<vector_reference<vector<T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_reference<vector<T> >::const_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_reference<vector<T> >::iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_reference<vector<T> >::const_reverse_iterator>::constraints ();
        MutableIndexedRandomAccess1DIteratorConcept<vector_reference<vector<T> >::reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_unary<vector<T>, scalar_identity<T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_unary<vector<T>, scalar_identity<T>  >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_unary<vector<T>, scalar_identity<T>  >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_binary<vector<T>, vector<T>, scalar_plus<T, T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary<vector<T>, vector<T>, scalar_plus<T, T> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary<vector<T>, vector<T>, scalar_plus<T, T> >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_binary_scalar1<T, vector<T>, scalar_multiplies<T, T>  > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar1<T, vector<T>, scalar_multiplies<T, T> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar1<T, vector<T>, scalar_multiplies<T, T> >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_binary_scalar2<vector<T>, scalar_value<T>, scalar_multiplies<T, T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar2<vector<T>, scalar_value<T>, scalar_multiplies<T, T> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar2<vector<T>, scalar_value<T>, scalar_multiplies<T, T> >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_binary_scalar1<scalar_value<T>, vector<T>, scalar_multiplies<T, T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar1<scalar_value<T>, vector<T>, scalar_multiplies<T, T> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar1<scalar_value<T>, vector<T>, scalar_multiplies<T, T> >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<vector_binary_scalar2<vector<T>, scalar_value<T>, scalar_multiplies<T, T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar2<vector<T>, scalar_value<T>, scalar_multiplies<T, T> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<vector_binary_scalar2<vector<T>, scalar_value<T>, scalar_multiplies<T, T> >::const_reverse_iterator>::constraints ();

        ScalarExpressionConcept<vector_scalar_unary<vector<T>, vector_sum<T> > >::constraints ();
        ScalarExpressionConcept<vector_scalar_unary<vector<T>, vector_norm_1<T> > >::constraints ();
        ScalarExpressionConcept<vector_scalar_unary<vector<T>, vector_norm_2<T> > >::constraints ();
        ScalarExpressionConcept<vector_scalar_unary<vector<T>, vector_norm_inf<T> > >::constraints ();

        ScalarExpressionConcept<vector_scalar_binary<vector<T>, vector<T>, vector_inner_prod<T, T, T> > >::constraints ();
#endif

        // Matrix Expressions
#if defined (INTERNAL_EXPRESSION) || defined (INTERNAL_MATRIX_EXPRESSION)
        MatrixExpressionConcept<matrix_reference<matrix<T> > >::constraints ();
        MutableMatrixExpressionConcept<matrix_reference<matrix<T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_reference<matrix<T> >::const_iterator1,
                                             matrix_reference<matrix<T> >::const_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_reference<matrix<T> >::iterator1,
                                                    matrix_reference<matrix<T> >::iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_reference<matrix<T> >::const_reverse_iterator1,
                                             matrix_reference<matrix<T> >::const_reverse_iterator2>::constraints ();
        MutableIndexedRandomAccess2DIteratorConcept<matrix_reference<matrix<T> >::reverse_iterator1,
                                                    matrix_reference<matrix<T> >::reverse_iterator2>::constraints ();

        MatrixExpressionConcept<vector_matrix_binary<vector<T>, vector<T>, scalar_multiplies<T, T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<vector_matrix_binary<vector<T>, vector<T>, scalar_multiplies<T, T> >::const_iterator1,
                                             vector_matrix_binary<vector<T>, vector<T>, scalar_multiplies<T, T> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<vector_matrix_binary<vector<T>, vector<T>, scalar_multiplies<T, T> >::const_reverse_iterator1,
                                             vector_matrix_binary<vector<T>, vector<T>, scalar_multiplies<T, T> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_unary1<matrix<T>, scalar_identity<T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_unary1<matrix<T>, scalar_identity<T> >::const_iterator1,
                                             matrix_unary1<matrix<T>, scalar_identity<T> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_unary1<matrix<T>, scalar_identity<T> >::const_reverse_iterator1,
                                             matrix_unary1<matrix<T>, scalar_identity<T> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_unary2<matrix<T>, scalar_identity<T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_unary2<matrix<T>, scalar_identity<T> >::const_iterator1,
                                             matrix_unary2<matrix<T>, scalar_identity<T> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_unary2<matrix<T>, scalar_identity<T> >::const_reverse_iterator1,
                                             matrix_unary2<matrix<T>, scalar_identity<T> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_binary<matrix<T>, matrix<T>, scalar_plus<T, T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary<matrix<T>, matrix<T>, scalar_plus<T, T> >::const_iterator1,
                                             matrix_binary<matrix<T>, matrix<T>, scalar_plus<T, T> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary<matrix<T>, matrix<T>, scalar_plus<T, T> >::const_reverse_iterator1,
                                             matrix_binary<matrix<T>, matrix<T>, scalar_plus<T, T> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_binary_scalar1<T, matrix<T>, scalar_multiplies<T, T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar1<T, matrix<T>, scalar_multiplies<T, T> >::const_iterator1,
                                             matrix_binary_scalar1<T, matrix<T>, scalar_multiplies<T, T> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar1<T, matrix<T>, scalar_multiplies<T, T> >::const_reverse_iterator1,
                                             matrix_binary_scalar1<T, matrix<T>, scalar_multiplies<T, T> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_binary_scalar2<matrix<T>, T, scalar_multiplies<T, T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar2<matrix<T>, T, scalar_multiplies<T, T> >::const_iterator1,
                                             matrix_binary_scalar2<matrix<T>, T, scalar_multiplies<T, T> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar2<matrix<T>, T, scalar_multiplies<T, T> >::const_reverse_iterator1,
                                             matrix_binary_scalar2<matrix<T>, T, scalar_multiplies<T, T> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_binary_scalar1<scalar_value<T>, matrix<T>, scalar_multiplies<T, T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar1<scalar_value<T>, matrix<T>, scalar_multiplies<T, T> >::const_iterator1,
                                             matrix_binary_scalar1<scalar_value<T>, matrix<T>, scalar_multiplies<T, T> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar1<scalar_value<T>, matrix<T>, scalar_multiplies<T, T> >::const_reverse_iterator1,
                                             matrix_binary_scalar1<scalar_value<T>, matrix<T>, scalar_multiplies<T, T> >::const_reverse_iterator2>::constraints ();

        MatrixExpressionConcept<matrix_binary_scalar2<matrix<T>, scalar_value<T>, scalar_multiplies<T, T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar2<matrix<T>, scalar_value<T>, scalar_multiplies<T, T> >::const_iterator1,
                                             matrix_binary_scalar2<matrix<T>, scalar_value<T>, scalar_multiplies<T, T> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_binary_scalar2<matrix<T>, scalar_value<T>, scalar_multiplies<T, T> >::const_reverse_iterator1,
                                             matrix_binary_scalar2<matrix<T>, scalar_value<T>, scalar_multiplies<T, T> >::const_reverse_iterator2>::constraints ();

        VectorExpressionConcept<matrix_vector_binary1<matrix<T>, vector<T>, matrix_vector_prod1<T, T, T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_binary1<matrix<T>, vector<T>, matrix_vector_prod1<T, T, T> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_binary1<matrix<T>, vector<T>, matrix_vector_prod1<T, T, T> >::const_reverse_iterator>::constraints ();

        VectorExpressionConcept<matrix_vector_binary2<vector<T>, matrix<T>, matrix_vector_prod2<T, T, T> > >::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_binary2<vector<T>, matrix<T>, matrix_vector_prod2<T, T, T> >::const_iterator>::constraints ();
        IndexedRandomAccess1DIteratorConcept<matrix_vector_binary2<vector<T>, matrix<T>, matrix_vector_prod2<T, T, T> >::const_reverse_iterator>::constraints ();

        MatrixExpressionConcept<matrix_matrix_binary<matrix<T>, matrix<T>, matrix_matrix_prod<T, T, T> > >::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_matrix_binary<matrix<T>, matrix<T>, matrix_matrix_prod<T, T, T> >::const_iterator1,
                                             matrix_matrix_binary<matrix<T>, matrix<T>, matrix_matrix_prod<T, T, T> >::const_iterator2>::constraints ();
        IndexedRandomAccess2DIteratorConcept<matrix_matrix_binary<matrix<T>, matrix<T>, matrix_matrix_prod<T, T, T> >::const_reverse_iterator1,
                                             matrix_matrix_binary<matrix<T>, matrix<T>, matrix_matrix_prod<T, T, T> >::const_reverse_iterator2>::constraints ();

        ScalarExpressionConcept<matrix_scalar_unary<matrix<T>, matrix_norm_1<T> > >::constraints ();
        ScalarExpressionConcept<matrix_scalar_unary<matrix<T>, matrix_norm_frobenius<T> > >::constraints ();
        ScalarExpressionConcept<matrix_scalar_unary<matrix<T>, matrix_norm_inf<T> > >::constraints ();
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
    }

}}}

#endif
