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

/** \file operation.hpp
 *  \brief This file contains some specialized products.
 */

#ifndef BOOST_UBLAS_OPERATION_H
#define BOOST_UBLAS_OPERATION_H

// axpy-based products
// Alexei Novakov had a lot of ideas to improve these. Thanks.
// Hendrik Kueck proposed some new kernel. Thanks again.

namespace boost { namespace numeric { namespace ublas {

    template<class V, class T1, class IA1, class TA1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const compressed_matrix<T1, row_major, 0, IA1, TA1> &e1,
               const vector_expression<E2> &e2,
               V &v, row_major_tag) {
        typedef typename V::size_type size_type;
        typedef typename V::value_type value_type;

        for (size_type i = 0; i < e1.size1 (); ++ i) {
            size_type begin = e1.index1_data () [i];
            size_type end = e1.index1_data () [i + 1];
            value_type t (v (i));
            for (size_type j = begin; j < end; ++ j)
                t += e1.value_data () [j] * e2 () (e1.index2_data () [j]);
            v (i) = t;
        }
        return v;
    }

    template<class V, class T1, class IA1, class TA1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const compressed_matrix<T1, column_major, 0, IA1, TA1> &e1,
               const vector_expression<E2> &e2,
               V &v, column_major_tag) {
        typedef typename V::size_type size_type;

        for (size_type j = 0; j < e1.size2 (); ++ j) {
            size_type begin = e1.index1_data () [j];
            size_type end = e1.index1_data () [j + 1];
            for (size_type i = begin; i < end; ++ i)
                v (e1.index2_data () [i]) += e1.value_data () [i] * e2 () (j);
        }
        return v;
    }

    // Dispatcher
    template<class V, class T1, class F1, class IA1, class TA1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const compressed_matrix<T1, F1, 0, IA1, TA1> &e1,
               const vector_expression<E2> &e2,
               V &v, bool init = true) {
        typedef typename V::value_type value_type;
        typedef typename F1::orientation_category orientation_category;

        if (init)
            v.assign (zero_vector<value_type> (e1.size1 ()));
#if BOOST_UBLAS_TYPE_CHECK
        vector<value_type> cv (v);
        typedef typename type_traits<value_type>::real_type real_type;
        real_type verrorbound (norm_1 (v) + norm_1 (e1) * norm_1 (e2));
        indexing_vector_assign (scalar_plus_assign<typename vector<value_type>::reference, value_type> (), cv, prod (e1, e2));
#endif
        axpy_prod (e1, e2, v, orientation_category ());
#if BOOST_UBLAS_TYPE_CHECK
        BOOST_UBLAS_CHECK (norm_1 (v - cv) <= 2 * std::numeric_limits<real_type>::epsilon () * verrorbound, internal_logic ());
#endif
        return v;
    }
    template<class V, class T1, class F, class IA1, class TA1, class E2>
    BOOST_UBLAS_INLINE
    V
    axpy_prod (const compressed_matrix<T1, F, 0, IA1, TA1> &e1,
               const vector_expression<E2> &e2) {
        typedef V vector_type;

        vector_type v (e1.size1 ());
        // FIXME: needed for c_matrix?!
        // return axpy_prod (e1, e2, v, false);
        return axpy_prod (e1, e2, v, true);
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const matrix_expression<E1> &e1,
               const vector_expression<E2> &e2,
               V &v, packed_random_access_iterator_tag, row_major_tag) {
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename V::size_type size_type;

        typename expression1_type::const_iterator1 it1 (e1 ().begin1 ());
        typename expression1_type::const_iterator1 it1_end (e1 ().end1 ());
        while (it1 != it1_end) {
            size_type index1 (it1.index1 ());
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename expression1_type::const_iterator2 it2 (it1.begin ());
            typename expression1_type::const_iterator2 it2_end (it1.end ());
#else
            typename expression1_type::const_iterator2 it2 (boost::numeric::ublas::begin (it1, iterator1_tag ()));
            typename expression1_type::const_iterator2 it2_end (boost::numeric::ublas::end (it1, iterator1_tag ()));
#endif
            while (it2 != it2_end) {
                v (index1) += *it2 * e2 () (it2.index2 ());
                ++ it2;
            }
            ++ it1;
        }
        return v;
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const matrix_expression<E1> &e1,
               const vector_expression<E2> &e2,
               V &v, packed_random_access_iterator_tag, column_major_tag) {
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename V::size_type size_type;

        typename expression1_type::const_iterator2 it2 (e1 ().begin2 ());
        typename expression1_type::const_iterator2 it2_end (e1 ().end2 ());
        while (it2 != it2_end) {
            size_type index2 (it2.index2 ());
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename expression1_type::const_iterator1 it1 (it2.begin ());
            typename expression1_type::const_iterator1 it1_end (it2.end ());
#else
            typename expression1_type::const_iterator1 it1 (boost::numeric::ublas::begin (it2, iterator2_tag ()));
            typename expression1_type::const_iterator1 it1_end (boost::numeric::ublas::end (it2, iterator2_tag ()));
#endif
            while (it1 != it1_end) {
                v (it1.index1 ()) += *it1 * e2 () (index2);
                ++ it1;
            }
            ++ it2;
        }
        return v;
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const matrix_expression<E1> &e1,
               const vector_expression<E2> &e2,
               V &v, sparse_bidirectional_iterator_tag) {
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename V::size_type size_type;

        typename expression2_type::const_iterator it (e2 ().begin ());
        typename expression2_type::const_iterator it_end (e2 ().end ());
        while (it != it_end) {
            v.plus_assign (column (e1 (), it.index ()) * *it);
            ++ it;
        }
        return v;
    }

    // Dispatcher
    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const matrix_expression<E1> &e1,
               const vector_expression<E2> &e2,
               V &v, packed_random_access_iterator_tag) {
        typedef typename E1::orientation_category orientation_category;
        return axpy_prod (e1, e2, v, packed_random_access_iterator_tag (), orientation_category ());
    }


  /** \brief computes <tt>v += A x</tt> or <tt>v = A x</tt> in an
          optimized fashion.

          \param e1 the matrix expression \c A
          \param e2 the vector expression \c x
          \param v  the result vector \c v
          \param init a boolean parameter

          <tt>axpy_prod(A, x, v, init)</tt> implements the well known
          axpy-product.  Setting \a init to \c true is equivalent to call
          <tt>v.clear()</tt> before <tt>axpy_prod</tt>. Currently \a init
          defaults to \c true, but this may change in the future.

          Up to now there are some specialisation for compressed
          matrices that give a large speed up compared to prod.
          
          \ingroup blas2

          \internal
          
          template parameters:
          \param V type of the result vector \c v
          \param E1 type of a matrix expression \c A
          \param E2 type of a vector expression \c x
  */
    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const matrix_expression<E1> &e1,
               const vector_expression<E2> &e2,
               V &v, bool init = true) {
        typedef typename V::value_type value_type;
        typedef typename E2::const_iterator::iterator_category iterator_category;

        if (init)
            v.assign (zero_vector<value_type> (e1 ().size1 ()));
#if BOOST_UBLAS_TYPE_CHECK
        vector<value_type> cv (v);
        typedef typename type_traits<value_type>::real_type real_type;
        real_type verrorbound (norm_1 (v) + norm_1 (e1) * norm_1 (e2));
        indexing_vector_assign (scalar_plus_assign<typename vector<value_type>::reference, value_type> (), cv, prod (e1, e2));
#endif
        axpy_prod (e1, e2, v, iterator_category ());
#if BOOST_UBLAS_TYPE_CHECK
        BOOST_UBLAS_CHECK (norm_1 (v - cv) <= 2 * std::numeric_limits<real_type>::epsilon () * verrorbound, internal_logic ());
#endif
        return v;
    }
    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V
    axpy_prod (const matrix_expression<E1> &e1,
               const vector_expression<E2> &e2) {
        typedef V vector_type;

        vector_type v (e1 ().size1 ());
        // FIXME: needed for c_matrix?!
        // return axpy_prod (e1, e2, v, false);
        return axpy_prod (e1, e2, v, true);
    }

    template<class V, class E1, class T2, class IA2, class TA2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const vector_expression<E1> &e1,
               const compressed_matrix<T2, column_major, 0, IA2, TA2> &e2,
               V &v, column_major_tag) {
        typedef typename V::size_type size_type;
        typedef typename V::value_type value_type;

        for (size_type j = 0; j < e2.size2 (); ++ j) {
            size_type begin = e2.index1_data () [j];
            size_type end = e2.index1_data () [j + 1];
            value_type t (v (j));
            for (size_type i = begin; i < end; ++ i)
                t += e2.value_data () [i] * e1 () (e2.index2_data () [i]);
            v (j) = t;
        }
        return v;
    }

    template<class V, class E1, class T2, class IA2, class TA2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const vector_expression<E1> &e1,
               const compressed_matrix<T2, row_major, 0, IA2, TA2> &e2,
               V &v, row_major_tag) {
        typedef typename V::size_type size_type;

        for (size_type i = 0; i < e2.size1 (); ++ i) {
            size_type begin = e2.index1_data () [i];
            size_type end = e2.index1_data () [i + 1];
            for (size_type j = begin; j < end; ++ j)
                v (e2.index2_data () [j]) += e2.value_data () [j] * e1 () (i);
        }
        return v;
    }

    // Dispatcher
    template<class V, class E1, class T2, class F2, class IA2, class TA2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const vector_expression<E1> &e1,
               const compressed_matrix<T2, F2, 0, IA2, TA2> &e2,
               V &v, bool init = true) {
        typedef typename V::value_type value_type;
        typedef typename F2::orientation_category orientation_category;

        if (init)
            v.assign (zero_vector<value_type> (e2 ().size2 ()));
#if BOOST_UBLAS_TYPE_CHECK
        vector<value_type> cv (v);
        typedef typename type_traits<value_type>::real_type real_type;
        real_type verrorbound (norm_1 (v) + norm_1 (e1) * norm_1 (e2));
        indexing_vector_assign (scalar_plus_assign<typename vector<value_type>::reference, value_type> (), cv, prod (e1, e2));
#endif
        axpy_prod (e1, e2, v, orientation_category ());
#if BOOST_UBLAS_TYPE_CHECK
        BOOST_UBLAS_CHECK (norm_1 (v - cv) <= 2 * std::numeric_limits<real_type>::epsilon () * verrorbound, internal_logic ());
#endif
        return v;
    }
    template<class V, class E1, class T2, class F2, class IA2, class TA2>
    BOOST_UBLAS_INLINE
    V
    axpy_prod (const vector_expression<E1> &e1,
               const compressed_matrix<T2, F2, 0, IA2, TA2> &e2) {
        typedef V vector_type;

        vector_type v (e2 ().size2 ());
        // FIXME: needed for c_matrix?!
        // return axpy_prod (e1, e2, v, false);
        return axpy_prod (e1, e2, v, true);
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const vector_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               V &v, packed_random_access_iterator_tag, column_major_tag) {
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename V::size_type size_type;

        typename expression2_type::const_iterator2 it2 (e2 ().begin2 ());
        typename expression2_type::const_iterator2 it2_end (e2 ().end2 ());
        while (it2 != it2_end) {
            size_type index2 (it2.index2 ());
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename expression2_type::const_iterator1 it1 (it2.begin ());
            typename expression2_type::const_iterator1 it1_end (it2.end ());
#else
            typename expression2_type::const_iterator1 it1 (boost::numeric::ublas::begin (it2, iterator2_tag ()));
            typename expression2_type::const_iterator1 it1_end (boost::numeric::ublas::end (it2, iterator2_tag ()));
#endif
            while (it1 != it1_end) {
                v (index2) += *it1 * e1 () (it1.index1 ());
                ++ it1;
            }
            ++ it2;
        }
        return v;
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const vector_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               V &v, packed_random_access_iterator_tag, row_major_tag) {
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename V::size_type size_type;

        typename expression2_type::const_iterator1 it1 (e2 ().begin1 ());
        typename expression2_type::const_iterator1 it1_end (e2 ().end1 ());
        while (it1 != it1_end) {
            size_type index1 (it1.index1 ());
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename expression2_type::const_iterator2 it2 (it1.begin ());
            typename expression2_type::const_iterator2 it2_end (it1.end ());
#else
            typename expression2_type::const_iterator2 it2 (boost::numeric::ublas::begin (it1, iterator1_tag ()));
            typename expression2_type::const_iterator2 it2_end (boost::numeric::ublas::end (it1, iterator1_tag ()));
#endif
            while (it2 != it2_end) {
                v (it2.index2 ()) += *it2 * e1 () (index1);
                ++ it2;
            }
            ++ it1;
        }
        return v;
    }

    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const vector_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               V &v, sparse_bidirectional_iterator_tag) {
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename V::size_type size_type;

        typename expression1_type::const_iterator it (e1 ().begin ());
        typename expression1_type::const_iterator it_end (e1 ().end ());
        while (it != it_end) {
            v.plus_assign (*it * row (e2 (), it.index ()));
            ++ it;
        }
        return v;
    }

    // Dispatcher
    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const vector_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               V &v, packed_random_access_iterator_tag) {
        typedef typename E2::orientation_category orientation_category;
        return axpy_prod (e1, e2, v, packed_random_access_iterator_tag (), orientation_category ());
    }


  /** \brief computes <tt>v += A<sup>T</sup> x</tt> or <tt>v = A<sup>T</sup> x</tt> in an
          optimized fashion.

          \param e1 the vector expression \c x
          \param e2 the matrix expression \c A
          \param v  the result vector \c v
          \param init a boolean parameter

          <tt>axpy_prod(x, A, v, init)</tt> implements the well known
          axpy-product.  Setting \a init to \c true is equivalent to call
          <tt>v.clear()</tt> before <tt>axpy_prod</tt>. Currently \a init
          defaults to \c true, but this may change in the future.

          Up to now there are some specialisation for compressed
          matrices that give a large speed up compared to prod.
          
          \ingroup blas2

          \internal
          
          template parameters:
          \param V type of the result vector \c v
          \param E1 type of a vector expression \c x
          \param E2 type of a matrix expression \c A
  */
    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V &
    axpy_prod (const vector_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               V &v, bool init = true) {
        typedef typename V::value_type value_type;
        typedef typename E1::const_iterator::iterator_category iterator_category;

        if (init)
            v.assign (zero_vector<value_type> (e2 ().size2 ()));
#if BOOST_UBLAS_TYPE_CHECK
        vector<value_type> cv (v);
        typedef typename type_traits<value_type>::real_type real_type;
        real_type verrorbound (norm_1 (v) + norm_1 (e1) * norm_1 (e2));
        indexing_vector_assign (scalar_plus_assign<typename vector<value_type>::reference, value_type> (), cv, prod (e1, e2));
#endif
        axpy_prod (e1, e2, v, iterator_category ());
#if BOOST_UBLAS_TYPE_CHECK
        BOOST_UBLAS_CHECK (norm_1 (v - cv) <= 2 * std::numeric_limits<real_type>::epsilon () * verrorbound, internal_logic ());
#endif
        return v;
    }
    template<class V, class E1, class E2>
    BOOST_UBLAS_INLINE
    V
    axpy_prod (const vector_expression<E1> &e1,
               const matrix_expression<E2> &e2) {
        typedef V vector_type;

        vector_type v (e2 ().size2 ());
        // FIXME: needed for c_matrix?!
        // return axpy_prod (e1, e2, v, false);
        return axpy_prod (e1, e2, v, true);
    }

    template<class M, class E1, class E2, class F>
    BOOST_UBLAS_INLINE
    M &
    axpy_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               M &m, F,
               dense_proxy_tag, row_major_tag) {
        typedef M matrix_type;
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;

#if BOOST_UBLAS_TYPE_CHECK
        matrix<value_type, row_major> cm (m);
        typedef typename type_traits<value_type>::real_type real_type;
        real_type merrorbound (norm_1 (m) + norm_1 (e1) * norm_1 (e2));
        indexing_matrix_assign (scalar_plus_assign<typename matrix<value_type, row_major>::reference, value_type> (), cm, prod (e1, e2), row_major_tag ());
#endif
        size_type size1 (e1 ().size1 ());
        size_type size2 (e1 ().size2 ());
        for (size_type i = 0; i < size1; ++ i)
            for (size_type j = 0; j < size2; ++ j)
                row (m, i).plus_assign (e1 () (i, j) * row (e2 (), j));
#if BOOST_UBLAS_TYPE_CHECK
        BOOST_UBLAS_CHECK (norm_1 (m - cm) <= 2 * std::numeric_limits<real_type>::epsilon () * merrorbound, internal_logic ());
#endif
        return m;
    }
    template<class M, class E1, class E2, class F>
    BOOST_UBLAS_INLINE
    M &
    axpy_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               M &m, F,
               sparse_proxy_tag, row_major_tag) {
        typedef M matrix_type;
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;
        typedef F functor_type;

#if BOOST_UBLAS_TYPE_CHECK
        matrix<value_type, row_major> cm (m);
        typedef typename type_traits<value_type>::real_type real_type;
        real_type merrorbound (norm_1 (m) + norm_1 (e1) * norm_1 (e2));
        indexing_matrix_assign (scalar_plus_assign<typename matrix<value_type, row_major>::reference, value_type> (), cm, prod (e1, e2), row_major_tag ());
#endif
        typename expression1_type::const_iterator1 it1 (e1 ().begin1 ());
        typename expression1_type::const_iterator1 it1_end (e1 ().end1 ());
        while (it1 != it1_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename expression1_type::const_iterator2 it2 (it1.begin ());
            typename expression1_type::const_iterator2 it2_end (it1.end ());
#else
            typename expression1_type::const_iterator2 it2 (boost::numeric::ublas::begin (it1, iterator1_tag ()));
            typename expression1_type::const_iterator2 it2_end (boost::numeric::ublas::end (it1, iterator1_tag ()));
#endif
            while (it2 != it2_end) {
                // row (m, it1.index1 ()).plus_assign (*it2 * row (e2 (), it2.index2 ()));
                matrix_row<expression2_type> mr (e2 (), it2.index2 ());
                typename matrix_row<expression2_type>::const_iterator itr (mr.begin ());
                typename matrix_row<expression2_type>::const_iterator itr_end (mr.end ());
                while (itr != itr_end) {
                    if (functor_type ().other (it1.index1 (), itr.index ()))
                        m (it1.index1 (), itr.index ()) += *it2 * *itr;
                    ++ itr;
                }
                ++ it2;
            }
            ++ it1;
        }
#if BOOST_UBLAS_TYPE_CHECK
        BOOST_UBLAS_CHECK (norm_1 (m - cm) <= 2 * std::numeric_limits<real_type>::epsilon () * merrorbound, internal_logic ());
#endif
        return m;
    }

    template<class M, class E1, class E2, class F>
    BOOST_UBLAS_INLINE
    M &
    axpy_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               M &m, F,
               dense_proxy_tag, column_major_tag) {
        typedef M matrix_type;
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;

#if BOOST_UBLAS_TYPE_CHECK
        matrix<value_type, column_major> cm (m);
        typedef typename type_traits<value_type>::real_type real_type;
        real_type merrorbound (norm_1 (m) + norm_1 (e1) * norm_1 (e2));
        indexing_matrix_assign (scalar_plus_assign<typename matrix<value_type, column_major>::reference, value_type> (), cm, prod (e1, e2), column_major_tag ());
#endif
        size_type size1 (e2 ().size1 ());
        size_type size2 (e2 ().size2 ());
        for (size_type j = 0; j < size2; ++ j)
            for (size_type i = 0; i < size1; ++ i)
                column (m, j).plus_assign (e2 () (i, j) * column (e1 (), i));
#if BOOST_UBLAS_TYPE_CHECK
        BOOST_UBLAS_CHECK (norm_1 (m - cm) <= 2 * std::numeric_limits<real_type>::epsilon () * merrorbound, internal_logic ());
#endif
        return m;
    }
    template<class M, class E1, class E2, class F>
    BOOST_UBLAS_INLINE
    M &
    axpy_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               M &m, F,
               sparse_proxy_tag, column_major_tag) {
        typedef M matrix_type;
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;
        typedef F functor_type;

#if BOOST_UBLAS_TYPE_CHECK
        matrix<value_type, column_major> cm (m);
        typedef typename type_traits<value_type>::real_type real_type;
        real_type merrorbound (norm_1 (m) + norm_1 (e1) * norm_1 (e2));
        indexing_matrix_assign (scalar_plus_assign<typename matrix<value_type, column_major>::reference, value_type> (), cm, prod (e1, e2), column_major_tag ());
#endif
        typename expression2_type::const_iterator2 it2 (e2 ().begin2 ());
        typename expression2_type::const_iterator2 it2_end (e2 ().end2 ());
        while (it2 != it2_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename expression2_type::const_iterator1 it1 (it2.begin ());
            typename expression2_type::const_iterator1 it1_end (it2.end ());
#else
            typename expression2_type::const_iterator1 it1 (boost::numeric::ublas::begin (it2, iterator2_tag ()));
            typename expression2_type::const_iterator1 it1_end (boost::numeric::ublas::end (it2, iterator2_tag ()));
#endif
            while (it1 != it1_end) {
                // column (m, it2.index2 ()).plus_assign (*it1 * column (e1 (), it1.index1 ()));
                matrix_column<expression1_type> mc (e1 (), it1.index1 ());
                typename matrix_column<expression1_type>::const_iterator itc (mc.begin ());
                typename matrix_column<expression1_type>::const_iterator itc_end (mc.end ());
                while (itc != itc_end) {
                    if (functor_type ().other (itc.index (), it2.index2 ()))
                        m (itc.index (), it2.index2 ()) += *it1 * *itc;
                    ++ itc;
                }
                ++ it1;
            }
            ++ it2;
        }
#if BOOST_UBLAS_TYPE_CHECK
        BOOST_UBLAS_CHECK (norm_1 (m - cm) <= 2 * std::numeric_limits<real_type>::epsilon () * merrorbound, internal_logic ());
#endif
        return m;
    }

    // Dispatcher
    template<class M, class E1, class E2, class F>
    BOOST_UBLAS_INLINE
    M &
    axpy_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               M &m, F, bool init = true) {
        typedef typename M::value_type value_type;
        typedef typename M::storage_category storage_category;
        typedef typename M::orientation_category orientation_category;
        typedef F functor_type;

        if (init)
            m.assign (zero_matrix<value_type> (e1 ().size1 (), e2 ().size2 ()));
        return axpy_prod (e1, e2, m, functor_type (), storage_category (), orientation_category ());
    }
    template<class M, class E1, class E2, class F>
    BOOST_UBLAS_INLINE
    M
    axpy_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               F) {
        typedef M matrix_type;
        typedef F functor_type;

        matrix_type m (e1 ().size1 (), e2 ().size2 ());
        // FIXME: needed for c_matrix?!
        // return axpy_prod (e1, e2, m, functor_type (), false);
        return axpy_prod (e1, e2, m, functor_type (), true);
    }

  /** \brief computes <tt>M += A X</tt> or <tt>M = A X</tt> in an
          optimized fashion.

          \param e1 the matrix expression \c A
          \param e2 the matrix expression \c X
          \param m  the result matrix \c M
          \param init a boolean parameter

          <tt>axpy_prod(A, X, M, init)</tt> implements the well known
          axpy-product.  Setting \a init to \c true is equivalent to call
          <tt>M.clear()</tt> before <tt>axpy_prod</tt>. Currently \a init
          defaults to \c true, but this may change in the future.

          Up to now there are no specialisations.
          
          \ingroup blas3

          \internal
          
          template parameters:
          \param M type of the result matrix \c M
          \param E1 type of a matrix expression \c A
          \param E2 type of a matrix expression \c X
  */
    template<class M, class E1, class E2>
    BOOST_UBLAS_INLINE
    M &
    axpy_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2,
               M &m, bool init = true) {
        typedef typename M::value_type value_type;
        typedef typename M::storage_category storage_category;
        typedef typename M::orientation_category orientation_category;

        if (init)
            m.assign (zero_matrix<value_type> (e1 ().size1 (), e2 ().size2 ()));
        return axpy_prod (e1, e2, m, full (), storage_category (), orientation_category ());
    }
    template<class M, class E1, class E2>
    BOOST_UBLAS_INLINE
    M
    axpy_prod (const matrix_expression<E1> &e1,
               const matrix_expression<E2> &e2) {
        typedef M matrix_type;

        matrix_type m (e1 ().size1 (), e2 ().size2 ());
        // FIXME: needed for c_matrix?!
        // return axpy_prod (e1, e2, m, full (), false);
        return axpy_prod (e1, e2, m, full (), true);
    }

    template<class M, class E1, class E2, class F>
    BOOST_UBLAS_INLINE
    M &
    opb_prod (const matrix_expression<E1> &e1,
              const matrix_expression<E2> &e2,
              M &m, F,
              dense_proxy_tag, row_major_tag) {
        typedef M matrix_type;
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;

#if BOOST_UBLAS_TYPE_CHECK
        matrix<value_type, row_major> cm (m);
        typedef typename type_traits<value_type>::real_type real_type;
        real_type merrorbound (norm_1 (m) + norm_1 (e1) * norm_1 (e2));
        indexing_matrix_assign (scalar_plus_assign<typename matrix<value_type, row_major>::reference, value_type> (), cm, prod (e1, e2), row_major_tag ());
#endif
        size_type size (BOOST_UBLAS_SAME (e1 ().size2 (), e2 ().size1 ()));
        for (size_type k = 0; k < size; ++ k) {
            vector<value_type> ce1 (column (e1 (), k));
            vector<value_type> re2 (row (e2 (), k));
            m.plus_assign (outer_prod (ce1, re2));
        }
#if BOOST_UBLAS_TYPE_CHECK
        BOOST_UBLAS_CHECK (norm_1 (m - cm) <= 2 * std::numeric_limits<real_type>::epsilon () * merrorbound, internal_logic ());
#endif
        return m;
    }

    template<class M, class E1, class E2, class F>
    BOOST_UBLAS_INLINE
    M &
    opb_prod (const matrix_expression<E1> &e1,
              const matrix_expression<E2> &e2,
              M &m, F,
              dense_proxy_tag, column_major_tag) {
        typedef M matrix_type;
        typedef const E1 expression1_type;
        typedef const E2 expression2_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;

#if BOOST_UBLAS_TYPE_CHECK
        matrix<value_type, column_major> cm (m);
        typedef typename type_traits<value_type>::real_type real_type;
        real_type merrorbound (norm_1 (m) + norm_1 (e1) * norm_1 (e2));
        indexing_matrix_assign (scalar_plus_assign<typename matrix<value_type, column_major>::reference, value_type> (), cm, prod (e1, e2), column_major_tag ());
#endif
        size_type size (BOOST_UBLAS_SAME (e1 ().size2 (), e2 ().size1 ()));
        for (size_type k = 0; k < size; ++ k) {
            vector<value_type> ce1 (column (e1 (), k));
            vector<value_type> re2 (row (e2 (), k));
            m.plus_assign (outer_prod (ce1, re2));
        }
#if BOOST_UBLAS_TYPE_CHECK
        BOOST_UBLAS_CHECK (norm_1 (m - cm) <= 2 * std::numeric_limits<real_type>::epsilon () * merrorbound, internal_logic ());
#endif
        return m;
    }

    // Dispatcher
    template<class M, class E1, class E2, class F>
    BOOST_UBLAS_INLINE
    M &
    opb_prod (const matrix_expression<E1> &e1,
              const matrix_expression<E2> &e2,
              M &m, F, bool init = true) {
        typedef typename M::value_type value_type;
        typedef typename M::storage_category storage_category;
        typedef typename M::orientation_category orientation_category;
        typedef F functor_type;

        if (init)
            m.assign (zero_matrix<value_type> (e1 ().size1 (), e2 ().size2 ()));
        return opb_prod (e1, e2, m, functor_type (), storage_category (), orientation_category ());
    }
    template<class M, class E1, class E2, class F>
    BOOST_UBLAS_INLINE
    M
    opb_prod (const matrix_expression<E1> &e1,
              const matrix_expression<E2> &e2,
              F) {
        typedef M matrix_type;
        typedef F functor_type;

        matrix_type m (e1 ().size1 (), e2 ().size2 ());
        // FIXME: needed for c_matrix?!
        // return opb_prod (e1, e2, m, functor_type (), false);
        return opb_prod (e1, e2, m, functor_type (), true);
    }

  /** \brief computes <tt>M += A X</tt> or <tt>M = A X</tt> in an
          optimized fashion.

          \param e1 the matrix expression \c A
          \param e2 the matrix expression \c X
          \param m  the result matrix \c M
          \param init a boolean parameter

          <tt>opb_prod(A, X, M, init)</tt> implements the well known
          axpy-product. Setting \a init to \c true is equivalent to call
          <tt>M.clear()</tt> before <tt>opb_prod</tt>. Currently \a init
          defaults to \c true, but this may change in the future.

          This function may give a speedup if \c A has less columns than
          rows, because the product is computed as a sum of outer
          products.
          
          \ingroup blas3

          \internal
          
          template parameters:
          \param M type of the result matrix \c M
          \param E1 type of a matrix expression \c A
          \param E2 type of a matrix expression \c X
  */
    template<class M, class E1, class E2>
    BOOST_UBLAS_INLINE
    M &
    opb_prod (const matrix_expression<E1> &e1,
              const matrix_expression<E2> &e2,
              M &m, bool init = true) {
        typedef typename M::value_type value_type;
        typedef typename M::storage_category storage_category;
        typedef typename M::orientation_category orientation_category;

        if (init)
            m.assign (zero_matrix<value_type> (e1 ().size1 (), e2 ().size2 ()));
        return opb_prod (e1, e2, m, full (), storage_category (), orientation_category ());
    }
    template<class M, class E1, class E2>
    BOOST_UBLAS_INLINE
    M
    opb_prod (const matrix_expression<E1> &e1,
              const matrix_expression<E2> &e2) {
        typedef M matrix_type;

        matrix_type m (e1 ().size1 (), e2 ().size2 ());
        // FIXME: needed for c_matrix?!
        // return opb_prod (e1, e2, m, full (), false);
        return opb_prod (e1, e2, m, full (), true);
    }

}}}

#endif
