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

#ifndef BOOST_UBLAS_MATRIX_ASSIGN_H
#define BOOST_UBLAS_MATRIX_ASSIGN_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

// Iterators based on ideas of Jeremy Siek

namespace boost { namespace numeric { namespace ublas {

    template<class E1, class E2>
    BOOST_UBLAS_INLINE
    bool equals (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
        typedef BOOST_UBLAS_TYPENAME type_traits<BOOST_UBLAS_TYPENAME promote_traits<BOOST_UBLAS_TYPENAME E1::value_type,
                                                                                     BOOST_UBLAS_TYPENAME E2::value_type>::promote_type>::real_type real_type;
#ifndef __GNUC__
        return norm_inf (e1 - e2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
               std::max<real_type> (std::max<real_type> (norm_inf (e1),
                                                         norm_inf (e2)),
                                    BOOST_UBLAS_TYPE_CHECK_MIN);
#else
        // GCC 3.1, oops?!
        return norm_inf (e1 - e2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
               (std::max) (real_type ((std::max) (real_type (norm_inf (e1)), real_type (norm_inf (e2)))),
                         real_type (BOOST_UBLAS_TYPE_CHECK_MIN));
#endif
    }


    template<class M, class E, class F>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void make_conformant (M &m, const matrix_expression<E> &e, row_major_tag, F) {
        BOOST_UBLAS_CHECK (m.size1 () == e ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (m.size2 () == e ().size2 (), bad_size ());
        typedef F functor_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
        // FIXME unbounded_array with push_back maybe better
        std::vector<std::pair<size_type, size_type> > index;
        typename M::iterator1 it1 (m.begin1 ());
        typename M::iterator1 it1_end (m.end1 ());
        typename E::const_iterator1 it1e (e ().begin1 ());
        typename E::const_iterator1 it1e_end (e ().end1 ());
        while (it1 != it1_end && it1e != it1e_end) {
            difference_type compare = it1.index1 () - it1e.index1 ();
            if (compare == 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename M::iterator2 it2 (it1.begin ());
                typename M::iterator2 it2_end (it1.end ());
                typename E::const_iterator2 it2e (it1e.begin ());
                typename E::const_iterator2 it2e_end (it1e.end ());
#else
                typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
                typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
                typename E::const_iterator2 it2e (begin (it1e, iterator1_tag ()));
                typename E::const_iterator2 it2e_end (end (it1e, iterator1_tag ()));
#endif
                if (it2 != it2_end && it2e != it2e_end) {
                    size_type it2_index = it2.index2 (), it2e_index = it2e.index2 ();
                    while (true) {
                        difference_type compare = it2_index - it2e_index;
                        if (compare == 0) {
                            ++ it2, ++ it2e;
                            if (it2 != it2_end && it2e != it2e_end) {
                                it2_index = it2.index2 ();
                                it2e_index = it2e.index2 ();
                            } else
                                break;
                        } else if (compare < 0) {
                            increment (it2, it2_end, - compare);
                            if (it2 != it2_end)
                                it2_index = it2.index2 ();
                            else
                                break;
                        } else if (compare > 0) {
                            if (functor_type::other (it2e.index1 (), it2e.index2 ()))
                                if (*it2e != value_type (0))
                                    index.push_back (std::pair<size_type, size_type> (it2e.index1 (), it2e.index2 ()));
                            ++ it2e;
                            if (it2e != it2e_end)
                                it2e_index = it2e.index2 ();
                            else
                                break;
                        }
                    }
                }
                while (it2e != it2e_end) {
                    if (functor_type::other (it2e.index1 (), it2e.index2 ()))
                        if (*it2e != value_type (0))
                            index.push_back (std::pair<size_type, size_type> (it2e.index1 (), it2e.index2 ()));
                    ++ it2e;
                }
                ++ it1, ++ it1e;
            } else if (compare < 0) {
                increment (it1, it1_end, - compare);
            } else if (compare > 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename E::const_iterator2 it2e (it1e.begin ());
                typename E::const_iterator2 it2e_end (it1e.end ());
#else
                typename E::const_iterator2 it2e (begin (it1e, iterator1_tag ()));
                typename E::const_iterator2 it2e_end (end (it1e, iterator1_tag ()));
#endif
                while (it2e != it2e_end) {
                    if (functor_type::other (it2e.index1 (), it2e.index2 ()))
                        if (*it2e != value_type (0))
                            index.push_back (std::pair<size_type, size_type> (it2e.index1 (), it2e.index2 ()));
                    ++ it2e;
                }
                ++ it1e;
            }
        }
        while (it1e != it1e_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename E::const_iterator2 it2e (it1e.begin ());
            typename E::const_iterator2 it2e_end (it1e.end ());
#else
            typename E::const_iterator2 it2e (begin (it1e, iterator1_tag ()));
            typename E::const_iterator2 it2e_end (end (it1e, iterator1_tag ()));
#endif
            while (it2e != it2e_end) {
                if (functor_type::other (it2e.index1 (), it2e.index2 ()))
                    if (*it2e != value_type (0))
                        index.push_back (std::pair<size_type, size_type> (it2e.index1 (), it2e.index2 ()));
                ++ it2e;
            }
            ++ it1e;
        }
        for (size_type k = 0; k < index.size (); ++ k)
            m (index [k].first, index [k].second) = value_type (0);
    }
    template<class M, class E, class F>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void make_conformant (M &m, const matrix_expression<E> &e, column_major_tag, F) {
        BOOST_UBLAS_CHECK (m.size1 () == e ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (m.size2 () == e ().size2 (), bad_size ());
        typedef F functor_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
        std::vector<std::pair<size_type, size_type> > index;
        typename M::iterator2 it2 (m.begin2 ());
        typename M::iterator2 it2_end (m.end2 ());
        typename E::const_iterator2 it2e (e ().begin2 ());
        typename E::const_iterator2 it2e_end (e ().end2 ());
        while (it2 != it2_end && it2e != it2e_end) {
            difference_type compare = it2.index2 () - it2e.index2 ();
            if (compare == 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename M::iterator1 it1 (it2.begin ());
                typename M::iterator1 it1_end (it2.end ());
                typename E::const_iterator1 it1e (it2e.begin ());
                typename E::const_iterator1 it1e_end (it2e.end ());
#else
                typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
                typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
                typename E::const_iterator1 it1e (begin (it2e, iterator2_tag ()));
                typename E::const_iterator1 it1e_end (end (it2e, iterator2_tag ()));
#endif
                if (it1 != it1_end && it1e != it1e_end) {
                    size_type it1_index = it1.index1 (), it1e_index = it1e.index1 ();
                    while (true) {
                        difference_type compare = it1_index - it1e_index;
                        if (compare == 0) {
                            ++ it1, ++ it1e;
                            if (it1 != it1_end && it1e != it1e_end) {
                                it1_index = it1.index1 ();
                                it1e_index = it1e.index1 ();
                            } else
                                break;
                        } else if (compare < 0) {
                            increment (it1, it1_end, - compare);
                            if (it1 != it1_end)
                                it1_index = it1.index1 ();
                            else
                                break;
                        } else if (compare > 0) {
                            if (functor_type::other (it1e.index1 (), it1e.index2 ()))
                                if (*it1e != value_type (0))
                                    index.push_back (std::pair<size_type, size_type> (it1e.index1 (), it1e.index2 ()));
                            ++ it1e;
                            if (it1e != it1e_end)
                                it1e_index = it1e.index1 ();
                            else
                                break;
                        }
                    }
                }
                while (it1e != it1e_end) {
                    if (functor_type::other (it1e.index1 (), it1e.index2 ()))
                        if (*it1e != value_type (0))
                            index.push_back (std::pair<size_type, size_type> (it1e.index1 (), it1e.index2 ()));
                    ++ it1e;
                }
                ++ it2, ++ it2e;
            } else if (compare < 0) {
                increment (it2, it2_end, - compare);
            } else if (compare > 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename E::const_iterator1 it1e (it2e.begin ());
                typename E::const_iterator1 it1e_end (it2e.end ());
#else
                typename E::const_iterator1 it1e (begin (it2e, iterator2_tag ()));
                typename E::const_iterator1 it1e_end (end (it2e, iterator2_tag ()));
#endif
                while (it1e != it1e_end) {
                    if (functor_type::other (it1e.index1 (), it1e.index2 ()))
                        if (*it1e != value_type (0))
                            index.push_back (std::pair<size_type, size_type> (it1e.index1 (), it1e.index2 ()));
                    ++ it1e;
                }
                ++ it2e;
            }
        }
        while (it2e != it2e_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename E::const_iterator1 it1e (it2e.begin ());
            typename E::const_iterator1 it1e_end (it2e.end ());
#else
            typename E::const_iterator1 it1e (begin (it2e, iterator2_tag ()));
            typename E::const_iterator1 it1e_end (end (it2e, iterator2_tag ()));
#endif
            while (it1e != it1e_end) {
                if (functor_type::other (it1e.index1 (), it1e.index2 ()))
                    if (*it1e != value_type (0))
                        index.push_back (std::pair<size_type, size_type> (it1e.index1 (), it1e.index2 ()));
                ++ it1e;
            }
            ++ it2e;
        }
        for (size_type k = 0; k < index.size (); ++ k)
            m (index [k].first, index [k].second) = value_type (0);
    }

    // Iterating row major case
    template<class F, class M, class T>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void iterating_matrix_assign_scalar (F, M &m, const T &t, row_major_tag) {
        typedef F functor_type;
        typedef typename M::difference_type difference_type;
        difference_type size1 (m.size1 ());
        difference_type size2 (m.size2 ());
        typename M::iterator1 it1 (m.begin1 ());
        BOOST_UBLAS_CHECK (size2 == 0 || m.end1 () - it1 == size1, bad_size ());
        while (-- size1 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator2 it2 (it1.begin ());
#else
            typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
#endif
            BOOST_UBLAS_CHECK (it1.end () - it2 == size2, bad_size ());
            difference_type temp_size2 (size2);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            while (-- temp_size2 >= 0)
                functor_type::apply (*it2, t), ++ it2;
#else
            DD (temp_size2, 4, r, (functor_type::apply (*it2, t), ++ it2));
#endif
            ++ it1;
        }
    }
    // Iterating column major case
    template<class F, class M, class T>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void iterating_matrix_assign_scalar (F, M &m, const T &t, column_major_tag) {
        typedef F functor_type;
        typedef typename M::difference_type difference_type;
        difference_type size2 (m.size2 ());
        difference_type size1 (m.size1 ());
        typename M::iterator2 it2 (m.begin2 ());
        BOOST_UBLAS_CHECK (size1 == 0 || m.end2 () - it2 == size2, bad_size ());
        while (-- size2 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator1 it1 (it2.begin ());
#else
            typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
#endif
            BOOST_UBLAS_CHECK (it2.end () - it1 == size1, bad_size ());
            difference_type temp_size1 (size1);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            while (-- temp_size1 >= 0)
                functor_type::apply (*it1, t), ++ it1;
#else
            DD (temp_size1, 4, r, (functor_type::apply (*it1, t), ++ it1));
#endif
            ++ it2;
        }
    }
    // Indexing row major case
    template<class F, class M, class T>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void indexing_matrix_assign_scalar (F, M &m, const T &t, row_major_tag) {
        typedef F functor_type;
        typedef typename M::size_type size_type;
        size_type size1 (m.size1 ());
        size_type size2 (m.size2 ());
        for (size_type i = 0; i < size1; ++ i) {
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            for (size_type j = 0; j < size2; ++ j)
                functor_type::apply (m (i, j), t);
#else
            size_type j (0);
            DD (size2, 4, r, (functor_type::apply (m (i, j), t), ++ j));
#endif
        }
    }
    // Indexing column major case
    template<class F, class M, class T>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void indexing_matrix_assign_scalar (F, M &m, const T &t, column_major_tag) {
        typedef F functor_type;
        typedef typename M::size_type size_type;
        size_type size2 (m.size2 ());
        size_type size1 (m.size1 ());
        for (size_type j = 0; j < size2; ++ j) {
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            for (size_type i = 0; i < size1; ++ i)
                functor_type::apply (m (i, j), t);
#else
            size_type i (0);
            DD (size1, 4, r, (functor_type::apply (m (i, j), t), ++ i));
#endif
        }
    }

    // Dense (proxy) case
    template<class F, class M, class T, class C>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign_scalar (F, M &m, const T &t, dense_proxy_tag, C) {
        typedef F functor_type;
        typedef C orientation_category;
#ifdef BOOST_UBLAS_USE_INDEXING
        indexing_matrix_assign_scalar (functor_type (), m, t, orientation_category ());
#elif BOOST_UBLAS_USE_ITERATING
        iterating_matrix_assign_scalar (functor_type (), m, t, orientation_category ());
#else
        typedef typename M::size_type size_type;
        size_type size1 (m.size1 ());
        size_type size2 (m.size2 ());
        if (size1 >= BOOST_UBLAS_ITERATOR_THRESHOLD &&
            size2 >= BOOST_UBLAS_ITERATOR_THRESHOLD)
            iterating_matrix_assign_scalar (functor_type (), m, t, orientation_category ());
        else
            indexing_matrix_assign_scalar (functor_type (), m, t, orientation_category ());
#endif
    }
    // Packed (proxy) row major case
    template<class F, class M, class T>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign_scalar (F, M &m, const T &t, packed_proxy_tag, row_major_tag) {
        typedef F functor_type;
        typedef typename M::difference_type difference_type;
        typename M::iterator1 it1 (m.begin1 ());
        difference_type size1 (m.end1 () - it1);
        while (-- size1 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator2 it2 (it1.begin ());
            difference_type size2 (it1.end () - it2);
#else
            typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
            difference_type size2 (end (it1, iterator1_tag ()) - it2);
#endif
            while (-- size2 >= 0)
                functor_type::apply (*it2, t), ++ it2;
            ++ it1;
        }
    }
    // Packed (proxy) column major case
    template<class F, class M, class T>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign_scalar (F, M &m, const T &t, packed_proxy_tag, column_major_tag) {
        typedef F functor_type;
        typedef typename M::difference_type difference_type;
        typename M::iterator2 it2 (m.begin2 ());
        difference_type size2 (m.end2 () - it2);
        while (-- size2 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator1 it1 (it2.begin ());
            difference_type size1 (it2.end () - it1);
#else
            typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
            difference_type size1 (end (it2, iterator2_tag ()) - it1);
#endif
            while (-- size1 >= 0)
                functor_type::apply (*it1, t), ++ it1;
            ++ it2;
        }
    }
    // Sparse (proxy) row major case
    template<class F, class M, class T>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign_scalar (F, M &m, const T &t, sparse_proxy_tag, row_major_tag) {
        typedef F functor_type;
        typename M::iterator1 it1 (m.begin1 ());
        typename M::iterator1 it1_end (m.end1 ());
        while (it1 != it1_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator2 it2 (it1.begin ());
            typename M::iterator2 it2_end (it1.end ());
#else
            typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
            typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
#endif
            while (it2 != it2_end)
                functor_type::apply (*it2, t), ++ it2;
            ++ it1;
        }
    }
    // Sparse (proxy) column major case
    template<class F, class M, class T>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign_scalar (F, M &m, const T &t, sparse_proxy_tag, column_major_tag) {
        typedef F functor_type;
        typename M::iterator2 it2 (m.begin2 ());
        typename M::iterator2 it2_end (m.end2 ());
        while (it2 != it2_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator1 it1 (it2.begin ());
            typename M::iterator1 it1_end (it2.end ());
#else
            typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
            typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
#endif
            while (it1 != it1_end)
                functor_type::apply (*it1, t), ++ it1;
            ++ it2;
        }
    }

    // Dispatcher
    template<class F, class M, class T>
    BOOST_UBLAS_INLINE
    void matrix_assign_scalar (F, M &m, const T &t) {
        typedef F functor_type;
        typedef typename M::storage_category storage_category;
        typedef typename M::orientation_category orientation_category;
        matrix_assign_scalar (functor_type (), m, t, storage_category (), orientation_category ());
    }

    template<class LS, class A, class RI1, class RI2>
    struct matrix_assign_traits {
        typedef LS storage_category;
    };

    template<>
    struct matrix_assign_traits<dense_tag, assign_tag, packed_random_access_iterator_tag, packed_random_access_iterator_tag> {
        typedef packed_tag storage_category;
    };
    template<>
    struct matrix_assign_traits<dense_tag, computed_assign_tag, packed_random_access_iterator_tag, packed_random_access_iterator_tag> {
        typedef packed_tag storage_category;
    };
    template<>
    struct matrix_assign_traits<dense_tag, assign_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_tag storage_category;
    };
    template<>
    struct matrix_assign_traits<dense_tag, computed_assign_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    template<>
    struct matrix_assign_traits<dense_proxy_tag, assign_tag, packed_random_access_iterator_tag, packed_random_access_iterator_tag> {
        typedef packed_proxy_tag storage_category;
    };
    template<>
    struct matrix_assign_traits<dense_proxy_tag, computed_assign_tag, packed_random_access_iterator_tag, packed_random_access_iterator_tag> {
        typedef packed_proxy_tag storage_category;
    };
    template<>
    struct matrix_assign_traits<dense_proxy_tag, assign_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };
    template<>
    struct matrix_assign_traits<dense_proxy_tag, computed_assign_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    template<>
    struct matrix_assign_traits<packed_tag, assign_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_tag storage_category;
    };
    template<>
    struct matrix_assign_traits<packed_tag, computed_assign_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    template<>
    struct matrix_assign_traits<packed_proxy_tag, assign_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };
    template<>
    struct matrix_assign_traits<packed_proxy_tag, computed_assign_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    template<>
    struct matrix_assign_traits<sparse_tag, computed_assign_tag, dense_random_access_iterator_tag, dense_random_access_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };
    template<>
    struct matrix_assign_traits<sparse_tag, computed_assign_tag, packed_random_access_iterator_tag, packed_random_access_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };
    template<>
    struct matrix_assign_traits<sparse_tag, computed_assign_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    // Iterating row major case
    template<class F, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void iterating_matrix_assign (F, M &m, const matrix_expression<E> &e, row_major_tag) {
        typedef F functor_type;
        typedef typename M::difference_type difference_type;
        difference_type size1 (BOOST_UBLAS_SAME (m.size1 (), e ().size1 ()));
        difference_type size2 (BOOST_UBLAS_SAME (m.size2 (), e ().size2 ()));
        typename M::iterator1 it1 (m.begin1 ());
        BOOST_UBLAS_CHECK (size2 == 0 || m.end1 () - it1 == size1, bad_size ());
        typename E::const_iterator1 it1e (e ().begin1 ());
        BOOST_UBLAS_CHECK (size2 == 0 || e ().end1 () - it1e == size1, bad_size ());
        while (-- size1 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator2 it2 (it1.begin ());
            typename E::const_iterator2 it2e (it1e.begin ());
#else
            typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
            typename E::const_iterator2 it2e (begin (it1e, iterator1_tag ()));
#endif
            BOOST_UBLAS_CHECK (it1.end () - it2 == size2, bad_size ());
            BOOST_UBLAS_CHECK (it1e.end () - it2e == size2, bad_size ());
            difference_type temp_size2 (size2);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            while (-- temp_size2 >= 0)
                functor_type::apply (*it2, *it2e), ++ it2, ++ it2e;
#else
            DD (temp_size2, 2, r, (functor_type::apply (*it2, *it2e), ++ it2, ++ it2e));
#endif
            ++ it1, ++ it1e;
        }
    }
    // Iterating column major case
    template<class F, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void iterating_matrix_assign (F, M &m, const matrix_expression<E> &e, column_major_tag) {
        typedef F functor_type;
        typedef typename M::difference_type difference_type;
        difference_type size2 (BOOST_UBLAS_SAME (m.size2 (), e ().size2 ()));
        difference_type size1 (BOOST_UBLAS_SAME (m.size1 (), e ().size1 ()));
        typename M::iterator2 it2 (m.begin2 ());
        BOOST_UBLAS_CHECK (size1 == 0 || m.end2 () - it2 == size2, bad_size ());
        typename E::const_iterator2 it2e (e ().begin2 ());
        BOOST_UBLAS_CHECK (size1 == 0 || e ().end2 () - it2e == size2, bad_size ());
        while (-- size2 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator1 it1 (it2.begin ());
            typename E::const_iterator1 it1e (it2e.begin ());
#else
            typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
            typename E::const_iterator1 it1e (begin (it2e, iterator2_tag ()));
#endif
            BOOST_UBLAS_CHECK (it2.end () - it1 == size1, bad_size ());
            BOOST_UBLAS_CHECK (it2e.end () - it1e == size1, bad_size ());
            difference_type temp_size1 (size1);
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            while (-- temp_size1 >= 0)
                functor_type::apply (*it1, *it1e), ++ it1, ++ it1e;
#else
            DD (temp_size1, 2, r, (functor_type::apply (*it1, *it1e), ++ it1, ++ it1e));
#endif
            ++ it2, ++ it2e;
        }
    }
    // Indexing row major case
    template<class F, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void indexing_matrix_assign (F, M &m, const matrix_expression<E> &e, row_major_tag) {
        typedef F functor_type;
        typedef typename M::size_type size_type;
        size_type size1 (BOOST_UBLAS_SAME (m.size1 (), e ().size1 ()));
        size_type size2 (BOOST_UBLAS_SAME (m.size2 (), e ().size2 ()));
        for (size_type i = 0; i < size1; ++ i) {
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            for (size_type j = 0; j < size2; ++ j)
                functor_type::apply (m (i, j), e () (i, j));
#else
            size_type j (0);
            DD (size2, 2, r, (functor_type::apply (m (i, j), e () (i, j)), ++ j));
#endif
        }
    }
    // Indexing column major case
    template<class F, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void indexing_matrix_assign (F, M &m, const matrix_expression<E> &e, column_major_tag) {
        typedef F functor_type;
        typedef typename M::size_type size_type;
        size_type size2 (BOOST_UBLAS_SAME (m.size2 (), e ().size2 ()));
        size_type size1 (BOOST_UBLAS_SAME (m.size1 (), e ().size1 ()));
        for (size_type j = 0; j < size2; ++ j) {
#ifndef BOOST_UBLAS_USE_DUFF_DEVICE
            for (size_type i = 0; i < size1; ++ i)
                functor_type::apply (m (i, j), e () (i, j));
#else
            size_type i (0);
            DD (size1, 2, r, (functor_type::apply (m (i, j), e () (i, j)), ++ i));
#endif
        }
    }

    // Dense (proxy) case
    template<class F, class M, class E, class C>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign (F, full, M &m, const matrix_expression<E> &e, dense_proxy_tag, C) {
        typedef F functor_type;
        typedef C orientation_category;
#ifdef BOOST_UBLAS_USE_INDEXING
        indexing_matrix_assign (functor_type (), m, e, orientation_category ());
#elif BOOST_UBLAS_USE_ITERATING
        iterating_matrix_assign (functor_type (), m, e, orientation_category ());
#else
        typedef typename M::difference_type difference_type;
        size_type size1 (BOOST_UBLAS_SAME (m.size1 (), e ().size1 ()));
        size_type size2 (BOOST_UBLAS_SAME (m.size2 (), e ().size2 ()));
        if (size1 >= BOOST_UBLAS_ITERATOR_THRESHOLD &&
            size2 >= BOOST_UBLAS_ITERATOR_THRESHOLD)
            iterating_matrix_assign (functor_type (), m, e, orientation_category ());
        else
            indexing_matrix_assign (functor_type (), m, e, orientation_category ());
#endif
    }
    // Packed (proxy) row major case
    template<class F1, class F2, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign (F1, F2, M &m, const matrix_expression<E> &e, packed_proxy_tag, row_major_tag) {
        BOOST_UBLAS_CHECK (m.size1 () == e ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (m.size2 () == e ().size2 (), bad_size ());
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
#if BOOST_UBLAS_TYPE_CHECK
        matrix<value_type, row_major> cm (m.size1 (), m.size2 ());
#ifndef BOOST_UBLAS_NO_ELEMENT_PROXIES
        indexing_matrix_assign (scalar_assign<typename matrix<value_type, row_major>::reference, value_type> (), cm, m, row_major_tag ());
        indexing_matrix_assign (functor1_type::template make_debug_functor<typename matrix<value_type, row_major>::reference, value_type> (), cm, e, row_major_tag ());
#else
        indexing_matrix_assign (scalar_assign<value_type&, value_type> (), cm, m, row_major_tag ());
        indexing_matrix_assign (functor1_type (), cm, e, row_major_tag ());
#endif
#endif
        typename M::iterator1 it1 (m.begin1 ());
        typename M::iterator1 it1_end (m.end1 ());
        typename E::const_iterator1 it1e (e ().begin1 ());
        typename E::const_iterator1 it1e_end (e ().end1 ());
        difference_type it1_size (it1_end - it1);
        difference_type it1e_size (it1e_end - it1e);
        difference_type diff1 (0);
        if (it1_size > 0 && it1e_size > 0)
            diff1 = it1.index1 () - it1e.index1 ();
        if (diff1 != 0) {
            difference_type size1 = (std::min) (diff1, it1e_size);
            if (size1 > 0) {
                it1e += size1;
                it1e_size -= size1;
                diff1 -= size1;
            }
            size1 = (std::min) (- diff1, it1_size);
            if (size1 > 0) {
                it1_size -= size1;
                if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
                    while (-- size1 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                        typename M::iterator2 it2 (it1.begin ());
                        typename M::iterator2 it2_end (it1.end ());
#else
                        typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
                        typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
#endif
                        difference_type size2 (it2_end - it2);
                        while (-- size2 >= 0)
                            functor1_type::apply (*it2, value_type (0)), ++ it2;
                        ++ it1;
                    }
                } else {
                    it1 += size1;
                }
                diff1 += size1;
            }
        }
        difference_type size1 ((std::min) (it1_size, it1e_size));
        it1_size -= size1;
        it1e_size -= size1;
        while (-- size1 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator2 it2 (it1.begin ());
            typename M::iterator2 it2_end (it1.end ());
            typename E::const_iterator2 it2e (it1e.begin ());
            typename E::const_iterator2 it2e_end (it1e.end ());
#else
            typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
            typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
            typename E::const_iterator2 it2e (begin (it1e, iterator1_tag ()));
            typename E::const_iterator2 it2e_end (end (it1e, iterator1_tag ()));
#endif
            difference_type it2_size (it2_end - it2);
            difference_type it2e_size (it2e_end - it2e);
            difference_type diff2 (0);
            if (it2_size > 0 && it2e_size > 0) {
                diff2 = it2.index2 () - it2e.index2 ();
                difference_type size2 = (std::min) (diff2, it2e_size);
                if (size2 > 0) {
                    it2e += size2;
                    it2e_size -= size2;
                    diff2 -= size2;
                }
                size2 = (std::min) (- diff2, it2_size);
                if (size2 > 0) {
                    it2_size -= size2;
                    if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
                        while (-- size2 >= 0)
                            functor1_type::apply (*it2, value_type (0)), ++ it2;
                    } else {
                        it2 += size2;
                    }
                    diff2 += size2;
                }
            }
            difference_type size2 ((std::min) (it2_size, it2e_size));
            it2_size -= size2;
            it2e_size -= size2;
            while (-- size2 >= 0)
                functor1_type::apply (*it2, *it2e), ++ it2, ++ it2e;
            size2 = it2_size;
            if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
                while (-- size2 >= 0)
                    functor1_type::apply (*it2, value_type (0)), ++ it2;
            } else {
                it2 += size2;
            }
            ++ it1, ++ it1e;
        }
        size1 = it1_size;
        if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
            while (-- size1 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename M::iterator2 it2 (it1.begin ());
                typename M::iterator2 it2_end (it1.end ());
#else
                typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
                typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
#endif
                difference_type size2 (it2_end - it2);
                while (-- size2 >= 0)
                    functor1_type::apply (*it2, value_type (0)), ++ it2;
                ++ it1;
            }
        } else {
            it1 += size1;
        }
#if BOOST_UBLAS_TYPE_CHECK
        if (! disable_type_check<bool>::value)
            BOOST_UBLAS_CHECK (equals (m, cm), external_logic ());
#endif
    }
    // Packed (proxy) column major case
    template<class F1, class F2, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign (F1, F2, M &m, const matrix_expression<E> &e, packed_proxy_tag, column_major_tag) {
        BOOST_UBLAS_CHECK (m.size2 () == e ().size2 (), bad_size ());
        BOOST_UBLAS_CHECK (m.size1 () == e ().size1 (), bad_size ());
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
#if BOOST_UBLAS_TYPE_CHECK
        matrix<value_type, column_major> cm (m.size1 (), m.size2 ());
#ifndef BOOST_UBLAS_NO_ELEMENT_PROXIES
        indexing_matrix_assign (scalar_assign<typename matrix<value_type, column_major>::reference, value_type> (), cm, m, column_major_tag ());
        indexing_matrix_assign (functor1_type::template make_debug_functor<typename matrix<value_type, column_major>::reference, value_type> (), cm, e, column_major_tag ());
#else
        indexing_matrix_assign (scalar_assign<value_type&, value_type> (), cm, m, column_major_tag ());
        indexing_matrix_assign (functor1_type (), cm, e, column_major_tag ());
#endif
#endif
        typename M::iterator2 it2 (m.begin2 ());
        typename M::iterator2 it2_end (m.end2 ());
        typename E::const_iterator2 it2e (e ().begin2 ());
        typename E::const_iterator2 it2e_end (e ().end2 ());
        difference_type it2_size (it2_end - it2);
        difference_type it2e_size (it2e_end - it2e);
        difference_type diff2 (0);
        if (it2_size > 0 && it2e_size > 0)
            diff2 = it2.index2 () - it2e.index2 ();
        if (diff2 != 0) {
            difference_type size2 = (std::min) (diff2, it2e_size);
            if (size2 > 0) {
                it2e += size2;
                it2e_size -= size2;
                diff2 -= size2;
            }
            size2 = (std::min) (- diff2, it2_size);
            if (size2 > 0) {
                it2_size -= size2;
                if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
                    while (-- size2 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                        typename M::iterator1 it1 (it2.begin ());
                        typename M::iterator1 it1_end (it2.end ());
#else
                        typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
                        typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
#endif
                        difference_type size1 (it1_end - it1);
                        while (-- size1 >= 0)
                            functor1_type::apply (*it1, value_type (0)), ++ it1;
                        ++ it2;
                    }
                } else {
                    it2 += size2;
                }
                diff2 += size2;
            }
        }
        difference_type size2 ((std::min) (it2_size, it2e_size));
        it2_size -= size2;
        it2e_size -= size2;
        while (-- size2 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator1 it1 (it2.begin ());
            typename M::iterator1 it1_end (it2.end ());
            typename E::const_iterator1 it1e (it2e.begin ());
            typename E::const_iterator1 it1e_end (it2e.end ());
#else
            typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
            typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
            typename E::const_iterator1 it1e (begin (it2e, iterator2_tag ()));
            typename E::const_iterator1 it1e_end (end (it2e, iterator2_tag ()));
#endif
            difference_type it1_size (it1_end - it1);
            difference_type it1e_size (it1e_end - it1e);
            difference_type diff1 (0);
            if (it1_size > 0 && it1e_size > 0) {
                diff1 = it1.index1 () - it1e.index1 ();
                difference_type size1 = (std::min) (diff1, it1e_size);
                if (size1 > 0) {
                    it1e += size1;
                    it1e_size -= size1;
                    diff1 -= size1;
                }
                size1 = (std::min) (- diff1, it1_size);
                if (size1 > 0) {
                    it1_size -= size1;
                    if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
                        while (-- size1 >= 0)
                            functor1_type::apply (*it1, value_type (0)), ++ it1;
                    } else {
                        it1 += size1;
                    }
                    diff1 += size1;
                }
            }
            difference_type size1 ((std::min) (it1_size, it1e_size));
            it1_size -= size1;
            it1e_size -= size1;
            while (-- size1 >= 0)
                functor1_type::apply (*it1, *it1e), ++ it1, ++ it1e;
            size1 = it1_size;
            if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
                while (-- size1 >= 0)
                    functor1_type::apply (*it1, value_type (0)), ++ it1;
            } else {
                it1 += size1;
            }
            ++ it2, ++ it2e;
        }
        size2 = it2_size;
        if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
            while (-- size2 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename M::iterator1 it1 (it2.begin ());
                typename M::iterator1 it1_end (it2.end ());
#else
                typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
                typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
#endif
                difference_type size1 (it1_end - it1);
                while (-- size1 >= 0)
                    functor1_type::apply (*it1, value_type (0)), ++ it1;
                ++ it2;
            }
        } else {
            it2 += size2;
        }
#if BOOST_UBLAS_TYPE_CHECK
        if (! disable_type_check<bool>::value)
            BOOST_UBLAS_CHECK (equals (m, cm), external_logic ());
#endif
    }
    // Sparse row major case
    template<class VT, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME M::reference, VT>, full, M &m, const matrix_expression<E> &e, sparse_tag, row_major_tag) {
        BOOST_UBLAS_CHECK (m.size1 () == e ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (m.size2 () == e ().size2 (), bad_size ());
        typedef typename M::value_type value_type;
#if BOOST_UBLAS_TYPE_CHECK
        // DEPRECATED Sparse type has no numeric constraints
#endif
        m.clear ();
        typename E::const_iterator1 it1e (e ().begin1 ());
        typename E::const_iterator1 it1e_end (e ().end1 ());
        while (it1e != it1e_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename E::const_iterator2 it2e (it1e.begin ());
            typename E::const_iterator2 it2e_end (it1e.end ());
#else
            typename E::const_iterator2 it2e (begin (it1e, iterator1_tag ()));
            typename E::const_iterator2 it2e_end (end (it1e, iterator1_tag ()));
#endif
            while (it2e != it2e_end) {
                value_type t (*it2e);
                if (t != value_type (0))
                    m.insert (it2e.index1 (), it2e.index2 (), t);
                ++ it2e;
            }
            ++ it1e;
        }
    }
    // Sparse column major case
    template<class VT, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME M::reference, VT>, full, M &m, const matrix_expression<E> &e, sparse_tag, column_major_tag) {
        BOOST_UBLAS_CHECK (m.size1 () == e ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (m.size2 () == e ().size2 (), bad_size ());
        typedef typename M::value_type value_type;
#if BOOST_UBLAS_TYPE_CHECK
        // DEPRECATED Sparse type has no numeric constraints
#endif
        m.clear ();
        typename E::const_iterator2 it2e (e ().begin2 ());
        typename E::const_iterator2 it2e_end (e ().end2 ());
        while (it2e != it2e_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename E::const_iterator1 it1e (it2e.begin ());
            typename E::const_iterator1 it1e_end (it2e.end ());
#else
            typename E::const_iterator1 it1e (begin (it2e, iterator2_tag ()));
            typename E::const_iterator1 it1e_end (end (it2e, iterator2_tag ()));
#endif
            while (it1e != it1e_end) {
                value_type t (*it1e);
                if (t != value_type (0))
                    m.insert (it1e.index1 (), it1e.index2 (), t);
                ++ it1e;
            }
            ++ it2e;
        }
    }
    // Sparse proxy row major case
    template<class F1, class F2, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign (F1, F2, M &m, const matrix_expression<E> &e, sparse_proxy_tag, row_major_tag) {
        BOOST_UBLAS_CHECK (m.size1 () == e ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (m.size2 () == e ().size2 (), bad_size ());
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
#if BOOST_UBLAS_TYPE_CHECK
        matrix<value_type, row_major> cm (m.size1 (), m.size2 ());
#ifndef BOOST_UBLAS_NO_ELEMENT_PROXIES
        indexing_matrix_assign (scalar_assign<typename matrix<value_type, row_major>::reference, value_type> (), cm, m, row_major_tag ());
        indexing_matrix_assign (functor1_type::template make_debug_functor<typename matrix<value_type, row_major>::reference, value_type> (), cm, e, row_major_tag ());
#else
        indexing_matrix_assign (scalar_assign<value_type&, value_type> (), cm, m, row_major_tag ());
        indexing_matrix_assign (functor1_type (), cm, e, row_major_tag ());
#endif
#endif
#ifdef BOOST_UBLAS_NON_CONFORMANT_PROXIES
        make_conformant (m, e, row_major_tag (), functor2_type ());
#endif
        typename M::iterator1 it1 (m.begin1 ());
        typename M::iterator1 it1_end (m.end1 ());
        typename E::const_iterator1 it1e (e ().begin1 ());
        typename E::const_iterator1 it1e_end (e ().end1 ());
        while (it1 != it1_end && it1e != it1e_end) {
            difference_type compare = it1.index1 () - it1e.index1 ();
            if (compare == 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename M::iterator2 it2 (it1.begin ());
                typename M::iterator2 it2_end (it1.end ());
                typename E::const_iterator2 it2e (it1e.begin ());
                typename E::const_iterator2 it2e_end (it1e.end ());
#else
                typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
                typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
                typename E::const_iterator2 it2e (begin (it1e, iterator1_tag ()));
                typename E::const_iterator2 it2e_end (end (it1e, iterator1_tag ()));
#endif
                if (it2 != it2_end && it2e != it2e_end) {
                    size_type it2_index = it2.index2 (), it2e_index = it2e.index2 ();
                    while (true) {
                        difference_type compare = it2_index - it2e_index;
                        if (compare == 0) {
                            functor1_type::apply (*it2, *it2e);
                            ++ it2, ++ it2e;
                            if (it2 != it2_end && it2e != it2e_end) {
                                it2_index = it2.index2 ();
                                it2e_index = it2e.index2 ();
                            } else
                                break;
                        } else if (compare < 0) {
                            if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
                                functor1_type::apply (*it2, value_type (0));
                                ++ it2;
                            } else
                                increment (it2, it2_end, - compare);
                            if (it2 != it2_end)
                                it2_index = it2.index2 ();
                            else
                                break;
                        } else if (compare > 0) {
                            increment (it2e, it2e_end, compare);
                            if (it2e != it2e_end)
                                it2e_index = it2e.index2 ();
                            else
                                break;
                        }
                    }
                }
                if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
                    while (it2 != it2_end) {
                        functor1_type::apply (*it2, value_type (0));
                        ++ it2;
                    }
                } else {
                    it2 = it2_end;
                }
                ++ it1, ++ it1e;
            } else if (compare < 0) {
                if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                    typename M::iterator2 it2 (it1.begin ());
                    typename M::iterator2 it2_end (it1.end ());
#else
                    typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
                    typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
#endif
                    while (it2 != it2_end) {
                        functor1_type::apply (*it2, value_type (0));
                        ++ it2;
                    }
                    ++ it1;
                } else {
                    increment (it1, it1_end, - compare);
                }
            } else if (compare > 0) {
                increment (it1e, it1e_end, compare);
            }
        }
        if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
            while (it1 != it1_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename M::iterator2 it2 (it1.begin ());
                typename M::iterator2 it2_end (it1.end ());
#else
                typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
                typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
#endif
                while (it2 != it2_end) {
                    functor1_type::apply (*it2, value_type (0));
                    ++ it2;
                }
                ++ it1;
            }
        } else {
            it1 = it1_end;
        }
#if BOOST_UBLAS_TYPE_CHECK
        if (! disable_type_check<bool>::value)
            BOOST_UBLAS_CHECK (equals (m, cm), external_logic ());
#endif
    }
    // Sparse proxy column major case
    template<class F1, class F2, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_assign (F1, F2, M &m, const matrix_expression<E> &e, sparse_proxy_tag, column_major_tag) {
        BOOST_UBLAS_CHECK (m.size1 () == e ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (m.size2 () == e ().size2 (), bad_size ());
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
#if BOOST_UBLAS_TYPE_CHECK
        matrix<value_type, column_major> cm (m.size1 (), m.size2 ());
#ifndef BOOST_UBLAS_NO_ELEMENT_PROXIES
        indexing_matrix_assign (scalar_assign<typename matrix<value_type, column_major>::reference, value_type> (), cm, m, column_major_tag ());
        indexing_matrix_assign (functor1_type::template make_debug_functor<typename matrix<value_type, column_major>::reference, value_type> (), cm, e, column_major_tag ());
#else
        indexing_matrix_assign (scalar_assign<value_type&, value_type> (), cm, m, column_major_tag ());
        indexing_matrix_assign (functor1_type (), cm, e, column_major_tag ());
#endif
#endif
#ifdef BOOST_UBLAS_NON_CONFORMANT_PROXIES
        make_conformant (m, e, column_major_tag (), functor2_type ());
#endif
        typename M::iterator2 it2 (m.begin2 ());
        typename M::iterator2 it2_end (m.end2 ());
        typename E::const_iterator2 it2e (e ().begin2 ());
        typename E::const_iterator2 it2e_end (e ().end2 ());
        while (it2 != it2_end && it2e != it2e_end) {
            difference_type compare = it2.index2 () - it2e.index2 ();
            if (compare == 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename M::iterator1 it1 (it2.begin ());
                typename M::iterator1 it1_end (it2.end ());
                typename E::const_iterator1 it1e (it2e.begin ());
                typename E::const_iterator1 it1e_end (it2e.end ());
#else
                typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
                typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
                typename E::const_iterator1 it1e (begin (it2e, iterator2_tag ()));
                typename E::const_iterator1 it1e_end (end (it2e, iterator2_tag ()));
#endif
                if (it1 != it1_end && it1e != it1e_end) {
                    size_type it1_index = it1.index1 (), it1e_index = it1e.index1 ();
                    while (true) {
                        difference_type compare = it1_index - it1e_index;
                        if (compare == 0) {
                            functor1_type::apply (*it1, *it1e);
                            ++ it1, ++ it1e;
                            if (it1 != it1_end && it1e != it1e_end) {
                                it1_index = it1.index1 ();
                                it1e_index = it1e.index1 ();
                            } else
                                break;
                        } else if (compare < 0) {
                            if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
                                functor1_type::apply (*it1, value_type (0));
                                ++ it1;
                            } else
                                increment (it1, it1_end, - compare);
                            if (it1 != it1_end)
                                it1_index = it1.index1 ();
                            else
                                break;
                        } else if (compare > 0) {
                            increment (it1e, it1e_end, compare);
                            if (it1e != it1e_end)
                                it1e_index = it1e.index1 ();
                            else
                                break;
                        }
                    }
                }
                if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
                    while (it1 != it1_end) {
                        functor1_type::apply (*it1, value_type (0));
                        ++ it1;
                    }
                } else {
                    it1 = it1_end;
                }
                ++ it2, ++ it2e;
            } else if (compare < 0) {
                if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                    typename M::iterator1 it1 (it2.begin ());
                    typename M::iterator1 it1_end (it2.end ());
#else
                    typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
                    typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
#endif
                    while (it1 != it1_end) {
                        functor1_type::apply (*it1, value_type (0));
                        ++ it1;
                    }
                    ++ it2;
                } else {
                    increment (it2, it2_end, - compare);
                }
            } else if (compare > 0) {
                increment (it2e, it2e_end, compare);
            }
        }
        if (boost::is_same<BOOST_UBLAS_TYPENAME functor1_type::assign_category, assign_tag>::value) {
            while (it2 != it2_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename M::iterator1 it1 (it2.begin ());
                typename M::iterator1 it1_end (it2.end ());
#else
                typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
                typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
#endif
                while (it1 != it1_end) {
                    functor1_type::apply (*it1, value_type (0));
                    ++ it1;
                }
                ++ it2;
            }
        } else {
            it2 = it2_end;
        }
#if BOOST_UBLAS_TYPE_CHECK
        if (! disable_type_check<bool>::value)
            BOOST_UBLAS_CHECK (equals (m, cm), external_logic ());
#endif
    }

    // Dispatcher
    template<class F, class M, class E>
    BOOST_UBLAS_INLINE
    void matrix_assign (F, M &m, const matrix_expression<E> &e) {
        typedef F functor_type;
        typedef typename matrix_assign_traits<BOOST_UBLAS_TYPENAME M::storage_category,
                                              BOOST_UBLAS_TYPENAME F::assign_category,
                                              BOOST_UBLAS_TYPENAME E::const_iterator1::iterator_category,
                                              BOOST_UBLAS_TYPENAME E::const_iterator2::iterator_category>::storage_category storage_category;
        // give preference to expressions orientation if known
        typedef typename boost::mpl::if_<boost::is_same<BOOST_UBLAS_TYPENAME E::orientation_category, unknown_orientation_tag>,
                                          BOOST_UBLAS_TYPENAME M::orientation_category ,
                                          BOOST_UBLAS_TYPENAME E::orientation_category >::type orientation_category;
        matrix_assign (functor_type (), full (), m, e, storage_category (), orientation_category ());
    }
    template<class F1, class F2, class M, class E>
    BOOST_UBLAS_INLINE
    void matrix_assign (F1, F2, M &m, const matrix_expression<E> &e) {
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef typename matrix_assign_traits<BOOST_UBLAS_TYPENAME M::storage_category,
                                              BOOST_UBLAS_TYPENAME F1::assign_category,
                                              BOOST_UBLAS_TYPENAME E::const_iterator1::iterator_category,
                                              BOOST_UBLAS_TYPENAME E::const_iterator2::iterator_category>::storage_category storage_category;
        // give preference to expressions orientation if known
        typedef typename boost::mpl::if_<boost::is_same<BOOST_UBLAS_TYPENAME E::orientation_category, unknown_orientation_tag>,
                                          BOOST_UBLAS_TYPENAME M::orientation_category ,
                                          BOOST_UBLAS_TYPENAME E::orientation_category >::type orientation_category;
        matrix_assign (functor1_type (), functor2_type (), m, e, storage_category (), orientation_category ());
    }

    template<class LS, class RI1, class RI2>
    struct matrix_swap_traits {
        typedef LS storage_category;
    };

    template<>
    struct matrix_swap_traits<dense_proxy_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    template<>
    struct matrix_swap_traits<packed_proxy_tag, sparse_bidirectional_iterator_tag, sparse_bidirectional_iterator_tag> {
        typedef sparse_proxy_tag storage_category;
    };

    // Dense (proxy) row major case
    template<class F, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_swap (F, full, M &m, matrix_expression<E> &e, dense_proxy_tag, row_major_tag) {
        typedef F functor_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typename M::iterator1 it1 (m.begin1 ());
        typename E::iterator1 it1e (e ().begin1 ());
        difference_type size1 (BOOST_UBLAS_SAME (m.size1 (), size_type (e ().end1 () - it1e)));
        while (-- size1 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator2 it2 (it1.begin ());
            typename E::iterator2 it2e (it1e.begin ());
            difference_type size2 (BOOST_UBLAS_SAME (m.size2 (), size_type (it1e.end () - it2e)));
#else
            typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
            typename E::iterator2 it2e (begin (it1e, iterator1_tag ()));
            difference_type size2 (BOOST_UBLAS_SAME (m.size2 (), size_type (end (it1e, iterator1_tag ()) - it2e)));
#endif
            while (-- size2 >= 0)
                functor_type::apply (*it2, *it2e), ++ it2, ++ it2e;
            ++ it1, ++ it1e;
        }
    }
    // Dense (proxy) column major case
    template<class F, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_swap (F, full, M &m, matrix_expression<E> &e, dense_proxy_tag, column_major_tag) {
        typedef F functor_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typename M::iterator2 it2 (m.begin2 ());
        typename E::iterator2 it2e (e ().begin2 ());
        difference_type size2 (BOOST_UBLAS_SAME (m.size2 (), size_type (e ().end2 () - it2e)));
        while (-- size2 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator1 it1 (it2.begin ());
            typename E::iterator1 it1e (it2e.begin ());
            difference_type size1 (BOOST_UBLAS_SAME (m.size1 (), size_type (it2e.end () - it1e)));
#else
            typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
            typename E::iterator1 it1e (begin (it2e, iterator2_tag ()));
            difference_type size1 (BOOST_UBLAS_SAME (m.size1 (), size_type (end (it2e, iterator2_tag ()) - it1e)));
#endif
            while (-- size1 >= 0)
                functor_type::apply (*it1, *it1e), ++ it1, ++ it1e;
            ++ it2, ++ it2e;
        }
    }
    // Packed (proxy) row major case
    template<class F1, class F2, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_swap (F1, F2, M &m, matrix_expression<E> &e, packed_proxy_tag, row_major_tag) {
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typename M::iterator1 it1 (m.begin1 ());
        typename E::iterator1 it1e (e ().begin1 ());
        difference_type size1 (BOOST_UBLAS_SAME (m.end1 () - it1, e ().end1 () - it1e));
        while (-- size1 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator2 it2 (it1.begin ());
            typename E::iterator2 it2e (it1e.begin ());
            difference_type size2 (BOOST_UBLAS_SAME (it1.end () - it2, it1e.end () - it2e));
#else
            typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
            typename E::iterator2 it2e (begin (it1e, iterator1_tag ()));
            difference_type size2 (BOOST_UBLAS_SAME (end (it1, iterator1_tag ()) - it2, end (it1e, iterator1_tag ()) - it2e));
#endif
            while (-- size2 >= 0)
                functor1_type::apply (*it2, *it2e), ++ it2, ++ it2e;
            ++ it1, ++ it1e;
        }
    }
    // Packed (proxy) column major case
    template<class F1, class F2, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_swap (F1, F2, M &m, matrix_expression<E> &e, packed_proxy_tag, column_major_tag) {
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typename M::iterator2 it2 (m.begin2 ());
        typename E::iterator2 it2e (e ().begin2 ());
        difference_type size2 (BOOST_UBLAS_SAME (m.end2 () - it2, e ().end2 () - it2e));
        while (-- size2 >= 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator1 it1 (it2.begin ());
            typename E::iterator1 it1e (it2e.begin ());
            difference_type size1 (BOOST_UBLAS_SAME (it2.end () - it1, it2e.end () - it1e));
#else
            typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
            typename E::iterator1 it1e (begin (it2e, iterator2_tag ()));
            difference_type size1 (BOOST_UBLAS_SAME (end (it2, iterator2_tag ()) - it1, end (it2e, iterator2_tag ()) - it1e));
#endif
            while (-- size1 >= 0)
                functor1_type::apply (*it1, *it1e), ++ it1, ++ it1e;
            ++ it2, ++ it2e;
        }
    }
    // Sparse (proxy) row major case
    template<class F1, class F2, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_swap (F1, F2, M &m, matrix_expression<E> &e, sparse_proxy_tag, row_major_tag) {
        BOOST_UBLAS_CHECK (m.size1 () == e ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (m.size2 () == e ().size2 (), bad_size ());
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
#ifdef BOOST_UBLAS_NON_CONFORMANT_PROXIES
        make_conformant (m, e, row_major_tag (), functor2_type ());
        make_conformant (e (), m, row_major_tag (), functor2_type ());
#endif
        typename M::iterator1 it1 (m.begin1 ());
        typename M::iterator1 it1_end (m.end1 ());
        typename E::iterator1 it1e (e ().begin1 ());
        typename E::iterator1 it1e_end (e ().end1 ());
        while (it1 != it1_end && it1e != it1e_end) {
            difference_type compare = it1.index1 () - it1e.index1 ();
            if (compare == 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename M::iterator2 it2 (it1.begin ());
                typename M::iterator2 it2_end (it1.end ());
                typename E::iterator2 it2e (it1e.begin ());
                typename E::iterator2 it2e_end (it1e.end ());
#else
                typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
                typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
                typename E::iterator2 it2e (begin (it1e, iterator1_tag ()));
                typename E::iterator2 it2e_end (end (it1e, iterator1_tag ()));
#endif
                if (it2 != it2_end && it2e != it2e_end) {
                    size_type it2_index = it2.index2 (), it2e_index = it2e.index2 ();
                    while (true) {
                        difference_type compare = it2_index - it2e_index;
                        if (compare == 0) {
                            functor1_type::apply (*it2, *it2e);
                            ++ it2, ++ it2e;
                            if (it2 != it2_end && it2e != it2e_end) {
                                it2_index = it2.index2 ();
                                it2e_index = it2e.index2 ();
                            } else
                                break;
                        } else if (compare < 0) {
                            increment (it2, it2_end, - compare);
                            if (it2 != it2_end)
                                it2_index = it2.index2 ();
                            else
                                break;
                        } else if (compare > 0) {
                            increment (it2e, it2e_end, compare);
                            if (it2e != it2e_end)
                                it2e_index = it2e.index2 ();
                            else
                                break;
                        }
                    }
                }
#if BOOST_UBLAS_TYPE_CHECK
                increment (it2e, it2e_end);
                increment (it2, it2_end);
#endif
                ++ it1, ++ it1e;
            } else if (compare < 0) {
#if BOOST_UBLAS_TYPE_CHECK
                while (it1.index1 () < it1e.index1 ()) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                    typename M::iterator2 it2 (it1.begin ());
                    typename M::iterator2 it2_end (it1.end ());
#else
                    typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
                    typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
#endif
                    increment (it2, it2_end);
                    ++ it1;
                }
#else
                increment (it1, it1_end, - compare);
#endif
            } else if (compare > 0) {
#if BOOST_UBLAS_TYPE_CHECK
                while (it1e.index1 () < it1.index1 ()) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                    typename E::iterator2 it2e (it1e.begin ());
                    typename E::iterator2 it2e_end (it1e.end ());
#else
                    typename E::iterator2 it2e (begin (it1e, iterator1_tag ()));
                    typename E::iterator2 it2e_end (end (it1e, iterator1_tag ()));
#endif
                    increment (it2e, it2e_end);
                    ++ it1e;
                }
#else
                increment (it1e, it1e_end, compare);
#endif
            }
        }
#if BOOST_UBLAS_TYPE_CHECK
        while (it1e != it1e_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename E::iterator2 it2e (it1e.begin ());
            typename E::iterator2 it2e_end (it1e.end ());
#else
            typename E::iterator2 it2e (begin (it1e, iterator1_tag ()));
            typename E::iterator2 it2e_end (end (it1e, iterator1_tag ()));
#endif
            increment (it2e, it2e_end);
            ++ it1e;
        }
        while (it1 != it1_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator2 it2 (it1.begin ());
            typename M::iterator2 it2_end (it1.end ());
#else
            typename M::iterator2 it2 (begin (it1, iterator1_tag ()));
            typename M::iterator2 it2_end (end (it1, iterator1_tag ()));
#endif
            increment (it2, it2_end);
            ++ it1;
        }
#endif
    }
    // Sparse (proxy) column major case
    template<class F1, class F2, class M, class E>
    // This function seems to be big. So we do not let the compiler inline it.
    // BOOST_UBLAS_INLINE
    void matrix_swap (F1, F2, M &m, matrix_expression<E> &e, sparse_proxy_tag, column_major_tag) {
        BOOST_UBLAS_CHECK (m.size1 () == e ().size1 (), bad_size ());
        BOOST_UBLAS_CHECK (m.size2 () == e ().size2 (), bad_size ());
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
#ifdef BOOST_UBLAS_NON_CONFORMANT_PROXIES
        make_conformant (m, e, column_major_tag (), functor2_type ());
        make_conformant (e (), m, column_major_tag (), functor2_type ());
#endif
        typename M::iterator2 it2 (m.begin2 ());
        typename M::iterator2 it2_end (m.end2 ());
        typename E::iterator2 it2e (e ().begin2 ());
        typename E::iterator2 it2e_end (e ().end2 ());
        while (it2 != it2_end && it2e != it2e_end) {
            difference_type compare = it2.index2 () - it2e.index2 ();
            if (compare == 0) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                typename M::iterator1 it1 (it2.begin ());
                typename M::iterator1 it1_end (it2.end ());
                typename E::iterator1 it1e (it2e.begin ());
                typename E::iterator1 it1e_end (it2e.end ());
#else
                typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
                typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
                typename E::iterator1 it1e (begin (it2e, iterator2_tag ()));
                typename E::iterator1 it1e_end (end (it2e, iterator2_tag ()));
#endif
                if (it1 != it1_end && it1e != it1e_end) {
                    size_type it1_index = it1.index1 (), it1e_index = it1e.index1 ();
                    while (true) {
                        difference_type compare = it1_index - it1e_index;
                        if (compare == 0) {
                            functor1_type::apply (*it1, *it1e);
                            ++ it1, ++ it1e;
                            if (it1 != it1_end && it1e != it1e_end) {
                                it1_index = it1.index1 ();
                                it1e_index = it1e.index1 ();
                            } else
                                break;
                        }  else if (compare < 0) {
                            increment (it1, it1_end, - compare);
                            if (it1 != it1_end)
                                it1_index = it1.index1 ();
                            else
                                break;
                        } else if (compare > 0) {
                            increment (it1e, it1e_end, compare);
                            if (it1e != it1e_end)
                                it1e_index = it1e.index1 ();
                            else
                                break;
                        }
                    }
                }
#if BOOST_UBLAS_TYPE_CHECK
                increment (it1e, it1e_end);
                increment (it1, it1_end);
#endif
                ++ it2, ++ it2e;
            } else if (compare < 0) {
#if BOOST_UBLAS_TYPE_CHECK
                while (it2.index2 () < it2e.index2 ()) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                    typename M::iterator1 it1 (it2.begin ());
                    typename M::iterator1 it1_end (it2.end ());
#else
                    typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
                    typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
#endif
                    increment (it1, it1_end);
                    ++ it2;
                }
#else
                increment (it2, it2_end, - compare);
#endif
            } else if (compare > 0) {
#if BOOST_UBLAS_TYPE_CHECK
                while (it2e.index2 () < it2.index2 ()) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
                    typename E::iterator1 it1e (it2e.begin ());
                    typename E::iterator1 it1e_end (it2e.end ());
#else
                    typename E::iterator1 it1e (begin (it2e, iterator2_tag ()));
                    typename E::iterator1 it1e_end (end (it2e, iterator2_tag ()));
#endif
                    increment (it1e, it1e_end);
                    ++ it2e;
                }
#else
                increment (it2e, it2e_end, compare);
#endif
            }
        }
#if BOOST_UBLAS_TYPE_CHECK
        while (it2e != it2e_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename E::iterator1 it1e (it2e.begin ());
            typename E::iterator1 it1e_end (it2e.end ());
#else
            typename E::iterator1 it1e (begin (it2e, iterator2_tag ()));
            typename E::iterator1 it1e_end (end (it2e, iterator2_tag ()));
#endif
            increment (it1e, it1e_end);
            ++ it2e;
        }
        while (it2 != it2_end) {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            typename M::iterator1 it1 (it2.begin ());
            typename M::iterator1 it1_end (it2.end ());
#else
            typename M::iterator1 it1 (begin (it2, iterator2_tag ()));
            typename M::iterator1 it1_end (end (it2, iterator2_tag ()));
#endif
            increment (it1, it1_end);
            ++ it2;
        }
#endif
    }

    // Dispatcher
    template<class F, class M, class E>
    BOOST_UBLAS_INLINE
    void matrix_swap (F, M &m, matrix_expression<E> &e) {
        typedef F functor_type;
        typedef typename matrix_swap_traits<BOOST_UBLAS_TYPENAME M::storage_category,
                                            BOOST_UBLAS_TYPENAME E::const_iterator1::iterator_category,
                                            BOOST_UBLAS_TYPENAME E::const_iterator2::iterator_category>::storage_category storage_category;
        // FIXME: can't we improve the dispatch here?
        // typedef typename E::orientation_category orientation_category;
        typedef typename M::orientation_category orientation_category;
        matrix_swap (functor_type (), full (), m, e, storage_category (), orientation_category ());
    }
    template<class F1, class F2, class M, class E>
    BOOST_UBLAS_INLINE
    void matrix_swap (F1, F2, M &m, matrix_expression<E> &e) {
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef typename matrix_swap_traits<BOOST_UBLAS_TYPENAME M::storage_category,
                                            BOOST_UBLAS_TYPENAME E::const_iterator1::iterator_category,
                                            BOOST_UBLAS_TYPENAME E::const_iterator2::iterator_category>::storage_category storage_category;
        // FIXME: can't we improve the dispatch here?
        // typedef typename E::orientation_category orientation_category;
        typedef typename M::orientation_category orientation_category;
        matrix_swap (functor1_type (), functor2_type (), m, e, storage_category (), orientation_category ());
    }

}}}

#endif
