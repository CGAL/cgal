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

#ifndef BOOST_UBLAS_FWD_H
#define BOOST_UBLAS_FWD_H

namespace boost {

    // Borrowed from Dave Abraham's noncopyable.
    // I believe this should be part of utility.hpp one day...
    class nonassignable {
    protected:
        nonassignable(){}
        ~nonassignable(){}
    private:  // emphasize the following members are private
        const nonassignable& operator=( const nonassignable& );
    }; // nonassignable

}

// Forward declarations
namespace boost { namespace numeric { namespace ublas {

    struct concrete_tag {};
    struct abstract_tag {};

    template<class T>
    class unbounded_array;

    class range;
    class slice;
    template<class A = unbounded_array<std::size_t> >
    class indirect_array;

    template<class I, class T>
    class map_array;

    struct vector_tag {};

    template<class E>
    struct vector_expression;
    template<class E>
    class vector_reference;

    struct matrix_tag {};

    template<class E>
    struct matrix_expression;
    template<class E>
    class matrix_reference;

    template<class E>
    class vector_range;
    template<class E>
    class vector_slice;
    template<class E, class IA = indirect_array<> >
    class vector_indirect;

    template<class E>
    class matrix_row;
    template<class E>
    class matrix_column;
    template<class E>
    class matrix_range;
    template<class E>
    class matrix_slice;
    template<class E, class IA = indirect_array<> >
    class matrix_indirect;

    template<class T, class A = unbounded_array<T> >
    class vector;
    template<class T, std::size_t N>
    class bounded_vector;

    template<class T>
    class unit_vector;

    template<class T>
    class scalar_vector;

    template<class T, std::size_t N>
    class c_vector;

    template<class T, class A = map_array<std::size_t, T> >
    class sparse_vector;

    template<class T, std::size_t IB = 0, class IA = unbounded_array<std::size_t>, class TA = unbounded_array<T> >
    class compressed_vector;

    template<class T, std::size_t IB = 0, class IA = unbounded_array<std::size_t>, class TA = unbounded_array<T> >
    class coordinate_vector;

    struct unknown_orientation_tag {};

    struct row_major_tag {};
    struct row_major;

    struct column_major_tag {};
    struct column_major;

    template<class T, class F = row_major, class A = unbounded_array<T> >
    class matrix;
    template<class T, std::size_t M, std::size_t N, class F = row_major>
    class bounded_matrix;

    template<class T>
    class identity_matrix;

    template<class T>
    class scalar_matrix;

    template<class T, std::size_t M, std::size_t N>
    class c_matrix;

    template<class T, class F = row_major, class A = unbounded_array<unbounded_array<T> > >
    class vector_of_vector;

    template<class T, class F = row_major, class A = unbounded_array<T> >
    class banded_matrix;
    template<class T, class F = row_major, class A = unbounded_array<T> >
    class diagonal_matrix;

    struct lower_tag {};
    struct lower;

    struct upper_tag {};
    struct upper;

    struct unit_lower_tag: public lower_tag {};
    struct unit_lower;

    struct unit_upper_tag: public upper_tag {};
    struct unit_upper;

    template<class T, class F1 = lower, class F2 = row_major, class A = unbounded_array<T> >
    class triangular_matrix;

    template<class M, class F = lower>
    class triangular_adaptor;

    template<class T, class F1 = lower, class F2 = row_major, class A = unbounded_array<T> >
    class symmetric_matrix;

    template<class M, class F = lower>
    class symmetric_adaptor;

    template<class T, class F1 = lower, class F2 = row_major, class A = unbounded_array<T> >
    class hermitian_matrix;

    template<class M, class F = lower>
    class hermitian_adaptor;

    template<class T, class F = row_major, class A = map_array<std::size_t, T> >
    class sparse_matrix;

    template<class T, class F = row_major, class A = map_array<std::size_t, map_array<std::size_t, T> > >
    class sparse_vector_of_sparse_vector;

    template<class T, class F = row_major, std::size_t IB = 0, class IA = unbounded_array<std::size_t>, class TA = unbounded_array<T> >
    class compressed_matrix;

    template<class T, class F = row_major, std::size_t IB = 0, class IA = unbounded_array<std::size_t>, class TA = unbounded_array<T> >
    class coordinate_matrix;

    // Some syntactic sugar. I'd like to drop it ;-)
    template<class V>
    typename V::size_type num_elements (const V &v) {
        return v.size ();
    }
    template<class M>
    typename M::size_type num_rows (const M &m) {
        return m.size1 ();
    }
    template<class M>
    typename M::size_type num_columns (const M &m) {
        return m.size2 ();
    }
    template<class MV>
    typename MV::size_type num_non_zeros (const MV &mv) {
        return mv.non_zeros ();
    }

}}}

#endif




























