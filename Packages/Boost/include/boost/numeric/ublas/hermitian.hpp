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

#ifndef BOOST_UBLAS_HERMITIAN_H
#define BOOST_UBLAS_HERMITIAN_H

#include <boost/numeric/ublas/config.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Iterators based on ideas of Jeremy Siek
// Hermitian matrices are square. Thanks to Peter Schmitteckert for spotting this.

namespace boost { namespace numeric { namespace ublas {

    template<class M>
    bool is_hermitian (const M &m) {
        typedef typename M::size_type size_type;

        if (m.size1 () != m.size2 ())
            return false;
        size_type size = BOOST_UBLAS_SAME (m.size1 (), m.size2 ());
        for (size_type i = 0; i < size; ++ i) {
            for (size_type j = i; j < size; ++ j) {
                if (m (i, j) != conj (m (j, i)))
                    return false;
            }
        }
        return true;
    }

#ifdef BOOST_UBLAS_STRICT_HERMITIAN

    template<class M>
    class hermitian_matrix_element:
       public container_reference<M> {
    public:
        typedef M matrix_type;
        typedef typename M::size_type size_type;
        typedef typename M::value_type value_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef value_type *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        hermitian_matrix_element (matrix_type &m, size_type i, size_type j, value_type d):
            container_reference<matrix_type> (m), i_ (i), j_ (j), d_ (d), dirty_ (false) {}
        BOOST_UBLAS_INLINE
        hermitian_matrix_element (const hermitian_matrix_element &p):
            container_reference<matrix_type> (p), i_ (p.i_), d_ (p.d_), dirty_ (p.dirty_) {}
        BOOST_UBLAS_INLINE
        ~hermitian_matrix_element () {
            if (dirty_)
                ((*this) ()).at (i_, j_, d_);
        }

        // Assignment
        BOOST_UBLAS_INLINE
        hermitian_matrix_element &operator = (const hermitian_matrix_element &p) {
            // Overide the implict copy assignment
            d_ = p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        hermitian_matrix_element &operator = (const D &d) {
            d_ = d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        hermitian_matrix_element &operator += (const D &d) {
            d_ += d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        hermitian_matrix_element &operator -= (const D &d) {
            d_ -= d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        hermitian_matrix_element &operator *= (const D &d) {
            d_ *= d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        hermitian_matrix_element &operator /= (const D &d) {
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
        void swap (hermitian_matrix_element p) {
            if (this != &p) {
                dirty_ = true;
                p.dirty_ = true;
                std::swap (d_, p.d_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (hermitian_matrix_element p1, hermitian_matrix_element p2) {
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
    struct type_traits<hermitian_matrix_element<M> > {
        typedef typename M::value_type element_type;
        typedef type_traits<hermitian_matrix_element<M> > self_type;
        typedef typename type_traits<element_type>::value_type value_type;
        typedef typename type_traits<element_type>::const_reference const_reference;
        typedef hermitian_matrix_element<M> reference;
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
    struct promote_traits<hermitian_matrix_element<M1>, T2> {
        typedef typename promote_traits<typename hermitian_matrix_element<M1>::value_type, T2>::promote_type promote_type;
    };
    template<class T1, class M2>
    struct promote_traits<T1, hermitian_matrix_element<M2> > {
        typedef typename promote_traits<T1, typename hermitian_matrix_element<M2>::value_type>::promote_type promote_type;
    };
    template<class M1, class M2>
    struct promote_traits<hermitian_matrix_element<M1>, hermitian_matrix_element<M2> > {
        typedef typename promote_traits<typename hermitian_matrix_element<M1>::value_type,
                                        typename hermitian_matrix_element<M2>::value_type>::promote_type promote_type;
    };

#endif

    // Array based hermitian matrix class
    template<class T, class F1, class F2, class A>
    class hermitian_matrix:
        public matrix_expression<hermitian_matrix<T, F1, F2, A> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<hermitian_matrix<T, F1, F2, A> >::operator ();
#endif
        typedef typename A::size_type size_type;
        typedef typename A::difference_type difference_type;
        typedef T value_type;
        // FIXME: no better way to not return the address of a temporary?
        // typedef const T &const_reference;
        typedef const T const_reference;
#ifndef BOOST_UBLAS_STRICT_HERMITIAN
        typedef T &reference;
#else
        typedef hermitian_matrix_element<hermitian_matrix<T, F1, F2, A> > reference;
#endif
        typedef A array_type;
    private:
        typedef T &true_reference;
        typedef T *pointer;
        typedef F1 functor1_type;
        typedef F2 functor2_type;
        typedef hermitian_matrix<T, F1, F2, A> self_type;
    public:
#ifndef BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
        typedef const matrix_const_reference<const self_type> const_closure_type;
#else
        typedef const matrix_reference<const self_type> const_closure_type;
#endif
        typedef matrix_reference<self_type> closure_type;
        typedef vector<T, A> vector_temporary_type;
        typedef matrix<T, F2, A> matrix_temporary_type;  // general sub-matrix
        typedef packed_tag storage_category;
        typedef typename F1::packed_category packed_category;
        typedef typename F2::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        hermitian_matrix ():
            matrix_expression<self_type> (),
            size_ (0), data_ (0) {}
        BOOST_UBLAS_INLINE
        hermitian_matrix (size_type size):
            matrix_expression<self_type> (),
            size_ (BOOST_UBLAS_SAME (size, size)), data_ (functor1_type::packed_size (size, size)) {
        }
        BOOST_UBLAS_INLINE
        hermitian_matrix (size_type size1, size_type size2):
            matrix_expression<self_type> (),
            size_ (BOOST_UBLAS_SAME (size1, size2)), data_ (functor1_type::packed_size (size1, size2)) {
        }
        BOOST_UBLAS_INLINE
        hermitian_matrix (size_type size, const array_type &data):
            matrix_expression<self_type> (),
            size_ (size), data_ (data) {}
        BOOST_UBLAS_INLINE
        hermitian_matrix (const hermitian_matrix &m):
            matrix_expression<self_type> (),
            size_ (m.size_), data_ (m.data_) {}
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_matrix (const matrix_expression<AE> &ae):
            matrix_expression<self_type> (),
            size_ (BOOST_UBLAS_SAME (ae ().size1 (), ae ().size2 ())),
            data_ (functor1_type::packed_size (size_, size_)) {
            matrix_assign (scalar_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
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
        const array_type &data () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        array_type &data () {
            return data_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size, bool preserve = true) {
            size_ = size;
            if (preserve) {
                self_type temporary (size_, size_);
                // FIXME use matrix_resize_preserve on conformant compilers
                // detail::matrix_resize_preserve<functor_type> (*this, temporary, size_, size_);
                assign_temporary (temporary);
            }
            else
                data ().resize (functor1_type::packed_size (size_, size_));
        }
        BOOST_UBLAS_INLINE
        void resize (size_type size1, size_type size2, bool preserve = true) {
            resize (BOOST_UBLAS_SAME (size1, size2), preserve);
        }
        BOOST_UBLAS_INLINE
        void resize_packed_preserve (size_type size) {
            size_ = BOOST_UBLAS_SAME (size, size);
            data ().resize (functor1_type::packed_size (size_, size_), value_type (0));
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference at_element (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            BOOST_UBLAS_CHECK (j < size_, bad_index ());
            // if (i == j)
            //    return type_traits<value_type>::real (data () [functor1_type::element (functor2_type (), i, size_, i, size_)]);
            // else
            if (functor1_type::other (i, j))
                return data () [functor1_type::element (functor2_type (), i, size_, j, size_)];
            else
                return type_traits<value_type>::conj (data () [functor1_type::element (functor2_type (), j, size_, i, size_)]);
        }
        BOOST_UBLAS_INLINE
        true_reference at_element (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            BOOST_UBLAS_CHECK (j < size_, bad_index ());
            if (functor1_type::other (i, j))
                return data () [functor1_type::element (functor2_type (), i, size_, j, size_)];
            else {
                external_logic ().raise ();
                return conj_ = type_traits<value_type>::conj (data () [functor1_type::element (functor2_type (), j, size_, i, size_)]);
            }
        }
        BOOST_UBLAS_INLINE
        void at (size_type i, size_type j, value_type t) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            BOOST_UBLAS_CHECK (j < size_, bad_index ());
            // if (i == j)
            //    data () [functor1_type::element (functor2_type (), i, size_, i, size_)] = type_traits<value_type>::real (t);
            // else
            if (functor1_type::other (i, j))
                data () [functor1_type::element (functor2_type (), i, size_, j, size_)] = t;
            else
                data () [functor1_type::element (functor2_type (), j, size_, i, size_)] = type_traits<value_type>::conj (t);
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
        if (functor1_type::other (i, j))
            return reference (*this, i, j, data () [functor1_type::element (functor2_type (), i, size_, j, size_)]);
        else
            return reference (*this, i, j, type_traits<value_type>::conj (data () [functor1_type::element (functor2_type (), j, size_, i, size_)]));
#endif
        }

        // Assignment
        BOOST_UBLAS_INLINE
        hermitian_matrix &operator = (const hermitian_matrix &m) {
            size_ = m.size_;
            data () = m.data ();
            return *this;
        }
        BOOST_UBLAS_INLINE
        hermitian_matrix &assign_temporary (hermitian_matrix &m) {
            swap (m);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_matrix &operator = (const matrix_expression<AE> &ae) {
            // return assign (self_type (ae));
            self_type temporary (ae);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_matrix &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_matrix& operator += (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this + ae));
            self_type temporary (*this + ae);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_matrix &plus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_plus_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_matrix& operator -= (const matrix_expression<AE> &ae) {
            // return assign (self_type (*this - ae));
            self_type temporary (*this - ae);
            return assign_temporary (temporary);
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_matrix &minus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_minus_assign<true_reference, BOOST_UBLAS_TYPENAME AE::value_type> (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        hermitian_matrix& operator *= (const AT &at) {
            // Multiplication is only allowed for real scalars,
            // otherwise the resulting matrix isn't hermitian.
            // Thanks to Peter Schmitteckert for spotting this.
            BOOST_UBLAS_CHECK (type_traits<value_type>::imag (at) == 0, non_real ());
            matrix_assign_scalar (scalar_multiplies_assign<true_reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        hermitian_matrix& operator /= (const AT &at) {
            // Multiplication is only allowed for real scalars,
            // otherwise the resulting matrix isn't hermitian.
            // Thanks to Peter Schmitteckert for spotting this.
            BOOST_UBLAS_CHECK (type_traits<value_type>::imag (at) == 0, non_real ());
            matrix_assign_scalar (scalar_divides_assign<true_reference, AT> (), *this, at);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (hermitian_matrix &m) {
            if (this != &m) {
                std::swap (size_, m.size_);
                data ().swap (m.data ());
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (hermitian_matrix &m1, hermitian_matrix &m2) {
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
            if (functor1_type::other (i, j)) {
                size_type k = functor1_type::element (functor2_type (), i, size_, j, size_);
                BOOST_UBLAS_CHECK (type_traits<value_type>::equals (data () [k], value_type (0)) ||
                                   type_traits<value_type>::equals (data () [k], t), bad_index ());
                // data ().insert (data ().begin () + k, t);
                data () [k] = t;
            } else {
                size_type k = functor1_type::element (functor2_type (), j, size_, i, size_);
                BOOST_UBLAS_CHECK (type_traits<value_type>::equals (data () [k], value_type (0)) ||
                                   type_traits<value_type>::equals (data () [k], type_traits<value_type>::conj (t)), bad_index ());
                // data ().insert (data ().begin () + k, type_traits<value_type>::conj (t));
                data () [k] = type_traits<value_type>::conj (t);
            }
        }
        BOOST_UBLAS_INLINE
        void erase (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            BOOST_UBLAS_CHECK (j < size_, bad_index ());
            if (functor1_type::other (i, j)) {
                size_type k = functor1_type::element (functor2_type (), i, size_, j, size_);
                // data ().erase (data ().begin () + k);
                data () [k] = value_type (0);
            } else {
                size_type k = functor1_type::element (functor2_type (), j, size_, i, size_);
                // data ().erase (data ().begin () + k);
                data () [k] = value_type (0);
            }
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
            public container_const_reference<hermitian_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator1, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename hermitian_matrix::value_type value_type;
            typedef typename hermitian_matrix::difference_type difference_type;
            typedef typename hermitian_matrix::const_reference reference;
            typedef const typename hermitian_matrix::pointer pointer;
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
                return (*this) ().at_element (it1_, it2_);
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
            public container_reference<hermitian_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator1, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename hermitian_matrix::value_type value_type;
            typedef typename hermitian_matrix::difference_type difference_type;
            typedef typename hermitian_matrix::true_reference reference;
            typedef typename hermitian_matrix::pointer pointer;
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
                return (*this) ().at_element (it1_, it2_);
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
            public container_const_reference<hermitian_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               const_iterator2, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
            typedef const_reference reference;
#else
            typedef typename hermitian_matrix::value_type value_type;
            typedef typename hermitian_matrix::difference_type difference_type;
            typedef typename hermitian_matrix::const_reference reference;
            typedef const typename hermitian_matrix::pointer pointer;
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
                return (*this) ().at_element (it1_, it2_);
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
            public container_reference<hermitian_matrix>,
            public random_access_iterator_base<packed_random_access_iterator_tag,
                                               iterator2, value_type> {
        public:
            typedef packed_random_access_iterator_tag iterator_category;
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename hermitian_matrix::value_type value_type;
            typedef typename hermitian_matrix::difference_type difference_type;
            typedef typename hermitian_matrix::true_reference reference;
            typedef typename hermitian_matrix::pointer pointer;
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
                return (*this) ().at_element (it1_, it2_);
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
        static value_type conj_;
    };

    template<class T, class F1, class F2, class A>
    typename hermitian_matrix<T, F1, F2, A>::value_type hermitian_matrix<T, F1, F2, A>::conj_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME hermitian_matrix<T, F1, F2, A>::value_type ()
#endif
    ;

    // Hermitian matrix adaptor class
    template<class M, class F>
    class hermitian_adaptor:
        public matrix_expression<hermitian_adaptor<M, F> > {
    public:
#ifndef BOOST_UBLAS_NO_PROXY_SHORTCUTS
        BOOST_UBLAS_USING matrix_expression<hermitian_adaptor<M, F> >::operator ();
#endif
        typedef const M const_matrix_type;
        typedef M matrix_type;
        typedef F functor_type;
        typedef typename M::size_type size_type;
        typedef typename M::difference_type difference_type;
        typedef typename M::value_type value_type;
#ifndef BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
        // FIXME: no better way to not return the address of a temporary?
        // typedef typename M::const_reference const_reference;
        typedef typename M::value_type const_reference;
#ifndef BOOST_UBLAS_STRICT_HERMITIAN
        typedef typename M::reference reference;
#else
        typedef hermitian_matrix_element<hermitian_adaptor<M, F> > reference;
#endif
#else
        typedef typename M::value_type const_reference;
#ifndef BOOST_UBLAS_STRICT_HERMITIAN
        typedef typename boost::mpl::if_<boost::is_const<M>,
                                          typename M::value_type,
                                          typename M::reference>::type reference;
#else
        typedef typename boost::mpl::if_<boost::is_const<M>,
                                          typename M::value_type,
                                          hermitian_matrix_element<const hermitian_adaptor<M, F> > >::type reference;
#endif
#endif
#ifndef BOOST_UBLAS_CT_PROXY_CLOSURE_TYPEDEFS
        typedef typename M::closure_type matrix_closure_type;
#else
        typedef typename boost::mpl::if_<boost::is_const<M>,
                                          typename M::const_closure_type,
                                          typename M::closure_type>::type matrix_closure_type;
#endif
    private:
        typedef hermitian_adaptor<M, F> self_type;
    public:
        typedef const self_type const_closure_type;
        typedef self_type closure_type;
        typedef typename M::vector_temporary_type vector_temporary_type;
        typedef typename M::matrix_temporary_type matrix_temporary_type;
        typedef typename storage_restrict_traits<typename M::storage_category,
                                                 packed_proxy_tag>::storage_category storage_category;
        typedef typename F::packed_category packed_category;
        typedef typename M::orientation_category orientation_category;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        hermitian_adaptor ():
            matrix_expression<self_type> (),
            data_ (nil_) {
            BOOST_UBLAS_CHECK (data_.size1 () == data_.size2 (), bad_size ());
        }
        BOOST_UBLAS_INLINE
        hermitian_adaptor (matrix_type &data):
            matrix_expression<self_type> (),
            data_ (data) {
            BOOST_UBLAS_CHECK (data_.size1 () == data_.size2 (), bad_size ());
        }
        BOOST_UBLAS_INLINE
        hermitian_adaptor (const hermitian_adaptor &m):
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

        // Element access
#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
            // if (i == j)
            //     return type_traits<value_type>::real (data () (i, i));
            // else
            if (functor_type::other (i, j))
                return data () (i, j);
            else
                return type_traits<value_type>::conj (data () (j, i));
        }
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
#ifndef BOOST_UBLAS_STRICT_HERMITIAN
            if (functor_type::other (i, j))
                return data () (i, j);
            else {
                external_logic ().raise ();
                return conj_ = type_traits<value_type>::conj (data () (j, i));
            }
#else
            if (functor_type::other (i, j))
                return reference (*this, i, j, data () (i, j));
            else
                return reference (*this, i, j, type_traits<value_type>::conj (data () (j, i)));
#endif
        }
        BOOST_UBLAS_INLINE
        void at (size_type i, size_type j, value_type t) {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
            // if (i == j)
            //     data () (i, i) = type_traits<value_type>::real (t);
            // else
            if (functor_type::other (i, j))
                data () (i, j) = t;
            else
                data () (j, i) = type_traits<value_type>::conj (t);
        }
#else
        BOOST_UBLAS_INLINE
        reference operator () (size_type i, size_type j) const {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
#ifndef BOOST_UBLAS_STRICT_HERMITIAN
            if (functor_type::other (i, j))
                return data () (i, j);
            else {
                external_logic ().raise ();
                return conj_ = type_traits<value_type>::conj (data () (j, i));
            }
#else
            if (functor_type::other (i, j))
                return reference (*this, i, j, data () (i, j));
            else
                return reference (*this, i, j, type_traits<value_type>::conj (data () (j, i)));
#endif
        }
        BOOST_UBLAS_INLINE
        void at (size_type i, size_type j, value_type t) const {
            BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
            BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
            // if (i == j)
            //     data () (i, i) = type_traits<value_type>::real (t);
            // else
            if (functor_type::other (i, j))
                data () (i, j) = t;
            else
                data () (j, i) = type_traits<value_type>::conj (t);
        }
#endif

        // Assignment
        BOOST_UBLAS_INLINE
        hermitian_adaptor &operator = (const hermitian_adaptor &m) {
            matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, value_type> (), functor_type (), *this, m);
            return *this;
        }
        BOOST_UBLAS_INLINE
        hermitian_adaptor &assign_temporary (hermitian_adaptor &m) {
            *this = m;
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_adaptor &operator = (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, value_type> (), functor_type (), *this, matrix<value_type> (ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_adaptor &assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), functor_type (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_adaptor& operator += (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, value_type> (), functor_type (), *this, matrix<value_type> (*this + ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_adaptor &plus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_plus_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), functor_type (), *this, ae);
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_adaptor& operator -= (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, value_type> (), functor_type (), *this, matrix<value_type> (*this - ae));
            return *this;
        }
        template<class AE>
        BOOST_UBLAS_INLINE
        hermitian_adaptor &minus_assign (const matrix_expression<AE> &ae) {
            matrix_assign (scalar_minus_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, BOOST_UBLAS_TYPENAME AE::value_type> (), functor_type (), *this, ae);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        hermitian_adaptor& operator *= (const AT &at) {
            // Multiplication is only allowed for real scalars,
            // otherwise the resulting matrix isn't hermitian.
            // Thanks to Peter Schmitteckert for spotting this.
            BOOST_UBLAS_CHECK (type_traits<value_type>::imag (at) == 0, non_real ());
            matrix_assign_scalar (scalar_multiplies_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, AT> (), *this, at);
            return *this;
        }
        template<class AT>
        BOOST_UBLAS_INLINE
        hermitian_adaptor& operator /= (const AT &at) {
            // Multiplication is only allowed for real scalars,
            // otherwise the resulting matrix isn't hermitian.
            // Thanks to Peter Schmitteckert for spotting this.
            BOOST_UBLAS_CHECK (type_traits<value_type>::imag (at) == 0, non_real ());
            matrix_assign_scalar (scalar_divides_assign<BOOST_UBLAS_TYPENAME iterator1_type::reference, AT> (), *this, at);
            return *this;
        }

        // Closure comparison
        BOOST_UBLAS_INLINE
        bool same_closure (const hermitian_adaptor &ha) const {
            return (*this).data ().same_closure (ha.data ());
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (hermitian_adaptor &m) {
            if (this != &m)
                matrix_swap (scalar_swap<BOOST_UBLAS_TYPENAME iterator1_type::reference, BOOST_UBLAS_TYPENAME iterator1_type::reference> (), functor_type (), *this, m);
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (hermitian_adaptor &m1, hermitian_adaptor &m2) {
            m1.swap (m2);
        }
#endif

        // Iterator types
    private:
        // Use matrix iterator
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
            public container_const_reference<hermitian_adaptor>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               const_iterator1, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename const_iterator1_type::iterator_category,
                                                      dense_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename const_iterator1_type::value_type value_type;
            typedef typename const_iterator1_type::difference_type difference_type;
            // FIXME: no better way to not return the address of a temporary?
            // typedef typename const_iterator1_type::reference reference;
            typedef typename const_iterator1_type::value_type reference;
            typedef typename const_iterator1_type::pointer pointer;
#else
            typedef typename iterator_restrict_traits<typename M::const_iterator1::iterator_category,
                                                      dense_random_access_iterator_tag>::iterator_category iterator_category;
            // FIXME: no better way to not return the address of a temporary?
            // typedef const_reference reference;
            typedef value_type reference;
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
            // FIXME cannot compile
            //  iterator1 does not have these members!
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
                    size_type d = (std::min) (n, it1_end_ - it1_);
                    it1_ += d;
                    n -= d;
                    if (n > 0 || (end_ == 1 && it1_ == it1_end_)) {
                        BOOST_UBLAS_CHECK (end_ == 1, external_logic ());
                        d = (std::min) (n, it2_end_ - it2_begin_);
                        it2_ = it2_begin_ + d;
                        n -= d;
                        current_ = 1;
                    }
                } else /* if (current_ == 1) */ {
                    size_type d = (std::min) (n, it2_end_ - it2_);
                    it2_ += d;
                    n -= d;
                    if (n > 0 || (end_ == 0 && it2_ == it2_end_)) {
                        BOOST_UBLAS_CHECK (end_ == 0, external_logic ());
                        d = (std::min) (n, it1_end_ - it1_begin_);
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
                    size_type d = (std::min) (n, it1_ - it1_begin_);
                    it1_ -= d;
                    n -= d;
                    if (n > 0) {
                        BOOST_UBLAS_CHECK (end_ == 1, external_logic ());
                        d = (std::min) (n, it2_end_ - it2_begin_);
                        it2_ = it2_end_ - d;
                        n -= d;
                        current_ = 1;
                    }
                } else /* if (current_ == 1) */ {
                    size_type d = (std::min) (n, it2_ - it2_begin_);
                    it2_ -= d;
                    n -= d;
                    if (n > 0) {
                        BOOST_UBLAS_CHECK (end_ == 0, external_logic ());
                        d = (std::min) (n, it1_end_ - it1_begin_);
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
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    BOOST_UBLAS_CHECK (it1_ != it1_end_, internal_logic ());
                    if (functor_type::other (index1 (), index2 ()))
                        return *it1_;
                    else
                        return type_traits<value_type>::conj (*it1_);
                } else /* if (current_ == 1) */ {
                    BOOST_UBLAS_CHECK (it2_ != it2_end_, internal_logic ());
                    if (functor_type::other (index1 (), index2 ()))
                        return *it2_;
                    else
                        return type_traits<value_type>::conj (*it2_);
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
            public container_reference<hermitian_adaptor>,
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
            public container_const_reference<hermitian_adaptor>,
            public random_access_iterator_base<dense_random_access_iterator_tag,
                                               const_iterator2, value_type> {
        public:
#ifndef BOOST_MSVC_STD_ITERATOR
            typedef typename iterator_restrict_traits<typename const_iterator2_type::iterator_category,
                                                      dense_random_access_iterator_tag>::iterator_category iterator_category;
            typedef typename const_iterator2_type::value_type value_type;
            typedef typename const_iterator2_type::difference_type difference_type;
            // FIXME: no better way to not return the address of a temporary?
            // typedef typename const_iterator2_type::reference reference;
            typedef typename const_iterator2_type::value_type reference;
            typedef typename const_iterator2_type::pointer pointer;
#else
            typedef typename iterator_restrict_traits<typename M::const_iterator2::iterator_category,
                                                      dense_random_access_iterator_tag>::iterator_category iterator_category;
            // FIXME: no better way to not return the address of a temporary?
            // typedef const_reference reference;
            typedef value_type reference;
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
            // FIXME cannot compiler
            //  iterator2 does not have these members!
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
                    size_type d = (std::min) (n, it1_end_ - it1_);
                    it1_ += d;
                    n -= d;
                    if (n > 0 || (end_ == 1 && it1_ == it1_end_)) {
                        BOOST_UBLAS_CHECK (end_ == 1, external_logic ());
                        d = (std::min) (n, it2_end_ - it2_begin_);
                        it2_ = it2_begin_ + d;
                        n -= d;
                        current_ = 1;
                    }
                } else /* if (current_ == 1) */ {
                    size_type d = (std::min) (n, it2_end_ - it2_);
                    it2_ += d;
                    n -= d;
                    if (n > 0 || (end_ == 0 && it2_ == it2_end_)) {
                        BOOST_UBLAS_CHECK (end_ == 0, external_logic ());
                        d = (std::min) (n, it1_end_ - it1_begin_);
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
                    size_type d = (std::min) (n, it1_ - it1_begin_);
                    it1_ -= d;
                    n -= d;
                    if (n > 0) {
                        BOOST_UBLAS_CHECK (end_ == 1, external_logic ());
                        d = (std::min) (n, it2_end_ - it2_begin_);
                        it2_ = it2_end_ - d;
                        n -= d;
                        current_ = 1;
                    }
                } else /* if (current_ == 1) */ {
                    size_type d = (std::min) (n, it2_ - it2_begin_);
                    it2_ -= d;
                    n -= d;
                    if (n > 0) {
                        BOOST_UBLAS_CHECK (end_ == 0, external_logic ());
                        d = (std::min) (n, it1_end_ - it1_begin_);
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
            const_reference operator * () const {
                BOOST_UBLAS_CHECK (current_ == 0 || current_ == 1, internal_logic ());
                if (current_ == 0) {
                    BOOST_UBLAS_CHECK (it1_ != it1_end_, internal_logic ());
                    if (functor_type::other (index1 (), index2 ()))
                        return *it1_;
                    else
                        return type_traits<value_type>::conj (*it1_);
                } else /* if (current_ == 1) */ {
                    BOOST_UBLAS_CHECK (it2_ != it2_end_, internal_logic ());
                    if (functor_type::other (index1 (), index2 ()))
                        return *it2_;
                    else
                        return type_traits<value_type>::conj (*it2_);
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
            public container_reference<hermitian_adaptor>,
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
        static value_type conj_;
    };

    template<class M, class F>
    typename hermitian_adaptor<M, F>::matrix_type hermitian_adaptor<M, F>::nil_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME hermitian_adaptor<M, F>::matrix_type ()
#endif
    ;
    template<class M, class F>
    typename hermitian_adaptor<M, F>::value_type hermitian_adaptor<M, F>::conj_
#ifdef BOOST_UBLAS_STATIC_OLD_INIT
        = BOOST_UBLAS_TYPENAME hermitian_adaptor<M, F>::value_type ()
#endif
    ;

}}}

#endif
