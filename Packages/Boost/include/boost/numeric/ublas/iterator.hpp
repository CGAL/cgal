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

#ifndef BOOST_UBLAS_ITERATOR_H
#define BOOST_UBLAS_ITERATOR_H

#include <iterator>

#include <boost/numeric/ublas/exception.hpp>

// Using older GCC the following is missing:
//
// namespace std {
//
//    template <class C, class T, class D = std::ptrdiff_t, class P = T *, class R = T &>
//    struct iterator {
//        typedef C iterator_category;
//        typedef T value_type;
//        typedef D difference_type;
//        typedef P pointer;
//        typedef R reference;
//    };
//
// }
//
// We therefore include the following header
#include <boost/iterator.hpp>
// and use namespace boost instead of std.

namespace boost { namespace numeric { namespace ublas {

  /** \brief Base class of all proxy classes that contain
   *       a (redirectable) reference to an immutable object.
   *
   *       \param C the type of the container referred to
   */
    template<class C>
    class container_const_reference:
        private nonassignable {
    public:
        typedef C container_type;

        BOOST_UBLAS_INLINE
        container_const_reference ():
            c_ (0) {}
        BOOST_UBLAS_INLINE
        container_const_reference (const container_type &c):
            c_ (&c) {}

        BOOST_UBLAS_INLINE
        const container_type &operator () () const {
            return *c_;
        }

        BOOST_UBLAS_INLINE
        container_const_reference &assign (const container_type *c) {
            c_ = c;
            return *this;
        }
        
        // Closure comparison
        BOOST_UBLAS_INLINE
        bool same_closure (const container_const_reference &cr) const {
            return c_ == cr.c_;
        }

    private:
        const container_type *c_;
    };

  /** \brief Base class of all proxy classes that contain
   *         a (redirectable) reference to a mutable object.
   *
   * \param C the type of the container referred to
   */
    template<class C>
    class container_reference:
        private nonassignable {
    public:
        typedef C container_type;

        BOOST_UBLAS_INLINE
        container_reference ():
            c_ (0) {}
        BOOST_UBLAS_INLINE
        container_reference (container_type &c):
            c_ (&c) {}

        BOOST_UBLAS_INLINE
        container_type &operator () () const {
           return *c_;
        }

        BOOST_UBLAS_INLINE
        container_reference &assign (container_type *c) {
            c_ = c;
            return *this;
        }

        // Closure comparison
        BOOST_UBLAS_INLINE
        bool same_closure (const container_reference &cr) const {
            return c_ == cr.c_;
        }

    private:
        container_type *c_;
    };

  /** \brief Base class of all forward iterators.
   * 
   *  \param IC the iterator category
   *  \param I the derived iterator type
   *  \param T the value type
   * 
   * The forward iterator can only proceed in one direction
   * via the post increment operator.
   */
    template<class IC, class I, class T>
    struct forward_iterator_base:
        public boost::iterator<IC, T> {
        typedef I derived_iterator_type;
        typedef T derived_value_type;

        // Arithmetic
        BOOST_UBLAS_INLINE
        derived_iterator_type operator ++ (int) {
            derived_iterator_type &d (*static_cast<const derived_iterator_type *> (this));
            derived_iterator_type tmp (d);
            ++ d;
            return tmp;
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend derived_iterator_type operator ++ (derived_iterator_type &d, int) {
            derived_iterator_type tmp (d);
            ++ d;
            return tmp;
        }
#endif

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator != (const derived_iterator_type &it) const {
            const derived_iterator_type *d = static_cast<const derived_iterator_type *> (this);
            return ! (*d == it);
        }
    };

#ifdef BOOST_UBLAS_NO_MEMBER_FRIENDS
    template<class IC, class I, class T>
    BOOST_UBLAS_INLINE
    typename forward_iterator_base<IC, I, T>::derived_iterator_type operator ++ (forward_iterator_base<IC, I, T> &it, int) {
        typedef BOOST_UBLAS_TYPENAME forward_iterator_base<IC, I, T>::derived_iterator_type derived_iterator_type;
        derived_iterator_type &d (static_cast<derived_iterator_type &> (it));
        derived_iterator_type tmp (d);
        ++ d;
        return tmp;
    }
#endif

  /** \brief Base class of all bidirectional iterators.
   *
   * \param IC the iterator category
   * \param I the derived iterator type
   * \param T the value type
   *
   * The bidirectional iterator can proceed in both directions
   * via the post increment and post decrement operator.
   */
    template<class IC, class I, class T>
    struct bidirectional_iterator_base:
        public boost::iterator<IC, T> {
        typedef I derived_iterator_type;
        typedef T derived_value_type;

        // Arithmetic
        BOOST_UBLAS_INLINE
        derived_iterator_type operator ++ (int) {
            derived_iterator_type &d (*static_cast<const derived_iterator_type *> (this));
            derived_iterator_type tmp (d);
            ++ d;
            return tmp;
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend derived_iterator_type operator ++ (derived_iterator_type &d, int) {
            derived_iterator_type tmp (d);
            ++ d;
            return tmp;
        }
#endif
        BOOST_UBLAS_INLINE
        derived_iterator_type operator -- (int) {
            derived_iterator_type &d (*static_cast<const derived_iterator_type *> (this));
            derived_iterator_type tmp (d);
            -- d;
            return tmp;
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend derived_iterator_type operator -- (derived_iterator_type &d, int) {
            derived_iterator_type tmp (d);
            -- d;
            return tmp;
        }
#endif

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator != (const derived_iterator_type &it) const {
            const derived_iterator_type *d = static_cast<const derived_iterator_type *> (this);
            return ! (*d == it);
        }
    };

#ifdef BOOST_UBLAS_NO_MEMBER_FRIENDS
    template<class IC, class I, class T>
    BOOST_UBLAS_INLINE
    typename bidirectional_iterator_base<IC, I, T>::derived_iterator_type operator ++ (bidirectional_iterator_base<IC, I, T> &it, int) {
        typedef BOOST_UBLAS_TYPENAME bidirectional_iterator_base<IC, I, T>::derived_iterator_type derived_iterator_type;
        derived_iterator_type &d (static_cast<derived_iterator_type &> (it));
        derived_iterator_type tmp (d);
        ++ d;
        return tmp;
    }
    template<class IC, class I, class T>
    BOOST_UBLAS_INLINE
    typename bidirectional_iterator_base<IC, I, T>::derived_iterator_type operator -- (bidirectional_iterator_base<IC, I, T> &it, int) {
        typedef BOOST_UBLAS_TYPENAME bidirectional_iterator_base<IC, I, T>::derived_iterator_type derived_iterator_type;
        derived_iterator_type &d (static_cast<derived_iterator_type &> (it));
        derived_iterator_type tmp (d);
        -- d;
        return tmp;
    }
#endif

  /** \brief Base class of all random access iterators.
   *
   * \param IC the iterator category
   * \param I the derived iterator type
   * \param T the value type
   * \param D the difference type, default: std::ptrdiff_t
   *
   * The random access iterator can proceed in both directions
   * via the post increment/decrement operator or in larger steps
   * via the +, - and +=, -= operators. The random access iterator
   * is LessThan Comparable.
   */
    template<class IC, class I, class T, class D = std::ptrdiff_t>
    // ISSUE the default here seems rather dangerous as it can easlly be (silently) incorrect
    struct random_access_iterator_base:
        public boost::iterator<IC, T> {
        typedef I derived_iterator_type;
        typedef T derived_value_type;
        typedef D derived_difference_type;
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef D difference_type;
#endif

        /*
         *  FIXME Need to explicitly pass derived_refernce_type as otherwise I undefined type or foward declared
        typedef BOOST_UBLAS_TYPENAME derived_iterator_type::reference derived_reference_type;
        // Indexed element
        BOOST_UBLAS_INLINE
        derived_reference_type operator [] (derived_difference_type n) {
            return *(*this + n);
        }
        */

        // Arithmetic
        BOOST_UBLAS_INLINE
        derived_iterator_type operator ++ (int) {
            derived_iterator_type &d (*static_cast<derived_iterator_type *> (this));
            derived_iterator_type tmp (d);
            ++ d;
            return tmp;
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend derived_iterator_type operator ++ (derived_iterator_type &d, int) {
            derived_iterator_type tmp (d);
            ++ d;
            return tmp;
        }
#endif
        BOOST_UBLAS_INLINE
        derived_iterator_type operator -- (int) {
            derived_iterator_type &d (*static_cast<derived_iterator_type *> (this));
            derived_iterator_type tmp (d);
            -- d;
            return tmp;
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend derived_iterator_type operator -- (derived_iterator_type &d, int) {
            derived_iterator_type tmp (d);
            -- d;
            return tmp;
        }
#endif
        BOOST_UBLAS_INLINE
        derived_iterator_type operator + (derived_difference_type n) const {
            derived_iterator_type tmp (*static_cast<const derived_iterator_type *> (this));
            return tmp += n;
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend derived_iterator_type operator + (const derived_iterator_type &d, derived_difference_type n) {
            derived_iterator_type tmp (d);
            return tmp += n;
        }
        BOOST_UBLAS_INLINE
        friend derived_iterator_type operator + (derived_difference_type n, const derived_iterator_type &d) {
            derived_iterator_type tmp (d);
            return tmp += n;
        }
#endif
        BOOST_UBLAS_INLINE
        derived_iterator_type operator - (derived_difference_type n) const {
            derived_iterator_type tmp (*static_cast<const derived_iterator_type *> (this));
            return tmp -= n;
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend derived_iterator_type operator - (const derived_iterator_type &d, derived_difference_type n) {
            derived_iterator_type tmp (d);
            return tmp -= n;
        }
#endif

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator != (const derived_iterator_type &it) const {
            const derived_iterator_type *d = static_cast<const derived_iterator_type *> (this);
            return ! (*d == it);
        }
        BOOST_UBLAS_INLINE
        bool operator <= (const derived_iterator_type &it) const {
            const derived_iterator_type *d = static_cast<const derived_iterator_type *> (this);
            return ! (it < *d);
        }
        BOOST_UBLAS_INLINE
        bool operator >= (const derived_iterator_type &it) const {
            const derived_iterator_type *d = static_cast<const derived_iterator_type *> (this);
            return ! (*d < it);
        }
        BOOST_UBLAS_INLINE
        bool operator > (const derived_iterator_type &it) const {
            const derived_iterator_type *d = static_cast<const derived_iterator_type *> (this);
            return it < *d;
        }
    };

#ifdef BOOST_UBLAS_NO_MEMBER_FRIENDS
    template<class IC, class I, class T>
    BOOST_UBLAS_INLINE
    typename random_access_iterator_base<IC, I, T>::derived_iterator_type operator ++ (random_access_iterator_base<IC, I, T> &it, int) {
        typedef BOOST_UBLAS_TYPENAME random_access_iterator_base<IC, I, T>::derived_iterator_type derived_iterator_type;
        derived_iterator_type &d (static_cast<derived_iterator_type &> (it));
        derived_iterator_type tmp (d);
        ++ d;
        return tmp;
    }
    template<class IC, class I, class T>
    BOOST_UBLAS_INLINE
    typename random_access_iterator_base<IC, I, T>::derived_iterator_type operator -- (random_access_iterator_base<IC, I, T> &it, int) {
        typedef BOOST_UBLAS_TYPENAME random_access_iterator_base<IC, I, T>::derived_iterator_type derived_iterator_type;
        derived_iterator_type &d (static_cast<derived_iterator_type &> (it));
        derived_iterator_type tmp (d);
        -- d;
        return tmp;
    }
    template<class IC, class I, class T>
    BOOST_UBLAS_INLINE
    typename random_access_iterator_base<IC, I, T>::derived_iterator_type operator + (const random_access_iterator_base<IC, I, T> &it, std::ptrdiff_t n) {
        typedef BOOST_UBLAS_TYPENAME random_access_iterator_base<IC, I, T>::derived_iterator_type derived_iterator_type;
        derived_iterator_type tmp (static_cast<const derived_iterator_type &> (it));
        return tmp += n;
    }
    template<class IC, class I, class T>
    BOOST_UBLAS_INLINE
    typename random_access_iterator_base<IC, I, T>::derived_iterator_type operator + (std::ptrdiff_t n, const random_access_iterator_base<IC, I, T> &it) {
        typedef BOOST_UBLAS_TYPENAME random_access_iterator_base<IC, I, T>::derived_iterator_type derived_iterator_type;
        derived_iterator_type tmp (static_cast<const derived_iterator_type &> (it));
        return tmp += n;
    }
    template<class IC, class I, class T>
    BOOST_UBLAS_INLINE
    typename random_access_iterator_base<IC, I, T>::derived_iterator_type operator - (const random_access_iterator_base<IC, I, T> &it, std::ptrdiff_t n) {
        typedef BOOST_UBLAS_TYPENAME random_access_iterator_base<IC, I, T>::derived_iterator_type derived_iterator_type;
        derived_iterator_type tmp (static_cast<const derived_iterator_type &> (it));
        return tmp -= n;
    }
#endif

#ifdef BOOST_MSVC_STD_ITERATOR

  /** \brief Base class of all reverse iterators. (MSVC version)
   *
   * \param I the derived iterator type
   * \param T the value type
   * \param R the reference type
   *
   * The reverse iterator implements a bidirectional iterator
   * reversing the elements of the underlying iterator. It
   * implements most operators of a random access iterator.
   *
   * uBLAS extension: it.index()
   */

    // Renamed this class from reverse_iterator to get
    // typedef reverse_iterator<...> reverse_iterator
    // working. Thanks to Gabriel Dos Reis for explaining this.
    template <class I, class T, class R>
    class reverse_iterator_base:
        public std::reverse_bidirectional_iterator<I, T, R> {
    public:
        typedef typename I::container_type container_type;
        typedef typename container_type::size_type size_type;
        typedef typename I::difference_type difference_type;
        typedef I iterator_type;
        typedef T value_type;
        typedef R reference;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        reverse_iterator_base ():
            std::reverse_bidirectional_iterator<iterator_type, value_type, reference> () {}
        BOOST_UBLAS_INLINE
        reverse_iterator_base (const iterator_type &it):
            std::reverse_bidirectional_iterator<iterator_type, value_type, reference> (it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        reverse_iterator_base &operator += (difference_type n) {
            // Comeau recommends...
            return *this = this->base () - n;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base &operator -= (difference_type n) {
            // Comeau recommends...
            return *this = this->base () + n;
        }

        BOOST_UBLAS_INLINE
        const container_type &operator () () const {
            // Comeau recommends...
            return this->base () ();
        }

        BOOST_UBLAS_INLINE
        size_type index () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).index ();
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator < (const reverse_iterator_base &it) const {
            return ! (this->base () < it.base ());
        }
        BOOST_UBLAS_INLINE
        bool operator <= (const reverse_iterator_base &it) const {
            return ! (this->base () <= it.base ());
        }
        BOOST_UBLAS_INLINE
        bool operator >= (const reverse_iterator_base &it) const {
            return ! (this->base () >= it.base ());
        }
        BOOST_UBLAS_INLINE
        bool operator > (const reverse_iterator_base &it) const {
            return ! (this->base () > it.base ());
        }
    };

    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    reverse_iterator_base<I, T, R> operator + (const reverse_iterator_base<I, T, R> &it, std::ptrdiff_t n) {
        reverse_iterator_base<I, T, R> tmp (it);
        return tmp += n;
    }
    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    reverse_iterator_base<I, T, R> operator + (std::ptrdiff_t n, const reverse_iterator_base<I, T, R> &it) {
        reverse_iterator_base<I, T, R> tmp (it);
        return tmp += n;
    }
    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    reverse_iterator_base<I, T, R> operator - (const reverse_iterator_base<I, T, R> &it, std::ptrdiff_t n) {
        reverse_iterator_base<I, T, R> tmp (it);
        return tmp -= n;
    }
    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    std::ptrdiff_t operator - (const reverse_iterator_base<I, T, R> &it1, const reverse_iterator_base<I, T, R> &it2) {
        return it2.base () - it1.base ();
    }

  /** \brief 1st base class of all matrix reverse iterators. (MSVC version)
   *
   * \param I the derived iterator type
   * \param T the value type
   * \param R the reference type
   *
   * The reverse iterator implements a bidirectional iterator
   * reversing the elements of the underlying iterator. It
   * implements most operators of a random access iterator.
   *
   * uBLAS extension: it.index1(), it.index2() and access to
   * the dual iterator via begin(), end(), rbegin(), rend()
   */

    // Renamed this class from reverse_iterator1 to get
    // typedef reverse_iterator1<...> reverse_iterator1
    // working. Thanks to Gabriel Dos Reis for explaining this.
    template <class I, class T, class R>
    class reverse_iterator_base1:
        public std::reverse_bidirectional_iterator<I, T, R> {
    public:
        typedef typename I::container_type container_type;
        typedef typename container_type::size_type size_type;
        typedef typename I::difference_type difference_type;
        typedef I iterator_type;
        typedef T value_type;
        typedef R reference;
        typedef typename I::dual_iterator_type dual_iterator_type;
        typedef typename I::dual_reverse_iterator_type dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 ():
            std::reverse_bidirectional_iterator<iterator_type, value_type, reference> () {}
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 (const iterator_type &it):
            std::reverse_bidirectional_iterator<iterator_type, value_type, reference> (it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 &operator += (difference_type n) {
            // Comeau recommends...
            return *this = this->base () - n;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 &operator -= (difference_type n) {
            // Comeau recommends...
            return *this = this->base () + n;
        }

        BOOST_UBLAS_INLINE
        const container_type &operator () () const {
            // Comeau recommends...
            return this->base () ();
        }

        BOOST_UBLAS_INLINE
        size_type index1 () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).index1 ();
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).index2 ();
        }

        BOOST_UBLAS_INLINE
        dual_iterator_type begin () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).begin ();
        }
        BOOST_UBLAS_INLINE
        dual_iterator_type end () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).end ();
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rbegin () const {
            return dual_reverse_iterator_type (end ());
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rend () const {
            return dual_reverse_iterator_type (begin ());
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator < (const reverse_iterator_base1 &it) const {
            return ! (this->base () < it.base ());
        }
        BOOST_UBLAS_INLINE
        bool operator <= (const reverse_iterator_base1 &it) const {
            return ! (this->base () <= it.base ());
        }
        BOOST_UBLAS_INLINE
        bool operator >= (const reverse_iterator_base1 &it) const {
            return ! (this->base () >= it.base ());
        }
        BOOST_UBLAS_INLINE
        bool operator > (const reverse_iterator_base1 &it) const {
            return ! (this->base () > it.base ());
        }
    };

    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    reverse_iterator_base1<I, T, R> operator + (const reverse_iterator_base1<I, T, R> &it, std::ptrdiff_t n) {
        reverse_iterator_base1<I, T, R> tmp (it);
        return tmp += n;
    }
    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    reverse_iterator_base1<I, T, R> operator + (std::ptrdiff_t n, const reverse_iterator_base1<I, T, R> &it) {
        reverse_iterator_base1<I, T, R> tmp (it);
        return tmp += n;
    }
    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    reverse_iterator_base1<I, T, R> operator - (const reverse_iterator_base1<I, T, R> &it, std::ptrdiff_t n) {
        reverse_iterator_base1<I, T, R> tmp (it);
        return tmp -= n;
    }
    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    std::ptrdiff_t operator - (const reverse_iterator_base1<I, T, R> &it1, const reverse_iterator_base1<I, T, R> &it2) {
        return it2.base () - it1.base ();
    }

  /** \brief 2nd base class of all matrix reverse iterators. (MSVC version)
   *
   * \param I the derived iterator type
   * \param T the value type
   * \param R the reference type
   *
   * The reverse iterator implements a bidirectional iterator
   * reversing the elements of the underlying iterator. It
   * implements most operators of a random access iterator.
   *
   * uBLAS extension: it.index1(), it.index2() and access to
   * the dual iterator via begin(), end(), rbegin(), rend()
   *
   * Note: This class is _identical_ to reverse_iterator_base1
   */

    // Renamed this class from reverse_iterator2 to get
    // typedef reverse_iterator2<...> reverse_iterator2
    // working. Thanks to Gabriel Dos Reis for explaining this.
    template <class I, class T, class R>
    class reverse_iterator_base2:
        public std::reverse_bidirectional_iterator<I, T, R> {
    public:
        typedef typename I::container_type container_type;
        typedef typename container_type::size_type size_type;
        typedef typename I::difference_type difference_type;
        typedef I iterator_type;
        typedef T value_type;
        typedef R reference;
        typedef typename I::dual_iterator_type dual_iterator_type;
        typedef typename I::dual_reverse_iterator_type dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 ():
            std::reverse_bidirectional_iterator<iterator_type, value_type, reference> () {}
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 (const iterator_type &it):
            std::reverse_bidirectional_iterator<iterator_type, value_type, reference> (it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 &operator += (difference_type n) {
            // Comeau recommends...
            return *this = this->base () - n;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 &operator -= (difference_type n) {
            // Comeau recommends...
            return *this = this->base () + n;
        }

        BOOST_UBLAS_INLINE
        const container_type &operator () () const {
            // Comeau recommends...
            return this->base () ();
        }

        BOOST_UBLAS_INLINE
        size_type index1 () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).index1 ();
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).index2 ();
        }

        BOOST_UBLAS_INLINE
        dual_iterator_type begin () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).begin ();
        }
        BOOST_UBLAS_INLINE
        dual_iterator_type end () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).end ();
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rbegin () const {
            return dual_reverse_iterator_type (end ());
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rend () const {
            return dual_reverse_iterator_type (begin ());
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator < (const reverse_iterator_base2 &it) const {
            return ! (this->base () < it.base ());
        }
        BOOST_UBLAS_INLINE
        bool operator <= (const reverse_iterator_base2 &it) const {
            return ! (this->base () <= it.base ());
        }
        BOOST_UBLAS_INLINE
        bool operator >= (const reverse_iterator_base2 &it) const {
            return ! (this->base () >= it.base ());
        }
        BOOST_UBLAS_INLINE
        bool operator > (const reverse_iterator_base2 &it) const {
            return ! (this->base () > it.base ());
        }
    };

    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    reverse_iterator_base2<I, T, R> operator + (const reverse_iterator_base2<I, T, R> &it, std::ptrdiff_t n) {
        reverse_iterator_base2<I, T, R> tmp (it);
        return tmp += n;
    }
    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    reverse_iterator_base2<I, T, R> operator + (std::ptrdiff_t n, const reverse_iterator_base2<I, T, R> &it) {
        reverse_iterator_base2<I, T, R> tmp (it);
        return tmp += n;
    }
    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    reverse_iterator_base2<I, T, R> operator - (const reverse_iterator_base2<I, T, R> &it, std::ptrdiff_t n) {
        reverse_iterator_base2<I, T, R> tmp (it);
        return tmp -= n;
    }
    template<class I, class T, class R>
    BOOST_UBLAS_INLINE
    std::ptrdiff_t operator - (const reverse_iterator_base2<I, T, R> &it1, const reverse_iterator_base2<I, T, R> &it2) {
        return it2.base () - it1.base ();
    }

#else

  /** \brief Base class of all reverse iterators. (non-MSVC version)
   *
   * \param I the derived iterator type
   * \param T the value type
   * \param R the reference type
   *
   * The reverse iterator implements a bidirectional iterator
   * reversing the elements of the underlying iterator. It
   * implements most operators of a random access iterator.
   *
   * uBLAS extension: it.index()
   */

    // Renamed this class from reverse_iterator to get
    // typedef reverse_iterator<...> reverse_iterator
    // working. Thanks to Gabriel Dos Reis for explaining this.
    template <class I>
    class reverse_iterator_base:
        public std::reverse_iterator<I> {
    public:
        typedef typename I::container_type container_type;
        typedef typename container_type::size_type size_type;
        typedef typename I::difference_type difference_type;
        typedef I iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        reverse_iterator_base ():
            std::reverse_iterator<iterator_type> () {}
        BOOST_UBLAS_INLINE
        reverse_iterator_base (const iterator_type &it):
            std::reverse_iterator<iterator_type> (it) {}

#ifndef BOOST_UBLAS_NO_REVERSE_ITERATOR_OVERLOADS
        // Arithmetic
        BOOST_UBLAS_INLINE
        reverse_iterator_base &operator ++ () {
            // Comeau recommends...
            return *this = -- this->base ();
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base operator ++ (int) {
            // Comeau recommends...
            reverse_iterator_base tmp (*this);
            *this = -- this->base ();
            return tmp;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base &operator -- () {
            // Comeau recommends...
            return *this = ++ this->base ();
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base operator -- (int) {
            // Comeau recommends...
            reverse_iterator_base tmp (*this);
            *this = ++ this->base ();
            return tmp;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base &operator += (difference_type n) {
            // Comeau recommends...
            return *this = this->base () - n;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base &operator -= (difference_type n) {
            // Comeau recommends...
            return *this = this->base () + n;
        }
#endif

#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend reverse_iterator_base operator + (const reverse_iterator_base &it, difference_type n) {
            reverse_iterator_base tmp (it);
            return tmp += n;
        }
        BOOST_UBLAS_INLINE
        friend reverse_iterator_base operator + (difference_type n, const reverse_iterator_base &it) {
            reverse_iterator_base tmp (it);
            return tmp += n;
        }
        BOOST_UBLAS_INLINE
        friend reverse_iterator_base operator - (const reverse_iterator_base &it, difference_type n) {
            reverse_iterator_base tmp (it);
            return tmp -= n;
        }
        BOOST_UBLAS_INLINE
        friend difference_type operator - (const reverse_iterator_base &it1, const reverse_iterator_base &it2) {
            return it2.base () - it1.base ();
        }
#endif

        BOOST_UBLAS_INLINE
        const container_type &operator () () const {
            return this->base () ();
        }

        BOOST_UBLAS_INLINE
        size_type index () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).index ();
        }
    };

#ifdef BOOST_UBLAS_NO_MEMBER_FRIENDS
    template<class I>
    BOOST_UBLAS_INLINE
    reverse_iterator_base<I> operator + (const reverse_iterator_base<I> &it, std::ptrdiff_t n) {
        reverse_iterator_base<I> tmp (it);
        return tmp += n;
    }
    template<class I>
    BOOST_UBLAS_INLINE
    reverse_iterator_base<I> operator + (std::ptrdiff_t n, const reverse_iterator_base<I> &it) {
        reverse_iterator_base<I> tmp (it);
        return tmp += n;
    }
    template<class I>
    BOOST_UBLAS_INLINE
    reverse_iterator_base<I> operator - (const reverse_iterator_base<I> &it, std::ptrdiff_t n) {
        reverse_iterator_base<I> tmp (it);
        return tmp -= n;
    }
    template<class I>
    BOOST_UBLAS_INLINE
    std::ptrdiff_t operator - (const reverse_iterator_base<I> &it1, const reverse_iterator_base<I> &it2) {
        return it2.base () - it1.base ();
    }
#endif

  /** \brief 1st base class of all matrix reverse iterators. (non-MSVC version)
   *
   * \param I the derived iterator type
   *
   * The reverse iterator implements a bidirectional iterator
   * reversing the elements of the underlying iterator. It
   * implements most operators of a random access iterator.
   *
   * uBLAS extension: it.index1(), it.index2() and access to
   * the dual iterator via begin(), end(), rbegin(), rend()
   */

    // Renamed this class from reverse_iterator1 to get
    // typedef reverse_iterator1<...> reverse_iterator1
    // working. Thanks to Gabriel Dos Reis for explaining this.
    template <class I>
    class reverse_iterator_base1:
        public std::reverse_iterator<I> {
    public:
        typedef typename I::container_type container_type;
        typedef typename container_type::size_type size_type;
        typedef typename I::difference_type difference_type;
        typedef I iterator_type;
        typedef typename I::dual_iterator_type dual_iterator_type;
        typedef typename I::dual_reverse_iterator_type dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 ():
            std::reverse_iterator<iterator_type> () {}
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 (const iterator_type &it):
            std::reverse_iterator<iterator_type> (it) {}

#ifndef BOOST_UBLAS_NO_REVERSE_ITERATOR_OVERLOADS
        // Arithmetic
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 &operator ++ () {
            // Comeau recommends...
            return *this = -- this->base ();
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 operator ++ (int) {
            // Comeau recommends...
            reverse_iterator_base1 tmp (*this);
            *this = -- this->base ();
            return tmp;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 &operator -- () {
            // Comeau recommends...
            return *this = ++ this->base ();
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 operator -- (int) {
            // Comeau recommends...
            reverse_iterator_base1 tmp (*this);
            *this = ++ this->base ();
            return tmp;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 &operator += (difference_type n) {
            // Comeau recommends...
            return *this = this->base () - n;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base1 &operator -= (difference_type n) {
            // Comeau recommends...
            return *this = this->base () + n;
        }
#endif

#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend reverse_iterator_base1 operator + (const reverse_iterator_base1 &it, difference_type n) {
            reverse_iterator_base1 tmp (it);
            return tmp += n;
        }
        BOOST_UBLAS_INLINE
        friend reverse_iterator_base1 operator + (difference_type n, const reverse_iterator_base1 &it) {
            reverse_iterator_base1 tmp (it);
            return tmp += n;
        }
        BOOST_UBLAS_INLINE
        friend reverse_iterator_base1 operator - (const reverse_iterator_base1 &it, difference_type n) {
            reverse_iterator_base1 tmp (it);
            return tmp -= n;
        }
        BOOST_UBLAS_INLINE
        friend difference_type operator - (const reverse_iterator_base1 &it1, const reverse_iterator_base1 &it2) {
            return it2.base () - it1.base ();
        }
#endif

        BOOST_UBLAS_INLINE
        const container_type &operator () () const {
            // Comeau recommends...
            return this->base () ();
        }

        BOOST_UBLAS_INLINE
        size_type index1 () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).index1 ();
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).index2 ();
        }

        BOOST_UBLAS_INLINE
        dual_iterator_type begin () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).begin ();
        }
        BOOST_UBLAS_INLINE
        dual_iterator_type end () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).end ();
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rbegin () const {
            return dual_reverse_iterator_type (end ());
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rend () const {
            return dual_reverse_iterator_type (begin ());
        }
    };

#ifdef BOOST_UBLAS_NO_MEMBER_FRIENDS
    template<class I>
    BOOST_UBLAS_INLINE
    reverse_iterator_base1<I> operator + (const reverse_iterator_base1<I> &it, std::ptrdiff_t n) {
        reverse_iterator_base1<I> tmp (it);
        return tmp += n;
    }
    template<class I>
    BOOST_UBLAS_INLINE
    reverse_iterator_base1<I> operator + (std::ptrdiff_t n, const reverse_iterator_base1<I> &it) {
        reverse_iterator_base1<I> tmp (it);
        return tmp += n;
    }
    template<class I>
    BOOST_UBLAS_INLINE
    reverse_iterator_base1<I> operator - (const reverse_iterator_base1<I> &it, std::ptrdiff_t n) {
        reverse_iterator_base1<I> tmp (it);
        return tmp -= n;
    }
    template<class I>
    BOOST_UBLAS_INLINE
    std::ptrdiff_t operator - (const reverse_iterator_base1<I> &it1, const reverse_iterator_base1<I> &it2) {
        return it2.base () - it1.base ();
    }
#endif

  /** \brief 2nd base class of all matrix reverse iterators. (non-MSVC version)
   *
   * \param I the derived iterator type
   *
   * The reverse iterator implements a bidirectional iterator
   * reversing the elements of the underlying iterator. It
   * implements most operators of a random access iterator.
   *
   * uBLAS extension: it.index1(), it.index2() and access to
   * the dual iterator via begin(), end(), rbegin(), rend()
   *
   * Note: this type is _identical_ to reverse_iterator_base1
   */

    // Renamed this class from reverse_iterator2 to get
    // typedef reverse_iterator2<...> reverse_iterator2
    // working. Thanks to Gabriel Dos Reis for explaining this.
    template <class I>
    class reverse_iterator_base2:
        public std::reverse_iterator<I> {
    public:
        typedef typename I::container_type container_type;
        typedef typename container_type::size_type size_type;
        typedef typename I::difference_type difference_type;
        typedef I iterator_type;
        typedef typename I::dual_iterator_type dual_iterator_type;
        typedef typename I::dual_reverse_iterator_type dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 ():
            std::reverse_iterator<iterator_type> () {}
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 (const iterator_type &it):
            std::reverse_iterator<iterator_type> (it) {}

#ifndef BOOST_UBLAS_NO_REVERSE_ITERATOR_OVERLOADS
        // Arithmetic
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 &operator ++ () {
            // Comeau recommends...
            return *this = -- this->base ();
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 operator ++ (int) {
            // Comeau recommends...
            reverse_iterator_base2 tmp (*this);
            *this = -- this->base ();
            return tmp;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 &operator -- () {
            // Comeau recommends...
            return *this = ++ this->base ();
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 operator -- (int) {
            // Comeau recommends...
            reverse_iterator_base2 tmp (*this);
            *this = ++ this->base ();
            return tmp;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 &operator += (difference_type n) {
            // Comeau recommends...
            return *this = this->base () - n;
        }
        BOOST_UBLAS_INLINE
        reverse_iterator_base2 &operator -= (difference_type n) {
            // Comeau recommends...
            return *this = this->base () + n;
        }
#endif

#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend reverse_iterator_base2 operator + (const reverse_iterator_base2 &it, difference_type n) {
            reverse_iterator_base2 tmp (it);
            return tmp += n;
        }
        BOOST_UBLAS_INLINE
        friend reverse_iterator_base2 operator + (difference_type n, const reverse_iterator_base2 &it) {
            reverse_iterator_base2 tmp (it);
            return tmp += n;
        }
        BOOST_UBLAS_INLINE
        friend reverse_iterator_base2 operator - (const reverse_iterator_base2 &it, difference_type n) {
            reverse_iterator_base2 tmp (it);
            return tmp -= n;
        }
        BOOST_UBLAS_INLINE
        friend difference_type operator - (const reverse_iterator_base2 &it1, const reverse_iterator_base2 &it2) {
            return it2.base () - it1.base ();
        }
#endif

        BOOST_UBLAS_INLINE
        const container_type &operator () () const {
            // Comeau recommends...
            return this->base () ();
        }

        BOOST_UBLAS_INLINE
        size_type index1 () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).index1 ();
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).index2 ();
        }

        BOOST_UBLAS_INLINE
        dual_iterator_type begin () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).begin ();
        }
        BOOST_UBLAS_INLINE
        dual_iterator_type end () const {
            // Comeau recommends...
            iterator_type tmp (this->base ());
            return (-- tmp).end ();
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rbegin () const {
            return dual_reverse_iterator_type (end ());
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rend () const {
            return dual_reverse_iterator_type (begin ());
        }
    };

#ifdef BOOST_UBLAS_NO_MEMBER_FRIENDS
    template<class I>
    BOOST_UBLAS_INLINE
    reverse_iterator_base2<I> operator + (const reverse_iterator_base2<I> &it, std::ptrdiff_t n) {
        reverse_iterator_base2<I> tmp (it);
        return tmp += n;
    }
    template<class I>
    BOOST_UBLAS_INLINE
    reverse_iterator_base2<I> operator + (std::ptrdiff_t n, const reverse_iterator_base2<I> &it) {
        reverse_iterator_base2<I> tmp (it);
        return tmp += n;
    }
    template<class I>
    BOOST_UBLAS_INLINE
    reverse_iterator_base2<I> operator - (const reverse_iterator_base2<I> &it, std::ptrdiff_t n) {
        reverse_iterator_base2<I> tmp (it);
        return tmp -= n;
    }
    template<class I>
    BOOST_UBLAS_INLINE
    std::ptrdiff_t operator - (const reverse_iterator_base2<I> &it1, const reverse_iterator_base2<I> &it2) {
        return it2.base () - it1.base ();
    }
#endif

#endif

  /** \brief A class implementing an indexed random access iterator.
   *
   * \param C the mutable container type
   * \param IC the iterator category
   *
   * This class implements a random access iterator. The current 
   * position is stored as the unsigned integer it_ and the
   * values are accessed via operator()(it_) of the container.
   *
   * uBLAS extension: index()
   */

    template<class C, class IC>
    class indexed_iterator:
        public container_reference<C>,
        public random_access_iterator_base<IC,
                                           indexed_iterator<C, IC>,
                                           typename C::value_type,
                                           typename C::difference_type> {
    public:
        typedef C container_type;
        typedef IC iterator_category;
        typedef typename container_type::size_type size_type;
        typedef typename container_type::difference_type difference_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::reference reference;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        indexed_iterator ():
            container_reference<container_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        indexed_iterator (container_type &c, size_type it):
            container_reference<container_type> (c), it_ (it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        indexed_iterator &operator ++ () {
            ++ it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_iterator &operator -- () {
            -- it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_iterator &operator += (difference_type n) {
            it_ += n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_iterator &operator -= (difference_type n) {
            it_ -= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        difference_type operator - (const indexed_iterator &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            return it_ - it.it_;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        reference operator * () const {
            BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
            return (*this) () (it_);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (difference_type n) const {
            return *((*this) + n);
        }

        // Index
        BOOST_UBLAS_INLINE
        size_type index () const {
            return it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        indexed_iterator &operator = (const indexed_iterator &it) {
            // FIX: ICC needs full qualification?!
            // assign (&it ());
            container_reference<C>::assign (&it ());
            it_ = it.it_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const indexed_iterator &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            return it_ == it.it_;
        }
        BOOST_UBLAS_INLINE
        bool operator < (const indexed_iterator &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            return it_ < it.it_;
        }

    private:
        size_type it_;
    };

#ifdef BOOST_MSVC_STD_ITERATOR
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator<C, I> operator ++ (const indexed_iterator<C, I> &it, int) {
        indexed_iterator<C, I> tmp (it);
        ++ tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator<C, I> operator -- (const indexed_iterator<C, I> &it, int) {
        indexed_iterator<C, I> tmp (it);
        -- tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator<C, I> operator + (const indexed_iterator<C, I> &it, std::ptrdiff_t n) {
        indexed_iterator<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator<C, I> operator + (std::ptrdiff_t n, const indexed_iterator<C, I> &it) {
        indexed_iterator<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator<C, I> operator - (const indexed_iterator<C, I> &it, std::ptrdiff_t n) {
        indexed_iterator<C, I> tmp (it);
        return tmp -= n;
    }
#endif

  /** \brief A class implementing an indexed random access iterator.
   *
   * \param C the mutable container type
   * \param IC the iterator category
   *
   * This class implements a random access iterator. The current 
   * position is stored as the unsigned integer \c it_ and the
   * values are accessed via \c operator()(it_) of the container.
   *
   * uBLAS extension: \c index()
   *
   * Note: there is an automatic conversion from 
   * \c indexed_iterator to \c indexed_const_iterator
   */

    template<class C, class IC>
    class indexed_const_iterator:
        public container_const_reference<C>,
        public random_access_iterator_base<IC,
                                           indexed_const_iterator<C, IC>,
                                           typename C::value_type,
                                           typename C::difference_type> {
    public:
        typedef C container_type;
        typedef IC iterator_category;
        typedef typename container_type::size_type size_type;
        typedef typename container_type::difference_type difference_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::const_reference reference;
        typedef indexed_iterator<container_type, iterator_category> iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        indexed_const_iterator ():
            container_const_reference<container_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        indexed_const_iterator (const container_type &c, size_type it):
            container_const_reference<container_type> (c), it_ (it) {}
        BOOST_UBLAS_INLINE 
        indexed_const_iterator (const iterator_type &it):
            container_const_reference<container_type> (it ()), it_ (it.index ()) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        indexed_const_iterator &operator ++ () {
            ++ it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_const_iterator &operator -- () {
            -- it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_const_iterator &operator += (difference_type n) {
            it_ += n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_const_iterator &operator -= (difference_type n) {
            it_ -= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        difference_type operator - (const indexed_const_iterator &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            return it_ - it.it_;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        reference operator * () const {
            BOOST_UBLAS_CHECK (index () < (*this) ().size (), bad_index ());
            return (*this) () (it_);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (difference_type n) const {
            return *((*this) + n);
        }

        // Index
        BOOST_UBLAS_INLINE
        size_type index () const {
            return it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        indexed_const_iterator &operator = (const indexed_const_iterator &it) {
            // FIX: ICC needs full qualification?!
            // assign (&it ());
            container_const_reference<C>::assign (&it ());
            it_ = it.it_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const indexed_const_iterator &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            return it_ == it.it_;
        }
        BOOST_UBLAS_INLINE
        bool operator < (const indexed_const_iterator &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            return it_ < it.it_;
        }

    private:
        size_type it_;

        friend class indexed_iterator<container_type, iterator_category>;
    };

#ifdef BOOST_MSVC_STD_ITERATOR
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator<C, I> operator ++ (const indexed_const_iterator<C, I> &it, int) {
        indexed_const_iterator<C, I> tmp (it);
        ++ tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator<C, I> operator -- (const indexed_const_iterator<C, I> &it, int) {
        indexed_const_iterator<C, I> tmp (it);
        -- tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator<C, I> operator + (const indexed_const_iterator<C, I> &it, std::ptrdiff_t n) {
        indexed_const_iterator<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator<C, I> operator + (std::ptrdiff_t n, const indexed_const_iterator<C, I> &it) {
        indexed_const_iterator<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator<C, I> operator - (const indexed_const_iterator<C, I> &it, std::ptrdiff_t n) {
        indexed_const_iterator<C, I> tmp (it);
        return tmp -= n;
    }
#endif

    template<class C, class IC>
    class indexed_iterator2;

  /** \brief A class implementing an indexed random access iterator 
   * of a matrix.
   *
   * \param C the mutable container type
   * \param IC the iterator category
   *
   * This class implements a random access iterator. The current
   * position is stored as two unsigned integers \c it1_ and \c it2_
   * and the values are accessed via \c operator()(it1_, it2_) of the
   * container. The iterator changes the first index.
   *
   * uBLAS extension: \c index1(), \c index2() and access to the
   * dual iterator via \c begin(), \c end(), \c rbegin() and \c rend()
   *
   * Note: The container has to support the \code find2(rank, i, j) \endcode 
   * method
   */

    template<class C, class IC>
    class indexed_iterator1:
        public container_reference<C>, 
        public random_access_iterator_base<IC,
                                           indexed_iterator1<C, IC>, 
                                           typename C::value_type,
                                           typename C::reference> {
    public:
        typedef C container_type;
        typedef IC iterator_category;
        typedef typename container_type::size_type size_type;
        typedef typename container_type::difference_type difference_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::reference reference;
        typedef indexed_iterator2<container_type, iterator_category> dual_iterator_type;
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base2<dual_iterator_type, value_type, reference> dual_reverse_iterator_type;
#else
        typedef reverse_iterator_base2<dual_iterator_type> dual_reverse_iterator_type;
#endif

        // Construction and destruction
        BOOST_UBLAS_INLINE
        indexed_iterator1 ():
            container_reference<container_type> (), it1_ (), it2_ () {}
        BOOST_UBLAS_INLINE 
        indexed_iterator1 (container_type &c, size_type it1, size_type it2):
            container_reference<container_type> (c), it1_ (it1), it2_ (it2) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        indexed_iterator1 &operator ++ () {
            ++ it1_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_iterator1 &operator -- () {
            -- it1_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_iterator1 &operator += (difference_type n) {
            it1_ += n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_iterator1 &operator -= (difference_type n) {
            it1_ -= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        difference_type operator - (const indexed_iterator1 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
            return it1_ - it.it1_;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        reference operator * () const {
            BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
            BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
            return (*this) () (it1_, it2_);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (difference_type n) const {
            return *((*this) + n);
        }

        // Index
        BOOST_UBLAS_INLINE
        size_type index1 () const {
            return it1_;
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
            return it2_;
        }

        BOOST_UBLAS_INLINE
        dual_iterator_type begin () const {
            return (*this) ().find2 (1, index1 (), 0); 
        }
        BOOST_UBLAS_INLINE
        dual_iterator_type end () const {
            return (*this) ().find2 (1, index1 (), (*this) ().size2 ());
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rbegin () const {
            return dual_reverse_iterator_type (end ());
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rend () const {
            return dual_reverse_iterator_type (begin ());
        }

        // Assignment
        BOOST_UBLAS_INLINE
        indexed_iterator1 &operator = (const indexed_iterator1 &it) {
            // FIX: ICC needs full qualification?!
            // assign (&it ());
            container_reference<C>::assign (&it ());
            it1_ = it.it1_;
            it2_ = it.it2_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const indexed_iterator1 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
            return it1_ == it.it1_;
        }
        BOOST_UBLAS_INLINE
        bool operator < (const indexed_iterator1 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
            return it1_ < it.it1_;
        }

    private:
        size_type it1_;
        size_type it2_;
    };

#ifdef BOOST_MSVC_STD_ITERATOR
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator1<C, I> operator ++ (const indexed_iterator1<C, I> &it, int) {
        indexed_iterator1<C, I> tmp (it);
        ++ tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator1<C, I> operator -- (const indexed_iterator1<C, I> &it, int) {
        indexed_iterator1<C, I> tmp (it);
        -- tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator1<C, I> operator + (const indexed_iterator1<C, I> &it, std::ptrdiff_t n) {
        indexed_iterator1<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator1<C, I> operator + (std::ptrdiff_t n, const indexed_iterator1<C, I> &it) {
        indexed_iterator1<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator1<C, I> operator - (const indexed_iterator1<C, I> &it, std::ptrdiff_t n) {
        indexed_iterator1<C, I> tmp (it);
        return tmp -= n;
    }
#endif

    template<class C, class IC>
    class indexed_const_iterator2;

  /** \brief A class implementing an indexed random access iterator 
   * of a matrix.
   *
   * \param C the (immutable) container type
   * \param IC the iterator category
   *
   * This class implements a random access iterator. The current
   * position is stored as two unsigned integers \c it1_ and \c it2_
   * and the values are accessed via \c operator()(it1_, it2_) of the
   * container. The iterator changes the first index.
   *
   * uBLAS extension: \c index1(), \c index2() and access to the
   * dual iterator via \c begin(), \c end(), \c rbegin() and \c rend()
   *
   * Note 1: The container has to support the find2(rank, i, j) method
   *
   * Note 2: there is an automatic conversion from 
   * \c indexed_iterator1 to \c indexed_const_iterator1
   */

    template<class C, class IC>
    class indexed_const_iterator1:
        public container_const_reference<C>, 
        public random_access_iterator_base<IC,
                                           indexed_const_iterator1<C, IC>, 
                                           typename C::value_type,
                                           typename C::const_reference> {
    public:
        typedef C container_type;
        typedef IC iterator_category;
        typedef typename container_type::size_type size_type;
        typedef typename container_type::difference_type difference_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::const_reference reference;
        typedef indexed_iterator1<container_type, iterator_category> iterator_type;
        typedef indexed_const_iterator2<container_type, iterator_category> dual_iterator_type;
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base2<dual_iterator_type, value_type, reference> dual_reverse_iterator_type;
#else
        typedef reverse_iterator_base2<dual_iterator_type> dual_reverse_iterator_type;
#endif

        // Construction and destruction
        BOOST_UBLAS_INLINE
        indexed_const_iterator1 ():
            container_const_reference<container_type> (), it1_ (), it2_ () {}
        BOOST_UBLAS_INLINE
        indexed_const_iterator1 (const container_type &c, size_type it1, size_type it2):
            container_const_reference<container_type> (c), it1_ (it1), it2_ (it2) {}
        BOOST_UBLAS_INLINE 
        indexed_const_iterator1 (const iterator_type &it):
            container_const_reference<container_type> (it ()), it1_ (it.index1 ()), it2_ (it.index2 ()) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        indexed_const_iterator1 &operator ++ () {
            ++ it1_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_const_iterator1 &operator -- () {
            -- it1_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_const_iterator1 &operator += (difference_type n) {
            it1_ += n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_const_iterator1 &operator -= (difference_type n) {
            it1_ -= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        difference_type operator - (const indexed_const_iterator1 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
            return it1_ - it.it1_;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        reference operator * () const {
            BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
            BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
            return (*this) () (it1_, it2_);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (difference_type n) const {
            return *((*this) + n);
        }

        // Index
        BOOST_UBLAS_INLINE
        size_type index1 () const {
            return it1_;
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
            return it2_;
        }

        BOOST_UBLAS_INLINE
        dual_iterator_type begin () const {
            return (*this) ().find2 (1, index1 (), 0); 
        }
        BOOST_UBLAS_INLINE
        dual_iterator_type end () const {
            return (*this) ().find2 (1, index1 (), (*this) ().size2 ()); 
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rbegin () const {
            return dual_reverse_iterator_type (end ()); 
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rend () const {
            return dual_reverse_iterator_type (begin ()); 
        }

        // Assignment
        BOOST_UBLAS_INLINE
        indexed_const_iterator1 &operator = (const indexed_const_iterator1 &it) {
            // FIX: ICC needs full qualification?!
            // assign (&it ());
            container_const_reference<C>::assign (&it ());
            it1_ = it.it1_;
            it2_ = it.it2_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const indexed_const_iterator1 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
            return it1_ == it.it1_;
        }
        BOOST_UBLAS_INLINE
        bool operator < (const indexed_const_iterator1 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
            return it1_ < it.it1_;
        }

    private:
        size_type it1_;
        size_type it2_;

        friend class indexed_iterator1<container_type, iterator_category>;
    };

#ifdef BOOST_MSVC_STD_ITERATOR
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator1<C, I> operator ++ (const indexed_const_iterator1<C, I> &it, int) {
        indexed_const_iterator1<C, I> tmp (it);
        ++ tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator1<C, I> operator -- (const indexed_const_iterator1<C, I> &it, int) {
        indexed_const_iterator1<C, I> tmp (it);
        -- tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator1<C, I> operator + (const indexed_const_iterator1<C, I> &it, std::ptrdiff_t n) {
        indexed_const_iterator1<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator1<C, I> operator + (std::ptrdiff_t n, const indexed_const_iterator1<C, I> &it) {
        indexed_const_iterator1<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator1<C, I> operator - (const indexed_const_iterator1<C, I> &it, std::ptrdiff_t n) {
        indexed_const_iterator1<C, I> tmp (it);
        return tmp -= n;
    }
#endif

  /** \brief A class implementing an indexed random access iterator 
   * of a matrix.
   *
   * \param C the mutable container type
   * \param IC the iterator category
   *
   * This class implements a random access iterator. The current
   * position is stored as two unsigned integers \c it1_ and \c it2_
   * and the values are accessed via \c operator()(it1_, it2_) of the
   * container. The iterator changes the second index.
   *
   * uBLAS extension: \c index1(), \c index2() and access to the
   * dual iterator via \c begin(), \c end(), \c rbegin() and \c rend()
   *
   * Note: The container has to support the find1(rank, i, j) method
   */
    template<class C, class IC>
    class indexed_iterator2:
        public container_reference<C>, 
        public random_access_iterator_base<IC,
                                           indexed_iterator2<C, IC>, 
                                           typename C::value_type,
                                           typename C::reference> {
    public:
        typedef C container_type;
        typedef IC iterator_category;
        typedef typename container_type::size_type size_type;
        typedef typename container_type::difference_type difference_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::reference reference;
        typedef indexed_iterator1<container_type, iterator_category> dual_iterator_type;
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<dual_iterator_type, value_type, reference> dual_reverse_iterator_type;
#else
        typedef reverse_iterator_base1<dual_iterator_type> dual_reverse_iterator_type;
#endif

        // Construction and destruction
        BOOST_UBLAS_INLINE
        indexed_iterator2 ():
            container_reference<container_type> (), it1_ (), it2_ () {}
        BOOST_UBLAS_INLINE
        indexed_iterator2 (container_type &c, size_type it1, size_type it2):
            container_reference<container_type> (c), it1_ (it1), it2_ (it2) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        indexed_iterator2 &operator ++ () {
            ++ it2_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_iterator2 &operator -- () {
            -- it2_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_iterator2 &operator += (difference_type n) {
            it2_ += n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_iterator2 &operator -= (difference_type n) {
            it2_ -= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        difference_type operator - (const indexed_iterator2 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
            return it2_ - it.it2_;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        reference operator * () const {
            BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
            BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
            return (*this) () (it1_, it2_);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (difference_type n) const {
            return *((*this) + n);
        }

        // Index
        BOOST_UBLAS_INLINE
        size_type index1 () const {
            return it1_;
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
            return it2_;
        }

        BOOST_UBLAS_INLINE
        dual_iterator_type begin () const {
            return (*this) ().find1 (1, 0, index2 ());
        }
        BOOST_UBLAS_INLINE
        dual_iterator_type end () const {
            return (*this) ().find1 (1, (*this) ().size1 (), index2 ());
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rbegin () const {
            return dual_reverse_iterator_type (end ());
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rend () const {
            return dual_reverse_iterator_type (begin ());
        }

        // Assignment
        BOOST_UBLAS_INLINE
        indexed_iterator2 &operator = (const indexed_iterator2 &it) {
            // FIX: ICC needs full qualification?!
            // assign (&it ());
            container_reference<C>::assign (&it ());
            it1_ = it.it1_;
            it2_ = it.it2_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const indexed_iterator2 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
            return it2_ == it.it2_;
        }
        BOOST_UBLAS_INLINE
        bool operator < (const indexed_iterator2 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
            return it2_ < it.it2_;
        }

    private:
        size_type it1_;
        size_type it2_;
    };

#ifdef BOOST_MSVC_STD_ITERATOR
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator2<C, I> operator ++ (const indexed_iterator2<C, I> &it, int) {
        indexed_iterator2<C, I> tmp (it);
        ++ tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator2<C, I> operator -- (const indexed_iterator2<C, I> &it, int) {
        indexed_iterator2<C, I> tmp (it);
        -- tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator2<C, I> operator + (const indexed_iterator2<C, I> &it, std::ptrdiff_t n) {
        indexed_iterator2<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator2<C, I> operator + (std::ptrdiff_t n, const indexed_iterator2<C, I> &it) {
        indexed_iterator2<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_iterator2<C, I> operator - (const indexed_iterator2<C, I> &it, std::ptrdiff_t n) {
        indexed_iterator2<C, I> tmp (it);
        return tmp -= n;
    }
#endif

  /** \brief A class implementing an indexed random access iterator 
   * of a matrix.
   *
   * \param C the (immutable) container type
   * \param IC the iterator category
   *
   * This class implements a random access iterator. The current
   * position is stored as two unsigned integers \c it1_ and \c it2_
   * and the values are accessed via \c operator()(it1_, it2_) of the
   * container. The iterator changes the second index.
   *
   * uBLAS extension: \c index1(), \c index2() and access to the
   * dual iterator via \c begin(), \c end(), \c rbegin() and \c rend()
   *
   * Note 1: The container has to support the \c find2(rank, i, j) method
   *
   * Note 2: there is an automatic conversion from 
   * \c indexed_iterator2 to \c indexed_const_iterator2
   */

    template<class C, class IC>
    class indexed_const_iterator2:
        public container_const_reference<C>,
        public random_access_iterator_base<IC,
                                           indexed_const_iterator2<C, IC>,
                                           typename C::value_type,
                                           typename C::const_reference> {
    public:
        typedef C container_type;
        typedef IC iterator_category;
        typedef typename container_type::size_type size_type;
        typedef typename container_type::difference_type difference_type;
        typedef typename container_type::value_type value_type;
        typedef typename container_type::const_reference reference;
        typedef indexed_iterator2<container_type, iterator_category> iterator_type;
        typedef indexed_const_iterator1<container_type, iterator_category> dual_iterator_type;
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef reverse_iterator_base1<dual_iterator_type, value_type, reference> dual_reverse_iterator_type;
#else
        typedef reverse_iterator_base1<dual_iterator_type> dual_reverse_iterator_type;
#endif

        // Construction and destruction
        BOOST_UBLAS_INLINE
        indexed_const_iterator2 ():
            container_const_reference<container_type> (), it1_ (), it2_ () {}
        BOOST_UBLAS_INLINE
        indexed_const_iterator2 (const container_type &c, size_type it1, size_type it2):
            container_const_reference<container_type> (c), it1_ (it1), it2_ (it2) {}
        BOOST_UBLAS_INLINE
        indexed_const_iterator2 (const iterator_type &it):
            container_const_reference<container_type> (it ()), it1_ (it.index1 ()), it2_ (it.index2 ()) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        indexed_const_iterator2 &operator ++ () {
            ++ it2_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_const_iterator2 &operator -- () {
            -- it2_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_const_iterator2 &operator += (difference_type n) {
            it2_ += n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        indexed_const_iterator2 &operator -= (difference_type n) {
            it2_ -= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        difference_type operator - (const indexed_const_iterator2 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
            return it2_ - it.it2_;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        reference operator * () const {
            BOOST_UBLAS_CHECK (index1 () < (*this) ().size1 (), bad_index ());
            BOOST_UBLAS_CHECK (index2 () < (*this) ().size2 (), bad_index ());
            return (*this) () (it1_, it2_);
        }
        BOOST_UBLAS_INLINE
        reference operator [] (difference_type n) const {
            return *((*this) + n);
        }

        // Index
        BOOST_UBLAS_INLINE
        size_type index1 () const {
            return it1_;
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
            return it2_;
        }

        BOOST_UBLAS_INLINE
        dual_iterator_type begin () const {
            return (*this) ().find1 (1, 0, index2 ());
        }
        BOOST_UBLAS_INLINE
        dual_iterator_type end () const {
            return (*this) ().find1 (1, (*this) ().size1 (), index2 ());
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rbegin () const {
            return dual_reverse_iterator_type (end ());
        }
        BOOST_UBLAS_INLINE
        dual_reverse_iterator_type rend () const {
            return dual_reverse_iterator_type (begin ());
        }

        // Assignment
        BOOST_UBLAS_INLINE
        indexed_const_iterator2 &operator = (const indexed_const_iterator2 &it) {
            // FIX: ICC needs full qualification?!
            // assign (&it ());
            container_const_reference<C>::assign (&it ());
            it1_ = it.it1_;
            it2_ = it.it2_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const indexed_const_iterator2 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
            return it2_ == it.it2_;
        }
        BOOST_UBLAS_INLINE
        bool operator < (const indexed_const_iterator2 &it) const {
            BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
            BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
            return it2_ < it.it2_;
        }

    private:
        size_type it1_;
        size_type it2_;

        friend class indexed_iterator2<container_type, iterator_category>;
    };

#ifdef BOOST_MSVC_STD_ITERATOR
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator2<C, I> operator ++ (const indexed_const_iterator2<C, I> &it, int) {
        indexed_const_iterator2<C, I> tmp (it);
        ++ tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator2<C, I> operator -- (const indexed_const_iterator2<C, I> &it, int) {
        indexed_const_iterator2<C, I> tmp (it);
        -- tmp;
        return tmp;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator2<C, I> operator + (const indexed_const_iterator2<C, I> &it, std::ptrdiff_t n) {
        indexed_const_iterator2<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator2<C, I> operator + (std::ptrdiff_t n, const indexed_const_iterator2<C, I> &it) {
        indexed_const_iterator2<C, I> tmp (it);
        return tmp += n;
    }
    template<class C, class I>
    BOOST_UBLAS_INLINE
    indexed_const_iterator2<C, I> operator - (const indexed_const_iterator2<C, I> &it, std::ptrdiff_t n) {
        indexed_const_iterator2<C, I> tmp (it);
        return tmp -= n;
    }
#endif

}}}

#endif
