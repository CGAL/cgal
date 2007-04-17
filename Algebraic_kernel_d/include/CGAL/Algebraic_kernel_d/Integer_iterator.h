// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file LiS/Integer_iterator.h
    \brief declares class \c LiS::Integer_iterator
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_INTEGER_ITERATOR_H
#define CGAL_ALGEBRAIC_KERNEL_D_INTEGER_ITERATOR_H 1

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*! \ingroup LiS_iterator
    \brief Iterator for the ordered sequence of all integers

Let \c T be a number type representing the integers (or at least a
subset of them). An object \c ii of type \c Integer_iterator<T>
stores an integer \c i of type T and represents its position in the
ordered sequence of all integers. This means that dereferencing
\c ii yields the value \c i , and moving \c ii forwards or backwards
by some offset amounts to addition or subtraction of that offset
from \c i .

This awkward way of encapsulating the integers is useful as a basis
for further iterator plumbing, like constructing an iterator range
representing successive function values by combining a
\c LiS::Integer_iterator with a \c LiS::Mapping_iterator.

Dereferencing an \c Integer_iterator yields a value, not a
reference. Hence it is not possible to write through an \c
Integer_iterator. This also implies that \c Integer_iterator
is just a model of an InputIterator and not of a more advanced
iterator category, even though \c Integer_iterator supports all
standard iterator operations such as \c += or \c [] .

Each iterator operation is only available if the underlying
integer type \c T supports the corresponding operation.
 */
template <class T>
class Integer_iterator {
    T count;
public:
    //! this class itself
    typedef Integer_iterator Self;
    //! iterator category
    typedef std::input_iterator_tag iterator_category;
    //! value type \c T
    typedef T value_type;
    //! type of a pointer to \c T
    typedef value_type* pointer;
    //! type of a reference to \c T
    typedef value_type& reference;
    //! type for offsets: \c T again
    typedef T difference_type;

    //! default constructor (points to the number \c T(0))
    Integer_iterator() : count(0) {}
    //! constructor with specific number \c t
    explicit Integer_iterator(T t) : count(t) {}

    //! return number
    value_type operator * () const { return count; }
    //! return number plus offset \c n
    value_type operator [] (difference_type n) const { return count + n; }
    //! subtraction of iterators, yielding an offset
    difference_type operator - (Integer_iterator it) const {
        return count - it.count;
    }

    //! equality of numbers
    bool operator == (Integer_iterator it) const { return count == it.count; }
    //! inequality of numbers
    bool operator != (Integer_iterator it) const { return count != it.count; }
    //! comparison of numbers
    bool operator <  (Integer_iterator it) const { return count <  it.count; }
    //! comparison of numbers
    bool operator <= (Integer_iterator it) const { return count <= it.count; }
    //! comparison of numbers
    bool operator >= (Integer_iterator it) const { return count >= it.count; }
    //! comparison of numbers
    bool operator >  (Integer_iterator it) const { return count >  it.count; }

    //! preincrement
    Integer_iterator operator ++ () { ++count; return *this; }
    //! postincrement
    Integer_iterator operator ++ (int) {
        Self tmp = *this;
        count++;
        return tmp;
    }
    //! compound addition
    Integer_iterator operator += (difference_type n) {
        count += n; return *this;
    }
    //! addition
    Integer_iterator operator + (difference_type n) const {
        return Integer_iterator(count + n);
    }

    //! predecrement
    Integer_iterator operator -- () { --count; return *this; }
    //! postdecrement
    Integer_iterator operator -- (int) {
        Self tmp = *this;
        count--;
        return tmp;
    }
    //! compound subtraction
    Integer_iterator operator -= (difference_type n) {
        count -= n; return *this;
    }
    //! subtraction iterator minus offet, yielding iterator
    Integer_iterator operator - (difference_type n) const {
        return Integer_iterator(count - n);
    }
};

/*! \relates LiS::Integer_iterator
 *  \brief addition
 *
 */
template <class T> inline
Integer_iterator<T> operator + (
    typename Integer_iterator<T>::difference_type n,
    Integer_iterator<T> it
) {
        return it + n;
}

/*! \relates LiS::Integer_iterator
 *  \brief returns \c Integer_iterator constructed from \c t
 *
 */
template <class T> inline
Integer_iterator<T> integer_iterator(T t) { return Integer_iterator<T>(t); }

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_INTEGER_ITERATOR_H
// EOF
