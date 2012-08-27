// Copyright (c) 2001-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>

#ifndef CGAL_FLATTENING_ITERATOR_H
#define CGAL_FLATTENING_ITERATOR_H 1

#include <CGAL/basic.h>
#include <CGAL/Nested_iterator.h>

#include <iterator>
#include <functional>

// LiS2CGAL check whether Nested_iterator in CGAL works in STL_extensions

namespace CGAL {

/*! Suppose you have an iterator range and suppose further that
 *  the values in that iterator range are containers and thus
 *  define iterator ranges themselves (accessible through
 *  some \c begin and \c end functions).
 *  \e Flattening the iterator range means to turn it into one
 *  iterator range whose values are the values of the containers
 *  in the original range.
 *  For example, a range of lists (1,2),(3),(4,5,6) can be
 *  flattened into a range 1,2,3,4,5,6
 *
 *  The class templates \c Flattening_iterator and
 *  \c Flattening_const_iterator implement this flattening.
 *  So far, only the \c const version has been implemented,
 *  i.e. the values in the flattened range are read-only.
 *
 *  But what about an iterator range whose values are lists of
 *  vectors of vectors of ... and so on?
 *  Yes, it is possible to flatten such nested containers, too,
 *  by recursive application of \c Flattening*_iterator.
 *  The class templates \c Recursive_flattening and
 *  \c Recursive_const_flattening offer the necessary
 *  typedefs and conversions.
 *  So far, only the \c const version has been implemented.
 */

template <int level_, class InputIterator>
class Recursive_const_flattening;

template <class InputIterator>
class Recursive_const_flattening<0, InputIterator> {
public:
    typedef Recursive_const_flattening Self;
    static const int level = 0;
    typedef InputIterator  Input_iterator;

    typedef Input_iterator Recursive_flattening_iterator;

    struct Flatten {
        typedef Recursive_flattening_iterator result_type;
        typedef Input_iterator argument_type;

        Recursive_flattening_iterator
        operator () (Input_iterator /*end*/, Input_iterator it) {
            return it;
        }
    };
};

template <class InputIterator>
class Recursive_const_flattening<1, InputIterator> {
public:
    typedef Recursive_const_flattening Self;
    static const int level = 1;
    typedef InputIterator Input_iterator;

private:
    struct Nested_iterator_traits
    {
      typedef Input_iterator Base_iterator;
      typedef typename std::iterator_traits<Input_iterator>::value_type::const_iterator
              Iterator;

      Iterator begin(Input_iterator it) const { return it->begin(); }
      Iterator end  (Input_iterator it) const { return it->end();   }
    };

public:

    typedef CGAL::Nested_iterator< Input_iterator, Nested_iterator_traits >
            Recursive_flattening_iterator;

    struct Flatten {
        typedef Recursive_flattening_iterator result_type;
        typedef Input_iterator                argument_type;

        Recursive_flattening_iterator
        operator () (Input_iterator end, Input_iterator it) {
            return Recursive_flattening_iterator(end,it);
        }
    };
};

/*! \ingroup LiS_Flattening_iterator
    \brief Recursive application of \c Flattening_const_iterator

    An instance \c Recursive_const_flattening<level,InputIterator>
    of this class template contains a typedef \c Recursive_flattening_iterator.
    This is a \c level -fold nested instance of
    \c Flattening_const_iterator, where \c level can be any
    non-negative integer. At each nesting level, begin and end
    functors are supplied to \c Flattening_const_iterator
    which assume that their arguments have a
    \c const_iterator typedef and \c begin() and \c end()
    member functions in the style of the STL.

    The functor <tt>Flatten()(it)</tt> converts an \c InputIterator \c it
    to a \c Recursive_flattening_iterator. Converting each endpoint
    of an iterator range [first,beyond) yields an iterator range
    in which \c level levels of packing into containers have been
    unpacked. (In particular, the case \c level==0 is permissible
    and a no-op.)

    To get a \c const_iterator out of a container \c c of type \c C
    which is not \c const itself, one can use syntax like
    <tt>const_cast<const C&>(c).begin())</tt>. This may sometimes
    be necessary in conjunction with \c Recursive_const_flattening.
 */
template< int level_, class InputIterator >
class Recursive_const_flattening {
public:
    //! this instance itself
    typedef Recursive_const_flattening Self;
    //! this instance's first template argument
    static const int level = level_;
    //! this instance's second template argument
    typedef InputIterator Input_iterator;

    typedef Recursive_const_flattening<
            level-1,
            typename std::iterator_traits<Input_iterator>::value_type::const_iterator
        > Nested_self;

private:
    struct Nested_iterator_traits
    {
      typedef Input_iterator
              Base_iterator;
      typedef typename Nested_self::Recursive_flattening_iterator
              Iterator;

      Iterator begin(Input_iterator it) const { return typename Nested_self::Flatten()(it->end(),it->begin()); }
      Iterator end  (Input_iterator it) const { return typename Nested_self::Flatten()(it->end(),it->end());   }
    };

public:
    typedef CGAL::Nested_iterator< Input_iterator, Nested_iterator_traits >
            Recursive_flattening_iterator;

    //! conversion functor (model of STL concept \c AdaptableUnaryFunction )
    struct Flatten {
        //! result type
        typedef Recursive_flattening_iterator result_type;
        //! argument type
        typedef Input_iterator                argument_type;

        //! conversion functor call
        Recursive_flattening_iterator
        operator () (Input_iterator end,Input_iterator it) {
            return Recursive_flattening_iterator(end,it);
        }
    };
};

#ifdef DOXYGEN_RUNNING
/*! \ingroup LiS_Flattening_iterator
    \brief (unimplemented)

    This class template is unimplemented.
    Only \c Recursive_const_flattening is available at this point.
 */
template <int level_, class InputIterator>
class Recursive_flattening { };
#endif // DOXYGEN_RUNNING

/*
 * Part 2: Helper functions
 */

/*! \relates Recursive_const_flattening
 *  \brief map \c it to
 *  <tt>Recursive_const_flattening<level, InputIterator>::Flatten()(it)</tt>
 *
 *  See \c Recursive_const_flattening for explanations.
 *  This function just exists to save typing via overloading resolution.
 *  You have to specify \c level explicitly in any case.
 *  For \c recursive_const_flattener<1,...>() , there is the
 *  shorthand \c const_flattener<...>() .
 */
template <int level, class InputIterator> inline
typename Recursive_const_flattening<level, InputIterator>
    ::Recursive_flattening_iterator
recursive_const_flattener(InputIterator end,InputIterator it) {
    return typename Recursive_const_flattening<level, InputIterator>::Flatten()(end,it);
}

/*! \relates Flattening_const_iterator
 *  \brief map \c it to
 *  <tt>Recursive_const_flattening<1, InputIterator>::Flatten()(it)</tt>
 *
 *  A function call <tt>fi = const_flattener(it)</tt> converts the iterator
 *   \c it to an instance of \c Flattening_const_iterator (see ibid.).
 *  The template arguments are chosen as follows:
 *    - \c InputIterator1 is \c InputIterator
 *    - \c InputIterator2 is \c InputIterator::value_type::const_iterator
 *    - \c UnaryFunction1/2 are set to functors that call \c begin()
 *       and \c end() member functions.
 *
 *  This function helps to save typing and to avoid manual specification
 *  of \c UnaryFunction1/2.
 */
template <class InputIterator> inline
typename Recursive_const_flattening<1, InputIterator>
    ::Recursive_flattening_iterator
const_flattener(InputIterator end,InputIterator it) {
    return typename Recursive_const_flattening<1, InputIterator>::Flatten()(end,it);
}

#ifdef DOXYGEN_RUNNING
/*! \relates Recursive_flattener
    \brief (unimplemented)

    This function template is unimplemented.
    Only \c recursive_const_flattener is available at this point.
 */
template <int level, class InputIterator> inline
typename Recursive_flattening<level, InputIterator>
    ::Recursive_flattening_iterator
recursive_flattener(InputIterator it);

/*! \relates Flattening_iteratos
    \brief (unimplemented)

    This function template is unimplemented.
    Only \c const_flattener is available at this point.
 */
template <class InputIterator> inline
typename Recursive_flattening<1, InputIterator>
    ::Recursive_flattening_iterator
flattener(InputIterator it);
#endif // DOXYGEN_RUNNING


} //namespace CGAL

#endif // LiS_FLATTENING_ITERATOR_H
