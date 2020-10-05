// Copyright (c) 2001-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//

#ifndef CGAL_CACHE_H
#define CGAL_CACHE_H 1

#include <CGAL/basic.h>
#include <CGAL/function_objects.h>

#include <map>

namespace CGAL {

/*! \brief The Cache serves as a constructor for an object of type Output from
 * a object of type Input.
 *
 * The Cache internally uses the class std::map<Input_, Output_,Compare_>.
 * The Cache works as follows. First the \c input is canonicalized by the
 * Canonicalizer and the  result is used as the key of the cache. If the
 * respective object is not in the cache yet, it is constructed by the Creator
 * functor from the canonicalized \c input and added to the map.\n
 *
 * The Compare functor serves as the key comparison function, which must
 * define a total ordering on the Input. \n
 *
 * The default for the Creator_ functor is Creator_1<Input_,Output_>.\n
 * The default for the Canonicalizer_ functor is Creator_1<Input_,
 * Input_>  \n
 * The default for the Compare_ functor is std::less<Input_>. \n
 *
 */
template < class Input_,
           class Output_,
           class Creator_ = Creator_1<Input_, Output_>,
           class Canonicalizer_ = Creator_1<Input_, Input_>,
           class Compare_ = std::less<Input_> >
class Cache{
public:
    //! The Input type of the Cache
    typedef Input_ Input;

    //! The Output type of the Cache
    typedef Output_ Output;

    /*! \brief The Creator functor of the Cache
     *  This functor is used to create the output from the input.
     *  If Output has a constructor from input you can use the default.*/
    typedef Creator_ Creator;

    /*! \brief The Canonicalizer functor of the Cache
     *  This functor is used to canonicalize the input in advance. \n
     *  If the Input already is canonicalized you can use the default.\n */
    typedef Canonicalizer_ Canonicalizer;

    /*! \brief The Compare functor used by the Cache */
    typedef Compare_ Compare;

    typedef Cache<Input,Output,Creator,Canonicalizer,Compare> Self;
private:
    typedef std::map<Input,Output,Compare> Map;
    Map map;
public:
    typedef Map _Rep_type;
    //! Iterator type
    typedef typename _Rep_type::iterator Iterator;
    //! Const_iterator type
    typedef typename _Rep_type::const_iterator Const_iterator;
    //! Reverse_iterator type
    typedef typename _Rep_type::reverse_iterator Reverse_iterator;
    //! Const_reverse_iterator type
    typedef typename _Rep_type::const_reverse_iterator Const_reverse_iterator;
    //! Size_type type
    typedef typename _Rep_type::size_type Size_type;
public:
    //! default constructor with empty table
    Cache() : map() {};

    /*! \brief Returns the respective object of type Output.
     *
     *  If the object is not in the cache, it is constructed form \c key and
     *  added to the cache.
     */
    Output operator () (const Input& input) {
        Canonicalizer canonicalize;
        Input key = canonicalize(input);
        typename Map::iterator it = map.find(key);
        if (it == map.end()) {
            Creator create;
            Output out = create(key);
            map.insert(it,typename Map::value_type(key,out));
            return out;
        } else {
            return (*it).second;
        }
    }

    //! Clears the cache contents
    void clear() { map.clear(); }

    //! Returns whether the cache is empty.
    bool is_empty() const { return map.empty(); }

    //! Returns the current number of different objects in the cache
    Size_type size() { return map.size(); }

    //! Returns the largest possible size of the cache.
    Size_type max_size() const { return map.max_size(); }

    //! Returns an Iterator pointing to the beginning of the cache.
    Iterator begin() { return map.begin(); }

    //! Returns an Iterator pointing to the end of the cache.
    Iterator end() { return map.end(); }

    //!        Returns a Const_iterator pointing to the beginning of the cache.
    Const_iterator begin() const { return map.begin(); }

    //! Returns a Const_iterator pointing to the end of the cache.
    Const_iterator end() const { return map.end(); }

    /*!        \brief Returns a Reverse_iterator pointing to the beginning of the
     *  reversed cache.
     */
    Reverse_iterator rbegin() { return map.rbegin(); }

    /*!        \brief Returns a Reverse_iterator pointing to the end of the reversed
     *  cache.
     */
    Reverse_iterator rend() { return map.rend(); }

    /*!        \brief Returns a Const_reverse_iterator pointing to the beginning of
     *  the reversed cache.
     */
    Const_reverse_iterator rbegin() const { return map.rbegin(); }

    /*!        \brief Returns a Const_reverse_iterator pointing to the end of the
     *  reversed cache.
     */
    Const_reverse_iterator rend() const { return map.rend(); }

};

} //namespace CGAL

#endif
