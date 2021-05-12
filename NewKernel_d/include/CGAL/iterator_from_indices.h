// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_ITERATOR_FROM_INDICES_H
#define CGAL_ITERATOR_FROM_INDICES_H
#include <CGAL/config.h>
#include <boost/iterator/iterator_facade.hpp>
namespace CGAL {
template <class Ref_>
struct Default_coordinate_access {
        typedef Ref_ result_type;
        template<class T> Ref_ operator()(T const& t, std::ptrdiff_t i)const{
                return t[i];
        }
};

//TODO: default type for Value_: typename same_cv<Container_,typename remove_cv<Container_>::type::value_type>::type
template <class Container_, class Value_, class Ref_=
        decltype(std::declval<Container_>()[0])
        , class Coord_access = Default_coordinate_access<Ref_>
        >
class Iterator_from_indices
: public boost::iterator_facade<Iterator_from_indices<Container_,Value_,Ref_,Coord_access>,
        Value_, std::bidirectional_iterator_tag, Ref_>
{
        friend class boost::iterator_core_access;
        //FIXME: use int to save space
        //TODO: use a tuple to save space when Coord_access is empty
        typedef std::ptrdiff_t index_t;
        Container_* cont;
        index_t index;
        CGAL_NO_UNIQUE_ADDRESS Coord_access ca;
        void increment(){ ++index; }
        void decrement(){ --index; }
        void advance(std::ptrdiff_t n){ index+=n; }
        ptrdiff_t distance_to(Iterator_from_indices const& other)const{
                return other.index-index;
        }
        bool equal(Iterator_from_indices const& other)const{
                return index==other.index;
        }
        Ref_ dereference()const{
                //FIXME: use the functor properly
                //Uh, and what did I mean by that?
                return ca(*cont,index);
        }
        public:
        Iterator_from_indices(){}
        Iterator_from_indices(Container_& cont_,std::size_t n)
                : cont(&cont_), index(n) {}
        template<class T>
        Iterator_from_indices(Container_& cont_,std::size_t n,T const&t)
                : cont(&cont_), index(n), ca(t) {}
};
}
#endif // CGAL_ITERATOR_FROM_INDICES_H
