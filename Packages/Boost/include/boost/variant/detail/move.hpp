//-----------------------------------------------------------------------------
// boost variant/detail/move.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002-2003
// Eric Friedman
//
// See below original copyright by Andrei Alexandrescu.
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_VARIANT_DETAIL_MOVE_HPP
#define BOOST_VARIANT_DETAIL_MOVE_HPP

#include <iterator> // for iterator_traits
#include <new> // for placement new

#include "boost/config.hpp"
#include "boost/detail/workaround.hpp"
#include "boost/mpl/if.hpp"
#include "boost/type_traits/is_base_and_derived.hpp"

namespace boost {
namespace detail { namespace variant {

//////////////////////////////////////////////////////////////////////////
// forward declares
//
// NOTE: Incomplete until (if?) Boost.Move becomes part of Boost.
//
template <typename Deriving> class moveable;
template <typename T>        class move_source;
template <typename T>        class move_return;

namespace detail {

// (detail) moveable_tag
//
// Concrete type from which moveable<T> derives.
//
// TODO: Move into moveable_fwd.hpp and define has_move_constructor.
//
template <typename Deriving>
struct moveable_tag
{
};

} // namespace detail

//////////////////////////////////////////////////////////////////////////
// function template move
//
// Takes a T& and returns, if T derives moveable<T>, a move_source<T> for
// the object; else, returns the T&.
//

namespace detail {

// (detail) class template move_type
//
// Metafunction that, given moveable T, provides move_source<T>, else T&.
//
template <typename T>
struct move_type
{
public: // metafunction result

    typedef typename mpl::if_<
          is_base_and_derived<detail::moveable_tag<T>, T>
        , move_source<T>
        , T&
        >::type type;

};

} // namespace detail

template <typename T>
inline
    typename detail::move_type<T>::type
move(T& source)
{
    typedef typename detail::move_type<T>::type
        move_t;

    return move_t(source);
}

//////////////////////////////////////////////////////////////////////////
// class template return_t
//
// Metafunction that, given moveable T, provides move_return<T>, else T.
//
template <typename T>
struct return_t
{
public: // metafunction result

    typedef typename mpl::if_<
          is_base_and_derived<moveable<T>, T>
        , move_return<T>
        , T
        >::type type;

};

//////////////////////////////////////////////////////////////////////////
// function template move_swap
//
// Swaps using Koenig lookup but falls back to move-swap for primitive
// types and on non-conforming compilers.
//

#if   defined(BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP)   \
 ||   BOOST_WORKAROUND(__GNUC__, BOOST_TESTED_AT(2))

// [Indicate that move_swap by overload is disabled...]
#define BOOST_NO_MOVE_SWAP_BY_OVERLOAD

// [...and provide straight swap-by-move implementation:]
template <typename T>
inline void move_swap(T& lhs, T& rhs)
{
    T tmp( boost::detail::variant::move(lhs) );
    lhs = boost::detail::variant::move(rhs);
    rhs = boost::detail::variant::move(tmp);
}

#else// !workaround

namespace detail { namespace move_swap {

template <typename T>
inline void swap(T& lhs, T& rhs)
{
    T tmp( boost::detail::variant::move(lhs) );
    lhs = boost::detail::variant::move(rhs);
    rhs = boost::detail::variant::move(tmp);
}

}} // namespace detail::move_swap

template <typename T>
inline void move_swap(T& lhs, T& rhs)
{
    using detail::move_swap::swap;

    swap(lhs, rhs);
}

#endif // workaround

}} // namespace detail::variant
} // namespace boost

#endif // BOOST_VARIANT_DETAIL_MOVE_HPP


/* This file derivative of MoJO. Much thanks to Andrei for his initial work.
 * See <http://www.cuj.com/experts/2102/alexandr.htm> for information on MOJO.

 * Original copyright -- on mojo.h -- follows:

////////////////////////////////////////////////////////////////////////////////
// MOJO: MOving Joint Objects
// Copyright (c) 2002 by Andrei Alexandrescu
//
// Created by Andrei Alexandrescu
//
// Permission to use, copy, modify, distribute and sell this software for any 
//     purpose is hereby granted without fee, provided that the above copyright 
//     notice appear in all copies and that both that copyright notice and this 
//     permission notice appear in supporting documentation.
// The author makes no representations about the suitability of this software 
//     for any purpose. It is provided "as is" 
//     without express or implied warranty.
////////////////////////////////////////////////////////////////////////////////

*/
