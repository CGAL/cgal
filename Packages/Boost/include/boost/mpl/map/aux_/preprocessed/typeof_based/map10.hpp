
// Copyright Aleksey Gurtovoy 2000-2004
// Copyright David Abrahams 2003-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//

// Preprocessed version of "boost/mpl/map/map10.hpp" header
// -- DO NOT modify by hand!

namespace boost { namespace mpl {

template<
      typename P0
    >
struct map1
    : m_item<
          typename P0::first
        , typename P0::second
        , map0<  >
        >
{
};

template<
      typename P0, typename P1
    >
struct map2
    : m_item<
          typename P1::first
        , typename P1::second
        , map1<P0>
        >
{
};

template<
      typename P0, typename P1, typename P2
    >
struct map3
    : m_item<
          typename P2::first
        , typename P2::second
        , map2< P0,P1 >
        >
{
};

template<
      typename P0, typename P1, typename P2, typename P3
    >
struct map4
    : m_item<
          typename P3::first
        , typename P3::second
        , map3< P0,P1,P2 >
        >
{
};

template<
      typename P0, typename P1, typename P2, typename P3, typename P4
    >
struct map5
    : m_item<
          typename P4::first
        , typename P4::second
        , map4< P0,P1,P2,P3 >
        >
{
};

template<
      typename P0, typename P1, typename P2, typename P3, typename P4
    , typename P5
    >
struct map6
    : m_item<
          typename P5::first
        , typename P5::second
        , map5< P0,P1,P2,P3,P4 >
        >
{
};

template<
      typename P0, typename P1, typename P2, typename P3, typename P4
    , typename P5, typename P6
    >
struct map7
    : m_item<
          typename P6::first
        , typename P6::second
        , map6< P0,P1,P2,P3,P4,P5 >
        >
{
};

template<
      typename P0, typename P1, typename P2, typename P3, typename P4
    , typename P5, typename P6, typename P7
    >
struct map8
    : m_item<
          typename P7::first
        , typename P7::second
        , map7< P0,P1,P2,P3,P4,P5,P6 >
        >
{
};

template<
      typename P0, typename P1, typename P2, typename P3, typename P4
    , typename P5, typename P6, typename P7, typename P8
    >
struct map9
    : m_item<
          typename P8::first
        , typename P8::second
        , map8< P0,P1,P2,P3,P4,P5,P6,P7 >
        >
{
};

template<
      typename P0, typename P1, typename P2, typename P3, typename P4
    , typename P5, typename P6, typename P7, typename P8, typename P9
    >
struct map10
    : m_item<
          typename P9::first
        , typename P9::second
        , map9< P0,P1,P2,P3,P4,P5,P6,P7,P8 >
        >
{
};

}}
