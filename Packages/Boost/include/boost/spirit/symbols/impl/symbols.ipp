/*=============================================================================
    Copyright (c) 2001-2003 Joel de Guzman
    http://spirit.sourceforge.net/

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#ifndef BOOST_SPIRIT_SYMBOLS_IPP
#define BOOST_SPIRIT_SYMBOLS_IPP

///////////////////////////////////////////////////////////////////////////////
#include <boost/spirit/symbols/impl/tst.ipp>

// MSVC: void warning about the use of 'this' pointer in constructors
#if defined(BOOST_MSVC)
#pragma warning(disable : 4355)
#endif

///////////////////////////////////////////////////////////////////////////////
namespace boost { namespace spirit {

///////////////////////////////////////////////////////////////////////////////
//
//  symbols class implementation
//
///////////////////////////////////////////////////////////////////////////////
template <typename T, typename CharT, typename SetT>
inline symbols<T, CharT, SetT>::symbols()
: SetT()
, add(*this)
{
}

//////////////////////////////////
template <typename T, typename CharT, typename SetT>
symbols<T, CharT, SetT>::symbols(symbols const& other)
: SetT(other)
, parser<symbols<T, CharT, SetT> >()
, add(*this)
{
}

//////////////////////////////////
template <typename T, typename CharT, typename SetT>
inline symbols<T, CharT, SetT>::~symbols()
{}

//////////////////////////////////
template <typename T, typename CharT, typename SetT>
inline symbols<T, CharT, SetT>&
symbols<T, CharT, SetT>::operator=(symbols const& other)
{
    SetT::operator=(other);
    return *this;
}

//////////////////////////////////
template <typename T, typename CharT, typename SetT>
inline symbol_inserter<T, SetT> const&
symbols<T, CharT, SetT>::operator=(CharT const* str)
{
    return add, str;
}

///////////////////////////////////////////////////////////////////////////////
//
//  Symbol table utilities
//
///////////////////////////////////////////////////////////////////////////////
template <typename T, typename CharT, typename SetT>
inline T*
find(symbols<T, CharT, SetT> const& table, CharT const* sym)
{
    CharT const* last = sym;
    while (*last)
        last++;
    scanner<CharT const *> scan(sym, last);
    T* result = table.find(scan);
    return scan.at_end()? result: 0;
}

//////////////////////////////////
template <typename T, typename CharT, typename SetT>
inline T*
add(symbols<T, CharT, SetT>& table, CharT const* sym, T const& data)
{
    CharT const* last = sym;
    while (*last)
        last++;
    scanner<CharT const *> scan(sym, last);
    if (table.find(scan))
        return 0;               // symbol already contained in symbol table
    table.add(sym, last, data);
    return table.find(scan);    // refind the inserted symbol
}

///////////////////////////////////////////////////////////////////////////////
}} // namespace boost::spirit

#if defined(BOOST_MSVC)
#pragma warning(default : 4355)
#endif

#endif
