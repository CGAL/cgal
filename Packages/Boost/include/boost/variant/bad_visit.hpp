//-----------------------------------------------------------------------------
// boost variant/bad_visit.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002-2003
// Eric Friedman
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_VARIANT_BAD_VISIT_HPP
#define BOOST_VARIANT_BAD_VISIT_HPP

#include <exception>

namespace boost {

//////////////////////////////////////////////////////////////////////////
// class bad_visit
//
// Exception thrown when a visitation attempt via apply_visitor fails due
// to invalid visited subtype or contents.
//
struct bad_visit
    : std::exception
{
public: // std::exception interface

    virtual const char * what() const throw()
    {
        return "boost::bad_visit: "
               "failed visitation using boost::apply_visitor";
    }

};

} // namespace boost

#endif // BOOST_VARIANT_BAD_VISIT_HPP
