//  (C) Copyright Gennadiy Rozental 2004.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : common code used by any agent serving as XML printer
// ***************************************************************************

#ifndef BOOST_TEST_XML_PRINTER_HPP_071894GER
#define BOOST_TEST_XML_PRINTER_HPP_071894GER

// Boost.Test
#include <boost/test/detail/basic_cstring/basic_cstring.hpp>
#include <boost/test/detail/fixed_mapping.hpp>

// BOOST
#include <boost/config.hpp>

// STL
#include <iostream>

namespace boost {

namespace unit_test {

namespace ut_detail {

// ************************************************************************** //
// **************                  xml_printer                 ************** //
// ************************************************************************** //

class xml_printer {
    static inline std::ostream& print_escaped( std::ostream& where_to, const_string value )
    {
        static fixed_mapping<char,char const*> char_type(
            '<' , "lt",
            '>' , "gt",
            '&' , "amp",
            '\'', "apos" ,
            '"' , "quot",

            0
            );

        for( const_string::iterator it = value.begin(); it != value.end(); ++it ) {
            char const* ref = char_type[*it];

            if( ref )
                where_to << '&' << ref << ';';
            else
                where_to << *it;
        }

        return where_to;
    }
protected:
    static inline std::ostream&    print_attr_value( std::ostream& where_to, const_string value )
    {
        where_to << "=\"";
        return print_escaped( where_to, value ) << '"';
    }

    template<typename T>
    static inline std::ostream&    print_attr_value( std::ostream& where_to, T const& value )
    {
        return where_to << "=\"" << value << '"';
    }

    static inline std::ostream&    print_pcdata( std::ostream& where_to, const_string value )
    {
        return print_escaped( where_to, value );
    }
};

} // namespace ut_detail

} // namespace unit_test

} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/11/20 10:52:21  spion
//  Initial revision
//
//  Revision 1.2  2004/08/04 02:50:27  rogeeff
//  darwin workarounds
//
//  Revision 1.1  2004/07/19 12:22:15  rogeeff
//  shared xml printer utils
//
// ***************************************************************************

#endif // BOOST_TEST_XML_PRINTER_HPP_071894GER
