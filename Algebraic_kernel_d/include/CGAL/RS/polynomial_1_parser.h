// Copyright (c) 2007-2010 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>
//
// (the first version of this file was written by Elias Tsigaridas)

#ifndef CGAL_RS_PARSER_1_H
#define CGAL_RS_PARSER_1_H

#include <iostream>
#include <string>
#include <algorithm>

// if boost version is 1.38 or newer, we use the new version of spirit
#include <boost/version.hpp>
#if BOOST_VERSION >= 103800
#include <boost/spirit/include/classic_core.hpp>
#define CGAL_BOOST_SPIRIT boost::spirit::classic
#else
#include <boost/spirit/core.hpp>
#define CGAL_BOOST_SPIRIT boost::spirit
#endif

namespace CGAL{

// The semantics.
// In this class we store the current pair and the whole result;
struct Semantics
{
    typedef std::pair< int, std::string >    pair_t;

    std::string                     variable;
    pair_t                          current;
    std::vector< pair_t >           result;
    std::string                     sign;

    Semantics( const std::string&  var) : variable( var) { }

    void set_variable( const std::string& var)
        {
            variable = var;
        }

    std::string get_variable() const
        {
            return variable;
        }

    bool is_first() const
        {
            return current == pair_t(-1, "F");
        }

    void add()
        {
            result.push_back( current);
            current = pair_t(-1, "F");
        }

    void init()
        {
            result.clear();
            current = pair_t(-1, "F");
            sign = "";
        }

};

//-------------    Semantic actions -----------------------------------

// Add a new monomial at the result
struct AddMonomial
{
    void operator()( char const* first, char const* last) const
        {
            if (!sem.is_first()) {
                sem.add();
            }
        }

    AddMonomial( Semantics& sem_) : sem( sem_) { }
    Semantics& sem;
};

// Set the coefficient of the current monomial
struct SetCoeff
{
    void operator() ( char const* first, char const* last) const
        {
            sem.current.second = std::string( first, last);
            if (sem.current.first == -1) {
                sem.current.first = 0;
            }
        }

    SetCoeff( Semantics& sem_) : sem( sem_) { }
    Semantics& sem;
};

// Set the exponent of the current monomial
struct SetExp
{
    void
    operator()( const char* first, const char* last) const
        {
            sem.current.first = 1;
            if (sem.current.second == "F") sem.current.second = "1";
        }

    void operator()( unsigned num) const
        {
            sem.current.first = num;
            if (sem.current.second == "F") sem.current.second = "1";
        }

    SetExp( Semantics& sem_) : sem( sem_) { }
    Semantics& sem;
};

// Set the sign of the current monomial
struct SetSign
{
    void
    operator()( char const* first, char const* last) const
        {
            std::string str( first, last);
            sem.sign = str;
        }

    SetSign( Semantics& sem_) : sem( sem_) { }
    Semantics& sem;
};

// Adjust the sign of the coefficient of the current monomial
struct AdjustCoeff
{
    void
    operator()( char const* first, char const* last) const
        {
            if (sem.sign == "-") {
                if (sem.current.second[0] == '-') {
                    sem.current.second[0] = '+';
                } else if (sem.current.second[0] == '+') {
                    sem.current.second[0] = '-';
                } else {
                    sem.current.second = '-' + sem.current.second;
                }
            }
        }

    AdjustCoeff( Semantics& sem_) : sem( sem_) { }
    Semantics& sem;
};

//----------------  Grammar -------------------------------------

// The grammar of the polynomial parser.
// The coefficient of the polynomial can also be rationals.
struct polynomial_p : public CGAL_BOOST_SPIRIT::grammar<polynomial_p>
{

    template < typename ScannerT >
    struct definition
    {
        definition(polynomial_p const& self)
            {
                CGAL_BOOST_SPIRIT::chlit<> lpar('(');
                CGAL_BOOST_SPIRIT::chlit<> rpar(')');

                CGAL_BOOST_SPIRIT::chlit<> plus('+');
                CGAL_BOOST_SPIRIT::chlit<> minus('-');
                CGAL_BOOST_SPIRIT::chlit<> mul('*');

                CGAL_BOOST_SPIRIT::strlit<>
                        varX(self.sem.get_variable().c_str());

                poly = (
                        (monomial | smonomial) [ AddMonomial( self.sem) ] >>
                        *( smonomial [ AddMonomial( self.sem) ] )
                       );

                monomial
                    = ( unumber [ SetCoeff( self.sem) ] >> !( mul >>  X ) )
                    | ( lpar >> number[ SetCoeff( self.sem) ] >> rpar >> !( mul >> X) )
                    | ( X >> mul >> number [ SetCoeff( self.sem) ] )
                    | ( X >> mul >> lpar >> number [ SetCoeff( self.sem) ] >> rpar )
                    | X;

                smonomial = ( ( plus | minus ) [ SetSign( self.sem) ] >>
                             monomial) [ AdjustCoeff( self.sem) ];

                X = ( (varX)[ SetExp( self.sem) ]  >>
                      !(
                        CGAL_BOOST_SPIRIT::ch_p('^') >>
                        ( (CGAL_BOOST_SPIRIT::uint_p) [ SetExp( self.sem) ]
                          | (lpar >>
                             CGAL_BOOST_SPIRIT::uint_p [ SetExp( self.sem) ] >>
                             rpar))
                      )
                    );

                unumber = +(CGAL_BOOST_SPIRIT::digit_p) >>
                          !(CGAL_BOOST_SPIRIT::ch_p("/") >>
                            +(CGAL_BOOST_SPIRIT::digit_p));
                number = ( (plus | minus) >> unumber );
            }

        CGAL_BOOST_SPIRIT::rule<ScannerT>
                poly, monomial, smonomial, X, unumber, number;

        CGAL_BOOST_SPIRIT::rule<ScannerT> const& start() const
            {
                return poly;
            }
    };

    polynomial_p( Semantics& sem_) : sem(sem_) { }

    Semantics& sem;
};

// The polynomial parser.
// The constructor takes the name of the variable.
struct Polynomial_parser_1
{
protected:
    // internal::Semantics                 sem;
    Semantics                           sem;
    // internal::polynomial_p              poly_p;
    polynomial_p                        poly_p;
    CGAL_BOOST_SPIRIT::parse_info<>     info;


// The default conversion is to int
    struct Converter
    {
        int operator()( const std::string str) const
            {
                return atoi( str.c_str());
            }
    };


public:
    Polynomial_parser_1( const std::string& var = "x")
        : sem( var), poly_p( sem)
        {
        }

    void parse( const std::string& in_str)
        {
            // Eliminate spaces
            std::string str(in_str);
            while (str.find(" ") < str.size() )
                {
                    size_t pos = str.find(" ");
                    str.erase( pos, 1);
                }

            sem.init();
            info = CGAL_BOOST_SPIRIT::parse(str.c_str(),
                                            poly_p,
                                            CGAL_BOOST_SPIRIT::space_p);
            std::sort( sem.result.begin(), sem.result.end());
        }

    bool is_correct() const
        {
            return info.full;
        }

    std::string error() const
        {
            if (!is_correct()) {
                return info.stop;
            }
            return "";
        }


    std::string get_variable() const
        {
            return sem.get_variable();
        }

    void set_variable( const std::string& var)
        {
            sem.set_variable( var);
        }



    template < typename OutputIterator,
               typename CONVERTER>
    OutputIterator
    result( OutputIterator oi, CONVERTER conv) const
        {
            //for (int i = 0; i < sem.result.size(); ++i) {
            //  std::cout << "(" << sem.result[i].first << ", " <<
            //    sem.result[i].second << ") ";
            //}
            //std::cout << std::endl;

            for (unsigned i = 0, j = 0; j < sem.result.size(); ++i) {
                if (sem.result[j].first == static_cast<int>( i)) {
                    *oi++ = conv( sem.result[j].second);
                    ++j;
                } else {
                    *oi++ = conv( "0");
                }

            }
            return oi;
        }


    template < typename OutputIterator >
    OutputIterator
    result( OutputIterator oi) const
        {
            return this->result( oi, Converter());
        }

};

template < typename T >
struct Convert_to
{
    T
    operator()( const std::string str) const {
        return static_cast<T>( atoi( str.c_str()));
    }
};

struct Convert_to_Gmpz
{
    Gmpz
    operator()( const std::string str) const {
        return Gmpz(str.c_str());
    }
};

} // namespace CGAL

#undef CGAL_BOOST_SPIRIT

#endif  // CGAL_RS_PARSER_1_H
