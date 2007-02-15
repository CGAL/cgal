// Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Elias Tsigaridas <Elias.Tsigaridas@loria.fr>

#ifndef AK_INTERNAL_PARSERS_H
#define AK_INTERNAL_PARSERS_H


//#include <AK/init.h>

#include <boost/spirit/core.hpp>
#include <iostream>
#include <string>
#include <algorithm>


CGAL_BEGIN_NAMESPACE
  
namespace spirit = boost::spirit;

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
struct polynomial_p : public spirit::grammar<polynomial_p>
{

    template < typename ScannerT >
    struct definition 
    {
	definition(polynomial_p const& self) 
	    {
		spirit::chlit<> lpar('(');
		spirit::chlit<> rpar(')');
      
		spirit::chlit<> plus('+');
		spirit::chlit<> minus('-');
		spirit::chlit<> mul('*');
      
		spirit::strlit<> varX( self.sem.get_variable().c_str());
      
		poly = (
		    (monomial | smonomial) [ AddMonomial( self.sem) ]
		    >> *( smonomial [ AddMonomial( self.sem) ] )
		    );
      
		monomial 
		    = ( unumber [ SetCoeff( self.sem) ] >> !( mul >>  X ) )
		    | ( lpar >> number[ SetCoeff( self.sem) ] >> rpar  >> !( mul >>  X) )
		    | ( X >> mul >> number [ SetCoeff( self.sem) ] )
		    | ( X >> mul >> lpar >> number [ SetCoeff( self.sem) ] >> rpar )
		    | X;
	    
		smonomial = ( ( plus | minus ) [ SetSign( self.sem) ] >> monomial) [ AdjustCoeff( self.sem) ];
	    
		X = ( (varX)[ SetExp( self.sem) ]  >> 
		      !(
			  spirit::ch_p('^') >> (  (spirit::uint_p) [ SetExp( self.sem) ] |
						  (lpar >> spirit::uint_p [ SetExp( self.sem) ] >> rpar) )
			  )
		    );

		unumber = +(spirit::digit_p) >> !( spirit::ch_p("/") >> +(spirit::digit_p) );
		number = ( (plus | minus) >> unumber );
	    }
    
	spirit::rule<ScannerT> poly, monomial, smonomial, X, unumber, number;

	spirit::rule<ScannerT> const& start() const
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
    /*internal::*/Semantics              sem;
    /*internal::*/polynomial_p           poly_p;
    boost::spirit::parse_info<>      info;
  

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
	    info = boost::spirit::parse(str.c_str(), poly_p, boost::spirit::space_p);
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
	    //  std::cout << "(" << sem.result[i].first << ", " << sem.result[i].second << ") ";
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

struct The_Convert_to
{
    Gmpz
    operator()( const std::string str) const {
	return Gmpz(str.c_str());
    }
};



// int
// main()
// {

//   // Initialize the parser. The variable is 'x'
//   Polynomial_parser_1 parser;
  
//   // The string that we will parse
//   std::string tstr("-  3 +  2*x -( -  23)*x^8 - (- 2328723)  *x^ 2 ");

//   // Parse the string
//   parser.parse( tstr);

//   // If the parse was succeful
//   if (parser.is_correct()) {

//     // The coefficients of the polynomial
//     typedef double NT;
//     std::vector< NT > Coeff;
    
//     // We get the coefficients from the parser.
//     // Notice that we use our convertor. If we didn't supply a convertor
//     // then the default convertor would be used (to ints)
//     parser.result( std::back_inserter( Coeff), MyConverter());
    
//     // Output the result
//     std::cout.precision( 20);
    
//     std::cout << "Coeff: ";
//     std::copy( Coeff.begin(), Coeff.end(), std::ostream_iterator<NT>( std::cout, " "));
//     std::cout << std::endl; 
//   } else {
//     // Something went wrong
//     std::cout << "Failure..." << std::endl; 
//     // The error was at...
//     std::cout << "at: " << parser.error() << std::endl; 
//   }
  
//   std::cout << "End of story..." << std::endl; 
  
//   return 0;
// }

    
CGAL_END_NAMESPACE

#endif // AK_INTERNAL_PARSERS_H  
