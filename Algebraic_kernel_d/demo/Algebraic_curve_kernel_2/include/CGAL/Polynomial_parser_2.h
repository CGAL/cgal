// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*  \brief Defines functor CGAL::Polynomial_parser_2
 *  
 *  Parser to read bivariate polynomials in MAPLE format
 */

#ifndef CGAL_POLYNOMIAL_PARSER_2_H
#define CGAL_POLYNOMIAL_PARSER_2_H

#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <CGAL/ipower.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

class Parser_exception {

public:
    Parser_exception(std::string msg) :
        _m_error(msg) {
    }
        
    std::string get_message() const
    { return _m_error; }

private:
    std::string _m_error;
 
};

// coefficient validity check
template <class Poly2>
struct _Default_checker {

    // coefficient type
    typedef typename CGAL::Polynomial_traits_d<Poly2>::Innermost_coefficient_type NT;

    bool operator()(const NT& x) const {
        return true;
    }
};

} // namespace CGALi

/*! 
 *  \brief this functor implements parsing of bivariate polynomials
 * 
 * input format (y is outermost variable), e.g.:
 * (y-1)^4 + (-1)*y^3 + (x + y)^2 - (2123234523*x^2 - 2*y*y*x + 3*x*132123)^3
 * (y + x - 3)^3 = x + y - 123 + y*x^2
 */
template < typename Poly2, 
           typename ValidyChecker = CGAL::CGALi::_Default_checker<Poly2> >
struct Polynomial_parser_2
{
    //!\name public typedefs
    //!@{

    //! this instance's template argument
    typedef Poly2 Poly_2;
    //! this instance's second template argument
    typedef ValidyChecker Validy_checker;
    
    //! an instance of univariate polynomial
    typedef typename Poly_2::NT Poly_1;
    //! an instance of innermost coefficient
    typedef typename Poly_1::NT NT;

    //!@}
public: 
    /// \name Public methods
    //!@{

    //! \brief functor invokation operator
    //!
    //! if \c max_exp_ > 0 it defines the maximal allowed exponent 
    //! for input polynomial
    bool operator()(const std::string& in, Poly_2& poly,
        unsigned max_exp_ = -1)
    {
        max_exp = max_exp_;
        try {
            // remove white spaces from the string: look for all possible
            // whitespace characters
            std::string s(in);
            const char *whitespace = " \f\n\r\t\v";
            unsigned int cnt = s.find_first_of(whitespace,0);
            while(cnt != std::string::npos){
                s.erase(cnt, 1);
                cnt = s.find_first_of(whitespace, cnt);
            }
                        
            //to handle the case when there is one '=' sign
            //we compute sides separetely and then subtract one from another
            unsigned int loc;
            std::string rhs;
            if((loc=s.find('=',0)) != std::string::npos) {
                rhs = s.substr(loc+1); // right-hand side
                s.erase(loc);
            }
            poly = get_poly(s);
            if(loc != std::string::npos) {
                Poly_2 tmp = get_poly(rhs);
                poly -= tmp;
            }
        } 
        catch(CGAL::CGALi::Parser_exception ex) {
            std::cerr << "Parser error: " << ex.get_message() << std::endl;
            return false;
        }
        return true;
    }
    
    //!@}
protected: 
    /// \name protected methods
    //!@{ 
    
    //! given a string \c cstr of length \c len starting with an open
    //! parentheses returns the place marking the end of the corresponding
    //! closing parentheses
    unsigned int match_parenth(std::istream& is)
    {   
        int count = 0;
        unsigned int pos = 0, start = is.tellg(); // save pos of the first '('
        do {
            if(is.eof())
            return -1; // illegal number of parentheses
            char ch = is.peek();
            if(ch == '(')
                count++;
            else if(ch == ')') 
                count--;
            is.seekg(1, std::ios_base::cur);
            pos++;
        } while(count != 0); //j is one more than the matching ')'
        is.seekg(start);
        // pos is the number of chars including ( and )
        return pos - 1u; // return position of the last ')'
    }

    //! constructs {x/y}^exp and returns the result as bivariate polynomial 
    //! \c res \c var encodes a variable (x or y) 
    Poly_2 construct_xy(char var, unsigned int exp)
    {
        if(exp == 0) // just return 1
            return Poly_2(Poly_1(NT(1)));
        
        if(var == 'x'||var == 'X') { // otherwise construct monomial
            typename Poly_1::Vector v(exp+1, NT(0));
            v[exp] = NT(1);
            Poly_1 tmp(v.begin(), v.end());
            return Poly_2(tmp); // y^0*x^exp
        } 
        // else ys
        Poly_1 zero(NT(0)), one(NT(1));
        typename Poly_2::Vector v(exp+1, zero);
        v[exp] = one; // y^exp*x^0
        Poly_2 tmp(v.begin(), v.end());
        return tmp;
    }

    void get_basic_term(std::istringstream& is, Poly_2& res)
    {
        char var = 'x', ch = is.peek();
        //std::cout << "getbasicterm: " << is.tellg() << std::endl;
    
        Poly_2 tmp;
        NT coeff;
        unsigned int which_case = 0, power = 1;
        if(isdigit(ch)) {
            is >> CGAL::iformat(coeff);
            if(!checker(coeff))
                throw CGAL::CGALi::Parser_exception(
                    "Coefficient validity check failed");
            which_case = 0; // number        
        
        } else if(ch == 'x' || ch == 'X'|| ch == 'y'|| ch == 'Y'){
            which_case = 1;
            var = is.get();
            
        } else if(ch =='(') {
        
            // pos is the index of closing parentheses relative 
            // to the opening ')'
            unsigned int pos = match_parenth(is);
            if(pos == -1u)
                throw CGAL::CGALi::Parser_exception(
                    "Parentheses do not match in basic term");
        
            unsigned int beg = (unsigned int)is.tellg() + 1u;
            std::string t = is.str().substr(beg, pos-1);
           is.seekg(pos+1, std::ios_base::cur);
           //std::cout << "polynomial in parenth: " << t << "\n";
           //printf("next char to be read: %c\n", is.peek());
            tmp = get_poly(t);
            which_case = 2;
        } else 
            throw CGAL::CGALi::Parser_exception("Error in parsing basic term");
        
        // adjust i to be pointed to the next char
        if(is.peek() == '^') {
            is.ignore(1); // ignore one character
            if(!isdigit(is.peek()))
                throw CGAL::CGALi::Parser_exception(
                    "Incorrect power for basic term");
            is >> CGAL::iformat(power);
            if(power >= max_exp)
                throw CGAL::CGALi::Parser_exception(
                    "Power is too large for basic term");
        } 
        switch(which_case) {
        case 0:
            coeff = CGAL::ipower(coeff, static_cast< long >(power));
            tmp = Poly_2(Poly_1(coeff));
            break;
        case 1:
            tmp = construct_xy(var, power);
            break;
        case 2: // control degree overflow
            int degree = CGAL::total_degree(tmp);
            if(degree * power >= max_exp) 
                throw CGAL::CGALi::Parser_exception(
                    "Power is too large for polynomial in basic  term ");
            tmp = CGAL::ipower(tmp, static_cast< long >(power));    
        }
        res = tmp;
        // returned index points to the next syntactic object
        //std::cout << "next read pos: " << is.tellg() << "\n";
        //std::cout << "getbasicterm result: " << res << "\n";
    }
    
    void get_term(std::istringstream& is, Poly_2& res)
    {
        //std::cout << "getterm: " << is.tellg() << "\n";
        if(is.eof()) {
            res=Poly_2(Poly_1(NT(0)));
            return;
        }
        Poly_2 mul;
        get_basic_term(is, mul);
        
        char ch = is.peek();
        //printf("getterm next char to be read: %c\n", ch);
        //while(ind < len && cstr[ind] != '+' && cstr[ind] != '-') {
        while(!is.eof() && ch != '+' && ch != '-') {
            //Recursively get the basic terms till we reach the end or see
            // a '+' or '-' sign.
            if(ch == '*') 
                is.ignore(1);
            Poly_2 tmp;
            get_basic_term(is, tmp);
            mul *= tmp;
            ch = is.peek();
            //printf("getterm next char to be read: %c\n", ch);
        }
        res = mul;
        //std::cout << "getterm result: " << res << "\n";
    }
    
    Poly_2 get_poly(std::string &s)
    {
       //std::cout << "getpoly: " << s << "\n";
        unsigned int len = s.length();
        if(len == 0) // zero polynomial
            return Poly_2(Poly_1(NT(0)));
               
        len = s.length();
        std::istringstream is(s);
        Poly_2 res;
        // res will be the polynomial in which we accumulate the
        // sum and difference of the different terms.
        if(is.peek() == '-') {
            is.ignore(1);
            get_term(is, res);
            res *= Poly_1(NT(-1)); // negate polynomial
        } else 
            get_term(is, res);
        
        //  here ind either points to +/- or to the end of string        
        while(!is.eof()) {
            
            Poly_2 tmp;
            char ch = is.get(); // must be +/-
            get_term(is, tmp);
            if(ch == '+')
                res += tmp;
            else if(ch == '-')
                res -= tmp;
            else
                throw CGAL::CGALi::Parser_exception(
                    "Illegal character while parsing polynomial");
        }
        //std::cout << "getpoly result: " << res << "\n";
        return res;
    }

protected:

    Validy_checker checker; 
         
    //! a maximal exponent allowed for input polynomial
    unsigned max_exp;
    
};  // Polynomial_parser_2

CGAL_END_NAMESPACE

#endif // CGALPOLYNOMIAL_PARSER_2_H


