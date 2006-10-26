// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================

#ifndef CGAL_NEEDS_PARENTHESES_AS_PRODUCT_H
#define  CGAL_NEEDS_PARENTHESES_AS_PRODUCT_H

#include <CGAL/IO/io.h>
 
CGAL_BEGIN_NAMESPACE

/*! 
 * oformat flag for parentheses if needed for a coefficient
 */
class Parens_as_product_tag {};

/*! \ingroup NiX_io_parens
 *  \brief Decides whether this number requires parentheses 
 *  in case it appears within a produkt.
 */
template <class NT>
struct Needs_parens_as_product{
    bool operator()(const NT& x){ return true; }
};

/*! \ingroup NiX_io_parens
 *  \brief Decides whether this number requires parentheses 
 *  in case it appears within a produkt.
 */
template <class NT>
inline bool needs_parens_as_product(const NT& x){
    typedef Needs_parens_as_product<NT> NPAP;
    return NPAP()(x);
}

/*! \ingroup NiX_io_parens
 *  \brief \c oformat \c Output_rep specialized for
 *  \c NiX::Parens_as_product_tag
 */
template <class T>
class Output_rep<T, Parens_as_product_tag> {
    const T& t;
public:
    Output_rep(const T& tt) : t(tt) {}
    std::ostream& operator () (std::ostream& out) const { 
        if ( ::CGAL::needs_parens_as_product(t)) {
            return out << "(" << oformat(t) << ")";
        } else {
            return out << oformat(t);
        }
    }
};


// built-in number types:
template <> struct Needs_parens_as_product<short>{
    bool operator()(const short& x){return x < short(0);} 
};
template <> struct Needs_parens_as_product<int>{
    bool operator()(const int& x){return x < int(0);} 
};
template <> struct Needs_parens_as_product<long>{
    bool operator()(const long& x){return x < long(0);} 
};
template <> struct Needs_parens_as_product<long long>{
    bool operator()(const long long& x){return x < (long long)(0);} 
};
template <> struct Needs_parens_as_product<float>{
    bool operator()(const float& x){return x < float(0);} 
};
template <> struct Needs_parens_as_product<double>{
    bool operator()(const double& x){return x < double(0);} 
};
template <> struct Needs_parens_as_product<long double>{
    bool operator()(const long double& x){return x < (long double)(0);} 
};

CGAL_END_NAMESPACE

#endif  //CGAL_NEEDS_PARENTHESES_AS_PRODUCT_H
