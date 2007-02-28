// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
// 
//
// Author(s)     :  
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_POLYNOMIAL_COERCION_TRAITS_H
#define CGAL_POLYNOMIAL_COERCION_TRAITS_H

CGAL_BEGIN_NAMESPACE

// COERCION_TRAITS BEGIN 

//Coercion_traits_polynomial-----------------------------------
// If there is a Polynomial_traits, valid for more than one Polynomial
// class this part should be adapted, using a Polynomial_traits 
// and the nesting_depth 
template <class A,class B>
class Coercion_traits_for_level<Polynomial<A>, Polynomial<B>, CTL_POLYNOMIAL >{
    typedef Coercion_traits<A,B> CT;            
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Polynomial<typename CT::Type> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const Polynomial<A>& poly) const { 
            typename CT::Cast cast; 
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                    ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
        Type operator()(const Polynomial<B>& poly) const {  
            typename CT::Cast cast;  
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                    ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
    }; 
};
        
template <class A,class B>
class Coercion_traits_for_level<Polynomial<A>,B ,CTL_POLYNOMIAL >{
    typedef Coercion_traits<A,B> CT;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    typedef Polynomial<typename CT::Type> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const Polynomial<A>& poly) const {
            typename CT::Cast cast;
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                       ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
        Type operator()(const B& x) const {
            typename CT::Cast cast;
            return Type(cast(x));
        } 
    };                                                        
}; 
template <class A,class B> 
class Coercion_traits_for_level<B,Polynomial<A>,CTL_POLYNOMIAL  >
    :public Coercion_traits_for_level<Polynomial<A>,B,CTL_POLYNOMIAL >
{};

// COERCION_TRAITS END

CGAL_END_NAMESPACE

#endif // CGAL_POLYNOMIAL_COERCION_TRAITS_H
