// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================

/*! \file CGAL/Algebraic_extension_traits.h
 *  \brief Defines traits class CGAL::Algebraic_extension_traits. 
*/

#ifndef CGAL_ALGEBRAIC_NUMBER_TRAITS_H
#define CGAL_ALGEBRAIC_NUMBER_TRAITS_H 1

CGAL_BEGIN_NAMESPACE

template< class T >
class Algebraic_extension_traits {
public:
    //! \name Typedefs 
    //! the number type for which this instance has been instantiated
    typedef T Type;
    //! standard number types are not extended
    typedef CGAL::Tag_false Is_extended;
  
    //! computes the factor which normalizes a number to be integral after 
    //  multiplication
    class Normalization_factor 
        : public Unary_function<Type,Type> {
    private:
        static Type 
        normalization_factor(const Type&,Integral_domain_without_division_tag){
            return Type(1);
        }
        static Type 
        normalization_factor(const Type& a, Field_tag){
            return Type(1)/a;
        }
    public:
        //! determine normalization factor
        Type operator () (const Type& a) {
            CGAL_precondition(a != Type(0));
            typedef typename Algebraic_structure_traits<Type>::Algebraic_category
                Tag;
            return normalization_factor(a, Tag());
        }
    };
    
    class Denominator_for_algebraic_integers 
        : public Unary_function<Type,Type> {
    public: 
        //! determine normalization factor
        Type operator () (const Type&) {
            return Type(1);
        }
        
        template <class InputIterator>
        Type operator () (InputIterator, InputIterator) {
            return Type(1);
        }
    };
};

CGAL_END_NAMESPACE

#endif // NiX_ALGEBRAIC_NUMBER_TRAITS_H
// EOF
