// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================


/*! \file CGAL/Algebraic_extension_traits.h
 *  \brief Defines traits class CGAL::Algebraic_extension_traits. 
*/

#ifndef CGAL_ALGEBRAIC_NUMBER_TRAITS_H
#define CGAL_ALGEBRAIC_NUMBER_TRAITS_H 1

#include <numeric> // for std::accumulate
#include <functional> // for std::unary_function
#include <CGAL/tags.h>
#include <CGAL/Algebraic_structure_traits.h>

namespace CGAL {

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
        : public std::unary_function<Type,Type> {
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
        : public std::unary_function<Type,Type> {
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

} //namespace CGAL

#endif // NiX_ALGEBRAIC_NUMBER_TRAITS_H
// EOF
