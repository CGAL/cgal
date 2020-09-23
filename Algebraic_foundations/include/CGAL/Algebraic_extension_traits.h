// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
        : public CGAL::cpp98::unary_function<Type,Type> {
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
        : public CGAL::cpp98::unary_function<Type,Type> {
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
