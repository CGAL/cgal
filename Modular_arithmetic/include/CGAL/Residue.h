// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Hemmer

/*! \file CGAL/Residue.h
    \brief Defines the class CGAL::Residue and CGAL::Modular_traits.
 
    Provides the \c CGAL::Modular_traits specialization for the build in number 
    types. 
*/

#ifndef CGAL_RESIDUE_H
#define CGAL_RESIDUE_H 1

#include <CGAL/basic.h>
#include <CGAL/Modular_arithmetic/Residue_type.h>
#include <CGAL/Coercion_traits.h>

namespace CGAL {


/*! \brief Specialization of CGAL::NT_traits for \c Residue, which is a model
 * of the \c Field concept. 
 * \ingroup CGAL_NT_traits_spec
 */
template <>
class Algebraic_structure_traits<Residue>
    : public Algebraic_structure_traits_base< Residue ,Field_tag >{
public: 
    typedef CGAL::Tag_true Is_exact; 
};

CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short,CGAL::Residue)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int  ,CGAL::Residue)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long ,CGAL::Residue)

} //namespace CGAL

#endif //#ifnedef CGAL_RESIDUE_H 1
 
