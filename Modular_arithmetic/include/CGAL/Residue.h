// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Hemmer

/*! \file CGAL/Residue.h
    \brief Defines the class CGAL::Residue and CGAL::Modular_traits.

    Provides the \c CGAL::Modular_traits specialization for the build in number
    types.
*/

#ifndef CGAL_RESIDUE_H
#define CGAL_RESIDUE_H 1

#include <CGAL/Modular_arithmetic/Residue_type.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/int.h>
#include <CGAL/Algebraic_structure_traits.h>

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

