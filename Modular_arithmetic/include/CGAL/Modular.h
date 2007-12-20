// Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// $URL:$
// $Id:$
//
// Author(s)     : Michael Hemmer

/*! \file CGAL/Modular.h
    \brief Defines the class CGAL::Modular and CGAL::Modular_traits.
 
    Provides the \c CGAL::Modular_traits specialization for the build in number 
    types. 
*/

#ifndef CGAL_MODULAR_H
#define CGAL_MODULAR_H 1

#include <CGAL/basic.h>
#include <CGAL/Modular_arithmetic/Modular_type.h>

CGAL_BEGIN_NAMESPACE

/*! \brief Specialization of CGAL::NT_traits for \c Modular, which is a model
 * of the \c Field concept. 
 * \ingroup CGAL_NT_traits_spec
 */
template <>
struct Algebraic_structure_traits<Modular>
    : public Algebraic_structure_traits_base< Modular ,Field_tag >{
    typedef CGAL::Tag_true Is_exact; 
};

CGAL_END_NAMESPACE

#endif //#ifnedef CGAL_MODULAR_H 1
 
