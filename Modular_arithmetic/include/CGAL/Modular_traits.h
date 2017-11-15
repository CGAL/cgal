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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Michael Hemmer

#ifndef CGAL_MODULAR_TRAITS_H
#define CGAL_MODULAR_TRAITS_H 1

#include <CGAL/basic.h>
#include <CGAL/Residue.h>
#include <CGAL/Modular_arithmetic/Residue_type.h>
#include <vector>


namespace CGAL { 


/*! \ingroup CGAL_Modular_traits_spec 
    \brief A model of concept ModularTraits. 
    
    This is the definition of general class template, 
    for unsupported types. Note that this support is optional. 
    \see CGAL_Modular_traits_spec for supported types. 
 */
 
template<class NT_>
class Modular_traits{
public: 
    typedef NT_ NT;
    typedef ::CGAL::Tag_false Is_modularizable;
    typedef ::CGAL::Null_functor Residue_type;
    typedef ::CGAL::Null_functor Modular_image;  
    typedef ::CGAL::Null_functor Modular_image_representative;    
};

template <class NT>
inline
typename CGAL::Modular_traits<NT>::Residue_type 
modular_image(const NT& x){
    typename CGAL::Modular_traits<NT>::Modular_image modular_image;
    return modular_image(x);
}

}///namespace CGAL
#endif //#ifnedef CGAL_MODULAR_TRAITS_H 1
 
