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

#ifndef CGAL_MODULAR_TRAITS_H
#define CGAL_MODULAR_TRAITS_H 1

#include <CGAL/tags.h>

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

