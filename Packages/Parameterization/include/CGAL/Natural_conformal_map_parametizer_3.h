// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_NATURAL_CONFORMAL_MAP_PARAMETIZER_3_H
#define CGAL_NATURAL_CONFORMAL_MAP_PARAMETIZER_3_H

#include <CGAL/Parametizer_traits_3.h>

CGAL_BEGIN_NAMESPACE


/// Class Natural_conformal_map_parametizer_3
/// implements Natural Conformal Map parameterization algorithm (Alliez et al).
/// No need to map the surface's border onto a convex polygon
/// but 1 to 1 mapping not guaranteed.
/// This is a conformal parameterization, i.e. it attempts to preserve angles.
///
/// Concept: Model of the ParametizerTraits_3 concept.
///
/// Design pattern:
/// ParametizerTraits_3 models are Strategies (see [GOF95]): they implement
/// a strategy of surface parameterization for models of MeshAdaptor_3.

/// NOT YET IMPLEMENTED

template <...>
class Natural_conformal_map_parametizer_3
    : public Parametizer_traits_3<...>
{
/// Public stuff
public:

/// Protected stuff
protected:

/// Private stuff
private:

};


CGAL_END_NAMESPACE

#endif //CGAL_NATURAL_CONFORMAL_MAP_PARAMETIZER_3_H