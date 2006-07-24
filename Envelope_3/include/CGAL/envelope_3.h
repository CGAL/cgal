// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Ron Wein     <wein@post.tau.ac.il>

#ifndef ENVELOPE_3_H
#define ENVELOPE_3_H

#include <CGAL/Envelope_pm_dcel.h>
#include <CGAL/Envelope_caching_traits_3.h>
#include <CGAL/Envelope_overlay_functor.h>

#include <iostream>

CGAL_BEGIN_NAMESPACE

template <class Traits3> class Envelope_diagram_2 :
  public Arrangement_2<Envelope_caching_traits_3<Traits>,
                       Envelope_pm_dcel<Envelope_caching_traits_3<Traits>,
                                        typename Envelope_caching_traits_3<Traits>::Xy_monotone_surface_3> >
{
};

template <class InputIterator, class Traits3>
void lower_envelope_3 (InputIterator begin, InputIterator end,
                       Envelope_diagram_2<Traits3>& min_diagram);

template <class InputIterator, class Traits3>
void upper_envelope_3 (InputIterator begin, InputIterator end,
                       Envelope_diagram_2<Traits3>& max_diagram);

CGAL_END_NAMESPACE

#endif
