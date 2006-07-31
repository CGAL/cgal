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

#include <CGAL/Envelope_divide_and_conquer_3.h>
#include <CGAL/Envelope_pm_dcel.h>
#include <CGAL/Envelope_caching_traits_3.h>
#include <CGAL/Envelope_overlay_functor.h>

#include <iostream>

CGAL_BEGIN_NAMESPACE

template <class Traits_3> class Envelope_diagram_2 :
  public Arrangement_2<Envelope_caching_traits_3<Traits_3>,
                       Envelope_pm_dcel<Envelope_caching_traits_3<Traits_3>,
                                        typename Envelope_caching_traits_3<Traits_3>::Xy_monotone_surface_3> >
{
  public:
  typedef Arrangement_2<Envelope_caching_traits_3<Traits_3>,
                       Envelope_pm_dcel<Envelope_caching_traits_3<Traits_3>,
                                        typename Envelope_caching_traits_3<Traits_3>::Xy_monotone_surface_3> > Base;

  Envelope_diagram_2()
  {}

  Envelope_diagram_2(Envelope_caching_traits_3<Traits_3>* tr): Base(tr)
  {}

  ~Envelope_diagram_2()
  {}
};

template <class InputIterator, class Traits_>
void lower_envelope_3 (InputIterator begin, InputIterator end,
                       Envelope_diagram_2<Traits_>& min_diagram)
{
  typedef Envelope_caching_traits_3<Traits_>  Traits_3;
  typedef Envelope_divide_and_conquer_3<Traits_3, typename Envelope_diagram_2<Traits_>::Base >
    Envelope_alg;
  Traits_3 traits(LOWER);
  Envelope_alg envelope_algorithm(&traits);
  envelope_algorithm.construct_lu_envelope(begin, end, min_diagram);
}

template <class InputIterator, class Traits_>
void upper_envelope_3 (InputIterator begin, InputIterator end,
                       Envelope_diagram_2<Traits_>& max_diagram)
{
  typedef Envelope_caching_traits_3<Traits_>  Traits_3;
  typedef Envelope_divide_and_conquer_3<Traits_3, typename Envelope_diagram_2<Traits_>::Base >
    Envelope_alg;
  Traits_3 traits(UPPER);
  Envelope_alg envelope_algorithm(&traits);
  envelope_algorithm.construct_lu_envelope(begin, end, max_diagram);
}

CGAL_END_NAMESPACE

#endif
