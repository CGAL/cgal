// Copyright(c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_DRAW_CT2_H
#define CGAL_DRAW_CT2_H

#include <CGAL/license/Triangulation_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_2/internal/In_domain.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>

namespace CGAL
{

  template<class InDomainPmap>
  CGAL::IO::Color color_of_face(Facet_const_handle fh, InDomainPmap ipm)
  { return get(ipm, fh)? CGAL::IO::yellow() : CGAL::IO::white(); }

  CGAL::IO::Color color_of_edge(Edge_const_handle eh)
  { return t2.is_constrained(*eh)? CGAL::IO::green() : CGAL::IO::black(); }


// Specialization of draw function.
#define CGAL_T2_TYPE CGAL::Constrained_triangulation_2<Gt, Tds, Itag>

template<class Gt, class Tds, class Itag, class InDomainPmap>
void draw(const CGAL_T2_TYPE& at2, InDomainPmap ipm)
{
}


template<class Gt, class Tds, class Itag>
void draw(const CGAL_T2_TYPE& at2)
{
  internal::In_domain<CGAL_T2_TYPE> in_domain;
  draw(at2, in_domain);
}

#undef CGAL_T2_TYPE

} // End namespace CGAL

#else

namespace CGAL {
// Specialization of draw function.
#define CGAL_T2_TYPE CGAL::Constrained_triangulation_2<Gt, Tds, Itag>

template<class Gt, class Tds, class Itag, class InDomainPmap>
void draw(const CGAL_T2_TYPE& ,
          InDomainPmap )
{}
#undef CGAL_T2_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_CT2_H
