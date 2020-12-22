// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_DEGENERACY_TESTER_BINDERS_H
#define CGAL_VORONOI_DIAGRAM_2_DEGENERACY_TESTER_BINDERS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=======================================================================

template<class VDA>
class Edge_rejector_binder
{
private:
  typedef typename VDA::Accessor::Edge_rejector  ER;

public:
  typedef typename ER::result_type               result_type;

  Edge_rejector_binder(const VDA* vda = nullptr) : vda_(vda) {}

  template<class A>
  bool operator()(const A& a) const {
    CGAL_precondition( vda_ != nullptr );
    return vda_->edge_rejector()(vda_->dual(), a);
  }

private:
  const VDA* vda_;
};

//=======================================================================

template<class VDA>
class Face_rejector_binder
{
private:
  typedef typename VDA::Accessor::Face_rejector  FR;

public:
  typedef typename FR::result_type               result_type;

  Face_rejector_binder(const VDA* vda = nullptr) : vda_(vda) {}

  template<class A>
  bool operator()(const A& a) const {
    CGAL_precondition( vda_ != nullptr );
    return vda_->face_rejector()(vda_->dual(), a);
  }

private:
  const VDA* vda_;
};

//=======================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_DEGENERACY_TESTER_BINDERS_H
