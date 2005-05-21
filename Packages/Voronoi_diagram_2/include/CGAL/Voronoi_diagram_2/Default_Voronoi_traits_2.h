// Copyright (c) 2005 University of Crete (Greece).
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
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_DEFAULT_VORONOI_TRAITS_2_H
#define CGAL_VORONOI_DIAGRAM_2_DEFAULT_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Cached_degeneracy_testers.h>
#include <CGAL/Handle_for_virtual.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class DG, class ET, class FT>
class Default_Voronoi_traits_2
{
 private:
  typedef Default_Voronoi_traits_2<DG,ET,FT>   Self;

 public:
  typedef DG  Dual_graph;
  typedef ET  Edge_degeneracy_tester;
  typedef FT  Face_degeneracy_tester;

  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


//=========================================================================
//=========================================================================

template<class DG, class ET, class FT>
class Default_cached_Voronoi_traits_2
{
 private:
  typedef ET  Edge_degeneracy_tester_base;
  typedef FT  Face_degeneracy_tester_base;

  typedef Default_cached_Voronoi_traits_2<DG,ET,FT>  Self;

 public:
  typedef DG           Dual_graph;

  typedef Cached_edge_degeneracy_tester<Edge_degeneracy_tester_base>
  Edge_degeneracy_tester;

  typedef Cached_face_degeneracy_tester<Face_degeneracy_tester_base,
					Edge_degeneracy_tester>
  Face_degeneracy_tester;


 public:
  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


//=========================================================================
//=========================================================================


template<class DG, class ET, class FT>
class Default_ref_counted_Voronoi_traits_2
{
 private:
  typedef ET  Edge_degeneracy_tester_base;
  typedef FT  Face_degeneracy_tester_base;

  typedef Default_ref_counted_Voronoi_traits_2<DG,ET,FT>  Self;

 public:
  typedef DG           Dual_graph;

  typedef Ref_counted_edge_degeneracy_tester<Edge_degeneracy_tester_base>
  Edge_degeneracy_tester;

  typedef Ref_counted_face_degeneracy_tester<Face_degeneracy_tester_base,
					     Edge_degeneracy_tester>
  Face_degeneracy_tester;

 public:
  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};

//=========================================================================
//=========================================================================


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_DEFAULT_VORONOI_TRAITS_2_H
