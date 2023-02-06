// Copyright (c) 2014  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando de Goes, Pierre Alliez, Ivo Vigan, Clément Jamin

#ifndef CGAL_RECONSTRUCTION_VERTEX_BASE_2_H_
#define CGAL_RECONSTRUCTION_VERTEX_BASE_2_H_

#include <CGAL/license/Optimal_transportation_reconstruction_2.h>


#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/OTR_2/Sample.h>

namespace CGAL {
namespace OTR_2 {

/// The Reconstruction_vertex_base_2 class is the default
/// vertex class of the Reconstruction_triangulation_2 class.
///
/// - Each vertex stores a Sample as well as the corresponding relocated point.
///
/// @param Traits_  Geometric traits class
/// @param Vb       Vertex base class, model of TriangulationVertexBase_2.
template < class Traits_, class Vb = Triangulation_vertex_base_2<Traits_> >
class Reconstruction_vertex_base_2 : public Vb
{

  /// \cond SKIP_IN_MANUAL

public:
  typedef Vb Base;
  typedef typename Traits_::FT        FT;
  typedef OTR_2::Sample<Traits_>      Sample_;
  typedef typename Traits_::Point_2   Point;
  typedef typename Base::Face_handle  Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Base::template Rebind_TDS<TDS2>::Other Vb2;
    typedef Reconstruction_vertex_base_2<Traits_,Vb2> Other;
  };

private:
  int       m_id;
  bool      m_pinned;
  int       m_sample;
  Point     m_relocated;
  FT        m_relevance;


public:
  Reconstruction_vertex_base_2()
  : Base(),
    m_id(-1),
    m_pinned(false),
    m_sample(-1),
    m_relevance(0)
{
}

  Reconstruction_vertex_base_2(const Point & p)
  : Base(p),
    m_id(-1),
    m_pinned(false),
    m_sample(-1),
    m_relevance(0)
  {
  }

  Reconstruction_vertex_base_2(Face_handle f)
  : Base(f),
    m_id(-1),
    m_pinned(false),
    m_sample(-1),
    m_relevance(0)
  {
  }

  Reconstruction_vertex_base_2(const Point & p, Face_handle f)
  : Base(p, f),
    m_id(-1),
    m_pinned(false),
    m_sample(-1),
    m_relevance(0)
  {
  }

  ~Reconstruction_vertex_base_2() { }

  int  id() const { return m_id; }
  int& id() { return m_id; }

  bool  pinned() const { return m_pinned; }
  bool& pinned() { return m_pinned; }

  FT relevance() const { return m_relevance; }
  void set_relevance(FT relevance) { m_relevance = relevance; }

  int sample() const { return m_sample; }
  void set_sample(int sample) { m_sample = sample; }

  const Point& relocated() const { return m_relocated; }
  Point& relocated() { return m_relocated; }

  bool  has_sample_assigned() const { return sample() != -1; }
};
//---------------STRUCT LESS VERTEX_HANDLE---------------------
template <class T>
struct less_Vertex_handle
{
  bool operator() (const T& a, const T& b) const
  {
    return (a->id() < b->id());
  }
};


/// \endcond

} } //end namespaces

#endif // CGAL_RECONSTRUCTION_VERTEX_BASE_2_H_
