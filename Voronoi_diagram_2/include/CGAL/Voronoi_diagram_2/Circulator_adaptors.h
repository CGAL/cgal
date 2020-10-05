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

#ifndef CGAL_VORONOI_DIAGRAM_2_CIRCULATOR_ADAPTORS_H
#define CGAL_VORONOI_DIAGRAM_2_CIRCULATOR_ADAPTORS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/circulator_bases.h>
#include <CGAL/Voronoi_diagram_2/Handle_adaptor.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================

template<class Arg, class Circulator>
class Circulator_from_halfedge_adaptor
  : public Bidirectional_circulator_base<Arg>
{
 private:
  typedef Bidirectional_circulator_base<Arg>                Base;
  typedef Circulator_from_halfedge_adaptor<Arg,Circulator>  Self;
  typedef Handle_adaptor<Arg>                               Handle;

 public:
  Circulator_from_halfedge_adaptor() : cur_() {}
  Circulator_from_halfedge_adaptor(const Arg& he) : cur_(he) {}

  operator Handle() const {
    return Handle(cur_);
  }

  Circulator& operator++() {
    Circulator* ptr = static_cast<Circulator*>(this);
    ptr->increment();
    return *ptr;
  }

  Circulator& operator--() {
    Circulator* ptr = static_cast<Circulator*>(this);
    ptr->decrement();
    return *ptr;
  }

  Circulator operator++(int) {
    Circulator tmp(*this);
    (static_cast<Circulator*>(this))->operator++();
    return tmp;
  }

  Circulator operator--(int) {
    Circulator tmp(*this);
    (static_cast<Circulator*>(this))->operator--();
    return tmp;
  }

  typename Base::pointer   operator->() { return &cur_; }
  typename Base::reference operator*() { return cur_; }

  bool operator==(const Circulator& other) const {
    return cur_ == other.cur_;

  }

  bool operator!=(const Circulator& other) const {
    return cur_ != other.cur_;
  }

 protected:
  Arg cur_;
};

//=========================================================================
//=========================================================================

template<class Halfedge>
class Halfedge_around_vertex_circulator_adaptor
  : public Circulator_from_halfedge_adaptor<Halfedge,
      Halfedge_around_vertex_circulator_adaptor<Halfedge> >
{
 private:
  typedef Halfedge                                         Arg;
  typedef Halfedge_around_vertex_circulator_adaptor<Arg>   Self;
  typedef Circulator_from_halfedge_adaptor<Arg,Self>       Base;

  friend class Circulator_from_halfedge_adaptor<Arg,Self>;

 public:
  Halfedge_around_vertex_circulator_adaptor() : Base() {}
  Halfedge_around_vertex_circulator_adaptor(const Arg& he)
    : Base(he) {}

 private:
  Halfedge_around_vertex_circulator_adaptor(const Base& b) : Base(b) {}

  void increment() {
    this->cur_ = *this->cur_.next()->opposite();
  }

  void decrement() {
    this->cur_ = *this->cur_.opposite()->previous();
  }
};

//=========================================================================
//=========================================================================

template<class Halfedge>
class Ccb_halfedge_circulator_adaptor
  : public Circulator_from_halfedge_adaptor<Halfedge,
      Ccb_halfedge_circulator_adaptor<Halfedge> >
{
 private:
  typedef Halfedge                                    Arg;
  typedef Ccb_halfedge_circulator_adaptor<Arg>        Self;
  typedef Circulator_from_halfedge_adaptor<Arg,Self>  Base;

  friend class Circulator_from_halfedge_adaptor<Arg,Self>;

 public:
  Ccb_halfedge_circulator_adaptor() : Base() {}
  Ccb_halfedge_circulator_adaptor(const Arg& he) : Base(he) {}

 private:
  Ccb_halfedge_circulator_adaptor(const Base& b) : Base(b) {}

  void increment() {
    this->cur_ = *this->cur_.next();
  }

  void decrement() {
    this->cur_ = *this->cur_.previous();
  }
};

//=========================================================================
//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL


#endif // CGAL_VORONOI_DIAGRAM_2_CIRCULATOR_ADAPTORS_H
