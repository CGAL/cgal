// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_EDGE_CRITERIA_3_H
#define CGAL_MESH_EDGE_CRITERIA_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_constant_domain_field_3.h>
#include <CGAL/Mesh_3/Is_mesh_domain_field_3.h>
#include <type_traits>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

  // Those two classes are designed to handle dynamic initialization of
  // Sizing_field type (using named parameters of make_mesh_3 for example)
  template < typename FT_, typename Point_, typename Index_ >
  class Sizing_field_interface
  {
  public:
    typedef FT_     FT;
    typedef Point_  Point_3;
    typedef Index_  Index;

    virtual ~Sizing_field_interface() {}

    virtual FT operator()(const Point_3& p,
                          const int dim,
                          const Index& index) const = 0;

    virtual Sizing_field_interface* clone() const = 0;
  };

  template < typename Sizing_field,
             typename FT,
             typename Point_3,
             typename Index>
  struct Sizing_field_container
    : public Sizing_field_interface < FT,
                                      Point_3,
                                      Index >
  {
    typedef Sizing_field_interface < FT,
                                     Point_3,
                                     Index > Base;

    typedef Sizing_field_container<Sizing_field, FT, Point_3, Index> Self;

  public:
    Sizing_field_container(const Sizing_field& s) : s_(s) {}
    virtual ~Sizing_field_container() {}

    virtual FT operator()(const Point_3& p,
                          const int dim,
                          const Index& index) const
    {
      return s_(p,dim,index);
    }

    virtual Base* clone() const
    {
      return new Self(*this);
    }

  private:
    Sizing_field s_;
  };

} // end namespace internal
} // end namespace Mesh_3

template < typename Tr >
class Mesh_edge_criteria_3
{
  typedef Mesh_edge_criteria_3 Self;

public:
  typedef typename Tr::Vertex::Index  Index;
  typedef typename Tr::Geom_traits    Gt;
  typedef typename Gt::FT             FT;
  typedef typename Tr::Bare_point     Point_3;

  /// Constructors
  Mesh_edge_criteria_3(const FT& value)
    : p_size_(new Mesh_3::internal::Sizing_field_container<
                Mesh_constant_domain_field_3<Gt,Index> ,
                FT,
                Point_3,
                Index>(value))
  {}

  // Nb: SFINAE to avoid wrong matches with built-in numerical types
  // as int.
  template < typename Sizing_field >
  Mesh_edge_criteria_3
  (
   const Sizing_field& size,
   typename std::enable_if<Mesh_3::Is_mesh_domain_field_3<Tr, Sizing_field>::value>::type* = 0
   )
  {
    p_size_ = new Mesh_3::internal::Sizing_field_container<Sizing_field,
                                                           FT,
                                                           Point_3,
                                                           Index>(size);
  }

  Mesh_edge_criteria_3(const Self& rhs)
    : p_size_(rhs.p_size_->clone()) {}

  /// Destructor
  ~Mesh_edge_criteria_3()
  {
    delete p_size_;
  }

  /// Returns size of tuple (p,dim,index)
  FT sizing_field(const Point_3& p, const int dim, const Index& index) const
  { return (*p_size_)(p,dim,index); }

private:
  typedef Mesh_3::internal::Sizing_field_interface<FT,Point_3,Index>
    Sizing_field_interface;

  // A pointer to Sizing_field_interface to handle dynamic wrapping of
  // real Sizing_field type
  Sizing_field_interface* p_size_;
};

} // end namespace CGAL

#endif // CGAL_MESH_EDGE_CRITERIA_3_H
