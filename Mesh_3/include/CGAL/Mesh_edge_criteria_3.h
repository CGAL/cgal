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


/*!
\ingroup PkgMesh3MeshClasses

The function object class `Mesh_edge_criteria_3` is a model of `MeshEdgeCriteria_3`. It
provides a bound for the size criterion.

\cgalModels{MeshEdgeCriteria_3}

\sa `MeshCriteria_3`
\sa `CGAL::Mesh_criteria_3<Tr>`
\sa `MeshDomainField_3`
*/
template < typename Tr >
class Mesh_edge_criteria_3
{
private:
  typedef Mesh_edge_criteria_3 Self;
  typedef typename Tr::Geom_traits GT;

public:

  /// \name Types
  /// @{
  /*!
  Numerical type.
  */
  typedef typename Tr::Geom_traits::FT  FT;
  typedef typename Tr::Vertex::Index    Index;
  typedef typename Tr::Bare_point       Point_3;

  /// @}


  /// \name Creation
  /// @{
  /*!
  * returns an object to serve as criteria for edges.
  *
  * \param length_bound is an upper bound
  * for the length of the edges which are used to discretize the curves.
  * \param min_length_bound is a desired lower bound
  * for the length of the edges which are used to discretize the curves.
  * Only edges that are longer than this bound will be refined. Using
  * this lower bound can be handy on some domains, but using it may
  * break all the surface topology guarantees of the meshing algorithm.
  * It is not guaranteed to be exactly respected in the output mesh.
  *
  * \note If one parameter is set to 0, then its corresponding criterion is ignored.
  */
  Mesh_edge_criteria_3(const FT& length_bound,
                       const FT& min_length_bound = 0)

    : p_size_(new Mesh_3::internal::Sizing_field_container<
                Mesh_constant_domain_field_3<GT,Index> ,
                FT,
                Point_3,
                Index>(length_bound))
    , min_length_bound_(min_length_bound)
  {}

  // Nb: SFINAE to avoid wrong matches with built-in numerical types
  // as int.

  /*!
  * returns an object to serve as criteria for edges.
  * The behavior and semantic of the argument are the same
  * as above, except that the `length_bound`
  * parameter is a functional instead of a constant.
  *
  * @tparam SizingField a model of `MeshDomainField_3`
  */
  template < typename SizingField >
  Mesh_edge_criteria_3
  (

   const SizingField& length_bound,
   const FT& min_length_bound = 0
#ifndef DOXYGEN_RUNNING
    , std::enable_if_t<Mesh_3::Is_mesh_domain_field_3<Tr, SizingField>::value>* = 0
#endif
   )
   : min_length_bound_(min_length_bound)
  {
    p_size_ = new Mesh_3::internal::Sizing_field_container<SizingField,
                                                           FT,
                                                           Point_3,
                                                           Index>(length_bound);
  }

  /// @}

#ifndef DOXYGEN_RUNNING
  Mesh_edge_criteria_3(const Self& rhs)
    : p_size_(rhs.p_size_->clone())
    , min_length_bound_(rhs.min_length_bound_)
  {}

  /// Destructor
  ~Mesh_edge_criteria_3()
  {
    delete p_size_;
  }

  /// Returns size of tuple (p,dim,index)
  FT sizing_field(const Point_3& p, const int dim, const Index& index) const
  { return (*p_size_)(p,dim,index); }

public:
  const FT& min_length_bound() const
  {
    return min_length_bound_;
  }
#endif

private:
  typedef Mesh_3::internal::Sizing_field_interface<FT,Point_3,Index>
    Sizing_field_interface;

  // A pointer to Sizing_field_interface to handle dynamic wrapping of
  // real Sizing_field type
  Sizing_field_interface* p_size_;
  const FT min_length_bound_;
};

} // end namespace CGAL

#endif // CGAL_MESH_EDGE_CRITERIA_3_H
