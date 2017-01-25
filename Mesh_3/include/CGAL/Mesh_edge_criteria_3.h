// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_EDGE_CRITERIA_3_H
#define CGAL_MESH_EDGE_CRITERIA_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_constant_domain_field_3.h>
#include <boost/type_traits.hpp>

namespace CGAL {

namespace internal {
namespace Mesh_3 {

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
  
  template < typename Sizing_field >
  struct Sizing_field_container
    : public Sizing_field_interface < typename Sizing_field::FT,
                                      typename Sizing_field::Point_3,
                                      typename Sizing_field::Index >
  {
    typedef Sizing_field_interface <
              typename Sizing_field::FT,
              typename Sizing_field::Point_3,
              typename Sizing_field::Index > Base;
    
    typedef Sizing_field_container<Sizing_field> Self;
    
  public:
    typedef typename Base::FT       FT;
    typedef typename Base::Point_3  Point_3;
    typedef typename Base::Index    Index;
    
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
  
}} // end namespace internal::Mesh_3
  
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
    : p_size_(new internal::Mesh_3::Sizing_field_container<
                Mesh_constant_domain_field_3<Gt,Index> >(value))
  {}
  
  // Nb: SFINAE (dummy) to avoid wrong matches with built-in numerical types
  // as int.
  template < typename Sizing_field >
  Mesh_edge_criteria_3(const Sizing_field& size,
                       typename Sizing_field::FT /*dummy*/ = 0 )
  {
    CGAL_static_assertion((boost::is_same<typename Sizing_field::FT,
                                          FT>::value));
    CGAL_static_assertion((boost::is_same<typename Sizing_field::Point_3,
                                          Point_3>::value));
    CGAL_static_assertion((boost::is_same<typename Sizing_field::Index,
                                          Index>::value));
                          
    p_size_ = new internal::Mesh_3::Sizing_field_container<Sizing_field>(size);
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
  typedef internal::Mesh_3::Sizing_field_interface<FT,Point_3,Index>
    Sizing_field_interface;
  
  // A pointer to Sizing_field_interface to handle dynamic wrapping of
  // real Sizing_field type
  Sizing_field_interface* p_size_;
};
  
} // end namespace CGAL

#endif // CGAL_MESH_EDGE_CRITERIA_3_H
