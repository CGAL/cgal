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

#ifndef CGAL_MESH_3_MESH_CONSTANT_DOMAIN_FIELD_3_H
#define CGAL_MESH_3_MESH_CONSTANT_DOMAIN_FIELD_3_H

#include <map>
#include <utility>

namespace CGAL {

template <typename Gt, typename Index_>
class Mesh_constant_domain_field_3
{
public:
  typedef typename Gt::FT         FT;
  typedef typename Gt::Point_3    Point_3;
  typedef Index_                  Index;
  
private:
  // Map to store field values
  typedef std::map<std::pair<int,Index>,FT> Values;
  
public:
  /// Constructor
  Mesh_constant_domain_field_3(const FT& d) : d_(d) {}

  /// Returns size
  FT operator()(const Point_3&, const int dim, const Index& index) const
  {
    typename Values::const_iterator it = values_.find(std::make_pair(dim,index));
    if ( it != values_.end() ) { return it->second; }
    
    return d_;
  }
  
  /// Sets size at any point of dimension \c dim and index \c index
  void set_size(const FT& size, const int dim, const Index& index)
  {
    values_.insert(std::make_pair(std::make_pair(dim,index),size));
  }

private:
  // default value
  FT d_;
  // Other values
  Values values_;
};

} // end namespace CGAL

#endif // CGAL_MESH_3_MESH_CONSTANT_DOMAIN_FIELD_3_H
