// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: https://scm.gforge.inria.fr/svn/cgal/branches/features/Mesh_3-experimental-GF/Mesh_3/include/CGAL/internal/Mesh_3/get_index.h $
// $Id: get_index.h 67573 2012-02-02 14:54:51Z lrineau $
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************

#ifndef CGAL_INTERNAL_MESH_3_GET_INDEX_3_H
#define CGAL_INTERNAL_MESH_3_GET_INDEX_3_H

#include <CGAL/license/Mesh_3.h>


#include <boost/type_traits/is_same.hpp>
#include <CGAL/Mesh_3/Has_features.h>
#include <CGAL/IO/io.h>

namespace CGAL {
namespace internal {
namespace Mesh_3 {
  
  
template <typename T, typename Boost_variant>
const T& get_index(const Boost_variant& x,
                   typename boost::disable_if<boost::is_same<T, Boost_variant> >::type * = 0)
{ return boost::get<T>(x); }

template <typename T>
const T& get_index(const T& x) { return x; }

template <typename Mesh_domain, 
          bool has_feature = Has_features<Mesh_domain>::value> 
struct Read_mesh_domain_index {
  // here we have has_feature==true

  typedef Mesh_domain MT; // was named "mesh traits" previously

  typename Mesh_domain::Index 
  operator()(int dimension, std::istream& is) const {
    switch(dimension) {
    case 0: 
      typename MT::Corner_index ci;
      if(is_ascii(is)) is >> ci;
      else CGAL::read(is, ci);
      return  ci;
      break;
    case 1: 
      typename MT::Curve_segment_index si;
      if(is_ascii(is)) is >> si;
      else CGAL::read(is, si);
      return  si;
      break;
    default:
      return Read_mesh_domain_index<Mesh_domain, false>()(dimension, is);
    }
  }
}; // end template partial specialization 
   // Read_mesh_domain_index<Mesh_domain, true>

template <typename Mesh_domain, 
          bool has_feature = Has_features<Mesh_domain>::value> 
struct Write_mesh_domain_index {
  // here we have has_feature==true

  typedef Mesh_domain MT; // was named "mesh traits" previously
  typedef typename MT::Corner_index Ci;
  typedef typename MT::Curve_segment_index  Si;

  void
  operator()(std::ostream& os, int dimension,
             const typename Mesh_domain::Index& index) const {
    switch(dimension) {
    case 0: {
      const Ci& ci = get_index<Ci>(index);
      if(is_ascii(os)) os << oformat(ci);
      else CGAL::write(os, ci);
    }
      break;
    case 1: {
      const Si& si = get_index<Si>(index);
      if(is_ascii(os)) os << oformat(si);
      else CGAL::write(os, si);
    }
      break;
    default:
      Write_mesh_domain_index<Mesh_domain, false>()(os, dimension, index);
    }
  }
}; // end template partial specialization 
   // Write_mesh_domain_index<Mesh_domain, true>

template <typename Mesh_domain>
struct Read_mesh_domain_index<Mesh_domain, false> {
  // here we have has_feature==false

  typedef Mesh_domain MT; // was named "mesh traits" previously

  typename Mesh_domain::Index 
  operator()(int dimension, std::istream& is) const {
    switch(dimension) {
    case 2: {
      typename MT::Surface_patch_index spi;
      if(is_ascii(is)) is >> iformat(spi);
      else CGAL::read(is, spi);
      return  spi;
    }
      break;
    default: {// 3
      typename MT::Subdomain_index di;
      if(is_ascii(is)) is >> iformat(di);
      else CGAL::read(is, di);
      return  di;
    }
      break;
    }
  }
}; // end template partial specialization 
   // Read_mesh_domain_index<Mesh_domain, false>

template <typename Mesh_domain>
struct Write_mesh_domain_index<Mesh_domain, false> {
  // here we have has_feature==false

  typedef Mesh_domain MT; // was named "mesh traits" previously
  typedef typename MT::Surface_patch_index Spi;
  typedef typename MT::Subdomain_index Di;

  void
  operator()(std::ostream& os, int dimension,
             const typename Mesh_domain::Index& index) const {
    switch(dimension) {
    case 2: {
      const Spi& spi = get_index<Spi>(index);
      if(is_ascii(os)) os << oformat(spi);
      else CGAL::write(os, spi);
    }
      break;
    default: {// 3
      const Di& di = get_index<Di>(index);
      if(is_ascii(os)) os << oformat(di);
      else CGAL::write(os, di);
    }
      break;
    }
  }
}; // end template partial specialization 
   // Write_mesh_domain_index<Mesh_domain, false>



}
}
}

#endif
