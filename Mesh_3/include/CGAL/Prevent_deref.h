  
// Copyright (c) 20012 GeometryFactory
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
// $URL: https://jtournois@geometryfactory.codebasehq.com/total/velocity.svn/Mesh_3-dec2011/CGAL/Prevent_deref.h $
// $Id: Prevent_deref.h 100 2012-05-29 14:42:38Z afabri $
//
//
// Author(s)     : Philipp Moeller
//
//******************************************************************************
// File Description : 
//******************************************************************************


#ifndef CGAL_PREVENT_DEREF_H
#define CGAL_PREVENT_DEREF_H

#include <boost/iterator_adaptors.hpp>

namespace CGAL {

template<typename I>
class Prevent_deref 
  : public boost::iterator_adaptor<
  Prevent_deref<I>
  , I // base
  , I // value
  >
{
public:
  typedef typename Prevent_deref::iterator_adaptor_ IA;
  typedef typename IA::reference reference;

  Prevent_deref() : Prevent_deref::iterator_adaptor_() {};
  Prevent_deref(const I& i) : Prevent_deref::iterator_adaptor_(i) {};
private:
  friend class boost::iterator_core_access;
  reference dereference() const { return const_cast<reference>(this->base_reference()); }
};

} // namespace CGAL

#endif // CGAL_PREVENT_DEREF_H
