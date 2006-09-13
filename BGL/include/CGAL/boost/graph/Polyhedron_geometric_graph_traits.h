// Copyright (c) 2006 Geometry Factory (France).
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Polyhedron_BGL.h $
// $Id: Polyhedron_BGL.h 32048 2006-06-23 13:59:36Z lsaboret $
// 
//
// Author(s): Andreas Fabri <andreas.fabri@geometryfactory.com>, Fernando Cacciola <fernando.cacciola@gmail.com>

#ifndef CGAL_BOOST_GRAPH_POLYHEDRON_GRAPH_TRAITS_H
#define CGAL_BOOST_GRAPH_POLYHEDRON_GRAPH_TRAITS_H

#include <CGAL/boost/graph/Geometric_graph_traits.h>

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
#  define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS
#else
#  define CGAL_HDS_PARAM_ class HDS
#endif


CGAL_BEGIN_NAMESPACE

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class geometric_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const > 
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;
  
public:
  
  typedef typename Polyhedron::Point_3 Point_3 ;
};

        
CGAL_END_NAMESPACE

#undef CGAL_HDS_

#endif // CGAL_BOOST_GRAPH_POLYHEDRON_GRAPH_TRAITS_H
