// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.4-branch/Triangulation_2/include/CGAL/Regular_Triangulation_face_base_sphere_2.h $
// $Id: Regular_Triangulation_face_base_sphere_2.h 41117 2007-12-07 13:14:48Z yvinec $
// 
//
// Author(s)     : Frederic Fichel, Mariette Yvinec

#ifndef CGAL_DELAUNAY_TRIANGULATION_FACE_BASE_SPHERE_2_H
#define CGAL_DELAUNAY_TRIANGULATION_FACE_BASE_SPHERE_2_H

#include <list>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_face_base_sphere_2.h>

namespace CGAL {


template <class Gt, class Fb = Triangulation_face_base_sphere_2<Gt> >
class Delaunay_triangulation_face_base_sphere_2
  :  public Fb
{
  typedef Fb                                            Fbase;
  typedef typename Fbase::Triangulation_data_structure  TDS;
public:
  typedef Gt                                   Geom_traits;
  typedef TDS                                  Triangulation_data_structure;
  typedef typename TDS::Vertex_handle          Vertex_handle;
  typedef typename TDS::Face_handle            Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other   Fb2;
    typedef Delaunay_triangulation_face_base_sphere_2<Gt,Fb2>             Other;
  };

  typedef std::list<Vertex_handle>             Vertex_list;

protected:
  Vertex_list vlist;
  unsigned char _in_conflict_flag;

public:
  void set_in_conflict_flag(unsigned char f) { _in_conflict_flag = f; }
  unsigned char get_in_conflict_flag() const { return _in_conflict_flag; }

  Delaunay_triangulation_face_base_sphere_2()
  : Fbase(),  vlist()
  {
    set_in_conflict_flag(0);
  }

 Delaunay_triangulation_face_base_sphere_2(Vertex_handle v0, 
				    Vertex_handle v1, 
				    Vertex_handle v2)
  : Fbase(v0,v1,v2), vlist()
  {
    set_in_conflict_flag(0);
  }

  Delaunay_triangulation_face_base_sphere_2(Vertex_handle v0, 
				    Vertex_handle v1, 
				    Vertex_handle v2,
				    Face_handle n0, 
				    Face_handle n1, 
				    Face_handle n2)
  : Fbase(v0,v1,v2,n0,n1,n2),  vlist()
  { 
       set_in_conflict_flag(0);
  }

  ~Delaunay_triangulation_face_base_sphere_2()
  { 
    vlist.clear();
  }


  Vertex_list& vertex_list()
  {
    return vlist;
  }

	
	bool is_valid(bool  verbose  = false, int  level  = 0) const
	{return true;}
	


};

} //namespace CGAL 

#endif // CGAL_REGULAR_Triangulation_face_base_sphere_2_H

