// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mariette Yvinec

#ifndef CGAL_TRIANGULATION_DS_FACE_BASE_2_H
#define CGAL_TRIANGULATION_DS_FACE_BASE_2_H

#include <CGAL/license/TDS_2.h>


#include <CGAL/config.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Dummy_tds_2.h>

namespace CGAL {

template < typename TDS = void>
class Triangulation_ds_face_base_2
{
public:
  typedef TDS                          Triangulation_data_structure;
  typedef typename TDS::Vertex_handle  Vertex_handle;
  typedef typename TDS::Face_handle    Face_handle;

  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_ds_face_base_2<TDS2> Other; };

private:
  Vertex_handle V[3];
  Face_handle   N[3];

public:
  Triangulation_ds_face_base_2();
  Triangulation_ds_face_base_2(Vertex_handle v0,
                               Vertex_handle v1,
                               Vertex_handle v2);
  Triangulation_ds_face_base_2(Vertex_handle v0,
                               Vertex_handle v1,
                               Vertex_handle v2,
                               Face_handle n0,
                               Face_handle n1,
                               Face_handle n2);

  Vertex_handle vertex(int i) const;
  bool has_vertex(Vertex_handle v) const;
  bool has_vertex(Vertex_handle v, int& i) const ;
  int index(Vertex_handle v) const ;

  Face_handle neighbor(int i) const ;
  bool has_neighbor(Face_handle n) const;
  bool has_neighbor(Face_handle n, int& i) const;
  int index(Face_handle n) const;

  void set_vertex(int i, Vertex_handle v);
  void set_vertices();
  void set_vertices(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2);
  void set_neighbor(int i, Face_handle n) ;
  void set_neighbors();
  void set_neighbors(Face_handle n0, Face_handle n1, Face_handle n2);
  void reorient();
  void ccw_permute();
  void cw_permute();

  int dimension() const;
  //the following trivial is_valid to allow
  // the user of derived face base classes
  // to add their own purpose checking
  bool is_valid(bool /* verbose */ = false, int /* level */ = 0) const
  {return true;}

   // For use by Compact_container.
  void * for_compact_container() const {return N[0].for_compact_container(); }
  void for_compact_container(void* p) { N[0].for_compact_container(p);}


  static int ccw(int i) {return Triangulation_cw_ccw_2::ccw(i);}
  static int  cw(int i) {return Triangulation_cw_ccw_2::cw(i);}
};

template <class TDS>
Triangulation_ds_face_base_2<TDS> ::
Triangulation_ds_face_base_2()
{
  set_vertices();
  set_neighbors();
}

template <class TDS>
Triangulation_ds_face_base_2<TDS> ::
Triangulation_ds_face_base_2( Vertex_handle v0,
                              Vertex_handle v1,
                              Vertex_handle v2)
{
  set_vertices(v0, v1, v2);
  set_neighbors();
}

template <class TDS>
Triangulation_ds_face_base_2<TDS> ::
Triangulation_ds_face_base_2(Vertex_handle v0,
                             Vertex_handle v1,
                             Vertex_handle v2,
                             Face_handle n0,
                             Face_handle n1,
                             Face_handle n2)
{
  set_vertices(v0, v1, v2);
  set_neighbors(n0, n1, n2);
}


template <class TDS>
inline
typename Triangulation_ds_face_base_2<TDS>::Vertex_handle
Triangulation_ds_face_base_2<TDS>::
vertex(int i) const
{
  CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
  return V[i];
}


template <class TDS>
inline bool
Triangulation_ds_face_base_2<TDS> ::
has_vertex(Vertex_handle v) const
{
  return (V[0] == v) || (V[1] == v) || (V[2]== v);
}

template <class TDS>
inline bool
Triangulation_ds_face_base_2<TDS> ::
has_vertex(Vertex_handle v, int& i) const
{
  if (v == V[0]) {
    i = 0;
    return true;
  }
  if (v == V[1]) {
    i = 1;
    return true;
  }
  if (v == V[2]) {
    i = 2;
    return true;
  }
  return false;
}

template <class TDS>
inline int
Triangulation_ds_face_base_2<TDS> ::
index(Vertex_handle v) const
{
  if (v == V[0]) return 0;
  if (v == V[1]) return 1;
  CGAL_triangulation_assertion( v == V[2] );
  return 2;
}

template <class TDS>
inline
typename Triangulation_ds_face_base_2<TDS>::Face_handle
Triangulation_ds_face_base_2<TDS>::
neighbor(int i) const
{
  CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
  return N[i];
}

template <class TDS>
inline bool
Triangulation_ds_face_base_2<TDS> ::
has_neighbor(Face_handle n) const
{
  return (N[0] == n) || (N[1] == n) || (N[2] == n);
}


template <class TDS>
inline bool
Triangulation_ds_face_base_2<TDS> ::
has_neighbor(Face_handle n, int& i) const
{
  if(n == N[0]){
    i = 0;
    return true;
  }
  if(n == N[1]){
    i = 1;
    return true;
  }
  if(n == N[2]){
    i = 2;
    return true;
  }
  return false;
}



template <class TDS>
inline int
Triangulation_ds_face_base_2<TDS> ::
index(Face_handle n) const
{
  if (n == N[0]) return 0;
  if (n == N[1]) return 1;
  CGAL_triangulation_assertion( n == N[2] );
  return 2;
}

template <class TDS>
inline void
Triangulation_ds_face_base_2<TDS> ::
set_vertex(int i, Vertex_handle v)
{
  CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
  V[i] = v;
}

template <class TDS>
inline void
Triangulation_ds_face_base_2<TDS> ::
set_neighbor(int i, Face_handle n)
{
  CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition( this != &*n );
  N[i] = n;
}

template <class TDS>
inline void
Triangulation_ds_face_base_2<TDS> ::
set_vertices()
{
  V[0] = V[1] = V[2] = Vertex_handle();
}

template <class TDS>
inline void
Triangulation_ds_face_base_2<TDS> ::
set_vertices(Vertex_handle v0,  Vertex_handle v1, Vertex_handle v2)
{
  V[0] = v0;
  V[1] = v1;
  V[2] = v2;
}

template <class TDS>
inline void
Triangulation_ds_face_base_2<TDS> ::
set_neighbors()
{
  N[0] = N[1] = N[2] = Face_handle();
}

template <class TDS>
inline void
Triangulation_ds_face_base_2<TDS> ::
set_neighbors(Face_handle n0,Face_handle n1, Face_handle n2)
{
  CGAL_triangulation_precondition( this != n0.operator->() );
  CGAL_triangulation_precondition( this != n1.operator->() );
  CGAL_triangulation_precondition( this != n2.operator->() );
  N[0] = n0;
  N[1] = n1;
  N[2] = n2;
}

template <class TDS>
void
Triangulation_ds_face_base_2<TDS> ::
reorient()
{
  //exchange the vertices 0 and 1
  set_vertices (V[1],V[0],V[2]);
  set_neighbors(N[1],N[0],N[2]);
}

template <class TDS>
inline void
Triangulation_ds_face_base_2<TDS> ::
ccw_permute()
{
  set_vertices (V[2],V[0],V[1]);
  set_neighbors(N[2],N[0],N[1]);
}


template <class TDS>
inline void
Triangulation_ds_face_base_2<TDS> ::
cw_permute()
{
  set_vertices (V[1],V[2],V[0]);
  set_neighbors(N[1],N[2],N[0]);
}


template < class TDS>
inline  int
Triangulation_ds_face_base_2<TDS> ::
dimension() const
{
  if (V[2] != Vertex_handle()) {return 2;}
  else return( V[1] != Vertex_handle() ? 1 : 0);
}

template < class TDS >
inline
std::istream&
operator>>(std::istream &is, Triangulation_ds_face_base_2<TDS> &)
  // non combinatorial information. Default = nothing
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os, const Triangulation_ds_face_base_2<TDS> &)
  // non combinatorial information. Default = nothing
{
  return os;
}

// Specialisation for void.
template <>
class Triangulation_ds_face_base_2<void>
{
public:
  typedef Dummy_tds_2                       Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Face_handle     Face_handle;
  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_ds_face_base_2<TDS2> Other; };
};



} //namespace CGAL

#endif //CGAL_DS_TRIANGULATION_FACE_BASE_2_H
