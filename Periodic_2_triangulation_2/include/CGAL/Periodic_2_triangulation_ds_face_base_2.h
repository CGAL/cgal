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
// $URL$
// $Id$
// 
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_DS_FACE_BASE_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_DS_FACE_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Dummy_tds_2.h>

namespace CGAL { 

template < typename TDS = void>
class Periodic_2_triangulation_ds_face_base_2 
{
public:
  typedef TDS                          Triangulation_data_structure;
  typedef typename TDS::Vertex_handle  Vertex_handle;
  typedef typename TDS::Face_handle    Face_handle;
#ifdef CGAL_TDS2_DATA
  typedef typename TDS::Face_data      TDS_data;
#endif
  
  template <typename TDS2>
  struct Rebind_TDS { typedef Periodic_2_triangulation_ds_face_base_2<TDS2> Other; }; 
  
private:
  Vertex_handle V[3];
  Face_handle   N[3];
#ifdef CGAL_TDS2_DATA
  TDS_data      _tds_data;
#endif
  
  // 2 respective bits are the _offset in x and y
  // right to left: 
  // bit[0]-bit[1]: vertex(0),
  // bit[2]-bit[3]: vertex(1) and 
  // bit[4]-bit[5]: vertex(2)
  // Thus the underlying data type needs to have at least 6 bit,
  // which is true for an unsigned char.
  // bit[6]: Used to convert 9 sheeted covering to a 1 sheeted covering
  unsigned char _off;
  
public:
  Periodic_2_triangulation_ds_face_base_2();
  Periodic_2_triangulation_ds_face_base_2(Vertex_handle v0, 
                                          Vertex_handle v1, 
                                          Vertex_handle v2);
  Periodic_2_triangulation_ds_face_base_2(Vertex_handle v0, 
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
  void * & for_compact_container()     { return N[0].for_compact_container();}
  
  
  static int ccw(int i) {return Triangulation_cw_ccw_2::ccw(i);}
  static int  cw(int i) {return Triangulation_cw_ccw_2::cw(i);}
  
#ifdef CGAL_TDS2_DATA
  // TDS internal data access functions.
  TDS_data& tds_data()       { return _tds_data; }
  const TDS_data& tds_data() const { return _tds_data; }
#endif
  
  /// Periodic functions
  int offset(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i < 3 );
    return ((_off>>2*i)&3);
  }
  void set_offsets(unsigned int o0, unsigned int o1, unsigned int o2) {
    // 192=11000000
    _off = _off | 192;
    unsigned int off0[2] = {(o0>>1)&1, (o0&1)};
    unsigned int off1[2] = {(o1>>1)&1, (o1&1)};
    unsigned int off2[2] = {(o2>>1)&1, (o2&1)};
    for (int i=0; i<2; i++) {
      unsigned int _off0 = ( _off    &3);
      unsigned int _off1 = ((_off>>2)&3);
      unsigned int _off2 = ((_off>>4)&3);
      
      _off0 = ( (_off0<<1) + off0[i]);
      _off1 = ( (_off1<<1) + off1[i]);
      _off2 = ( (_off2<<1) + off2[i]);
      
      // 252=11111100
      // 243=11110011
      // 207=11001111
      _off = ((_off&252) | (_off0   ));
      _off = ((_off&243) | (_off1<<2));
      _off = ((_off&207) | (_off2<<4));
    }
  }
  
  void set_additional_flag(unsigned char b) {
    CGAL_assertion(b < 4);
    // 63=00111111
    _off = ((_off & 63) | (b<<6));
  }
  unsigned char get_additional_flag() {
    return (_off>>6);
  }

};

template <class TDS>
Periodic_2_triangulation_ds_face_base_2<TDS> ::
Periodic_2_triangulation_ds_face_base_2()
: _off(0)
{
  set_vertices();
  set_neighbors();
}

template <class TDS>
Periodic_2_triangulation_ds_face_base_2<TDS> ::
Periodic_2_triangulation_ds_face_base_2( Vertex_handle v0, 
                                        Vertex_handle v1, 
                                        Vertex_handle v2)
: _off(0)
{
  set_vertices(v0, v1, v2);
  set_neighbors();
}

template <class TDS>
Periodic_2_triangulation_ds_face_base_2<TDS> ::
Periodic_2_triangulation_ds_face_base_2(Vertex_handle v0, 
                                        Vertex_handle v1, 
                                        Vertex_handle v2,
                                        Face_handle n0, 
                                        Face_handle n1, 
                                        Face_handle n2)
: _off(0)
{
  set_vertices(v0, v1, v2);
  set_neighbors(n0, n1, n2);
}


template <class TDS>
inline 
typename Periodic_2_triangulation_ds_face_base_2<TDS>::Vertex_handle
Periodic_2_triangulation_ds_face_base_2<TDS>::
vertex(int i) const
{
  CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
  return V[i];
} 


template <class TDS>
inline bool
Periodic_2_triangulation_ds_face_base_2<TDS> ::
has_vertex(Vertex_handle v) const
{
  return (V[0] == v) || (V[1] == v) || (V[2]== v);
}

template <class TDS>
inline bool
Periodic_2_triangulation_ds_face_base_2<TDS> ::    
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
Periodic_2_triangulation_ds_face_base_2<TDS> :: 
index(Vertex_handle v) const
{
  if (v == V[0]) return 0;
  if (v == V[1]) return 1;
  CGAL_triangulation_assertion( v == V[2] );
  return 2;
}

template <class TDS>    
inline 
typename Periodic_2_triangulation_ds_face_base_2<TDS>::Face_handle 
Periodic_2_triangulation_ds_face_base_2<TDS>::
neighbor(int i) const
{
  CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
  return N[i];
}

template <class TDS>      
inline bool 
Periodic_2_triangulation_ds_face_base_2<TDS> ::
has_neighbor(Face_handle n) const
{
  return (N[0] == n) || (N[1] == n) || (N[2] == n);
}


template <class TDS>      
inline bool 
Periodic_2_triangulation_ds_face_base_2<TDS> ::
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
Periodic_2_triangulation_ds_face_base_2<TDS> ::
index(Face_handle n) const
{
  if (n == N[0]) return 0;
  if (n == N[1]) return 1;
  CGAL_triangulation_assertion( n == N[2] );
  return 2;
}

template <class TDS>      
inline void
Periodic_2_triangulation_ds_face_base_2<TDS> :: 
set_vertex(int i, Vertex_handle v)
{
  CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
  V[i] = v;
}

template <class TDS>      
inline void
Periodic_2_triangulation_ds_face_base_2<TDS> ::     
set_neighbor(int i, Face_handle n)
{
  CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition( this != &*n );
  N[i] = n;
}

template <class TDS>      
inline void
Periodic_2_triangulation_ds_face_base_2<TDS> :: 
set_vertices()
{
  V[0] = V[1] = V[2] = Vertex_handle();
}

template <class TDS>      
inline void
Periodic_2_triangulation_ds_face_base_2<TDS> ::     
set_vertices(Vertex_handle v0,  Vertex_handle v1, Vertex_handle v2)
{
  V[0] = v0;
  V[1] = v1;
  V[2] = v2;
}

template <class TDS>      
inline void
Periodic_2_triangulation_ds_face_base_2<TDS> :: 
set_neighbors()
{
  N[0] = N[1] = N[2] = Face_handle();
}

template <class TDS>      
inline void
Periodic_2_triangulation_ds_face_base_2<TDS> :: 
set_neighbors(Face_handle n0,Face_handle n1, Face_handle n2)
{
  CGAL_triangulation_precondition( this != &*n0 );
  CGAL_triangulation_precondition( this != &*n1 );
  CGAL_triangulation_precondition( this != &*n2 );
  N[0] = n0;
  N[1] = n1;
  N[2] = n2;
}

template <class TDS>
void
Periodic_2_triangulation_ds_face_base_2<TDS> :: 
reorient()
{
  //exchange the vertices 0 and 1
  set_vertices (V[1],V[0],V[2]);
  set_neighbors(N[1],N[0],N[2]);
}

template <class TDS>
inline void 
Periodic_2_triangulation_ds_face_base_2<TDS> ::
ccw_permute()
{
  set_vertices (V[2],V[0],V[1]);
  set_neighbors(N[2],N[0],N[1]);
}


template <class TDS>
inline void 
Periodic_2_triangulation_ds_face_base_2<TDS> ::
cw_permute()
{
  set_vertices (V[1],V[2],V[0]);
  set_neighbors(N[1],N[2],N[0]);
}


template < class TDS>
inline  int 
Periodic_2_triangulation_ds_face_base_2<TDS> ::
dimension() const
{
  if (V[2] != Vertex_handle()) {return 2;}
  else return( V[1] != Vertex_handle() ? 1 : 0);
}

template < class TDS >
inline
std::istream&
operator>>(std::istream &is, Periodic_2_triangulation_ds_face_base_2<TDS> &)
// non combinatorial information. Default = nothing
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os, const Periodic_2_triangulation_ds_face_base_2<TDS> &)
// non combinatorial information. Default = nothing
{
  return os;
}

// Specialisation for void.
template <>
class Periodic_2_triangulation_ds_face_base_2<void>
{
public:
  typedef Dummy_tds_2                       Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Face_handle     Face_handle;
  template <typename TDS2>
  struct Rebind_TDS { typedef Periodic_2_triangulation_ds_face_base_2<TDS2> Other; };
};



} //namespace CGAL 

#endif //CGAL_PERIODIC_2_TRIANGULATION_DS_FACE_BASE_2_H
