// ============================================================================
//
// Copyright (c) 1999,2000,2001,2002,2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : Triangulation/include/CGAL/Triangulation_ds_vertex_base_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_TRIANGULATION_DS_VERTEX_BASE_2_H
#define CGAL_TRIANGULATION_DS_VERTEX_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Dummy_tds_2.h>

CGAL_BEGIN_NAMESPACE

template < class TDS = void >
class Triangulation_ds_vertex_base_2 
{
  typedef typename TDS::Face_handle    Face_handle;
public:
  typedef TDS                          Triangulation_data_structure;

   // Borland seems to require it.
  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_ds_vertex_base_2<TDS2> Other; };

  Triangulation_ds_vertex_base_2 ()    : _f(NULL)  {}
  Triangulation_ds_vertex_base_2(Face_handle f)    :  _f(f)    {}

  Face_handle face() const { return _f;}
  void set_face(Face_handle f) { _f = f ;}

  //the following trivial is_valid to allow
  // the user of derived face base classes 
  // to add their own purpose checking
  bool is_valid(bool /*verbose*/=false, int /*level*/= 0) const
    {return face() != NULL;}

private:
  Face_handle _f;
  // FIXME : temporary cruft needed for DS_Container (sizeof() too small) :
  int i;
};

// Specialization for void.
template <>
class Triangulation_ds_vertex_base_2<void>
{
public:
  typedef Dummy_tds_2  Triangulation_data_structure;
  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_ds_vertex_base_2<TDS2> Other; };
};


template < class TDS >
inline
std::istream&
operator>>(std::istream &is, Triangulation_ds_vertex_base_2<TDS> &)
  // no combinatorial information.
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os, const Triangulation_ds_vertex_base_2<TDS> &)
  // no combinatorial information.
{
  return os;
}
CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_DS_VERTEX_BASE_2_H
