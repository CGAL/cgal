// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_vertex_base_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class GT > class Triangulation_vertex_base_3;
template < class GT >
std::istream& operator >> (std::istream&, Triangulation_vertex_base_3<GT>&);
  
template < class GT >
class Triangulation_vertex_base_3
{
  friend std::istream& operator >> CGAL_NULL_TMPL_ARGS
  (std::istream&, Triangulation_vertex_base_3<GT>&);

public:
  typedef typename GT::Point_3 Point;
  
  // CONSTRUCTORS
  
  inline 
  Triangulation_vertex_base_3()
    : _p(), _c(NULL)
  {}
  
  inline 
  Triangulation_vertex_base_3(const Point & p)
    :  _p(p), _c(NULL)
  {}
    
  inline 
  Triangulation_vertex_base_3(const Point & p, void* c)
    :  _p(p), _c(c)
  {}

  inline 
  Triangulation_vertex_base_3(void* c)
    :  _p(), _c(c)
  {}

  // ACCES 

  inline 
  const Point & point() const
  { return _p; }
    
  inline 
  void* cell() const
  { return _c; }

  // SETTING

  inline 
  void set_point(const Point & p)
  { _p = p; }
    
  inline 
  void set_cell(void* c)
  { _c = c; }

  // CHECKING

  // the following trivial is_valid allows
  // the user of derived cell base classes 
  // to add their own purpose checking
  bool is_valid(bool, int ) const
  { 
    //return true; 
    return( cell() != NULL );
  }


private:
  Point _p;
  void * _c;
  
};

template < class GT >
std::istream& operator >>
(std::istream& is, Triangulation_vertex_base_3<GT> & v)
  // non combinatorial information. Default = point
{
  typedef typename GT::Point_3 Point;
  //  Point p;
  //  is >> p;
  is >> v._p;
  //  v.set_point(p);
  return is;
}
template < class GT >
std::ostream& operator<<
(std::ostream& os, const Triangulation_vertex_base_3<GT> & v)
  // non combinatorial information. Default = point
{
  os << v.point();
  //  if(is_ascii(os)){
  //    os << endl;
  //  }
  return os;
}

template < class GT >
class Triangulation_vertex_base_pointer_3;
template < class GT >
std::istream& operator >> 
(std::istream&, Triangulation_vertex_base_pointer_3<GT>&);

template < class GT >
class Triangulation_vertex_base_pointer_3
{
  friend std::istream& operator >> CGAL_NULL_TMPL_ARGS
  (std::istream&, Triangulation_vertex_base_pointer_3<GT>&);

public:
  typedef typename GT::Point_3 Point;
  //typedef Point Point_3;
  
  // CONSTRUCTORS
  
  inline 
  Triangulation_vertex_base_pointer_3()
    : _p(NULL), _c(NULL)
  {}
  
  inline 
  Triangulation_vertex_base_pointer_3(const Point & p)
    :  _p(&p), _c(NULL)
  {
    //    const Point * q =  &p;
    //    _p = q;
  }
    
  inline 
  Triangulation_vertex_base_pointer_3(const Point & p, void* c)
    :  _p(p), _c(c)
  {}

  inline 
  Triangulation_vertex_base_pointer_3(void* c)
    :  _p(NULL), _c(c)
  {}

  // ACCES 

  inline 
  const Point & point() const
  { return *_p; }
    
  inline 
  void* cell() const
  { return _c; }

  // SETTING

//   inline 
//   void set_point(Point * p)
//   {     
//     const Point * q =  &p;
//     p = q; 
//   }
    
  inline 
  void set_point(const Point & p)
  {     
    //    const Point * q =  &p;
    _p = &p; 
  }
    
  inline 
  void set_cell(void* c)
  { _c = c; }

  // CHECKING

  // the following trivial is_valid allows
  // the user of derived cell base classes 
  // to add their own purpose checking
  bool is_valid(bool, int ) const
  { 
    //return true; 
    return( cell() != NULL );
  }


private:
  const Point * _p;
  void * _c;
  
};

template < class GT >
std::istream& operator >>
(std::istream& is, Triangulation_vertex_base_pointer_3<GT> & v)
  // non combinatorial information. Default = point
{ 
  typedef typename GT::Point_3 Point;
  Point q;
  is >> q;
  const Point * p = new Point(q);
  v._p = p;
  return is;
}
template < class GT >
std::ostream& operator<<
(std::ostream& os, const Triangulation_vertex_base_pointer_3<GT> & v)
  // non combinatorial information. Default = point
{
  os << v.point();
  //  if(is_ascii(os)){
  //    os << endl;
  //  }
  return os;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_VERTEX_BASE_3_H
