// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Direction_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DIRECTION_D_H
#define CGAL_CARTESIAN_DIRECTION_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/d_tuple.h>
#include <algorithm>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class DirectionCd CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle
{
public:
  typedef R_                                    R;
  typedef typename R::RT                        RT;
  typedef const RT*                             const_iterator;
  typedef RT*                                   iterator;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef DirectionCd<R CGAL_CTAG>              Self;
  typedef typename R::Vector_d                  Vector_d;
  typedef typename R::Aff_transformation_d      Aff_transformation_d;
#else
  typedef DirectionCd<R>                        Self;
  typedef typename R::Vector_d_base             Vector_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  DirectionCd(int d = 0);
  DirectionCd(const Self &d);
  DirectionCd(const Vector_d &v);
  ~DirectionCd();

  template < class InputIterator >
  DirectionCd(int dim, const InputIterator &first, const InputIterator &last)
  {
    CGAL_kernel_precondition( dim > 0);
    CGAL_kernel_precondition( last-first == dim );
    PTR = new _d_tuple<RT>(dim);
    std::copy_n(first,dim,begin());
  }

  Self&          operator=(const Self &d);

  bool           operator==(const Self &d) const;
  bool           operator!=(const Self &d) const;

  long           id() const;
  
  bool           is_degenerate() const;
  Vector_d       to_vector() const;
  Self           operator-() const;
  Self           transform(const Aff_transformation_d &t) const;

  RT             delta(int i) const;

  int            dimension() const { return ptr()->d; }
  const_iterator begin()     const { return ptr()->e; }
  const_iterator end()       const { return ptr()->e + dimension(); }

// protected:
  iterator       begin()           { return ptr()->e; }
  iterator       end()             { return ptr()->e + dimension(); }

private:
  const _d_tuple<RT>* ptr()  const { return (const _d_tuple<RT>*)PTR; }
  _d_tuple<RT>*       ptr()        { return (_d_tuple<RT>*)PTR; }
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Direction_d.C>
#endif 

#endif // CGAL_CARTESIAN_DIRECTION_D_H
