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
// file          : include/CGAL/Cartesian/Point_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================
 
#ifndef CGAL_CARTESIAN_POINT_D_H
#define CGAL_CARTESIAN_POINT_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/d_tuple.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
class DACd;

template < class R_ >
class PointCd CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
  typedef const FT*                             const_iterator;
  typedef FT*                                   iterator;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef PointCd<R,Cartesian_tag>              Self;
  typedef typename R::Vector_d                  Vector_d;
  typedef typename R::Direction_d               Direction_d;
  typedef typename R::Aff_transformation_d      Aff_transformation_d;
#else
  typedef PointCd<R>                            Self;
  typedef typename R::Vector_d_base             Vector_d;
  typedef typename R::Direction_d_base          Direction_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  friend class DACd<FT>;

  friend CGAL_FRIEND_INLINE
  Self operator- CGAL_NULL_TMPL_ARGS (const Self &p, const Vector_d &v);

  friend CGAL_FRIEND_INLINE
  Self operator+ CGAL_NULL_TMPL_ARGS (const Self &p, const Vector_d &v);

  PointCd (int dim= 0);
  PointCd (const Self &p);
  PointCd (const Vector_d &);
  PointCd (int dim, const Origin&);
  ~PointCd();

  template <class InputIterator>
  PointCd (int dim, InputIterator first, InputIterator last)
  {
    CGAL_kernel_precondition( dim >= 0);
    CGAL_kernel_precondition( (last-first == dim) || (last-first == dim+1) );
    PTR = new _d_tuple<FT>(dim);
    std::copy_n(first,dim,begin());
    if (last-first == dim+1)
      std::transform(begin(),end(),begin(),
                     bind2nd(divides<FT>(),*(first+dim)));
  }

  Self&          operator=(const Self &p);

  bool           operator==(const Self &p) const;
  bool           operator!=(const Self &p) const;
  bool           operator==(const Origin&) const;
  bool           operator!=(const Origin&) const;

  long           id() const;

  FT             homogeneous (int i) const;
  FT             cartesian (int i) const;
  FT             operator[] (int i) const;

  Self           transform(const Aff_transformation_d &t) const;

  const int      dimension() const { return ptr()->d; }
  const_iterator begin()     const { return ptr()->e; }
  const_iterator end()       const { return ptr()->e + dimension(); }

// protected:
  iterator       begin()           { return ptr()->e; }
  iterator       end()             { return ptr()->e + dimension(); }

private:
  const _d_tuple<FT>* ptr()  const { return (const _d_tuple<FT>*)PTR; }
  _d_tuple<FT>*       ptr()        { return (_d_tuple<FT>*)PTR; }
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Point_d.C>
#endif 

#endif // CGAL_CARTESIAN_POINT_D_H
