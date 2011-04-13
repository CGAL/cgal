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
// file          : include/CGAL/Cartesian/Simplex_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_SIMPLEX_D_H
#define CGAL_CARTESIAN_SIMPLEX_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/d_tuple.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class SimplexCd CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef SimplexCd<R CGAL_CTAG>                Self;
  typedef typename R::Point_d                   Point_d;
  typedef typename R::Plane_d                   Plane_d;
  typedef typename R::Aff_transformation_d      Aff_transformation_d;
#else
  typedef SimplexCd<R>                          Self;
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Plane_d_base              Plane_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif
  typedef const Point_d*                        const_iterator;
  typedef Point_d*                              iterator;

  SimplexCd();
  SimplexCd(const Self &t);
  template < class PointIterator >
  SimplexCd(const PointIterator &first, const PointIterator &last)
    { 
      int dimension = last-first-1;
      // kernel assertions: the following will be pruned if they're off
      CGAL_kernel_precondition( dimension>0 );
      PointIterator p;
      for (p=first; p!=last; ++p)
	CGAL_kernel_precondition( p->dimension() == dimension );
      // This is the code
      new_rep(dimension);
      std::copy(first,last,begin());
    }
  ~SimplexCd();

  Self &operator=(const Self &t);

  Point_d       vertex(int i) const;
  Point_d       operator[](int i) const;

  bool          operator==(const Self &t) const;
  bool          operator!=(const Self &t) const;
  long          id() const;

  int           dimension() const { return ptr()->d-1; }
  int           number_of_vertices() const { return ptr()->d; }

  // Bbox_d     bbox() const;
  Self       transform(const Aff_transformation_d &t) const;

  Orientation   orientation() const;
  Oriented_side oriented_side(const Point_d &p) const;
  Bounded_side  bounded_side(const Point_d &p) const;

  bool          has_on_positive_side(const Point_d &p) const;
  bool          has_on_negative_side(const Point_d &p) const;
  bool          has_on_bounded_side(const Point_d &p) const;
  bool          has_on_unbounded_side(const Point_d &p) const;
  bool          has_on_boundary(const Point_d &p) const;
  bool          has_on(const Point_d &p) const;

  bool          is_degenerate() const;

  // coordinate access
  const_iterator begin()     const { return ptr()->e; }
  const_iterator end()       const { return ptr()->e+ptr()->d; }

// protected:
  iterator       begin()           { return ptr()->e; }
  iterator       end()             { return ptr()->e+ptr()->d; }

private:
  void           new_rep(int dim)  { PTR = new _d_tuple<Point_d>(dim+1); } 
  _d_tuple< Point_d >* ptr() const { return (_d_tuple< Point_d >*)PTR; }
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Simplex_d.C>
#endif 

#endif
