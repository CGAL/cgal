// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_PREDICATES_ON_POINTS_3_H
#define CGAL_CARTESIAN_PREDICATES_ON_POINTS_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/predicates/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
Comparison_result
compare_lexicographically_xyz(const PointC3<R CGAL_CTAG> &p,
                              const PointC3<R CGAL_CTAG> &q)
{
  return compare_lexicographically_xyzC3(p.x(),p.y(),p.z(),
                                         q.x(),q.y(),q.z());
}

template < class R >
bool
lexicographically_xyz_smaller_or_equal(const PointC3<R CGAL_CTAG> &p,
                                       const PointC3<R CGAL_CTAG> &q)
{
  return compare_lexicographically_xyz(p,q) != LARGER;
 }

template < class R >
bool
lexicographically_xyz_smaller(const PointC3<R CGAL_CTAG> &p,
                              const PointC3<R CGAL_CTAG> &q)
{
  return compare_lexicographically_xyz(p,q) == SMALLER;
}

template < class R >
inline
bool
x_equal(const PointC3<R CGAL_CTAG> &p, const PointC3<R CGAL_CTAG> &q)
{
  return p.x() == q.x();
}

template < class R >
inline
bool
y_equal(const PointC3<R CGAL_CTAG> &p, const PointC3<R CGAL_CTAG> &q)
{
  return p.y() == q.y();
}

template < class R >
inline
bool
z_equal(const PointC3<R CGAL_CTAG> &p, const PointC3<R CGAL_CTAG> &q)
{
  return p.z() == q.z();
}

template < class R >
inline
bool
equal_xy(const PointC3<R CGAL_CTAG> &p, const PointC3<R CGAL_CTAG> &q)
{
  return p.x() == q.x() && p.y() == q.y();
}

template < class R >
inline
bool
equal_xyz(const PointC3<R CGAL_CTAG> &p, const PointC3<R CGAL_CTAG> &q)
{
  return p.x() == q.x() && p.y() == q.y() && p.z() == q.z();
}


template < class R >
inline
Comparison_result
compare_x(const PointC3<R CGAL_CTAG> &p, const PointC3<R CGAL_CTAG> &q)
{
  return CGAL::compare(p.x(), q.x());
}

template < class R >
inline
Comparison_result
compare_y(const PointC3<R CGAL_CTAG> &p, const PointC3<R CGAL_CTAG> &q)
{
 return CGAL::compare(p.y(), q.y());
}

template < class R >
inline
Comparison_result
compare_z(const PointC3<R CGAL_CTAG> &p, const PointC3<R CGAL_CTAG> &q)
{
 return CGAL::compare(p.z(), q.z());
}


template < class R >
inline
bool
strict_dominance(const PointC3<R CGAL_CTAG> &p,
		 const PointC3<R CGAL_CTAG> &q)
{
  return strict_dominanceC3(p.x(), p.y(), p.z(),
			    q.x(), q.y(), q.z());
}

template < class R >
inline
bool
dominance(const PointC3<R CGAL_CTAG> &p,
	  const PointC3<R CGAL_CTAG> &q)
{
  return dominanceC3(p.x(), p.y(), p.z(),
		     q.x(), q.y(), q.z());
}


template < class R >
inline
bool
collinear(const PointC3<R CGAL_CTAG> &p,
               const PointC3<R CGAL_CTAG> &q,
               const PointC3<R CGAL_CTAG> &r)
{
  return collinearC3(p.x(), p.y(), p.z(),
                     q.x(), q.y(), q.z(),
                     r.x(), r.y(), r.z());
}

template < class R >
inline
Orientation
orientation(const PointC3<R CGAL_CTAG> &p,
            const PointC3<R CGAL_CTAG> &q,
            const PointC3<R CGAL_CTAG> &r,
            const PointC3<R CGAL_CTAG> &s)
{
  return orientationC3(p.x(), p.y(), p.z(),
                       q.x(), q.y(), q.z(),
                       r.x(), r.y(), r.z(),
                       s.x(), s.y(), s.z());
}

template < class R >
inline
bool
coplanar(const PointC3<R CGAL_CTAG> &p,
         const PointC3<R CGAL_CTAG> &q,
         const PointC3<R CGAL_CTAG> &r,
         const PointC3<R CGAL_CTAG> &s)
{
  return orientation(p, q, r, s) == COPLANAR;
}

template < class R >
inline
Orientation
coplanar_orientation(const PointC3<R CGAL_CTAG> &q,
         const PointC3<R CGAL_CTAG> &r,
         const PointC3<R CGAL_CTAG> &s,
         const PointC3<R CGAL_CTAG> &p)
{
  // p,q,r,s supposed to be coplanar                                   
  // q,r,s supposed to be non collinear                                
  // tests whether p is on the same side of q,r as s                   
  // returns :                                                         
  // COLLINEAR if pqr collinear                                        
  // POSITIVE if qrp and qrs have the same orientation                 
  // NEGATIVE if qrp and qrs have opposite orientations       
  CGAL_kernel_exactness_precondition( ! collinear(q, r, s) );
  CGAL_kernel_exactness_precondition( coplanar(p, q, r, s) );
  return coplanar_orientationC3(q.x(), q.y(), q.z(),
                                r.x(), r.y(), r.z(),
                                s.x(), s.y(), s.z(),
                                p.x(), p.y(), p.z());
}

template < class R>
inline
bool
are_positive_oriented(const PointC3<R CGAL_CTAG>& p,
                      const PointC3<R CGAL_CTAG>& q,
                      const PointC3<R CGAL_CTAG>& r,
                      const PointC3<R CGAL_CTAG>& s)
{
  return orientation(p,q,r,s) == POSITIVE;
}

template < class R>
inline
bool
are_negative_oriented(const PointC3<R CGAL_CTAG>& p,
                      const PointC3<R CGAL_CTAG>& q,
                      const PointC3<R CGAL_CTAG>& r,
                      const PointC3<R CGAL_CTAG>& s)
{
  return orientation(p,q,r,s) == NEGATIVE;
}

template < class R >
inline
bool
are_ordered_along_line(const PointC3<R CGAL_CTAG> &p,
                       const PointC3<R CGAL_CTAG> &q,
                       const PointC3<R CGAL_CTAG> &r)
{
  return (collinear(p, q, r)) ? collinear_are_ordered_along_line(p, q, r)
                              : false;
}

template < class R >
inline
bool
collinear_are_ordered_along_line(const PointC3<R CGAL_CTAG> &p,
                                 const PointC3<R CGAL_CTAG> &q,
                                 const PointC3<R CGAL_CTAG> &r)
{
  CGAL_kernel_exactness_precondition( collinear(p, q, r) );
  return collinear_are_ordered_along_lineC3(p.x(),p.y(),p.z(),
                                            q.x(),q.y(),q.z(),
                                            r.x(),r.y(),r.z());
}

template < class R >
inline
bool
are_strictly_ordered_along_line(const PointC3<R CGAL_CTAG> &p,
                                const PointC3<R CGAL_CTAG> &q,
                                const PointC3<R CGAL_CTAG> &r)
{
  return (collinear(p, q, r))
         ? collinear_are_strictly_ordered_along_line(p, q, r)
         : false;
}

template < class R >
inline
bool
collinear_are_strictly_ordered_along_line(const PointC3<R CGAL_CTAG> &p,
                                          const PointC3<R CGAL_CTAG> &q,
                                          const PointC3<R CGAL_CTAG> &r)
{
  CGAL_kernel_exactness_precondition( collinear(p, q, r) );
  return collinear_are_strictly_ordered_along_lineC3(p.x(),p.y(),p.z(),
                                                     q.x(),q.y(),q.z(),
                                                     r.x(),r.y(),r.z());
}

template <class R >
Oriented_side
side_of_oriented_sphere(const PointC3<R CGAL_CTAG> &p,
                        const PointC3<R CGAL_CTAG> &q,
                        const PointC3<R CGAL_CTAG> &r,
                        const PointC3<R CGAL_CTAG> &s,
                        const PointC3<R CGAL_CTAG> &test)
{
  return side_of_oriented_sphereC3(p.x(),p.y(),p.z(),
                                   q.x(),q.y(),q.z(),
                                   r.x(),r.y(),r.z(),
                                   s.x(),s.y(),s.z(),
                                   test.x(),test.y(),test.z());
}

template <class R >
Bounded_side
side_of_bounded_sphere(const PointC3<R CGAL_CTAG> &p,
                       const PointC3<R CGAL_CTAG> &q,
                       const PointC3<R CGAL_CTAG> &r,
                       const PointC3<R CGAL_CTAG> &s,
                       const PointC3<R CGAL_CTAG> &test)
{
  return side_of_bounded_sphereC3(p.x(),p.y(),p.z(),
                                  q.x(),q.y(),q.z(),
                                  r.x(),r.y(),r.z(),
                                  s.x(),s.y(),s.z(),
                                  test.x(),test.y(),test.z());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_POINTS_3_H
