 // Copyright (c) 2009  GeometryFactory (France)
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
// Author(s)     : Laurent Rineau


#ifndef CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H
#define CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H

namespace CGAL {

namespace TriangulationProjectionTraitsCartesianFunctors {

template <class Traits>
class Coplanar_orientation_with_normal_3
{
  // private members
  const Traits& traits;
  typedef typename Traits::K K;
  typename K::Construct_vector_3 vector;

  // private type alias
  typedef typename Traits::Point_2 Point;
public:
  Coplanar_orientation_with_normal_3(const Traits& traits_)
    : traits(traits_),
      vector(traits.kernel().construct_vector_3_object())
  {
  }

  Orientation operator()(const Point& p,
			 const Point& q,
			 const Point& r) const
  {
    // same as orientation( p, q, r, p + traits.normal() )
    return enum_cast<Orientation>(CGAL_NTS sign(determinant(vector(p, q),
							    vector(p, r),
							    traits.normal())));
  }
}; // end class Coplanar_orientation_with_normal_3<Traits>

template <class Traits>
class Coplanar_side_of_oriented_circle_with_normal_3
{
  // private members
  const Traits& traits;
  typedef typename Traits::K K;
  typename K::Compute_scalar_product_3 scalar;
  typename K::Compute_squared_length_3 sq_norm;
  typename K::Construct_vector_3 vector;

  // private types aliases
  typedef typename Traits::Point_2 Point;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::FT FT;
public:
  Coplanar_side_of_oriented_circle_with_normal_3(const Traits& traits_)
    : traits(traits_),
      scalar(traits.kernel().compute_scalar_product_3_object()),
      sq_norm(traits.kernel().compute_squared_length_3_object()),
      vector(traits.kernel().construct_vector_3_object())
  {
  }

  Oriented_side operator()(const Point& p,
			   const Point& q,
			   const Point& r,
			   const Point& t) const
  {
    const Vector_3& u = traits.normal();
    const FT u2 = sq_norm(u);

//     const FT tpx = p.x() - t.x();
//     const FT tpy = p.y() - t.y();
//     const FT tpz = p.z() - t.z();
//     const FT tqx = q.x() - t.x();
//     const FT tqy = q.y() - t.y();
//     const FT tqz = q.z() - t.z();
//     const FT trx = r.x() - t.x();
//     const FT try = r.y() - t.y();
//     const FT trz = r.z() - t.z();

    const Vector_3 tp = vector(t, p);
    const Vector_3 tq = vector(t, q);
    const Vector_3 tr = vector(t, r);

    const FT tp2 = sq_norm(tp);
    const FT tq2 = sq_norm(tq);
    const FT tr2 = sq_norm(tr);

    const FT k_p = scalar(tp, u);
    const FT k_q = scalar(tq, u);
    const FT k_r = scalar(tr, u);

    return enum_cast<Orientation>
      (
       sign_of_determinant(tp.x(), tp.y(), tp.z(), (tp2 + k_p) * u2 - k_p * k_p,
			   tr.x(), tr.y(), tr.z(), (tr2 + k_r) * u2 - k_r * k_r,
			   tq.x(), tq.y(), tq.z(), (tq2 + k_q) * u2 - k_q * k_q,
  			    u.x(),  u.y(),  u.z(), u2 * u2)
       );
    // Note that q and r have been swapped in the determinant above.
  }
}; // end class Coplanar_side_of_oriented_circle_with_normal_3

template <class Traits>
class Coplanar_square_distance_to_line_with_normal_3
{
  // private members
  const Traits& traits;
  typedef typename Traits::K K;
  typename K::Compute_scalar_product_3 scalar;
  typename K::Compute_squared_length_3 sq_norm;
  typename K::Construct_cross_product_vector_3 cross;
  typename K::Construct_vector_3 vector;
  
  // private types aliases
  typedef typename Traits::Point_2 Point;
  typedef typename Traits::Line_2 Line;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::FT FT;

public:
  Coplanar_square_distance_to_line_with_normal_3(const Traits& traits_)
    : traits(traits_),
      scalar(traits.kernel().compute_scalar_product_3_object()),
      sq_norm(traits.kernel().compute_squared_length_3_object()),
      cross(traits.kernel().construct_cross_product_vector_3_object()),
      vector(traits.kernel().construct_vector_3_object())
  {}

  FT operator()(const Line& line, const Point& p)
  {
    const FT den = sq_norm(cross(traits.normal(),line.to_vector()));

    if(den == FT(0)) {
      // den == 0 if the line is vertical
      // In that case, the distance is the distance to the line
      return sq_norm(cross(vector(p, line.point()),
                           line.to_vector()))
        / sq_norm(line.to_vector());
    }
    else {
      const FT det = determinant(traits.normal(), 
                                 line.to_vector(), 
                                 vector(p, line.point()));
      return det * det / den;
    }
  }
}; // end class Coplanar_square_distance_to_line_with_normal_3

template <class Traits>
class Coplanar_intersect_3
{
  // private members
  const Traits& traits;
  typedef typename Traits::K K;
  typename K::Compute_scalar_product_3 scalar;
  typename K::Compute_squared_length_3 sq_norm;
  typename K::Construct_cross_product_vector_3 cross;
  typename K::Construct_vector_3 vector;
  typename K::Intersect_3 intersect;
  
  // private types aliases
  typedef typename Traits::Point_2 Point;
  typedef typename Traits::Line_2 Line;
  typedef typename Traits::Segment_2 Segment;
  typedef typename K::Plane_3 Plane_3;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::FT FT;
public:
  Coplanar_intersect_3(const Traits& traits_)
    : traits(traits_),
      scalar(traits.kernel().compute_scalar_product_3_object()),
      sq_norm(traits.kernel().compute_squared_length_3_object()),
      cross(traits.kernel().construct_cross_product_vector_3_object()),
      vector(traits.kernel().construct_vector_3_object()),
      intersect(traits.kernel().intersect_3_object())
  {}

  Object operator()(const Segment& s1, const Segment& s2)
  {
    const Vector_3 u1 = cross(s1.to_vector(), traits.normal());
    if(u1 == NULL_VECTOR)
      return intersect(s1.supporting_line(), s2);

    const Vector_3 u2 = cross(s2.to_vector(), traits.normal());
    if(u2 == NULL_VECTOR)
      return intersect(s1, s2.supporting_line());
    
    const Plane_3 plane_1(s1.source(), u1);
    const Plane_3 plane_2(s2.source(), u2);
    
    Object planes_intersection = intersect(plane_1, plane_2);
    if(planes_intersection.empty()) {
      std::cerr << "planes_intersection is empty\n";
      return planes_intersection;
    }
    if(const Line* line = object_cast<Line>(&planes_intersection))
    {
      const Point& pi = line->point(0);
      if(scalar(cross(traits.normal(), vector(s1.source(), pi)),
                cross(traits.normal(), vector(s1.target(), pi))) > FT(0)
         ||
         scalar(cross(traits.normal(), vector(s2.source(), pi)),
                cross(traits.normal(), vector(s2.target(), pi))) > FT(0) )
      {
        // the intersection of the lines is not inside the segments
        std::cerr << "intersection not inside\n";
        return Object();
      }
      else
      {
        // Let the plane passing through s1.source() and with normal
        // the cross product of s1.to_vector() and s2.to_vector(). That
        // plane should intersect *l, now.
        return intersect(*line, Plane_3(s1.source(), cross(s1.to_vector(),
                                                           s2.to_vector())));
      }
    }
    if(object_cast<Plane_3>(&planes_intersection))
    {
      std::cerr << "coplanar lines\n";
      CGAL_error();
      return Object();
    }
    return Object();
  }
}; // end class Coplanar_intersect_3

} // end namespace TriangulationProjectionTraitsCartesianFunctors

template <class R>
class  Intersect_xy_3
{
public:
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Segment_3 Segment_3;
  typedef typename R::Point_2   Point_2; 
  typedef typename R::Vector_2  Vector_2; 
  typedef typename R::Segment_2 Segment_2;
  
  typename R::FT x(const Point_3 &p) const { return p.x(); }
  typename R::FT y(const Point_3 &p) const { return p.y(); }

  Point_2 project(const Point_3& p)
  {
    return Point_2(x(p),y(p));
  } 

#if USE_ROBUST_INTERSECTIONS
  Object operator()(const Segment_3& s1, const Segment_3& s2)
  {
    typedef typename CGAL::Exact_predicates_exact_constructions_kernel::Point_2 ExactPoint_2;
    typedef typename CGAL::Exact_predicates_exact_constructions_kernel::Segment_2 ExactSegment_2;
    
    
    

    ExactPoint_2 s1p(s1.source().x(), s1.source().y());
    ExactPoint_2 t1p(s1.target().x(), s1.target().y());
    ExactSegment_2 s1_2(s1p, t1p);
    ExactSegment_2 s2_2(ExactPoint_2(s2.source().x(), s2.source().y()),
			ExactPoint_2(s2.target().x(), s2.target().y()));
    Object o = intersection(s1_2,s2_2);
    ExactPoint_2 pi;
    if(assign(pi,o)){
      double l1 = std::sqrt(to_double(squared_distance(s1p,t1p)));
      double l2 = std::sqrt(to_double(squared_distance(s1p,pi)));
      double ratio = l2/l1;
      double z = s1.source().z() + ratio * (s1.target().z() - s1.source().z());
      Point_3 res(to_double(pi.x()), to_double(pi.y()), z);
      return make_object(res);
    } else {
      std::cerr << "NOT YET IMPLEMENTED: Intersection is not a point" << std::endl;
      Point_3 res;
      return make_object(res);
    }
  }

#else 

  Object operator()(const Segment_3& s1, const Segment_3& s2)
  {
    Point_2 s1p = project(s1.source());
    Point_2 t1p = project(s1.target());
    Segment_2 s1_2(s1p, t1p);
    Segment_2 s2_2(project(s2.source()), project(s2.target()));
    Object o = intersection(s1_2,s2_2);
    Point_2 pi;
    if(assign(pi,o)){
      double l1 = std::sqrt(to_double(squared_distance(s1p,t1p)));
      double l2 = std::sqrt(to_double(squared_distance(s1p,pi)));
      double ratio = l2/l1;
      Point_3 p = s1.source() + ratio * (s1.target() - s1.source());
      Point_3 res(pi.x(), pi.y(), p.z());
      return make_object(res);
    } else {
      std::cerr << "NOT YET IMPLEMENTED: Intersection is not a point" << std::endl;
      return Object();
    }
  }
#endif
};

template <class R>
class Compare_distance_xy_3
{
public:
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Point_2   Point_2;   
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return p.x(); }
  typename R::FT y(const Point_3 &p) const { return p.y(); }

  Point_2 project(const Point_3& p)
  {
    return Point_2(x(p),y(p));
  }

  Comparison_result operator()(const Point_3& p,const Point_3& q,const Point_3& r)
  {
    Point_2 p2 = project(p);
    Point_2 q2 = project(q);
    Point_2 r2 = project(r);
    return compare_distance_to_point(p2,q2,r2);
  }
};


template <class R>
class Squared_distance_xy_3
{
public:
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Point_2   Point_2; 
  typedef typename R::Line_3    Line_3; 
  typedef typename R::Line_2    Line_2;
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return p.x(); }
  typename R::FT y(const Point_3 &p) const { return p.y(); }

  Point_2 project(const Point_3& p)
  {
    return Point_2(x(p),y(p));
  }

  RT operator()(const Line_3& l, const Point_3& p)
  {
    Point_2 p2 = project(p);
    Line_2 l2 = Line_2(project(l.point(0)), project(l.point(1)));
    return squared_distance(p2, l2);
  }
};

// template <class R>
// class Squared_radius_xy_3
// {
// public:
//   typedef typename R::Point_3   Point_3; 
//   typedef typename R::Point_2   Point_2; 
//   typedef typename R::Line_3    Line_3; 
//   typedef typename R::Line_2    Line_2;
//   typedef typename R::FT        RT;
//   typename R::FT x(const Point_3 &p) const { return p.x(); }
//   typename R::FT y(const Point_3 &p) const { return p.y(); }

//   Point_2 project(const Point_3& p)
//   {
//     return Point_2(x(p),y(p));
//   }

// #if USE_ROBUST_SQUARED_RADIUS
//   RT operator()(const Point_3& p, const Point_3& q, const Point_3& r)
//   {
//     typedef typename CGAL::Exact_predicates_exact_constructions_kernel::Point_2 ExactPoint_2;
//     typedef typename CGAL::Exact_predicates_exact_constructions_kernel::FT ExactFT;
//     typedef typename CGAL::Exact_predicates_exact_constructions_kernel::Compute_squared_radius_2 ExactSqRadius_2;

//     ExactPoint_2 pp(p.x(), p.y());
//     ExactPoint_2 qq(q.x(), q.y());
//     ExactPoint_2 rr(r.x(), r.y());

//     ExactFT exact_sq_radius = exact_comp_sqradius_2(pp, qq, rr);
//     return RT(CGAL::to_double(exact_sq_radius));
//   }

//   RT operator()(const Point_3& p, const Point_3& q)
//   {
//     typedef typename CGAL::Exact_predicates_exact_constructions_kernel::Point_2 ExactPoint_2;
//     typedef typename CGAL::Exact_predicates_exact_constructions_kernel::FT ExactFT;
//     typedef typename CGAL::Exact_predicates_exact_constructions_kernel::Compute_squared_radius_2 ExactSqRadius_2;

//     ExactPoint_2 pp(p.x(), p.y());
//     ExactPoint_2 qq(q.x(), q.y());

//     ExactFT exact_sq_radius = exact_comp_sqradius_2(pp, qq);
//     return RT(CGAL::to_double(exact_sq_radius));
//   }
// #else
//   RT operator()(const Point_3& p, const Point_3& q, const Point_3& r)
//   {
//     Point_2 p2 = project(p);
//     Point_2 q2 = project(q);
//     Point_2 r2 = project(r);
//     typename R::Compute_squared_radius_2 = R().compute_squared_radius_2_object();
//     return squared_radius(p2, q2, r2);
//   }

//   RT operator()(const Point_3& p, const Point_3& q)
//   {
//     Point_2 p2 = project(p);
//     Point_2 q2 = project(q);
//     typename R::Compute_squared_radius_2 = R().compute_squared_radius_2_object();
//     return squared_radius(p2, q2, r2);
//   }
// #endif
// };

namespace TODO {

template <class R>
class Side_of_oriented_circle_xy_3 
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT x(const Point &p) const { return p.x(); }
  typename R::FT y(const Point &p) const { return p.y(); }

  CGAL::Oriented_side operator() (const Point &p, 
				  const Point &q,
				  const Point &r, 
				  const Point &s) 
    {
      typename R::Side_of_oriented_circle_2 side_of_oriented_circle_2 = R().side_of_oriented_circle_2_object();
#ifdef CGAL_3_2

      return side_of_oriented_circle_2(p,
				       q,
				       r,
				       s);
#else
      typename R::Point_2 pp(x(p), y(p)), qq(x(q), y(q)), rr(x(r), y(r)), ss(x(s), y(s));
      return side_of_oriented_circle_2(pp,
				       qq,
				       rr,
				       ss);
#endif
    }
};

template <class R>
class Side_of_bounded_circle_xy_3
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT x(const Point &p) const { return p.x(); }
  typename R::FT y(const Point &p) const { return p.y(); }

  CGAL::Oriented_side operator() (const Point &p, 
				  const Point &q,
				  const Point &r, 
				  const Point &s) 
    {
      typename R::Side_of_bounded_circle_2 side_of_bounded_circle_2 = R().side_of_bounded_circle_2_object();
      typename R::Point_2 pp(x(p), y(p)), qq(x(q), y(q)), rr(x(r), y(r)), ss(x(s), y(s));
      return side_of_bounded_circle_2(pp,
				      qq,
				      rr,
				      ss);
    }
};
} // end namespace TODO

template < class Kernel >
class Triangulation_2_projection_traits_3
{
  typedef Triangulation_2_projection_traits_3<Kernel> Self;

  typename Kernel::Vector_3 n;
  Kernel k;

public:
  typedef typename Kernel::Vector_3 Vector_3;


  Triangulation_2_projection_traits_3(const Vector_3& n_, 
				      Kernel k_ = Kernel())
    : n(n_), k(k_)
  {}

  const Vector_3& normal() const
  {
    return n;
  }

  Kernel& kernel()
  {
    return k;
  }

  const Kernel& kernel() const
  {
    return k;
  }

  Triangulation_2_projection_traits_3(const Self& other)
    : n(other.n)
  {}

  Self& operator=(const Self& other)
  {
    if(this != &other) {
      n = other.n;
    }
    return *this;
  }

  typedef Kernel K;
  typedef typename K::FT          FT;
  typedef typename K::Point_3     Point_2;
  typedef typename K::Segment_3   Segment_2;
  typedef typename K::Triangle_3  Triangle_2;
  typedef typename K::Line_3      Line_2;

  // Maybe not a good choice
  typedef typename K::Less_xy_3            Less_x_2;
  typedef typename K::Less_z_3             Less_y_2;

//   class Compare_x_2 {
//     const Self& traits;
//   public:
//     Compare_x_2(const Self& traits_) : traits(traits_) {}

//     Comparison_result operator()(const typename K::Point_3& p,
// 				 const typename K::Point_3& q) {
//       if(determinant(p-ORIGIN, q-ORIGIN, traits.normal()) == FT(0))
// 	return EQUAL;
//       else
// 	return typename K::Compare_xy_3()(p,q);
//     }
//   };

//   class Compare_y_2 {
//     const Self& traits;
//   public:
//     Compare_y_2(const Self& traits_) : traits(traits_) {}

//     Comparison_result operator()(const typename K::Point_3& p,
// 				 const typename K::Point_3& q) {
//       if(collinear(p, q, p+traits.normal()))
// 	return EQUAL;
//       else
// 	return compare_z(p,q);
//     }
//   };

  typedef typename K::Compare_xy_3                           Compare_x_2;
  typedef typename K::Compare_z_3                            Compare_y_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
    Coplanar_orientation_with_normal_3<Self>                 Orientation_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
    Coplanar_side_of_oriented_circle_with_normal_3<Self>     Side_of_oriented_circle_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
  Coplanar_square_distance_to_line_with_normal_3<Self>       Compute_squared_distance_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
  Coplanar_intersect_3<Self>                                 Intersect_2;
//   typedef Projection_traits::Squared_distance_xy_3<K>         Compute_squared_distance_2;
//   typedef Projection_traits::Squared_radius_xy_3<K>           Compute_squared_radius_2;
//   typedef Projection_traits::Compare_distance_xy_3<K>         Compare_distance_2;
  typedef typename K::Construct_segment_3  Construct_segment_2;
  typedef typename K::Construct_line_3     Construct_line_2;
  typedef typename K::Construct_triangle_3 Construct_triangle_2;
//   typedef Intersect_xy_3<K>                Intersect_2;
  
  Less_x_2
  less_x_2_object() const
    { return Less_x_2();}

  Less_y_2
  less_y_2_object() const
    { return Less_y_2();}

  Compare_x_2
  compare_x_2_object() const
  {
    return Compare_x_2(// *this
		       );
  }

  Compare_y_2
  compare_y_2_object() const
  { 
    return Compare_y_2(// *this
		       );
  }
  
//   Compare_distance_2
//   compare_distance_2_object() const
//   {
//     return Compare_distance_2();
//   }
  
  Orientation_2 
  orientation_2_object() const
  {
    return Orientation_2(*this); 
  }

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  {
    return Side_of_oriented_circle_2(*this);
  }
  
  Compute_squared_distance_2
  compute_squared_distance_2_object() const
  {
    return Compute_squared_distance_2(*this);
  }
  
  Intersect_2
  intersect_2_object () const
  {
    return Intersect_2(*this);
  }

//   Compute_squared_distance_2
//   compute_squared_distance_2_object () const
//   {
//     return Compute_squared_distance_2();
//   }
  
  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}
  
  Construct_line_2  construct_line_2_object() const
    {return Construct_line_2();}
  
  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}

}; // end class Triangulation_2_projection_traits_3<Kernel>
  

} // end namespace CGAL

#endif // CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H
