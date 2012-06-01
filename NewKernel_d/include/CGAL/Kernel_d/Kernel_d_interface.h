#ifndef CGAL_KD_KERNEL_3_INTERFACE_H
#define CGAL_KD_KERNEL_3_INTERFACE_H

#include <CGAL/functor_tags.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/marcutils.h>
#include <CGAL/tuple.h>


namespace CGAL {
template <class Base_> struct Kernel_d_interface : public Base_ {
	typedef Base_ Base;
	typedef Kernel_d_interface<Base> Kernel;
	typedef typename Base::Flat_orientation Flat_orientation_d;
	typedef typename Base::Point Point_d;
	typedef typename Base::Vector Vector_d;
	//typedef typename Base::Segment Segment_d;
	typedef typename Base::template Functor<Compute_point_cartesian_coordinate_tag>::type Compute_coordinate_d;
	typedef typename Base::template Functor<Compare_lexicographically_tag>::type Compare_lexicographically_d;
	typedef typename Base::template Functor<Equal_points_tag>::type Equal_d;
	typedef typename Base::template Functor<Less_lexicographically_tag>::type Less_lexicographically_d;
	typedef typename Base::template Functor<Less_or_equal_lexicographically_tag>::type Less_or_equal_lexicographically_d;
	typedef typename Base::template Functor<Orientation_of_points_tag>::type Orientation_d;
	typedef typename Base::template Functor<Point_dimension_tag>::type Point_dimension_d;
	typedef typename Base::template Functor<Side_of_oriented_sphere_tag>::type Side_of_oriented_sphere_d;
	typedef typename Base::template Functor<Contained_in_affine_hull_tag>::type Contained_in_affine_hull_d;
	typedef typename Base::template Functor<Construct_flat_orientation_tag>::type Construct_flat_orientation_d;
	typedef typename Base::template Functor<In_flat_orientation_tag>::type In_flat_orientation_d;
	typedef typename Base::template Functor<In_flat_side_of_oriented_sphere_tag>::type In_flat_side_of_oriented_sphere_d;
	typedef typename Base::template Functor<Point_to_vector_tag>::type Point_to_vector_d;
	typedef typename Base::template Functor<Vector_to_point_tag>::type Vector_to_point_d;
	struct Component_accessor_d : private Store_kernel<Kernel_d_interface> {
	  CGAL_FUNCTOR_INIT_STORE(Component_accessor_d)
	  typedef typename Base::Point Point;
	  typedef typename Base::RT RT;
	  typedef typename Base::FT FT;
	  int dimension(Point const&p){
	    return this->kernel().point_dimension_d_object()(p);
	  }
	  FT cartesian(Point const&p, int i){
	    return this->kernel().compute_coordinate_d_object()(p);
	  }
	  RT homogeneous(Point const&p, int i){
	    // FIXME
	    return this->kernel().compute_coordinate_d_object()(p);
	  }
	};


	Compute_coordinate_d compute_coordinate_d_object()const{ return Compute_coordinate_d(*this); }
	Compare_lexicographically_d compare_lexicographically_d_object()const{ return Compare_lexicographically_d(*this); }
	Equal_d equal_d_object()const{ return Equal_d(*this); }
	Less_lexicographically_d less_lexicographically_d_object()const{ return Less_lexicographically_d(*this); }
	Less_or_equal_lexicographically_d less_or_equal_lexicographically_d_object()const{ return Less_or_equal_lexicographically_d(*this); }
	Orientation_d orientation_d_object()const{ return Orientation_d(*this); }
	Point_dimension_d point_dimension_d_object()const{ return Point_dimension_d(*this); }
	Side_of_oriented_sphere_d side_of_oriented_sphere_d_object()const{ return Side_of_oriented_sphere_d(*this); }
	Contained_in_affine_hull_d contained_in_affine_hull_d_object()const{ return Contained_in_affine_hull_d(*this); }
	Construct_flat_orientation_d construct_flat_orientation_d_object()const{ return Construct_flat_orientation_d(*this); }
	In_flat_orientation_d in_flat_orientation_d_object()const{ return In_flat_orientation_d(*this); }
	In_flat_side_of_oriented_sphere_d in_flat_side_of_oriented_sphere_d_object()const{ return In_flat_side_of_oriented_sphere_d(*this); }
	Point_to_vector_d point_to_vector_d_object()const{ return Point_to_vector_d(*this); }
	Vector_to_point_d vector_to_point_d_object()const{ return Vector_to_point_d(*this); }
	Component_accessor_d component_accessor_d_object()const{ return Component_accessor_d(*this); }
};
}

#endif
