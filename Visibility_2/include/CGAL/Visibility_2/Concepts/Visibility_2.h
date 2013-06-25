#ifndef CGAL_VISIBILITY_2_H
#define CGAL_VISIBILITY_2_H

#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

namespace CGAL {

namespace Visibility_2 {

template<class Arrangement_2> 
class Visibility_2 
{
public:
	typedef typename Arrangement_2::Geometry_traits_2		Geometry_traits_2;
	// Currently only consider with same type for both
	typedef Input_Arrangement_2								Arrangement_2;
	typedef Output_Arrangement_2							Arrangement_2;

	typedef typename Arrangement_2::Halfedge_const_handle	Halfedge_const_handle;
	typedef typename Arrangement_2::Face_const_handle		Face_const_handle;

	typedef typename Geometry_traits_2::Point_2				Point_2;

	// TODO: Add RegularizationTag

	Visibility_2(const Input_Arrangement_2 &arr/*, Regularization_tag r_t*/): p_arr(arr) {};

	bool is_attached() {
		return (p_arr != NULL);
	}

	void attach(const Input_Arrangement_2 &arr) {
		p_arr = arr;
	}

	void detach() {
		p_arr = NULL;
	}

	Input_Arrangement_2 arr() {
		return (*p_arr);
	}

	void visibility_region(const Point_2 &q, 
						   const Face_const_handle f,
						   Output_Arrangement_2 &out_arr
						   );
	void visibility_region(const Point_2 &q, 
						   const Halfedge_const_handle he,
						   Output_Arrangement_2 &out_arr
						   );
protected:
	Input_Arrangement_2 *p_arr;
};

} // namespace Visibility_2
} // namespace CGAL

#endif