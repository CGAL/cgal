#include <CGAL/basic.h>
#include <CGAL/enum.h>
using CGAL::Comparison_result;
template <typename Traits , typename Input_iterator>
CGAL::Orientation polygon_orientation( Input_iterator begin , Input_iterator end , const Traits& traits)
{
	const typename Traits::Compare_xy_2 cmp_xy = traits.compare_xy_2_object();
	const typename Traits::Compare_y_at_x_right_2 cmp_y_at_x_right = traits.compare_y_at_x_right_2_object();
	const typename Traits::Compare_endpoints_xy_2 cmp_endpoints = traits.compare_endpoints_xy_2_object();
	const typename Traits::Construct_min_vertex_2 ctr_min =	traits.construct_min_vertex_2_object();
	// Locate the smallest polygon vert ex .
	Input_iterator in_min , out_min;
	const typename Traits::Point_2* min_point = NULL;
	// Go over the x-monotone curves t hat compose the polygon boundary .
	Input_iterator curr = end ; --curr ;
	Input_iterator next = begin ;
	bool next_directed_right = (cmp_endpoints (* curr ) == CGAL::SMALLER) ;
	for ( ; next != end ; curr = next++) 
	{
		bool curr_directed_right = next_directed_right ;
		next_directed_right = (cmp_endpoints (*next ) == CGAL::SMALLER) ;
		// I t i s s u f f i c i e n t to look for a pair of curves with opposite d i r e c t i o n s .
		if ( curr_directed_right || !next_directed_right ) continue;
		const typename Traits::Point_2& point = ctr_min (* curr );
		Comparison_result res;
		if(( min_point == NULL) || (( res = cmp_xy(point,*min_point)) == CGAL::SMALLER))
		{
			in_min = curr ;
			out_min = next ;
			min_point = &point ;
			continue;
		}
		//Handle the degenerate case where 2 polygon v e r t i c e s coincide .
		if (res == CGAL::EQUAL)
		{
			Comparison_result res_in = cmp_y_at_x_right(* curr , *in_min , point );
			Comparison_result res_out = cmp_y_at_x_right(*next , *out_min , point);
			CGAL_assertion((res_in != CGAL::EQUAL) && ( res_out != CGAL::EQUAL));
			if(((res_in == CGAL::SMALLER) && ( res_out == CGAL::SMALLER)) || (((res_in == CGAL::SMALLER) && ( res_out == CGAL::LARGER)) && (cmp_y_at_x_right(* curr , *out_min , point ) == CGAL::SMALLER)) || ((( res_in == CGAL::LARGER) && ( res_out == CGAL::SMALLER)) && (cmp_y_at_x_right(*next , *in_min , point ) == CGAL::SMALLER)))
			{
				in_min = curr ;
				out_min = next ;
			}
		}
	}
	// Perform the comparison near the smallest vert ex .
	Comparison_result res = cmp_y_at_x_right(*in_min , *out_min , *min_point ) ;
	CGAL_assertion( res != CGAL::EQUAL) ;
	return ( res == CGAL::SMALLER) ? CGAL::CLOCKWISE : CGAL::COUNTERCLOCKWISE;
}
