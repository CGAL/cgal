// Author     : Nader Salman


// an object used in Gyroviz_cdt2 to store the vector
// of cdt2 into a vector of Gyroviz_triangle_with_cam

#ifndef _Gyroviz_triangle_with_cam_
#define _Gyroviz_triangle_with_cam_

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

template < class Triangulation >
class Gyroviz_triangle_with_cam
{

private:
	// Triangle_3
	typedef typename Triangulation::Triangle_3 Triangle_3;
	Triangle_3 triangle;


	// Camera number
	int camera_number;

public:

	Gyroviz_triangle_with_cam(){}
	Gyroviz_triangle_with_cam(Triangle_3 t, int cam_num):triangle(t),camera_number(cam_num){}

	// accessors
	const	Triangle_3&	get_triangle3() const { return	triangle; }
	const	int	    get_camera_number() const { return	camera_number;	}

	// modificators
	void	set_triangle3(Triangle_3 t)       { triangle = t; }
	void	set_camera_number(int camera_num) { camera_number =  camera_num;}


};

#endif // _Gyroviz_triangle_with_cam_