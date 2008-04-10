// Author     : Nader Salman


// an object used as info in the CGAL::Triangulation_vertex_base_with_info_2

#ifndef _Gyroviz_info_for_dt2_
#define _Gyroviz_info_for_dt2_

class	Gyroviz_info_for_dt2
{

private:
  // 3D Point
  Point_3 point;

  // Camera number
  int camera_number;

public:
  Gyroviz_info_for_dt2(){}
  Gyroviz_info_for_dt2(Point_3 p, int cam_num):point(p),camera_number(cam_num){}

  // accessors
  const	Point_3	get_point3()        const { return	point; }
  const	int	    get_camera_number() const { return	camera_number;	}

  // modificators
  void	set_point3(Point_3 point3) { point = point3; }
  void	set_camera_number(int camera_num) { camera_number =  camera_num;}
};

#endif //  _Gyroviz_info_for_dt2_
