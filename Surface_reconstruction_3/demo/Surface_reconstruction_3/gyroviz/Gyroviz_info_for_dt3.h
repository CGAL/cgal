// Author     : Nader Salman


// an object used as info in the CGAL::Triangulation_vertex_base_with_info_3
// to be revised !!

#ifndef _Gyroviz_info_for_dt3_
#define _Gyroviz_info_for_dt3_

#include <vector>

class	Gyroviz_info_for_dt3
{

protected:

  // List of images
  std::vector<Point_3> list_of_cameras;
  double greatest_camera_angle;

public:
  Gyroviz_info_for_dt3(){}
  
  Gyroviz_info_for_dt3(Point_3 position, const std::vector<Point_3>& v):list_of_cameras(v)
  {
	  // give a score to each vertex: the score will help us to detect 			  
	  // the outliers in the point cloud extracted by Voodoo
	  double vertex_score = 0;
	  Vector_3 v1, v2;
	  double v1_v2, n_v1, n_v2;
	  double intermediate_score;

	  for(int i=0; i<list_of_cameras.size()-1; ++i)
	  {
		  for(int j=i+1; j<list_of_cameras.size(); ++j)
		  {
			  v1 = list_of_cameras[i] - position;
			  v2 = list_of_cameras[j] - position;
			  n_v1  = sqrt(v1.squared_length());  // returns the length
			  n_v2  = sqrt(v2.squared_length()); 
			  v1_v2 = v1 * v2; // returns the scalar product
			  intermediate_score = acos(v1_v2/(n_v1*n_v2));

			  if(intermediate_score > vertex_score)
				  vertex_score = intermediate_score;
		  }
	  }
	  greatest_camera_angle =  vertex_score;
  }

  // accessors
  const	std::vector<Point_3>&	  get_list_of_cameras() const { return	list_of_cameras; }
  double	  get_greatest_camera_angle() const { return	greatest_camera_angle; }
};

#endif // _Gyroviz_info_for_dt3_