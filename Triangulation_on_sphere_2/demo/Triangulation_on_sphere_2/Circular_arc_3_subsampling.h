#ifndef CIRCULAR_ARC_3_SUBSAMPLING_H
#define CIRCULAR_ARC_3_SUBSAMPLING_H

#include <list>
#include <CGAL/basic_classes.h>
#include <CGAL/Linear_algebraCd.h>
#include <cmath>


#define MIN_EDGE_SIZE 0.05

template <class Kernel>
double get_theta( typename Kernel::Point_3& pt,
                  typename Kernel::Vector_3& V1,
                  typename Kernel::Vector_3& V2,
                  typename Kernel::Vector_3& V3)
{
  typedef CGAL::Linear_algebraCd<typename Kernel::FT> LA;
  std::list<typename LA::Vector> lvect;
  lvect.push_back( typename LA::Vector(V1.cartesian_begin(),V1.cartesian_end()) );
  lvect.push_back( typename LA::Vector(V2.cartesian_begin(),V2.cartesian_end()) );
  lvect.push_back( typename LA::Vector(V3.cartesian_begin(),V3.cartesian_end()) );

  typename LA::Matrix M(lvect.begin(),lvect.end());
  
  typename Kernel::FT det=1;
  typename LA::Matrix inverse=LA::inverse(M,det);

  typename LA::Vector pt_in_vect(pt.cartesian_begin(),pt.cartesian_end());
  typename Kernel::FT X=inverse.row(0)*pt_in_vect;
  typename Kernel::FT Y=inverse.row(1)*pt_in_vect;
  
  double angle=atan2(Y*sqrt(V2.squared_length()),X*sqrt(V1.squared_length()));
  
  return angle;
}

template <class Kernel>
inline typename Kernel::Point_3 compute_point(const typename Kernel::Point_3& center,
                                              double radius,double theta,
                                              const typename Kernel::Vector_3& b1,
                                              const typename Kernel::Vector_3& b2)
{
  return center + radius * cos(theta) * b1/sqrt(b1.squared_length()) + radius * sin(theta) * b2/sqrt(b2.squared_length());
}

template <class Kernel,class Output_iterator>
void subsample_circular_arc_3(const typename Kernel::Circle_3& circle,
                              double source,double target,
                              const typename Kernel::Vector_3& b1,const typename Kernel::Vector_3& b2,
                              Output_iterator out_pts,double min_edge_size)
{
  double radius=sqrt(circle.squared_radius());
  if (source > target) target+=2*CGAL_PI;
  double edge_len= (target-source) * radius;
  int nb_of_segments=static_cast<int>(floor(edge_len/min_edge_size));
  
  *out_pts++=compute_point<Kernel>(circle.center(),radius,source,b1,b2);
  double step_size=(target-source)/static_cast<double>(nb_of_segments);
  double current_theta=source;
  for (int i=0; i< nb_of_segments-1 ; ++i){
    current_theta+=step_size;
    *out_pts++=compute_point<Kernel>(circle.center(),radius,current_theta,b1,b2);  
  }
  *out_pts++=compute_point<Kernel>(circle.center(),radius,target,b1,b2);
}

//Subsample from source to target seen ccw from the side of plane pointed by its orthogonal_vector()
template <class Kernel,class Output_iterator>
void subsample_circle_3(const typename Kernel::Circle_3& circle,Output_iterator out_pts,double min_edge_size=MIN_EDGE_SIZE)
{
  typename Kernel::Vector_3 b1=circle.supporting_plane().base1();
  typename Kernel::Vector_3 b2=circle.supporting_plane().base2();
  
  subsample_circular_arc_3<Kernel>(circle,0,2*CGAL_PI,b1,b2,out_pts,min_edge_size);
}

//Subsample from source to target seen ccw from the side of plane pointed by its orthogonal_vector()
template <class Kernel,class Output_iterator>
void subsample_circular_arc_3(const typename Kernel::Circle_3& circle,
                              const typename Kernel::Plane_3& plane,
                              const typename Kernel::Point_3& source,
                              const typename Kernel::Point_3& target,
                              Output_iterator out_pts,
                              double min_edge_size=MIN_EDGE_SIZE)
{
  typename Kernel::Vector_3 b1=plane.base1();
  typename Kernel::Vector_3 b2=plane.base2();
  typename Kernel::Vector_3 b3=plane.orthogonal_vector();
  
  
  typename Kernel::Point_3 tmp=CGAL::ORIGIN + typename Kernel::Vector_3(circle.center(),source);
  double theta_source=get_theta<Kernel>(tmp,b1,b2,b3 );
  tmp=CGAL::ORIGIN + typename Kernel::Vector_3(circle.center(),target);
  double theta_target=get_theta<Kernel>(tmp,b1,b2,b3 );
  subsample_circular_arc_3<Kernel>(circle,theta_source,theta_target,b1,b2,out_pts,min_edge_size);
}



#endif //CIRCULAR_ARC_3_SUBSAMPLING_H
