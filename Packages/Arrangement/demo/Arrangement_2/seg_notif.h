





#ifndef SEG_NOTIF
#define SEG_NOTIF
#include"cgal_types1.h"



class Seg_notification : public Seg_arr_change_notification 
{
public:

	//default constructor
	Seg_notification(){}
	
  void add_edge(const  Traits::X_monotone_curve_2 & x_mon_curve,
                Planar_map::Halfedge_handle halfedge_handle, 
                bool  left_to_right , bool overlap = false)
  {
   
    std::cout << "add_edge" << std::endl;

	// we want to remove the const out of x_mon_curve in order to change its data
	Traits::X_monotone_curve_2 & my_x_mon_curve = const_cast<Traits::X_monotone_curve_2&>(x_mon_curve); 

	Curve_data d = my_x_mon_curve.get_data();
	d.halfedge_handle = halfedge_handle;
	my_x_mon_curve.set_data( d );    
    
  }

  void split_edge(Planar_map::Halfedge_handle  orig_edge , 
                  Planar_map::Halfedge_handle  new_edge ,
                  const Traits::X_monotone_curve_2 & c1,
                  const Traits::X_monotone_curve_2 & c2)
  {
    std::cout << "split_edge" << std::endl;

	

	// we want to remove the const out of c1 and c2 in order to change their data
    Traits::X_monotone_curve_2 & my_c1 = const_cast<Traits::X_monotone_curve_2&>(c1); 
    Traits::X_monotone_curve_2 & my_c2 = const_cast<Traits::X_monotone_curve_2&>(c2); 

	// updating the new curves

	Curve_data d1 = c1.get_data();
	d1.halfedge_handle = orig_edge;
	my_c1.set_data( d1 );   

	Curve_data d2 = c2.get_data();
	d2.halfedge_handle = new_edge;
	my_c2.set_data( d2 );   
	
	return;
  }

  void split_face(Planar_map::Face_handle  orig_face , 
                  Planar_map::Face_handle  new_face )
  {
    std::cout << "split_face" << std::endl;
	new_face->set_info( orig_face->info());  // the new face's color will be the same
											 //as the original face

  }

  void add_hole(Planar_map::Face_handle /* in_face */, 
                Planar_map::Halfedge_handle /* new_hole */)
  {
    std::cout << "add_hole" << std::endl;
  }

  
};

#endif //SEG_NOTIF

