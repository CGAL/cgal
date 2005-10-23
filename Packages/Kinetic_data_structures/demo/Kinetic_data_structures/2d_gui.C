#include <CGAL/KDS/Exact_simulation_traits_2.h>
#include <CGAL/KDS/Enclosing_box_2.h>


#ifdef CGAL_USE_QT
#include <CGAL/KDS/IO/Qt_widget_2.h>
#include <CGAL/KDS/IO/Qt_moving_points_2.h>
#endif

#include <CGAL/KDS/Insert_event.h>

/*!
  \file coin_check.cc A simple example using the qt GUI.
*/

int main(int argc, char*argv[]){

#ifdef CGAL_USE_QT
  CGAL_KDS_SET_LOG_LEVEL(CGAL::KDS::LOG_LOTS);

  typedef CGAL::KDS::Exact_simulation_traits_2 Traits;

  typedef CGAL::KDS::Qt_widget_2<Traits::Simulator> Gui;
  typedef CGAL::KDS::Qt_moving_points_2<Traits, Gui> Qt_moving_points;
  typedef CGAL::KDS::Insert_event<Traits::Moving_point_table> Insert_event;
  typedef Traits::Kinetic_kernel::Point_2 Moving_point;
  typedef CGAL::KDS::Enclosing_box_2<Traits> Box;
  Traits tr;

  Gui::Pointer qtsim=new Gui(argc, argv, tr.simulator_pointer());
  Qt_moving_points::Pointer qtmptp= new Qt_moving_points(qtsim, tr);
  Box::Pointer box= new Box(tr);
  
 
  Traits::Function_kernel::Construct_function cf= tr.function_kernel_object().construct_function_object();

  tr.simulator_pointer()->new_event(Traits::Simulator::Time(0), Insert_event( Moving_point(cf(0,2,4),
						     cf(2,-3,4)), tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Traits::Simulator::Time(.1), Insert_event( Moving_point(cf(3,2,4),
						      cf(1,-3,4)), tr.moving_point_table_pointer()));

  return qtsim->begin_event_loop();
#else
  bool warning_qt_is_not_available_no_2d_visualization;
  return 0;
#endif
}
