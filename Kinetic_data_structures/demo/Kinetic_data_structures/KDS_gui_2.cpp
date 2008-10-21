#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Enclosing_box_2.h>

#include <CGAL/Kinetic/IO/Qt_widget_2.h>
#include <CGAL/Kinetic/IO/Qt_moving_points_2.h>

#include <CGAL/Kinetic/Insert_event.h>

/*!
  \file coin_check.cc A simple example using the qt GUI.
*/

int main(int argc, char*argv[])
{

  CGAL_SET_LOG_LEVEL(CGAL::Log::LOTS);

  typedef CGAL::Kinetic::Exact_simulation_traits Traits;

  typedef CGAL::Kinetic::Qt_widget_2<Traits::Simulator> Gui;
  typedef CGAL::Kinetic::Qt_moving_points_2<Traits, Gui> Qt_moving_points;
  typedef CGAL::Kinetic::Insert_event<Traits::Active_points_2_table> Insert_event;
  typedef Traits::Kinetic_kernel::Point_2 Moving_point;
  typedef CGAL::Kinetic::Enclosing_box_2<Traits> Box;
  Traits tr(0,100000.0);

  Gui::Handle qtsim=new Gui(argc, argv, tr.simulator_handle());
  Qt_moving_points::Handle qtmptp= new Qt_moving_points(qtsim, tr);
  Box::Handle box= new Box(tr);

  Traits::Kinetic_kernel::Function_kernel::Construct_function cf= tr.kinetic_kernel_object().function_kernel_object().construct_function_object();

  tr.simulator_handle()->new_event(Traits::Simulator::Time(0), Insert_event( Moving_point(cf(0,2,4),
											  cf(2,-3,4)), tr.active_points_2_table_handle()));
  tr.simulator_handle()->new_event(Traits::Simulator::Time(.1),
				   Insert_event( Moving_point(cf(3,2,4),
							      cf(1,-3,4)),
						 tr.active_points_2_table_handle()));

  return qtsim->begin_event_loop();
}
