#include <CGAL/KDS/Regular_triangulation_3.h>
#include <CGAL/KDS/Regular_triangulation_exact_simulation_traits_3.h>
#include <CGAL/KDS/IO/Qt_moving_points_3.h>
#include <CGAL/KDS/IO/Qt_triangulation_3.h>
#include <CGAL/KDS/IO/Qt_widget_3.h>
#include <CGAL/KDS/Insert_event.h>




int main(int argc, char *argv[]){

  typedef CGAL::KDS::Regular_triangulation_exact_simulation_traits_3 Traits;
  typedef CGAL::KDS::Qt_widget_3<Traits::Simulator> Qt_gui;
  typedef CGAL::KDS::Qt_moving_points_3<Traits, Qt_gui> Qt_mpt;
  typedef Traits::Kinetic_kernel::Motion_function MF;
  typedef Traits::Kinetic_kernel::Weighted_point_3 MP;
  typedef Traits::Kinetic_kernel::Point_3 MPP;

  Traits tr;
  
  CGAL_KDS_SET_LOG_LEVEL(CGAL::KDS::LOG_LOTS);

  Qt_gui::Pointer qtsim= new Qt_gui(argc, argv, tr.simulator_pointer());
  Qt_mpt::Pointer qtmpt= new Qt_mpt(tr, qtsim);
  
  typedef CGAL::KDS::Regular_triangulation_3<Traits> KDel;
  typedef CGAL::KDS::Qt_triangulation_3<KDel, Qt_gui, Qt_mpt> CoinKDel;
  KDel::Pointer kdel= new KDel(tr);
  CoinKDel::Pointer cd= new CoinKDel(kdel, qtsim, qtmpt);
  
  typedef Traits::Simulator::Time Time;
  typedef CGAL::KDS::Insert_event<Traits::Moving_point_table> MOI;

  Traits::Function_kernel::Construct_function cf= tr.function_kernel_object().construct_function_object();

  tr.simulator_pointer()->new_event(Time(0.000000), MOI(MP(MPP(cf(10, 1),
							       cf(0, .1),  
							       cf(0)), cf(.1)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000001), MOI(MP(MPP(cf(10), 
							       cf(0),  
							       cf(2)), cf(.1)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000002), MOI(MP(MPP(cf(12), 
							       cf(0),  
							       cf(0)), cf(.1)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000003), MOI(MP(MPP(cf(10), 
							    cf(2),  
							       cf(0)), cf(.1)), tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000000035), MOI(MP(MPP(cf(1), cf(10),  cf(0)), cf(.1)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000004),MOI(MP(MPP(cf(10.5, .05), 
							      cf(.5),  
							      cf(.5)), cf(.1)), tr.moving_point_table_pointer()));
  
  /*tr.simulator_pointer()->new_event(Time(0.000005), MOI(MPP(cf(0,1,-2), 
					  cf(2,5),  
					  cf(2,5)), tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000006), MOI(MPP(cf(1,3,-3), 
					  cf(3,-3),  
					  cf(3,5)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000007),MOI( MPP(cf(2,1,-1), cf(5,2),  cf(2,5)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000008), MOI(MPP(cf(3,2,3), cf(3,6),  cf(-1,5)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000009), MOI(MPP(cf(5,2,3), cf(-3,6),  cf(4,5)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000010), MOI(MPP(cf(-5), cf(-4),  cf(-6,5)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000011), MOI(MPP(cf(7,2,5), cf(-1,6),  cf(6,5)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000012), MOI(MPP(cf(3.1,-1,7), cf(-4.5,6),  cf(-2,5)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->new_event(Time(0.000013), MOI(MPP(cf(-1,2,-5), cf(2.3,6),  cf(0,5)) , tr.moving_point_table_pointer()));
  tr.simulator_pointer()->end_time();*/

  
  CGAL_KDS_SET_LOG_LEVEL(CGAL::KDS::LOG_LOTS);
  tr.simulator_pointer()->set_current_event_number(4);
  kdel->set_has_certificates(true);

  return qtsim->begin_event_loop();
};
