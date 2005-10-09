#include <CGAL/KDS/Regular_triangulation_inexact_simulation_traits_3.h>
#include <CGAL/KDS/Insert_event.h>
#include <CGAL/KDS/Regular_triangulation_3.h>


int main(int, char *[]){

  typedef CGAL::KDS::Regular_triangulation_inexact_simulation_traits_3 Traits;
  typedef Traits::Kinetic_kernel::Point_3 MP;
  typedef Traits::Kinetic_kernel::Weighted_point_3 WMP;

  typedef CGAL::KDS::Regular_triangulation_3<Traits> KDel;

  Traits tr;
  KDel kdel(tr);


  typedef Traits::Simulator::Time Time;
  typedef CGAL::KDS::Insert_event<Traits::Moving_point_table> MOI;
  Traits::Simulator::Pointer sp= tr.simulator_pointer();

  Traits::Function_kernel::Construct_function cf= tr.function_kernel_object().construct_function_object();
  kdel.set_has_certificates(true);
   
  sp->new_event(Time(0.00000000), MOI(WMP(MP(cf(10),
					     cf(0),  
					     cf(0)), cf(1)), tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000001), MOI(WMP(MP(cf(10), 
					     cf(0),  
					     cf(2)), cf(1)), tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000002), MOI(WMP(MP(cf(12), cf(0), cf(0)), cf(1)),
				      tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000003), MOI(WMP(MP(cf(10), cf(2), cf(0)), cf(1)), 
				      tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000004),MOI( WMP(MP(cf(10.5, .05), 
					     cf(.5), cf(.5)), cf(1)),
				      tr.moving_point_table_pointer()));
  
  sp->new_event(Time(0.00000005), MOI(WMP(MP(cf(0,1,-2), 
					     cf(2,5),  
					     cf(2,5)), cf(1)),
				      tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000006), MOI(WMP(MP(cf(1,3,-3), 
					     cf(3,-3),  
					     cf(3,5)), cf(1)),
				      tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000007),
		MOI( WMP(MP(cf(2,1,-1), cf(5,2),  cf(2,5)), cf(1)), 
		     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000008), 
		MOI(WMP(MP(cf(3,2,3), cf(3,6),  cf(-1,5)), cf(1)), 
		    tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000009), 
		MOI(WMP(MP(cf(5,2,3), cf(-3,6),  cf(4,5)), cf(1)), 
		    tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000010), 
		MOI(WMP(MP(cf(-5), cf(-4),  cf(-6,5)), cf(1)), 
		    tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000011), 
		MOI(WMP(MP(cf(7,2,5), cf(-1,6),  cf(6,5)), cf(1)), 
		    tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000012), 
		MOI(WMP(MP(cf(3.1,-1,7), cf(-4.5,6),  cf(-2,5)), cf(1)), 
		    tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00000013), 
		MOI(WMP(MP(cf(-1,2,-5), cf(2.3,6),  cf(0,5)), cf(1)), 
		    tr.moving_point_table_pointer()));



  sp->set_current_event_number(10000);
  return EXIT_SUCCESS;
};
