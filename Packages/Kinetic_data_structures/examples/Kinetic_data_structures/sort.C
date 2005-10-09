#include <CGAL/KDS/basic.h>

#include <CGAL/KDS/Exact_simulation_traits_1.h>
#include <CGAL/KDS/Insert_event.h>
#include <CGAL/KDS/Sort.h>




int main(int, char *[]){

  typedef CGAL::KDS::Exact_simulation_traits_1 Traits;

  typedef CGAL::KDS::Insert_event<Traits::Moving_point_table> Insert_event;
  typedef Traits::Moving_point_table::Data Moving_point;
  typedef CGAL::KDS::Sort<Traits> Sort;
  typedef Traits::Simulator::Time Time;

  Traits tr;
  Sort sort(tr);
  Traits::Simulator::Pointer sp= tr.simulator_pointer();

  Traits::Function_kernel::Construct_function fc= tr.function_kernel_object().construct_function_object();


  sp->new_event(Time(0.00000), 
		Insert_event(Moving_point(fc(0,3.5)), 
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00001),
		Insert_event(Moving_point(fc(3,0)),
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00002), 
		Insert_event(Moving_point(fc(2,1,-1)), 
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00003), 
		Insert_event(Moving_point(fc(1,3,-3)), 
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00004), 
		Insert_event(Moving_point(fc(2.3,-1)),
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00005), 
		Insert_event(Moving_point(fc(6,3,-3)),
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00006), 
		Insert_event(Moving_point(fc(-9,3,-3)),
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00007), 
		Insert_event(Moving_point(fc(11,3,-3)),
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00008), 
		Insert_event(Moving_point(fc(5,3,-3)),
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00009), 
		Insert_event(Moving_point(fc(-2,3,-3, 9)), 
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00010), 
		Insert_event(Moving_point(fc(-6,3,-13, 94)),
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00011), 
		Insert_event(Moving_point(fc(-7,32,-3, 8)),
			     tr.moving_point_table_pointer()));
  sp->new_event(Time(0.00012), 
		Insert_event(Moving_point(fc(-9,1,-12)),
			     tr.moving_point_table_pointer()));
  
  while (sp->next_event_time() != sp->end_time()){
    sp->set_current_event_number(sp->current_event_number()+1);
  }

  return EXIT_SUCCESS;
};
