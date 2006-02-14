#include <CGAL/Kinetic/Regular_triangulation_3.h>
#include <CGAL/Kinetic/Regular_triangulation_exact_simulation_traits_3.h>

#ifdef CGAL_USE_SOQT
#include "include/SoQt_moving_points_3.h"
#include "include/SoQt_triangulation_3.h"
#include "include/SoQt_widget_3.h"
#include <CGAL/Kinetic/Insert_event.h>
#endif


int main(int argc, char *argv[])
{
#ifdef CGAL_USE_SOQT
    typedef CGAL::Kinetic::Regular_triangulation_exact_simulation_traits_3 Traits;
    typedef CGAL::Kinetic::SoQt_widget_3<Traits::Simulator> Qt_gui;
    typedef CGAL::Kinetic::SoQt_moving_points_3<Traits, Qt_gui> Qt_mpt;
    typedef Traits::Kinetic_kernel::Motion_function MF;
    typedef Traits::Kinetic_kernel::Weighted_point_3 MP;
    typedef Traits::Kinetic_kernel::Point_3 MPP;

    Traits tr;

    CGAL_KINETIC_SET_LOG_LEVEL(CGAL::Kinetic::LOG_LOTS);

    Qt_gui::Pointer qtsim= new Qt_gui(argc, argv, tr.simulator_pointer());
    Qt_mpt::Pointer qtmpt= new Qt_mpt(tr, qtsim);

    typedef CGAL::Kinetic::Regular_triangulation_3<Traits> KDel;
    typedef CGAL::Kinetic::SoQt_triangulation_3<KDel, Qt_gui, Qt_mpt> CoinKDel;
    KDel::Pointer kdel= new KDel(tr);
    CoinKDel::Pointer cd= new CoinKDel(kdel, qtsim, qtmpt);

    typedef Traits::Simulator::Time Time;
    typedef CGAL::Kinetic::Insert_event<Traits::Active_objects_table> MOI;

    Traits::Kinetic_kernel::Function_kernel::Construct_function cf= tr.kinetic_kernel_object().function_kernel_object().construct_function_object();

    tr.simulator_pointer()->new_event(Time(0.000000), MOI(MP(MPP(cf(10, 1),
								 cf(0, .1),
								 cf(0)), cf(.1)) , 
							  tr.active_objects_table_pointer()));
    tr.simulator_pointer()->new_event(Time(0.000001), MOI(MP(MPP(cf(10),
								 cf(0),
								 cf(2)), cf(.1)) ,
							  tr.active_objects_table_pointer()));
    tr.simulator_pointer()->new_event(Time(0.000002), MOI(MP(MPP(cf(12),
								 cf(0),
								 cf(0)), cf(.1)) ,
							  tr.active_objects_table_pointer()));
    tr.simulator_pointer()->new_event(Time(0.000003), MOI(MP(MPP(cf(10),
								 cf(2),
								 cf(0)), cf(.1)),
							  tr.active_objects_table_pointer()));
    tr.simulator_pointer()->new_event(Time(0.000000035), MOI(MP(MPP(cf(1), cf(10),  cf(0)), cf(.1)) , 
							     tr.active_objects_table_pointer()));
    tr.simulator_pointer()->new_event(Time(0.000004),MOI(MP(MPP(cf(10.5, .05),
								cf(.5),
								cf(.5)), cf(.1)),
							 tr.active_objects_table_pointer()));

    CGAL_KINETIC_SET_LOG_LEVEL(CGAL::Kinetic::LOG_LOTS);
    tr.simulator_pointer()->set_current_event_number(4);
    kdel->set_has_certificates(true);

    return qtsim->begin_event_loop();
#else
    bool regular_triangulation_3_compiled_without_CGAL_USE_SOQT_defined;
    return EXIT_FAILURE;
#endif
};
