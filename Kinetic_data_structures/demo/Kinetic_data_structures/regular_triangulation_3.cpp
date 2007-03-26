#include <CGAL/Kinetic/Regular_triangulation_3.h>
#include <CGAL/Kinetic/Regular_triangulation_exact_simulation_traits.h>
#include <fstream>
#ifdef CGAL_USE_COIN
#include "include/SoQt_moving_points_3.h"
#include "include/SoQt_triangulation_3.h"
#include "include/SoQt_widget_3.h"
#include <CGAL/Kinetic/Insert_event.h>
#endif


int main(int argc, char *argv[])
{
#ifdef CGAL_USE_COIN
    typedef CGAL::Kinetic::Regular_triangulation_exact_simulation_traits Traits;
    typedef CGAL::Kinetic::SoQt_widget_3<Traits::Simulator> Qt_gui;
    typedef CGAL::Kinetic::SoQt_moving_points_3<Traits, Qt_gui> Qt_mpt;
    typedef Traits::Kinetic_kernel::Motion_function MF;
    typedef Traits::Kinetic_kernel::Weighted_point_3 MP;
    typedef Traits::Kinetic_kernel::Point_3 MPP;

    Traits tr(0,100000);

    CGAL_KINETIC_SET_LOG_LEVEL(CGAL::Kinetic::LOG_LOTS);

    Qt_gui::Handle qtsim= new Qt_gui(argc, argv, tr.simulator_handle());
    Qt_mpt::Handle qtmpt= new Qt_mpt(tr, qtsim);

    typedef CGAL::Kinetic::Regular_triangulation_3<Traits> KDel;
    typedef CGAL::Kinetic::SoQt_triangulation_3<KDel, Qt_gui, Qt_mpt> CoinKDel;
    KDel::Handle kdel= new KDel(tr);
    CoinKDel::Handle cd= new CoinKDel(kdel, qtsim, qtmpt);

    if (argc==1) {

      typedef Traits::Simulator::Time Time;
      typedef CGAL::Kinetic::Insert_event<Traits::Active_points_3_table> MOI;

      Traits::Kinetic_kernel::Function_kernel::Construct_function cf= tr.kinetic_kernel_object().function_kernel_object().construct_function_object();

      tr.simulator_handle()->new_event(Time(0.000000), MOI(MP(MPP(cf(10, 1),
								  cf(0, .1),
								 cf(0)), cf(.1)) ,
							   tr.active_points_3_table_handle()));
      tr.simulator_handle()->new_event(Time(0.000001), MOI(MP(MPP(cf(10),
								  cf(0),
								 cf(2)), cf(.1)) ,
							   tr.active_points_3_table_handle()));
      tr.simulator_handle()->new_event(Time(0.000002), MOI(MP(MPP(cf(12),
								  cf(0),
								  cf(0)), cf(.1)) ,
							  tr.active_points_3_table_handle()));
      tr.simulator_handle()->new_event(Time(0.000003), MOI(MP(MPP(cf(10),
								  cf(2),
								  cf(0)), cf(.1)),
							   tr.active_points_3_table_handle()));
      tr.simulator_handle()->new_event(Time(0.000000035), MOI(MP(MPP(cf(1), cf(10),  cf(0)), cf(.1)) ,
							      tr.active_points_3_table_handle()));
      tr.simulator_handle()->new_event(Time(0.000004),MOI(MP(MPP(cf(10.5, .05),
								 cf(.5),
								 cf(.5)), cf(.1)),
							  tr.active_points_3_table_handle()));

      tr.simulator_handle()->set_current_event_number(4);
    } else {
      std::cout << "Reading from " << argv[1] << std::endl;
      std::ifstream in(argv[1]);
      in >> *tr.active_points_3_table_handle();
    }
    kdel->set_has_certificates(true);


    std::cout << "This program displays a 3D kinetic Delaunay triangulation.\n";
    std::cout << "What is displayed can be controlled by pushing the following keys over the Coin window and with the arrow tool selected.\n";
    std::cout << "Press 'h' to hide/show the convex hull.\n";
    std::cout << "Press 'f' to hide/show faces.\n";
    std::cout << "Press 's' to show spheres for points, 'p' to show points.\n";

    return qtsim->begin_event_loop();
 #else
    std::cout << "An install of Inventor and SoQt are required for this demo.  "
      "Please make sure they are installed and then compile "
      "using the makefile 'makefile.soqt'.\n"
      "They can be found at http://www.coin3d.org or as an rpm from "
      "your linux distribution (they are part of Fedora extras, for example).\n";
    return EXIT_FAILURE;
    #endif
};
