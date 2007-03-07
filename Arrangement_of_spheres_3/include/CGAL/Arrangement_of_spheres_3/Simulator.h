#ifndef CGAL_ARR_OF_S_3_H_SIM_TRIATS
#define CGAL_ARR_OF_S_3_H_SIM_TRIATS


#include <CGAL/Arrangement_of_spheres_3/Function_kernel.h>
#include <CGAL/Kinetic/IO/Qt_widget_2.h>
#include <CGAL/Kinetic/Two_list_pointer_event_queue.h>
#include <CGAL/Kinetic/Default_simulator.h>

typedef CGAL::Kinetic::Default_simulator<Function_kernel,
					 CGAL::Kinetic::Two_list_pointer_event_queue<Function_kernel, true> > Simulator;
typedef CGAL::Kinetic::Qt_widget_2<Simulator> Qt_gui;
#endif
