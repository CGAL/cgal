#ifndef CGAL_GPS_SEGEMENT_TRAITS_2_APURVA_H
#define CGAL_GPS_SEGEMENT_TRAITS_2_APURVA_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Gps_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <iostream>
//using namespace std;

namespace CGAL {

template <class Kernel_>
class Gps_segement_traits_2_apurva :  public Gps_traits_2<Arr_segment_traits_2<Kernel_> >
{
public:
  Gps_segement_traits_2_apurva<Kernel_>() : Gps_traits_2<Arr_segment_traits_2<Kernel_> >()
  {
    std::cout<<"using Apurva's gps_segement_traits_2 setup"<<std::endl;
  }

};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
