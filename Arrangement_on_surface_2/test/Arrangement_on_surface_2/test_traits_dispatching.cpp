#include <CGAL/basic.h>
#include <iostream>
#include <cassert>

#include <CGAL/Arrangement_2/Arr_traits_adaptor_2_dispatching.h>

#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>

int dispatch(CGAL::Arr_use_dummy_tag) {
  return 0;
}

int dispatch(CGAL::Arr_use_traits_tag) {
  return 1;
}

int main ()
{

  ////////////////
  // left-right //
  ////////////////

  // oblivious-oblivious
  // oblivious-open
  // oblivious-contracted
  // oblivious-closed

  // Arr_left_right_implementation_dispatch oblivious-oblivious
  {
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_oblivious_side_tag, CGAL::Arr_oblivious_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 0);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 0);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 0);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 0);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 0);
  }
  {
    // Arr_left_right_implementation_dispatch oblivious-open
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_oblivious_side_tag, CGAL::Arr_open_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 0);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 0);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 0);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch oblivious-contracted
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_oblivious_side_tag, CGAL::Arr_contracted_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 0);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 0);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch oblivious-closed
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_oblivious_side_tag, CGAL::Arr_closed_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 1);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 1);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }

  // open-oblivious
  // open-open
  // open-contracted
  // open-closed

  // Arr_left_right_implementation_dispatch open-oblivious
  {
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_open_side_tag, CGAL::Arr_oblivious_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 0);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 0);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 0);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch open-open
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_open_side_tag, CGAL::Arr_open_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 0);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 0);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 0);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch open-contracted
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_open_side_tag, CGAL::Arr_contracted_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 0);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 0);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch open-closed
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_open_side_tag, CGAL::Arr_closed_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 1);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 1);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }

  // contracted-oblivious
  // contracted-open
  // contracted-contracted
  // contracted-closed
  // Arr_left_right_implementation_dispatch contracted-oblivious
  {
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_contracted_side_tag, CGAL::Arr_oblivious_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 0);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 0);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch contracted-open
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_contracted_side_tag, CGAL::Arr_open_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 0);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 0);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch contracted-contracted
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_contracted_side_tag, CGAL::Arr_contracted_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 0);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 0);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch contracted-closed
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_contracted_side_tag, CGAL::Arr_closed_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 1);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 1);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }

  // closed-oblivious
  // closed-open
  // closed-contracted
  // closed-closed
  // Arr_left_right_implementation_dispatch closed-oblivious
  {
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_closed_side_tag, CGAL::Arr_oblivious_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 1);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 1);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch closed-open
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_closed_side_tag, CGAL::Arr_open_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 1);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 1);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch closed-contracted
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_closed_side_tag, CGAL::Arr_contracted_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 1);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 1);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_left_right_implementation_dispatch closed-closed
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_closed_side_tag, CGAL::Arr_closed_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 1);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 1);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 0);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 0);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 1);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }


  // identified-identified
  {
    // Arr_left_right_implementation_dispatch identified-identified
    typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    CGAL::Arr_oblivious_side_tag, CGAL::Arr_identified_side_tag > LR;

    typedef LR::Parameter_space_in_x_2_curve_end_tag Psx_2_curve_end;
    assert(dispatch(Psx_2_curve_end()) == 1);
    typedef LR::Parameter_space_in_x_2_curve_tag Psx_2_curve;
    assert(dispatch(Psx_2_curve()) == 0);
    typedef LR::Parameter_space_in_x_2_point_tag Psx_2_point;
    assert(dispatch(Psx_2_point()) == 0);

    typedef LR::Is_on_y_identification_2_curve_tag Ioyi_2_curve;
    assert(dispatch(Ioyi_2_curve()) == 1);
    typedef LR::Is_on_y_identification_2_point_tag Ioyi_2_point;
    assert(dispatch(Ioyi_2_point()) == 1);

    typedef LR::Compare_y_on_boundary_2_points_tag Cmp_y_ob_2_points;
    assert(dispatch(Cmp_y_ob_2_points()) == 1);

    typedef LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_y_nb_2_curve_ends;
    assert(dispatch(Cmp_y_nb_2_curve_ends()) == 1);
  }

  ////////////////
  // bottom-top //
  ////////////////

  // oblivious-oblivious
  // oblivious-open
  // oblivious-contracted
  // oblivious-closed
  {
    // Arr_bottom_top_implementation_dispatch oblivious-oblivious
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_oblivious_side_tag, CGAL::Arr_oblivious_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 0);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 0);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 0);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 0);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 0);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 0);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 0);
  }
  {
    // Arr_bottom_top_implementation_dispatch oblivious-open
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_oblivious_side_tag, CGAL::Arr_open_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 0);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 0);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 0);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 0);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 0);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 0);
  }
  {
    // Arr_bottom_top_implementation_dispatch oblivious-contracted
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_oblivious_side_tag, CGAL::Arr_contracted_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 0);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 0);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 1);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 1);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_bottom_top_implementation_dispatch oblivious-closed
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_oblivious_side_tag, CGAL::Arr_closed_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 1);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 1);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 0);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 0);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }

  // open-oblivious
  // open-open
  // open-contracted
  // open-closed
  {
    // Arr_bottom_top_implementation_dispatch open-oblivious
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_open_side_tag, CGAL::Arr_oblivious_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 0);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 0);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 0);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 0);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 0);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 0);
  }
  {
    // Arr_bottom_top_implementation_dispatch open-open
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_open_side_tag, CGAL::Arr_open_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 0);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 0);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 0);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 0);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 0);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 0);
  }
  {
    // Arr_bottom_top_implementation_dispatch open-contracted
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_open_side_tag, CGAL::Arr_contracted_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 0);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 0);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 1);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 1);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_bottom_top_implementation_dispatch open-closed
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_open_side_tag, CGAL::Arr_closed_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 1);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 1);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 0);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 0);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }

  // contracted-oblivious
  // contracted-open
  // contracted-contracted
  // contracted-closed
  {
    // Arr_bottom_top_implementation_dispatch contracted-oblivious
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_contracted_side_tag, CGAL::Arr_oblivious_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 0);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 0);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 1);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 1);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_bottom_top_implementation_dispatch contracted-open
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_contracted_side_tag, CGAL::Arr_open_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 0);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 0);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 1);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 1);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_bottom_top_implementation_dispatch contracted-contracted
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_contracted_side_tag, CGAL::Arr_contracted_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 0);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 0);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 1);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 1);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_bottom_top_implementation_dispatch contracted-closed
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_contracted_side_tag, CGAL::Arr_closed_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 1);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 1);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 1);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 1);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }

  // closed-oblivious
  // closed-open
  // closed-contracted
  // closed-closed
  {
    // Arr_bottom_top_implementation_dispatch closed-oblivious
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_closed_side_tag, CGAL::Arr_oblivious_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 1);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 1);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 0);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 0);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_bottom_top_implementation_dispatch closed-open
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_closed_side_tag, CGAL::Arr_open_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 1);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 1);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 0);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 0);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_bottom_top_implementation_dispatch closed-contracted
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_closed_side_tag, CGAL::Arr_contracted_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 1);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 1);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 1);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 1);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }
  {
    // Arr_bottom_top_implementation_dispatch closed-closed
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_closed_side_tag, CGAL::Arr_closed_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 1);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 1);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 0);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 0);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 1);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 0);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 0);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }

  // identified-identified
  {
    // Arr_bottom_top_implementation_dispatch identified-identified
    typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    CGAL::Arr_identified_side_tag, CGAL::Arr_identified_side_tag > BT;

    typedef BT::Parameter_space_in_y_2_curve_end_tag Psy_2_curve_end;
    assert(dispatch(Psy_2_curve_end()) == 1);
    typedef BT::Parameter_space_in_y_2_curve_tag Psy_2_curve;
    assert(dispatch(Psy_2_curve()) == 0);
    typedef BT::Parameter_space_in_y_2_point_tag Psy_2_point;
    assert(dispatch(Psy_2_point()) == 0);

    typedef BT::Is_on_x_identification_2_curve_tag Ioxi_2_curve;
    assert(dispatch(Ioxi_2_curve()) == 1);
    typedef BT::Is_on_x_identification_2_point_tag Ioxi_2_point;
    assert(dispatch(Ioxi_2_point()) == 1);

    typedef BT::Compare_x_on_boundary_2_points_tag Cmp_x_ob_2_points;
    assert(dispatch(Cmp_x_ob_2_points()) == 1);
    typedef BT::Compare_x_on_boundary_2_point_curve_end_tag Cmp_x_ob_2_point_curve_end;
    assert(dispatch(Cmp_x_ob_2_point_curve_end()) == 0);
    typedef BT::Compare_x_on_boundary_2_curve_ends_tag Cmp_x_ob_2_curve_ends;
    assert(dispatch(Cmp_x_ob_2_curve_ends()) == 0);
    typedef BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_x_nb_2_curve_ends;
    assert(dispatch(Cmp_x_nb_2_curve_ends()) == 1);
  }
  return EXIT_SUCCESS;

}
