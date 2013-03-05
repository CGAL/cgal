// Copyright (c) 2009 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent Rineau
// 

template <class R>
bool
_test_compare_dihedral_angle_3(const R& rep)
{
  typedef typename R::Point_3 Point_3;
  typedef typename R::Vector_3 Vector_3;
  typedef typename R::FT FT;
  typename R::Compare_dihedral_angle_3 compare_dih_angle
    = rep.compare_dihedral_angle_3_object();

  for(int theta1 = -170; theta1 <= 180; theta1+= 10) {
    const double angle1 = CGAL_PI*theta1/180.;
    Point_3 p1(0, 0, -1);
    Point_3 p2(0, 0, 1);
    Point_3 p3(1, 0, 0);
    Point_3 p4((int)(std::cos(angle1)*1000), (int)(std::sin(angle1)*1000), 0);
    if(theta1 == 180) {
      p4 = Point_3(-1, 0, 0);
    }
    for(int theta2 = -170; theta2 <= 180; theta2+= 10) {
      // if(theta1 == theta2) continue;
      const double angle2 = CGAL_PI*theta2/180.;
      Point_3 q1(0, 1, 0);
      Point_3 q2(0, -1, 0);
      Point_3 q3(0, 0, 1);
      Point_3 q4((int)(std::sin(angle2)*1000), 0, (int)(std::cos(angle2)*1000));
      if(theta2 == 180) {
        q4 = Point_3(0, 0, -1);
      }
      const FT approx_cosine2 = FT((int)(std::cos(angle2)*1000))/FT(1000);
      const typename R::Comparison_result 
        comp_result = compare_dih_angle(p1, p2, p3, p4, q1, q2, q3, q4),
        global_fct_call = CGAL::compare_dihedral_angle(p1, p2, p3, p4,
                                                       q1, q2, q3, q4),
        call_with_vectors = compare_dih_angle(p2 - p1,
                                              p3 - p1,
                                              p4 - p1,
                                              q2 - q1,
                                              q3 - q1,
                                              q4 - q1),
        call_fct_with_vectors = CGAL::compare_dihedral_angle(p2 - p1,
                                                             p3 - p1,
                                                             p4 - p1,
                                                             q2 - q1,
                                                             q3 - q1,
                                                             q4 - q1),
        theorical_result = CGAL::compare(std::abs(theta1),std::abs(theta2)),
        compare_with_cosine = compare_dih_angle(p1, p2, p3, p4,
                                                approx_cosine2),
        compare_with_cosine_with_vectors = compare_dih_angle(p2 - p1,
                                                             p3 - p1,
                                                             p4 - p1,
                                                             approx_cosine2),
        compare_with_cosine_fct = CGAL::compare_dihedral_angle(p1, p2, p3, p4,
                                                               approx_cosine2),
        compare_with_cosine_with_vectors_fct = 
        CGAL::compare_dihedral_angle(p2 - p1,
                                     p3 - p1,
                                     p4 - p1,
                                     approx_cosine2);
      if(comp_result != theorical_result || 
         comp_result != global_fct_call ||
         comp_result != call_with_vectors ||
         comp_result != call_fct_with_vectors ||
         ( (theta1 != theta2) && (theta1 != -theta2)
           && (
               comp_result != compare_with_cosine || 
               comp_result != compare_with_cosine_fct ||
               comp_result != compare_with_cosine_with_vectors ||
               comp_result != compare_with_cosine_with_vectors_fct
               ) 
           )
         )
      {
        std::cerr << "Error compare_dihedral_angle_3, with angles "
                  << theta1 << " and " << theta2 << std::endl;
        std::cerr << "Results are: "
                  << comp_result << " " << global_fct_call << " "
                  << call_with_vectors << " " << call_fct_with_vectors << " "
                  << compare_with_cosine << " " << compare_with_cosine_fct << " "
                  << compare_with_cosine_with_vectors << " " 
                  << compare_with_cosine_with_vectors_fct << "\n";
        const Vector_3 u1 = p2 - p1;
        const Vector_3 v1 = p3 - p1;
        const Vector_3 w1 = p4 - p1;
        const Vector_3 uv1 = cross_product(u1, v1);
        const Vector_3 uw1 = cross_product(u1, w1);
        const Vector_3 u2 = q2 - q1;
        const Vector_3 v2 = q3 - q1;
        const Vector_3 w2 = q4 - q1;
        const Vector_3 uv2 = cross_product(u2, v2);
        const Vector_3 uw2 = cross_product(u2, w2);
        std::cerr << "Squared cosine and signs are: \n"
                  << CGAL::to_double((uv1 * uw1) * (uv1 * uw1) / 
                                     ((uv1*uv1) * ( uw1*uw1))) << " " 
                  << CGAL::sign(uv1*uw1) << "\n"
                  << CGAL::to_double((uv2 * uw2) * (uv2 * uw2) /
                                     ((uv2*uv2) * ( uw2*uw2))) << " " 
                  << CGAL::sign(uv2*uw2) << "\n";
        return false;
      }
    } // end loop on theta2
    // if(CGAL::compare_dihedral_angle(p1, p2, p3, p4, FT(-1)/FT(2)) != 
    //    CGAL::compare(std::abs(theta1), 120) ||
    //    CGAL::compare_dihedral_angle(p1, p2, p3, p4, FT(1)/FT(2)) != 
    //    CGAL::compare(std::abs(theta1), 60))
    // {
    //   std::cerr << "Error compare_dihedral_angle_3, with angle "
    //             << theta1 << " and cosine" << std::endl;
    //   return false;
    // }
  } // end loop and theta1
  return true;
}

