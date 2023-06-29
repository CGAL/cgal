// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) :

// Partially supported by the IST Programme of the EU as a
// STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#include <set>

template <class SK>
inline
void test_make_monotone_an_already_monotone_arc(const typename SK::Circular_arc_3& arc,
                                                const typename SK::Sphere_3& ref_sphere)
{
  typename SK::Is_theta_monotone_3 is_t_mon=SK().is_theta_monotone_3_object(ref_sphere);
  assert (is_t_mon(arc));
  CGAL::is_theta_monotone(arc,ref_sphere);
  typename SK::Make_theta_monotone_3 mk_mono=SK().make_theta_monotone_3_object(ref_sphere);
  std::vector<typename SK::Circular_arc_3> arcs;
  mk_mono(arc,std::back_inserter(arcs));
  assert(arcs.size()==1);
  assert (arcs[0]==arc);
}

template <class SK>
inline
void test_non_monotone_arcs_decomposition(const typename SK::Circular_arc_3& arc,
                                          const typename SK::Sphere_3& ref_sphere,
                                          const std::vector< typename SK::Circular_arc_3>& expected)
{
  typename SK::Is_theta_monotone_3 is_t_mon=SK().is_theta_monotone_3_object(ref_sphere);
  assert( !is_t_mon(arc) );
  assert( !is_theta_monotone(arc,ref_sphere) );

  std::vector<typename SK::Circular_arc_3> res;
  typename SK::Make_theta_monotone_3 mk_mono=SK().make_theta_monotone_3_object(ref_sphere);
  mk_mono(arc,std::back_inserter(res));
  assert(expected.size()==res.size());
  std::set<unsigned> matches;
  for (typename std::vector<typename SK::Circular_arc_3>::iterator it=res.begin();it!=res.end();++it){
    assert(is_theta_monotone(*it,ref_sphere));
    bool found=false;
    for (unsigned i=0;i<expected.size();++i){
      if ( matches.find(i)!=matches.end() ) continue;
      if (*it==expected[i]){
        found=true;
        matches.insert(i);
      }
    }
    assert(found);
  }
}


//FUNCTIONS ONLY VALID FOR CIRCLES DEFINED BELOW
template <class SK>
void test_normal_circle_monotonicity(const typename SK::Circle_3& circle,
                                     const typename SK::Sphere_3& ref_sphere)
{
  typename SK::Is_theta_monotone_3 is_t_mon=SK().is_theta_monotone_3_object(ref_sphere);
  std::vector<std::variant<typename SK::Circle_3, std::pair<typename SK::Circular_arc_point_3, unsigned>>> vect_obj;
  typename SK::FT zcoord = CGAL::SphericalFunctors::extremal_points_z_coordinate<SK>(circle,ref_sphere);
  //create extremal points of circle
  SK().intersect_3_object()(circle,typename SK::Plane_3(0,0,1,-zcoord),std::back_inserter(vect_obj));
  assert(vect_obj.size()==2);
  typename SK::Circular_arc_point_3 extrems[2];
  extrems[0]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&vect_obj[0])->first;
  extrems[1]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&vect_obj[1])->first;
  //create non extremal points on circle
  vect_obj.clear();
  SK().intersect_3_object()(circle,typename SK::Plane_3(0,0,1,-zcoord-typename SK::FT(0.1)),std::back_inserter(vect_obj));
  assert(vect_obj.size()==2);
  typename SK::Circular_arc_point_3 other_pts[2];
  other_pts[0]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&vect_obj[0])->first;
  other_pts[1]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&vect_obj[1])->first;

  //Test Make_theta_monotone+[Ii]s_theta_monotone(_3) on monotone arcs
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,extrems[0],extrems[1]),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,extrems[1],extrems[0]),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,other_pts[0],other_pts[1]),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,extrems[0],other_pts[1]),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,extrems[0],other_pts[0]),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,other_pts[1],extrems[1]),ref_sphere);


  std::vector<typename SK::Circular_arc_3> expected;
  expected.push_back(typename SK::Circular_arc_3(circle,other_pts[1],extrems[1]));
  expected.push_back(typename SK::Circular_arc_3(circle,extrems[1],extrems[0]));
  expected.push_back(typename SK::Circular_arc_3(circle,extrems[0],other_pts[0]));
  test_non_monotone_arcs_decomposition<SK>( typename SK::Circular_arc_3(circle,other_pts[1],other_pts[0]),
                                        ref_sphere, expected
  );

  expected.clear();
  expected.push_back(typename SK::Circular_arc_3(circle,other_pts[1],extrems[1]));
  expected.push_back(typename SK::Circular_arc_3(circle,extrems[1],extrems[0]));
  test_non_monotone_arcs_decomposition<SK>( typename SK::Circular_arc_3(circle,other_pts[1],extrems[0]),
                                        ref_sphere, expected
  );

  //Make circle montone
  std::vector<typename SK::Circular_arc_3> res;
  typename SK::Make_theta_monotone_3 mk_mono=SK().make_theta_monotone_3_object(ref_sphere);
  mk_mono(circle,std::back_inserter(res));
  assert(res.size()==2);
  assert(res[0].source()==res[1].target() && res[0].target()==res[1].source());
  assert(is_theta_monotone(res[0],ref_sphere));
  assert(is_theta_monotone(res[1],ref_sphere));
}


template <class SK>
void test_threaded_circle_monotonicity(const typename SK::Circle_3& circle,
                                       const typename SK::Sphere_3& ref_sphere)
{
  typename SK::Is_theta_monotone_3 is_t_mon=SK().is_theta_monotone_3_object(ref_sphere);
  std::vector<std::variant<typename SK::Circle_3, std::pair<typename SK::Circular_arc_point_3,unsigned>>> vect_obj;
  SK().intersect_3_object()(circle,typename SK::Plane_3(1,0,0,-ref_sphere.center().x()),std::back_inserter(vect_obj));
  assert(vect_obj.size()==2);
  typename SK::Circular_arc_point_3 pts1[2];
  pts1[0]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&vect_obj[0])->first;
  pts1[1]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&vect_obj[1])->first;
  vect_obj.clear();
  SK().intersect_3_object()(circle,typename SK::Plane_3(1,1,0,-ref_sphere.center().x()-ref_sphere.center().y()),std::back_inserter(vect_obj));
  assert(vect_obj.size()==2);
  typename SK::Circular_arc_point_3 pts2[2];
  pts2[0]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&vect_obj[0])->first;
  pts2[1]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&vect_obj[1])->first;
  //assertions
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,pts1[0],pts1[1]),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,pts2[0],pts2[1]),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,pts1[0],pts2[1]),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,pts1[1],pts2[0]),ref_sphere);

  //Make circle montone
  std::vector<typename SK::Circular_arc_3> res;
  typename SK::Make_theta_monotone_3 mk_mono=SK().make_theta_monotone_3_object(ref_sphere);
  mk_mono(circle,std::back_inserter(res));
  assert(res.size()==1);
  assert(res[0].source()==res[0].target());
  assert(is_theta_monotone(res[0],ref_sphere));
}

template<class SK>
void test_polar_circle_monotonicity(const typename SK::Circle_3& circle,
                                    const typename SK::Sphere_3& ref_sphere)
{
  typename SK::Is_theta_monotone_3 is_t_mon=SK().is_theta_monotone_3_object(ref_sphere);
  std::vector<std::variant<typename SK::Circle_3, std::pair<typename SK::Circular_arc_point_3,unsigned>>> vect_obj;
  SK().intersect_3_object()(circle,typename SK::Plane_3(0,0,1,-circle.center().z()),std::back_inserter(vect_obj));
  assert(vect_obj.size()==2);
  typename SK::Circular_arc_point_3 pts[2];
  pts[0]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&vect_obj[0])->first;
  pts[1]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&vect_obj[1])->first;



  typename SK::Root_of_2 radius=CGAL::make_root_of_2(typename SK::FT(0),typename SK::FT(1),ref_sphere.squared_radius());
  bool npole=circle.center().z() > ref_sphere.center().z();
  typename SK::Circular_arc_point_3 pole=
  typename SK::Algebraic_kernel::Root_for_spheres_2_3( ref_sphere.center().x(),
                                                       ref_sphere.center().y(),
                                                       ref_sphere.center().z()+(npole?1:-1)*radius
  );


  if (npole){
    test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,pts[1],pts[0]),ref_sphere);
    std::vector<typename SK::Circular_arc_3> expected;
    expected.push_back(typename SK::Circular_arc_3(circle,pts[0],pole));
    expected.push_back(typename SK::Circular_arc_3(circle,pole,pts[1]));
    test_non_monotone_arcs_decomposition<SK>( typename SK::Circular_arc_3(circle,pts[0],pts[1]),
                                          ref_sphere, expected
    );
  }
  else{
    test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,pts[0],pts[1]),ref_sphere);
    std::vector<typename SK::Circular_arc_3> expected;
    expected.push_back(typename SK::Circular_arc_3(circle,pts[1],pole));
    expected.push_back(typename SK::Circular_arc_3(circle,pole,pts[0]));
    test_non_monotone_arcs_decomposition<SK>( typename SK::Circular_arc_3(circle,pts[1],pts[0]),
                                          ref_sphere, expected
    );
  }

  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,pole,pts[1]),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,pole,pts[0]),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,pts[0],pole),ref_sphere);
  test_make_monotone_an_already_monotone_arc<SK>(typename SK::Circular_arc_3(circle,pts[1],pole),ref_sphere);

  //Make circle montone
  std::vector<typename SK::Circular_arc_3> res;
  typename SK::Make_theta_monotone_3 mk_mono=SK().make_theta_monotone_3_object(ref_sphere);
  mk_mono(circle,std::back_inserter(res));
  assert(res.size()==1);
  assert(res[0].source()==res[0].target());
  assert(res[0].source()==pole);
  assert(is_theta_monotone(res[0],ref_sphere));
}

template <class SK>
inline
typename SK::Vector_3 get_vector_in_plane(unsigned i){
  double Pi_16=CGAL_PI / 16.;
  int s=i>15?-1:1;
  switch (i){
    case  0:
    case 16:
      return typename SK::Vector_3(s*1,0,0); //theta=0 and Pi
    case  4:
    case 20:
      return typename SK::Vector_3(s*1,s*1,0); //theta=Pi/4 and 5/4 Pi
    case  8:
    case 24:
      return typename SK::Vector_3(0,s*1,0); //theta=Pi/2 and 3/2 Pi
    case 12:
    case 28:
      return typename SK::Vector_3(-1*s,s*1,0); //theta=3/4 Pi and 7/4 Pi
    default:
      return typename SK::Vector_3(cos(i*Pi_16),sin(i*Pi_16),0);
  }
}


template <class SK>
inline
typename SK::Plane_3 get_meridians(unsigned i,const typename SK::Point_3& center){
  typedef typename SK::FT FT;
  double Pi_16=CGAL_PI / 16.;
  switch (i){
    case 0:
      return typename SK::Plane_3(0,1,0,-center.y()); //theta=0 and Pi
    case 4:
      return typename SK::Plane_3(1,-1,0,-center.x()+center.y()); //theta=Pi/4 and 5/4 Pi
    case 8:
      return typename SK::Plane_3(1,0,0,-center.x()); //theta=Pi/2 and 3/2 Pi
    case 12:
      return typename SK::Plane_3(1,1,0,-center.x()-center.y()); //theta=3/4 Pi and 7/4 Pi
    default:
      return typename SK::Plane_3(sin(i*Pi_16),-cos(i*Pi_16),0,-FT(sin(i*Pi_16))*center.x()+FT(cos(i*Pi_16))*center.y());
  }
}

unsigned get_num(unsigned plane_num,unsigned int_num){
  assert(int_num<2);
  if (plane_num<9)
    return plane_num + (int_num==0?16:0);
  return plane_num + (int_num==0?0:16);
}

template <class SK>
void fill_intersections(const typename SK::Circle_3& circle,
                        const typename SK::Sphere_3& ref,
                        typename SK::Circular_arc_point_3* inters,unsigned i)
{
  typename SK::Plane_3 meridians=get_meridians<SK>(i,ref.center());
  std::vector<std::variant<typename SK::Circle_3, std::pair<typename SK::Circular_arc_point_3,unsigned>>> objs;
  SK().intersect_3_object()(circle,meridians,std::back_inserter(objs));
  assert(objs.size()==2);
  inters[get_num(i,0)]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[0])->first;
  inters[get_num(i,1)]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;
}


template <class SK>
std::pair<typename SK::Vector_3,typename SK::Vector_3>
get_bounding_vectors(const typename SK::Circular_arc_point_3& point,
                     const typename SK::Sphere_3& ref_sphere)
{
  CGAL::Interval_nt<false> int_x=CGAL::to_interval(point.x());
  CGAL::Interval_nt<false> int_y=CGAL::to_interval(point.y());
  typename SK::Vector_3 all[4];

  all[0]=typename SK::Vector_3(int_x.inf(),int_y.inf(),0)-typename SK::Vector_3(ref_sphere.center().x(),ref_sphere.center().y(),0);
  all[1]=typename SK::Vector_3(int_x.inf(),int_y.sup(),0)-typename SK::Vector_3(ref_sphere.center().x(),ref_sphere.center().y(),0);
  all[2]=typename SK::Vector_3(int_x.sup(),int_y.inf(),0)-typename SK::Vector_3(ref_sphere.center().x(),ref_sphere.center().y(),0);
  all[3]=typename SK::Vector_3(int_x.sup(),int_y.sup(),0)-typename SK::Vector_3(ref_sphere.center().x(),ref_sphere.center().y(),0);
  int smaller=-1;
  int greater=-1;


  for (int i=0;i!=4;++i){
      CGAL::Comparison_result res=CGAL::compare_theta(all[i],point,ref_sphere);
    if(res==CGAL::EQUAL){
      std::cout << "here" << std::endl;
      continue;
    }
    if (res==CGAL::SMALLER)
      smaller=i;
    else
      greater=i;
  }
  assert(smaller!=-1 && greater!=-1);
  return std::make_pair(all[smaller],all[greater]);
}


template <class SK>
void test_extremal_points(const typename SK::Circle_3& circle,
                          const typename SK::Sphere_3& ref_sphere,
                          bool cut_by_M0)
{
  typename SK::Circular_arc_point_3 xtrms[2];
  CGAL::theta_extremal_points(circle,ref_sphere,xtrms);
  assert( xtrms[0]==CGAL::theta_extremal_point(circle,ref_sphere,true) );
  assert( xtrms[1]==CGAL::theta_extremal_point(circle,ref_sphere,false) );
  const typename SK::Circular_arc_point_3& inter1=cut_by_M0?xtrms[1]:xtrms[0];
  const typename SK::Circular_arc_point_3& inter2=cut_by_M0?xtrms[0]:xtrms[1];

  typename SK::Intersect_3 func=SK().intersect_3_object();
  std::vector<std::variant<typename SK::Circle_3, std::pair<typename SK::Circular_arc_point_3,unsigned>>> intersections;

  std::pair<typename SK::Vector_3,typename SK::Vector_3> vect_pair=get_bounding_vectors<SK>(inter1,ref_sphere);
  func(circle,typename SK::Plane_3(ref_sphere.center(),ref_sphere.center()+typename SK::Vector_3(0,0,1),ref_sphere.center()+vect_pair.first),std::back_inserter(intersections));
  assert(intersections.size()==0);
  func(circle,typename SK::Plane_3(ref_sphere.center(),ref_sphere.center()+typename SK::Vector_3(0,0,1),ref_sphere.center()+vect_pair.second),std::back_inserter(intersections));
  assert(intersections.size()==2);

  vect_pair=get_bounding_vectors<SK>(inter2,ref_sphere);
  func(circle,typename SK::Plane_3(ref_sphere.center(),ref_sphere.center()+typename SK::Vector_3(0,0,1),ref_sphere.center()+vect_pair.first),std::back_inserter(intersections));
  assert(intersections.size()==4);
  func(circle,typename SK::Plane_3(ref_sphere.center(),ref_sphere.center()+typename SK::Vector_3(0,0,1),ref_sphere.center()+vect_pair.second),std::back_inserter(intersections));
  assert(intersections.size()==4);
}

template <class SK>
void
test_functionalities_on_a_reference_sphere(const typename SK::Point_3& ref_sphere_center){
  typedef typename SK::FT FT;
  std::cout << "test functionalities on a sphere" << std::endl;
  //=============DATA=========================================================================
  typename SK::Sphere_3 ref_sphere(ref_sphere_center,1);
  typename SK::Root_of_2 radius=CGAL::make_root_of_2(typename SK::FT(0),typename SK::FT(1),ref_sphere.squared_radius());
  typename SK::Circular_arc_point_3 north_pole(
    typename SK::Algebraic_kernel::Root_for_spheres_2_3( ref_sphere.center().x(),
                                                         ref_sphere.center().y(),
                                                         ref_sphere.center().z()+radius)
  );
  typename SK::Circular_arc_point_3 south_pole(
    typename SK::Algebraic_kernel::Root_for_spheres_2_3( ref_sphere.center().x(),
                                                         ref_sphere.center().y(),
                                                         ref_sphere.center().z()-radius)
  );


  typename SK::Circle_3 great_threaded(ref_sphere,typename SK::Plane_3(0,0,1  ,-ref_sphere_center.z() ));
  typename SK::Circle_3 threaded      (ref_sphere,typename SK::Plane_3(0,0,1  ,FT(0.5)-ref_sphere_center.z() ));
  typename SK::Circle_3 bipolar       (ref_sphere,typename SK::Plane_3(0,1,0  ,-ref_sphere_center.y()));
  typename SK::Circle_3 normal1       (ref_sphere,typename SK::Plane_3(0,1,0  ,-FT(0.5)-ref_sphere_center.y()));
  typename SK::Circle_3 normal2       (ref_sphere,typename SK::Plane_3(0,1,FT(0.1),-FT(0.5)-ref_sphere_center.y()-FT(0.1)*ref_sphere_center.z()));
  typename SK::Circle_3 north_polar   (ref_sphere,typename SK::Plane_3(0,1,1  ,-FT(1)-ref_sphere_center.y()-ref_sphere_center.z()));
  typename SK::Circle_3 south_polar   (ref_sphere,typename SK::Plane_3(0,1,-1 ,-FT(1)-ref_sphere_center.y()+ref_sphere_center.z()));

//=============TEST CLASSIFY================================================================
  assert(CGAL::classify(great_threaded,ref_sphere)==CGAL::THREADED);
  assert(CGAL::classify(threaded,ref_sphere)==CGAL::THREADED);
  assert(CGAL::classify(bipolar,ref_sphere)==CGAL::BIPOLAR);
  assert(CGAL::classify(normal1,ref_sphere)==CGAL::NORMAL);
  assert(CGAL::classify(normal2,ref_sphere)==CGAL::NORMAL);
  assert(CGAL::classify(north_polar,ref_sphere)==CGAL::POLAR);
  assert(CGAL::classify(south_polar,ref_sphere)==CGAL::POLAR);
  std::cout << "Test classify OK" << std::endl;

//=============TEST IS_THETA_MONOTONE=======================================================
  test_normal_circle_monotonicity<SK>(normal1,ref_sphere);
  test_normal_circle_monotonicity<SK>(normal2,ref_sphere);
  test_threaded_circle_monotonicity<SK>(threaded,ref_sphere);
  test_threaded_circle_monotonicity<SK>(great_threaded,ref_sphere);
  test_polar_circle_monotonicity<SK>(north_polar,ref_sphere);
  test_polar_circle_monotonicity<SK>(south_polar,ref_sphere);
  std::cout << "Test monotonicity using functor and global function OK" << std::endl;
  std::cout << "Test Make_theta_monotone_3 OK" << std::endl;

//=============TEST COMPARE_THETA_3=========================================================
  {
    //compare y_extremal_points and theta_extremal_points
    typename SK::Compare_theta_3 cmp_theta=SK().compare_theta_3_object(ref_sphere);
    typename SK::Compare_theta_z_3 cmp_theta_z=SK().compare_theta_z_3_object(ref_sphere);

    typename SK::Circle_3 normal_cut_M0(ref_sphere,typename SK::Plane_3(1,0,FT(0.1)  ,-FT(0.5)-ref_sphere_center.x()-FT(0.1)*ref_sphere_center.z()));
    std::vector <typename SK::Circular_arc_point_3> y_xtrems;
    CGAL::y_extremal_points(normal_cut_M0.diametral_sphere(),std::back_inserter(y_xtrems));
    assert( cmp_theta(y_xtrems[0],y_xtrems[1])==CGAL::LARGER );

    typename SK::FT zcoord=CGAL::SphericalFunctors::extremal_points_z_coordinate<SK>(normal_cut_M0,ref_sphere);
    std::variant<typename SK::Circle_3, std::pair<typename SK::Circular_arc_point_3,unsigned>> objs[2];
    SK().intersect_3_object()(normal_cut_M0,typename SK::Plane_3(0,0,1,-zcoord),objs);
    typename SK::Circular_arc_point_3 xtrms[2];
    xtrms[0]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[0])->first;
    xtrms[1]=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;

    assert( cmp_theta(xtrms[0],xtrms[1])==CGAL::LARGER );
    assert( cmp_theta(xtrms[0],y_xtrems[0])==CGAL::SMALLER );
    assert( cmp_theta(xtrms[1],y_xtrems[1])==CGAL::LARGER);

    //cut by several planes
    typename SK::Circular_arc_point_3 inter_threaded[32];
    typename SK::Circular_arc_point_3 inter_great[32];
    for (unsigned i=0;i<16;++i){
      fill_intersections<SK>(threaded,ref_sphere,inter_threaded,i);
      fill_intersections<SK>(great_threaded,ref_sphere,inter_great,i);
    }

    for (unsigned i=0;i<31;++i){
      typename SK::Vector_3 vect_i=get_vector_in_plane<SK>(i);
      assert (cmp_theta(inter_threaded[i],inter_great[i])==CGAL::EQUAL );
      assert (cmp_theta(vect_i,vect_i)==CGAL::EQUAL );
      assert (cmp_theta_z(inter_threaded[i],inter_great[i])==CGAL::SMALLER );
      assert (CGAL::compare_theta(vect_i,vect_i)==CGAL::EQUAL );
      assert (CGAL::compare_theta(inter_threaded[i],inter_great[i],ref_sphere)==CGAL::EQUAL );
      assert (CGAL::compare_theta_z(inter_threaded[i],inter_great[i],ref_sphere)==CGAL::SMALLER );
      for (unsigned j=i+1;j<32;++j){
        typename SK::Vector_3 vect_j=get_vector_in_plane<SK>(j);
        assert( cmp_theta(inter_threaded[i],inter_great[j])==CGAL::SMALLER );
        assert( cmp_theta(inter_threaded[i],vect_j)==CGAL::SMALLER );
        assert( cmp_theta(vect_i,inter_great[j])==CGAL::SMALLER );
        assert( cmp_theta_z(inter_threaded[i],inter_great[j])==CGAL::SMALLER );
        assert( CGAL::compare_theta(inter_threaded[i],inter_great[j],ref_sphere)==CGAL::SMALLER );
        assert( CGAL::compare_theta(vect_i,inter_great[j],ref_sphere)==CGAL::SMALLER );
        assert( CGAL::compare_theta(inter_threaded[i],vect_j,ref_sphere)==CGAL::SMALLER );
        assert( CGAL::compare_theta_z(inter_threaded[i],inter_great[j],ref_sphere)==CGAL::SMALLER );
      }
    }
  }
  std::cout << "Test compare theta (z) global function and functor OK" << std::endl;

//=============TEST COMPARE_Z_AT_THETA_3====================================================
  {
    typename SK::Compare_z_at_theta_3 cmp_z_at_theta= SK().compare_z_at_theta_3_object(ref_sphere);

  //point vs arc
    //polar circle
    typename SK::Line_3 line_n(ref_sphere_center,north_polar.center()-ref_sphere_center);
    typename SK::Line_3 line_s(ref_sphere_center,south_polar.center()-ref_sphere_center);
    std::variant<typename SK::Circle_3, std::pair<typename SK::Circular_arc_point_3,unsigned>> objs[2];
    SK().intersect_3_object()(line_n,ref_sphere,objs);
    typename SK::Circular_arc_point_3 cn=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;
    SK().intersect_3_object()(line_s,ref_sphere,objs);
    typename SK::Circular_arc_point_3 cs=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;

    assert (cmp_z_at_theta(cn,typename SK::Circular_arc_3(north_polar,north_pole))==CGAL::LARGER );
    assert (cmp_z_at_theta(cs,typename SK::Circular_arc_3(north_polar,north_pole))==CGAL::SMALLER );
    assert (cmp_z_at_theta(cn,typename SK::Circular_arc_3(south_polar,south_pole))==CGAL::LARGER );
    assert (cmp_z_at_theta(cs,typename SK::Circular_arc_3(south_polar,south_pole))==CGAL::SMALLER );

    typename SK::Circular_arc_point_3 on_pt(ref_sphere_center.x(),ref_sphere_center.y()+radius,ref_sphere_center.z());
    assert (cmp_z_at_theta(south_pole,typename SK::Circular_arc_3(south_polar,south_pole))==CGAL::EQUAL);
    assert (cmp_z_at_theta(on_pt,typename SK::Circular_arc_3(south_polar,south_pole))==CGAL::EQUAL);
    assert (cmp_z_at_theta(north_pole,typename SK::Circular_arc_3(north_polar,north_pole))==CGAL::EQUAL);
    assert (cmp_z_at_theta(on_pt,typename SK::Circular_arc_3(north_polar,north_pole))==CGAL::EQUAL);

    //threaded circle
    assert (cmp_z_at_theta(cn,typename SK::Circular_arc_3(great_threaded))==CGAL::LARGER );
    assert (cmp_z_at_theta(cn,typename SK::Circular_arc_3(threaded))==CGAL::LARGER );
    assert (cmp_z_at_theta(north_pole,typename SK::Circular_arc_3(great_threaded))==CGAL::LARGER );
    assert (cmp_z_at_theta(north_pole,typename SK::Circular_arc_3(threaded))==CGAL::LARGER );
    assert (cmp_z_at_theta(cs,typename SK::Circular_arc_3(great_threaded))==CGAL::SMALLER );
    assert (cmp_z_at_theta(cs,typename SK::Circular_arc_3(threaded))==CGAL::SMALLER );
    assert (cmp_z_at_theta(south_pole,typename SK::Circular_arc_3(great_threaded))==CGAL::SMALLER );
    assert (cmp_z_at_theta(south_pole,typename SK::Circular_arc_3(threaded))==CGAL::SMALLER );
    assert (cmp_z_at_theta(on_pt,typename SK::Circular_arc_3(great_threaded))==CGAL::EQUAL);

    //normal circle
    typename SK::Line_3 line_1(ref_sphere_center,normal1.center()-ref_sphere_center);
    typename SK::Line_3 line_2(ref_sphere_center,normal2.center()-ref_sphere_center);
    SK().intersect_3_object()(line_1,ref_sphere,objs);
    typename SK::Circular_arc_point_3 c1=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;
    SK().intersect_3_object()(line_2,ref_sphere,objs);
    typename SK::Circular_arc_point_3 c2=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;

    typename SK::Circular_arc_point_3 xtrms1[2];
    typename SK::Circular_arc_point_3 xtrms2[2];
    CGAL::theta_extremal_points(normal1,ref_sphere,xtrms1);
    CGAL::theta_extremal_points(normal2,ref_sphere,xtrms2);

    assert ( cmp_z_at_theta(north_pole,typename SK::Circular_arc_3(normal1,xtrms1[0],xtrms1[1]))==CGAL::LARGER );
    assert ( cmp_z_at_theta(north_pole,typename SK::Circular_arc_3(normal1,xtrms1[1],xtrms1[0]))==CGAL::LARGER );
    assert ( cmp_z_at_theta(north_pole,typename SK::Circular_arc_3(normal2,xtrms2[0],xtrms2[1]))==CGAL::LARGER );
    assert ( cmp_z_at_theta(north_pole,typename SK::Circular_arc_3(normal2,xtrms2[1],xtrms2[0]))==CGAL::LARGER );
    assert ( cmp_z_at_theta(south_pole,typename SK::Circular_arc_3(normal1,xtrms1[0],xtrms1[1]))==CGAL::SMALLER );
    assert ( cmp_z_at_theta(south_pole,typename SK::Circular_arc_3(normal1,xtrms1[1],xtrms1[0]))==CGAL::SMALLER );
    assert ( cmp_z_at_theta(south_pole,typename SK::Circular_arc_3(normal2,xtrms2[0],xtrms2[1]))==CGAL::SMALLER );
    assert ( cmp_z_at_theta(south_pole,typename SK::Circular_arc_3(normal2,xtrms2[1],xtrms2[0]))==CGAL::SMALLER );

    assert ( cmp_z_at_theta(c1,typename SK::Circular_arc_3(normal1,xtrms1[0],xtrms1[1]))==CGAL::LARGER );
    assert ( cmp_z_at_theta(c1,typename SK::Circular_arc_3(normal1,xtrms1[1],xtrms1[0]))==CGAL::SMALLER );
    assert ( cmp_z_at_theta(c2,typename SK::Circular_arc_3(normal2,xtrms2[0],xtrms2[1]))==CGAL::LARGER );
    assert ( cmp_z_at_theta(c2,typename SK::Circular_arc_3(normal2,xtrms2[1],xtrms2[0]))==CGAL::SMALLER );
    assert ( cmp_z_at_theta(xtrms2[0],typename SK::Circular_arc_3(normal2,xtrms2[0],xtrms2[1]))==CGAL::EQUAL );
    assert ( cmp_z_at_theta(xtrms2[1],typename SK::Circular_arc_3(normal2,xtrms2[1],xtrms2[0]))==CGAL::EQUAL);
  //arc vs arc
    //polar vs polar
    assert ( cmp_z_at_theta(typename SK::Circular_arc_3(north_polar,north_pole),typename SK::Circular_arc_3(south_polar,south_pole),typename SK::Vector_3(0,1,0))==CGAL::EQUAL );

    //normal vs polar
    assert ( cmp_z_at_theta(typename SK::Circular_arc_3(north_polar,north_pole),
                            typename SK::Circular_arc_3(normal1,xtrms1[1],xtrms1[0]),
                            typename SK::Vector_3(0,1,0))==CGAL::SMALLER
    );
    assert ( cmp_z_at_theta(typename SK::Circular_arc_3(north_polar,north_pole),
                            typename SK::Circular_arc_3(normal1,xtrms1[0],xtrms1[1]),
                            typename SK::Vector_3(0,1,0))==CGAL::LARGER
    );

    //threaded  vs threaded
    assert ( cmp_z_at_theta(typename SK::Circular_arc_3(great_threaded),
                            typename SK::Circular_arc_3(threaded),
                            typename SK::Vector_3(typename SK::FT(0.25),-4,0))==CGAL::LARGER
    );

    //threaded vs polar
    assert ( cmp_z_at_theta(typename SK::Circular_arc_3(north_polar,north_pole),typename SK::Circular_arc_3(great_threaded),typename SK::Vector_3(0,1,0))==CGAL::EQUAL );
    assert ( cmp_z_at_theta(typename SK::Circular_arc_3(south_polar,south_pole),typename SK::Circular_arc_3(great_threaded),typename SK::Vector_3(0,1,0))==CGAL::EQUAL );
    assert ( cmp_z_at_theta(typename SK::Circular_arc_3(north_polar,north_pole),typename SK::Circular_arc_3(great_threaded),typename SK::Vector_3(-1.2,1.1,0))==CGAL::LARGER );
    assert ( cmp_z_at_theta(typename SK::Circular_arc_3(south_polar,south_pole),typename SK::Circular_arc_3(great_threaded),typename SK::Vector_3(-1.7,1.8,0))==CGAL::SMALLER );

    //normal vs normal
    assert ( cmp_z_at_theta(typename SK::Circular_arc_3(normal1,xtrms1[1],xtrms1[0]),
                            typename SK::Circular_arc_3(normal1,xtrms1[0],xtrms1[1]),
                            typename SK::Vector_3(FT(0.1),FT(0.8),0))==CGAL::LARGER
    );
    assert ( cmp_z_at_theta(typename SK::Circular_arc_3(normal2,xtrms2[0],xtrms2[1]),
                            typename SK::Circular_arc_3(normal1,xtrms1[0],xtrms1[1]),
                            typename SK::Vector_3(FT(0.1),FT(0.8),0))==CGAL::LARGER
    );
  }
  std::cout << "Test Compare_z_at_theta_3  OK" << std::endl;

//=============TEST COMPARE_Z_TO_RIGHT_3====================================================
  {
    typename SK::Compare_z_to_right_3 cmp_right= SK().compare_z_to_right_3_object(ref_sphere);
  //at intersection points global
    //normal vs normal
    typename SK::Circle_3 normal3 (ref_sphere,typename SK::Plane_3(0.08,1.1,0.9  ,-1-FT(0.08)*ref_sphere_center.x() -FT(1.1)*ref_sphere_center.y()-FT(0.9)*ref_sphere_center.z()));
    typename SK::Circle_3 normal4 (ref_sphere,typename SK::Plane_3(0.05,1.1,-0.9 ,-1-FT(0.05)*ref_sphere_center.x() -FT(1.1)*ref_sphere_center.y()+FT(0.9)*ref_sphere_center.z()));
    std::vector<std::variant<typename SK::Circle_3, std::pair<typename SK::Circular_arc_point_3,unsigned>>> objs;
    SK().intersect_3_object()(normal3,normal4,std::back_inserter(objs));
    typename SK::Circular_arc_point_3 int1=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;
    typename SK::Circular_arc_point_3 int2=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[0])->first;
    typename SK::Circular_arc_point_3 xtrms3[2];
    typename SK::Circular_arc_point_3 xtrms4[2];
    CGAL::theta_extremal_points(normal3,ref_sphere,xtrms3);
    CGAL::theta_extremal_points(normal4,ref_sphere,xtrms4);
    assert ( cmp_right(typename SK::Circular_arc_3(normal3,xtrms3[0],xtrms3[1]),typename SK::Circular_arc_3(normal4,xtrms4[1],xtrms4[0]),int1)==CGAL::SMALLER );
    assert ( cmp_right(typename SK::Circular_arc_3(normal3,xtrms3[0],xtrms3[1]),typename SK::Circular_arc_3(normal4,xtrms4[1],xtrms4[0]),int2)==CGAL::LARGER );
    //normal vs threaded
    objs.clear();
    SK().intersect_3_object()(threaded,normal4,std::back_inserter(objs));
    int1=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;
    int2=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[0])->first;
    assert ( cmp_right(typename SK::Circular_arc_3(threaded),typename SK::Circular_arc_3(normal4,xtrms4[1],xtrms4[0]),int1)==CGAL::SMALLER );
    assert ( cmp_right(typename SK::Circular_arc_3(threaded),typename SK::Circular_arc_3(normal4,xtrms4[1],xtrms4[0]),int2)==CGAL::LARGER );
    //normal vs polar
    objs.clear();
    typename SK::Circular_arc_point_3 xtrms2[2];
    CGAL::theta_extremal_points(normal2,ref_sphere,xtrms2);
    SK().intersect_3_object()(south_polar,normal2,std::back_inserter(objs));
    int1=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;
    int2=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[0])->first;
    assert ( cmp_right(typename SK::Circular_arc_3(south_polar,south_pole),typename SK::Circular_arc_3(normal2,xtrms2[0],xtrms2[1]),int1)==CGAL::LARGER );
    assert ( cmp_right(typename SK::Circular_arc_3(south_polar,south_pole),typename SK::Circular_arc_3(normal2,xtrms2[0],xtrms2[1]),int2)==CGAL::SMALLER );
    //polar vs polar
    typename SK::Circle_3 npolar(ref_sphere,typename SK::Plane_3(0,1,0.25,-FT(0.25)-1*ref_sphere_center.y()-FT(0.25)*ref_sphere_center.z()));
    typename SK::Circle_3 spolar(ref_sphere,typename SK::Plane_3(0,1,-0.25,-FT(0.25)-1*ref_sphere_center.y()+FT(0.25)*ref_sphere_center.z()));
    assert(CGAL::classify(npolar,ref_sphere)==CGAL::POLAR);
    assert(CGAL::classify(spolar,ref_sphere)==CGAL::POLAR);
    objs.clear();
    SK().intersect_3_object()(npolar,spolar,std::back_inserter(objs));
    int1=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;
    int2=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[0])->first;
    assert ( cmp_right(typename SK::Circular_arc_3(npolar,north_pole),typename SK::Circular_arc_3(spolar,south_pole),int1)==CGAL::SMALLER );
    assert ( cmp_right(typename SK::Circular_arc_3(npolar,north_pole),typename SK::Circular_arc_3(spolar,south_pole),int2)==CGAL::LARGER );
    //polar vs threaded
    objs.clear();
    SK().intersect_3_object()(south_polar,threaded,std::back_inserter(objs));
    int1=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;
    int2=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[0])->first;
    assert ( cmp_right(typename SK::Circular_arc_3(south_polar,south_pole),typename SK::Circular_arc_3(threaded),int1)==CGAL::LARGER );
    assert ( cmp_right(typename SK::Circular_arc_3(south_polar,south_pole),typename SK::Circular_arc_3(threaded),int2)==CGAL::SMALLER );
    //threaded vs threaded
    typename SK::Circle_3 threaded1 (ref_sphere,typename SK::Plane_3(0,1.1,0.9,-FT(1.1)*ref_sphere_center.y()-FT(0.9)*ref_sphere_center.z()));
    typename SK::Circle_3 threaded2 (ref_sphere,typename SK::Plane_3(0,1.1,-0.9,-FT(1.1)*ref_sphere_center.y()+FT(0.9)*ref_sphere_center.z()));
    objs.clear();
    SK().intersect_3_object()(threaded1,threaded2,std::back_inserter(objs));
    int1=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[1])->first;
    int2=std::get_if<std::pair<typename SK::Circular_arc_point_3,unsigned> >(&objs[0])->first;
    assert ( cmp_right(typename SK::Circular_arc_3(threaded1),typename SK::Circular_arc_3(threaded2),int1)==CGAL::SMALLER );
    assert ( cmp_right(typename SK::Circular_arc_3(threaded1),typename SK::Circular_arc_3(threaded2),int2)==CGAL::LARGER );
  //tangency tests global
    {
      typename SK::Circular_arc_point_3 tgt_pt1(ref_sphere_center.x(),ref_sphere_center.y()+radius,ref_sphere_center.z());
      //polar vs polar
      assert ( cmp_right(typename SK::Circular_arc_3(north_polar,north_pole),typename SK::Circular_arc_3(south_polar,south_pole),tgt_pt1)==CGAL::LARGER    );
      //polar vs threaded
      assert ( cmp_right(typename SK::Circular_arc_3(great_threaded),typename SK::Circular_arc_3(south_polar,south_pole),tgt_pt1)==CGAL::LARGER    );
      assert ( cmp_right(typename SK::Circular_arc_3(north_polar,north_pole),typename SK::Circular_arc_3(great_threaded),tgt_pt1)==CGAL::LARGER    );
      //polar vs normal
      typename SK::Circle_3 normal5(ref_sphere,typename SK::Plane_3(0,1,0.8,-1-1*ref_sphere_center.y()-FT(0.8)*ref_sphere_center.z()) );
      assert(CGAL::classify(normal5,ref_sphere)==CGAL::NORMAL);
      typename SK::Circular_arc_point_3 xtrms5[2];
      CGAL::theta_extremal_points(normal5,ref_sphere,xtrms5);
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[0],xtrms5[1]),typename SK::Circular_arc_3(north_polar,north_pole),tgt_pt1)==CGAL::LARGER );
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[0],xtrms5[1]),typename SK::Circular_arc_3(south_polar,south_pole),tgt_pt1)==CGAL::LARGER );
      //threaded vs threaded
      typename SK::Circle_3 threaded3 (ref_sphere,typename SK::Plane_3(0,0.1,1,-FT(0.1)-FT(0.1)*ref_sphere_center.y()-1*ref_sphere_center.z()) );
      typename SK::Circle_3 threaded4 (ref_sphere,typename SK::Plane_3(0,0.1,-1,-FT(0.1)-FT(0.1)*ref_sphere_center.y()+1*ref_sphere_center.z()));
      assert(CGAL::classify(threaded3,ref_sphere)==CGAL::THREADED);
      assert(CGAL::classify(threaded4,ref_sphere)==CGAL::THREADED);
      assert ( cmp_right(typename SK::Circular_arc_3(threaded3),typename SK::Circular_arc_3(threaded4),tgt_pt1)==CGAL::LARGER );
      //threaded vs normal
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[0],xtrms5[1]),typename SK::Circular_arc_3(threaded3),tgt_pt1)==CGAL::LARGER );
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[0],xtrms5[1]),typename SK::Circular_arc_3(threaded4),tgt_pt1)==CGAL::LARGER );
      //normal vs normal
      typename SK::Circle_3 normal6(ref_sphere,typename SK::Plane_3(0,1,-0.8,-1-1*ref_sphere_center.y()+FT(0.8)*ref_sphere_center.z()) );
      assert(CGAL::classify(normal6,ref_sphere)==CGAL::NORMAL);
      typename SK::Circular_arc_point_3 xtrms6[2];
      CGAL::theta_extremal_points(normal6,ref_sphere,xtrms6);
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[0],xtrms5[1]),typename SK::Circular_arc_3(normal6,xtrms6[1],xtrms6[0]),tgt_pt1)==CGAL::LARGER );
    }
  //tangency tests at theta=0
    {
      typename SK::Circular_arc_point_3 tgt_pt1(ref_sphere_center.x()+radius,ref_sphere_center.y(),ref_sphere_center.z());
      typename SK::Circle_3 north_polar0  (ref_sphere,typename SK::Plane_3(1,0,1  ,-1-ref_sphere_center.x()-ref_sphere_center.z()));
      typename SK::Circle_3 south_polar0  (ref_sphere,typename SK::Plane_3(1,0,-1 ,-1-ref_sphere_center.x()+ref_sphere_center.z()));
      //polar vs polar
      assert ( cmp_right(typename SK::Circular_arc_3(north_polar0,north_pole),typename SK::Circular_arc_3(south_polar0,south_pole),tgt_pt1)==CGAL::LARGER    );
      //polar vs threaded
      assert ( cmp_right(typename SK::Circular_arc_3(great_threaded),typename SK::Circular_arc_3(south_polar0,south_pole),tgt_pt1)==CGAL::LARGER    );
      assert ( cmp_right(typename SK::Circular_arc_3(north_polar0,north_pole),typename SK::Circular_arc_3(great_threaded),tgt_pt1)==CGAL::LARGER    );
      //polar vs normal
      typename SK::Circle_3 normal5(ref_sphere,typename SK::Plane_3(1,0,0.8,-1-1*ref_sphere_center.x()-FT(0.8)*ref_sphere_center.z()) );
      assert(CGAL::classify(normal5,ref_sphere)==CGAL::NORMAL);
      typename SK::Circular_arc_point_3 xtrms5[2];
      CGAL::theta_extremal_points(normal5,ref_sphere,xtrms5);
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),typename SK::Circular_arc_3(north_polar0,north_pole),tgt_pt1)==CGAL::LARGER );
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),typename SK::Circular_arc_3(south_polar0,south_pole),tgt_pt1)==CGAL::LARGER );
      //threaded vs threaded
      typename SK::Circle_3 threaded3 (ref_sphere,typename SK::Plane_3(0.1,0,1,-FT(0.1)-FT(0.1)*ref_sphere_center.x()-1*ref_sphere_center.z()) );
      typename SK::Circle_3 threaded4 (ref_sphere,typename SK::Plane_3(0.1,0,-1,-FT(0.1)-FT(0.1)*ref_sphere_center.x()+1*ref_sphere_center.z()));
      assert(CGAL::classify(threaded3,ref_sphere)==CGAL::THREADED);
      assert(CGAL::classify(threaded4,ref_sphere)==CGAL::THREADED);
      assert ( cmp_right(typename SK::Circular_arc_3(threaded3),typename SK::Circular_arc_3(threaded4),tgt_pt1)==CGAL::LARGER );
      //threaded vs normal
      assert( CGAL::SphericalFunctors::is_upper_arc<SK>(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),ref_sphere)==false );
      assert ( cmp_right(typename SK::Circular_arc_3(threaded3),typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),tgt_pt1)==CGAL::SMALLER );
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),typename SK::Circular_arc_3(threaded3),tgt_pt1)==CGAL::LARGER );
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),typename SK::Circular_arc_3(threaded4),tgt_pt1)==CGAL::LARGER );
      //normal vs normal
      typename SK::Circle_3 normal6(ref_sphere,typename SK::Plane_3(1,0,-0.8,-1-1*ref_sphere_center.x()+FT(0.8)*ref_sphere_center.z()) );
      assert(CGAL::classify(normal6,ref_sphere)==CGAL::NORMAL);
      typename SK::Circular_arc_point_3 xtrms6[2];
      CGAL::theta_extremal_points(normal6,ref_sphere,xtrms6);
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),typename SK::Circular_arc_3(normal6,xtrms6[0],xtrms6[1]),tgt_pt1)==CGAL::LARGER );
    }
  //at intersections at theta=0
    {
      typename SK::Circular_arc_point_3 int_pt1(ref_sphere_center.x()+radius,ref_sphere_center.y(),ref_sphere_center.z());
      typename SK::Circle_3 north_polar_at_0  (ref_sphere,typename SK::Plane_3(1,1,1  ,-1-ref_sphere_center.x()-ref_sphere_center.y()-ref_sphere_center.z()));
      typename SK::Circle_3 south_polar_at_0  (ref_sphere,typename SK::Plane_3(1,1,-1 ,-1-ref_sphere_center.x()-ref_sphere_center.y()+ref_sphere_center.z()));
      assert(CGAL::classify(north_polar_at_0,ref_sphere)==CGAL::POLAR);
      assert(CGAL::classify(south_polar_at_0,ref_sphere)==CGAL::POLAR);
      //polar vs polar
      assert ( cmp_right(typename SK::Circular_arc_3(north_polar_at_0,north_pole),typename SK::Circular_arc_3(south_polar_at_0,south_pole),int_pt1)==CGAL::SMALLER    );
      //polar vs threaded
      assert ( cmp_right(typename SK::Circular_arc_3(great_threaded),typename SK::Circular_arc_3(south_polar_at_0,south_pole),int_pt1)==CGAL::SMALLER    );
      assert ( cmp_right(typename SK::Circular_arc_3(north_polar_at_0,north_pole),typename SK::Circular_arc_3(great_threaded),int_pt1)==CGAL::SMALLER    );
      //polar vs normal
      typename SK::Circle_3 normal5(ref_sphere,typename SK::Plane_3(1,1,0.8,-1-ref_sphere_center.x()-ref_sphere_center.y()-FT(0.8)*ref_sphere_center.z()) );
      assert(CGAL::classify(normal5,ref_sphere)==CGAL::NORMAL);
      typename SK::Circular_arc_point_3 xtrms5[2];
      CGAL::theta_extremal_points(normal5,ref_sphere,xtrms5);
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),typename SK::Circular_arc_3(north_polar_at_0,north_pole),int_pt1)==CGAL::SMALLER );
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),typename SK::Circular_arc_3(south_polar_at_0,south_pole),int_pt1)==CGAL::SMALLER );
      //threaded vs threaded
      typename SK::Circle_3 threaded3 (ref_sphere,typename SK::Plane_3(0,1,1,-1*ref_sphere_center.y()-1*ref_sphere_center.z()) );
      assert(CGAL::classify(threaded3,ref_sphere)==CGAL::THREADED);
      assert ( cmp_right(typename SK::Circular_arc_3(threaded3),typename SK::Circular_arc_3(great_threaded),int_pt1)==CGAL::SMALLER );
      //~ //threaded vs normal
      assert( CGAL::SphericalFunctors::is_upper_arc<SK>(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),ref_sphere)==false );
      assert ( cmp_right(typename SK::Circular_arc_3(threaded3),typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),int_pt1)==CGAL::LARGER );
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),typename SK::Circular_arc_3(threaded3),int_pt1)==CGAL::SMALLER );
      //~ //normal vs normal
      typename SK::Circle_3 normal6(ref_sphere,typename SK::Plane_3(1,-1,0.8,-1-ref_sphere_center.x()+ref_sphere_center.y()-FT(0.8)*ref_sphere_center.z()) );
      assert(CGAL::classify(normal6,ref_sphere)==CGAL::NORMAL);
      typename SK::Circular_arc_point_3 xtrms6[2];
      CGAL::theta_extremal_points(normal6,ref_sphere,xtrms6);
      assert ( cmp_right(typename SK::Circular_arc_3(normal5,xtrms5[1],xtrms5[0]),typename SK::Circular_arc_3(normal6,xtrms6[1],xtrms6[0]),int_pt1)==CGAL::SMALLER );
    }
  }
  std::cout << "Test Compare_z_to_right_3 OK" << std::endl;

//=============TEST THETA_EXTREMAL_POINT(S)=================================================
  {
    test_extremal_points<SK>(normal1,ref_sphere,false);
    test_extremal_points<SK>(normal2,ref_sphere,false);
    test_extremal_points<SK>(typename SK::Circle_3(ref_sphere,typename SK::Plane_3(1,0.1,0.1,-FT(0.5)-ref_sphere_center.x()-FT(0.1)*ref_sphere_center.y()-FT(0.1)*ref_sphere_center.z())),ref_sphere,true);
    test_extremal_points<SK>(typename SK::Circle_3(ref_sphere,typename SK::Plane_3(-1,0.1,0.1,-FT(0.5)+ref_sphere_center.x()-FT(0.1)*ref_sphere_center.y()-FT(0.1)*ref_sphere_center.z())),ref_sphere,false);

    typename SK::Circle_3 circle_shq(ref_sphere,typename SK::Plane_3(-1,-0.25,0.1,-FT(1.03)+ref_sphere_center.x()+FT(0.25)*ref_sphere_center.y()-FT(0.1)*ref_sphere_center.z()));
    assert ( CGAL::SphericalFunctors::half_quadrant<SK>(CGAL::theta_extremal_point(circle_shq,ref_sphere,false),ref_sphere) ==
             CGAL::SphericalFunctors::half_quadrant<SK>(CGAL::theta_extremal_point(circle_shq,ref_sphere,true),ref_sphere)
    );
    test_extremal_points<SK>(circle_shq,ref_sphere,false);
  }
  std::cout << "Test theta_extremal_point(s) OK" << std::endl;
}
