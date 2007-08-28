#include <predicates.h>

int main(int, char *[]) {

  initialize_rational_points();
 
  
  {
    std::cout << "\n\nTesting compare_point_to_circle_rule_c" << std::endl;
    Summary sum;
    for (int j=0; j< 50; ++j) {
      P pt= random_point();
      std::cout << "Point is " << pt << " 1/10000 1" << std::endl;
      SP t= random_sp(pt);
      CI rule_coordinate(plane_coordinate(j%2));
      std::pair<S,S> rp=random_sr_point(pt, rule_coordinate);
      std::vector<S> spheres;
      spheres.push_back(rp.first);
      spheres.push_back(rp.second);
      CGAL::Bbox_3 bb=rp.first.bbox() + rp.second.bbox();
      Traits tr(spheres.begin(), spheres.end());
      for (unsigned int i=0; i< 100; ++i) {
	P qpt(pt.x(), random_coordinate(bb.ymin(), bb.ymax()), random_coordinate(bb.zmin(), bb.zmax()));
	SP qsp= random_sp(qpt);
	CGAL::Comparison_result cr0= tr.compare_point_to_circle_rule_c(qsp,
								       K(1), K(0),
								       rule_coordinate,
								       plane_coordinate(0));
	CGAL::Comparison_result cr1= tr.compare_point_to_circle_rule_c(qsp,
								       K(1), K(0),
								       rule_coordinate,
								       plane_coordinate(1));
	sum(cr0);
	sum(cr1);
	CGAL::Comparison_result ccr0=CGAL::compare(qpt[plane_coordinate(0).index()], 
						   pt[plane_coordinate(0).index()]);
	CGAL::Comparison_result ccr1=CGAL::compare(qpt[plane_coordinate(1).index()], 
						   pt[plane_coordinate(1).index()]);
	if (cr0 != ccr0 || cr1 != ccr1) {
	  std::cout << "Error: " << j%2 << " ";
	  std::cout << qpt << " 1/10000 1" << std::endl;
	  std::cout << cr0 << " " << cr1 << " " << ccr0 << " " << ccr1 << std::endl;
	  tr.compare_point_to_circle_rule_c(qsp,
					    K(1), K(0),
					    plane_coordinate(i%2),
					    plane_coordinate(0));
	  tr.compare_point_to_circle_rule_c(qsp,
					    K(1), K(0),
					    plane_coordinate(i%2),
					    plane_coordinate(1));
	}
	CGAL_assertion(cr0==ccr0);
	CGAL_assertion(cr1==ccr1);

      }
      
      for (unsigned int i=0; i< 10; ++i) {
	SP sp= random_sp(pt);
	CGAL::Comparison_result cr0= tr.compare_point_to_circle_rule_c(sp,
								       K(1), K(0),
								       rule_coordinate,
								       plane_coordinate(0));
	CGAL::Comparison_result cr1= tr.compare_point_to_circle_rule_c(sp,
								       K(1), K(0),
								       rule_coordinate,
								       plane_coordinate(1));
	sum(cr0);
	sum(cr1);
	CGAL_assertion(cr0 == CGAL::EQUAL);
	CGAL_assertion(cr1==CGAL::EQUAL);
      }
     
    }
  }
  {
    Summary sum;
    std::cout << "\n\nTesting compare_point_to_circle_circle_c" << std::endl;
    for (int j=0; j< 50; ++j) {
      P pt= random_point();
      std::cout << "Point is " << pt << " 1/10000 1" << std::endl;
      SP t= random_sp(pt);
      std::pair<S,S> rp=random_ss_point(pt);
      std::cout << "Spheres are " << rp.first << ": " << rp.second << std::endl;
      std::vector<S> spheres;
      spheres.push_back(rp.first);
      spheres.push_back(rp.second);
      CGAL::Bbox_3 bb= rp.first.bbox()+rp.second.bbox();
      Traits tr(spheres.begin(), spheres.end());
      for (unsigned int i=0; i< 100; ++i) {
	P qpt(pt.x(), random_coordinate(bb.ymin(), bb.ymax()), 
	      random_coordinate(bb.zmin(), bb.zmax()));
	SP qsp= random_sp(qpt);
	CGAL::Comparison_result cr0= tr.compare_point_to_circle_circle_c(qsp,
									 K(0), K(1),
									 
									 plane_coordinate(0));
	CGAL::Comparison_result cr1= tr.compare_point_to_circle_circle_c(qsp,
									 K(0), K(1),
									 
									 plane_coordinate(1));
	sum(cr0);
	sum(cr1);
	CGAL::Comparison_result ccr0=CGAL::compare(qpt[plane_coordinate(0).index()], 
						   pt[plane_coordinate(0).index()]);
	CGAL::Comparison_result ccr1=CGAL::compare(qpt[plane_coordinate(1).index()], 
						   pt[plane_coordinate(1).index()]);
	if (cr0 != ccr0 || cr1 != ccr1) {
	  std::cout << "Error: ";
	  std::cout << qpt << " 1/10000 1" << std::endl;
	  std::cout << cr0 << " " << cr1 << " " << ccr0 << " " << ccr1 << std::endl;
	  tr.compare_point_to_circle_circle_c(qsp,
					      K(0), K(1),
					      plane_coordinate(0));
	  tr.compare_point_to_circle_circle_c(qsp,
					      K(0), K(1),
					      plane_coordinate(1));
	}
	CGAL_assertion(cr0==ccr0);
	CGAL_assertion(cr1==ccr1);
      }


      for (unsigned int i=0; i< 10; ++i) {
	SP sp= random_sp(pt);
	CGAL::Comparison_result cr0= tr.compare_point_to_circle_circle_c(sp,
									 K(0), K(1),
									 plane_coordinate(0));
	CGAL::Comparison_result cr1= tr.compare_point_to_circle_circle_c(sp,
									 K(0), K(1),
									 plane_coordinate(1));
	sum(cr0);
	sum(cr1);
	CGAL_assertion(cr0 == CGAL::EQUAL);
	CGAL_assertion(cr1==CGAL::EQUAL);
      }
    }
    
   
  }





 
  
  {
    Summary sum;
    std::cout << "\n\nTesting compare_center_to_circle_circle_c" << std::endl;
    for (int j=0; j< 50; ++j) {
      P pt= random_point();
      SP t= random_sp(pt);
      std::cout << "Point is " << pt << " 1/10000 1" << std::endl;
      std::pair<S,S> rp=random_ss_point(pt);
      std::cout << "Spheres are " << rp.first << ": " << rp.second << std::endl;
      std::vector<S> spheres;
      spheres.push_back(rp.first);
      spheres.push_back(rp.second);
      for (unsigned int i=0; i< 100; ++i) {
	S qs=random_sphere(pt[sweep.index()], rp.first.bbox()+rp.second.bbox());// random_point(), random_radius());
	spheres.push_back(qs);
      }
      spheres.push_back(S(pt, random_coordinate()));
      Traits tr(spheres.begin(), spheres.end());
      for (unsigned int i=0; i< 101; ++i) {
	{
	  CGAL::Comparison_result cr0= tr.compare_center_to_circle_circle_c(K(i+2),
									    t,
									    K(0), K(1),
									    plane_coordinate(0));
	  CGAL::Comparison_result cr1= tr.compare_center_to_circle_circle_c(K(i+2),
									    t,
									    K(0), K(1),
									    plane_coordinate(1));
	  sum(cr0);
	  sum(cr1);
	  CGAL::Comparison_result ccr0=CGAL::compare(spheres[i+2].center()[plane_coordinate(0).index()], 
						     pt[plane_coordinate(0).index()]);
	  CGAL::Comparison_result ccr1=CGAL::compare(spheres[i+2].center()[plane_coordinate(1).index()], 
						     pt[plane_coordinate(1).index()]);
	  if (cr0 != ccr0 || cr1 != ccr1) {
	    std::cout << "Error: " << pt.x() << std::endl;
	    std::cout << cr0 << " " << cr1 << " " << ccr0 << " " << ccr1 << std::endl;
	    std::cout <<  tr.debug_separating_plane(K(0), K(1)).oriented_side(pt) 
		      << " " << tr.debug_separating_plane(K(1), K(0)).oriented_side(pt)  << std::endl;
	    std::cout << spheres[0] << std::endl << std::endl;
	    std::cout << spheres[1] << std::endl << std::endl;
	    std::cout << spheres[i+2] << std::endl << std::endl;
	    std::cout << pt.x() << " " << spheres[i+2].center().y() << " "
		      << spheres[i+2].center().z() << " 1/10000 1"<< std::endl << std::endl;
	    std::cout << "Point is " << pt << " 1/10000 1" << std::endl << std::endl;
	    std::cout << "p " << tr.debug_separating_plane(K(0), K(1)) << std::endl << std::endl;
	    std::cout << "p " << tr.debug_equipower_plane(K(0), K(1)) << std::endl << std::endl;
	    tr.compare_center_to_circle_circle_c(K(i+2),
						 t,
						 K(0), K(1),
						 plane_coordinate(0));
	    tr.compare_center_to_circle_circle_c(K(i+2),
						 t,
						 K(0), K(1),
						 plane_coordinate(1));
	  }
	  CGAL_assertion(cr0==ccr0);
	  CGAL_assertion(cr1==ccr1);
	  if (i==100) {
	    CGAL_assertion(cr0 == CGAL::EQUAL);
	    CGAL_assertion(cr1 == CGAL::EQUAL);
	  }
	}
      }
    }
    
  }
 
  {
    Summary sum;
    Traits tr;
    std::vector<std::pair<P, SP> > pts;
    for (unsigned int i=0; i< 50; ++i) {
      P p= random_point();
      pts.push_back(std::make_pair(p, random_sp(p)));
      pts.push_back(std::make_pair(p, random_sp(p)));
    }

    std::cout << "\n\nTesting compare_points_c" << std::endl;
    for (unsigned int i=0; i< pts.size(); ++i) {
      for (unsigned int j=0; j<i; ++j) {
	for (unsigned int k=0; k<3; ++k) {
	  CGAL_AOS3_INTERNAL_NS::Coordinate_index ci(k);
	  CGAL::Comparison_result cr=tr.compare_points_c(pts[i].second, pts[j].second, ci);
	  CGAL_assertion(cr
			 == CGAL::compare(pts[i].first[ci.index()],
					  pts[j].first[ci.index()]));
	  sum(cr);
	}
      }
    }
  }
 
 
  
  {
    Summary sum;
    std::cout << "\n\nTesting compare_center_to_circle_rule_c" << std::endl;
    for (int j=0; j< 50; ++j) {
      P pt= random_point();
      SP t= random_sp(pt);
      CI rule_coordinate(plane_coordinate(j%2));
      std::pair<S,S> rp=random_sr_point(pt, rule_coordinate);
      std::vector<S> spheres;
      spheres.push_back(rp.first);
      spheres.push_back(rp.second);
      for (unsigned int i=0; i< 100; ++i) {
	S qs=random_sphere(pt.x(), rp.first.bbox()+rp.second.bbox());
	spheres.push_back(qs);
      }
      for (unsigned int i=0; i< 10; ++i ) {
	if (i%2 == 0) {
	  spheres.push_back(S(P(pt.x(), pt.y(), random_coordinate(rp.second.bbox().zmin(),
								  rp.second.bbox().zmax())),
			      2));
	} else {
	  spheres.push_back(S(P(pt.x(), random_coordinate(rp.first.bbox().ymin(),
							  rp.first.bbox().ymax()), pt.z()),
			      2));
	}
      }
      Traits tr(spheres.begin(), spheres.end());
      for (unsigned int i=0; i< 110; ++i) {
	CGAL::Comparison_result cr0= tr.compare_center_to_circle_rule_c(K(i+2),
									t,
									K(1), K(0),
									rule_coordinate,
									plane_coordinate(0));
	CGAL::Comparison_result cr1= tr.compare_center_to_circle_rule_c(K(i+2),
									t,
									K(1), K(0),
									rule_coordinate,
									plane_coordinate(1));
	sum(cr0);
	sum(cr1);
	CGAL::Comparison_result ccr0=CGAL::compare(spheres[i+2].center()[plane_coordinate(0).index()], 
						   pt[plane_coordinate(0).index()]);
	CGAL::Comparison_result ccr1=CGAL::compare(spheres[i+2].center()[plane_coordinate(1).index()], 
						   pt[plane_coordinate(1).index()]);
	if (cr0 != ccr0 || cr1 != ccr1) {
	  std::cout << "Error: " << j%2 << " ";
	  std::cout << spheres[i+2] << std::endl;
	  tr.compare_center_to_circle_rule_c(K(i+2),
					     t,
					     K(1), K(0),
					     plane_coordinate(i%2),
					     plane_coordinate(0));
	  tr.compare_center_to_circle_rule_c(K(i+2),
					     t,
					     K(1), K(0),
					     plane_coordinate(i%2),
					     plane_coordinate(1));
	}
	CGAL_assertion(cr0==ccr0);
	CGAL_assertion(cr1==ccr1);
      }
    }
  }
 





  return 0;
}
