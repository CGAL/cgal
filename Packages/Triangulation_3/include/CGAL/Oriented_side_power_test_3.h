#ifndef CGAL_ORIENTED_SIDE_POWER_TEST_3_H
#define CGAL_ORIENTED_SIDE_POWER_TEST_3_H

#define square(x) x*x
#define DEBUG_POWER_TEST_3 0

template < class FT,class Weight >
CGAL_Oriented_side
CGAL_power_test_3 (  const CGAL_Weighted_point_3<
		           CGAL_Point_3<CGAL_Cartesian<FT> >, Weight> &a,
		     const CGAL_Weighted_point_3<
		           CGAL_Point_3<CGAL_Cartesian<FT> >, Weight> &b,
		     const CGAL_Weighted_point_3<
		           CGAL_Point_3<CGAL_Cartesian<FT> >, Weight> &c,
		     const CGAL_Weighted_point_3<
		           CGAL_Point_3<CGAL_Cartesian<FT> >, Weight> &d,
		     const CGAL_Weighted_point_3<
		           CGAL_Point_3<CGAL_Cartesian<FT> >, Weight> &test)
{
  // comput the sign of
  // |1                            1         1        1       1     |
  // |ax                           bx        cx       dx      testx |
  // |ay                                                            |
  // |az                                                            |
  // |ax^2 + ay^2 + az^2 - ar^2                                     |

  //cout << "power test en dimension 3\n";

  FT FT0(0);

  FT dax = a.point().x() - test.point().x();
  FT day = a.point().y() - test.point().y();
  FT daz = a.point().z() - test.point().z();
  FT das = square(dax) + square(day) + square(daz) -
    a.weight()  + test.weight() ;
  
  FT dbx = b.point().x() - test.point().x();
  FT dby = b.point().y() - test.point().y();
  FT dbz = b.point().z() - test.point().z();
  FT dbs = square(dbx) + square(dby) + square(dbz) -
    b.weight()  + test.weight() ;
  
  FT dcx = c.point().x() - test.point().x();
  FT dcy = c.point().y() - test.point().y();
  FT dcz = c.point().z() - test.point().z();
  FT dcs = square(dcx) + square(dcy) + square(dcz) -
    c.weight()  + test.weight() ;
  
  FT ddx = d.point().x() - test.point().x();
  FT ddy = d.point().y() - test.point().y();
  FT ddz = d.point().z() - test.point().z();
  FT dds = square(ddx) + square(ddy) + square(ddz) -
    d.weight()  + test.weight() ;
  
  FT det = CGAL_det4x4_by_formula(dax,day,daz,das,
				  dbx,dby,dbz,dbs,
				  dcx,dcy,dcz,dcs,
				  ddx,ddy,ddz,dds);
  
  return (det > FT0) ? CGAL_ON_NEGATIVE_SIDE  : 
    ((det == FT0) ? CGAL_ON_ORIENTED_BOUNDARY :
     CGAL_ON_POSITIVE_SIDE);
}


// power test for 2 dimensions triangulation (p,q,r,test : coplanar)

template <class FT,class Weight >
CGAL_Oriented_side
CGAL_power_test_3( const CGAL_Weighted_point_3<
		           CGAL_Point_3<CGAL_Cartesian<FT> >, Weight> &p,
		     const CGAL_Weighted_point_3<
		           CGAL_Point_3<CGAL_Cartesian<FT> >, Weight> &q,
		     const CGAL_Weighted_point_3<
		           CGAL_Point_3<CGAL_Cartesian<FT> >, Weight> &r,
		     const CGAL_Weighted_point_3<
		           CGAL_Point_3<CGAL_Cartesian<FT> >, Weight> &test)
{

	FT FT0(0);
	FT FT1(1);
	const FT &px = p.x();
	const FT &py = p.y();
	const FT &pz = p.z();
	const FT &pw = p.weight();
	const FT &qx = q.x();
	const FT &qy = q.y();
	const FT &qz = q.z();
	const FT &qw = q.weight();
	const FT &rx = r.x();
	const FT &ry = r.y();
	const FT &rz = r.z();
	const FT &rw = r.weight();
	const FT &tx = test.x();
	const FT &ty = test.y();
	const FT &tz = test.z();
	const FT &tw = test.weight();
	FT det;
	FT det_inf ;

	if (DEBUG_POWER_TEST_3) {
	  cout << "power test en dimension 2\n";
	}
	
	CGAL_kernel_precondition( ! CGAL_collinear(p,q,r) );
	
	// -----------------------------------------------------------------
	det = CGAL_det4x4_by_formula(px, py, px*px + py*py + pz*pz -pw, FT1,
				     qx, qy, qx*qx + qy*qy + qz*qz -qw, FT1,
				     rx, ry, rx*rx + ry*ry + rz*rz -rw, FT1,
				     tx, ty, tx*tx + ty*ty + tz*tz -tw, FT1);
	if (DEBUG_POWER_TEST_3) {
	  cout << "det = projection suivant z : " << det << endl;
	}

	if (det != 0) {
	  det_inf = CGAL_det4x4_by_formula(px, py, px*px + py*py + pz*pz -pw, FT1,
					   qx, qy, qx*qx + qy*qy + qz*qz -qw, FT1,
					   rx, ry, rx*rx + ry*ry + rz*rz -rw, FT1,
					   FT0, FT0, FT1, FT0);
	  return (det*det_inf > FT0) ? CGAL_ON_NEGATIVE_SIDE  :
	    CGAL_ON_POSITIVE_SIDE;
	}



	// -----------------------------------------------------------------
	det = CGAL_det4x4_by_formula(px, pz, px*px + py*py + pz*pz -pw, FT1,
				     qx, qz, qx*qx + qy*qy + qz*qz -qw, FT1,
				     rx, rz, rx*rx + ry*ry + rz*rz -rw, FT1,
				     tx, tz, tx*tx + ty*ty + tz*tz -tw, FT1);
	if (DEBUG_POWER_TEST_3) {
	  cout << "det = projection suivant y : " << det << endl;
	}

	if (det != 0) {
	  det_inf = CGAL_det4x4_by_formula(px, pz, px*px + py*py + pz*pz -pw, FT1,
					   qx, qz, qx*qx + qy*qy + qz*qz -qw, FT1,
					   rx, rz, rx*rx + ry*ry + rz*rz -rw, FT1,
					   FT0, FT0, FT1, FT0);
	  return (det*det_inf > FT0) ? CGAL_ON_NEGATIVE_SIDE  :
	    CGAL_ON_POSITIVE_SIDE;
	}

	// -----------------------------------------------------------------
	det = CGAL_det4x4_by_formula(py, pz, px*px + py*py + pz*pz -pw, FT1,
				     qy, qz, qx*qx + qy*qy + qz*qz -qw, FT1,
				     ry, rz, rx*rx + ry*ry + rz*rz -rw, FT1,
				     ty, tz, tx*tx + ty*ty + tz*tz -tw, FT1);
	if (DEBUG_POWER_TEST_3) {
	  cout << "det = projection suivant x : " << det << endl;
	}
	
	det_inf = CGAL_det4x4_by_formula(py, pz, px*px + py*py + pz*pz -pw, FT1,
					 qy, qz, qx*qx + qy*qy + qz*qz -qw, FT1,
					 ry, rz, rx*rx + ry*ry + rz*rz -rw, FT1,
					 FT0, FT0, FT1, FT0);
	
	return (det*det_inf > FT0) ? CGAL_ON_NEGATIVE_SIDE  : 
	  ((det==FT0) ? CGAL_ON_ORIENTED_BOUNDARY :
	   CGAL_ON_POSITIVE_SIDE);
	
}

// power test for 1 dimension triangulation

template <class FT,class Weight >
CGAL_Oriented_side
CGAL_power_test_3(const CGAL_Weighted_point_3<CGAL_Point_3<CGAL_Cartesian<FT> >, 
		  Weight> &p,
		  const CGAL_Weighted_point_3<CGAL_Point_3<CGAL_Cartesian<FT> >, 
		  Weight> &q,
		  const CGAL_Weighted_point_3<CGAL_Point_3<CGAL_Cartesian<FT> >, 
		  Weight> &test)
{
	FT FT0(0);
	FT FT1(1);
	FT det;
	const FT &px = p.x();
	const FT &py = p.y();
	const FT &pz = p.z();
	const FT &pw = p.weight();
	const FT &qx = q.x();
	const FT &qy = q.y();
	const FT &qz = q.z();
	const FT &qw = q.weight();
	const FT &tx = test.x();
	const FT &ty = test.y();
	const FT &tz = test.z();
	const FT &tw = test.weight();

	if (DEBUG_POWER_TEST_3) {
	  cout << "power test en dimension 1\n";
	}

	CGAL_kernel_precondition( CGAL_collinear(p.point(),q.point(),test.point()) );
	CGAL_kernel_precondition( p.point() != q.point() );
	if (px!=qx)
	{	
		det = CGAL_det3x3_by_formula(px, px*px + py*py + pz*pz -pw, FT1,
		                             qx, qx*qx + qy*qy + qz*qz -qw, FT1,
		                             tx, tx*tx + ty*ty + tz*tz -tw, FT1);
		if (px<qx)
			det=-det;
	}
	else if (py!=qy)
	{
		det = CGAL_det3x3_by_formula(py, px*px + py*py + pz*pz -pw, FT1,
		                             qy, qx*qx + qy*qy + qz*qz -qw, FT1,
		                             ty, tx*tx + ty*ty + tz*tz -tw, FT1);
		if (py<qy)
			det=-det;
	}
	else
	{
		det = CGAL_det3x3_by_formula(pz, px*px + py*py + pz*pz -pw, FT1,
		                             qz, qx*qx + qy*qy + qz*qz -qw, FT1,
		                             tz, tx*tx + ty*ty + tz*tz -tw, FT1);
		if (pz<qz)
			det=-det;
	}
	return (det<FT0) ? CGAL_ON_NEGATIVE_SIDE
	                 : ((det==FT0) ? CGAL_ON_ORIENTED_BOUNDARY 
				: CGAL_ON_POSITIVE_SIDE);
}
#endif
