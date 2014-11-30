






class ArrConstructionFromTwoPointsTraits_2 {
public:

	typedef ArrTraits::X_monotone_curve_2				X_monotone_curve_2;
	typedef ArrTraits::Curve_2							Curve_2;
	typedef ArrTraits::Point_2 							Point_2; 

	/*!
	returns an \f$ x\f$-monotone curve connecting `p1` and `p2` (i.e., the 
	two input points are its endpoints). 
	*/ 
	X_monotone_curve_2 operator() ( Point_2 p1, Point_2 p2);

	/*!
	returns a curve connecting `p1` and `p2` (i.e., the 
	two input points are its endpoints). 
	*/ 
	Curve_2 operator() ( Point_2 p1, Point_2 p2);

}