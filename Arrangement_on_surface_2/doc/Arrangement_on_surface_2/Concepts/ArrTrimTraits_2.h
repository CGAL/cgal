/*








*/

class ArrTrimTraits_2 {
public:

	typedef ArrTraits::X_monotone_curve_2				X_monotone_curve_2;
	typedef ArrTraits::Point_2 							Point_2; 

	/*!
	returns an \f$ x\f$-monotone trimmed version \f$ x\f$-monotone 
	connecting `p1` and `p2` (i.e., the 
	two input points are its endpoints). 
	*/ 
	X_monotone_curve_2 operator() ( X_monotone_curve_2 xcv,
		 						    Point_2 p1, Point_2 p2);


};


