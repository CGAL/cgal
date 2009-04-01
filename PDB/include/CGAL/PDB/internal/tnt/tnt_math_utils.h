#ifndef CGAL_TNT_MATH_UTILS_H
#define CGAL_TNT_MATH_UTILS_H
#include <CGAL/PDB/basic.h>
/* needed for fabs, sqrt() below */
#include <cmath>


namespace CGAL { namespace PDB { namespace TNT {

/**
	@returns hypotenuse of real (non-complex) scalars a and b by 
	avoiding underflow/overflow
	using (a * sqrt( 1 + (b/a) * (b/a))), rather than
	sqrt(a*a + b*b).
*/
template <class Real>
Real hypot(const Real &a, const Real &b)
{
	
	if (a== 0)
		return abs(b);
	else
	{
		Real c = b/a;
		//return fabs(a) * sqrt(1 + c*c);
		return abs(a) * sqrt(1 + c*c);
	}
}

}}}

#endif
/* MATH_UTILS_H */
