#ifndef myeigen_h
#define myeigen_h
#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#include <Eigen/Core>

namespace Eigen {

	template<> struct NumTraits<CGAL::Gmpq>
	{
		typedef CGAL::Gmpq Real;
		typedef CGAL::Gmpq NonInteger;
		typedef CGAL::Gmpq Nested;

		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned,
			ReadCost = 10,
			AddCost = 100,
			MulCost = 1000
		};
	};

	template<> struct NumTraits<CGAL::Interval_nt_advanced>
	{
		typedef CGAL::Interval_nt_advanced Real;
		typedef CGAL::Interval_nt_advanced NonInteger;
		typedef CGAL::Interval_nt_advanced Nested;

		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned,
			ReadCost = 2,
			AddCost = 2,
			MulCost = 10
		};
	};

}

#endif
