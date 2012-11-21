#include <fenv.h>
#include <iostream>
#include <limits>

int modes[4] = { FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD, FE_TONEAREST };

std::string str (int mode)
{
	switch (mode)
		{
			case FE_TOWARDZERO: return "FE_TOWARDZERO";
			case FE_UPWARD:     return "FE_UPWARD";
			case FE_DOWNWARD:   return "FE_DOWNWARD";
			case FE_TONEAREST:  return "FE_TONEAREST";
			default:            throw __LINE__;
		}
}

int fetestround ()
{
	volatile double eps = std::numeric_limits<double>::denorm_min();

	double x = -1.0;
	double y =  1.0;
	volatile double x_plus_eps  = x + eps;
	volatile double y_minus_eps = y - eps;

	if ((x == x_plus_eps) && (y == y_minus_eps))
		return FE_TONEAREST;
	if (y == y_minus_eps)
		return FE_UPWARD;
	if (x == x_plus_eps)
		return FE_DOWNWARD;
	return FE_TOWARDZERO;
}

int main (int argc, char* argv[])
{
	int errors = 0;

	int mode_get = fegetround();
	std::cout << "fegetround() = " << str(mode_get) << "   " << std::endl;

	for (int i=0; i<4; i++)
		{
			int mode_set = modes[i];
			fesetround (mode_set);
			std::cout << "fesetround (" << str(mode_set) << ")" << std::endl;
			
			int mode_get = fegetround();
			std::cout << "fegetround() = " << str(mode_get) << "   ";
			bool ok_get = mode_get == mode_set;
			if (!ok_get)
				errors++;
			std::cout << (ok_get ? "(ok)" : "(error)") << std::endl;

			int mode_test = fetestround();
			std::cout << "fetestround() = " << str(mode_test) << "   ";
			bool ok_test = mode_test == mode_set;
			if (!ok_test)
				errors++;
			std::cout << (ok_test ? "(ok)" : "(error)") << std::endl;
		}

	return errors;
}
