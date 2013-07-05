// Sylvain Pion, 2008.

#include <CGAL/Uncertain.h>
#include <cassert>
#include <iostream>
#include <CGAL/enum.h>
#include <CGAL/exceptions.h>
#include <CGAL/use.h>

// The assert() macro, e.g. on Leopard, requires an int (due to the use of builtin_expect()).
// So I cast to bool manually (needed when the automatic conversion is switched off).
#define bool_assert(X) assert(bool(X))



// generic test, for both enums and bool
template < typename T >
void test()
{
	std::cout << "Testing Uncertain<T>\n";

	typedef CGAL::Uncertain<T> U;

	// Nested types.
	typename U::value_type t = T();
	typename U::Uncertain_conversion_exception e = CGAL::Uncertain_conversion_exception("");

	// Constructors, assignment.
	const T t0 = static_cast<T>(0);
	const T t2 = CGAL::indeterminate<T>();
	const U u = U::indeterminate();
	const U uuu = CGAL::indeterminate<U>();
	const U uu;
	const U v = t0;
	U w = U (T(), T());
	U v2 = u;
        CGAL_USE(v2);
	w = (w = t0);

	// Various functions
	t = u.inf();
	t = u.sup();

	bool_assert(! u.is_certain());
	bool_assert( uu.is_certain());
	bool_assert(  v.is_certain());
	bool_assert(t0 == v.make_certain());

#ifndef CGAL_NO_UNCERTAIN_CONVERSION_OPERATOR
	T conv = v;
	bool_assert(conv == t0);
#endif

	U indet = U::indeterminate();
        CGAL_USE(indet);

	t = CGAL::inf(u);
	t = CGAL::sup(u);

	bool_assert(! CGAL::is_certain(u));
	bool_assert(  CGAL::is_certain(v));
	bool_assert(  CGAL::is_certain(t0));
	bool_assert(  CGAL::is_indeterminate(u));
	bool_assert(  CGAL::is_indeterminate(uuu));
	bool_assert(! CGAL::is_indeterminate(v));
	bool_assert(! CGAL::is_indeterminate(t0));
	bool_assert(t0 == CGAL::get_certain(t0));
	bool_assert(t0 == CGAL::get_certain(v));
	bool_assert(t0 == CGAL::make_certain(t0));
	bool_assert(t0 == CGAL::make_certain(v));

	// Exceptions
	bool ok = true;

        // It is important to let that protection with CGAL_assertion_code,
        // because the following test checks that
        // CGAL::get_certain(Uncertain) does throw an assertion, *but if
        // CGAL_NDEBUG is set*.
        // -- Laurent Rineau, 2013/03/26
	CGAL_assertion_code( ok = false );
	try { CGAL::get_certain(u); }
	catch (CGAL::Assertion_exception) { ok = true; }
	bool_assert(ok);

	ok = false;
	try {
	  CGAL::make_certain(u);
#ifndef CGAL_NO_UNCERTAIN_CONVERSION_OPERATOR
	  T t = u;
	  CGAL_USE(t);
#endif
	}
	catch (CGAL::Uncertain_conversion_exception) { ok = true; }
	bool_assert(ok);

	U u2 = CGAL::make_uncertain(u);
	U u3 = CGAL::make_uncertain(T());
        CGAL_USE(u2);
        CGAL_USE(u3);

	// Operators
	bool_assert( v == v );
	bool_assert( v == t0 );
	bool_assert( t0 == v );
	bool_assert( ! (v != v) );
	bool_assert( ! (v != t0) );
	bool_assert( ! (t0 != v) );

	CGAL_USE(t);
	CGAL_USE(t2);
}


// test only for enums
template < typename T >
void test_enum()
{
	std::cout << "Testing Uncertain<enum>\n";

	test<T>(); // generic test

	typedef CGAL::Uncertain<T> U;

	T n = static_cast<T>(-1);
	T z = static_cast<T>(0);
	T p = static_cast<T>(1);

	U indet = U::indeterminate();

	// <, <=, >, >=
	bool_assert(U(n) < U(z));
	bool_assert(U(z) < U(p));
	bool_assert(U(n) < z);
	bool_assert(n < U(z));
	bool_assert(! CGAL::is_certain(indet < z));

	bool_assert(U(n) <= U(z));
	bool_assert(U(z) <= U(p));
	bool_assert(U(z) <= U(z));
	bool_assert(U(n) <= z);
	bool_assert(U(n) <= n);
	bool_assert(n <= U(z));
	bool_assert(z <= U(z));
	bool_assert(! CGAL::is_certain(indet <= z));
	bool_assert(! CGAL::is_certain(indet <= indet));

	bool_assert(U(z) > U(n));
	bool_assert(U(p) > U(z));
	bool_assert(U(z) > n);
	bool_assert(z > U(n));
	bool_assert(! CGAL::is_certain(z > indet));

	bool_assert(U(z) >= U(n));
	bool_assert(U(p) >= U(z));
	bool_assert(U(z) >= U(z));
	bool_assert(U(z) >= n);
	bool_assert(U(n) >= n);
	bool_assert(z >= U(n));
	bool_assert(z >= U(z));
	bool_assert(! CGAL::is_certain(z >= indet));
	bool_assert(! CGAL::is_certain(indet >= indet));

	bool_assert(-p == n);
	bool_assert(-n == p);
	bool_assert(-z == z);
}

// test only for enums with multiplication operator
template < typename T >
void test_mult_enum()
{
	std::cout << "Testing Uncertain<enum with operator*>\n";

	test_enum<T>(); // generic enum test.

	typedef CGAL::Uncertain<T> U;

	T n = static_cast<T>(-1);
	T z = static_cast<T>(0);
	T p = static_cast<T>(1);

	U indet = U::indeterminate();

	bool_assert(n*z == z);
	bool_assert(z*n == z);
	bool_assert(z*z == z);
	bool_assert(p*z == z);
	bool_assert(z*p == z);
	bool_assert(n*n == p);
	bool_assert(p*p == p);
	bool_assert(n*p == n);
	bool_assert(indet*z == z);
	bool_assert(z*indet == z);
	bool_assert(CGAL::is_indeterminate(p*indet));
	bool_assert(CGAL::is_indeterminate(n*indet));
}


template <typename T>
void error_on_bool(T) {}

void error_on_bool(bool) { abort(); }


// test only for bool
void test_bool()
{
	std::cout << "Testing Uncertain<bool>\n";

	test<bool>(); // generic test

	typedef CGAL::Uncertain<bool> U;

	U utrue = true;
	U ufalse = false;
	U indet = U::indeterminate();

	bool_assert(utrue);
	bool_assert(!ufalse);
	bool_assert(ufalse == !utrue);
	bool_assert(utrue != !utrue);
	bool_assert(false == !utrue);
	bool_assert(true != !utrue);

	bool_assert(utrue | utrue);
	bool_assert(utrue | ufalse);
	bool_assert(ufalse | utrue);
	bool_assert(! (ufalse | ufalse));
	bool_assert(utrue | true);
	bool_assert(utrue | false);
	bool_assert(ufalse | true);
	bool_assert(! (ufalse | false));
	bool_assert(true | utrue);
	bool_assert(true | ufalse);
	bool_assert(false | utrue);
	bool_assert(! (false | ufalse));

	bool_assert(utrue & utrue);
	bool_assert(! (utrue & ufalse));
	bool_assert(! (ufalse & utrue));
	bool_assert(! (ufalse & ufalse));
	bool_assert(utrue & true);
	bool_assert(! (utrue & false));
	bool_assert(! (ufalse & true));
	bool_assert(! (ufalse & false));
	bool_assert(true & utrue);
	bool_assert(! (true & ufalse));
	bool_assert(! (false & utrue));
	bool_assert(! (false & ufalse));

        bool_assert(!CGAL_AND(true, false));
        bool_assert(!CGAL_AND(false, true));
        bool_assert(CGAL_AND(true, true));
        bool_assert(!CGAL_AND(false, false));
        bool_assert(!CGAL_AND(utrue, ufalse));
        bool_assert(!CGAL_AND(ufalse, utrue));
        bool_assert(CGAL_AND(utrue, utrue));
        bool_assert(!CGAL_AND(ufalse, ufalse));

        bool_assert(CGAL_AND_3(true, true, true));
        bool_assert(!CGAL_AND_3(true, false, true));

        bool_assert(CGAL_OR(true, false));
        bool_assert(CGAL_OR(false, true));
        bool_assert(CGAL_OR(true, true));
        bool_assert(!CGAL_OR(false, false));
        bool_assert(CGAL_OR(utrue, ufalse));
        bool_assert(CGAL_OR(ufalse, utrue));
        bool_assert(CGAL_OR(utrue, utrue));
        bool_assert(!CGAL_OR(ufalse, ufalse));

        bool_assert(CGAL_OR_3(false, false, true));
        bool_assert(!CGAL_OR_3(false, false, false));

        try {
          bool_assert( CGAL_AND(utrue, utrue));
          bool_assert(!CGAL_AND(ufalse, ufalse));
          bool_assert(!CGAL_AND(ufalse, utrue));
          bool_assert(!CGAL_AND(utrue, ufalse));
          bool_assert(indet.is_same(CGAL_AND(indet, indet)));
          bool_assert(indet.is_same(CGAL_AND(utrue, indet)));
          bool_assert(indet.is_same(CGAL_AND(indet, utrue)));
          bool_assert(!CGAL_AND(ufalse, indet));
          bool_assert(!CGAL_AND(indet, ufalse));

          bool_assert( CGAL_OR(utrue, utrue));
          bool_assert(!CGAL_OR(ufalse, ufalse));
          bool_assert( CGAL_OR(ufalse, utrue));
          bool_assert( CGAL_OR(utrue, ufalse));
          bool_assert(indet.is_same(CGAL_OR(indet, indet)));
          bool_assert( CGAL_OR(utrue, indet));
          bool_assert( CGAL_OR(indet, utrue));
          bool_assert(indet.is_same(CGAL_OR(ufalse, indet)));
          bool_assert(indet.is_same(CGAL_OR(indet, ufalse)));

        } catch (CGAL::Uncertain_conversion_exception) {
#ifndef CGAL_CFG_NO_STATEMENT_EXPRESSIONS
         std::abort();
#endif
        }

	// Test exceptions
	bool ok = false;
	try { bool b = indet; CGAL_USE(b); }
	catch (CGAL::Uncertain_conversion_exception) { ok = true; }
	bool_assert(ok);
	// The following must throw.
	ok = false;
	try { U u = indet && utrue; u = indet || ufalse; }
	catch (CGAL::Uncertain_conversion_exception) { ok = true; }
	bool_assert(ok);
	// The following must not throw.
	try { bool b = utrue; b = ufalse; CGAL_USE(b); }
	catch (CGAL::Uncertain_conversion_exception) { bool_assert(false); }

	// certainly, possibly
	bool_assert(CGAL::certainly(true));
	bool_assert(CGAL::certainly(utrue));
	bool_assert(!CGAL::certainly(indet));
	bool_assert(!CGAL::certainly(false));
	bool_assert(!CGAL::certainly(ufalse));

	bool_assert(CGAL::possibly(true));
	bool_assert(CGAL::possibly(utrue));
	bool_assert(CGAL::possibly(indet));
	bool_assert(!CGAL::possibly(false));
	bool_assert(!CGAL::possibly(ufalse));

	bool_assert(!CGAL::certainly_not(true));
	bool_assert(!CGAL::certainly_not(utrue));
	bool_assert(!CGAL::certainly_not(indet));
	bool_assert(CGAL::certainly_not(false));
	bool_assert(CGAL::certainly_not(ufalse));

	bool_assert(!CGAL::possibly_not(true));
	bool_assert(!CGAL::possibly_not(utrue));
	bool_assert(CGAL::possibly_not(indet));
	bool_assert(CGAL::possibly_not(false));
	bool_assert(CGAL::possibly_not(ufalse));

	// conversion
	error_on_bool(0 ? CGAL::make_uncertain(true) : indet);
	// error_on_bool(0 ? true : indet); // fails due to the 2 possible conversions
}

void test_enum_cast()
{
	typedef CGAL::Uncertain<CGAL::Sign> Us;
	typedef CGAL::Uncertain<CGAL::Bounded_side> Ub;
	typedef CGAL::Uncertain<CGAL::Angle> Ua;

	Us s;
	Ub b = CGAL::enum_cast<CGAL::Bounded_side>(s);
	Ua a = CGAL::enum_cast<CGAL::Angle>(s);
        CGAL_USE(a);
	s = CGAL::enum_cast<CGAL::Sign>(b);
}

int main()
{
	test_bool();
	test_mult_enum<CGAL::Sign>();
	test_mult_enum<CGAL::Comparison_result>();
	test_mult_enum<CGAL::Orientation>();
	test_mult_enum<CGAL::Oriented_side>();
	test_enum<CGAL::Bounded_side>();
	test_enum<CGAL::Angle>();
	test_enum_cast();

	return 0;
}
