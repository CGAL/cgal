// Sylvain Pion, 2008.

#include <CGAL/Uncertain.h>
#include <cassert>
#include <iostream>
#include <CGAL/enum.h>
#include <CGAL/exceptions.h>

// "unused variable" warning killer
template <typename T>
void use(T) {}

// generic test , for both enums and bool
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
	const U u;
	const U v = t0;
	U w = U (T(), T());
	U v2 = u;
	w = (w = t0);

	// Various functions
	t = u.inf();
	t = u.sup();

	assert(! u.is_certain());
	assert(  v.is_certain());
	assert(t0 == v.make_certain());

	T conv = v;
	assert(conv == t0);

	U indet = U::indeterminate();

	t = CGAL::inf(u);
	t = CGAL::sup(u);

	assert(! CGAL::is_certain(u));
	assert(  CGAL::is_certain(v));
	assert(  CGAL::is_certain(t0));
	assert(  CGAL::is_indeterminate(u));
	assert(! CGAL::is_indeterminate(v));
	assert(! CGAL::is_indeterminate(t0));
	assert(t0 == CGAL::get_certain(t0));
	assert(t0 == CGAL::get_certain(v));
	assert(t0 == CGAL::make_certain(t0));
	assert(t0 == CGAL::make_certain(v));

	// Exceptions
	bool ok = true;
	CGAL_assertion_code( ok = false );
	try { CGAL::get_certain(u); }
	catch (CGAL::Assertion_exception) { ok = true; }
	assert(ok);
	ok = false;
	try { CGAL::make_certain(u); T t = u; use(t); }
	catch (CGAL::Uncertain_conversion_exception) { ok = true; }
	assert(ok);

	U u2 = CGAL::make_uncertain(u);
	U u3 = CGAL::make_uncertain(T());

	// Operators
	assert( v == v );
	assert( v == t0 );
	assert( t0 == v );
	assert( ! (v != v) );
	assert( ! (v != t0) );
	assert( ! (t0 != v) );

	use(t);
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

	U indet;

	// <, <=, >, >=
	assert(U(n) < U(z));
	assert(U(z) < U(p));
	assert(U(n) < z);
	assert(n < U(z));
	assert(! CGAL::is_certain(indet < z));

	assert(U(n) <= U(z));
	assert(U(z) <= U(p));
	assert(U(z) <= U(z));
	assert(U(n) <= z);
	assert(U(n) <= n);
	assert(n <= U(z));
	assert(z <= U(z));
	assert(! CGAL::is_certain(indet <= z));
	assert(! CGAL::is_certain(indet <= indet));

	assert(U(z) > U(n));
	assert(U(p) > U(z));
	assert(U(z) > n);
	assert(z > U(n));
	assert(! CGAL::is_certain(z > indet));

	assert(U(z) >= U(n));
	assert(U(p) >= U(z));
	assert(U(z) >= U(z));
	assert(U(z) >= n);
	assert(U(n) >= n);
	assert(z >= U(n));
	assert(z >= U(z));
	assert(! CGAL::is_certain(z >= indet));
	assert(! CGAL::is_certain(indet >= indet));

	assert(-p == n);
	assert(-n == p);
	assert(-z == z);
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

	U indet;

	assert(n*z == z);
	assert(z*n == z);
	assert(z*z == z);
	assert(p*z == z);
	assert(z*p == z);
	assert(n*n == p);
	assert(p*p == p);
	assert(n*p == n);
	assert(indet*z == z);
	assert(z*indet == z);
	assert(CGAL::is_indeterminate(p*indet));
	assert(CGAL::is_indeterminate(n*indet));
}


// test only for bool
void test_bool()
{
	std::cout << "Testing Uncertain<bool>\n";

	test<bool>(); // generic test

	typedef CGAL::Uncertain<bool> U;

	U utrue = true;
	U ufalse = false;
	U indet;

	assert(utrue);
	assert(!ufalse);
	assert(ufalse == !utrue);
	assert(utrue != !utrue);
	assert(false == !utrue);
	assert(true != !utrue);

	assert(utrue | utrue);
	assert(utrue | ufalse);
	assert(ufalse | utrue);
	assert(! (ufalse | ufalse));
	assert(utrue | true);
	assert(utrue | false);
	assert(ufalse | true);
	assert(! (ufalse | false));
	assert(true | utrue);
	assert(true | ufalse);
	assert(false | utrue);
	assert(! (false | ufalse));

	assert(utrue & utrue);
	assert(! (utrue & ufalse));
	assert(! (ufalse & utrue));
	assert(! (ufalse & ufalse));
	assert(utrue & true);
	assert(! (utrue & false));
	assert(! (ufalse & true));
	assert(! (ufalse & false));
	assert(true & utrue);
	assert(! (true & ufalse));
	assert(! (false & utrue));
	assert(! (false & ufalse));

	// Test exceptions
	bool ok = false;
	try { bool b = indet; use(b); }
	catch (CGAL::Uncertain_conversion_exception) { ok = true; }
	assert(ok);
	// The following must throw.
	ok = false;
	try { U u = indet && utrue; u = indet || ufalse; }
	catch (CGAL::Uncertain_conversion_exception) { ok = true; }
	assert(ok);
	// The following must not throw.
	try { bool b = utrue; b = ufalse; }
	catch (CGAL::Uncertain_conversion_exception) { assert(false); }

	// certainly, possibly
	assert(CGAL::certainly(true));
	assert(CGAL::certainly(utrue));
	assert(!CGAL::certainly(indet));
	assert(!CGAL::certainly(false));
	assert(!CGAL::certainly(ufalse));

	assert(CGAL::possibly(true));
	assert(CGAL::possibly(utrue));
	assert(CGAL::possibly(indet));
	assert(!CGAL::possibly(false));
	assert(!CGAL::possibly(ufalse));
}

void test_enum_cast()
{
	typedef CGAL::Uncertain<CGAL::Sign> Us;
	typedef CGAL::Uncertain<CGAL::Bounded_side> Ub;
	typedef CGAL::Uncertain<CGAL::Angle> Ua;

	Us s;
	Ub b = CGAL::enum_cast<CGAL::Bounded_side>(s);
	Ua a = CGAL::enum_cast<CGAL::Angle>(s);
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
