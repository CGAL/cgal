// Test program for the Boost Concept Checking Library.
//
// Sylvain Pion.

#include <boost/concept_check.hpp>
#include <boost/concept_archetype.hpp>

using namespace boost;

// Remove warnings.
template < typename T >
inline void use(const T&) {}

// The concept checking class for the Pipo concept.
template < typename Flute >
struct PipoConcept
{
  void constraints()
  {
     // There are 2 requirements.
     // The first is that Pipo refines the "Convertible to int" concept.
     function_requires< ConvertibleConcept< Flute, int > >();
     // The second is that it has a const member function pipo()
     // returning something convertible to int.
     const_constraints(flute);
  }

  void const_constraints(const Flute &f)
  {
     i = f.pipo();
  }

  Flute flute;
  int i;
};


// An archetype class for the concept.
class pipo_archetype
  : public convertible_to_archetype<int>
{
  struct proxy : public convertible_to_archetype<int> {};

public:

  const proxy & pipo() const { return static_object<proxy>::get(); }

};


// Define some generic functions using this concept.
namespace Internal {

// Internal function that actually requires the requirements.
template < typename T >
void
function_with_long_name_that_I_don_t_want_to_see_in_the_ERROR_message
    (const T&t)
{
  int i = t.pipo();
  use(i);
}

} // namespace Internal


// Here is the documented function that checks the requirements
// before calling the internal function.
template < typename T >
void
my_function(const T & t)
{
  function_requires< PipoConcept<T> >();

  Internal::
  function_with_long_name_that_I_don_t_want_to_see_in_the_ERROR_message(t);
}


// A class that matches the concept and a bit more.
struct A
{
  short pipo() const { return 0; }
  operator int() const { return 0; }
};


int main()
{
  A a;
  my_function(a);

  // "int" does not match the concept.
  // Uncomment the line below to see the effect.
  int i;
  // my_function(i);
  use(i);

  // Now try the archetype (without execution).
  if (false) {
    const pipo_archetype & pa = static_object<pipo_archetype>::get();
    my_function(pa);
  }

  return 0;
}
