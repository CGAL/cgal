#define CGAL_DEBUG
#include <CGAL/assertions.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/exceptions.h>

template <typename Check>
bool raises_exception(Check check)
{
  try
  {
    check();
  }
  catch(const CGAL::Failure_exception&)
  {
    return true;
  }
  return false;
}

int main()
{
  CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
#if defined(CGAL_EXPECT_EXPENSIVE_CHECKS) || defined(CGAL_CHECK_EXPENSIVE)
  const bool expected = true;
#else
  const bool expected = false;
#endif

  const bool assertion = raises_exception([] { CGAL_expensive_assertion(false); });
  const bool precondition = raises_exception([] { CGAL_expensive_precondition(false); });
  const bool postcondition = raises_exception([] { CGAL_expensive_postcondition(false); });
  return assertion == expected && precondition == expected && postcondition == expected ? 0 : 1;
}
