#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <CGAL/tags.h>
#include <cassert>
#include <type_traits>

using std::is_same;
using CGAL::Tag_true;
using CGAL::Tag_false;

int main() {
  CGAL_static_assertion(is_same< Tag_true::value_type, bool >::value);
  CGAL_static_assertion(is_same< Tag_true::type, Tag_true >::value);
  CGAL_static_assertion(Tag_true::value== true);
  assert( Tag_true() == true );

  CGAL_static_assertion(is_same< Tag_false::value_type, bool >::value);
  CGAL_static_assertion(is_same< Tag_false::type, Tag_false >::value);
  CGAL_static_assertion(Tag_false::value == false);
  assert( Tag_false() == false );

  return 0;
}
