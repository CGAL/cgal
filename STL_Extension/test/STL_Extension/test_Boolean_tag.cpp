#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <CGAL/tags.h>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <cassert>

using boost::is_same;
using CGAL::Tag_true;
using CGAL::Tag_false;

int main() {
  BOOST_MPL_ASSERT(( is_same< Tag_true::value_type, bool > ));
  BOOST_MPL_ASSERT(( is_same< Tag_true::type, Tag_true > ));
  BOOST_MPL_ASSERT_RELATION( (Tag_true::value), ==, true );
  assert( Tag_true() == true );

  BOOST_MPL_ASSERT(( is_same< Tag_false::value_type, bool > ));
  BOOST_MPL_ASSERT(( is_same< Tag_false::type, Tag_false > ));
  BOOST_MPL_ASSERT_RELATION( (Tag_false::value), ==, false );
  assert( Tag_false() == false );

  return 0;
}
