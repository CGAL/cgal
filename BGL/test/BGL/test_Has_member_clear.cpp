#include <CGAL/assertions.h>
#include <CGAL/boost/graph/internal/Has_member_clear.h>

struct with_clear {
  void clear() {}
};

struct wo_clear { };

struct with_clear_but_args {
  void clear(int) {}
};

struct with_clear_const {
  void clear() const {}
};


int main()
{
  using namespace CGAL::internal;
  static_assert(Has_member_clear<with_clear>::value);

  static_assert(!Has_member_clear<wo_clear>::value);

  static_assert(!Has_member_clear<with_clear_but_args>::value);

  static_assert(Has_member_clear<with_clear_const>::value);

  return 0;
}
