#include <CGAL/assertions.h>
#include <CGAL/boost/graph/internal/Has_member_clear.h>

struct with_clear {
  void clear() {}
};

struct wo_clear { };

struct with_clear_but_args { 
  void clear(int) {}
};

struct with_clear_but_const { 
  void clear() const {}
};


int main()
{
  using namespace CGAL::internal;
  CGAL_static_assertion(Has_member_clear<with_clear>::value);
  CGAL_static_assertion(!Has_member_clear<wo_clear>::value);
  CGAL_static_assertion(!Has_member_clear<with_clear_but_args>::value);
  CGAL_static_assertion(!Has_member_clear<with_clear_but_const>::value);
  return 0;
}
