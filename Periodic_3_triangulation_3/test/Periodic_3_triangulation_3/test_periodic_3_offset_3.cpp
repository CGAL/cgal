#include <iostream>
#include <fstream>

#include <CGAL/Timer.h>
#include <CGAL/Periodic_3_offset_3.h>

typedef CGAL::Periodic_3_offset_3 Offset;

#include <CGAL/_test_cls_periodic_3_offset_3.h>

int main(int, char**)
{
  CGAL::Timer timer;
  timer.start();
  _test_cls_periodic_3_offset_3( Offset() );

  std::cerr << timer.time() << " sec." << std::endl;
  return 0;
}
