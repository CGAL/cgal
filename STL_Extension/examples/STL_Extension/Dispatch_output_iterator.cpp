#include <vector>
#include <CGAL/iterator.h>
#include <boost/variant.hpp>
#include <boost/optional.hpp>

int main()
{
  std::vector<int> a;
  std::vector<double> b;
  std::vector<char> c;

  typedef CGAL::Dispatch_output_iterator< 
    CGAL::cpp0x::tuple<int, double, char>,
    CGAL::cpp0x::tuple<std::back_insert_iterator< std::vector<int> >,
                       std::back_insert_iterator< std::vector<double> >,
                       std::back_insert_iterator< std::vector<char> > 
                       > > Dispatch;

  Dispatch disp = CGAL::dispatch_output<int, double, char>(
    std::back_inserter(a),
    std::back_inserter(b),
    std::back_inserter(c));

  typedef boost::variant<int, double, char> var;
  var va = 23; var vb = 4.2; var vc = 'x';

  // goes to a
  *disp++ = va; 
  // goes to b
  *disp++ = vb; 
  // goes to c
  *disp++ = vc; 
  // goes to a
  *disp++ = 42;

  return 0;
}
