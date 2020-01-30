#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/squared_distance_2.h>

typedef CGAL::Simple_cartesian<double> SC;
typedef CGAL::Simple_cartesian<CGAL::Gmpq>  EC;

template <typename Kernel>
int fct() {

	typedef Kernel::Segment_2 Segment_2;
  const Segment_2 segi = {
    { -4.0380854964382, -1.9947196614192 },
    { 10.43442091460618, -0.5886833953492263 } };
  const Segment_2 segj = {
    { -11.5138934277993, -2.721011070186227 },
    { -8.822747585009402, -2.459560251317805 } };

  const auto dist = CGAL::squared_distance(segi, segj);
  std::cout << "#dist: " << dist << std::endl;

  return 0;
}

int main()
{
	fct<SC>();
	fct<CGAL::Epick>();
	fct<EC>();
	return 0;
}
