#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/iterator.h>
#include <CGAL/Origin.h>

#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel    EPECK;
typedef EPECK::Point_2                                       Point_2;
typedef EPECK::Compute_x_2                                   Compute_x_2;
typedef EPECK::Construct_midpoint_2                          Construct_midpoint_2;
typedef EPECK::Collinear_2                                   Collinear_2;

typedef std::vector<Point_2>                                 Point_vector;
typedef Point_vector::iterator                               PV_it;
typedef Point_vector::const_iterator                         PV_cit;

void test_join_input_iterator_1()
{
  Point_vector pv;
  pv.push_back(Point_2(0.1, 0.1));

  typedef CGAL::Join_input_iterator_1<PV_it, Compute_x_2>    Join;

  Join join;
  Join join_bis(join);
  join_bis = join;

  Join join_ter(pv.begin());
  Join join_quater(pv.begin(), Compute_x_2());
  assert(join_ter == join_quater);

  assert(join_ter.current_iterator1() == pv.begin());
  assert(*join_ter == 0.1); // calls Compute_x_2
  assert(join_ter[0] == 0.1); // calls Compute_x_2

  Join join_quinquies(pv.end(), Compute_x_2());
  assert(join_ter != join_quinquies);
  assert(join_ter < join_quinquies);

  ++join_ter;
  assert(join_ter == join_quinquies);
  --join_ter;
  assert(join_ter.current_iterator1() == pv.begin());

  assert((join_ter++).current_iterator1() == pv.begin());
  assert(join_ter == join_quinquies);
  assert((join_ter--).current_iterator1() == pv.end());
  assert(join_ter.current_iterator1() == pv.begin());

  join_ter += 1;
  assert(join_ter == join_quinquies);
  join_ter -= 1;
  assert(join_ter.current_iterator1() == pv.begin());

  join_ter = join_ter + 1;
  assert(join_ter == join_quinquies);
  join_ter = join_ter - 1;
  assert(join_ter.current_iterator1() == pv.begin());
}

void test_join_input_iterator_2()
{
  Point_vector pv;
  pv.push_back(Point_2(-1, -1));
  pv.push_back(Point_2(1, 1));

  typedef CGAL::Join_input_iterator_2<PV_cit, PV_cit, Construct_midpoint_2>    Join;

  PV_cit first = pv.begin(), second = ++(pv.begin());

  Join join;
  Join join_bis(join);
  join_bis = join;

  Join join_ter(first, second);
  Join join_quater(first, second, Construct_midpoint_2());
  assert(join_ter == join_quater);

  assert(join_ter.current_iterator1() == first);
  assert(join_ter.current_iterator2() == second);
  assert(*join_ter == CGAL::ORIGIN); // calls Construct_midpoint_2
  assert(join_ter[0] == CGAL::ORIGIN); // calls Construct_midpoint_2

  Join join_quinquies(second, pv.end(), Construct_midpoint_2());
  assert(join_ter != join_quinquies);
  assert(join_ter < join_quinquies);

  ++join_ter;
  assert(join_ter == join_quinquies);
  --join_ter;
  assert(join_ter.current_iterator1() == first);
  assert(join_ter.current_iterator2() == second);

  assert((join_ter++).current_iterator1() == first);
  assert(join_ter == join_quinquies);
  assert((join_ter--).current_iterator1() == second);
  assert(join_ter.current_iterator1() == first);
  assert(join_ter.current_iterator2() == second);

  join_ter += 1;
  assert(join_ter == join_quinquies);
  join_ter -= 1;
  assert(join_ter.current_iterator1() == first);
  assert(join_ter.current_iterator2() == second);

  join_ter = join_ter + 1;
  assert(join_ter == join_quinquies);
  join_ter = join_ter - 1;
  assert(join_ter.current_iterator1() == first);
  assert(join_ter.current_iterator2() == second);
}

void test_join_input_iterator_3()
{
  Point_vector pv;
  pv.push_back(Point_2(-0.1, -0.1));
  pv.push_back(Point_2(CGAL::ORIGIN));
  pv.push_back(Point_2(0.1, 0.1));
  pv.push_back(Point_2(0.1, 0));

  typedef CGAL::Join_input_iterator_3<PV_cit, PV_cit, PV_cit, Collinear_2>    Join;

  PV_cit first = pv.begin(), second = ++(pv.begin()),
         third = ++(++(pv.begin())), fourth = --(pv.end());

  Join join;
  Join join_bis(join);
  join_bis = join;

  Join join_ter(first, second, third);
  Join join_quater(first, second, third, Collinear_2());
  assert(join_ter == join_quater);

  assert(join_ter.current_iterator1() == first);
  assert(join_ter.current_iterator2() == second);
  assert(join_ter.current_iterator3() == third);
  assert(*join_ter); // calls Collinear_2
  assert(join_ter[0]); // calls Collinear_2
  assert(!(join_ter[1])); // calls Collinear_2

  const Join join_quinquies(second, third, fourth, Collinear_2());
  assert(join_ter != join_quinquies);
  assert(join_ter < join_quinquies);

  ++join_ter;
  assert(join_ter == join_quinquies);
  --join_ter;
  assert(join_ter.current_iterator1() == first);
  assert(join_ter.current_iterator2() == second);
  assert(join_ter.current_iterator3() == third);

  assert((join_ter++).current_iterator1() == first);
  assert(join_ter == join_quinquies);
  assert((join_ter--).current_iterator1() == second);
  assert(join_ter.current_iterator1() == first);
  assert(join_ter.current_iterator2() == second);
  assert(join_ter.current_iterator3() == third);

  join_ter += 1;
  assert(join_ter == join_quinquies);
  join_ter -= 1;
  assert(join_ter.current_iterator1() == first);
  assert(join_ter.current_iterator2() == second);
  assert(join_ter.current_iterator3() == third);

  join_ter = join_ter + 1;
  assert(join_ter == join_quinquies);
  join_ter = join_ter - 1;
  assert(join_ter.current_iterator1() == first);
  assert(join_ter.current_iterator2() == second);
  assert(join_ter.current_iterator3() == third);
}

int main()
{
  test_join_input_iterator_1();
  test_join_input_iterator_2();
  test_join_input_iterator_3();

  return EXIT_SUCCESS;
}
