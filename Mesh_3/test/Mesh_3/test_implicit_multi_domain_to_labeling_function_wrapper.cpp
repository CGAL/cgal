#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K_e_i;
typedef K_e_i::Point_3 Point_3;
typedef double (Function) (const Point_3&);

typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function*> Labeling_function;
typedef Labeling_function::Function_vector Function_vector;


double cube_function_1 (const Point_3& p)
{
  if( p.x() > 0 && p.x() < 2 &&
      p.y() > 0 && p.y() < 2 &&
      p.z() > 0 && p.z() < 2 )
    return -1.;
  return 1.;
}

double cube_function_2 (const Point_3& p)
{
  if( p.x() > 1 && p.x() < 3 &&
      p.y() > 1 && p.y() < 3 &&
      p.z() > 1 && p.z() < 3 )
    return -1.;
  return 1.;
}


void test_constructor_vfunc ()
{
  Function_vector vfunc;
  vfunc.push_back(cube_function_1);
  vfunc.push_back(cube_function_2);

  Labeling_function lfunc(vfunc);

  Point_3 p1(0.5, 0.5, 0.5);
  Point_3 p2(2.5, 2.5, 2.5);
  Point_3 p3(1.5, 1.5, 1.5);
  Point_3 p_out(4., 4., 4.);

  Labeling_function::return_type rp1 = lfunc(p1);
  Labeling_function::return_type rp2 = lfunc(p2);
  Labeling_function::return_type rp3 = lfunc(p3);
  Labeling_function::return_type rp_out = lfunc(p_out);

  assert(rp1 != 0);
  assert(rp2 != 0);
  assert(rp3 != 0);
  assert(rp_out == 0);

  assert(rp1 != rp2);
  assert(rp1 != rp3);
  assert(rp2 != rp3);
}

void test_constructor_vfunc_vvpos ()
{
  Function_vector vfunc;
  vfunc.push_back(cube_function_1);
  vfunc.push_back(cube_function_2);

  std::vector<std::vector<CGAL::Sign> > vvpos;
  std::vector<CGAL::Sign> vpos1;
  vpos1.push_back(CGAL::NEGATIVE);
  vpos1.push_back(CGAL::POSITIVE);
  vvpos.push_back(vpos1);
  std::vector<CGAL::Sign> vpos2;
  vpos2.push_back(CGAL::POSITIVE);
  vpos2.push_back(CGAL::NEGATIVE);
  vvpos.push_back(vpos2);
  std::vector<CGAL::Sign> vpos3;
  vpos3.push_back(CGAL::NEGATIVE);
  vpos3.push_back(CGAL::NEGATIVE);
  vvpos.push_back(vpos3);

  Labeling_function lfunc(vfunc, vvpos);

  Point_3 p1(0.5, 0.5, 0.5);
  Point_3 p2(2.5, 2.5, 2.5);
  Point_3 p3(1.5, 1.5, 1.5);
  Point_3 p_out(4., 4., 4.);

  Labeling_function::return_type rp1 = lfunc(p1);
  Labeling_function::return_type rp2 = lfunc(p2);
  Labeling_function::return_type rp3 = lfunc(p3);
  Labeling_function::return_type rp_out = lfunc(p_out);

  assert(rp1 != 0);
  assert(rp2 != 0);
  assert(rp3 != 0);
  assert(rp_out == 0);

  assert(rp1 != rp2);
  assert(rp1 != rp3);
  assert(rp2 != rp3);
}

void test_constructor_vfunc_vstr ()
{
  Function_vector vfunc;
  vfunc.push_back(cube_function_1);
  vfunc.push_back(cube_function_2);

  std::vector<std::string> vstr;
  vstr.push_back("-+");
  vstr.push_back("+-");
  vstr.push_back("--");

  Labeling_function lfunc(vfunc, vstr);

  Point_3 p1(0.5, 0.5, 0.5);
  Point_3 p2(2.5, 2.5, 2.5);
  Point_3 p3(1.5, 1.5, 1.5);
  Point_3 p_out(4., 4., 4.);

  Labeling_function::return_type rp1 = lfunc(p1);
  Labeling_function::return_type rp2 = lfunc(p2);
  Labeling_function::return_type rp3 = lfunc(p3);
  Labeling_function::return_type rp_out = lfunc(p_out);

  assert(rp1 != 0);
  assert(rp2 != 0);
  assert(rp3 != 0);
  assert(rp_out == 0);

  assert(rp1 != rp2);
  assert(rp1 != rp3);
  assert(rp2 != rp3);
}

int main ()
{
  test_constructor_vfunc();
  test_constructor_vfunc_vvpos();
  test_constructor_vfunc_vstr();
  return 0;
}
