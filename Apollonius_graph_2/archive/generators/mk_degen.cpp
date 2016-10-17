#include <CGAL/basic.h>
#include <CGAL/Apollonius_graph_2/make_degenerate.h>

#include <iostream>
#include <iomanip>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Apollonius_graph_filtered_traits_2.h>

typedef CGAL::Simple_cartesian<double> CK;
typedef CGAL::Simple_cartesian<CGAL::MP_Float> EK;

typedef CGAL::Field_with_sqrt_tag  CM;
typedef CGAL::Integral_domain_without_division_tag        EM;
typedef CGAL::Apollonius_graph_filtered_traits_2<CK,CM,EK,EM> Traits;

typedef Traits::Site_2     Site;
typedef std::vector<Site>  Site_vector;

int main()
{
  Site_vector input, output;
  Site site;

  std::back_insert_iterator<Site_vector> oit(output);

  while (std::cin >> site) { input.push_back(site); }

  oit = CGAL::make_degenerate(input.begin(), input.end(), oit, Traits());

  std::cout << std::setprecision(17);
  for (std::vector<Site>::iterator it = output.begin();
       it != output.end(); ++it)
    {
      std::cout << (*it) << std::endl;
    }

  return 0;
}


