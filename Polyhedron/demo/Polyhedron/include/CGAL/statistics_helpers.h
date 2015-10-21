#ifndef POLYHEDRON_DEMO_STATISTICS_HELPERS_H
#define POLYHEDRON_DEMO_STATISTICS_HELPERS_H

#include <cmath>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>

#include <map>
#include <boost/property_map/property_map.hpp>


using namespace boost::accumulators;

template<typename Polyhedron>
void angles(Polyhedron* poly, double& mini, double& maxi, double& ave)
{
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;

  double rad_to_deg = 180. / CGAL_PI;

  accumulator_set< double,
    features< tag::min, tag::max, tag::mean > > acc;

  typename boost::property_map<Polyhedron, CGAL::vertex_point_t>::type
    vpmap = get(CGAL::vertex_point, *poly);
  BOOST_FOREACH(halfedge_descriptor h, halfedges(*poly))
  {
    typename Kernel::Point_3 a = get(vpmap, source(h, *poly));
    typename Kernel::Point_3 b = get(vpmap, target(h, *poly));
    typename Kernel::Point_3 c = get(vpmap, target(next(h, *poly), *poly));

    typename Kernel::Vector_3 ba(b, a);
    typename Kernel::Vector_3 bc(b, c);
    double cos_angle = (ba * bc)
      / std::sqrt(ba.squared_length() * bc.squared_length());

    acc(std::acos(cos_angle) * rad_to_deg);
  }

  mini = extract_result< tag::min >(acc);
  maxi = extract_result< tag::max >(acc);
  ave = extract_result< tag::mean >(acc);
}

#endif // POLYHEDRON_DEMO_STATISTICS_HELPERS_H

