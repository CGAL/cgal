#ifndef POLYHEDRON_DEMO_STATISTICS_HELPERS_H
#define POLYHEDRON_DEMO_STATISTICS_HELPERS_H

#include <cmath>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <CGAL/squared_distance_3_0.h>
#include <map>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Polygon_mesh_processing/repair.h>


template<typename Polyhedron>
void angles(Polyhedron* poly, double& mini, double& maxi, double& ave)
{
  using namespace boost::accumulators;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;

  double rad_to_deg = 180. / CGAL_PI;

  accumulator_set< double,
    features< tag::min, tag::max, tag::mean > > acc;

  typename boost::property_map<Polyhedron, CGAL::vertex_point_t>::type
    vpmap = get(CGAL::vertex_point, *poly);
  BOOST_FOREACH(halfedge_descriptor h, halfedges(*poly))
  {
    if (face(h, *poly) == boost::graph_traits<Polyhedron>::null_face())
      continue;

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

template<typename Polyhedron>
void edges_length(Polyhedron* poly,
  double& mini, double& maxi, double& mean, double& mid,
  unsigned int& nb_degen)
{
  using namespace boost::accumulators;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;

  accumulator_set< double,
    features< tag::min, tag::max, tag::mean , tag::median> > acc;

  typename boost::property_map<Polyhedron, CGAL::vertex_point_t>::type
    vpmap = get(CGAL::vertex_point, *poly);
  nb_degen = 0;
  BOOST_FOREACH(edge_descriptor e, edges(*poly))
  {
    halfedge_descriptor h = halfedge(e, *poly);
    typename Kernel::Point_3 a = get(vpmap, source(h, *poly));
    typename Kernel::Point_3 b = get(vpmap, target(h, *poly));
    acc(CGAL::sqrt(CGAL::squared_distance(a, b)));

    if (a == b) ++nb_degen;
  }

  mini = extract_result< tag::min >(acc);
  maxi = extract_result< tag::max >(acc);
  mean = extract_result< tag::mean >(acc);
  mid =  extract_result< tag::median >(acc);
}

template<typename Polyhedron, typename VPmap>
unsigned int nb_degenerate_faces(Polyhedron* poly, VPmap vpmap)
{
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

  unsigned int nb = 0;
  BOOST_FOREACH(face_descriptor f, faces(*poly))
  {
    if (CGAL::Polygon_mesh_processing::is_degenerated(f, *poly, vpmap, Kernel()))
      ++nb;
  }
  return nb;
}

template<typename Polyhedron>
unsigned int nb_holes(Polyhedron* poly)
{
  //gets the number of holes
  //if is_closed is false, then there are borders (= holes)
  int n(0);
  int i = 0;

  // initialization : keep the original ids in memory and set them to 0
  std::vector<std::size_t> ids;
  for (typename Polyhedron::Halfedge_iterator it = poly->halfedges_begin();
      it != poly->halfedges_end(); ++it)
  {
    ids.push_back(it->id());
    it->id() = 0;
  }

  //if a border halfedge is found, increment the number of hole and set all the ids of the hole's border halfedges to 1 to prevent
  // the algorithm from counting them several times.
  for (typename Polyhedron::Halfedge_iterator it = poly->halfedges_begin();
       it != poly->halfedges_end();
       ++it)
  {
    if (it->is_border() && it->id() == 0){
      n++;
      typename Polyhedron::Halfedge_around_facet_circulator hf_around_facet = it->facet_begin();
      do {
        CGAL_assertion(hf_around_facet->id() == 0);
        hf_around_facet->id() = 1;
      } while (++hf_around_facet != it->facet_begin());
    }
  }
  //reset the ids to their initial value
  for (typename Polyhedron::Halfedge_iterator it = poly->halfedges_begin();
      it != poly->halfedges_end(); ++it)
  {
    it->id() = ids[i++];
  }
  return n;
}

#endif // POLYHEDRON_DEMO_STATISTICS_HELPERS_H

