#ifndef POLYHEDRON_DEMO_STATISTICS_HELPERS_H
#define POLYHEDRON_DEMO_STATISTICS_HELPERS_H

#include <CGAL/Distance_3/Point_3_Point_3.h>
#include <CGAL/Distance_3/Point_3_Line_3.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/property_map/property_map.hpp>

#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <vector>

template<typename Set>
struct Angles_test_with_set
{
  const Set& set;
  typedef typename Set::value_type value_type;
  Angles_test_with_set(const Set& set):set(set) {}

  bool operator()(const value_type& f)
  {
    return (set.find(f) == set.end());
  }
};

struct Angles_test
{
  Angles_test() {}
  template<typename T>
  bool operator()(const T&)
  {
    return false;
  }
};

template<typename Mesh, typename Tester>
void compute_angles(Mesh* poly,Tester tester , double& mini, double& maxi, double& ave)
{
  using namespace boost::accumulators;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPMap;
  typedef typename CGAL::Kernel_traits< typename boost::property_traits<VPMap>::value_type >::Kernel Traits;
  double rad_to_deg = 180. / CGAL_PI;

  accumulator_set< double,
    features< tag::min, tag::max, tag::mean > > acc;

  VPMap vpmap = get(CGAL::vertex_point, *poly);
  for(halfedge_descriptor h : halfedges(*poly))
  {
    face_descriptor f = face(h, *poly);
    if (f == boost::graph_traits<Mesh>::null_face()
        || tester(f))
      continue;

    typename Traits::Point_3 a = get(vpmap, source(h, *poly));
    typename Traits::Point_3 b = get(vpmap, target(h, *poly));
    typename Traits::Point_3 c = get(vpmap, target(next(h, *poly), *poly));

    typename Traits::Vector_3 ba(b, a);
    typename Traits::Vector_3 bc(b, c);
    double cos_angle = (ba * bc)
      / std::sqrt(ba.squared_length() * bc.squared_length());

    acc(std::acos(cos_angle) * rad_to_deg);
  }

  mini = extract_result< tag::min >(acc);
  maxi = extract_result< tag::max >(acc);
  ave = extract_result< tag::mean >(acc);
}

template<typename Mesh>
void angles(Mesh* poly, double& mini, double& maxi, double& ave)
{
  compute_angles(poly, Angles_test(), mini, maxi, ave);
}

template<typename Mesh, typename Face_set>
void angles(Mesh* poly, const Face_set& faces, double& mini, double& maxi, double& ave)
{
  Angles_test_with_set<Face_set> tester(faces);
  compute_angles(poly, tester, mini, maxi, ave);
}



template<typename Mesh, typename Edge_range>
void edges_length(Mesh* poly, const Edge_range& range,
  double& mini, double& maxi, double& mean, double& mid,
  unsigned int& nb_degen)
{
  using namespace boost::accumulators;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPMap;
  typedef typename boost::property_traits<VPMap>::value_type  Point;

  accumulator_set< double,
    features< tag::min, tag::max, tag::mean , tag::median> > acc;

  VPMap vpmap = get(CGAL::vertex_point, *poly);
  nb_degen = 0;
  for(edge_descriptor e : range)
  {
    halfedge_descriptor h = halfedge(e, *poly);
    Point a = get(vpmap, source(h, *poly));
    Point b = get(vpmap, target(h, *poly));
    acc(CGAL::sqrt(CGAL::squared_distance(a, b)));

    if (a == b) ++nb_degen;
  }

  mini = extract_result< tag::min >(acc);
  maxi = extract_result< tag::max >(acc);
  mean = extract_result< tag::mean >(acc);
  mid =  extract_result< tag::median >(acc);
}

template<typename Mesh>
void edges_length(Mesh* poly,
  double& mini, double& maxi, double& mean, double& mid,
  unsigned int& nb_degen)
{
  edges_length(poly, edges(*poly), mini, maxi, mean, mid, nb_degen);
}
template<typename Mesh>
unsigned int nb_degenerate_faces(Mesh* poly)
{
  typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;

  std::vector<face_descriptor> degenerate_faces;
  CGAL::Polygon_mesh_processing::degenerate_faces(*poly, std::back_inserter(degenerate_faces));

  return static_cast<unsigned int>(degenerate_faces.size());
}

template<typename Mesh, typename IDMap>
unsigned int nb_holes(Mesh* poly, IDMap idmap)
{

  //gets the number of holes
  //if is_closed is false, then there are borders (= holes)
  int n(0);

  // initialization :set all ids to 0 in vector
  std::vector<std::size_t> ids;
  ids.resize(num_halfedges(*poly));
  for (std::size_t i=0; i< ids.size(); ++i )
  {
    ids[i] = 0;
  }

  //if a border halfedge is found, increment the number of hole and set all the ids of the hole's border halfedges in the vector to 1 to prevent
  // the algorithm from counting them several times.
  for (typename boost::graph_traits<Mesh>::halfedge_iterator it = halfedges(*poly).begin();
       it != halfedges(*poly).end();
       ++it)
  {
    typename boost::graph_traits<Mesh>::halfedge_descriptor hd(*it);
    if (is_border(hd, *poly) && ids[get(idmap, hd)] == 0){
      n++;
      CGAL::Halfedge_around_face_circulator<Mesh> hf_around_facet(hd, *poly), done(hf_around_facet);
      do {
        CGAL_assertion(ids[get(idmap, *hf_around_facet)] == 0);
        ids[get(idmap, *hf_around_facet)] = 1;
      } while (++hf_around_facet != done);
    }
  }
  //reset the ids to their initial value
  //for (typename Mesh::Halfedge_iterator it = poly->halfedges_begin();
  //    it != poly->halfedges_end(); ++it)
  //{
  //  it->id() = ids[i++];
  //}
  return n;
}
template<typename Mesh>
unsigned int nb_holes(Mesh* poly)
{
  typedef typename boost::property_map<Mesh, boost::halfedge_index_t>::type IDMap;
  IDMap idmap = get(boost::halfedge_index, *poly);
  return nb_holes(poly, idmap);
}
template<typename Mesh, typename Face_range>
void faces_area(Mesh* poly, const Face_range& range,
  double& mini, double& maxi, double& mean, double& mid)
{
  using namespace boost::accumulators;
  typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPMap;
  typedef typename boost::property_traits<VPMap>::value_type Point;

  accumulator_set< double,
    features< tag::min, tag::max, tag::mean , tag::median> > acc;

  VPMap vpmap = get(CGAL::vertex_point, *poly);
  for(face_descriptor f : range)
  {
    halfedge_descriptor h = halfedge(f, *poly);
    Point a = get(vpmap, target(h, *poly));
    Point b = get(vpmap, target(next(h, *poly), *poly));
    Point c = get(vpmap, target(next(next(h, *poly), *poly), *poly));
    acc(CGAL::sqrt(CGAL::squared_area(a,b,c)));
  }

  mini = extract_result< tag::min >(acc);
  maxi = extract_result< tag::max >(acc);
  mean = extract_result< tag::mean >(acc);
  mid =  extract_result< tag::median >(acc);
}

template<typename Mesh>
void faces_area(Mesh* poly,
  double& mini, double& maxi, double& mean, double& mid)
{
  faces_area(poly, faces(*poly), mini, maxi, mean, mid);
}

template<typename Mesh, typename Face_range>
void faces_aspect_ratio(Mesh* poly, const Face_range& range,
  double& min_altitude, double& mini, double& maxi, double& mean)
{
  using namespace boost::accumulators;
  typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPMap;
  typedef typename CGAL::Kernel_traits< typename boost::property_traits<VPMap>::value_type >::Kernel Traits;


  accumulator_set< double,
    features< tag::min, tag::max, tag::mean> > acc;

  min_altitude = std::numeric_limits<double>::infinity();
  typename boost::property_map<Mesh, CGAL::vertex_point_t>::type
    vpmap = get(CGAL::vertex_point, *poly);
  for(face_descriptor f : range)
  {
    halfedge_descriptor h = halfedge(f, *poly);
    typename Traits::Point_3 points[3];
    points[0] = get(vpmap, target(h, *poly));
    points[1] = get(vpmap, target(next(h, *poly), *poly));
    points[2] = get(vpmap, target(next(next(h, *poly), *poly), *poly));
    //Compute smallest altitude
    double min_alt = std::numeric_limits<double>::infinity();
    double longest_edge = 0;
    for(int i=0; i<3; ++i)
    {
      double alt = CGAL::sqrt(CGAL::squared_distance(points[(0+i)%3], typename Traits::Line_3(points[(1+i)%3], points[(2+i)%3])));
      double edge =  CGAL::sqrt(CGAL::squared_distance(points[(1+i)%3], points[(2+i)%3]));
      if(alt < min_alt) { min_alt = alt; }
      if(edge > longest_edge) { longest_edge = edge; }
    }
    //compute aspect-ratio
    acc(longest_edge/min_alt);

    if(min_alt < min_altitude) { min_altitude = min_alt; }
  }

  mini = extract_result< tag::min >(acc);
  maxi = extract_result< tag::max >(acc);
  mean = extract_result< tag::mean >(acc);
}


template<typename Mesh>
void faces_aspect_ratio(Mesh* poly,
  double& min_altitude, double& mini, double& maxi, double& mean)
{
  faces_aspect_ratio(poly, faces(*poly), min_altitude,  mini, maxi, mean);
}
#endif // POLYHEDRON_DEMO_STATISTICS_HELPERS_H

