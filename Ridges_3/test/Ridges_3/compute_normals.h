#ifndef COMPUTE_NORMALS_H
#define COMPUTE_NORMALS_H

#include <CGAL/boost/graph/helpers.h>

#include <boost/foreach.hpp>

template <typename TriangleMesh, typename FaceVectorMap, typename Kernel>
const typename Kernel::Vector_3 
computeFacetsAverageUnitNormal(const TriangleMesh& tm,
                               typename boost::graph_traits<TriangleMesh>::vertex_descriptor v,
                               FaceVectorMap fvm,
                               const Kernel& )
{
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h;
  typename boost::graph_traits<TriangleMesh>::face_descriptor f;
  typename Kernel::Vector_3 sum(0., 0., 0.), n;

  CGAL::Halfedge_around_target_circulator<TriangleMesh> hedgeb(halfedge(v,tm),tm), hedgee = hedgeb;

  do
    {
      h = *hedgeb;
      if (is_border_edge(h,tm))
	{
	  hedgeb++;
	  continue;
	}

      f =  face(h,tm);
      n = get(fvm,f);
      sum = (sum + n);
      hedgeb++;
    }
  while (hedgeb != hedgee);
  sum = sum / std::sqrt(sum * sum);
  return sum;
}

template <typename TriangleMesh, typename FaceVectorMap,typename Kernel>
void compute_facets_normals(const TriangleMesh& tm,
                            FaceVectorMap fvm,
                            const Kernel& )
{
  typedef typename boost::property_traits<FaceVectorMap>::value_type Vector_3;

  typedef typename boost::property_map<TriangleMesh,CGAL::vertex_point_t>::const_type VPM;
  VPM vpm = get(CGAL::vertex_point,tm);
  BOOST_FOREACH(typename boost::graph_traits<TriangleMesh>::face_descriptor f, faces(tm)){
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h = halfedge(f,tm);
    Vector_3 normal = 
      CGAL::cross_product(get(vpm, target(h,tm)) -
			    get(vpm, target(opposite(h,tm),tm)),
                          get(vpm, target(next(h,tm),tm)) -
			    get(vpm, target(opposite(h,tm),tm)));
      put(fvm, f, normal / CGAL::sqrt(normal * normal));
  }
}



#endif
