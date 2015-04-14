
#ifndef CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H

#include <CGAL/boost/graph/Euler_operations.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <map>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

  template<typename PolygonMesh
         , typename VertexPointMap
         , typename GeomTraits
  >
  class Incremental_remesher
  {
    typedef PolygonMesh PM;
    typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PM>::edge_descriptor     edge_descriptor;
    typedef typename boost::graph_traits<PM>::vertex_descriptor   vertex_descriptor;

    typedef typename GeomTraits::Point_3 Point;

  public:
    Incremental_remesher(PolygonMesh& pmesh
                       , VertexPointMap& vpmap)
      : mesh_(pmesh)
      , vpmap_(vpmap)
    {}

    void split_long_edges(const double& high)
    {
      //collect long edges
      double sq_high = high*high;
      std::map<halfedge_descriptor, double/*squared length*/> long_edges;
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        double sqlen = sqlength(e);
        if(sqlen > sq_high)
          long_edges[halfedge(e, mesh_)] = sqlen;
      }

      //split long edges
      while (!long_edges.empty())
      {
        typename std::map<halfedge_descriptor, double>::iterator eit = long_edges.begin();
        edge_descriptor e = eit->first;
        double sqlen = eit->second;
        long_edges.erase(eit);
        Point refinement_point = this->midpoint(e);

        //split edge
        bool is_border = (face(halfedge(e, mesh_), mesh_) == boost::graph_traits<PM>::null_face());
        halfedge_descriptor hnew = (!is_border)
          ? CGAL::Euler::split_edge(halfedge(e, mesh_), mesh_)
          : CGAL::Euler::split_edge(opposite(halfedge(e, mesh_), mesh_), mesh_);

        vertex_descriptor vnew = target(hnew, mesh_);
        vpmap_[vnew] = refinement_point;

        //check sub-edges
        double sqlen_new = 0.25 * sqlen;
        if (sqlen_new > sq_high)
        {
          //if it was more than twice the "long" threshold, insert them
          long_edges[hnew] = sqlen_new;
          long_edges[next(hnew, mesh_)] = sqlen_new;
        }

        //insert new edges to keep triangular faces, and update long_edges
        if (face(hnew, mesh_) != boost::graph_traits<PM>::null_face())
        {
          halfedge_descriptor hnew2 = CGAL::Euler::split_face(hnew,
                                                              next(next(hnew, mesh_), mesh_),
                                                              mesh_);
          double sql = sqlength(hnew2);
          if (sql > sq_high)
            long_edges[hnew2] = sql;
        }
        //do it again on the other side if we're not on boundary
        halfedge_descriptor hnew_opp = opposite(hnew, mesh_);
        if (face(hnew_opp, mesh_) != boost::graph_traits<PM>::null_face())
        {
          halfedge_descriptor hnew2 = CGAL::Euler::split_face(prev(hnew_opp, mesh_),
                                                              next(hnew_opp, mesh_),
                                                              mesh_);
          double sql = sqlength(hnew2);
          if (sql > sq_high)
            long_edges[hnew2] = sql;
        }
      }
    }
    void collapse_short_edges(const double& low, const double& high)
    {
      ;
    }
    void equalize_valences()
    {
      ;
    }
    void tangential_relaxation()
    {
      ;
    }
    void project_to_surface()
    {
      ;
    }

  private:
    double sqlength(const halfedge_descriptor& h) const
    {
      vertex_descriptor v1 = target(h, mesh_);
      vertex_descriptor v2 = source(h, mesh_);
      return CGAL::squared_distance(vpmap_[v1], vpmap_[v2]);
    }
    double sqlength(const edge_descriptor& e) const
    {
      return sqlength(halfedge(e, mesh_));
    }

    Point midpoint(const edge_descriptor& e) const
    {
      Point p1 = vpmap_[target(halfedge(e, mesh_), mesh_)];
      Point p2 = vpmap_[source(halfedge(e, mesh_), mesh_)];
      return CGAL::midpoint(p1, p2);
    }

  private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;

  };//end class Incremenal_remesher
}//end namespace internal
}//end namesapce Polygon_mesh_processing
}//end namesapce CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
