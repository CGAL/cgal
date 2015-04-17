
#ifndef CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H

#include <CGAL/boost/graph/Euler_operations.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>

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
      typedef boost::bimap<
        boost::bimaps::set_of<halfedge_descriptor>,
        boost::bimaps::multiset_of<double, std::greater<double> > >  Boost_bimap;
      typedef typename Boost_bimap::value_type                       long_edge;

      std::cout << "Split long edges (" << high << ")...";
      double sq_high = high*high;

      //collect long edges
      Boost_bimap long_edges;
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        double sqlen = sqlength(e);
        if(sqlen > sq_high)
          long_edges.insert(long_edge(halfedge(e, mesh_), sqlen));
      }

      //split long edges
      while (!long_edges.empty())
      {
        typename Boost_bimap::right_map::iterator eit = long_edges.right.begin();
        halfedge_descriptor he = eit->second;
        double sqlen = eit->first;
        long_edges.right.erase(eit);
        Point refinement_point = this->midpoint(he);

        //split edge
        halfedge_descriptor hnew = (!is_border(he, mesh_))
          ? CGAL::Euler::split_edge(he, mesh_)
          : CGAL::Euler::split_edge(opposite(he, mesh_), mesh_);

        vertex_descriptor vnew = target(hnew, mesh_);
        vpmap_[vnew] = refinement_point;

        //check sub-edges
        double sqlen_new = 0.25 * sqlen;
        if (sqlen_new > sq_high)
        {
          //if it was more than twice the "long" threshold, insert them
          long_edges.insert(long_edge(hnew, sqlen_new));
          long_edges.insert(long_edge(next(hnew, mesh_), sqlen_new));
        }

        //insert new edges to keep triangular faces, and update long_edges
        if (!is_border(hnew, mesh_))
        {
          halfedge_descriptor hnew2 = CGAL::Euler::split_face(hnew,
                                                              next(next(hnew, mesh_), mesh_),
                                                              mesh_);
          double sql = sqlength(hnew2);
          if (sql > sq_high)
            long_edges.insert(long_edge(hnew2, sql));
        }
        //do it again on the other side if we're not on boundary
        halfedge_descriptor hnew_opp = opposite(hnew, mesh_);
        if (!is_border(hnew_opp, mesh_))
        {
          halfedge_descriptor hnew2 = CGAL::Euler::split_face(prev(hnew_opp, mesh_),
                                                              next(hnew_opp, mesh_),
                                                              mesh_);
          double sql = sqlength(hnew2);
          if (sql > sq_high)
            long_edges.insert(long_edge(hnew2, sql));
        }
      }
      std::cout << " done." << std::endl;
#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("1-edge_split.off");
#endif
    }

    void collapse_short_edges(const double& low, const double& high)
    {
      typedef boost::bimap<
        boost::bimaps::set_of<halfedge_descriptor>,
        boost::bimaps::multiset_of<double, std::less<double> > >  Boost_bimap;
      typedef typename Boost_bimap::value_type                    short_edge;

      std::cout << "Collapse short edges (" << low << ", " << high << ")...";
      double sq_low = low*low;
      double sq_high = high*high;

      Boost_bimap short_edges;
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        double sqlen = sqlength(e);
        if (sqlen < sq_low)
          short_edges.insert(short_edge(halfedge(e, mesh_), sqlen));
      }

      unsigned int nb_collapses = 0;
      while (!short_edges.empty())
      {
        //the edge with shortest length
        typename Boost_bimap::right_map::iterator eit = short_edges.right.begin();
        halfedge_descriptor he = eit->second;
        double sqlen = eit->first;
        short_edges.right.erase(eit);

        //let's try to collapse he into vb
        vertex_descriptor va = target(he, mesh_);
        vertex_descriptor vb = source(he, mesh_);

        //handle the boundary case : an edge incident to boundary can be collapsed,
        //but only if the boundary vertex is kept, so re-insert opposite(he)
        //to collapse it
        if (is_border(va, mesh_))
        {
          if (!is_border(vb, mesh_))
            short_edges.insert(short_edge(opposite(he, mesh_), sqlen));
          continue;
        }

        if (degree(va, mesh_) < 3
          || degree(vb, mesh_) < 3
          || !CGAL::Euler::satisfies_link_condition(he, mesh_))//necessary to collapse
          continue;

        //check that collapse would not create an edge with length > high
        //iterate on vertices va_i of the one-ring of va
        bool collapse_ok = true;
        BOOST_FOREACH(halfedge_descriptor ha, halfedges_around_target(va, mesh_))
        {
          vertex_descriptor va_i = source(ha, mesh_);
          if (sqlength(vb, va_i) > sq_high)
          {
            collapse_ok = false;
            break;
          }
        }
        //if it is allowed, perform the collapse
        if (collapse_ok)
        {
          //"collapse va into vb along e"
          // remove edges incident to va and vb, because their lengths will change
          BOOST_FOREACH(halfedge_descriptor ha, halfedges_around_target(va, mesh_))
          {
            short_edges.left.erase(ha);
            short_edges.left.erase(opposite(ha, mesh_));
          }
          BOOST_FOREACH(halfedge_descriptor hb, halfedges_around_target(vb, mesh_))
          {
            short_edges.left.erase(hb);
            short_edges.left.erase(opposite(hb, mesh_));
          }

          CGAL_assertion_code(
            halfedge_descriptor en = next(he, mesh_);
            halfedge_descriptor enp = next(opposite(he, mesh_), mesh_);
            );

          //perform collapse
          Point target_point = vpmap_[vb];
          vertex_descriptor vkept = CGAL::Euler::collapse_edge(edge(he, mesh_), mesh_);
          vpmap_[vkept] = target_point;
          ++nb_collapses;

          CGAL_assertion(source(en, mesh_) == source(enp, mesh_));
          CGAL_expensive_assertion(is_triangle_mesh(mesh_));

          //insert new/remaining short edges
          BOOST_FOREACH(halfedge_descriptor ht, halfedges_around_target(vkept, mesh_))
          {
            double sqlen = sqlength(ht);
            if (sqlen < sq_low)
              short_edges.insert(short_edge(ht, sqlen));
          }

          std::cout << ".";
          std::cout.flush();
        }
      }
      std::cout << " done (" << nb_collapses << " collapses)." << std::endl;
#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("2-edge_collapse.off");
#endif
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
    double sqlength(const vertex_descriptor& v1,
                    const vertex_descriptor& v2) const
    {
      return CGAL::squared_distance(vpmap_[v1], vpmap_[v2]);
    }

    double sqlength(const halfedge_descriptor& h) const
    {
      vertex_descriptor v1 = target(h, mesh_);
      vertex_descriptor v2 = source(h, mesh_);
      return sqlength(v1, v2);
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

    void dump(const char* filename) const
    {
      std::ofstream out(filename);
      out << mesh_;
      out.close();
    }

  private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;

  };//end class Incremenal_remesher
}//end namespace internal
}//end namesapce Polygon_mesh_processing
}//end namesapce CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
