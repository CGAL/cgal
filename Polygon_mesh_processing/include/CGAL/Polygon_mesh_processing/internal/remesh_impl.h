
#ifndef CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/get_border.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/range.hpp>

#include <map>
#include <list>
#include <vector>
#include <iterator>

namespace PMP = CGAL::Polygon_mesh_processing;


namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

  enum Halfedge_status {
    PATCH,       //h and hopp belong to the patch to be remeshed
    PATCH_BORDER,//h belongs to the patch, hopp is MESH
    MESH,        //h and hopp belong to the mesh, not the patch
    MESH_BORDER  //h belongs to the mesh, face(hopp, pmesh) == null_face()
  };

  template<typename PolygonMesh
         , typename FaceRange
         , typename VertexPointMap
         , typename GeomTraits
  >
  class Incremental_remesher
  {
    typedef PolygonMesh PM;
    typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PM>::edge_descriptor     edge_descriptor;
    typedef typename boost::graph_traits<PM>::vertex_descriptor   vertex_descriptor;
    typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;

    typedef typename GeomTraits::Point_3    Point;
    typedef typename GeomTraits::Vector_3   Vector_3;
    typedef typename GeomTraits::Plane_3    Plane_3;
    typedef typename GeomTraits::Triangle_3 Triangle_3;

    typedef std::list<Triangle_3>                                  Triangle_list;
    typedef typename Triangle_list::iterator                       Tr_iterator;
    typedef CGAL::AABB_triangle_primitive<GeomTraits, Tr_iterator> Primitive;
    typedef CGAL::AABB_traits<GeomTraits, Primitive>               Traits;
    typedef CGAL::AABB_tree<Traits>                                AABB_tree;

  public:
    Incremental_remesher(PolygonMesh& pmesh
                       , FaceRange face_range
                       , VertexPointMap& vpmap)
      : mesh_(pmesh)
      , vpmap_(vpmap)
      , own_tree_(true)
      , input_triangles_()
      , patch_(boost::begin(face_range), boost::end(face_range))
      , halfedge_status_map_()
    {
      CGAL_assertion(CGAL::is_triangle_mesh(mesh_));

      //build AABB tree of input surface
      //todo : add a constructor with aabb_tree as parameter
      //todo : do we really need to keep this copy of the input surface?
      BOOST_FOREACH(face_descriptor f, faces(mesh_))
      {
        halfedge_descriptor h = halfedge(f, mesh_);
        vertex_descriptor v1 = target(h, mesh_);
        vertex_descriptor v2 = target(next(h, mesh_), mesh_);
        vertex_descriptor v3 = target(next(next(h, mesh_), mesh_), mesh_);
        input_triangles_.push_back(Triangle_3(vpmap[v1], vpmap[v2], vpmap[v3]));
      }
      tree_ptr_ = new AABB_tree(input_triangles_.begin(), input_triangles_.end());

      tag_halfedges_status(face_range);
    }

    ~Incremental_remesher()
    {
      if (own_tree_)
        delete tree_ptr_;
    }
    
    // PMP book :
    // "visits all edges of the mesh
    //if an edge is longer than the given threshold `high`, the edge
    //is split at its midpoint and the two adjacent triangles are bisected (2-4 split)"
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
      unsigned int nb_splits = 0;
      while (!long_edges.empty())
      {
        //the edge with longest length
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

        ++nb_splits;
      }
      std::cout << " done ("<< nb_splits << " splits)." << std::endl;

#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("1-edge_split.off");
#endif
    }

    // PMP book :
    // "collapses and thus removes all edges that are shorter than a
    // threshold `low`. [...] testing before each collapse whether the collapse
    // would produce an edge that is longer than `high`"
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

        //handle the boundary case :
        //a boundary edge can be collapsed,
        //and an edge incident to boundary can be collapsed,
        //but only if the boundary vertex is kept, so re-insert opposite(he)
        //to collapse it
        if (!is_border_edge(he, mesh_))
        {
          if (is_border(va, mesh_))
          {
            if (!is_border(vb, mesh_))
              short_edges.insert(short_edge(opposite(he, mesh_), sqlen));
            continue;
          }
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
        }
      }
      std::cout << " done (" << nb_collapses << " collapses)." << std::endl;

#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("2-edge_collapse.off");
#endif
    }

    // PMP book :
    // "equalizes the vertex valences by flipping edges.
    // The target valence is 6 and 4 for interior and boundary vertices, resp.
    // The algo. tentatively flips each edge `e` and checks whether the deviation
    // to the target valences decreases. If not, the edge is flipped back"
    void equalize_valences()
    {
      std::cout << "Equalize valences...";
      unsigned int nb_flips = 0;
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        if (!is_flip_allowed(e))
          continue;

        halfedge_descriptor he = halfedge(e, mesh_);
        vertex_descriptor va = source(he, mesh_);
        vertex_descriptor vb = target(he, mesh_);
        vertex_descriptor vc = target(next(he, mesh_), mesh_);
        vertex_descriptor vd = target(next(opposite(he, mesh_), mesh_), mesh_);

        int deviation_pre = CGAL::abs(valence(va) - target_valence(va))
                          + CGAL::abs(valence(vb) - target_valence(vb))
                          + CGAL::abs(valence(vc) - target_valence(vc))
                          + CGAL::abs(valence(vd) - target_valence(vd));

        CGAL::Euler::flip_edge(he, mesh_);
        ++nb_flips;

        CGAL_assertion(
             (vc == target(he, mesh_) && vd == source(he, mesh_))
          || (vd == target(he, mesh_) && vc == source(he, mesh_)));

        int deviation_post = CGAL::abs(valence(va) - target_valence(va))
                          + CGAL::abs(valence(vb) - target_valence(vb))
                          + CGAL::abs(valence(vc) - target_valence(vc))
                          + CGAL::abs(valence(vd) - target_valence(vd));

        if (deviation_pre < deviation_post)
        {
          CGAL::Euler::flip_edge(he, mesh_);
          --nb_flips;
          CGAL_assertion(
               (va == source(he, mesh_) && vb == target(he, mesh_))
            || (vb == source(he, mesh_) && va == target(he, mesh_)));
        }
      }
      std::cout << "done. ("<< nb_flips << " flips)" << std::endl;

#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("3-edge_flips.off");
#endif
    }

    // PMP book :
    // "applies an iterative smoothing filter to the mesh.
    // The vertex movement has to be constrained to the vertex tangent plane [...]
    // smoothing algorithm with uniform Laplacian weights"
    void tangential_relaxation()
    {
      //todo : move border vertices along 1-dimensional features
      std::cout << "Tangential relaxation...";

      //todo : use boost::vector_property_map to improve computing time
      typedef std::map<vertex_descriptor, Vector_3> VNormalsMap;
      VNormalsMap vnormals;
      boost::associative_property_map<VNormalsMap> propmap_normals(vnormals);

      PMP::compute_vertex_normals(mesh_,
                                  propmap_normals,
                                  PMP::parameters::vertex_point_map(vpmap_).
                                  geom_traits(GeomTraits()));

      // at each vertex, compute barycenter of neighbors
      std::map<vertex_descriptor, Point> barycenters;
      BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
      {
        if (is_border(v, mesh_))
          continue;
        Vector_3 move = CGAL::NULL_VECTOR;
        unsigned int star_size = 0;
        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, mesh_))
        {
          move = move + Vector_3(vpmap_[v], vpmap_[source(h, mesh_)]);
          ++star_size;
        }
        move = (1. / (double)star_size) * move;
        barycenters[v] = vpmap_[v] + move;
      }

      // compute moves
      std::map<vertex_descriptor, Point> new_locations;
      BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
      {
        if (is_border(v, mesh_))
          continue;
        Vector_3 nv = boost::get(propmap_normals, v);
        Point qv = barycenters[v];
        new_locations[v] = qv + (nv * Vector_3(qv, vpmap_[v])) * nv;
      }

      // perform moves
      typedef typename std::map<vertex_descriptor, Point>::value_type VP_pair;
      BOOST_FOREACH(const VP_pair& vp, new_locations)
      {
        vpmap_[vp.first] = new_locations[vp.first];
      }

      CGAL_assertion(is_valid(mesh_));
      CGAL_assertion(is_triangle_mesh(mesh_));

      std::cout << "done." << std::endl;

#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("4-relaxation.off");
#endif
    }


    // PMP book :
    // "maps the vertices back to the surface"
    void project_to_surface()
    {
      //todo : handle the case of boundary vertices
      std::cout << "Project to surface...";

      BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
      {
        if (is_border(v, mesh_))
          continue;
        vpmap_[v] = tree_ptr_->closest_point(vpmap_[v]);
      }

      CGAL_assertion(is_valid(mesh_));
      CGAL_assertion(is_triangle_mesh(mesh_));

      std::cout << "done." << std::endl;

#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("5-project.off");
#endif
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

    int valence(const vertex_descriptor& v) const
    {
      return degree(v, mesh_);
    }

    int target_valence(const vertex_descriptor& v) const
    {
      return (is_border(v, mesh_)) ? 4 : 6;
    }

    bool is_flip_allowed(const edge_descriptor& e) const
    {
      if (is_border(e, mesh_))
        return false;//we can't flip border edges

      halfedge_descriptor he = halfedge(e, mesh_);
      Point p1 = vpmap_[target(he, mesh_)];
      Point p2 = vpmap_[source(he, mesh_)];
      Vector_3 normal(p1, p2);

      //construct planes passing through p1 and p2,
      // and orthogonal to e
      Plane_3 plane1(p1, normal);
      Plane_3 plane2(p2, normal);

      CGAL_assertion(//check orientation is consistent
        plane1.orthogonal_vector() * plane2.orthogonal_vector() > 0.);

      //get third points of triangles shared by e
      Point p3 = vpmap_[target(next(he, mesh_), mesh_)];
      Point p4 = vpmap_[target(next(opposite(he, mesh_), mesh_), mesh_)];

      //check whether p3 and p4 are between plane1 and plane2
      return (plane1.oriented_side(p3) != plane2.oriented_side(p3))
         &&  (plane1.oriented_side(p4) != plane2.oriented_side(p4));
    }

    template<typename FaceRange>
    void tag_halfedges_status(FaceRange face_range)
    {
      //tag PATCH,       //h and hopp belong to the patch to be remeshed
      BOOST_FOREACH(face_descriptor f, face_range)
      {
        BOOST_FOREACH(halfedge_descriptor h,
                      halfedges_around_face(halfedge(f, mesh_), mesh_))
        {
          halfedge_status_map_[h] = PATCH;
        }
      }

      //tag PATCH_BORDER,//h belongs to the patch, hopp doesn't
      std::vector<halfedge_descriptor> border_halfedges;
      PMP::get_border(mesh_, face_range, std::back_inserter(border_halfedges));
      BOOST_FOREACH(halfedge_descriptor h, border_halfedges)
      {
        halfedge_status_map_[h] = PATCH_BORDER;
      }

      //tag MESH,        //h and hopp belong to the mesh, not the patch
      //tag MESH_BORDER  //h belongs to the mesh, face(hopp, pmesh) == null_face()
      BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh_))
      {
        //being part of the border of the mesh is predominant
        if (is_border(h, mesh_))
          halfedge_status_map_[h] = MESH_BORDER; //erase previous value if exists
        else
        {
          //h is not border, does not belong to patch, nor to patch border
          if (halfedge_status_map_.find(h) == halfedge_status_map_.end())
            halfedge_status_map_[h] = MESH;
        }
      }
    }

    bool is_on_patch(const halfedge_descriptor& h) const
    {
      return halfedge_status_map_[h] == PATCH;
    }

    bool is_on_patch_border(const halfedge_descriptor& h) const
    {
      return halfedge_status_map_[h] == PATCH_BORDER;
    }
    bool is_on_patch_border(const edge_descriptor& e) const
    {
      return is_on_patch_border(halfedge(e,mesh_))
          || is_on_patch_border(opposite(halfedge(e, mesh_), mesh_));
    }

    bool is_on_border(const halfedge_descriptor& h) const
    {
      CGAL_assertion(is_border(h, mesh_) == halfedge_status_map_[h]);
      return halfedge_status_map_[h] == MESH_BORDER;
    }
    bool is_on_border(const edge_descriptor& e) const
    {
      return is_on_border(halfedge(e, mesh_))
          || is_on_border(opposite(halfedge(e, mesh_), mesh_));
    }

  private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    const AABB_tree* tree_ptr_;
    bool own_tree_;
    Triangle_list input_triangles_;
    std::vector<face_descriptor> patch_;
    std::map<halfedge_descriptor, Halfedge_status> halfedge_status_map_;

  };//end class Incremental_remesher
}//end namespace internal
}//end namespace Polygon_mesh_processing
}//end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
