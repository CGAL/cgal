
#ifndef CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/get_border.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

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

#ifdef CGAL_PMP_REMESHING_DEBUG
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#endif

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

  // A property map 
  template <typename PM, typename FaceRange>
  struct Border_constraint_pmap
  {
    typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PM>::edge_descriptor edge_descriptor;

    std::map<edge_descriptor, bool> border_edges;
    const PM& pmesh_;
  public:
    Border_constraint_pmap(const PM& pmesh, const FaceRange& faces)
      : pmesh_(pmesh)
    {
      std::vector<halfedge_descriptor> border;
      PMP::get_border(pmesh_, faces, std::back_inserter(border));

      BOOST_FOREACH(edge_descriptor e, edges(pmesh_))
        border_edges.insert(std::make_pair(e, false));

      BOOST_FOREACH(halfedge_descriptor h, border)
        border_edges[edge(h, pmesh_)] = true;
    }

    friend bool get(const Border_constraint_pmap<PM, FaceRange>& map,
                    const edge_descriptor& e)
    {
      CGAL_assertion(!map.border_edges.empty());
      typename std::map<edge_descriptor, bool>::const_iterator it
        = map.border_edges.find(e);

      CGAL_assertion(it != map.border_edges.end());
      return it->second;
    }

  };

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
                       , VertexPointMap& vpmap
                       , const bool protect_constraints)
      : mesh_(pmesh)
      , vpmap_(vpmap)
      , own_tree_(true)
      , input_triangles_()
      , halfedge_status_map_()
      , protect_constraints_(protect_constraints)
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
        input_triangles_.push_back(
          Triangle_3(get(vpmap, v1), get(vpmap, v2), get(vpmap, v3)));
      }
      tree_ptr_ = new AABB_tree(input_triangles_.begin(), input_triangles_.end());
    }

    ~Incremental_remesher()
    {
      if (own_tree_)
        delete tree_ptr_;
    }
    
    template<typename FaceRange
           , typename EdgeIsConstrainedMap>
    void init_faces_remeshing(const FaceRange& face_range
                            , const EdgeIsConstrainedMap& ecmap)
    {
      tag_halfedges_status(face_range, ecmap);
    }

    // split edges of edge_range that have their length > high
    template<typename EdgeRange>
    void split_long_edges(const EdgeRange& edge_range,
                          const double& high)
    {
      typedef boost::bimap<
        boost::bimaps::set_of<halfedge_descriptor>,
        boost::bimaps::multiset_of<double, std::greater<double> > >  Boost_bimap;
      typedef typename Boost_bimap::value_type                       long_edge;

      std::cout << "Split long edges (" << high << ")...";
      double sq_high = high*high;

      //collect long edges
      Boost_bimap long_edges;
      BOOST_FOREACH(edge_descriptor e, edge_range)
      {
        double sqlen = sqlength(e);
        if (sqlen > sq_high)
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

        //split edge
        Point refinement_point = this->midpoint(he);
        halfedge_descriptor hnew = CGAL::Euler::split_edge(he, mesh_);
        CGAL_assertion(he == next(hnew, mesh_));
        ++nb_splits;

        //move refinement point
        vertex_descriptor vnew = target(hnew, mesh_);
        put(vpmap_, vnew, refinement_point);

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
          halfedge_descriptor hnew2
            = CGAL::Euler::split_face(hnew, next(next(hnew, mesh_), mesh_), mesh_);
        }

        //do it again on the other side if we're not on boundary
        halfedge_descriptor hnew_opp = opposite(hnew, mesh_);
        if (!is_border(hnew_opp, mesh_))
        {
          halfedge_descriptor hnew2
            = CGAL::Euler::split_face(prev(hnew_opp, mesh_), next(hnew_opp, mesh_), mesh_);
        }
      }
      std::cout << " done (" << nb_splits << " splits)." << std::endl;
#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("0-border_split.off");
#endif
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
      std::cout.flush(); 
      double sq_high = high*high;

      //collect long edges
      Boost_bimap long_edges;
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        if (!is_split_allowed(e))
          continue;
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

        if (protect_constraints_ && !is_longest_on_faces(edge(he, mesh_)))
          continue;

        //split edge
        Point refinement_point = this->midpoint(he);
        halfedge_descriptor hnew = CGAL::Euler::split_edge(he, mesh_);
        CGAL_assertion(he == next(hnew, mesh_));
        ++nb_splits;

        //move refinement point
        vertex_descriptor vnew = target(hnew, mesh_);
        put(vpmap_, vnew, refinement_point);

        //after splitting
        halfedge_descriptor hnew_opp = opposite(hnew, mesh_);
        halfedge_added(hnew, status(he));
        halfedge_added(hnew_opp, status(opposite(he, mesh_)));

        //check sub-edges
        double sqlen_new = 0.25 * sqlen;
        if (sqlen_new > sq_high)
        {
          //if it was more than twice the "long" threshold, insert them
          long_edges.insert(long_edge(hnew,              sqlen_new));
          long_edges.insert(long_edge(next(hnew, mesh_), sqlen_new));
        }

        //insert new edges to keep triangular faces, and update long_edges
        if (!is_on_border(hnew))
        {
          halfedge_descriptor hnew2 = CGAL::Euler::split_face(hnew,
                                                              next(next(hnew, mesh_), mesh_),
                                                              mesh_);
          Halfedge_status snew = (is_on_patch(hnew) || is_on_patch_border(hnew))
            ? PATCH
            : MESH;
          halfedge_added(hnew2,                  snew);
          halfedge_added(opposite(hnew2, mesh_), snew);

          if (snew == PATCH)
          {
            double sql = sqlength(hnew2);
            if (sql > sq_high)
              long_edges.insert(long_edge(hnew2, sql));
          }
        }

        //do it again on the other side if we're not on boundary
        if (!is_on_border(hnew_opp))
        {
          halfedge_descriptor hnew2 = CGAL::Euler::split_face(prev(hnew_opp, mesh_),
                                                              next(hnew_opp, mesh_),
                                                              mesh_);
          Halfedge_status snew = (is_on_patch(hnew_opp) || is_on_patch_border(hnew_opp))
             ? PATCH
            : MESH;
          halfedge_added(hnew2,                  snew);
          halfedge_added(opposite(hnew2, mesh_), snew);

          if (snew == PATCH)
          {
            double sql = sqlength(hnew2);
            if (sql > sq_high)
              long_edges.insert(long_edge(hnew2, sql));
          }
        }
      }
      std::cout << " done ("<< nb_splits << " splits)." << std::endl;

#ifdef CGAL_PMP_REMESHING_DEBUG
      CGAL_expensive_assertion(is_triangle_mesh(mesh_));
      CGAL_assertion(halfedge_status_map_.size() == nb_valid_halfedges());
      debug_status_map();
      debug_patch_border();
      debug_mesh_border();
#endif

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
      std::cout.flush();
      double sq_low = low*low;
      double sq_high = high*high;

      Boost_bimap short_edges;
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        if (!is_collapse_allowed(e))
          continue;
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

        //handle the boundary case :
        //a PATCH_BORDER edge can be collapsed,
        //and an edge incident to PATCH_BORDER can be collapsed,
        //but only if the boundary vertex is kept,
        //so re-insert opposite(he) to collapse it
        if (!is_on_patch(he))
        {
          CGAL_assertion(!protect_constraints_);//is_collapse_allowed returned false
          if (is_on_border(he) || is_on_mesh(he))
          {
            he = opposite(he, mesh_); //he now is PATCH_BORDER
            CGAL_assertion(is_on_patch_border(he));
          }
        }//end if(not on PATCH)

        //let's try to collapse he into vb
        vertex_descriptor va = source(he, mesh_);
        vertex_descriptor vb = target(he, mesh_);

        if (is_on_patch_border(va) && !is_on_patch_border(vb))
        {
          he = opposite(he, mesh_);
          va = source(he, mesh_);
          vb = target(he, mesh_);
          CGAL_assertion(is_on_patch_border(vb) && !is_on_patch_border(va));
        }
        else if (is_on_patch(va) && is_on_patch(vb))
        {
          if(!collapse_does_not_invert_face(he))
          {
            if (collapse_does_not_invert_face(opposite(he, mesh_)))
            {
              he = opposite(he, mesh_);
              va = source(he, mesh_);
              vb = target(he, mesh_);
            }
            else
              continue;//both directions invert a face
          }
          CGAL_assertion(collapse_does_not_invert_face(he));
        }

        CGAL_assertion(is_collapse_allowed(edge(he, mesh_)));

        if (degree(va, mesh_) < 3
          || degree(vb, mesh_) < 3
          || !CGAL::Euler::satisfies_link_condition(edge(he, mesh_), mesh_))//necessary to collapse
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

          //before collapse
          halfedge_descriptor ep_p  = prev(opposite(he, mesh_), mesh_);
          halfedge_descriptor epo_p = opposite(ep_p, mesh_);
          halfedge_descriptor en    = next(he, mesh_);
          halfedge_descriptor en_p  = next(opposite(he, mesh_), mesh_);
          Halfedge_status s_ep_p    = status(ep_p);
          Halfedge_status s_en_p    = status(en_p);
          Halfedge_status s_epo_p   = status(epo_p);
          Halfedge_status s_ep      = status(prev(he, mesh_));
          Halfedge_status s_epo     = status(opposite(prev(he, mesh_), mesh_));

          bool mesh_border_case = is_on_border(opposite(he, mesh_));
          if (!mesh_border_case)
            halfedge_and_opp_removed(prev(opposite(he, mesh_), mesh_));
          halfedge_and_opp_removed(he);
          halfedge_and_opp_removed(prev(he, mesh_));

          //perform collapse
          Point target_point = get(vpmap_, vb);

          vertex_descriptor vkept = CGAL::Euler::collapse_edge(edge(he, mesh_), mesh_);
          put(vpmap_, vkept, target_point);
          ++nb_collapses;

#ifdef CGAL_PMP_REMESHING_DEBUG
          debug_normals(vkept);
          //PMP::remove_degenerate_faces(mesh_/*todo : add named parameters*/);
          //debug_self_intersections(vkept);
#endif

          // merge halfedge_status to keep the more important on both sides
          merge_status(en, s_epo, s_ep);
          if (!mesh_border_case)
            merge_status(en_p, s_epo_p, s_ep_p);

#ifdef CGAL_PMP_REMESHING_DEBUG
          unsigned int nbb = nb_valid_halfedges();
          CGAL_assertion(nbb == halfedge_status_map_.size());
          CGAL_assertion(source(en, mesh_) == source(en_p, mesh_));
          debug_status_map();
          debug_patch_border();
#endif

          //insert new/remaining short edges
          BOOST_FOREACH(halfedge_descriptor ht, halfedges_around_target(vkept, mesh_))
          {
            if (!is_collapse_allowed(edge(ht, mesh_)))
              continue;
            double sqlen = sqlength(ht);
            if (sqlen < sq_low)
              short_edges.insert(short_edge(ht, sqlen));
          }
        }//end if(collapse_ok)
      }

      PMP::remove_degenerate_faces(mesh_/*todo : add named parameters*/);

      std::cout << " done (" << nb_collapses << " collapses)." << std::endl;

#ifdef CGAL_DUMP_REMESHING_STEPS
      dump("2-edge_collapse.off");
#endif

#ifdef CGAL_PMP_REMESHING_DEBUG
      CGAL_assertion(nb_valid_halfedges() == halfedge_status_map_.size());
      CGAL_expensive_assertion(is_triangle_mesh(mesh_));
      debug_status_map();
      debug_patch_border();
      debug_mesh_border();
      debug_self_intersections();
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
      std::cout.flush(); 
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

        CGAL_assertion_code(Halfedge_status s1 = status(he));
        CGAL_assertion_code(Halfedge_status s1o = status(opposite(he, mesh_)));

        CGAL::Euler::flip_edge(he, mesh_);
        ++nb_flips;

        CGAL_assertion_code(Halfedge_status s2 = status(he));
        CGAL_assertion_code(Halfedge_status s2o = status(opposite(he, mesh_)));
        CGAL_assertion(s1 == s2   && s1 == PATCH);
        CGAL_assertion(s1o == s2o && s1o == PATCH);
        CGAL_assertion(nb_valid_halfedges() == halfedge_status_map_.size());
        CGAL_assertion(!is_border(he, mesh_));

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

          CGAL_assertion_code(Halfedge_status s3 = status(he));
          CGAL_assertion(s1 == s3);
          CGAL_assertion(!is_border(he, mesh_));
          CGAL_assertion(
               (va == source(he, mesh_) && vb == target(he, mesh_))
            || (vb == source(he, mesh_) && va == target(he, mesh_)));
        }
      }

      std::setprecision(17);
      dump("after-edge-flips.off");

      PMP::remove_degenerate_faces(mesh_
        , PMP::parameters::vertex_point_map(vpmap_)
        .geom_traits(GeomTraits()));

      std::cout << "done. ("<< nb_flips << " flips)" << std::endl;

#ifdef CGAL_PMP_REMESHING_DEBUG
      CGAL_assertion(nb_valid_halfedges() == halfedge_status_map_.size());
      debug_status_map();
#endif

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
      std::cout.flush();

      //todo : use boost::vector_property_map to improve computing time
      typedef std::map<vertex_descriptor, Vector_3> VNormalsMap;
      VNormalsMap vnormals;
      boost::associative_property_map<VNormalsMap> propmap_normals(vnormals);
      BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
      {
        if (!is_on_patch(v))
          continue;
        Vector_3 vn = PMP::compute_vertex_normal(v, mesh_
                            , PMP::parameters::vertex_point_map(vpmap_)
                            .geom_traits(GeomTraits()));
        put(propmap_normals, v, vn);
      }

      // at each vertex, compute barycenter of neighbors
      std::map<vertex_descriptor, Point> barycenters;
      BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
      {
        if (!is_on_patch(v))
          continue;
        Vector_3 move = CGAL::NULL_VECTOR;
        unsigned int star_size = 0;
        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, mesh_))
        {
          move = move + Vector_3(get(vpmap_, v), get(vpmap_, source(h, mesh_)));
          ++star_size;
        }
        CGAL_assertion(star_size > 0);
        move = (1. / (double)star_size) * move;

        barycenters[v] = get(vpmap_, v) + move;
      }

      // compute moves
      typedef typename std::map<vertex_descriptor, Point>::value_type VP_pair;
      std::map<vertex_descriptor, Point> new_locations;
      BOOST_FOREACH(const VP_pair& vp, barycenters)
      {
        vertex_descriptor v = vp.first;
        Point pv = get(vpmap_, v);
        Vector_3 nv = boost::get(propmap_normals, v);
        Point qv = vp.second; //barycenter at v

        new_locations[v] = qv + (nv * Vector_3(qv, pv)) * nv;
      }

      // perform moves
      BOOST_FOREACH(const VP_pair& vp, new_locations)
      {
        put(vpmap_, vp.first, vp.second);
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
      std::cout.flush();

      BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
      {
        if (!is_on_patch(v))
          continue;
        put(vpmap_, v, tree_ptr_->closest_point(get(vpmap_, v)));
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
      return CGAL::squared_distance(get(vpmap_, v1), get(vpmap_, v2));
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

    Point midpoint(const halfedge_descriptor& he) const
    {
      Point p1 = get(vpmap_, target(he, mesh_));
      Point p2 = get(vpmap_, source(he, mesh_));
      return CGAL::midpoint(p1, p2);
    }

    void dump(const char* filename) const
    {
      std::ofstream out(filename);
//      out << mesh_;
      out.close();
    }

    int valence(const vertex_descriptor& v) const
    {
      return static_cast<int>(degree(v, mesh_));
    }

    int target_valence(const vertex_descriptor& v) const
    {
      return (is_border(v, mesh_)) ? 4 : 6;
    }

    bool is_flip_allowed(const edge_descriptor& e) const
    {
      //only the patch, and non-border edges are
      //allowed to be flipped
      halfedge_descriptor he = halfedge(e, mesh_);
      if (!is_on_patch(he))
        return false;

      Point p1 = get(vpmap_, target(he, mesh_));
      Point p2 = get(vpmap_, source(he, mesh_));
      Vector_3 normal(p1, p2);

      //construct planes passing through p1 and p2,
      // and orthogonal to e
      Plane_3 plane1(p1, normal);
      Plane_3 plane2(p2, normal);

      CGAL_assertion(//check orientation is consistent
        plane1.orthogonal_vector() * plane2.orthogonal_vector() > 0.);

      //get third points of triangles shared by e
      Point p3 = get(vpmap_, target(next(he, mesh_), mesh_));
      Point p4 = get(vpmap_, target(next(opposite(he, mesh_), mesh_), mesh_));

      //check whether p3 and p4 are between plane1 and plane2
      return (plane1.oriented_side(p3) != plane2.oriented_side(p3))
         &&  (plane1.oriented_side(p4) != plane2.oriented_side(p4));
    }

    bool is_longest_on_faces(const edge_descriptor& e) const
    {
      halfedge_descriptor h = halfedge(e, mesh_);
      halfedge_descriptor hopp = opposite(h, mesh_);

      //check whether h is the longest edge in its associated face
      //overwise refinement will go for an endless loop
      double sqh = sqlength(h);
      return sqh >= sqlength(next(h, mesh_))
          && sqh >= sqlength(next(next(h, mesh_), mesh_))
          //do the same for hopp
          && sqh >= sqlength(next(hopp, mesh_))
          && sqh >= sqlength(next(next(hopp, mesh_), mesh_));
    }

    bool is_split_allowed(const edge_descriptor& e) const
    {
      halfedge_descriptor h = halfedge(e, mesh_);
      halfedge_descriptor hopp = opposite(h, mesh_);

      if (protect_constraints_)
      {
        return is_on_patch(h); //PATCH are the only splittable edges
      }
      else //allow splitting constraints
      {
        if (is_on_mesh(h) && is_on_mesh(hopp))
          return false;
        else if (is_on_mesh(h) && is_on_border(hopp))
          return false;
        else if (is_on_mesh(hopp) && is_on_border(h))
          return false;
        else
          return true;
      }
    }

    bool is_collapse_allowed(const edge_descriptor& e) const
    {
      halfedge_descriptor he = halfedge(e, mesh_);
      halfedge_descriptor hopp = opposite(he, mesh_);

      if (is_on_patch(he)) //hopp is also on patch
        return true;
      else if (is_on_patch_border(he) || is_on_patch_border(hopp))
        return !protect_constraints_;//allowed only when no protection
      else
        return false;
    }

    bool collapse_does_not_invert_face(const halfedge_descriptor& h) const
    {
      vertex_descriptor vs = source(h, mesh_);
      vertex_descriptor vt = target(h, mesh_);
      
      //backup source point
      Point ps = get(vpmap_, vs);
      //move source at target
      put(vpmap_, vs, get(vpmap_, vt));

      //collect normals to faces around vs AND vt
      //vertices are at the same location, but connectivity is still be same,
      //with plenty of degenerate triangles (which are common to both stars)
      std::vector<Vector_3> normals;
      BOOST_FOREACH(halfedge_descriptor hd,
                    halfedges_around_target(h, mesh_))
      {
        Vector_3 n = compute_normal(face(hd, mesh_));
        if (n != CGAL::NULL_VECTOR)
          normals.push_back(n);
      }
      BOOST_FOREACH(halfedge_descriptor hd,
                    halfedges_around_target(opposite(h, mesh_), mesh_))
      {
        Vector_3 n = compute_normal(face(hd, mesh_));
        if (n != CGAL::NULL_VECTOR)
          normals.push_back(n);
      }

      //check all normals have same orientation
      for(std::size_t i = 1; i < normals.size(); ++i)/*start at 1 on purpose*/
      {
        if (normals[i-1] * normals[i] <= 0.)
        {
          //restore position
          put(vpmap_, vs, ps);
          return false;
        }
      }
      //restore position
      put(vpmap_, vs, ps);
      return true;
    }

    Vector_3 compute_normal(const face_descriptor& f) const
    {
      halfedge_descriptor hd = halfedge(f, mesh_);
      typename GeomTraits::Triangle_3
        tr(get(vpmap_, target(hd, mesh_)),
           get(vpmap_, target(next(hd, mesh_), mesh_)),
           get(vpmap_, target(next(next(hd, mesh_), mesh_), mesh_)));

      if (tr.is_degenerate())
        return CGAL::NULL_VECTOR;
      else
        return PMP::compute_face_normal(f, mesh_);
    }

    template<typename FaceRange, typename EdgeIsConstrainedMap>
    void tag_halfedges_status(const FaceRange& face_range
                            , const EdgeIsConstrainedMap& ecmap)
    {
      //tag MESH,        //h and hopp belong to the mesh, not the patch
      //tag MESH_BORDER  //h belongs to the mesh, face(hopp, pmesh) == null_face()
      BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh_))
      {
        //being part of the border of the mesh is predominant
        if (is_border(h, mesh_))
          halfedge_status_map_[h] = MESH_BORDER; //erase previous value if exists
        else
          halfedge_status_map_[h] = MESH;
      }

      //tag PATCH,       //h and hopp belong to the patch to be remeshed
      BOOST_FOREACH(face_descriptor f, face_range)
      {
        BOOST_FOREACH(halfedge_descriptor h,
                      halfedges_around_face(halfedge(f, mesh_), mesh_))
        {
          halfedge_status_map_[h] = PATCH;
        }
      }

      //override the border of PATCH
      //tag PATCH_BORDER,//h belongs to the patch, hopp doesn't
      BOOST_FOREACH(edge_descriptor e, edges(mesh_))
      {
        if (get(ecmap, e))
        {
          //deal with h and hopp for borders that are sharp edges to be preserved
          halfedge_descriptor h = halfedge(e, mesh_);
          if (halfedge_status_map_[h] == PATCH)
            halfedge_status_map_[h] = PATCH_BORDER;
          halfedge_descriptor hopp = opposite(halfedge(e, mesh_), mesh_);
          if (halfedge_status_map_[hopp] == PATCH)
            halfedge_status_map_[hopp] = PATCH_BORDER;
        }
      }

#ifdef CGAL_PMP_REMESHING_DEBUG
      CGAL_assertion(halfedge_status_map_.size() == nb_valid_halfedges());
      debug_patch_border();
#endif
    }

    Halfedge_status status(const halfedge_descriptor& h) const
    {
      typename std::map < halfedge_descriptor, Halfedge_status >::const_iterator
        it = halfedge_status_map_.find(h);
      CGAL_assertion(it != halfedge_status_map_.end());
      return it->second;
    }

    void merge_status(const halfedge_descriptor& en,
                      const Halfedge_status& s_epo,
                      const Halfedge_status& s_ep)
    {
      CGAL_assertion(halfedge_status_map_.find(en) != halfedge_status_map_.end());

      //get missing data
      halfedge_descriptor eno = opposite(en, mesh_);
      Halfedge_status s_en = status(en);
      Halfedge_status s_eno = status(eno);

      if(s_epo == MESH_BORDER
        || s_ep == MESH_BORDER
        || s_epo == PATCH_BORDER
        || s_ep == PATCH_BORDER)
      {
        halfedge_status_map_[en]  = s_epo;
        halfedge_status_map_[eno] = s_ep;
      }
      // else keep current status for en and eno
    }

    bool is_on_patch(const halfedge_descriptor& h) const
    {
      bool res =(status(h) == PATCH);
      CGAL_assertion(res == (status(opposite(h, mesh_)) == PATCH));
      return res;
    }

    bool is_on_patch(const face_descriptor& f) const
    {
      BOOST_FOREACH(halfedge_descriptor h,
                    halfedges_around_face(halfedge(f, mesh_), mesh_))
      {
        if (is_on_patch(h) || is_on_patch_border(h))
          return true;
      }
      return false;
    }

    bool is_on_patch(const vertex_descriptor& v) const
    {
      BOOST_FOREACH(halfedge_descriptor h,
                    halfedges_around_target(v, mesh_))
      {
        if (!is_on_patch(h))
          return false;
      }
      return true;
    }

    bool is_on_patch_border(const halfedge_descriptor& h) const
    {
      bool res = (status(h) == PATCH_BORDER);
      if (res)
      {
        CGAL_assertion_code(Halfedge_status hs = status(opposite(h, mesh_)));
        CGAL_assertion(hs == MESH_BORDER
                    || hs == MESH
                    || hs == PATCH_BORDER);//when 2 incident patches are remeshed
      }
      return res;
    }
    bool is_on_patch_border(const edge_descriptor& e) const
    {
      return is_on_patch_border(halfedge(e,mesh_))
          || is_on_patch_border(opposite(halfedge(e, mesh_), mesh_));
    }
    bool is_on_patch_border(const vertex_descriptor& v) const
    {
      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, mesh_))
      {
        if (is_on_patch_border(h) || is_on_patch_border(opposite(h, mesh_)))
          return true;
      }
      return false;
    }

    bool is_on_border(const halfedge_descriptor& h) const
    {
      bool res = (status(h) == MESH_BORDER);
      CGAL_assertion(res == is_border(h, mesh_));
      CGAL_assertion(res == is_border(next(h, mesh_), mesh_));
      return res;
    }

    bool is_on_border(const edge_descriptor& e) const
    {
      return is_on_border(halfedge(e, mesh_))
          || is_on_border(opposite(halfedge(e, mesh_), mesh_));
    }

    bool is_on_mesh(const halfedge_descriptor& h) const
    {
      return status(h) == MESH;
    }

    void halfedge_added(const halfedge_descriptor& h,
                        const Halfedge_status& s)
    {
      halfedge_status_map_.insert(std::make_pair(h, s));
    }

    void halfedge_and_opp_removed(const halfedge_descriptor& h)
    {
      halfedge_status_map_.erase(h);
      halfedge_status_map_.erase(opposite(h, mesh_));
    }

    unsigned int nb_valid_halfedges() const
    {
      unsigned int nb = 0;
      BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh_))
        ++nb;
      return nb;
    }

    void debug_status_map() const
    {
      typedef typename std::map<halfedge_descriptor, Halfedge_status>::value_type
        HD_pair;
      BOOST_FOREACH(const HD_pair& hs, halfedge_status_map_)
      {
        bool b1 = is_on_patch(hs.first);
        bool b2 = is_on_patch_border(hs.first);
        bool b3 = is_on_mesh(hs.first);
        bool b4 = is_on_border(hs.first);
      }
    }

    void debug_patch_border() const
    {
      std::map<vertex_descriptor, unsigned int> patch_border;
      typedef typename std::map<halfedge_descriptor, Halfedge_status>::value_type
        HD_pair;
      BOOST_FOREACH(const HD_pair& hs, halfedge_status_map_)
      {
        if (is_on_patch_border(hs.first))
        {
          if (patch_border.find(target(hs.first, mesh_)) != patch_border.end())
            patch_border[target(hs.first, mesh_)]++;
          else
            patch_border[target(hs.first, mesh_)] = 1;

          if (patch_border.find(source(hs.first, mesh_)) != patch_border.end())
            patch_border[source(hs.first, mesh_)]++;
          else
            patch_border[source(hs.first, mesh_)] = 1;
        }
      }
      //check we found each vertex exactly twice
      typedef typename std::map<vertex_descriptor, unsigned int>::value_type
        V_pair;
      BOOST_FOREACH(const V_pair& v, patch_border)
        CGAL_assertion(v.second == 2);
    }

    void debug_self_intersections() const
    {
      std::cout << "Test self intersections...";
      std::vector<std::pair<face_descriptor, face_descriptor> > facets;
      PMP::does_self_intersect(
        mesh_,
        std::back_inserter(facets),
        PMP::parameters::vertex_point_map(vpmap_));
      CGAL_assertion(facets.empty());
      std::cout << "done." << std::endl;
    }

    void debug_self_intersections(const vertex_descriptor& v) const
    {
      std::cout << "Test self intersections...";
      std::vector<std::pair<face_descriptor, face_descriptor> > facets;
      PMP::does_self_intersect(
        faces_around_target(halfedge(v, mesh_), mesh_),
        mesh_,
        std::back_inserter(facets),
        PMP::parameters::vertex_point_map(vpmap_));
      CGAL_assertion(facets.empty());
      std::cout << "done." << std::endl;
    }

    void debug_normals(const vertex_descriptor& v) const
    {
      if (!is_on_patch(v))
        return;//not much to say if we are on a boundary/sharp edge

      std::vector<Vector_3> normals;
      BOOST_FOREACH(halfedge_descriptor hd,
                    halfedges_around_target(halfedge(v, mesh_), mesh_))
      {
        Vector_3 n = compute_normal(face(hd, mesh_));
        if (n != CGAL::NULL_VECTOR)
          normals.push_back(n);
      }
      //check all normals have same orientation
      for (std::size_t i = 1; i < normals.size(); ++i)/*start at 1 on purpose*/
        CGAL_assertion(normals[i - 1] * normals[i] > 0.);
    }

    void debug_mesh_border() const
    {
      std::map<vertex_descriptor, unsigned int> mesh_border;
      typedef typename std::map<halfedge_descriptor, Halfedge_status>::value_type
        HD_pair;
      BOOST_FOREACH(const HD_pair& hs, halfedge_status_map_)
      {
        if (is_on_border(hs.first))
        {
          if (mesh_border.find(target(hs.first, mesh_)) != mesh_border.end())
            mesh_border[target(hs.first, mesh_)]++;
          else
            mesh_border[target(hs.first, mesh_)] = 1;

          if (mesh_border.find(source(hs.first, mesh_)) != mesh_border.end())
            mesh_border[source(hs.first, mesh_)]++;
          else
            mesh_border[source(hs.first, mesh_)] = 1;
        }
      }
      //check we found each vertex exactly twice
      typedef typename std::map<vertex_descriptor, unsigned int>::value_type
        V_pair;
      BOOST_FOREACH(const V_pair& v, mesh_border)
        CGAL_assertion(v.second == 2);
    }

  private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    const AABB_tree* tree_ptr_;
    bool own_tree_;
    Triangle_list input_triangles_;
    std::map<halfedge_descriptor, Halfedge_status> halfedge_status_map_;
    bool protect_constraints_;

  };//end class Incremental_remesher
}//end namespace internal
}//end namespace Polygon_mesh_processing
}//end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
