//#define CGAL_PMP_REMESHING_VERBOSE
//#define CGAL_PMP_REMESHING_DEBUG
//#define CGAL_PMP_REMESHING_VERY_VERBOSE
//#define CGAL_PMP_REMESHING_VERBOSE_PROGRESS

#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>


#include "Scene_surface_mesh_item.h"

#include "Scene_polyhedron_selection_item.h"

#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/utility.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_set.hpp>
#include <CGAL/property_map.h>

#include <QElapsedTimer>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QDialog>
#include <QtPlugin>
#include <QMessageBox>

#include <vector>
#include <algorithm>
#include <queue>
#include <sstream>
#include <cmath>

#ifdef CGAL_LINKED_WITH_TBB
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/partitioner.h"
#endif

#include "ui_Isotropic_remeshing_dialog.h"


typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef Scene_facegraph_item::Face_graph FaceGraph;
typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;


// give a halfedge and a target edge length, put in `out` points
// which the edge equally spaced such that splitting the edge
// using the sequence of points make the edges shorter than
// `target_length`
template <class TriangleMesh, class PointPMap, class PointOutputIterator>
PointOutputIterator
sample_edge(
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
  TriangleMesh& triangle_mesh,
  double target_length,
  const PointPMap& pmap,
  PointOutputIterator out)
{
  typedef typename boost::property_traits<PointPMap>::value_type Point_3;
  typedef typename CGAL::Kernel_traits<Point_3>::Kernel::Vector_3 Vector_3;
  typename boost::property_traits<PointPMap>::reference src=get(pmap, source(hd,triangle_mesh) );
  typename boost::property_traits<PointPMap>::reference tgt=get(pmap, target(hd,triangle_mesh) );

  double length = std::sqrt( CGAL::squared_distance(src, tgt) );
  if ( length <= target_length ) return out;

  double nb_points = std::floor( length / target_length );
  Vector_3 unit = (tgt-src) / (nb_points+1);

  for(double i=0; i<nb_points; ++i)
    *out++=src+unit*(i+1);

  return out;
}

// given a set of points that are expected to be on an edge, split
// that edge and retriangulate the face incident to the edge
// Points are sorted so that they are sorted from the source to the target
// of the edge (the sequence does not contains edge endpoints)
template <class TriangleMesh, class PointPMap, class PointRange, class EdgeOutputIterator>
EdgeOutputIterator
split_identical_edges(
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
  TriangleMesh& tm,
  const PointPMap& pmap,
  const PointRange& points,
  EdgeOutputIterator out)
{
  typedef typename PointRange::value_type Point_3;
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  for(const Point_3& p : points)
  {
    // split the edge
    halfedge_descriptor new_hd=CGAL::Euler::split_edge(hd,tm);
    // set the vertex point
    put(pmap, target(new_hd, tm), p);
    *out++=edge(new_hd, tm);
  }
  *out++=edge(hd, tm);
  return out;
}

// HedgeRange is expected to be a range with value type being
// std::pair<halfedge_descriptor, TriangleMesh*>
// Given a set of halfedges representing different edges
// but with identical endpoints, and a target edge length
// we split all edges identically so that subedges are
// or length <= length
template <class HedgeRange, class Edges_to_protect>
void split_long_duplicated_edge(const HedgeRange& hedge_range,
                                double target_length,
                                Edges_to_protect& edges_to_protect)
{
  typedef typename HedgeRange::value_type Pair;
  typedef typename Pair::first_type halfedge_descriptor;
  typedef typename boost::remove_pointer<
    typename Pair::second_type>::type TriangleMesh;
  typedef typename boost::property_map<TriangleMesh,
    CGAL::vertex_point_t>::type PointPMap;
  typedef typename boost::property_traits<PointPMap>::value_type Point_3;

  if (hedge_range.empty()) return;

  const Pair& p = *hedge_range.begin();
  PointPMap pmap = get(boost::vertex_point, *p.second);

  std::vector<Point_3> points;
  halfedge_descriptor hd = p.first;

  // collect points to be add inside the edges
  sample_edge(hd, *p.second, target_length, pmap, std::back_inserter(points) );

  CGAL_assertion_code(Point_3 src = get(pmap, source(hd, *p.second));)
  CGAL_assertion_code(Point_3 tgt = get(pmap, target(hd, *p.second));)

  // split the edges and collect faces to triangulate
  for(const Pair& h_and_p : hedge_range)
  {
    halfedge_descriptor hc=h_and_p.first;
    TriangleMesh* polyc = h_and_p.second;

    PointPMap pmap_2 = get(boost::vertex_point, *polyc);
    //make sure halfedge are consistently oriented
    CGAL_assertion( get(pmap_2, source(hc, *polyc)) == src );
    CGAL_assertion( get(pmap_2, target(hc, *polyc)) == tgt );

    typedef typename Edges_to_protect::value_type::second_type Edge_set;
    Edge_set& edge_set = edges_to_protect[polyc];

    // now split the halfedge and incident faces
    split_identical_edges(hc,*polyc,pmap_2, points,
                          std::inserter( edge_set, edge_set.begin()));
  }
}

using namespace CGAL::Three;
class Polyhedron_demo_isotropic_remeshing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "isotropic_remeshing_plugin.json")

  typedef boost::graph_traits<FaceGraph>::edge_descriptor edge_descriptor;
  typedef boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;

  typedef boost::unordered_set<edge_descriptor>    Edge_set;
  typedef Scene_polyhedron_selection_item::Is_constrained_map<Edge_set> Edge_constrained_pmap;

  struct Visitor
  {
    typedef typename Scene_polyhedron_selection_item::Selection_set_facet Container;
    Container& faces;

    Visitor(Container& container)
      : faces(container)
    {}

    void before_subface_creations(face_descriptor fd)
    {
      Container::iterator it = faces.find(fd);
      faces.erase(it);
    }
    void after_subface_created(face_descriptor fd)
    {
      faces.insert(fd);
    }
    void after_subface_creations(){}
  };

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionIsotropicRemeshing_ = new QAction("Isotropic Remeshing", mw);
    actionIsotropicRemeshing_->setProperty("subMenuName", "Polygon Mesh Processing");
    if (actionIsotropicRemeshing_) {
      connect(actionIsotropicRemeshing_, SIGNAL(triggered()),
        this, SLOT(isotropic_remeshing()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionIsotropicRemeshing_;
  }

  bool applicable(QAction*) const
  {
    if (scene->selectionIndices().size() == 1)
    {
    return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()))
    || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
    }

    bool ok(true), found_poly(false);

    Q_FOREACH(int index, scene->selectionIndices())
    {
      if (!qobject_cast<Scene_facegraph_item*>(scene->item(index)))
        ok = false;
      else
        found_poly=true;
    }
    return ok && found_poly;
  }

  typedef boost::property_map<FaceGraph, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;

  void detect_and_split_duplicates(std::vector<Scene_facegraph_item*>& selection,
                                   std::map<FaceGraph*,Edge_set>& edges_to_protect,
                                   double target_length)
  {
    typedef EPICK::Point_3 Point_3;
    typedef std::pair<Point_3,Point_3> Segment_3;

    typedef std::map< Segment_3,
                      std::vector< std::pair<halfedge_descriptor, FaceGraph*> > > MapType;
    typedef boost::property_map<FaceGraph,
      CGAL::vertex_point_t>::type PointPMap;
    MapType duplicated_edges;

    for(Scene_facegraph_item* poly_item : selection){
      FaceGraph& pmesh = *poly_item->polyhedron();
      PointPMap pmap = get(boost::vertex_point, pmesh);
      for(edge_descriptor ed : edges(pmesh)){
        halfedge_descriptor hd = halfedge(ed,pmesh);
        Point_3 p = get(pmap, source(hd,pmesh)), q = get(pmap, target(hd,pmesh));
        Segment_3 s = CGAL::make_sorted_pair(p,q);
        if (s.first==q) hd=opposite(hd,pmesh); // make sure the halfedges are consistently oriented

        duplicated_edges[s].push_back( std::make_pair(hd,&pmesh) );
      }
    }

    // consistently split duplicate edges and triangulate incident faces
    typedef std::pair<face_descriptor, FaceGraph*> Face_and_poly;
    std::set< Face_and_poly > faces_to_triangulate;
    for(const MapType::value_type& p : duplicated_edges)
      if (p.second.size()>1){
        //collect faces to retriangulate
        typedef std::pair<halfedge_descriptor, FaceGraph*> Pair_type;
        for(const Pair_type& h_and_p : p.second)
        {
          halfedge_descriptor hc=h_and_p.first;
          FaceGraph* polyc = h_and_p.second;

          if ( !is_border(hc, *polyc) )
            faces_to_triangulate.insert( Face_and_poly(face(hc,*polyc), polyc) );
          if ( !is_border(opposite(hc, *polyc), *polyc) )
            faces_to_triangulate.insert(
              Face_and_poly(face(opposite(hc, *polyc),*polyc), polyc) );
        }
        // split the edges
        split_long_duplicated_edge(p.second, target_length, edges_to_protect);
      }
    // now retriangulate
    namespace PMP=CGAL::Polygon_mesh_processing;
    for(Face_and_poly f_and_p : faces_to_triangulate)
      PMP::triangulate_face(f_and_p.first, *f_and_p.second);
  }

  void do_split_edges(Scene_polyhedron_selection_item* selection_item,
                      SMesh& pmesh,
                      double target_length)
  {
    std::vector<edge_descriptor> p_edges;
    for(edge_descriptor e : edges(pmesh))
    {
      if(get(selection_item->constrained_edges_pmap(), e))
        p_edges.push_back(e);
    }
    for(face_descriptor f : selection_item->selected_facets)
    {
      for(halfedge_descriptor he : halfedges_around_face(halfedge(f, pmesh), pmesh))
      {
        if (selection_item->selected_facets.find(face(opposite(he, pmesh), pmesh))
            == selection_item->selected_facets.end())
          p_edges.push_back(edge(he, pmesh));
      }
    }
    if (!p_edges.empty())
      CGAL::Polygon_mesh_processing::split_long_edges(
            p_edges
            , target_length
            , *selection_item->polyhedron()
            , PMP::parameters::geom_traits(EPICK())
            .edge_is_constrained_map(selection_item->constrained_edges_pmap()));
    else
      std::cout << "No selected or boundary edges to be split" << std::endl;
  }

public Q_SLOTS:
  void isotropic_remeshing()
  {
    if (scene->selectionIndices().size() > 1)
    {
      isotropic_remeshing_of_several_polyhedra();
      return;
    }
    const Scene_interface::Item_id index = scene->mainSelectionIndex();

    Scene_facegraph_item* poly_item =
      qobject_cast<Scene_facegraph_item*>(scene->item(index));

    Scene_polyhedron_selection_item* selection_item =
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

    if (poly_item || selection_item)
    {
      if(selection_item && selection_item->selected_edges.empty())
      {
        QMessageBox::warning(mw, "Empty Edges", "There are no selected edges. Aborting.");
        return;
      }
      // Create dialog box
      QDialog dialog(mw);
      Ui::Isotropic_remeshing_dialog ui
        = remeshing_dialog(&dialog, poly_item, selection_item);

      // Get values
      int i = dialog.exec();
      if (i == QDialog::Rejected)
      {
        std::cout << "Remeshing aborted" << std::endl;
        return;
      }
      bool edges_only = ui.splitEdgesOnly_checkbox->isChecked();
      bool preserve_duplicates = ui.preserveDuplicates_checkbox->isChecked();
      double target_length = ui.edgeLength_dspinbox->value();
      unsigned int nb_iter = ui.nbIterations_spinbox->value();
      unsigned int nb_smooth = ui.nbSmoothing_spinbox->value();
      bool protect = ui.protect_checkbox->isChecked();
      bool smooth_features = ui.smooth1D_checkbox->isChecked();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QElapsedTimer time;
      time.start();

      typedef boost::graph_traits<FaceGraph>::edge_descriptor edge_descriptor;
      typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;

      FaceGraph& pmesh = (poly_item != NULL)
        ? *poly_item->polyhedron()
        : *selection_item->polyhedron();

     Patch_id_pmap fpmap = get(CGAL::face_patch_id_t<int>(), pmesh);
     bool fpmap_valid = false;
     {
       for(face_descriptor f : faces(pmesh))
       {
         if (get(fpmap, f) != 1)
         {
           fpmap_valid = true;
           break;/*1 is the default value for both Surface_mesh and Polyhedron*/
         }
       }
     }

      if (selection_item)
      {
        if (edges_only)
        {
          do_split_edges(selection_item, pmesh, target_length);
        }
        else //not edges_only
        {
            if(protect &&
               !CGAL::Polygon_mesh_processing::internal::constraints_are_short_enough(
                 *selection_item->polyhedron(),
                 selection_item->constrained_edges_pmap(),
                 get(CGAL::vertex_point, *selection_item->polyhedron()),
                 CGAL::Constant_property_map<face_descriptor, std::size_t>(1),
                 4. / 3. * target_length))
            {
              QApplication::restoreOverrideCursor();
              //If facets are selected, splitting edges will add facets that won't be selected, and it will mess up the rest.
              //If there is only edges, it will work fine because new edges are dealt with in the code, so we can directly
              //split and continue.
              // Possibility todo: check if the barycenter of a new face is inside an old selected face to
              //select it again.
              if(!selection_item->selected_facets.empty())
              {
                QMessageBox::warning(mw, tr("Error"),
                                      tr("Isotropic remeshing : protect_constraints cannot be set to"
                                         " true with constraints larger than 4/3 * target_edge_length."
                                         " Aborting."));
                return;
              }
              else if(QMessageBox::question(mw, tr("Error"),
                                            tr("Isotropic remeshing : protect_constraints cannot be set to"
                                               " true with constraints larger than 4/3 * target_edge_length."
                                               " Do you wish to split the constrained edges ?")) !=
                      QMessageBox::Yes)
              {
                return;
              }
              else
              {
                do_split_edges(selection_item, pmesh, target_length);
              }
            }

            if (selection_item->selected_facets.empty() && !selection_item->isEmpty())
            {
              if (!CGAL::is_triangle_mesh(pmesh))
              {
                QApplication::restoreOverrideCursor();
                if (QMessageBox::Ok ==
                    QMessageBox::question(mw, tr("Error - Triangulate Faces?"),
                      tr("The input mesh is not a triangulated surface mesh.\n"
                         "Do you wish to triangulate faces first, or cancel remeshing ?"),
                      (QMessageBox::Ok | QMessageBox::Cancel),
                      QMessageBox::Ok))
                {
                  QApplication::setOverrideCursor(Qt::WaitCursor);
                  CGAL::Polygon_mesh_processing::triangulate_faces(pmesh);
                }
                else
                {
                  return;
                }
              }

              if (fpmap_valid)
                CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(*selection_item->polyhedron())
                   , target_length
                   , *selection_item->polyhedron()
                   , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
                   .protect_constraints(protect)
                   .edge_is_constrained_map(selection_item->constrained_edges_pmap())
                   .relax_constraints(smooth_features)
                   .number_of_relaxation_steps(nb_smooth)
                   .vertex_is_constrained_map(selection_item->constrained_vertices_pmap())
                   .face_patch_map(fpmap));
              else
                CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(*selection_item->polyhedron())
                   , target_length
                   , *selection_item->polyhedron()
                   , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
                   .protect_constraints(protect)
                   .edge_is_constrained_map(selection_item->constrained_edges_pmap())
                   .relax_constraints(smooth_features)
                   .number_of_relaxation_steps(nb_smooth)
                   .vertex_is_constrained_map(selection_item->constrained_vertices_pmap())
                                                                   );
            }
            else //selected_facets not empty
            {
              for (auto f : selection_item->selected_facets)
              {
                if (!CGAL::is_triangle(halfedge(f, pmesh), pmesh))
                {
                  QApplication::restoreOverrideCursor();
                  if(QMessageBox::Ok ==
                     QMessageBox::question(mw, tr("Error - Triangulate Faces?"),
                       tr("The input faces selected for remeshing are not all triangle faces.\n"
                          "Do you wish to triangulate faces first, or cancel remeshing ?"),
                       (QMessageBox::Ok | QMessageBox::Cancel),
                       QMessageBox::Ok))
                  {
                    Visitor visitor(selection_item->selected_facets);
                    CGAL::Polygon_mesh_processing::triangulate_faces(selection_item->selected_facets,
                      pmesh,
                      CGAL::Polygon_mesh_processing::parameters::visitor(visitor));
                    break;
                  }
                  else
                  {
                    return;
                  }
                }
              }

              if (fpmap_valid)
                CGAL::Polygon_mesh_processing::isotropic_remeshing(selection_item->selected_facets
                  , target_length
                  , *selection_item->polyhedron()
                  , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
                  .protect_constraints(protect)
                  .edge_is_constrained_map(selection_item->constrained_edges_pmap())
                  .relax_constraints(smooth_features)
                  .number_of_relaxation_steps(nb_smooth)
                  .vertex_is_constrained_map(selection_item->constrained_vertices_pmap())
                  .face_patch_map(fpmap));
              else
                CGAL::Polygon_mesh_processing::isotropic_remeshing(selection_item->selected_facets
                  , target_length
                  , *selection_item->polyhedron()
                  , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
                  .protect_constraints(protect)
                  .edge_is_constrained_map(selection_item->constrained_edges_pmap())
                  .relax_constraints(smooth_features)
                  .number_of_relaxation_steps(nb_smooth)
                  .vertex_is_constrained_map(selection_item->constrained_vertices_pmap()));
            }
        }

        selection_item->polyhedron_item()->setColor(
              selection_item->polyhedron_item()->color());
        if(fpmap_valid)
        {
          selection_item->polyhedron_item()->setItemIsMulticolor(true);
          selection_item->polyhedron_item()->computeItemColorVectorAutomatically(true);
        }
        else
        {
          selection_item->polyhedron_item()->setItemIsMulticolor(false);
        }
        selection_item->setKeepSelectionValid(Scene_polyhedron_selection_item::Edge);
        selection_item->polyhedron_item()->invalidateOpenGLBuffers();
        Q_EMIT selection_item->polyhedron_item()->itemChanged();
        selection_item->invalidateOpenGLBuffers();
        selection_item->setKeepSelectionValid(Scene_polyhedron_selection_item::None);
      }
      else if (poly_item)
      {
        boost::property_map<FaceGraph, CGAL::edge_is_feature_t>::type eif
          = get(CGAL::edge_is_feature, pmesh);
        if (edges_only)
        {
          std::vector<edge_descriptor> edges_to_split;
          for(edge_descriptor e : edges(pmesh))
          {
            if( is_border(e, pmesh) || get(eif, e) )
              edges_to_split.push_back(e);
          }

          if (!edges_to_split.empty())
          {
            if (fpmap_valid)
              CGAL::Polygon_mesh_processing::split_long_edges(
                edges_to_split
                , target_length
                , pmesh
                , PMP::parameters::geom_traits(EPICK())
                . edge_is_constrained_map(eif)
                . face_patch_map(fpmap));
            else
              CGAL::Polygon_mesh_processing::split_long_edges(
                edges_to_split
                , target_length
                , pmesh
                , PMP::parameters::geom_traits(EPICK())
                . edge_is_constrained_map(eif));
          }
          else
            std::cout << "No border to be split" << std::endl;
        }
        else
        {
          // tricks to use the function detect_and_split_duplicates
          // that uses several poly items
          std::map<FaceGraph*, Edge_set > edges_to_protect_map;
          std::vector<Scene_facegraph_item*> poly_items(1, poly_item);
          Edge_set& edges_to_protect = edges_to_protect_map[poly_item->polyhedron()];
          if (preserve_duplicates)
          {
            detect_and_split_duplicates(poly_items, edges_to_protect_map, target_length);
          }
          Scene_polyhedron_selection_item::Is_constrained_map<Edge_set> ecm(&edges_to_protect);
          for(edge_descriptor e : edges(pmesh))
          {
            if (eif[e])
              edges_to_protect.insert(e);
          }

          if(protect &&
             !CGAL::Polygon_mesh_processing::internal::constraints_are_short_enough(
               pmesh,
               ecm,
               get(CGAL::vertex_point, pmesh),
               CGAL::Constant_property_map<face_descriptor, std::size_t>(1),
               4. / 3. * target_length))
          {
            QApplication::restoreOverrideCursor();
            QMessageBox::warning(mw, tr("Error"),
                                 tr("Isotropic remeshing : protect_constraints cannot be set to"
                                    " true with constraints larger than 4/3 * target_edge_length."
                                    " Aborting."));
            return;
          }

          if (!CGAL::is_triangle_mesh(pmesh))
          {
            QApplication::restoreOverrideCursor();
            if (QMessageBox::Ok ==
                QMessageBox::question(mw, tr("Error - Triangulate Faces?"),
                  tr("The input mesh is not a triangulated surface mesh.\n"
                     "Do you wish to triangulate faces first, or cancel remeshing ?"),
                  (QMessageBox::Ok | QMessageBox::Cancel), QMessageBox::Ok))
            {
              QApplication::setOverrideCursor(Qt::WaitCursor);
              CGAL::Polygon_mesh_processing::triangulate_faces(pmesh);
            }
            else
            {
              return;
            }
          }

          if (fpmap_valid)
            CGAL::Polygon_mesh_processing::isotropic_remeshing(
                 faces(*poly_item->polyhedron())
               , target_length
               , *poly_item->polyhedron()
               , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
               .protect_constraints(protect)
               .number_of_relaxation_steps(nb_smooth)
               .edge_is_constrained_map(ecm)
               .relax_constraints(smooth_features)
               .face_patch_map(fpmap));
          else
            CGAL::Polygon_mesh_processing::isotropic_remeshing(
                 faces(*poly_item->polyhedron())
               , target_length
               , *poly_item->polyhedron()
               , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
               .protect_constraints(protect)
               .number_of_relaxation_steps(nb_smooth)
               .edge_is_constrained_map(ecm)
               .relax_constraints(smooth_features));

          //recollect sharp edges
          for(edge_descriptor e : edges(pmesh))
            eif[e] = false;
          for(edge_descriptor e : edges_to_protect)
            eif[e] = true;
        }
        if (fpmap_valid)
        {
          poly_item->setItemIsMulticolor(true);
          poly_item->show_feature_edges(true);
        }
        else
          poly_item->setItemIsMulticolor(false);

        poly_item->invalidateOpenGLBuffers();

        Q_EMIT poly_item->itemChanged();
      }
      else{
        std::cout << "Can't remesh that type of thing" << std::endl;
      }
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
    }
    // default cursor
    QApplication::restoreOverrideCursor();
  }

  void isotropic_remeshing_of_several_polyhedra()
  {
    // Remeshing parameters
    bool edges_only = false, preserve_duplicates = false;
    double target_length = 0.;
    unsigned int nb_iter = 1;
    bool protect = false;
    bool smooth_features = true;

    std::vector<Scene_facegraph_item*> selection;
    for(int index : scene->selectionIndices())
    {
      Scene_facegraph_item* poly_item =
        qobject_cast<Scene_facegraph_item*>(scene->item(index));

      if (poly_item == NULL)
      {
        std::cout << scene->item(index)->name().data()
          << " is not a FaceGraph, remeshing skipped\n";
        continue;
      }
      else
      {
        selection.push_back(poly_item);

        if (target_length == 0.)//parameters have not been set yet
        {
        QDialog dialog(mw);
        Ui::Isotropic_remeshing_dialog ui = remeshing_dialog(&dialog, poly_item);
        ui.objectName->setText(QString::number(scene->selectionIndices().size())
          .append(QString(" items to be remeshed")));
        int i = dialog.exec();
        if (i == QDialog::Rejected)
        {
          std::cout << "Remeshing aborted" << std::endl;
          return;
        }

        edges_only = ui.splitEdgesOnly_checkbox->isChecked();
        preserve_duplicates = ui.preserveDuplicates_checkbox->isChecked();
        target_length = ui.edgeLength_dspinbox->value();
        nb_iter = ui.nbIterations_spinbox->value();
        protect = ui.protect_checkbox->isChecked();
        smooth_features = ui.smooth1D_checkbox->isChecked();
        }
      }
    }

    if(target_length == 0.)//parameters have not been set
    {                      // i.e. no item is a polyhedron
      std::cout << "Remeshing aborted" << std::endl;
      return;
    }


    //check non-triangulated surfaces
    for (Scene_facegraph_item* poly_item : selection)
    {
      if (!CGAL::is_triangle_mesh(*poly_item->polyhedron()))
      {
        if (QMessageBox::Ok == QMessageBox::question(mw,
              tr("Error - Triangulate Faces?"),
              tr("The input mesh ").append(poly_item->name())
               .append(tr(" is not a triangulated surface mesh.\n"
                "Do you wish to triangulate faces first, or cancel remeshing ?")),
              (QMessageBox::Ok | QMessageBox::Cancel),
              QMessageBox::Ok))
        {
          QApplication::setOverrideCursor(Qt::WaitCursor);
          CGAL::Polygon_mesh_processing::triangulate_faces(*poly_item->polyhedron());
          QApplication::restoreOverrideCursor();
        }
        else
        {
          QApplication::restoreOverrideCursor();
          return;
        }
      }
    }

    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);
    int total_time = 0;

    std::map<FaceGraph*,Edge_set > edges_to_protect;

    if(preserve_duplicates)
      detect_and_split_duplicates(selection, edges_to_protect, target_length);

#ifdef CGAL_LINKED_WITH_TBB
    QElapsedTimer time;
    time.start();

      tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, selection.size()),
        Remesh_polyhedron_item_for_parallel_for<Remesh_polyhedron_item>(
                                                                        selection, edges_to_protect, edges_only, target_length, nb_iter, protect, smooth_features));

    total_time = time.elapsed();

#else

    Remesh_polyhedron_item remesher(edges_only,
      target_length, nb_iter, protect, smooth_features);
    for(Scene_facegraph_item* poly_item : selection)
    {
      QElapsedTimer time;
      time.start();

      remesher(poly_item, edges_to_protect[poly_item->polyhedron()]);

      total_time += time.elapsed();
      std::cout << "Remeshing of " << poly_item->name().data()
                << " done in " << time.elapsed() << " ms" << std::endl;
    }
#endif
    std::cout << "Remeshing of all selected items done in "
      << total_time << " ms" << std::endl;

    for(Scene_facegraph_item* poly_item : selection)
    {
      //destroys the patch_id_map for the Surface_mesh_item to avoid assertions.
      poly_item->resetColors();
      poly_item->invalidateOpenGLBuffers();
      Q_EMIT poly_item->itemChanged();
    }

    // default cursor
    QApplication::restoreOverrideCursor();
  }

private:
  Scene_interface *scene;
  QMainWindow* mw;
  struct Remesh_polyhedron_item
  {
    typedef boost::graph_traits<FaceGraph>::edge_descriptor     edge_descriptor;
    typedef boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
    typedef boost::graph_traits<FaceGraph>::face_descriptor     face_descriptor;

    bool edges_only_;
    double target_length_;
    unsigned int nb_iter_;
    bool protect_;
    bool smooth_features_;

  protected:
    void remesh(Scene_facegraph_item* poly_item,
                Edge_set& edges_to_protect) const
    {
      //fill face_index property map

      if (edges_only_)
      {
        std::vector<halfedge_descriptor> border;
        CGAL::Polygon_mesh_processing::border_halfedges(
          faces(*poly_item->polyhedron())
          , *poly_item->polyhedron()
          , std::back_inserter(border));
        std::vector<edge_descriptor> border_edges;
        for(halfedge_descriptor h : border)
          border_edges.push_back(edge(h, *poly_item->polyhedron()));

        CGAL::Polygon_mesh_processing::split_long_edges(
            border_edges
          , target_length_
          , *poly_item->polyhedron());
      }
      else
      {
        std::cout << "Isotropic remeshing of "
          << poly_item->name().toStdString() << " started..." << std::endl;
        Scene_polyhedron_selection_item::Is_constrained_map<Edge_set> ecm(&edges_to_protect);
        CGAL::Polygon_mesh_processing::isotropic_remeshing(
            faces(*poly_item->polyhedron())
          , target_length_
          , *poly_item->polyhedron()
          , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter_)
          .protect_constraints(protect_)
          .edge_is_constrained_map(ecm)
          .face_patch_map(get(CGAL::face_patch_id_t<int>(), *poly_item->polyhedron()))
          .relax_constraints(smooth_features_));
        std::cout << "Isotropic remeshing of "
          << poly_item->name().toStdString() << " done." << std::endl;
      }
    }

  public:
    Remesh_polyhedron_item(
      const bool edges_only,
      const double target_length,
      const unsigned int nb_iter,
      const bool protect,
      const bool smooth_features)
      : edges_only_(edges_only)
      , target_length_(target_length)
      , nb_iter_(nb_iter)
      , protect_(protect)
      , smooth_features_(smooth_features)
    {}

    Remesh_polyhedron_item(const Remesh_polyhedron_item& remesh)
      : edges_only_(remesh.edges_only_)
      , target_length_(remesh.target_length_)
      , nb_iter_(remesh.nb_iter_)
      , protect_(remesh.protect_)
      , smooth_features_(remesh.smooth_features_)
    {}

    void operator()(Scene_facegraph_item* poly_item,
                    Edge_set& edges_to_protect) const
    {
      remesh(poly_item, edges_to_protect);
    }
  };

#ifdef CGAL_LINKED_WITH_TBB
  template<typename RemeshFunctor>
  struct Remesh_polyhedron_item_for_parallel_for
    : RemeshFunctor
  {
    const std::vector<Scene_facegraph_item*>& selection_;
    std::map<FaceGraph*,Edge_set >& edges_to_protect_;

  public:
    // Constructor
    Remesh_polyhedron_item_for_parallel_for(
      const std::vector<Scene_facegraph_item*>& selection,
      std::map<FaceGraph*,Edge_set >& edges_to_protect,
      const bool edges_only,
      const double target_length,
      const unsigned int nb_iter,
      const bool protect,
      const bool smooth_features)
      : RemeshFunctor(edges_only, target_length, nb_iter, protect, smooth_features)
      , selection_(selection), edges_to_protect_(edges_to_protect)
    {}

    // Constructor
    Remesh_polyhedron_item_for_parallel_for(
      const Remesh_polyhedron_item_for_parallel_for &remesh)
      : RemeshFunctor(remesh)
      , selection_(remesh.selection_)
      , edges_to_protect_(remesh.edges_to_protect_)
    {}

    // operator()
    void operator()(const tbb::blocked_range<size_t>& r) const
    {
      for (size_t i = r.begin(); i != r.end(); ++i)
        RemeshFunctor::remesh(selection_[i], edges_to_protect_[selection_[i]->polyhedron()]);
    }
  };
#endif

  Ui::Isotropic_remeshing_dialog
  remeshing_dialog(QDialog* dialog,
                   Scene_facegraph_item* poly_item,
                   Scene_polyhedron_selection_item* selection_item = NULL)
  {
    Ui::Isotropic_remeshing_dialog ui;
    ui.setupUi(dialog);
    connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

    //connect checkbox to spinbox
    connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
            ui.nbIterations_spinbox, SLOT(setDisabled(bool)));
    connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
            ui.protect_checkbox, SLOT(setDisabled(bool)));
    connect(ui.protect_checkbox, SIGNAL(toggled(bool)),
            ui.smooth1D_checkbox, SLOT(setDisabled(bool)));
    connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
            ui.smooth1D_checkbox, SLOT(setDisabled(bool)));
    connect(ui.preserveDuplicates_checkbox, SIGNAL(toggled(bool)),
            ui.protect_checkbox, SLOT(setChecked(bool)));
    connect(ui.preserveDuplicates_checkbox, SIGNAL(toggled(bool)),
            ui.protect_checkbox, SLOT(setDisabled(bool)));

    //Set default parameters
    Scene_interface::Bbox bbox = poly_item != NULL ? poly_item->bbox()
      : (selection_item != NULL ? selection_item->bbox()
        : scene->bbox());
    ui.objectName->setText(poly_item != NULL ? poly_item->name()
      : (selection_item != NULL ? selection_item->name()
        : QString("Remeshing parameters")));

    ui.objectNameSize->setText(
      tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
      .arg(bbox.xmax()-bbox.xmin(), 0, 'g', 3)
      .arg(bbox.ymax()-bbox.ymin(), 0, 'g', 3)
      .arg(bbox.zmax()-bbox.zmin(), 0, 'g', 3));

    double diago_length = CGAL::sqrt((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin())
                                   + (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin())
                                   + (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()));


    ui.edgeLength_dspinbox->setValue(0.05 * diago_length);

    std::ostringstream oss;
    oss << "Diagonal length of the Bbox of the selection to remesh is ";
    oss << diago_length << "." << std::endl;
    oss << "Default is 5% of it" << std::endl;
    ui.edgeLength_dspinbox->setToolTip(QString::fromStdString(oss.str()));

    ui.nbIterations_spinbox->setSingleStep(1);
    ui.nbIterations_spinbox->setRange(1/*min*/, 1000/*max*/);
    ui.nbIterations_spinbox->setValue(1);

    ui.protect_checkbox->setChecked(false);
    ui.smooth1D_checkbox->setChecked(true);

    if (NULL != selection_item)
    {
      //do not preserve duplicates in selection mode
      ui.preserveDuplicates_checkbox->setDisabled(true);
      ui.preserveDuplicates_checkbox->setChecked(false);
    }

    return ui;
  }


private:
  QAction* actionIsotropicRemeshing_;

}; // end Polyhedron_demo_isotropic_remeshing_plugin

#include "Isotropic_remeshing_plugin.moc"
