#ifndef SURFACE_MESH_ITEM_CLASSIFICATION_H
#define SURFACE_MESH_ITEM_CLASSIFICATION_H

//#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#define CGAL_CLASSIFICATION_VERBOSE

#include <CGAL/Three/Scene_item.h>

#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Item_classification_base.h"
#include "Kernel_type.h"

#include <CGAL/Classification.h>


#include <iostream>

typedef CGAL::Parallel_if_available_tag Concurrency_tag;

class Surface_mesh_item_classification : public Item_classification_base
{
public:

  typedef SMesh Mesh;
  typedef Kernel::Point_3 Point;
  typedef Scene_polyhedron_selection_item::Selection_set_facet Selection;
  typedef Mesh::Face_range Face_range;
  typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typedef CGAL::Identity_property_map<face_descriptor> Face_map;

  typedef CGAL::Classification::Mesh_neighborhood<Mesh> Neighborhood;

  typedef CGAL::Classification::Face_descriptor_to_center_of_mass_map<Mesh> Face_center_map;
  typedef CGAL::Classification::Face_descriptor_to_face_descriptor_with_bbox_map<Mesh> Face_descriptor_with_bbox_map;
  typedef CGAL::Classification::Mesh_feature_generator<Kernel, Mesh, Face_center_map> Generator;

public:

  Surface_mesh_item_classification(Scene_surface_mesh_item* mesh);
  ~Surface_mesh_item_classification();

  void backup_existing_colors_and_add_new();

  CGAL::Three::Scene_item* item() { return m_mesh; }
  void erase_item() { m_mesh = NULL; }

  void compute_features (std::size_t nb_scales, float voxel_size);

  void add_selection_to_training_set (std::size_t label)
  {
    if (m_selection == NULL)
      return;

    for (Selection::iterator it = m_selection->selected_facets.begin();
         it != m_selection->selected_facets.end(); ++ it)
    {
      m_classif[*it] = label;
      m_training[*it] = label;
    }
    m_selection->clear_all();

    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_set (std::size_t label)
  {
    for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
      if (m_training[fd] == label)
        m_training[fd] = std::size_t(-1);
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }

  void reset_training_set_of_selection ()
  {
    for (Selection::iterator it = m_selection->selected_facets.begin();
         it != m_selection->selected_facets.end(); ++ it)
    {
      m_training[*it] = -1;
      m_classif[*it] = -1;
    }
    m_selection->clear_all();

    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }

  void reset_training_sets()
  {
    for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
      m_training[fd] = std::size_t(-1);
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void select_random_region()
  {
    // TODO
    std::cerr << "Warning: operation not yet available for meshes." << std::endl;
  }
  void validate_selection ()
  {
    if (m_selection == NULL)
      return;

    for (Selection::iterator it = m_selection->selected_facets.begin();
         it != m_selection->selected_facets.end(); ++ it)
      m_training[*it] = m_classif[*it];
    m_selection->clear_all();

    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void train(int classifier, const QMultipleInputDialog& dialog);
  bool run (int method, int classifier, std::size_t subdivisions, double smoothing);

  void update_color() { change_color (m_index_color); }
  void change_color (int index, float* vmin = NULL, float* vmax = NULL);
  CGAL::Three::Scene_item* generate_one_item (const char* /* name */,
                                              int /* label */) const
  {
    // TODO
    std::cerr << "Warning: operation not yet available for meshes." << std::endl;
    return NULL;
  }
  void generate_one_item_per_label(std::vector<CGAL::Three::Scene_item*>&,
                                   const char*) const
  {
    // TODO
    std::cerr << "Warning: operation not yet available for meshes." << std::endl;
  }

  void set_selection_item (Scene_polyhedron_selection_item* selection)
  {
    m_selection = selection;
  }

protected:

  template <typename Classifier>
  bool run (int method, const Classifier& classifier,
            std::size_t subdivisions, double smoothing)
  {
    std::vector<std::size_t> indices(m_mesh->polyhedron()->faces().size(), -1);
    Face_center_map fc_map (m_mesh->polyhedron());

    if (method == 0)
      CGAL::Classification::classify<Concurrency_tag> (m_mesh->polyhedron()->faces(),
                                                       m_labels, classifier,
                                                       indices);
    else if (method == 1)
      CGAL::Classification::classify_with_local_smoothing<Concurrency_tag>
        (m_mesh->polyhedron()->faces(), Face_map(), m_labels, classifier,
         m_generator->neighborhood().n_ring_neighbor_query(2),
         indices);
    else if (method == 2)
    {
      // Fix: need bbox of face

      CGAL::Classification::classify_with_graphcut<Concurrency_tag>
        (m_mesh->polyhedron()->faces(), Face_descriptor_with_bbox_map(m_mesh->polyhedron()),
         m_labels, classifier,
         m_generator->neighborhood().n_ring_neighbor_query(1),
         smoothing, subdivisions, indices);
    }

    std::vector<std::size_t> ground_truth(num_faces(*(m_mesh->polyhedron())), std::size_t(-1));
    for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
    {
      m_classif[fd] = indices[fd];
      ground_truth[fd] = m_training[fd];
    }

    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);

    std::cerr << "Precision, recall, F1 scores and IoU:" << std::endl;

    CGAL::Classification::Evaluation eval (m_labels, ground_truth, indices);

    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      {
        std::cerr << " * " << m_labels[i]->name() << ": "
                  << eval.precision(m_labels[i]) << " ; "
                  << eval.recall(m_labels[i]) << " ; "
                  << eval.f1_score(m_labels[i]) << " ; "
                  << eval.intersection_over_union(m_labels[i]) << std::endl;
      }

    std::cerr << "Accuracy = " << eval.accuracy() << std::endl
              << "Mean F1 score = " << eval.mean_f1_score() << std::endl
              << "Mean IoU = " << eval.mean_intersection_over_union() << std::endl;


    return true;
  }

  Scene_surface_mesh_item* m_mesh;
  Scene_polyhedron_selection_item* m_selection;
  Mesh::Property_map<face_descriptor, std::size_t> m_training;
  Mesh::Property_map<face_descriptor, std::size_t> m_classif;
  Mesh::Property_map<face_descriptor, CGAL::Color> m_color;
  Mesh::Property_map<face_descriptor, CGAL::Color> m_real_color;

  std::vector<std::vector<float> > m_label_probabilities;

  Generator* m_generator;
  int m_index_color;
};




#endif // SURFACE_MESH_ITEM_CLASSIFICATION_H
