#include "Surface_mesh_item_classification.h"
#include "Color_ramp.h"

#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Three/Viewer_interface.h>

#include <QLineEdit>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>

Surface_mesh_item_classification::Surface_mesh_item_classification(Scene_surface_mesh_item* mesh)
  : m_mesh (mesh),
    m_selection (NULL),
    m_generator (NULL)
{
  m_index_color = 1;

  backup_existing_colors_and_add_new();
  m_training = m_mesh->polyhedron()->add_property_map<face_descriptor, std::size_t>("f:training", std::size_t(-1)).first;
  m_classif = m_mesh->polyhedron()->add_property_map<face_descriptor, std::size_t>("f:label", std::size_t(-1)).first;

  m_labels.add("ground");
  m_labels.add("vegetation");
  m_labels.add("roof");
  m_labels.add("facade");

  m_sowf = new Sum_of_weighted_features (m_labels, m_features);
  m_ethz = NULL;
#ifdef CGAL_LINKED_WITH_OPENCV
  m_random_forest = NULL;
#endif
#ifdef CGAL_LINKED_WITH_TENSORFLOW
  m_neural_network = NULL;
#endif
}


Surface_mesh_item_classification::~Surface_mesh_item_classification()
{
  if (m_sowf != NULL)
    delete m_sowf;
  if (m_ethz != NULL)
    delete m_ethz;
#ifdef CGAL_LINKED_WITH_OPENCV
  if (m_random_forest != NULL)
    delete m_random_forest;
#endif
#ifdef CGAL_LINKED_WITH_TENSORFLOW
  if (m_neural_network != NULL)
    delete m_neural_network;
#endif
  if (m_generator != NULL)
    delete m_generator;
}

void Surface_mesh_item_classification::backup_existing_colors_and_add_new()
{
  bool has_colors = false;
  boost::tie (m_color, has_colors) = m_mesh->polyhedron()->property_map<face_descriptor, CGAL::Color>("f:color");
  if (has_colors)
  {
    m_real_color
      = m_mesh->polyhedron()->add_property_map<face_descriptor, CGAL::Color>("f:real_color").first;
    for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
    {
      m_real_color[fd] = m_color[fd];
      m_color[fd] = CGAL::Color(128, 128, 128);
    }
  }
  else
    m_color =
      m_mesh->polyhedron()->add_property_map<face_descriptor, CGAL::Color>("f:color", CGAL::Color(128,128,128)).first;
}

void Surface_mesh_item_classification::change_color (int index, float* vmin, float* vmax)
{
  m_index_color = index;
  int index_color = index;
  if (index == 0 && m_real_color == Mesh::Property_map<face_descriptor, CGAL::Color>())
    index_color = -1;

  static Color_ramp ramp;
  ramp.build_rainbow();

  if (index_color == -1) // item color
  {
    for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
      m_color[fd] = CGAL::Color(128,128,128);
  }
  else if (index_color == 0) // real colors
  {
    for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
      m_color[fd] = m_real_color[fd];
  }
  else if (index_color == 1) // classif
  {
    for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
    {
      QColor color (128, 128, 128);
      std::size_t c = m_classif[fd];

      if (c != std::size_t(-1))
        color = label_qcolor (m_labels[c]);

      m_color[fd] = CGAL::Color(color.red(), color.green(), color.blue());
    }
  }
  else if (index_color == 2) // training
  {
    for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
    {
      QColor color (128, 128, 128);
      std::size_t c = m_training[fd];
      std::size_t c2 = m_classif[fd];

      if (c != std::size_t(-1))
        color = label_qcolor(m_labels[c]);

      float div = 1;
      if (c != c2)
        div = 2;
      m_color[fd] = CGAL::Color(color.red() / div,
                                color.green() / div,
                                color.blue() / div);
    }
  }
  else
  {
    std::size_t corrected_index = index_color - 3;
    if (corrected_index < m_labels.size()) // Display label probabilities
    {
      if (m_label_probabilities.size() <= corrected_index ||
          m_label_probabilities[corrected_index].size() != num_faces(*(m_mesh->polyhedron())))
      {
        for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
        {
          m_color[fd] = CGAL::Color((unsigned char)(128),
                                    (unsigned char)(128),
                                    (unsigned char)(128));
        }
      }
      else
      {
        for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
        {
          float v = (std::max) (0.f, (std::min)(1.f, m_label_probabilities[corrected_index][fd]));
          m_color[fd] = CGAL::Color((unsigned char)(ramp.r(v) * 255),
                                    (unsigned char)(ramp.g(v) * 255),
                                    (unsigned char)(ramp.b(v) * 255));
        }
      }
    }
    else
    {
      corrected_index -= m_labels.size();
      if (corrected_index >= m_features.size())
      {
        std::cerr << "Error: trying to access feature " << corrected_index << " out of " << m_features.size() << std::endl;
        return;
      }

      Feature_handle feature = m_features[corrected_index];

      float min = (std::numeric_limits<float>::max)();
      float max = -(std::numeric_limits<float>::max)();

      if (vmin != NULL && vmax != NULL
          && *vmin != std::numeric_limits<float>::infinity()
          && *vmax != std::numeric_limits<float>::infinity())
      {
        min = *vmin;
        max = *vmax;
      }
      else
      {
        for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
        {
          if (feature->value(fd) > max)
            max = feature->value(fd);
          if (feature->value(fd) < min)
            min = feature->value(fd);
        }
      }

      for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
      {
        float v = (feature->value(fd) - min) / (max - min);
        if (v < 0.f) v = 0.f;
        if (v > 1.f) v = 1.f;

        m_color[fd] = CGAL::Color((unsigned char)(ramp.r(v) * 255),
                                  (unsigned char)(ramp.g(v) * 255),
                                  (unsigned char)(ramp.b(v) * 255));
      }

      if (vmin != NULL && vmax != NULL)
      {
        *vmin = min;
        *vmax = max;
      }
    }
  }
}

void Surface_mesh_item_classification::compute_features (std::size_t nb_scales, float voxel_size)
{
  std::cerr << "Computing features with " << nb_scales << " scale(s) and ";
  if (voxel_size == -1)
    std::cerr << "automatic voxel size" << std::endl;
  else
    std::cerr << "voxel size = " << voxel_size << std::endl;

  m_features.clear();

  if (m_generator != NULL)
    delete m_generator;

  Face_center_map fc_map (m_mesh->polyhedron());

  m_generator = new Generator (*(m_mesh->polyhedron()), fc_map, nb_scales, voxel_size);

#ifdef CGAL_LINKED_WITH_TBB
  m_features.begin_parallel_additions();
#endif

  m_generator->generate_point_based_features(m_features);
  m_generator->generate_face_based_features(m_features);

#ifdef CGAL_LINKED_WITH_TBB
  m_features.end_parallel_additions();
#endif

  delete m_sowf;
  m_sowf = new Sum_of_weighted_features (m_labels, m_features);
  if (m_ethz != NULL)
  {
    delete m_ethz;
    m_ethz = NULL;
  }
#ifdef CGAL_LINKED_WITH_OPENCV
  if (m_random_forest != NULL)
  {
    delete m_random_forest;
    m_random_forest = NULL;
  }
#endif
#ifdef CGAL_LINKED_WITH_TENSORFLOW
  if (m_neural_network != NULL)
  {
    delete m_neural_network;
    m_neural_network = NULL;
  }
#endif
  std::cerr << "Features = " << m_features.size() << std::endl;
}

void Surface_mesh_item_classification::train (int classifier, const QMultipleInputDialog& dialog)
{
  if (m_features.size() == 0)
  {
    std::cerr << "Error: features not computed" << std::endl;
    return;
  }

  m_label_probabilities.clear();
  m_label_probabilities.resize (m_labels.size());
  for (std::size_t i = 0; i < m_label_probabilities.size(); ++ i)
    m_label_probabilities[i].resize (num_faces(*(m_mesh->polyhedron())));

  std::vector<std::size_t> training (num_faces(*(m_mesh->polyhedron())), std::size_t(-1));
  std::vector<std::size_t> indices (num_faces(*(m_mesh->polyhedron())), std::size_t(-1));

  std::vector<std::size_t> nb_label (m_labels.size(), 0);
  std::size_t nb_total = 0;

  for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
  {
    training[fd] = m_training[fd];
    if (training[fd] != std::size_t(-1))
    {
      nb_label[training[fd]] ++;
      ++ nb_total;
    }
  }

  std::cerr << nb_total << " face(s) used for training ("
            << 100. * (nb_total / double(m_mesh->polyhedron()->faces().size())) << "% of the total):" << std::endl;
  for (std::size_t i = 0; i < m_labels.size(); ++ i)
    std::cerr << " * " << m_labels[i]->name() << ": " << nb_label[i] << " face(s)" << std::endl;

  if (classifier == CGAL_CLASSIFICATION_SOWF_NUMBER)
  {
    m_sowf->train<Concurrency_tag>(training, dialog.get<QSpinBox>("trials")->value());
    CGAL::Classification::classify<Concurrency_tag> (m_mesh->polyhedron()->faces(),
                                                     m_labels, *m_sowf,
                                                     indices, m_label_probabilities);
  }
  else if (classifier == CGAL_CLASSIFICATION_ETHZ_NUMBER)
  {
    if (m_ethz != NULL)
      delete m_ethz;
    m_ethz = new ETHZ_random_forest (m_labels, m_features);
    m_ethz->train<Concurrency_tag>(training, true,
                                   dialog.get<QSpinBox>("num_trees")->value(),
                                   dialog.get<QSpinBox>("max_depth")->value());
    CGAL::Classification::classify<Concurrency_tag> (m_mesh->polyhedron()->faces(),
                                                     m_labels, *m_ethz,
                                                     indices, m_label_probabilities);
  }
  else if (classifier == CGAL_CLASSIFICATION_OPENCV_NUMBER)
  {
#ifdef CGAL_LINKED_WITH_OPENCV
    if (m_random_forest != NULL)
      delete m_random_forest;
    m_random_forest = new Random_forest (m_labels, m_features,
                                         dialog.get<QSpinBox>("max_depth")->value(), 5, 15,
                                         dialog.get<QSpinBox>("num_trees")->value());
    m_random_forest->train (training);

    CGAL::Classification::classify<Concurrency_tag> (m_mesh->polyhedron()->faces(),
                                                     m_labels, *m_random_forest,
                                                     indices, m_label_probabilities);
#endif
  }
  else if (classifier == CGAL_CLASSIFICATION_TENSORFLOW_NUMBER)
  {
#ifdef CGAL_LINKED_WITH_TENSORFLOW
    if (m_neural_network != NULL)
    {
      if (m_neural_network->initialized())
      {
        if (dialog.get<QCheckBox>("restart")->isChecked())
        {
          delete m_neural_network;
          m_neural_network = new Neural_network (m_labels, m_features);
        }
      }
      else
      {
        delete m_neural_network;
        m_neural_network = new Neural_network (m_labels, m_features);
      }
    }
    else
      m_neural_network = new Neural_network (m_labels, m_features);

    std::vector<std::size_t> hidden_layers;

    std::string hl_input = dialog.get<QLineEdit>("hidden_layers")->text().toStdString();
    if (hl_input != "")
    {
      std::istringstream iss(hl_input);
      int s;
      while (iss >> s)
        hidden_layers.push_back (std::size_t(s));
    }

    m_neural_network->train (training,
                             dialog.get<QCheckBox>("restart")->isChecked(),
                             dialog.get<QSpinBox>("trials")->value(),
                             dialog.get<DoubleEdit>("learning_rate")->value(),
                             dialog.get<QSpinBox>("batch_size")->value(),
                             hidden_layers);

    CGAL::Classification::classify<Concurrency_tag> (m_mesh->polyhedron()->faces(),
                                                     m_labels, *m_neural_network,
                                                     indices, m_label_probabilities);
#endif
  }

  for(face_descriptor fd : faces(*(m_mesh->polyhedron())))
    m_classif[fd] = indices[fd];

  if (m_index_color == 1 || m_index_color == 2)
     change_color (m_index_color);
}

bool Surface_mesh_item_classification::run (int method, int classifier,
                                            std::size_t subdivisions, double smoothing)
{
  if (m_features.size() == 0)
  {
    std::cerr << "Error: features not computed" << std::endl;
    return false;
  }

  if (classifier == CGAL_CLASSIFICATION_SOWF_NUMBER)
    run (method, *m_sowf, subdivisions, smoothing);
  else if (classifier == CGAL_CLASSIFICATION_ETHZ_NUMBER)
  {
    if (m_ethz == NULL)
    {
      std::cerr << "Error: ETHZ Random Forest must be trained or have a configuration loaded first" << std::endl;
      return false;
    }
    run (method, *m_ethz, subdivisions, smoothing);
  }
  else if (classifier == CGAL_CLASSIFICATION_OPENCV_NUMBER)
  {
#ifdef CGAL_LINKED_WITH_OPENCV
    if (m_random_forest == NULL)
    {
      std::cerr << "Error: OpenCV Random Forest must be trained or have a configuration loaded first" << std::endl;
      return false;
    }
    run (method, *m_random_forest, subdivisions, smoothing);
#endif
  }
  else if (classifier == CGAL_CLASSIFICATION_TENSORFLOW_NUMBER)
  {
#ifdef CGAL_LINKED_WITH_TENSORFLOW
    if (m_neural_network == NULL)
    {
      std::cerr << "Error: TensorFlow Neural Network must be trained or have a configuration loaded first" << std::endl;
      return false;
    }
    run (method, *m_neural_network, subdivisions, smoothing);
#endif
  }

  return true;
}
