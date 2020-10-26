// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_LABEL_SET_H
#define CGAL_CLASSIFICATION_LABEL_SET_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Label.h>

#include <CGAL/Random.h>

#include <vector>

namespace CGAL {

namespace Classification {

/*!
\ingroup PkgClassificationLabel

\brief Set of `Label` used as input by classification
algorithms.

*/
class Label_set
{
  using Base = std::vector<Label_handle>;

  CGAL::Random m_random;
  Base m_labels;

public:

#ifdef DOXYGEN_RUNNING
  using const_iterator = unspecified_type; ///< A random access constant iterator with value type `Label_handle`.
  using iterator = unspecified_type; ///< A random access iterator with value type `Label_handle`.
#else
  using const_iterator = std::vector<Label_handle>::const_iterator;
  using iterator = std::vector<Label_handle>::iterator;
#endif

  /// \name Constructors
  /// @{

  Label_set() { }

  /*!
    \brief constructs a label set from a set of label names.
  */
  Label_set(std::initializer_list<const char*> labels)
  {
    for (const char* l : labels)
      add (l);
  }

  /// @}

  /// \name Modifications
  /// @{

  /*!
    \brief adds a label.

    \note Names, standard indices and colors are not used for
    identification: two labels in the same set can have the same name,
    standard index or color, but not the same handle. Each call to
    `add()` generates a new distinct label.

    \param name name of the label.

    \param color used to represent the label.

    \param standard_index standard index of the classification label
    (i.e. index in the ASPRS standard).

    \return a handle to the newly added label.
  */
  Label_handle add (const char* name,
                    CGAL::Color color,
                    std::size_t standard_index = -1)
  {
    Label_handle out = std::make_shared<Classification::Label>
      (name, m_labels.size(), standard_index, color);
    m_labels.push_back (out);
    return out;
  }


  /*!
    \brief adds a label with default standard index and color.

    This functions tries to map label names to standard ASPRS labels
    and automatically picks the `standard_index` and `color` of the
    label:

    - `"unassigned"` is given standard index 2 and color `(0, 0, 0)`
    - `"ground"` is given standard index 2 and color `(186, 189, 182)`
    - `"low_vegetation"` is given standard index 3 and color `(78, 154, 6)`
    - `"medium_vegetation"` is given standard index 4 and color `(138, 226, 52)`
    - `"high_vegetation"` is given standard index 5 and color `(204, 255, 201)`
    - `"building"` is given standard index 6 and color `(245, 121, 0)`
    - `"noise"` is given standard index 7 and color `(128, 0, 0)`
    - `"reserved"` is given standard index 8 and color `(233, 185, 110)`
    - `"water"` is given standard index 9 and color `(114, 159, 207)`
    - `"rail"` is given standard index 10 and color `(136, 46, 25)`
    - `"road_surface"` is given standard index 11 and color `(56, 56, 56)`
    - `"reserved_2"` is given standard index 12 and color `(193, 138, 51)`
    - `"wire_guard"` is given standard index 13 and color `(37, 61, 136)`
    - `"wire_conductor"` is given standard index 14 and color `(173, 127, 168)`
    - `"transmission_tower"` is given standard index 15 and color `(136, 138, 133)`
    - `"wire_connect"` is given standard index 16 and color `(145, 64, 236)`
    - `"bridge_deck"` is given standard index 17 and color `(213, 93, 93)`
    - `"high_noise"` is given standard index 18 and color `(255, 0, 0)`

    If the name is not found, the label is given standard index
    `std::size_t(-1)` and a random color.

    \note Names are not used for identification: two labels in the
    same set can have the same name but not the same handle. Each call
    to `add()` generates a new distinct label.

    \param name name of the label.

    \return a handle to the newly added label.
  */
  Label_handle add (const char* name)
  {
    static std::unordered_map<std::string, std::pair<std::size_t, CGAL::Color> > init_map;
    if (init_map.empty())
    {
      init_map.insert (std::make_pair ("unassigned",
                                       std::make_pair (2, CGAL::Color (0, 0, 0))));
      init_map.insert (std::make_pair ("ground",
                                       std::make_pair (2, CGAL::Color (186, 189, 182))));
      init_map.insert (std::make_pair ("low_vegetation",
                                       std::make_pair (3, CGAL::Color (78, 154, 6))));
      init_map.insert (std::make_pair ("medium_vegetation",
                                       std::make_pair (4, CGAL::Color (138, 226, 52))));
      init_map.insert (std::make_pair ("high_vegetation",
                                       std::make_pair (5, CGAL::Color (204, 255, 201))));
      init_map.insert (std::make_pair ("building",
                                       std::make_pair (6, CGAL::Color (245, 121, 0))));
      init_map.insert (std::make_pair ("noise",
                                       std::make_pair (7, CGAL::Color (128, 0, 0))));
      init_map.insert (std::make_pair ("reserved",
                                       std::make_pair (8, CGAL::Color (233, 185, 110))));
      init_map.insert (std::make_pair ("water",
                                       std::make_pair (9, CGAL::Color (114, 159, 207))));
      init_map.insert (std::make_pair ("rail",
                                       std::make_pair (10, CGAL::Color (136, 46, 25))));
      init_map.insert (std::make_pair ("road_surface",
                                       std::make_pair (11, CGAL::Color (56, 56, 56))));
      init_map.insert (std::make_pair ("reserved_2",
                                       std::make_pair (12, CGAL::Color (193, 138, 51))));
      init_map.insert (std::make_pair ("wire_guard",
                                       std::make_pair (13, CGAL::Color (37, 61, 136))));
      init_map.insert (std::make_pair ("wire_conductor",
                                       std::make_pair (14, CGAL::Color (173, 127, 168))));
      init_map.insert (std::make_pair ("wire_conduct",
                                       std::make_pair (14, CGAL::Color (173, 127, 168))));
      init_map.insert (std::make_pair ("transmission_tower",
                                       std::make_pair (15, CGAL::Color (136, 138, 133))));
      init_map.insert (std::make_pair ("trans_tower",
                                       std::make_pair (15, CGAL::Color (136, 138, 133))));
      init_map.insert (std::make_pair ("wire_connect",
                                       std::make_pair (16, CGAL::Color (145, 64, 236))));
      init_map.insert (std::make_pair ("bridge_deck",
                                       std::make_pair (17, CGAL::Color (213, 93, 93))));
      init_map.insert (std::make_pair ("high_noise",
                                       std::make_pair (18, CGAL::Color (255, 0, 0))));

      // Undocumented additions
      init_map.insert (std::make_pair ("low_veget",
                                       std::make_pair (3, CGAL::Color (78, 154, 6))));
      init_map.insert (std::make_pair ("medium_veget",
                                       std::make_pair (4, CGAL::Color (138, 226, 52))));
      init_map.insert (std::make_pair ("vegetation",
                                       std::make_pair (4, CGAL::Color (138, 226, 52))));
      init_map.insert (std::make_pair ("high_veget",
                                       std::make_pair (5, CGAL::Color (204, 255, 201))));
      init_map.insert (std::make_pair ("roof",
                                       std::make_pair (6, CGAL::Color (245, 121, 0))));
      init_map.insert (std::make_pair ("facade",
                                       std::make_pair (-1, CGAL::Color (77, 131, 186))));
    }

    std::string sname (name);
    auto found = init_map.find (sname);
    if (found == init_map.end())
      return add (name,
                  CGAL::Color ((unsigned char)(m_random.get_int(64, 192)),
                               (unsigned char)(m_random.get_int(64, 192)),
                               (unsigned char)(m_random.get_int(64, 192))));

    // else
    return add (name, found->second.second, found->second.first);
  }
  /// \endcond

  /*!
    \brief removes a label.

    \param label the handle to the label that must be removed.

    \return `true` if the label was correctly removed,
    `false` if its handle was not found.
  */
  bool remove (Label_handle label)
  {
    if (label->index() >= m_labels.size()
        || m_labels[label->index()] != label)
      return false;

    for (std::size_t i = label->index() + 1; i < m_labels.size(); ++ i)
      m_labels[i]->m_index --;
    m_labels.erase (m_labels.begin() + label->index());

    return true;
  }

  /*!
    \brief removes all labels.
  */
  void clear ()
  {
    m_labels.clear();
  }

  /// @}


  /// \name Access
  /// @{

  const_iterator begin() const { return m_labels.begin(); }
  iterator begin() { return m_labels.begin(); }
  const_iterator end() const { return m_labels.end(); }
  iterator end() { return m_labels.end(); }

  /*!
    \brief returns how many labels are defined.
  */
  std::size_t size () const
  {
    return m_labels.size();
  }

  /*!
    \brief returns the \f$i^{th}\f$ label.
  */
  Label_handle operator[] (std::size_t i) const
  {
    return m_labels[i];
  }

  /// @}

  /// \name Validity
  /// @{

  /*!
    \brief checks the validity of the ground truth with respect to the
    label set.

    \param ground_truth range of label indices. This function checks
    that all these indices are either -1 (for unclassified) or a valid
    index of one of the labels. If at least one of the indices is out
    of range, this function returns `false`, otherwise it returns
    `true`.

    \param verbose if set to `true`, the number of inliers of each
    label, the number of unclassified items and the potential number
    of out-of-range items are displayed. Otherwise, this function does
    not display anything.
  */
  template <typename LabelIndexRange>
  bool is_valid_ground_truth (const LabelIndexRange& ground_truth,
                              bool verbose = false) const
  {
    std::vector<std::size_t> nb_inliers (m_labels.size() + 2, 0);
    std::size_t total = 0;

    for (const auto& gt : ground_truth)
    {
      int g = int(gt);
      if (g == -1)
        ++ nb_inliers[m_labels.size()];
      else if (g >= int(m_labels.size()))
      {
        ++ nb_inliers[m_labels.size() + 1];
        if (!verbose)
          break;
      }
      else
        ++ nb_inliers[std::size_t(gt)];
      ++ total;
    }

    bool valid = (nb_inliers[m_labels.size() + 1] == 0);

    if (verbose)
    {
      std::cout << "Ground truth is " << (valid ? "valid" : "invalid") << ":" << std::endl;
      std::cout << " * " << nb_inliers[m_labels.size()] << " unclassified item(s) ("
                << 100. * (nb_inliers[m_labels.size()] / double(total)) << "%)" << std::endl;
      for (std::size_t i = 0; i < m_labels.size(); ++ i)
        std::cout << " * " << nb_inliers[i] << " " << m_labels[i]->name() << " inlier(s) ("
                  << 100. * (nb_inliers[i] / double(total)) << "%)" << std::endl;
      if (!valid)
        std::cout << " * " << nb_inliers[m_labels.size() + 1] << " item(s) with out-of-range index ("
                  << 100. * (nb_inliers[m_labels.size() + 1] / double(total)) << "%)" << std::endl;
    }

    return valid;
  }

  /// @}

};



} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_LABEL_SET_H
