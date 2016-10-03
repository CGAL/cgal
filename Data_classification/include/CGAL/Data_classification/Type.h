#ifndef CGAL_DATA_CLASSIFICATION_TYPE_H
#define CGAL_DATA_CLASSIFICATION_TYPE_H

namespace CGAL {

namespace Data_classification {

/*!
\ingroup PkgDataClassification

\brief Definition of a classification type based on a set of
relationship with attributes.

A classification type is used to segment the input data set. Usual
classification types are ground, vegetation and buildings (but others
can be defined).

*/
class Type
{
public:
  
  enum Attribute_effect /// Defines the effect of the values of an attribute on the classification type.
    {
      FAVORED_ATT = 0, ///< High values of the attribute favor this type
      NEUTRAL_ATT = 1, ///< The attribute has no effect on this type
      PENALIZED_ATT = 2 ///< Low values of the attribute favor this type
    };

private:
  std::string m_id;
  std::map<Attribute_handle, Attribute_effect> m_attribute_effects;
  std::vector<std::size_t> m_training_set;

public:

  /*! 
    \param id The name of the classification type
    (e.g. vegetation).
  */ 
  Type (std::string id) : m_id (id) { }

  /*! 
    \brief Sets how an attribute affects the classification type.

    \param att Attribute whose effect on the classification type will be set

    \param effect The effect the attribute will have on the classification type

  */ 
  void set_attribute_effect (Attribute_handle att, Attribute_effect effect)
  {
    m_attribute_effects[att] = effect;
  }

  /*!
    \brief Gets the effects of an attribute on the classification type.

    \param att Attribute

    \return The effect of the attribute on the classification type.
   */
  Attribute_effect attribute_effect (Attribute_handle att) 
  {
    std::map<Attribute_handle, Attribute_effect>::iterator
      search = m_attribute_effects.find (att);
    return (search == m_attribute_effects.end () ? NEUTRAL_ATT : search->second);
  }

  /*!
    \brief Gets the ID of the classification type.

    \return The ID of the classification type.
  */
  const std::string& id() const { return m_id; }
  
  /*!
    \brief Set of input point indices used as inlier of this
    classification type for training.

    \return The training set.
  */
  std::vector<std::size_t>& training_set() { return m_training_set; }

  /// \cond SKIP_IN_MANUAL
  void info()
  {
    std::cerr << "Attribute " << m_id << ": ";
    for (std::map<Attribute_handle, Attribute_effect>::iterator it = m_attribute_effects.begin();
         it != m_attribute_effects.end(); ++ it)
      {
        if (it->second == NEUTRAL_ATT)
          continue;
        
        std::cerr << it->first;
        if (it->second == FAVORED_ATT) std::cerr << " (favored), ";
        else if (it->second == PENALIZED_ATT) std::cerr << " (penalized), ";
      }
    std::cerr << std::endl;
  }
  /// \endcond

};

typedef boost::shared_ptr<Type> Type_handle;

} // namespace Data_classification

} // namespace CGAL

#endif // CGAL_DATA_CLASSIFICATION_TYPE_H
