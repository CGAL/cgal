#ifndef CGAL_CLASSIFICATION_TYPE_H
#define CGAL_CLASSIFICATION_TYPE_H

namespace CGAL {

namespace Classification {

/*!
\ingroup PkgClassification

\brief Classification type.

A type is what is sought after when performing classification on an
input set (for example: vegetation, ground, etc.). It is defined
as a set of relationship with classification attributes.

*/
class Type
{
public:
  
  enum Attribute_effect /// Defines how an attribute affects this type.
    {
      FAVORED_ATT = 0, ///< High values of the attribute favor this type
      NEUTRAL_ATT = 1, ///< The attribute has no effect on this type
      PENALIZED_ATT = 2 ///< Low values of the attribute favor this type
    };

private:
  /// \cond SKIP_IN_MANUAL
  std::string m_id;
  std::map<Attribute_handle, Attribute_effect> m_attribute_effects;
  /// \endcond

public:

  /*! 
    \param id The name of the classification type
    (e.g. vegetation).
  */ 
  Type (std::string id) : m_id (id) { }

  /*! 
    \brief Sets how an attribute affects the classification type.

    \param att Attribute whose effect on the classification type is set

    \param effect The effect the attribute has on the classification type

  */ 
  void set_attribute_effect (Attribute_handle att, Attribute_effect effect)
  {
    m_attribute_effects[att] = effect;
  }

  /*!
    \brief Returns the effect of attribute `att` on the classification type.
   */
  Attribute_effect attribute_effect (Attribute_handle att) 
  {
    std::map<Attribute_handle, Attribute_effect>::iterator
      search = m_attribute_effects.find (att);
    return (search == m_attribute_effects.end () ? NEUTRAL_ATT : search->second);
  }

  /*!
    \brief Returns the ID of the classification type.
  */
  const std::string& id() const { return m_id; }
  
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

#ifdef DOXYGEN_RUNNING
/*!
  \ingroup PkgClassification

  \brief Handle to a classification `Type`.

  \cgalModels Handle
*/
class Type_handle { };
#else
typedef boost::shared_ptr<Type> Type_handle;
#endif

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_TYPE_H
