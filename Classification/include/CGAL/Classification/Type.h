#ifndef CGAL_CLASSIFICATION_TYPE_H
#define CGAL_CLASSIFICATION_TYPE_H

namespace CGAL {

namespace Classification {

/*!
\ingroup PkgClassification

\brief %Classification type (for example: vegetation, ground, etc.)
defined as a set of relationship with classification attributes.

*/
class Type
{
public:
  
private:
  /// \cond SKIP_IN_MANUAL
  std::string m_name;
  std::map<Attribute_handle, Attribute_effect> m_attribute_effects;
  /// \endcond

public:

  /*! 
    \param name The name of the classification type
    (e.g. vegetation).
  */ 
  Type (std::string name) : m_name (name) { }

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
    return (search == m_attribute_effects.end () ? NEUTRAL : search->second);
  }

  const std::string& name() const { return m_name; }
  
  /// \cond SKIP_IN_MANUAL
  void info()
  {
    std::cerr << "Attribute " << m_name << ": ";
    for (std::map<Attribute_handle, Attribute_effect>::iterator it = m_attribute_effects.begin();
         it != m_attribute_effects.end(); ++ it)
      {
        if (it->second == NEUTRAL)
          continue;
        
        std::cerr << it->first;
        if (it->second == FAVORING) std::cerr << " (favored), ";
        else if (it->second == PENALIZING) std::cerr << " (penalized), ";
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
