//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
// Copyright (C) 2014 GeometryFactory
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//


#ifndef CGAL_SURFACE_MESH_PROPERTY_H
#define CGAL_SURFACE_MESH_PROPERTY_H

#include <vector>
#include <string>
#include <algorithm>
#include <typeinfo>

#include <boost/property_map/property_map.hpp>

#include <CGAL/assertions.h>

namespace CGAL {

#ifndef DOXYGEN_RUNNING

/// \addtogroup PkgSurface_mesh
///
/// @{

/// @cond CGAL_DOCUMENT_INTERNALS
class Base_property_array
{
public:

    /// Default constructor
    Base_property_array(const std::string& name) : name_(name) {}

    /// Destructor.
    virtual ~Base_property_array() {}

    /// Reserve memory for n elements.
    virtual void reserve(size_t n) = 0;

    /// Resize storage to hold n elements.
    virtual void resize(size_t n) = 0;

    /// Free unused memory.
    virtual void shrink_to_fit() = 0;

    /// Extend the number of elements by one.
    virtual void push_back() = 0;

    /// Let two elements swap their storage place.
    virtual void swap(size_t i0, size_t i1) = 0;

    /// Return a deep copy of self.
    virtual Base_property_array* clone () const = 0;

    /// Return the type_info of the property
    virtual const std::type_info& type() = 0;

    /// Return the name of the property
    const std::string& name() const { return name_; }


protected:

    std::string name_;
};

  /// @endcond


//== CLASS DEFINITION =========================================================

/// @cond CGAL_DOCUMENT_INTERNALS

template <class T>
class Property_array : public Base_property_array
{
public:

    typedef T                                       value_type;
    typedef std::vector<value_type>                 vector_type;
    typedef typename vector_type::reference         reference;
    typedef typename vector_type::const_reference   const_reference;

    Property_array(const std::string& name, T t=T()) : Base_property_array(name), value_(t) {}

public: // virtual interface of Base_property_array

    virtual void reserve(size_t n)
    {
        data_.reserve(n);
    }

    virtual void resize(size_t n)
    {
        data_.resize(n, value_);
    }

    virtual void push_back()
    {
        data_.push_back(value_);
    }

    virtual void shrink_to_fit()
    {
        vector_type(data_).swap(data_);
    }

    virtual void swap(size_t i0, size_t i1)
    {
        T d(data_[i0]);
        data_[i0]=data_[i1];
        data_[i1]=d;
    }

    virtual Base_property_array* clone() const
    {
        Property_array<T>* p = new Property_array<T>(name_, value_);
        p->data_ = data_;
        return p;
    }

    virtual const std::type_info& type() { return typeid(T); }


public:

    /// Get pointer to array (does not work for T==bool)
    const T* data() const
    {
        return &data_[0];
    }

    /// Access the i'th element. No range check is performed!
    reference operator[](int _idx)
    {
        CGAL_assertion( size_t(_idx) < data_.size() );
        return data_[_idx];
    }

    /// Const access to the i'th element. No range check is performed!
    const_reference operator[](int _idx) const
    {
        CGAL_assertion( size_t(_idx) < data_.size());
        return data_[_idx];
    }



private:
    vector_type data_;
    value_type  value_;
};

  /// @endcond

// specialization for bool properties
template <>
inline const bool*
Property_array<bool>::data() const
{
    CGAL_assertion(false);
    return NULL;
}



//== CLASS DEFINITION =========================================================

/// @cond CGAL_DOCUMENT_INTERNALS

template<typename>
class Property_container;
/// @endcond 




//== CLASS DEFINITION =========================================================
/// @cond CGAL_DOCUMENT_INTERNALS

template <class, class>
class Property_map;

template<typename Key>
class Property_container
{
public:

    // default constructor
    Property_container() : size_(0) {}

    // destructor (deletes all property arrays)
    virtual ~Property_container() { clear(); }

    // copy constructor: performs deep copy of property arrays
    Property_container(const Property_container& _rhs) { operator=(_rhs); }

    // assignment: performs deep copy of property arrays
    Property_container& operator=(const Property_container& _rhs)
    {
        if (this != &_rhs)
        {
            clear();
            parrays_.resize(_rhs.n_properties());
            size_ = _rhs.size();
            for (unsigned int i=0; i<parrays_.size(); ++i)
                parrays_[i] = _rhs.parrays_[i]->clone();
        }
        return *this;
    }

    // returns the current size of the property arrays
    size_t size() const { return size_; }

    // returns the number of property arrays
    size_t n_properties() const { return parrays_.size(); }

    // returns a vector of all property names
    std::vector<std::string> properties() const
    {
        std::vector<std::string> names;
        for (unsigned int i=0; i<parrays_.size(); ++i)
            names.push_back(parrays_[i]->name());
        return names;
    }

    // add a property with name \c name and default value \c t
    template <class T> 
    std::pair<Property_map<Key, T>,bool>
    add(const std::string& name, const T t=T())
    {
        // if a property with this name already exists, return it
        for (unsigned int i=0; i<parrays_.size(); ++i)
        {
            if (parrays_[i]->name() == name)
            {
              return std::make_pair(Property_map<Key, T>(dynamic_cast<Property_array<T>*>(parrays_[i])), true);
            }
        }

        // otherwise add the property
        Property_array<T>* p = new Property_array<T>(name, t);
        p->resize(size_);
        parrays_.push_back(p);
        return std::make_pair(Property_map<Key, T>(p),false);
    }


    // get a property by its name. returns invalid property if it does not exist.
    template <class T> Property_map<Key, T> get(const std::string& name) const
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            if (parrays_[i]->name() == name)
                return Property_map<Key, T>(dynamic_cast<Property_array<T>*>(parrays_[i]));
        return Property_map<Key, T>();
    }


    // returns a property if it exists, otherwise it creates it first.
    template <class T> Property_map<Key, T> get_or_add(const std::string& name, const T t=T())
    {
        Property_map<Key, T> p = get<T>(name);
        if (!p) p = add<T>(name, t);
        return p;
    }


    // get the type of property by its name. returns typeid(void) if it does not exist.
    const std::type_info& get_type(const std::string& name)
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            if (parrays_[i]->name() == name)
                return parrays_[i]->type();
        return typeid(void);
    }


    // delete a property
    template <class T> void remove(Property_map<Key, T>& h)
    {
        std::vector<Base_property_array*>::iterator it=parrays_.begin(), end=parrays_.end();
        for (; it!=end; ++it)
        {
            if (*it == h.parray_)
            {
                delete *it;
                parrays_.erase(it);
                h.reset();
                break;
            }
        }
    }


    // delete all properties
    void clear()
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            delete parrays_[i];
        parrays_.clear();
        size_ = 0;
    }


    // reserve memory for n entries in all arrays
    void reserve(size_t n) const
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->reserve(n);
    }

    // resize all arrays to size n
    void resize(size_t n)
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->resize(n);
        size_ = n;
    }

    // free unused space in all arrays
    void shrink_to_fit() const
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->shrink_to_fit();
    }

    // add a new element to each vector
    void push_back()
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->push_back();
        ++size_;
    }

    // swap elements i0 and i1 in all arrays
    void swap(size_t i0, size_t i1) const
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->swap(i0, i1);
    }


private:
    std::vector<Base_property_array*>  parrays_;
    size_t  size_;
};

  /// @endcond

///
/// `Property_map` enables to attach properties to the simplices of a 
///  surface mesh.
/// 
/// @tparam Key The key type of the property map. It must be a model of `Index`.
/// @tparam Value The value type of the property.
///
/// \cgalModels `LvaluePropertyMap`
///
template <class I, class T>
class Property_map
/// @cond CGAL_DOCUMENT_INTERNALS
  : public boost::put_get_helper< 
           typename Property_array<T>::reference,
           Property_map< I, T > >
/// @endcond
{
    typedef void (Property_map::*bool_type)() const;
    void this_type_does_not_support_comparisons() const {}
public:
    typedef I key_type;
    typedef T value_type;
    typedef boost::lvalue_property_map_tag category;

#ifndef DOXYGEN_RUNNING

    typedef typename Property_array<T>::reference reference;

    typedef typename Property_array<T>::const_reference const_reference;
#else 
    /// A reference to the value type of the property.
  typedef unspecified_type reference;

    /// A const reference to the value type of the property.
  typedef unspecified_type const_reference;
#endif

#ifndef DOXYGEN_RUNNING
    friend class Property_container<I>;

    template <typename K>  friend class Surface_mesh;
#endif

public:
/// @cond CGAL_DOCUMENT_INTERNALS
    Property_map(Property_array<T>* p=NULL) : parray_(p) {}

    void reset()
    {
        parray_ = NULL;
    }
  /// @endcond 

public:
    /// \name Accessing Properties
    //@{
#ifdef DOXYGEN_RUNNING
    /// Conversion to a Boolean. It is \c true when the property map
    /// can be used, and \c false otherwise.  
  operator bool () const;
#else
    operator bool_type() const {
        return parray_ != NULL ?
            &Property_map::this_type_does_not_support_comparisons : 0;
    }
#endif
    /// Access the property associated with the key \c i.
    reference operator[](const I& i)
    {
      CGAL_assertion(parray_ != NULL);
      return (*parray_)[i.idx()];
    }

    /// Access the property associated with the key \c i.
    reference operator[](const I& i) const
    {
      CGAL_assertion(parray_ != NULL);
      return (*parray_)[i.idx()];
    }

    /// Allows access to the underlying storage of the property. This
    /// is useful when the key associated with the properties is
    /// unimportant and only the properties are of interest
    /// (e.g. rendering).
    ///
    /// \returns a pointer to the underlying storage of the property.
    const T* data() const
    {
      CGAL_assertion(parray_ != NULL);
      return parray_->data();
    }

    //@}
private:

    Property_array<T>& array()
    {
        CGAL_assertion(parray_ != NULL);
        return *parray_;
    }

    const Property_array<T>& array() const
    {
        CGAL_assertion(parray_ != NULL);
        return *parray_;
    }

    Property_array<T>* parray_;
};



///@}
#endif // DOXYGEN_RUNNING  
} // CGAL

//=============================================================================
#endif // CGAL_SURFACE_MESH_PROPERTY_H
//=============================================================================
