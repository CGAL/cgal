// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
// Copyright (C) 2014 GeometryFactory
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//


#ifndef CGAL_SURFACE_MESH_PROPERTY_H
#define CGAL_SURFACE_MESH_PROPERTY_H

#include <CGAL/license/Surface_mesh.h>

#ifndef DOXYGEN_RUNNING

#include <vector>
#include <string>
#include <algorithm>
#include <typeinfo>

#include <CGAL/property_map.h>
#include <CGAL/assertions.h>

namespace CGAL {

namespace Properties {

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

    /// Reset element to default value
    virtual void reset(size_t idx) = 0;

    virtual bool transfer(const Base_property_array& other) = 0;
    virtual bool transfer(const Base_property_array& other, std::size_t from, std::size_t to) = 0;

    /// Let two elements swap their storage place.
    virtual void swap(size_t i0, size_t i1) = 0;

    /// Return a deep copy of self.
    virtual Base_property_array* clone () const = 0;

    /// Return a empty copy of self.
    virtual Base_property_array* empty_clone () const = 0;

    /// Return the type_info of the property
    virtual const std::type_info& type() const = 0;

    /// Return the name of the property
    const std::string& name() const { return name_; }

    bool is_same (const Base_property_array& other)
    {
      return (name() == other.name() && type() == other.type());
    }

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
    typedef typename vector_type::iterator          iterator;
    typedef typename vector_type::const_iterator    const_iterator;

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

    virtual void reset(size_t idx)
    {
        data_[idx] = value_;
    }

    bool transfer(const Base_property_array& other)
    {
      const Property_array<T>* pa = dynamic_cast<const Property_array*>(&other);
      if(pa != nullptr){
        std::copy((*pa).data_.begin(), (*pa).data_.end(), data_.end()-(*pa).data_.size());
        return true;
      }
      return false;
    }

    bool transfer(const Base_property_array& other, std::size_t from, std::size_t to)
    {
      const Property_array<T>* pa = dynamic_cast<const Property_array*>(&other);
      if (pa != nullptr)
      {
        data_[to] = (*pa)[from];
        return true;
      }

      return false;
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
        Property_array<T>* p = new Property_array<T>(this->name_, this->value_);
        p->data_ = data_;
        return p;
    }

    virtual Base_property_array* empty_clone() const
    {
        Property_array<T>* p = new Property_array<T>(this->name_, this->value_);
        return p;
    }

    virtual const std::type_info& type() const { return typeid(T); }


public:

    /// Get pointer to array (does not work for T==bool)
    const T* data() const
    {
        return &data_[0];
    }

    /// Access the i'th element. No range check is performed!
    reference operator[](std::size_t _idx)
    {
        CGAL_assertion( _idx < data_.size() );
        return data_[_idx];
    }

    /// Const access to the i'th element. No range check is performed!
    const_reference operator[](std::size_t _idx) const
    {
        CGAL_assertion( _idx < data_.size());
        return data_[_idx];
    }

    iterator begin() { return data_.begin(); }
    iterator end() { return data_.end(); }
    const_iterator begin() const { return data_.begin(); }
    const_iterator end() const { return data_.end(); }

private:
    vector_type data_;
    value_type  value_;
};


  /// @endcond

//== CLASS DEFINITION =========================================================

/// @cond CGAL_DOCUMENT_INTERNALS

template<typename, typename>
class Property_container;
/// @endcond




//== CLASS DEFINITION =========================================================
/// @cond CGAL_DOCUMENT_INTERNALS

template<typename Ref_class, typename Key>
class Property_container
{
public:

    // default constructor
    Property_container() : size_(0), capacity_(0) {}

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
            capacity_ = _rhs.capacity();
            for (std::size_t i=0; i<parrays_.size(); ++i)
                parrays_[i] = _rhs.parrays_[i]->clone();
        }
        return *this;
    }

    void transfer(const Property_container& _rhs)
    {
      for(std::size_t i=0; i<parrays_.size(); ++i){
        for (std::size_t j=0; j<_rhs.parrays_.size(); ++j){
          if(parrays_[i]->is_same (*(_rhs.parrays_[j]))){
            parrays_[i]->transfer(* _rhs.parrays_[j]);
            break;
          }
        }
      }
    }

    // Copy properties that don't already exist from another container
    void copy_properties (const Property_container& _rhs)
    {
      for (std::size_t i = 0; i < _rhs.parrays_.size(); ++ i)
      {
        bool property_already_exists = false;
        for (std::size_t j = 0; j < parrays_.size(); ++ j)
          if (_rhs.parrays_[i]->is_same (*(parrays_[j])))
          {
            property_already_exists = true;
            break;
          }

        if (property_already_exists)
          continue;

        parrays_.push_back (_rhs.parrays_[i]->empty_clone());
        parrays_.back()->reserve(capacity_);
        parrays_.back()->resize(size_);
      }
    }

    // Transfer one element with all properties
    // WARNING: properties must be the same in the two containers
    bool transfer(const Property_container& _rhs, std::size_t from, std::size_t to)
    {
      bool out = true;
      for(std::size_t i=0; i<parrays_.size(); ++i)
        if (!(parrays_[i]->transfer(* _rhs.parrays_[i], from, to)))
          out = false;
      return out;
    }

    // returns the current size of the property arrays
    size_t size() const { return size_; }

    // returns the current capacity of the property arrays
    size_t capacity() const { return capacity_; }

    // returns the number of property arrays
    size_t n_properties() const { return parrays_.size(); }

    // returns a vector of all property names
    std::vector<std::string> properties() const
    {
        std::vector<std::string> names;
        for (std::size_t i=0; i<parrays_.size(); ++i)
            names.push_back(parrays_[i]->name());
        return names;
    }

    template <typename T>
    struct Get_pmap_type {
      typedef typename Ref_class::template Get_property_map<Key, T>::type type;
    };

    template <class T>
    std::pair<typename Get_pmap_type<T>::type, bool>
    get(const std::string& name, std::size_t i) const
    {
      typedef typename Ref_class::template Get_property_map<Key, T>::type Pmap;
      if (parrays_[i]->name() == name)
        {
          if (Property_array<T>* array = dynamic_cast<Property_array<T>*>(parrays_[i]))
            return std::make_pair (Pmap(array), true);
        }
      return std::make_pair(Pmap(), false);
    }

    // add a property with name \c name and default value \c t
    template <class T>
    std::pair<typename Get_pmap_type<T>::type, bool>
    add(const std::string& name, const T t=T())
    {
        typedef typename Ref_class::template Get_property_map<Key, T>::type Pmap;
        for (std::size_t i=0; i<parrays_.size(); ++i)
        {
            std::pair<Pmap, bool> out = get<T>(name, i);
            if (out.second)
              {
                out.second = false;
                return out;
              }
        }

        // otherwise add the property
        Property_array<T>* p = new Property_array<T>(name, t);
        p->reserve(capacity_);
        p->resize(size_);
        parrays_.push_back(p);
        return std::make_pair(Pmap(p), true);
    }


    // get a property by its name. returns invalid property if it does not exist.
    template <class T>
    std::pair<typename Get_pmap_type<T>::type, bool>
    get(const std::string& name) const
    {
        typedef typename Ref_class::template Get_property_map<Key, T>::type Pmap;
        for (std::size_t i=0; i<parrays_.size(); ++i)
          {
            std::pair<Pmap, bool> out = get<T>(name, i);
            if (out.second)
              return out;
          }
        return std::make_pair(Pmap(), false);
    }


    // returns a property if it exists, otherwise it creates it first.
    template <class T>
    typename Get_pmap_type<T>::type
    get_or_add(const std::string& name, const T t=T())
    {
      typename Ref_class::template Get_property_map<Key, T>::type p;
      bool b;
      boost::tie(p,b)= get<T>(name);
        if (!b) p = add<T>(name, t).first;
        return p;
    }


    // get the type of property by its name. returns typeid(void) if it does not exist.
    const std::type_info&
    get_type(const std::string& name) const
    {
        for (std::size_t i=0; i<parrays_.size(); ++i)
            if (parrays_[i]->name() == name)
                return parrays_[i]->type();
        return typeid(void);
    }


    // delete a property
    template <class T>
    bool
    remove(typename Get_pmap_type<T>::type& h)
    {
        typename std::vector<Base_property_array*>::iterator it=parrays_.begin(), end=parrays_.end();
        for (; it!=end; ++it)
        {
            if (*it == h.parray_)
            {
                delete *it;
                parrays_.erase(it);
                h.reset();
                return true;
            }
        }
        return false;
    }


    // delete all properties
    void clear()
    {
        for (std::size_t i=0; i<parrays_.size(); ++i)
            delete parrays_[i];
        parrays_.clear();
        size_ = 0;
    }


    // reserve memory for n entries in all arrays
    void reserve(size_t n)
    {
        for (std::size_t i=0; i<parrays_.size(); ++i)
            parrays_[i]->reserve(n);
        capacity_ = (std::max)(n, capacity_);
    }

    // resize all arrays to size n
    void resize(size_t n)
    {
        for (std::size_t i=0; i<parrays_.size(); ++i)
            parrays_[i]->resize(n);
        size_ = n;
    }

    // resize the vector of properties to n, deleting all other properties
    void resize_property_array(size_t n)
    {
        if (parrays_.size()<=n)
          return;
        for (std::size_t i=n; i<parrays_.size(); ++i)
            delete parrays_[i];
        parrays_.resize(n);
    }

    // free unused space in all arrays
    void shrink_to_fit()
    {
        for (std::size_t i=0; i<parrays_.size(); ++i)
            parrays_[i]->shrink_to_fit();
        capacity_ = size_;
    }

    // add a new element to each vector
    void push_back()
    {
        for (std::size_t i=0; i<parrays_.size(); ++i)
            parrays_[i]->push_back();
        ++size_;
        capacity_ = ((std::max)(size_, capacity_));
    }

    // reset element to its default property values
    void reset(size_t idx)
    {
        for (std::size_t i=0; i<parrays_.size(); ++i)
            parrays_[i]->reset(idx);
    }

    // swap elements i0 and i1 in all arrays
    void swap(size_t i0, size_t i1) const
    {
        for (std::size_t i=0; i<parrays_.size(); ++i)
            parrays_[i]->swap(i0, i1);
    }

    // swap content with other Property_container
    void swap (Property_container& other)
    {
      this->parrays_.swap (other.parrays_);
      std::swap(this->size_, other.size_);
    }

private:
    std::vector<Base_property_array*>  parrays_;
    size_t  size_;
    size_t  capacity_;
};

  /// @endcond

#ifndef DOXYGEN_RUNNING
///
///
/// `Property_map` enables to attach properties to the simplices of a
///  surface mesh.
///
/// @tparam Key The key type of the property map. It must be a model of `Index`.
/// @tparam Value The value type of the property.
///
/// \cgalModels `LvaluePropertyMap`
///
template <class I, class T, class CRTP_derived_class>
class Property_map_base
/// @cond CGAL_DOCUMENT_INTERNALS
  : public boost::put_get_helper<
           typename Property_array<T>::reference,
           CRTP_derived_class>
/// @endcond
{
    typedef void (Property_map_base::*bool_type)() const;
    void this_type_does_not_support_comparisons() const {}
public:
    typedef I key_type;
    typedef T value_type;
    typedef boost::lvalue_property_map_tag category;

#ifndef DOXYGEN_RUNNING

    typedef typename Property_array<T>::reference reference;
    typedef typename Property_array<T>::const_reference const_reference;
    typedef typename Property_array<T>::iterator iterator;
    typedef typename Property_array<T>::const_iterator const_iterator;
#else
    /// A reference to the value type of the property.
  typedef unspecified_type reference;

    /// A const reference to the value type of the property.
  typedef unspecified_type const_reference;
#endif

#ifndef DOXYGEN_RUNNING
    template <typename Ref_class, typename Key>
    friend class Property_container;
#endif

public:
/// @cond CGAL_DOCUMENT_INTERNALS
    Property_map_base(Property_array<T>* p=nullptr) : parray_(p) {}

    void reset()
    {
        parray_ = nullptr;
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
        return parray_ != nullptr ?
            &Property_map_base::this_type_does_not_support_comparisons : 0;
    }
#endif
    /// Access the property associated with the key \c i.
    reference operator[](const I& i)
    {
      CGAL_assertion(parray_ != nullptr);
      return (*parray_)[i];
    }

    /// Access the property associated with the key \c i.
    reference operator[](const I& i) const
    {
      CGAL_assertion(parray_ != nullptr);
      return (*parray_)[i];
    }

    iterator begin() { return parray_->begin(); }
    iterator end() { return parray_->end(); }
    const_iterator begin() const { return parray_->begin(); }
    const_iterator end() const { return parray_->end(); }

    bool transfer (const Property_map_base& other)
    {
      return parray_->transfer(*(other.parray_));
    }

    bool transfer (const Property_map_base& other, std::size_t from, std::size_t to)
    {
      return parray_->transfer(*(other.parray_), from, to);
    }

    /// Allows access to the underlying storage of the property. This
    /// is useful when the key associated with the properties is
    /// unimportant and only the properties are of interest
    /// (e.g. rendering).
    ///
    /// \returns a pointer to the underlying storage of the property.
    const T* data() const
    {
      CGAL_assertion(parray_ != nullptr);
      return parray_->data();
    }

    //@}
#ifndef CGAL_TEST_SURFACE_MESH
private:
#endif

    Property_array<T>& array()
    {
        CGAL_assertion(parray_ != nullptr);
        return *parray_;
    }

    const Property_array<T>& array() const
    {
        CGAL_assertion(parray_ != nullptr);
        return *parray_;
    }

    Property_array<T>* parray_;
};

#endif // DOXYGEN_RUNNING

///@}

} // Properties

} // CGAL

#endif // DOXYGEN_RUNNING

//=============================================================================
#endif // CGAL_SURFACE_MESH_PROPERTY_H
//=============================================================================
