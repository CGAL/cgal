// Copyright (C) 2014,2025 GeometryFactory
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_SURFACE_MESH_PROPERTY_H
#define CGAL_SURFACE_MESH_PROPERTY_H

#ifndef DOXYGEN_RUNNING

#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/tags.h>

#include <algorithm>
#include <atomic>
#include <optional>
#include <string>
#include <typeinfo>
#include <vector>

#ifdef CGAL_LINKED_WITH_TBB
#  include <tbb/concurrent_vector.h>
#endif

namespace CGAL {

namespace Properties {

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

    /// Grow the storage by n elements. Returns the index of the first new element.
    virtual size_t grow_by(size_t n) = 0;

    /// Free unused memory.
    virtual void shrink_to_fit() = 0;

    /// Extend the number of elements by one.
    virtual size_t push_back() = 0;

    /// Remove the last element.
    virtual void pop_back() = 0;

    /// Reset element to default value
    virtual void reset(size_t idx) = 0;

    virtual size_t capacity() const = 0;

    virtual bool transfer(const Base_property_array& other) = 0;
    virtual bool transfer(const Base_property_array& other, std::size_t from, std::size_t to) = 0;

    /// Let two elements swap their storage place.
    virtual void swap(size_t i0, size_t i1) = 0;

    /// Return a deep copy of self.
    virtual Base_property_array* clone () const = 0;

    /// Return an empty copy of self.
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

template <class T, class ConcurrencyTag = Sequential_tag>
class Property_array : public Base_property_array
{
public:

    typedef T                                       value_type;
    static constexpr bool is_parallel = std::is_convertible_v<ConcurrencyTag, Parallel_tag>;

#ifdef CGAL_LINKED_WITH_TBB
    using vector_type = std::conditional_t<is_parallel,
                                           tbb::concurrent_vector<value_type>,
                                           std::vector<value_type>>;
#else
    static_assert(false == is_parallel, "CGAL was not compiled with TBB support. "
                  "Please recompile CGAL with TBB support to use parallel features.");
    using vector_type = std::vector<value_type>;
#endif

    typedef typename vector_type::reference         reference;
    typedef typename vector_type::const_reference   const_reference;
    typedef typename vector_type::iterator          iterator;
    typedef typename vector_type::const_iterator    const_iterator;

    Property_array(const std::string& name, T t=T()) : Base_property_array(name), value_(t) {}

public: // virtual interface of Base_property_array

    size_t capacity() const final
    {
        return data_.capacity();
    }

    void reserve(size_t n) final
    {
        data_.reserve(n);
    }

    void resize(size_t n) final
    {
        data_.resize(n, value_);
    }

    size_t grow_by(size_t n) final
    {
      if constexpr(is_parallel) {
        return data_.grow_by(n, value_) - data_.begin();
      } else {
        const size_t old_size = data_.size();
        data_.resize(old_size + n, value_);
        return old_size;
      }
    }

    size_t push_back() final
    {
      if constexpr(is_parallel) {
        return data_.push_back(value_) - data_.begin();
      } else {
        data_.push_back(value_);
        return data_.size() - 1;
      }
    }

    void pop_back() final
    {
      if constexpr(is_parallel){
      }else{
        data_.pop_back();
      }
    }

    void reset(size_t idx) final
    {
        data_[idx] = value_;
    }

    bool transfer(const Base_property_array& other) final
    {
      const auto* pa = dynamic_cast<const Property_array*>(&other);
      if(pa != nullptr){
        std::copy((*pa).data_.begin(), (*pa).data_.end(), data_.end()-(*pa).data_.size());
        return true;
      }
      return false;
    }

    bool transfer(const Base_property_array& other, std::size_t from, std::size_t to) final
    {
      const auto* pa = dynamic_cast<const Property_array*>(&other);
      if (pa != nullptr)
      {
        data_[to] = (*pa)[from];
        return true;
      }

      return false;
    }

    void shrink_to_fit() final
    {
        vector_type(data_).swap(data_);
    }

    void swap(size_t i0, size_t i1) final
    {
        T d(data_[i0]);
        data_[i0]=data_[i1];
        data_[i1]=d;
    }

    Base_property_array* clone() const final
    {
        auto* p = new Property_array<T,ConcurrencyTag>(this->name_, this->value_);
        p->data_ = data_;
        return p;
    }

    Base_property_array* empty_clone() const final
    {
        auto* p = new Property_array<T,ConcurrencyTag>(this->name_, this->value_);
        return p;
    }

    const std::type_info& type() const final { return typeid(T); }


public:

    /// Get pointer to array (does not work for T==bool)
    const T* data() const
    {
        return &data_[0];
    }

    /// Access the i-th element. No range check is performed!
    reference operator[](std::size_t _idx)
    {
        CGAL_assertion( _idx < data_.size() );
        return data_[_idx];
    }

    /// Const access to the i-th element. No range check is performed!
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

template<typename, typename, typename ConcurrencyTag = Sequential_tag>
class Property_container;
/// @endcond




//== CLASS DEFINITION =========================================================
/// @cond CGAL_DOCUMENT_INTERNALS

template<typename Ref_class, typename Key, typename ConcurrencyTag>
class Property_container
{
public:
    static constexpr bool is_parallel = std::is_convertible_v<ConcurrencyTag, Parallel_tag>;

    // default constructor
    Property_container() = default;

    // destructor (deletes all property arrays)
    virtual ~Property_container() { clear(); }

    // copy constructor: performs deep copy of property arrays
    Property_container(const Property_container& _rhs) { operator=(_rhs); }

    Property_container(Property_container&& c) noexcept
    {
      c.swap(*this);
    }

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

    Property_container& operator=(Property_container&& c) noexcept
    {
      Property_container tmp(std::move(c));
      tmp.swap(*this);
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
    size_t size() const {
      if constexpr(is_parallel) {
        return size_.load(std::memory_order_acquire);
      } else {
       return size_;
      }
    }

    // returns the current capacity of the property arrays
    //size_t capacity() const { return capacity_; }

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
    std::optional<typename Get_pmap_type<T>::type>
    get(const std::string& name, std::size_t i) const
    {
      typedef typename Ref_class::template Get_property_map<Key, T>::type Pmap;
      if (parrays_[i]->name() == name)
      {
        if (auto* array = dynamic_cast<Property_array<T, ConcurrencyTag>*>(parrays_[i]))
          return std::optional(Pmap(array));
      }
      return std::nullopt;
    }

    // add a property with name \c name and default value \c t
    template <class T>
    std::pair<typename Get_pmap_type<T>::type, bool>
    add(const std::string& name, const T t=T())
    {
        typedef typename Ref_class::template Get_property_map<Key, T>::type Pmap;
        for (std::size_t i=0; i<parrays_.size(); ++i)
        {
            std::optional<Pmap> out = get<T>(name, i);
            if (out.has_value())
              return std::make_pair(*out, false);
        }

        // otherwise add the property
        auto* p = new Property_array<T, ConcurrencyTag>(name, t);
        p->reserve(capacity_);
        p->resize(size_);
        parrays_.push_back(p);
        return std::make_pair(Pmap(p), true);
    }


    // get a property by its name. Returns std::nullopt when it doesn't exist
    template <class T>
    std::optional<typename Get_pmap_type<T>::type>
    get(const std::string& name) const
    {
        typedef typename Ref_class::template Get_property_map<Key, T>::type Pmap;
        for (std::size_t i=0; i<parrays_.size(); ++i)
          {
            std::optional<Pmap> out = get<T>(name, i);
            if (out.has_value())
              return out;
          }
        return std::nullopt;
    }


    // returns a property if it exists, otherwise it creates it first.
    template <class T>
    typename Get_pmap_type<T>::type
    get_or_add(const std::string& name, const T t=T())
    {
      std::optional<typename Get_pmap_type<T>::type> out = get<T>(name);
      if (out.has_value())
        return out.value();
      else
        return add<T>(name, t).first;
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

    size_t capacity() const
    {
        if(parrays_.empty())
            return 0;
        return parrays_[0]->capacity();
    }

    size_t grow_by(size_t s) {
      if constexpr(is_parallel) {
        size_t pos = parrays_[0]->grow_by(s);
        for (std::size_t i=1, end = parrays_.size(); i<end; ++i) {
            parrays_[i]->grow_by(s);
        }
        size_.fetch_add(s, std::memory_order_relaxed);
        return pos;
      } else {
        size_t old_size = size_;
        size_ += s;
        for (std::size_t i=0, end = parrays_.size(); i<end; ++i) {
            parrays_[i]->resize(size_);
        }
        return old_size;
      }
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

    // increment the size of the property arrays by one
    // and return the old size
    size_t increment_size()
    {
      if constexpr(is_parallel) {
        auto old_size = size_.fetch_add(1, std::memory_order_relaxed);
        return old_size;
      } else {
        return size_++;
      }
    }

    void decrement_size()
    {
      if constexpr(is_parallel) {
        size_.fetch_sub(1, std::memory_order_relaxed);
      } else {
        --size_;
      }
    }

    // add a new element to each vector
    size_t push_back()
    {
      size_t pos = parrays_[0]->push_back();
      for (std::size_t i=1, end = parrays_.size(); i<end; ++i) {
          parrays_[i]->push_back();
      }
      auto old_size = increment_size();
      capacity_ = ((std::max)(old_size + 1, capacity_));
      return pos;
    }

    void pop_back()
    {
      for (std::size_t i=0; i<parrays_.size(); ++i) {
        parrays_[i]->pop_back();
      }
      decrement_size();
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
      size_t old_size = this->size_;
      this->size_ = other.size_;
      other.size_ = old_size;
      std::swap(this->capacity_, other.capacity_);
    }

private:
    std::vector<Base_property_array*>  parrays_;
    std::conditional_t<is_parallel,
                      std::atomic<size_t>,
                      size_t> size_ = 0;

    size_t  capacity_ = 0;
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
/// \cgalModels{LvaluePropertyMap}
///
template <class I, class T, class CRTP_derived_class, class ConcurrencyTag = Sequential_tag>
class Property_map_base
/// @cond CGAL_DOCUMENT_INTERNALS
  : public boost::put_get_helper<
           typename Property_array<T,ConcurrencyTag>::reference,
           CRTP_derived_class>
/// @endcond
{
public:
  using key_type = I;
  using value_type = T;
  using category = boost::lvalue_property_map_tag;
  using Property_array_type = Property_array<T,ConcurrencyTag>;

#ifndef DOXYGEN_RUNNING

  using reference = typename Property_array_type::reference;
  using const_reference = typename Property_array_type::const_reference;
  using iterator = typename Property_array_type::iterator;
  using const_iterator = typename Property_array_type::const_iterator;
#else
    /// A reference to the value type of the property.
  typedef unspecified_type reference;

    /// A const reference to the value type of the property.
  typedef unspecified_type const_reference;
#endif

#ifndef DOXYGEN_RUNNING
    template <typename Ref_class, typename Key, typename Tag>
    friend class Property_container;
#endif

public:
/// @cond CGAL_DOCUMENT_INTERNALS
    Property_map_base(Property_array_type* p=nullptr) : parray_(p) {}

    Property_map_base(Property_map_base&& pm) noexcept
      : parray_(std::exchange(pm.parray_, nullptr))
    {}

    Property_map_base(const Property_map_base& pm)
      : parray_(pm.parray_)
    {}

    Property_map_base& operator=(const Property_map_base& pm)
    {
      parray_ = pm.parray_;
      return *this;
    }

    Property_map_base& operator=(Property_map_base&& pm) noexcept
    {
      parray_ = std::exchange(pm.parray_, nullptr);
      return *this;
    }

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
    explicit operator bool() const {
        return parray_ != nullptr;
    }
#endif

    bool operator==(const Property_map_base& pm) const {
      return parray_ == pm.parray_;
    }

    bool operator!=(const Property_map_base& pm) const {
      return parray_ != pm.parray_;
    }

    /// Access the property associated with the key \c i.
    reference operator[](const I& i)
    {
      CGAL_assertion(parray_ != nullptr);
      return (*parray_)[i.id()];
    }

    /// Access the property associated with the key \c i.
    reference operator[](const I& i) const
    {
      CGAL_assertion(parray_ != nullptr);
      return (*parray_)[i.id()];
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

    void swap(Property_map_base& pm2)
    {
      std::swap(parray_, pm2.parray_);
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

    Property_array<T,ConcurrencyTag>& array()
    {
        CGAL_assertion(parray_ != nullptr);
        return *parray_;
    }

    const Property_array<T,ConcurrencyTag>& array() const
    {
        CGAL_assertion(parray_ != nullptr);
        return *parray_;
    }

    Property_array<T,ConcurrencyTag>* parray_;
};

#endif // DOXYGEN_RUNNING

} // Properties

} // CGAL

#endif // DOXYGEN_RUNNING

//=============================================================================
#endif // CGAL_SURFACE_MESH_PROPERTY_H
//=============================================================================
