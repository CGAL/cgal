//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public License
// as published by the Free Software Foundation, version 2.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//=============================================================================


#ifndef SURFACE_MESH_PROPERTY_H
#define SURFACE_MESH_PROPERTY_H


//== INCLUDES =================================================================


#include <vector>
#include <string>
#include <algorithm>
#include <typeinfo>


//== CLASS DEFINITION =========================================================


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
    virtual void free_memory() = 0;

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



//== CLASS DEFINITION =========================================================


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

    virtual void free_memory()
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
        assert( size_t(_idx) < data_.size() );
        return data_[_idx];
    }

    /// Const access to the i'th element. No range check is performed!
    const_reference operator[](int _idx) const
    {
        assert( size_t(_idx) < data_.size());
        return data_[_idx];
    }



private:
    vector_type data_;
    value_type  value_;
};


// specialization for bool properties
template <>
inline const bool*
Property_array<bool>::data() const
{
    assert(false);
    return NULL;
}



//== CLASS DEFINITION =========================================================


template <class T>
class Property
{
public:

    typedef typename Property_array<T>::reference reference;
    typedef typename Property_array<T>::const_reference const_reference;

    friend class Property_container;
    friend class Surface_mesh;


public:

    Property(Property_array<T>* p=NULL) : parray_(p) {}

    void reset()
    {
        parray_ = NULL;
    }

    operator bool() const
    {
        return parray_ != NULL;
    }

    reference operator[](int i)
    {
        assert(parray_ != NULL);
        return (*parray_)[i];
    }

    const_reference operator[](int i) const
    {
        assert(parray_ != NULL);
        return (*parray_)[i];
    }

    const T* data() const
    {
        assert(parray_ != NULL);
        return parray_->data();
    }


private:

    Property_array<T>& array()
    {
        assert(parray_ != NULL);
        return *parray_;
    }

    const Property_array<T>& array() const
    {
        assert(parray_ != NULL);
        return *parray_;
    }


private:
    Property_array<T>* parray_;
};



//== CLASS DEFINITION =========================================================


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
    template <class T> Property<T> add(const std::string& name, const T t=T())
    {
        // if a property with this name already exists, return an invalid property
        for (unsigned int i=0; i<parrays_.size(); ++i)
        {
            if (parrays_[i]->name() == name)
            {
                std::cerr << "[Property_container] A property with name \""
                          << name << "\" already exists. Returning invalid property.\n";
                return Property<T>();
            }
        }

        // otherwise add the property
        Property_array<T>* p = new Property_array<T>(name, t);
        p->resize(size_);
        parrays_.push_back(p);
        return Property<T>(p);
    }


    // get a property by its name. returns invalid property if it does not exist.
    template <class T> Property<T> get(const std::string& name) const
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            if (parrays_[i]->name() == name)
                return Property<T>(dynamic_cast<Property_array<T>*>(parrays_[i]));
        return Property<T>();
    }


    // returns a property if it exists, otherwise it creates it first.
    template <class T> Property<T> get_or_add(const std::string& name, const T t=T())
    {
        Property<T> p = get<T>(name);
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
    template <class T> void remove(Property<T>& h)
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
    void free_memory() const
    {
        for (unsigned int i=0; i<parrays_.size(); ++i)
            parrays_[i]->free_memory();
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


//=============================================================================
#endif // SURFACE_MESH_PROPERTY_H
//=============================================================================
