// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
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
// $URL$
// $Id$
// 
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
#ifndef CGAL_OPEN_HASH_H
#define CGAL_OPEN_HASH_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the Open_hash and the Hash_map class-templates.
 */

#include <vector>
#include <list>

namespace CGAL {

/*! \class
 * A default equality functor.
 * The parameter Type must support the equality operator. 
 */
template <class Type>
class Default_is_equal_functor
{
public:

  bool operator() (const Type& t1, const Type& t2) const
  {
    return (t1 == t2);
  }
};

/*! \class
 * An implementation of an open-addressing hash.
 */
template <class Key_,
          class Hash_functor_,
          class EqualKey_ = Default_is_equal_functor<Key_> >
class Open_hash
{
public:

  // Type definitions (for STL compatibility):
  typedef Key_               value_type;
  typedef Key_               key_type;
  typedef Hash_functor_      hasher;
  typedef EqualKey_          key_equal;
  typedef value_type*        pointer;
  typedef value_type&        reference;
  typedef const value_type*  const_pointer;
  typedef const value_type&  const_reference;
  typedef size_t             size_type;
  typedef size_t             difference_type;

protected:

  enum
  {
    DEFAULT_NUMBER_OF_BUCKETS = 1000
  };

  typedef std::list<Key_>                  Bucket;
  typedef typename Bucket::iterator        Bucket_iterator;
  typedef typename Bucket::const_iterator  Bucket_const_iterator;
  
  size_type            n_buckets;   // Number of buckets in the hash.
  size_type            n_objects;   // Number of objects stored in the hash.
  std::vector<Bucket>  buckets;     // The buckets (linked lists of objects).
  hasher               hash_func;   // The hashing functor.
  key_equal            equal_func;  // The equality functor.

private:

  // Copy constructor and assignment operator - not supported.
  Open_hash (const Open_hash& );
  Open_hash& operator= (const Open_hash& );

public:

  /// \name Construction and destruction methods.
  //@{

  /*! Default constructor. */
  Open_hash () :
    n_objects (0),
    buckets (DEFAULT_NUMBER_OF_BUCKETS),
    hash_func (),
    equal_func ()
  {
    n_buckets = buckets.size();
  }

  /*! Constructor with an initial number of buckets. */
  Open_hash (size_type n,
             hasher hash_f = hasher(),
             key_equal equal_f = key_equal()) :
    n_objects (0),
    buckets (n),
    hash_func (hash_f),
    equal_func (equal_f)
  {
    n_buckets = buckets.size();
  }

  /*! Destructor. */
  virtual ~Open_hash ()
  {}
  //@}

  /// \number Access the hash properties.
  //@{

  /*! Get the number of buckets in the hash. */
  size_type bucket_count () const
  {
    return (n_buckets);
  }

  /*! Get the number of objects stored in the hash. */
  size_type size () const
  {
    return (n_objects);
  }

  /*! Check whether the hash is empty. */
  bool empty () const
  {
    return (n_objects == 0);
  }

  /*! Get the hasher object used by the hash. */ 
  const hasher& hash_funct () const
  {
    return (hash_func);
  }

  /*! Get the key_equal object used by the hash. */
  const key_equal& key_eq () const
  {
    return (equal_func);
  }
  //@}

  /*! \class
   * An iterator for the open-addressing hash.
   */
  class iterator
  {
    friend class Open_hash<Key_, Hash_functor_, EqualKey_>;

  private:
  
    Open_hash<Key_, Hash_functor_, EqualKey_>  *p_hash;  // The scanned hash.
    int                           index;  // The index of the current bucket.
    Bucket_iterator               iter;   // Iterator for the current bucket.

    /*! Private constructor. */
    iterator (const Open_hash<Key_, Hash_functor_, EqualKey_>& hash,
              size_type ind) :
      p_hash (const_cast<Open_hash<Key_,Hash_functor_,EqualKey_>* >(&hash))
    {
      // Find the next non-empty bucket and get an iterator for it.
      index = p_hash->_next_non_empty_bucket (static_cast<int> (ind));

      if (index < static_cast<int>(p_hash->n_buckets))
        iter = (p_hash->buckets[index]).begin();
    }

    /*! Private constructor. */
    iterator (const Open_hash<Key_, Hash_functor_, EqualKey_>& hash,
              size_type ind,
              Bucket_iterator it) :
      p_hash (const_cast<Open_hash<Key_,Hash_functor_,EqualKey_>* >(&hash)),
      index (static_cast<int> (ind)),
      iter (it)
    {}

  public:

    /*! Default constructor. */
    iterator () :
      p_hash (NULL),
      index (0)
    {}

    /*! Equality operator. */
    bool operator== (const iterator& it)
    {
      if (p_hash != it.p_hash)
        return (false);
      
      if (p_hash == NULL)
        return (true);

      if ((index < 0 || index >= static_cast<int>(p_hash->n_buckets)) &&
          (it.index < 0 || it.index >= static_cast<int>(p_hash->n_buckets)))
        return (true);

      return (index == it.index && iter == it.iter);
    }

    /*! Inequality operator. */
    bool operator!= (const iterator& it)
    {
      return (! (*this == it));
    }

    /*! Get the current object the iterator points on. */
    const_reference operator* () const
    {
      CGAL_precondition (p_hash != NULL &&
                         index >= 0 && index < static_cast<int>(p_hash->n_buckets));

      return (*iter);
    }

    /*! Get a pointer for the current object the iterator points on. */
    const_pointer operator-> () const
    {
      CGAL_precondition (p_hash != NULL &&
                         index >= 0 && index < static_cast<int>(p_hash->n_buckets));

      return (&(*iter));
    }

    /* Increment the iterator (prefix notation). */
    iterator& operator++ ()
    {
      // Check end conditions.
      CGAL_precondition (p_hash != NULL);

      if (index >= static_cast<int>(p_hash->n_buckets))
        return (*this);

      if (index < 0)
      {
        // Find the first non-empty bucket and get an iterator for it.
        index = p_hash->_next_non_empty_bucket (0);

        if (index < static_cast<int>(p_hash->n_buckets))
          iter = (p_hash->buckets[index]).begin();

        return (*this);
      }
      
      // Try to increment the current list iterator.
      ++iter;

      if (iter == (p_hash->buckets[index]).end())
      {
        // Find the next non-empty bucket and get an iterator for it.
        index = p_hash->_next_non_empty_bucket (index + 1);

        if (index < static_cast<int>(p_hash->n_buckets))
          iter = (p_hash->buckets[index]).begin();
      }
        
      return (*this);
    }

    /* Increment the iterator (postfix notation). */
    iterator operator++ (int )
    {
      iterator   temp = *this;

      ++(*this);
      return (temp);
    }

    /* Decrement the iterator (prefix notation). */
    iterator& operator-- ()
    {
      // Check end conditions.
      CGAL_precondition (p_hash != NULL);

      if (index >= static_cast<int>(p_hash->n_buckets))
      {
        // Find the last non-empty bucket and get an iterator for it.
        index = p_hash->_prev_non_empty_bucket 
                           (static_cast<int>(p_hash->n_buckets) - 1);

        if (index >= 0)
        {
          iter = (p_hash->buckets[index]).end();
          --iter;
        }

        return (*this);
      }

      if (index < 0)
        return (*this);
      
      // Try to decrement the current list iterator.
      if (iter != (p_hash->buckets[index]).begin())
      {
        --iter;
      }
      else
      {
        // Find the nprevious non-empty bucket and get an iterator for it.
        index = p_hash->_prev_non_empty_bucket (index - 1);

        if (index >= 0)
        {
          iter = (p_hash->buckets[index]).end();
          --iter;
        }
      }
        
      return (*this);
    }

    /* Decrement the iterator (postfix notation). */
    iterator operator-- (int )
    {
      iterator   temp = *this;

      --(*this);
      return (temp);
    }

  };

  friend class iterator;

  /// \name Scan the objects in the hash.
  //@{

  /*! Get an iterator for the first object in the hash. */
  iterator begin () const
  {
    return (iterator (*this, 0));
  }

  /*! Get a past-the end iterator for the objects. */
  iterator end () const
  {
    return (iterator (*this, n_buckets));
  }

  //@}

  /// \name Hash operations.
  //@{

  /*!
   * Check if an object with the given key exists in the map.
   * \param ket The key to look for.
   * \return An iterator for an objects having the given key,
   *         or a past-the-end iterator if no such object was found.
   */
  iterator find (const key_type& key) const
  {
    // Look for the object in the approriate bucket.
    const size_type           index = hash_func (key) % n_buckets;
    Bucket&                   my_bucket = const_cast<Bucket&>(buckets[index]);
    Bucket_iterator           bucket_it = _find_in_bucket (index, key);

    if (bucket_it == my_bucket.end())
      // The object was not found.
      return (end());

    // Create an iterator pointing to the object.
    return (iterator (*this, index, bucket_it));
  }

  /*!
   * Insert an object into the hash.
   * \param val The object to insert.
   * \return A pair comprised of an iterator for the inserted object,
   *         and a flag indicating whether the insertion tool place
   *         (false if the object had already existed in the hash).
   */
  std::pair<iterator,bool> insert (const value_type& val)
  {
    // Look for the object in the approriate bucket.
    const size_type           index = hash_func (val) % n_buckets;
    Bucket_iterator           bucket_it = _find_in_bucket (index, val);
    std::pair<iterator,bool>  res;

    // Check if the object already exists.
    if (bucket_it != buckets[index].end())
    {
      // In case the object already exists:
      res.first = iterator (*this, index, bucket_it);
      res.second = false;
    }
    else
    {
      // Insert a new object to the approriate bucket.
      buckets[index].push_front (val);

      n_objects++;

      res.first = iterator (*this, index, buckets[index].begin());
      res.second = true;
    }
      
    return (res);
  }

  /*!
   * Erase an object given an iterator pointing to it.
   * \param pos An iterator pointing at the object to be erased.
   */
  void erase (iterator pos)
  {
    // Erase the object from the approriate bucket.
    CGAL_precondition (pos.p_hash == this);

    const int        index = pos.index;
    Bucket_iterator  bucket_it = pos.iter;

    CGAL_assertion (index >= 0 && index < static_cast<int>(n_buckets));
    
    buckets[index].erase (bucket_it);
    n_objects--;

    return;
  }

  /*!
   * Erase an object with the given key.
   * \param key The key to delete.
   * \return Whether an object has been found and erased.
   */
  bool erase (const key_type& key)
  {
    // Look for the object in the approriate bucket.
    const size_type           index = hash_func (key) % n_buckets;
    Bucket_iterator           bucket_it = _find_in_bucket (index, key);

    if (bucket_it == buckets[index].end())
      // The object was not found.
      return (false);

    // Erase the object from the approriate bucket.
    buckets[index].erase (bucket_it);
    n_objects--;

    return (true);
  }

  /*!
   * Clear the contents of the hash.
   */
  void clear ()
  {
    // Clear all buckets.
    typename std::vector<Bucket>::iterator  it;

    for (it = buckets.begin(); it != buckets.end(); ++it)
      (*it).clear();

    n_objects = 0;
    return;
  }
  //@}

  /// \name Re-arranging the hash structure.
  //@{

  /*!
   * Resize the hash (takes time proporsional to the number of objects
   * stored in the hash).
   * \param n The new number of buckets.
   */
  void resize (size_type n)
  {
    rehash (n, hash_func);
  };

  /*! Rehash the structure.
   * \param n The new number of buckets.
   * \param hash_f The new hashing functor.
   */
  void rehash (size_type n, hasher hash_f)
  {
    // Set the hashing functor.
    hash_func = hash_f;

    // Rehash the structure.
    std::vector<Bucket>                     new_buckets (n);
    typename std::vector<Bucket>::iterator  it;
    Bucket_iterator                         bucket_it;
    size_type                               index;

    for (it = buckets.begin(); it != buckets.end(); ++it)
    {
      for (bucket_it = (*it).begin(); bucket_it != (*it).end(); ++bucket_it)
      {
        const value_type&  val = *bucket_it;
        
        index = hash_func (val) % n;
        new_buckets[index].push_back (val);
      }
    }

    // Set the new buckets.
    buckets = new_buckets;
    n_buckets = n;

    return;
  }
  //@}

private:

  /*! Get the index of the previous non-empty bucket. */
  int _prev_non_empty_bucket (int index) const
  {
    while (index >= 0 && buckets[index].empty())
      index--;

    return (index);
  }

  /*! Get the index of the next non-empty bucket. */
  int _next_non_empty_bucket (int index) const
  {
    while (index < static_cast<int>(n_buckets) && buckets[index].empty())
      index++;

    return (index);
  }

  /*! Find the given object in the given bucket. */
  Bucket_iterator _find_in_bucket (std::size_t index,
                                   const value_type& val) const
  {
    Bucket&             my_bucket = const_cast<Bucket&>(buckets[index]);
    Bucket_iterator     iter = my_bucket.begin();

    while (iter != my_bucket.end() && ! equal_func (val, *iter))
      ++iter;

    return (iter);
  }
};

/*! \class
 * An implementation of a hash-based map.
 */
template <class Key_, class Data_,
          class Hash_functor_,
          class EqualKey_ = Default_is_equal_functor<Key_> >
class Hash_map
{
public:

  // Type definitions (for STL compatibility):
  typedef Key_                   key_type;
  typedef Data_                  data_type;
  typedef std::pair<Key_,Data_>  value_type;
  typedef Hash_functor_          hasher;
  typedef EqualKey_              key_equal;
  typedef value_type*            pointer;
  typedef value_type&            reference;
  typedef const value_type&      const_reference;
  typedef size_t                 size_type;
  typedef size_t                 difference_type;

protected:

  /*! \class
   * A functor for the value_type pair.
   */
  class Value_type_hash_functor
  {
  private:

    Hash_functor_       hash_func;

  public:

    Value_type_hash_functor () :
      hash_func ()
    {}

    Value_type_hash_functor (Hash_functor_ func) :
      hash_func (func)
    {}

    size_type operator() (const value_type& vt) const
    {
      return (hash_func (vt.first));
    }
  };

  /*! \class
   * A functor for comparing value_type pairs.
   */
  class Value_type_equal
  {
  private:

    EqualKey_          equal_func;

  public:

    Value_type_equal () :
      equal_func ()
    {}

    Value_type_equal (EqualKey_ func) :
      equal_func (func)
    {}

    bool operator() (const value_type& vt1, const value_type& vt2) const
    {
      return (equal_func (vt1.first, vt2.first));
    }
  };

  typedef Open_hash<value_type, 
                    Value_type_hash_functor, Value_type_equal>    Hash;

  Value_type_hash_functor            vt_hash_func;
  Value_type_equal                   vt_equal_func;
  Hash                               hash;

public:

  /// \name Construction and destruction methods.
  //@{

  /*! Default constructor. */
  Hash_map () :
    vt_hash_func (hasher()),
    vt_equal_func (key_equal()),
    hash ()
  {}

  /*! Constructor with an initial number of buckets. */
  Hash_map (size_type n,
            hasher hash_f = hasher(),
            key_equal equal_f = key_equal()) :
    vt_hash_func (hash_f),
    vt_equal_func (equal_f),
    hash (n, vt_hash_func, vt_equal_func)
  {}

  /*! Destructor. */
  virtual ~Hash_map ()
  {}
  //@}

  /// \number Access the hash properties.
  //@{

  /*! Get the number of buckets in the hash map. */
  size_type bucket_count () const
  {
    return (hash.bucket_count());
  }

  /*! Get the number of objects stored in the hash map. */
  size_type size () const
  {
    return (hash.size());
  }

  /*! Check whether the hash map is empty. */
  bool empty () const
  {
    return (hash.empty());
  }
  //@}

  typedef typename Hash::iterator    iterator;

  /// \name Scan the objects in the hash.
  //@{

  /*! Get an iterator for the first object in the hash map. */
  iterator begin () const
  {
    return (hash.begin());
  }

  /*! Get a past-the end iterator for the objects map. */
  iterator end () const
  {
    return (hash.end());
  }

  //@}

  /// \name Hash-map operations.
  //@{

  /*!
   * Check if an object with the given key exists in the map.
   * \param ket The key to look for.
   * \return An iterator for an objects having the given key,
   *         or a past-the-end iterator if no such object was found.
   */
  iterator find (const key_type& key) const
  {
    value_type    vt (key, data_type());
    return (hash.find (vt));
  }

  /*!
   * Returns a reference to the object that is associated with the given key.
   * If the hash map does not contain such an object, it inserts an entry
   * corresponding to the given key.
   * \param key The key.
   * \return The data associated with this key.
   */
  data_type& operator[] (const key_type& key)
  {
    value_type                  vt (key, data_type());
    std::pair<iterator, bool>   res = hash.insert (vt);

    return (const_cast<data_type&> ((*(res.first)).second));
  }

  /*!
   * Insert an object into the hash.
   * \param val The object to insert.
   * \return A pair comprised of an iterator for the inserted object,
   *         and a flag indicating whether the insertion tool place
   *         (false if the object had already existed in the hash).
   */
  std::pair<iterator,bool> insert (const value_type& val)
  {
    return (hash.insert (val));
  }

  /*!
   * Erase an object given an iterator pointing to it.
   * \param pos An iterator pointing at the object to be erased.
   */
  void erase (iterator pos)
  {
    hash.erase (pos);
    return;
  }

  /*!
   * Erase an object with the given key.
   * \param key The key to delete.
   * \return Whether an object has been found and erased.
   */
  bool erase (const key_type& key)
  {
    value_type    vt (key, data_type());
    return (hash.erase (vt));
  }

  /*!
   * Clear the contents of the hash map.
   */
  void clear ()
  {
    hash.clear();
  }
  //@}

  /// \name Re-arranging the hash structure.
  //@{

  /*!
   * Resize the hash map (takes time proporsional to the number of objects
   * stored in the hash).
   * \param n The new number of buckets.
   */
  void resize (size_type n)
  {
    hash.resize (n);
  };

  /*! Rehash the structure.
   * \param n The new number of buckets.
   * \param hash_f The new hashing functor.
   */
  void rehash (size_type n, hasher hash_f)
  {
    vt_hash_func = Value_type_hash_functor (hash_f);
    hash.rehash (n, vt_hash_func);
  }
  //@}

};

} //namespace CGAL

#endif
