// Copyright (c) 2005  Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_MULTISET_H
#define CGAL_MULTISET_H

#include <CGAL/config.h>
#include <CGAL/assertions.h>
#include <CGAL/multiset_assertions.h>
#include <CGAL/enum.h>
#include <CGAL/memory.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/Compact_container.h>
#include <iterator>

namespace CGAL {

/*!
 * Container class representing a red-black tree, which is a balanced binary
 * tree that has the following invariants:
 * 1. Each node has a color, which is either red or black.
 * 2. Each red node has two red children (if a child is missing, it is
 *    considered as a black node).
 * 3. The number of black nodes from every path from the tree root to a leaf
 *    is the same for all tree leaves (it is called the 'black height' of the
 *    tree).
 * Due to propeties 2-3, the height of a red-black tree containing n nodes
 * is bounded by 2*log_2(n).
 *
 * The Multiset template requires three template parmeters:
 * - The contained Type class represents the objects stored in the tree.
 *   It has to support the default constructor, the copy constructor and
 *   the assignment operator (operator=).
 * - Compare is a three-valued functor used to define the order of objects of
 *   class Type: It has to support an operator() that receives two objects from
 *   the Type class and returns SMALLER, EQUAL or LARGER, depending on the
 *   comparison result.
 *   In case the deafult parameter is supplied, the Type class has to support
 *   the less-than (<) and the equal (==) operators.
 * - The Allocator represents an allocator class. By default, it is the CGAL
 *   allocator.
 */

template <typename Type_,
          class Compare_ = CGAL::Compare<Type_>,
          class Allocator_ = CGAL_ALLOCATOR(int),
#ifdef CGAL_MULTISET_USE_COMPACT_CONTAINER_AS_DEFAULT
          typename UseCompactContainer = Tag_true>
#else
          typename UseCompactContainer = Tag_false>
#endif
class Multiset
{
public:

  // General type definitions:
  typedef Type_                                     Type;
  typedef Compare_                                  Compare;
  typedef Allocator_                                Allocator;
  typedef Multiset<Type, Compare, Allocator, UseCompactContainer>
                                                    Self;

  // Type definitions for STL compatibility.
  typedef Type                                      value_type;
  typedef Type                                      key_type;
  typedef value_type&                               reference;
  typedef const value_type&                         const_reference;
  typedef Compare                                   key_compare;
  typedef Compare                                   value_compare;
  typedef size_t                                    size_type;
  typedef std::ptrdiff_t                  difference_type;

protected:

  /*! \struct
   * Representation of a node in a red-black tree.
   */
  struct Node
  {

    enum Color
    {
      RED,                        // Regular red node in the tree.
      BLACK,                      // Regular black node in the tree.
      DUMMY_BEGIN,                // The dummy before-the-begin tree node.
      DUMMY_END                   // The dummy past-the-end tree node.
    };

    typedef char Node_color;

    Type         object;           // The stored object.
    Node_color   color;            // The color of the node.
    Node        *parentP;          // Points on the parent node.
    Node        *rightP;           // Points on the right child of the node.
    Node        *leftP;            // Points on the left child of the node.

    /*! Default constructor. */
    Node() :
      parentP(nullptr),
      rightP(nullptr),
      leftP(nullptr)
    {}

    /*!
     * Constructor of a red-black tree node.
     * \param _object The object stored in the node.
     * \param _color The color of the node.
     */
    Node (const Type& _object,  Node_color _color) :
      object(_object),
      color(_color),
      parentP(nullptr),
      rightP(nullptr),
      leftP(nullptr)
    {}

    /*!
     * Initialize the node.
     * \param _object The object stored in the node.
     * \param _color The color of the node.
     */
    void init (const Type& _object, Node_color _color)
    {
      object = _object;
      color = _color;
    }

    /*! Destructor. */
    ~Node ()
    {}

    /*!
     * Check if the node is valid (not a dummy node).
     */
    inline bool is_valid () const
    {
      return (color == BLACK || color == RED);
    }

    /*!
     * Check if the node is red.
     */
    inline bool is_red () const
    {
      return (color == RED);
    }

    /*!
     * Check if the node is black.
     * Note that dummy nodes are also considered to be black.
     */
    inline bool is_black () const
    {
      return (color != RED);
    }

    /*!
     * Get the previous node in the tree (according to the tree order).
     */
    Node* predecessor () const
    {
      // The DUMMY_BEGIN node has no predecessor.
      CGAL_multiset_assertion (color != DUMMY_BEGIN);

      Node        *predP;

      if (leftP != nullptr)
      {
        // If there is a left child, the predecessor is the maximal object in
        // the sub-tree spanned by this child.
        predP = leftP;
        while (predP->rightP != nullptr)
          predP = predP->rightP;
      }
      else
      {
        // Otherwise, go up the tree until reaching the parent from the right
        // direction (this also covers the case of the DUMMY_END node).
        const Node    *prevP = this;

        predP = parentP;
        while (predP != nullptr && prevP == predP->leftP)
        {
          prevP = predP;
          predP = predP->parentP;
        }
      }

      return (predP);
    }

    /*!
     * Get the next node in the tree (according to the tree order).
     */
    Node* successor () const
    {
      // The DUMMY_END node has no successor.
      CGAL_multiset_assertion (color != DUMMY_END);

      Node        *succP;

      if (rightP != nullptr)
      {
        // If there is a right child, the successor is the minimal object in
        // the sub-tree spanned by this child.
        succP = rightP;
        while (succP->leftP != nullptr)
          succP = succP->leftP;
      }
      else
      {
        // Otherwise, go up the tree until reaching the parent from the left
        // direction (this also covers the case of the DUMMY_BEGIN node).
        const Node    *prevP = this;

        succP = parentP;
        while (succP != nullptr && prevP == succP->rightP)
        {
          prevP = succP;
          succP = succP->parentP;
        }
      }

      return (succP);
    }

    void* for_compact_container() const
    {
      return parentP;
    }

    void for_compact_container (void * p)
    {
      reinterpret_cast<void*&>(parentP) = p;
    }

  };

  // Rebind the allocator to the Node type:
  class Default_node_allocator
  {
    typedef std::allocator_traits<Allocator> Allocator_traits;
    typedef typename Allocator_traits::template rebind_alloc<Node> Base;
    Base base;

  public:

    Node* allocate (const Node& n)
    {
      Node* new_node = base.allocate(1);
      std::allocator_traits<Base>::construct(base, new_node, n);
      return new_node;
    }

    void deallocate (Node* n)
    {
      std::allocator_traits<Base>::destroy(base, n);
      base.deallocate (n, 1);
    }
  };

  class CC_node_allocator
  {
    typedef Compact_container<Node> Base;
    Base base;

  public:

    Node* allocate (const Node& n)
    {
      Node* new_node = &*base.emplace(n);
      return new_node;
    }

    void deallocate (Node* n)
    {
      base.erase (base.iterator_to(*n));
    }
  };

  typedef typename std::conditional<UseCompactContainer::value,
                                    CC_node_allocator,
                                    Default_node_allocator>::type Node_allocator;

public:

  // Forward decleration:
  class const_iterator;

  /*! \class
   * An iterator for the red black tree.
   * The iterator also servers as a handle to a tree nodes.
   */
  class iterator
  {
    // Give the red-black tree class template access to the iterator's members.
    friend class Multiset<Type, Compare, Allocator, UseCompactContainer>;
    friend class const_iterator;

  public:

    // Type definitions:
    typedef std::bidirectional_iterator_tag  iterator_category;
    typedef Type                             value_type;
    typedef std::ptrdiff_t         difference_type;
    typedef size_t                           size_type;
    typedef value_type&                      reference;
    typedef value_type*                      pointer;

  private:

    Node    *nodeP;              // Points to a tree node.

    /*! Private constructor. */
    iterator (Node *_nodeP) :
      nodeP (_nodeP)
    {}

  public:

    /*! Deafult constructor. */
    iterator () :
      nodeP (nullptr)
    {}

    /*! Equality operator. */
    bool operator== (const iterator& iter) const
    {
      return (nodeP == iter.nodeP);
    }

    /*! Inequality operator. */
    bool operator!= (const iterator& iter) const
    {
      return (nodeP != iter.nodeP);
    }

    /*! Increment operator (prefix notation). */
    iterator& operator++ ()
    {
      CGAL_multiset_precondition (nodeP != nullptr);

      nodeP = nodeP->successor();
      return (*this);
    }

    /*! Increment operator (postfix notation). */
    iterator operator++ (int )
    {
      CGAL_multiset_precondition (nodeP != nullptr);

      iterator temp = *this;

      nodeP = nodeP->successor();
      return (temp);
    }

    /*! Decrement operator (prefix notation). */
    iterator& operator-- ()
    {
      CGAL_multiset_precondition (nodeP != nullptr);

      nodeP = nodeP->predecessor();
      return (*this);
    }

    /*! Decrement operator (postfix notation). */
    iterator operator-- (int )
    {
      CGAL_multiset_precondition (nodeP != nullptr);

      iterator temp = *this;

      nodeP = nodeP->predecessor();
      return (temp);
    }

    /*!
     * Get the object pointed by the iterator.
     */
    reference operator* () const
    {
      CGAL_multiset_precondition (nodeP != nullptr && nodeP->is_valid());

      return (nodeP->object);
    }

    /*!
     * Get the a pointer object pointed by the iterator.
     */
    pointer operator-> () const
    {
      CGAL_multiset_precondition (nodeP != nullptr && nodeP->is_valid());

      return (&(nodeP->object));
    }

  };
  friend class iterator;

  /*! \class
   * An const iterator for the red black tree.
   * The iterator also servers as a const handle to a tree nodes.
   */
  class const_iterator
  {
    // Give the red-black tree class template access to the iterator's members.
    friend class Multiset<Type, Compare, Allocator, UseCompactContainer>;

  public:

    // Type definitions:
    typedef std::bidirectional_iterator_tag  iterator_category;
    typedef Type                             value_type;
    typedef std::ptrdiff_t         difference_type;
    typedef size_t                           size_type;
    typedef const value_type&                reference;
    typedef const value_type*                pointer;

  private:

    const Node    *nodeP;     // Points to a tree node.

    /*! Private constructor. */
    const_iterator (const Node *_nodeP) :
      nodeP (_nodeP)
    {}

  public:

    /*! Deafult constructor. */
    const_iterator () :
      nodeP (nullptr)
    {}

    /*! Constructor from a mutable iterator. */
    const_iterator (const iterator& iter) :
      nodeP (iter.nodeP)
    {}

    /*! Equality operator. */
    bool operator== (const const_iterator& iter) const
    {
      return (nodeP == iter.nodeP);
    }

    /*! Inequality operator. */
    bool operator!= (const const_iterator& iter) const
    {
      return (nodeP != iter.nodeP);
    }

    /*! Increment operator (prefix notation). */
    const_iterator& operator++ ()
    {
      CGAL_multiset_precondition (nodeP != nullptr);

      nodeP = nodeP->successor();
      return (*this);
    }

    /*! Increment operator (postfix notation). */
    const_iterator operator++ (int )
    {
      CGAL_multiset_precondition (nodeP != nullptr);

      const_iterator temp = *this;

      nodeP = nodeP->successor();
      return (temp);
    }

    /*! Decrement operator (prefix notation). */
    const_iterator& operator-- ()
    {
      CGAL_multiset_precondition (nodeP != nullptr);

      nodeP = nodeP->predecessor();
      return (*this);
    }

    /*! Decrement operator (postfix notation). */
    const_iterator operator-- (int )
    {
      CGAL_multiset_precondition (nodeP != nullptr);

      const_iterator temp = *this;

      nodeP = nodeP->predecessor();
      return (temp);
    }

    /*!
     * Get the object pointed by the iterator.
     */
    reference operator* () const
    {
      CGAL_multiset_precondition (nodeP != nullptr && nodeP->is_valid());

      return (nodeP->object);
    }

    /*!
     * Get the a pointer object pointed by the iterator.
     */
    pointer operator-> () const
    {
      CGAL_multiset_precondition (nodeP != nullptr && nodeP->is_valid());

      return (&(nodeP->object));
    }

  };
  friend class const_iterator;

  // Define the reverse iterators:
  typedef std::reverse_iterator<iterator>         reverse_iterator;
  typedef std::reverse_iterator<const_iterator>   const_reverse_iterator;

protected:

  // Data members:
  Node            *rootP;        // Pointer to the tree root.
  size_t           iSize;        // Number of objects stored in the tree.
  size_t           iBlackHeight; // The black-height of the tree (the number
                                 // of black nodes from the root to each leaf).
  Compare          comp_f;       // A comparison functor.
  Node_allocator   node_alloc;   // Allocator for the tree nodes.
  Node             beginNode;    // A fictitious before-the-minimum node.
                                 // Its parent is the leftmost tree node.
  Node             endNode;      // A fictitious past-the-maximum node.
                                 // Its parent is the rightmost tree node.

public:

  /// \name Construction and destruction functions.
  //@{

  /*!
   * Default constructor. [takes O(1) operations]
   */
  Multiset ();

  /*!
   * Constructor with a comparison object. [takes O(1) operations]
   * \param comp A comparison object to be used by the tree.
   */
  Multiset (const Compare& comp);

  /*!
   * Copy constructor. [takes O(n) operations]
   * \param tree The copied tree.
   */
  Multiset (const Self& tree);

  /*!
   * Construct a tree that contains all objects in the given range.
   * [takes O(n log n) operations]
   * \param first An iterator for the first object in the range.
   * \param last A past-the-end iterator for the range.
   */
  template <class InputIterator>
  Multiset (InputIterator first, InputIterator last,
            const Compare& comp = Compare()) :
    rootP (nullptr),
    iSize (0),
    iBlackHeight (0),
    comp_f (comp)
  {
    // Mark the two fictitious nodes as dummies.
    beginNode.color = Node::DUMMY_BEGIN;
    endNode.color = Node::DUMMY_END;

    // Insert all objects to the tree.
    while (first != last)
    {
      insert (*first);
      ++first;
    }
  }

  /*!
   * Destructor. [takes O(n) operations]
   */
  virtual ~Multiset ();

  /*!
   * Assignment operator. [takes O(n) operations]
   * \param tree The copied tree.
   */
  Self& operator= (const Self& tree);

  /*!
   * Swap two trees. [takes O(1) operations]
   * \param tree The copied tree.
   */
  void swap (Self& tree);
  //@}

  /// \name Comparison operations.
  //@{

  /*!
   * Test two trees for equality. [takes O(n) operations]
   * \param tree The compared tree.
   */
  bool operator== (const Self& tree) const;

  /*!
   * Check if our tree is lexicographically smaller. [takes O(n) operations]
   * \param tree The compared tree.
   */
  bool operator< (const Self& tree) const;
  //@}

  /// \name Access functions.
  //@{

  /*!
   * Get the comparsion object used by the tree (non-const version).
   */
  inline Compare& key_comp ()
  {
    return (comp_f);
  }

  /*!
   * Get the comparsion object used by the tree (non-const version).
   */
  inline Compare& value_comp ()
  {
    return (comp_f);
  }


  /*!
   * Get the comparsion object used by the tree (const version).
   */
  inline const Compare& key_comp () const
  {
    return (comp_f);
  }

  /*!
   * Get the comparsion object used by the tree (const version).
   */
  inline const Compare& value_comp () const
  {
    return (comp_f);
  }

  /*!
   * Get an iterator for the minimum object in the tree (non-const version).
   */
  inline iterator begin ();

  /*!
   * Get a past-the-end iterator for the tree objects (non-const version).
   */
  inline iterator end ();

  /*!
   * Get an iterator for the minimum object in the tree (const version).
   */
  inline const_iterator begin () const;

  /*!
   * Get a past-the-end iterator for the tree objects (const version).
   */
  inline const_iterator end () const;

  /*!
   * Get a reverse iterator for the maxnimum object in the tree
   * (non-const version).
   */
  inline reverse_iterator rbegin ();

  /*!
   * Get a pre-the-begin reverse iterator for the tree objects
   * (non-const version).
   */
  inline reverse_iterator rend ();

  /*!
   * Get a reverse iterator for the maximum object in the tree (const version).
   */
  inline const_reverse_iterator rbegin () const;

  /*!
   * Get a pre-the-begin reverse iterator for the tree objects (const version).
   */
  inline const_reverse_iterator rend () const;

  /*!
   * Check whether the tree is empty.
   */
  inline bool empty () const
  {
    return (rootP == nullptr);
  }

  /*!
   * Get the size of the tree. [takes O(1) operations, unless the tree
   * was involved in a split operation, then it may take O(n) time.]
   * \return The number of objects stored in the tree.
   */
  size_t size () const;

  /*!
   * Get the maximal possible size (equivalent to size()).
   */
  size_t max_size () const
  {
    return (size());
  }
  //@}

  /// \name Insertion functions.

  /*!
   * Insert an object into the tree. [takes O(log n) operations]
   * \param object The object to be inserted.
   * \return An iterator pointing to the inserted object.
   */
  iterator insert (const Type& object);

  /*!
   * Insert a range of k objects into the tree. [takes O(k log n) operations]
   * \param first An iterator for the first object in the range.
   * \param last A past-the-end iterator for the range.
   */
  template <class InputIterator>
  void insert (InputIterator first, InputIterator last)
  {
    // Insert all objects to the tree.
    while (first != last)
    {
      insert (*first);
      ++first;
    }

    return;
  }

  /*!
   * Insert an object to the tree, with a given hint to its position.
   * [takes O(log n) operations at worst-case, but only O(1) amortized]
   * \param position A hint for the position of the object.
   * \param object The object to be inserted.
   * \return An iterator pointing to the inserted object.
   */
  iterator insert (iterator position,
                   const Type& object);

  /*!
   * Insert an object to the tree, as the successor the given object.
   * [takes O(log n) operations at worst-case, but only O(1) amortized]
   * \param position Points to the object after which the new object should
   *                 be inserted (or an invalid iterator to insert the object
   *                 as the tree minimum).
   * \param object The object to be inserted.
   * \pre The operation does not violate the tree properties.
   * \return An iterator pointing to the inserted object.
   */
  iterator insert_after (iterator position,
                         const Type& object);

  /*!
   * Insert an object to the tree, as the predecessor the given object.
   * [takes O(log n) operations at worst-case, but only O(1) amortized]
   * \param position Points to the object before which the new object should
   *                 be inserted (or an invalid iterator to insert the object
   *                 as the tree maximum).
   * \param object The object to be inserted.
   * \pre The operation does not violate the tree properties.
   * \return An iterator pointing to the inserted object.
   */
  iterator insert_before (iterator position,
                          const Type& object);

  /// \name Erasing functions.
  //@{

  /*!
   * Erase objects from the tree. [takes O(log n) operations]
   * \param object The object to be removed.
   * \return The number of objects removed from the tree.
   *          Note that all iterators to the erased objects become invalid.
   */
  size_t erase (const Type& object);

  /*!
   * Remove the object pointed by the given iterator.
   * [takes O(log n) operations at worst-case, but only O(1) amortized]
   * \param position An iterator pointing the object to be erased.
   * \pre The iterator must be a valid.
   *      Note that all iterators to the erased object become invalid.
   */
  void erase (iterator position);

  /*!
   * Clear the contents of the tree. [takes O(n) operations]
   */
  void clear ();

  //@}

  /// \name Search functions.
  //@{

  /*!
   * Search the tree for the given key (non-const version).
   * [takes O(log n) operations]
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \return A iterator pointing to the first equivalent object in the tree,
   *         or end() if no such object is found in the tree.
   */
  iterator find (const Type& object)
  {
    return (find (object, comp_f));
  }

  template <class Key, class CompareKey>
  iterator find (const Key& key,
                 const CompareKey& comp_key)
  {
    bool    is_equal;
    Node   *nodeP = _bound (LOWER_BOUND, key, comp_key, is_equal);

    if (_is_valid(nodeP) && is_equal)
      return (iterator (nodeP));
    else
      return (iterator (&endNode));
  }

  /*!
   * Search the tree for the given key (const version).
   * [takes O(log n) operations]
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \return A iterator pointing to the first equivalent object in the tree,
   *         or end() if no such object is found in the tree.
   */
  const_iterator find (const Type& object) const
  {
    return (find (object, comp_f));
  }

  template <class Key, class CompareKey>
  const_iterator find (const Key& key,
                       const CompareKey& comp_key) const
  {
    bool          is_equal;
    const Node   *nodeP = _bound (LOWER_BOUND, key, comp_key, is_equal);

    if (_is_valid (nodeP) && is_equal)
      return (const_iterator (nodeP));
    else
      return (const_iterator (&endNode));
  }

  /*!
   * Count the number of object in the tree equivalent to a given key.
   * [takes O(log n + d) operations]
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \return The number of equivalent objects.
   */
  size_type count (const Type& object) const
  {
    return (count (object, comp_f));
  }

  template <class Key, class CompareKey>
  size_type count (const Key& key,
                   const CompareKey& comp_key) const
  {
    // Get the first equivalent object (if any).
    size_t      n_objects = 0;
    bool        is_equal;
    const Node *nodeP = _bound (LOWER_BOUND, key, comp_key, is_equal);

    if (! is_equal)
      return (0);

    while (_is_valid (nodeP) &&
           comp_key (key, nodeP->object) == EQUAL)
    {
      n_objects++;

      // Proceed to the successor.
      nodeP = nodeP->successor();
    }

    return (n_objects);
  }

  /*!
   * Get the first element whose key is not less than a given key
   * (non-const version). [takes O(log n) operations]
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \return The lower bound of the key, or end() if the key is not found
   *         in the tree.
   */
  iterator lower_bound (const Type& object)
  {
    return (lower_bound (object, comp_f));
  }

  template <class Key, class CompareKey>
  iterator lower_bound (const Key& key,
                        const CompareKey& comp_key)
  {
    bool    is_equal;
    Node   *nodeP = _bound (LOWER_BOUND, key, comp_key, is_equal);

    if (_is_valid (nodeP))
      return (iterator (nodeP));
    else
      return (iterator (&endNode));
  }

  /*!
   * Get the first element whose key is not less than a given key
   * (non-const version). [takes O(log n) operations]
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \return The lower bound of the key, along with a flag indicating whether
   *         the bound equals the given key.
   */
  std::pair<iterator, bool> find_lower (const Type& object)
  {
    return (find_lower (object, comp_f));
  }

  template <class Key, class CompareKey>
  std::pair<iterator, bool> find_lower (const Key& key,
                                        const CompareKey& comp_key)
  {
    bool    is_equal;
    Node   *nodeP = _bound (LOWER_BOUND, key, comp_key, is_equal);

    if (_is_valid (nodeP))
      return (std::make_pair (iterator (nodeP), is_equal));
    else
      return (std::make_pair (iterator (&endNode), false));
  }

  /*!
   * Get the first element whose key is greater than a given key
   * (non-const version). [takes O(log n) operations]
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \return The upper bound of the key, or end() if the key is not found
   *         in the tree.
   */
  iterator upper_bound (const Type& object)
  {
    return (upper_bound (object, comp_f));
  }

  template <class Key, class CompareKey>
  iterator upper_bound (const Key& key,
                        const CompareKey& comp_key)
  {
    bool    is_equal;
    Node   *nodeP = _bound (UPPER_BOUND, key, comp_key, is_equal);

    if (_is_valid (nodeP))
      return (iterator (nodeP));
    else
      return (iterator (&endNode));
  }

  /*!
   * Get the first element whose key is not less than a given key
   * (const version). [takes O(log n) operations]
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \return The lower bound of the key, or end() if the key is not found
   *         in the tree.
   */
  const_iterator lower_bound (const Type& object) const
  {
    return (lower_bound (object, comp_f));
  }

  template <class Key, class CompareKey>
  const_iterator lower_bound (const Key& key,
                              const CompareKey& comp_key) const
  {
    bool          is_equal;
    const Node   *nodeP = _bound (LOWER_BOUND, key, comp_key, is_equal);

    if (_is_valid (nodeP))
      return (const_iterator (nodeP));
    else
      return (const_iterator (&endNode));
  }

  /*!
   * Get the first element whose key is not less than a given key
   * (const version). [takes O(log n) operations]
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \return The lower bound of the key, along with a flag indicating whether
   *         the bound equals the given key.
   */
  std::pair<const_iterator, bool> find_lower (const Type& object) const
  {
    return (find_lower (object, comp_f));
  }

  template <class Key, class CompareKey>
  std::pair<const_iterator, bool> find_lower (const Key& key,
                                              const CompareKey& comp_key) const
  {
    bool          is_equal;
    const Node   *nodeP = _bound (LOWER_BOUND, key, comp_key, is_equal);

    if (_is_valid (nodeP))
      return (std::make_pair (const_iterator (nodeP), is_equal));
    else
      return (std::make_pair (const_iterator (&endNode), false));
  }

  /*!
   * Get the first element whose key is greater than a given key
   * (const version). [takes O(log n) operations]
   * \param object The query object.
   * \return The upper bound of the key, or end() if the key is not found
   *         in the tree.
   */
  const_iterator upper_bound (const Type& object) const
  {
    return (upper_bound (object, comp_f));
  }

  template <class Key, class CompareKey>
  const_iterator upper_bound (const Key& key,
                              const CompareKey& comp_key) const
  {
    bool          is_equal;
    const Node   *nodeP = _bound (UPPER_BOUND, key, comp_key, is_equal);

    if (_is_valid (nodeP))
      return (const_iterator (nodeP));
    else
      return (const_iterator (&endNode));
  }

  /*!
   * Get the range of objects in the tree that are equivalent to a given key
   * (non-const version). [takes O(log n + d) operations]
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \return A pair of (lower_bound(key), upper_bound(key)).
   */
  std::pair<iterator, iterator> equal_range (const Type& object)
  {
    return (equal_range (object, comp_f));
  }

  template <class Key, class CompareKey>
  std::pair<iterator, iterator>
  equal_range (const Key& key,
               const CompareKey& comp_key)
  {
    // Get the first equivalent object (if any).
    const size_t   max_steps = static_cast<size_t> (1.5 * iBlackHeight);
    bool           is_equal;
    Node          *lowerP = _bound (LOWER_BOUND, key, comp_key, is_equal);
    Node          *upperP = lowerP;
    size_t         n_objects = 0;

    if (is_equal)
    {
      while (_is_valid (upperP) &&
             comp_key (key, upperP->object) == EQUAL)
      {
        n_objects++;
        if (n_objects >= max_steps)
        {
          // If we have more than log(n) objects in the range, locate the
          // upper bound directly.
          upperP = _bound (UPPER_BOUND, key, comp_key, is_equal);
          break;
        }

        // Proceed to the successor.
        upperP = upperP->successor();
      }
    }

    return (std::pair<iterator, iterator>
      ((_is_valid (lowerP)) ? iterator (lowerP) : iterator (&endNode),
       (_is_valid (upperP)) ? iterator (upperP) : iterator (&endNode)));
  }

  /*!
   * Get the range of objects in the tree that are equivalent to a given key
   * (const version). [takes O(log n + d) operations]
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \return A pair of (lower_bound(key), upper_bound(key)).
   */
  std::pair<const_iterator, const_iterator>
  equal_range (const Type& object) const
  {
    return (equal_range (object, comp_f));
  }

  template <class Key, class CompareKey>
  std::pair<const_iterator, const_iterator>
  equal_range (const Key& key,
               const CompareKey& comp_key) const
  {
    // Get the first equivalent object (if any).
    const size_t   max_steps = static_cast<size_t> (1.5 * iBlackHeight);
    bool           is_equal;
    const Node    *lowerP = _bound (LOWER_BOUND, key, comp_key, is_equal);
    const Node    *upperP = lowerP;
    size_t         n_objects = 0;

    if (is_equal)
    {
      while (_is_valid (upperP) &&
             comp_key (key, upperP->object) == EQUAL)
      {
        n_objects++;
        if (n_objects >= max_steps)
        {
          // If we have more than log(n) objects in the range, locate the
          // upper bound directly.
          upperP = _bound (UPPER_BOUND, key, comp_key, is_equal);
          break;
        }

        // Proceed to the successor.
        upperP = upperP->successor();
      }
    }

    return (std::pair<const_iterator, const_iterator>
      ((_is_valid (lowerP)) ? const_iterator(lowerP) :
                              const_iterator(&endNode),
       (_is_valid (upperP)) ? const_iterator(upperP) :
                              const_iterator(&endNode)));
  }
  //@}

  /// \name Special functions.
  //@{

  /*!
   * Replace the object pointed by a given iterator with another object.
   * [takes O(1) operations]
   * \param position An iterator pointing the object to be replaced.
   * \param object The new object.
   * \pre The given iterator is valid.
   *      The operation does not violate the tree properties.
   */
  void replace (iterator position,
                const Type& object);

  /*!
   * Swap the location two objects in the tree, given by their positions.
   * [takes O(1) operations]
   * \param pos1 An iterator pointing to the first object.
   * \param pos1 An iterator pointing to the second object.
   * \pre The two iterators are valid.
   *      The operation does not violate the tree properties.
   */
  void swap (iterator pos1, iterator pos2);

  /*!
   * Catenate the tree with a given tree, whose minimal object is not less
   * than the maximal object of this tree. [takes O(log n) operations]
   * The function clears the other given tree, but all its iterators remain
   * valid and can be used with the catenated tree.
   * \param tree The tree to catenate to out tree.
   * \pre The minimal object in the given tree is not less than the maximal
   *      objects.
   */
  void catenate (Self& tree);

  /*!
   * Split the tree such that all remaining objects are less than a given
   * key, and all objects greater then (or equal to) this key form
   * a new output tree. [takes O(log n) operations]
   * \param key The split key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \param tree Output: The tree that will eventually contain all objects
   *                     greater than the split object.
   * \pre The output tree is initially empty.
   */
  void split (const Type& object, Self& tree)
  {
    split (object, comp_f, tree);
    return;
  }

  template <class Key, class CompareKey>
  void split (const Key& key, const CompareKey& comp_key,
              Self& tree)
  {
    split (lower_bound (key, comp_key), tree);
    return;
  }

  /*!
   * Split the tree at a given position, such that it contains all objects
   * in the range [begin, position) and all objects in the range
   * [position, end) form a new output tree. [takes O(log n) operations]
   * \param position An iterator pointing at the split position.
   * \param tree Output: The output tree.
   * \pre The output tree is initially empty.
   */
  void split (iterator position, Self& tree);

  //@}

  /// \name Debugging utilities.
  //@{

  /*!
   * Check validation of the tree. [takes o(n) operations]
   * \return true iff the tree is valid.
   *
   */
  bool is_valid() const;

  /*!
   * Get the height of the tree. [takes O(n) operations]
   * \return The length of the longest path from the root to a leaf node.
   */
  size_t height () const;

  /*!
   * Get the black-height of the tree. [takes O(1) operations]
   * \return The number of black nodes from the root to each leaf node.
   */
  inline size_t black_height () const
  {
    return (iBlackHeight);
  }
  //@}

protected:

  /// \name Auxiliary predicates.
  //@{

  /*! Check whether a node is valid. */
  inline bool _is_valid (const Node *nodeP) const
  {
    return (nodeP != nullptr && nodeP->is_valid());
  }

  /*! Check whether a node is red. */
  inline bool _is_red (const Node *nodeP) const
  {
    return (nodeP != nullptr && nodeP->color == Node::RED);
  }

  /*! Check whether a node is black. */
  inline bool _is_black (const Node *nodeP) const
  {
    // Note that invalid nodes are considered ro be black as well.
    return (nodeP == nullptr || nodeP->color != Node::RED);
  }
  //@}

  /// \name Auxiliary operations.
  //@{

  /*!
   * Move the contents of one tree to another without actually duplicating
   * the nodes. This operation also clears the copied tree.
   * \param tree The tree to be copied (and eventuallu cleared).
   */
  void _shallow_assign (Self& tree);

  /*!
   * Clear the properties of the tree, without actually deallocating its
   * nodes.
   */
  void _shallow_clear ();

  /*!
   * Return the pointer to the node containing the first object that is not
   * less than (or greater than) the given key.
   * \param type The bound type (lower or upper).
   * \param key The query key.
   * \param comp_key A comparison functor for comparing keys and objects.
   * \param is_equal Output: In case of a lower bound, indicates whether the
   *                         returned node contains an object equivalent to
   *                         the given key.
   * \return A node that contains the first object that is not less than the
   *         given key (if type is LOWER_BOUND), or anode that contains the
   *         first object that is greater than the given key (if type is
   *         UPPER_BOUND) - or a nullptr node if no such nodes exist.
   */
  enum Bound_type {LOWER_BOUND, UPPER_BOUND};

  template <class Key, class CompareKey>
  Node* _bound (Bound_type type,
                const Key& key,
                const CompareKey& comp_key,
                bool& is_equal) const
  {
    // Initially mark that the key is not found in the tree.
    is_equal = false;

    if (rootP == nullptr)
      // The tree is empty:
      return (nullptr);

    Node                      *currentP = rootP;
    Node                      *prevP = currentP;
    Comparison_result          comp_res = EQUAL;

    while (_is_valid (currentP))
    {
      comp_res = comp_key (key, currentP->object);

      if (comp_res == EQUAL)
      {
        // The key is found in the tree:
        if (type == LOWER_BOUND)
        {
          is_equal = true;

          // Lower bound computation:
          // Go backward as long as we find equivalent objects.
          prevP = currentP->predecessor();
          while (_is_valid (prevP))
          {
            if (comp_key (key, prevP->object) != EQUAL)
              break;
            currentP = prevP;
            prevP = currentP->predecessor();
          }
        }
        else
        {
          // Upper bound computation:
          // Go backward until we encounter a non-equivalent objects.
          do
          {
            currentP = currentP->successor();
          } while (_is_valid (currentP) &&
                   comp_key (key, currentP->object) == EQUAL);
        }

        return (currentP);
      }
      else if (comp_res == SMALLER)
      {
        prevP = currentP;

        // Go down to the left child.
        currentP = currentP->leftP;
      }
      else // comp_res == LARGER
      {
        prevP = currentP;

        // Go down to the right child.
        currentP = currentP->rightP;
      }
    }

    // If we reached here, the object is not found in the tree. We check if the
    // last node we visited (given by prevP) contains an object greater than
    // the query object. If not, we return its successor.
    if (comp_res == SMALLER)
      return (prevP);
    else
      return (prevP->successor());
  }

  /*!
   * Remove the object stored in the given tree node.
   * \param nodeP The node storing the object to be removed from the tree.
   */
  void _remove_at (Node* nodeP);

  /*!
   * Swap the location two nodes in the tree.
   * \param node1_P The first node.
   * \param node2_P The second node.
   */
  void _swap (Node* node1_P, Node* node2_P);

  /*!
   * Swap the location two sibling nodes in the tree.
   * \param node1_P The first node.
   * \param node2_P The second node.
   * \pre The two nodes have a common parent.
   */
  void _swap_siblings (Node* node1_P, Node* node2_P);

  /*!
   * Calculate the height of the sub-tree spanned by the given node.
   * \param nodeP The sub-tree root.
   * \return The height of the sub-tree.
   */
  size_t _sub_height (const Node* nodeP) const;

  /*!
   * Check whether a sub-tree is valid.
   * \param nodeP The sub-tree root.
   * \param sub_size Output: The size of the sub-tree.
   * \param sub_bh Output: The black height of the sub-tree.
   */
  bool _sub_is_valid (const Node* nodeP,
                      size_t& sub_size,
                      size_t& sub_bh) const;

  /*!
   * Get the leftmost node in the sub-tree spanned by the given node.
   * \param nodeP The sub-tree root.
   * \return The sub-tree minimum.
   */
  Node* _sub_minimum (Node* nodeP) const;

  /*!
   * Get the rightmost node in the sub-tree spanned by the given node.
   * \param nodeP The sub-tree root.
   * \return The sub-tree maximum.
   */
  Node* _sub_maximum (Node* nodeP) const;

  /*!
   * Left-rotate the sub-tree spanned by the given node.
   * \param nodeP The sub-tree root.
   */
  void _rotate_left (Node* nodeP);

  /*!
   * Right-rotate the sub-tree spanned by the given node.
   * \param nodeP The sub-tree root.
   */
  void _rotate_right (Node* nodeP);

  /*!
   * Duplicate the entire sub-tree rooted at the given node.
   * \param nodeP The sub-tree root.
   * \return A pointer to the duplicated sub-tree root.
   */
  Node* _duplicate (const Node* nodeP);

  /*!
   * Destroy the entire sub-tree rooted at the given node.
   * \param nodeP The sub-tree root.
   */
  void _destroy (Node* nodeP);

  /*!
   * Fix-up the red-black tree properties after an insertion operation.
   * \param nodeP The node that has just been inserted to the tree.
   */
  void _insert_fixup (Node* nodeP);

  /*!
   * Fix-up the red-black tree properties after a removal operation.
   * \param nodeP The child of the node that has just been removed from
   *              the tree (may be a nullptr node).
   * \param parentP The parent node of nodeP (as nodeP may be a nullptr node,
   *                we have to specify its parent explicitly).
   */
  void _remove_fixup (Node* nodeP, Node* parentP);

  /*!
   * Allocate and initialize new tree node.
   * \param object The object stored in the node.
   * \param color The node color.
   * \return A pointer to the newly created node.
   */
  Node* _allocate_node (const Type& object,
                        typename Node::Node_color color);

  /*!
   * De-allocate a tree node.
   * \param nodeP The object node to be deallocated.
   */
  void _deallocate_node (Node* nodeP);
  //@}
};

//---------------------------------------------------------
// Default constructor.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
Multiset<Type, Compare, Allocator, UseCompactContainer>::Multiset () :
  rootP (nullptr),
  iSize (0),
  iBlackHeight (0),
  comp_f ()
{
  // Mark the two fictitious nodes as dummies.
  beginNode.color = Node::DUMMY_BEGIN;
  endNode.color = Node::DUMMY_END;
}

//---------------------------------------------------------
// Constructor with a pointer to comparison object.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
Multiset<Type, Compare, Allocator, UseCompactContainer>::Multiset (const Compare& comp) :
  rootP (nullptr),
  iSize (0),
  iBlackHeight (0),
  comp_f (comp)
{
  // Mark the two fictitious nodes as dummies.
  beginNode.color = Node::DUMMY_BEGIN;
  endNode.color = Node::DUMMY_END;
}

//---------------------------------------------------------
// Copy constructor.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
Multiset<Type, Compare, Allocator, UseCompactContainer>::Multiset (const Self& tree) :
  rootP (nullptr),
  iSize (tree.iSize),
  iBlackHeight (tree.iBlackHeight),
  comp_f (tree.comp_f)
{
  // Mark the two fictitious nodes as dummies.
  beginNode.color = Node::DUMMY_BEGIN;
  endNode.color = Node::DUMMY_END;

  // Copy all the copied tree's nodes recursively.
  if (tree.rootP != nullptr)
  {
    rootP = _duplicate (tree.rootP);

    // Set the dummy nodes.
    beginNode.parentP = _sub_minimum (rootP);
    beginNode.parentP->leftP = &beginNode;
    endNode.parentP = _sub_maximum (rootP);
    endNode.parentP->rightP = &endNode;
  }
  else
  {
    beginNode.parentP = nullptr;
    endNode.parentP = nullptr;
  }
}

//---------------------------------------------------------
// Destructor.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
Multiset<Type, Compare, Allocator, UseCompactContainer>::~Multiset ()
{
  if (UseCompactContainer::value)
    return;

  // Delete the entire tree recursively.
  if (rootP != nullptr)
    _destroy (rootP);

  rootP = nullptr;
  beginNode.parentP = nullptr;
  endNode.parentP = nullptr;
}

//---------------------------------------------------------
// Assignment operator.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
Multiset<Type, Compare, Allocator, UseCompactContainer>&
Multiset<Type, Compare, Allocator, UseCompactContainer>::operator= (const Self& tree)
{
  // Avoid self-assignment.
  if (this == &tree)
    return (*this);

  // Free all objects currently stored in the tree.
  clear();

  // Update the number of objects stored in the tree.
  iSize = tree.iSize;
  iBlackHeight = tree.iBlackHeight;

  // Copy all the copied tree's nodes recursively.
  if (tree.rootP != nullptr)
  {
    rootP = _duplicate (tree.rootP);

    // Set the dummy nodes.
    beginNode.parentP = _sub_minimum (rootP);
    beginNode.parentP->leftP = &beginNode;
    endNode.parentP = _sub_maximum (rootP);
    endNode.parentP->rightP = &endNode;
  }
  else
  {
    beginNode.parentP = nullptr;
    endNode.parentP = nullptr;
  }

  return (*this);
}

//---------------------------------------------------------
// Swap two trees (replace their contents).
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::swap (Self& tree)
{
  // Avoid self-swapping.
  if (this == &tree)
    return;

  // Replace the contents of the trees.
  Node            *tempP = rootP;
  rootP = tree.rootP;
  tree.rootP = tempP;

  size_t           iTemp = iSize;
  iSize = tree.iSize;
  tree.iSize = iTemp;

  iTemp = iBlackHeight;
  iBlackHeight = tree.iBlackHeight;
  tree.iBlackHeight = iTemp;

  // Update the fictitious begin and end nodes.
  tempP = beginNode.parentP;
  beginNode.parentP = tree.beginNode.parentP;
  if (beginNode.parentP != nullptr)
    beginNode.parentP->leftP = &beginNode;
  tree.beginNode.parentP = tempP;
  if (tree.beginNode.parentP != nullptr)
    tree.beginNode.parentP->leftP = &(tree.beginNode);

  tempP = endNode.parentP;
  endNode.parentP = tree.endNode.parentP;
  if (endNode.parentP != nullptr)
    endNode.parentP->rightP = &endNode;
  tree.endNode.parentP = tempP;
  if (tree.endNode.parentP != nullptr)
    tree.endNode.parentP->rightP = &(tree.endNode);

  return;
}

//---------------------------------------------------------
// Test two trees for equality.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
bool Multiset<Type,Compare,Allocator,UseCompactContainer>::operator== (const Self& tree) const
{
  // The sizes of the two trees must be the same.
  if (size() != tree.size())
    return (false);

  // Go over all elements in both tree and compare them pairwise.
  const_iterator   it1 = this->begin();
  const_iterator   it2 = tree.begin();

  while (it1 != this->end() && it2 != tree.end())
  {
    if (comp_f (*it1, *it2) != EQUAL)
      return (false);

    ++it1;
    ++it2;
  }

  // If we reached here, the two trees are equal.
  return (true);
}

//---------------------------------------------------------
// Check if our tree is lexicographically smaller that a given tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
bool Multiset<Type,Compare,Allocator,UseCompactContainer>::operator< (const Self& tree) const
{
  // Go over all elements in both tree and compare them pairwise.
  const_iterator   it1 = this->begin();
  const_iterator   it2 = tree.begin();

  while (it1 != this->end() && it2 != tree.end())
  {
    const Comparison_result   res = comp_f (*it1, *it2);

    if (res == SMALLER)
      return (true);
    if (res == LARGER)
      return (false);

    ++it1;
    ++it2;
  }

  // If we reached here, one tree is the prefix of the other tree. We now
  // check which tree contains more elements.
  if (it1 != this->end())
  {
    // Our tree contains the other tree as a prefix, so it is not smaller.
    return (false);
  }
  else if (it2 != tree.end())
  {
    // The other tree contains our tree as a prefix, so our tree is smaller.
    return (true);
  }

  // The two trees are equal:
  return (false);
}

//---------------------------------------------------------
// Get an iterator for the minimum object in the tree (non-const version).
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
inline typename Multiset<Type,Compare,Allocator,UseCompactContainer>::iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::begin ()
{
  if (beginNode.parentP != nullptr)
    return (iterator (beginNode.parentP));
  else
    return (iterator (&endNode));
}

//---------------------------------------------------------
// Get a past-the-end iterator for the tree objects (non-const version).
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
inline typename Multiset<Type, Compare,Allocator, UseCompactContainer>::iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::end ()
{
  return (iterator (&endNode));
}

//---------------------------------------------------------
// Get an iterator for the minimum object in the tree (const version).
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
inline typename Multiset<Type,Compare,Allocator,UseCompactContainer>::const_iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::begin () const
{
  if (beginNode.parentP != nullptr)
    return (const_iterator (beginNode.parentP));
  else
    return (const_iterator (&endNode));
}

//---------------------------------------------------------
// Get a past-the-end iterator for the tree objects (const version).
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
inline typename Multiset<Type,Compare,Allocator,UseCompactContainer>::const_iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::end () const
{
  return (const_iterator (&endNode));
}

//---------------------------------------------------------
// Get a reverse iterator for the maxnimum object in the tree
// (non-const version).
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
inline typename Multiset<Type,Compare,Allocator,UseCompactContainer>::reverse_iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::rbegin ()
{
  return (reverse_iterator (end()));
}

//---------------------------------------------------------
// Get a pre-the-begin reverse iterator for the tree objects
// (non-const version).
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
inline typename Multiset<Type,Compare,Allocator,UseCompactContainer>::reverse_iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::rend ()
{
  return (reverse_iterator (begin()));
}

//---------------------------------------------------------
// Get a reverse iterator for the maximum object in the tree (const version).
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
inline typename Multiset<Type,Compare,Allocator,UseCompactContainer>::const_reverse_iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::rbegin () const
{
  return (const_reverse_iterator (end()));
}

//---------------------------------------------------------
// Get a pre-the-begin reverse iterator for the tree objects (const version).
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
inline typename Multiset<Type,Compare,Allocator,UseCompactContainer>::const_reverse_iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::rend () const
{
  return (const_reverse_iterator (begin()));
}

//---------------------------------------------------------
// Get the size of the tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
size_t Multiset<Type, Compare, Allocator, UseCompactContainer>::size () const
{
  if (rootP == nullptr)
    // The tree is empty:
    return (0);
  else if (iSize > 0)
    return (iSize);

  // If we reached here, the tree is the result of a split operation and its
  // size is not known - compute it now.
  const Node    *nodeP = beginNode.parentP;
  size_t         iComputedSize = 0;

  while (nodeP != &endNode)
  {
    nodeP = nodeP->successor();
    iComputedSize++;
  }

  // Assign the computed size.
  Self  *myThis = const_cast<Self*> (this);
  myThis->iSize = iComputedSize;

  return (iComputedSize);
}

//---------------------------------------------------------
// Insert a new object to the tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
typename Multiset<Type, Compare, Allocator, UseCompactContainer>::iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::insert (const Type& object)
{
  if (rootP == nullptr)
  {
    // In case the tree is empty, assign a new rootP.
    // Notice that the root is always black.
    rootP = _allocate_node (object, Node::BLACK);

    iSize = 1;
    iBlackHeight = 1;

    // As the tree now contains a single node, it is both the tree
    // maximum and minimum.
    beginNode.parentP = rootP;
    rootP->leftP = &beginNode;
    endNode.parentP = rootP;
    rootP->rightP = &endNode;

    return (iterator (rootP));
  }

  // Find a place for the new object, and insert it as a red leaf.
  Node                      *currentP = rootP;
  Node                      *newNodeP = _allocate_node (object, Node::RED);

  Comparison_result         comp_res;
  bool                      is_leftmost = true;
  bool                      is_rightmost = true;

  while (_is_valid (currentP))
  {
    // Compare the inserted object with the object stored in the current node.
    comp_res = comp_f (object, currentP->object);

    if (comp_res == SMALLER)
    {
      is_rightmost = false;

      if (! _is_valid (currentP->leftP))
      {
        // Insert the new leaf as the left child of the current node.
        currentP->leftP = newNodeP;
        newNodeP->parentP = currentP;
        currentP = nullptr;        // In order to terminate the while loop.

        if (is_leftmost)
        {
          // Assign a new tree minimum.
          beginNode.parentP = newNodeP;
          newNodeP->leftP = &beginNode;
        }
      }
      else
      {
        // Go to the left sub-tree.
        currentP = currentP->leftP;
      }
    }
    else
    {
      is_leftmost = false;

      if (! _is_valid (currentP->rightP))
      {
        // Insert the new leaf as the right child of the current node.
        currentP->rightP = newNodeP;
        newNodeP->parentP = currentP;
        currentP = nullptr;        // In order to terminate the while loop.

        if (is_rightmost)
        {
          // Assign a new tree maximum.
          endNode.parentP = newNodeP;
          newNodeP->rightP = &endNode;
        }
      }
      else
      {
        // Go to the right sub-tree.
        currentP = currentP->rightP;
      }
    }
  }

  // Mark that a new node was added.
  if (iSize > 0)
    iSize++;

  // Fix up the tree properties.
  _insert_fixup (newNodeP);

  return (iterator(newNodeP));
}

//---------------------------------------------------------
// Insert an object to the tree, with a given hint to its position.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
typename Multiset<Type, Compare, Allocator, UseCompactContainer>::iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::insert (iterator position,
                                            const Type& object)
{
  Node  *nodeP = position.nodeP;

  CGAL_multiset_precondition (_is_valid (nodeP));

  // Compare the object to the one stored at the given node in order to decide
  // in which direction to proceed.
  const size_t       max_steps = static_cast<size_t> (1.5 * iBlackHeight);
  Comparison_result  res = comp_f(object, nodeP->object);
  bool               found_pos = true;
  size_t             k = 0;

  if (res == EQUAL)
    return (insert_after (position, object));

  if (res == SMALLER)
  {
    // Go back until locating an object not greater than the inserted one.
    Node  *predP = nodeP->predecessor();

    while (_is_valid (predP) &&
           comp_f (object, predP->object) == SMALLER)
    {
      k++;
      if (k > max_steps)
      {
        // In case the given position is too far away (more than log(n) steps)
        // from the true poisition of the object, break the loop.
        found_pos = false;
        break;
      }

      nodeP = predP;
      predP = nodeP->predecessor();
    }

    if (found_pos)
      return (insert_before (iterator (nodeP), object));
  }
  else
  {
    // Go forward until locating an object not less than the inserted one.
    Node  *succP = nodeP->successor();

    while (_is_valid (succP) &&
           comp_f (object, succP->object) == LARGER)
    {
      k++;
      if (k > max_steps)
      {
        // In case the given position is too far away (more than log(n) steps)
        // from the true poisition of the object, break the loop.
        found_pos = false;
        break;
      }

      nodeP = succP;
      succP = nodeP->successor();
    }

    if (found_pos)
      return (insert_after (iterator (nodeP), object));
  }

  // If the hint was too far than the actual position, perform a normal
  // insertion:
  return (insert (object));
}

//---------------------------------------------------------
// Insert a new object to the tree as the a successor of a given node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
typename Multiset<Type, Compare, Allocator, UseCompactContainer>::iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::insert_after (iterator position,
                                                  const Type& object)
{
  Node  *nodeP = position.nodeP;

  // In case we are given a nullptr node, object should be the tree minimum.
  CGAL_multiset_assertion (nodeP != &endNode);

  if (nodeP == &beginNode)
    nodeP = nullptr;

  if (rootP == nullptr)
  {
    // In case the tree is empty, make sure that we did not receive a valid
    // iterator.
    CGAL_multiset_precondition (nodeP == nullptr);

    // Assign a new root node. Notice that the root is always black.
    rootP = _allocate_node (object, Node::BLACK);

    iSize = 1;
    iBlackHeight = 1;

    // As the tree now contains a single node, it is both the tree
    // maximum and minimum:
    beginNode.parentP = rootP;
    rootP->leftP = &beginNode;
    endNode.parentP = rootP;
    rootP->rightP = &endNode;

    return (iterator (rootP));
  }

  // Insert the new object as a red leaf, being the successor of nodeP.
  Node        *parentP;
  Node        *newNodeP = _allocate_node (object, Node::RED);

  if (nodeP == nullptr)
  {
    // The new node should become the tree minimum: Place is as the left
    // child of the current minimal leaf.
    parentP = beginNode.parentP;

    CGAL_multiset_precondition (comp_f(object, parentP->object) != LARGER);

    parentP->leftP = newNodeP;

    // As we inserted a new tree minimum:
    beginNode.parentP = newNodeP;
    newNodeP->leftP = &beginNode;
  }
  else
  {
    // Make sure the insertion does not violate the tree order.
    CGAL_multiset_precondition_code (Node *_succP = nodeP->successor());
    CGAL_multiset_precondition (comp_f(object, nodeP->object) != SMALLER);
    CGAL_multiset_precondition (! _succP->is_valid() ||
                                comp_f(object, _succP->object) != LARGER);

    // In case given node has no right child, place the new node as its
    // right child. Otherwise, place it at the leftmost position at the
    // sub-tree rooted at its right side.
    if (! _is_valid (nodeP->rightP))
    {
      parentP = nodeP;
      parentP->rightP = newNodeP;
    }
    else
    {
      parentP = _sub_minimum (nodeP->rightP);
      parentP->leftP = newNodeP;
    }

    if (nodeP == endNode.parentP)
    {
      // As we inserted a new tree maximum:
      endNode.parentP = newNodeP;
      newNodeP->rightP = &endNode;
    }
  }

  newNodeP->parentP = parentP;

  // Mark that a new node was added.
  if (iSize > 0)
    iSize++;

  // Fix up the tree properties.
  _insert_fixup (newNodeP);

  return (iterator (newNodeP));
}

//---------------------------------------------------------
// Insert a new object to the tree as the a predecessor of a given node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
typename Multiset<Type, Compare, Allocator, UseCompactContainer>::iterator
Multiset<Type, Compare, Allocator, UseCompactContainer>::insert_before (iterator position,
                                                   const Type& object)
{
  Node  *nodeP = position.nodeP;

  // In case we are given a nullptr node, object should be the tree maximum.
  CGAL_multiset_assertion (nodeP != &beginNode);

  if (nodeP == &endNode)
    nodeP = nullptr;

  if (rootP == nullptr)
  {
    // In case the tree is empty, make sure that we did not receive a valid
    // iterator.
    CGAL_multiset_precondition (nodeP == nullptr);

    // Assign a new root node. Notice that the root is always black.
    rootP = _allocate_node(object, Node::BLACK);

    iSize = 1;
    iBlackHeight = 1;

    // As the tree now contains a single node, it is both the tree
    // maximum and minimum:
    beginNode.parentP = rootP;
    rootP->leftP = &beginNode;
    endNode.parentP = rootP;
    rootP->rightP = &endNode;

    return (iterator (rootP));
  }

  // Insert the new object as a red leaf, being the predecessor of nodeP.
  Node        *parentP;
  Node        *newNodeP = _allocate_node (object, Node::RED);

  if (nodeP == nullptr)
  {
    // The new node should become the tree maximum: Place is as the right
    // child of the current maximal leaf.
    parentP = endNode.parentP;

    CGAL_multiset_precondition (comp_f(object, parentP->object) != SMALLER);

    parentP->rightP = newNodeP;

    // As we inserted a new tree maximum:
    endNode.parentP = newNodeP;
    newNodeP->rightP = &endNode;
  }
  else
  {
    // Make sure the insertion does not violate the tree order.
    CGAL_multiset_precondition_code (Node *_predP = nodeP->predecessor());
    CGAL_multiset_precondition (comp_f(object, nodeP->object) != LARGER);
    CGAL_multiset_precondition (! _predP->is_valid() ||
                                comp_f(object, _predP->object) != SMALLER);

    // In case given node has no left child, place the new node as its
    // left child. Otherwise, place it at the rightmost position at the
    // sub-tree rooted at its left side.
    if (! _is_valid (nodeP->leftP))
    {
      parentP = nodeP;
      parentP->leftP = newNodeP;
    }
    else
    {
      parentP = _sub_maximum (nodeP->leftP);
      parentP->rightP = newNodeP;
    }

    if (nodeP == beginNode.parentP)
    {
      // As we inserted a new tree minimum:
      beginNode.parentP = newNodeP;
      newNodeP->leftP = &beginNode;
    }
  }

  newNodeP->parentP = parentP;

  // Mark that a new node was added.
  if (iSize > 0)
    iSize++;

  // Fix up the tree properties.
  _insert_fixup (newNodeP);

  return (iterator (newNodeP));
}

//---------------------------------------------------------
// Remove an object from the tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
size_t Multiset<Type, Compare, Allocator, UseCompactContainer>::erase (const Type& object)
{
  // Find the first node containing an object not less than the object to
  // be erased and from there look for objects equivalent to the given object.
  size_t      n_removed = 0;
  bool        is_equal;
  Node       *nodeP = _bound (LOWER_BOUND, object, comp_f, is_equal);
  Node       *succP;

  if (! is_equal)
    return (n_removed);

  while (_is_valid (nodeP) && comp_f (object, nodeP->object) == EQUAL)
  {
    // Keep a pointer to the successor node.
    succP = nodeP->successor();

    // Remove the current node.
    _remove_at (nodeP);
    n_removed++;

    // Proceed to the successor.
    nodeP = succP;
  }

  return (n_removed);
}

//---------------------------------------------------------
// Remove the object pointed by the given iterator.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::erase (iterator position)
{
  Node  *nodeP = position.nodeP;

  CGAL_multiset_precondition (_is_valid (nodeP));

  _remove_at (nodeP);
  return;
}

//---------------------------------------------------------
// Remove all objects from the tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::clear ()
{
  // Delete all the tree nodes recursively.
  if (rootP != nullptr)
    _destroy (rootP);

  rootP = nullptr;
  beginNode.parentP = nullptr;
  endNode.parentP = nullptr;

  // Mark that there are no more objects in the tree.
  iSize = 0;
  iBlackHeight = 0;

  return;
}

//---------------------------------------------------------
// Replace the object pointed by a given iterator with another object.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::replace (iterator position,
                                                  const Type& object)
{
  Node  *nodeP = position.nodeP;

  CGAL_multiset_precondition (_is_valid (nodeP));

  // Make sure the replacement does not violate the tree order.
  CGAL_multiset_precondition_code (Node *_succP = nodeP->successor());
  CGAL_multiset_precondition (_succP == nullptr ||
                              _succP->color == Node::DUMMY_END ||
                              comp_f(object, _succP->object) != LARGER);

  CGAL_multiset_precondition_code (Node *_predP = nodeP->predecessor());
  CGAL_multiset_precondition (_predP == nullptr ||
                              _predP->color == Node::DUMMY_BEGIN ||
                              comp_f(object, _predP->object) != SMALLER);

  // Replace the object at nodeP.
  nodeP->object = object;

  return;
}

//---------------------------------------------------------
// Swap the location two objects in the tree, given by their positions.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::swap (iterator pos1,
                                               iterator pos2)
{
  Node        *node1_P = pos1.nodeP;
  Node        *node2_P = pos2.nodeP;

  CGAL_multiset_precondition (_is_valid (node1_P));
  CGAL_multiset_precondition (_is_valid (node2_P));

  if (node1_P == node2_P)
    return;

  // Make sure the swap does not violate the tree order.
  CGAL_multiset_precondition_code (Node *_succ1_P = node1_P->successor());
  CGAL_multiset_precondition (! _is_valid (_succ1_P) ||
                              comp_f (node2_P->object,
                                      _succ1_P->object) != LARGER);

  CGAL_multiset_precondition_code (Node *_pred1_P = node1_P->predecessor());
  CGAL_multiset_precondition (! _is_valid (_pred1_P) ||
                              comp_f (node2_P->object,
                                      _pred1_P->object) != SMALLER);

  CGAL_multiset_precondition_code (Node *_succ2_P = node2_P->successor());
  CGAL_multiset_precondition (! _is_valid (_succ2_P) ||
                              comp_f (node1_P->object,
                                      _succ2_P->object) != LARGER);

  CGAL_multiset_precondition_code (Node *_pred2_P = node2_P->predecessor());
  CGAL_multiset_precondition (! _is_valid (_pred2_P) ||
                              comp_f (node1_P->object,
                                      _pred2_P->object) != SMALLER);

  // Perform the swap.
  if (node1_P->parentP == node2_P->parentP)
    _swap_siblings (node1_P, node2_P);
  else
    _swap (node1_P, node2_P);

  return;
}

//---------------------------------------------------------
// Check if the tree is a valid one.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
bool Multiset<Type, Compare, Allocator, UseCompactContainer>::is_valid () const
{
  if (rootP == nullptr)
  {
    // If there is no root, make sure that the tree is empty.
    if (iSize != 0 || iBlackHeight != 0)
      return (false);

    if (beginNode.parentP != nullptr || endNode.parentP != nullptr)
      return (false);

    return (true);
  }

  // Check the validity of the fictitious nodes.
  if (beginNode.parentP == nullptr || beginNode.parentP->leftP != &beginNode)
    return (false);

  if (endNode.parentP == nullptr || endNode.parentP->rightP != &endNode)
    return (false);

  // Check recursively whether the tree is valid.
  size_t           iComputedSize;
  size_t           iComputedBHeight;

  if (! _sub_is_valid (rootP, iComputedSize, iComputedBHeight))
    return (false);

  // Make sure the tree properties are correct.
  if (iSize != 0)
  {
    if (iSize != iComputedSize)
      return (false);
  }
  else
  {
    // Assign the computed size.
    Self  *myThis = const_cast<Self*> (this);
    myThis->iSize = iComputedSize;
  }

  if (iBlackHeight != iComputedBHeight)
    return (false);

  // If we reached here, the entire tree is valid.
  return (true);
}

//---------------------------------------------------------
// Get the height of the tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
size_t Multiset<Type, Compare, Allocator, UseCompactContainer>::height () const
{
  if (rootP == nullptr)
    // Empty tree.
    return (0);

  // Return the height of the root's sub-tree (the entire tree).
  return (_sub_height (rootP));
}

//---------------------------------------------------------
// Catenate the tree with another given tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::catenate (Self& tree)
{
  // Get the maximal node in our tree and the minimal node in the other tree.
  Node    *max1_P = endNode.parentP;
  Node    *min2_P = tree.beginNode.parentP;

  if (min2_P == nullptr)
  {
    // The other tree is empty - nothing to do.
    CGAL_multiset_assertion (tree.rootP == nullptr);
    return;
  }
  else if (max1_P == nullptr)
  {
    // Our tree is empty: Copy all other tree properties to our tree.
    CGAL_multiset_assertion (rootP == nullptr);

    _shallow_assign (tree);
    return;
  }

  // Make sure that the minimal object in the other tree is not less than the
  // maximal object in our tree.
  CGAL_multiset_precondition (comp_f (max1_P->object,
                                      min2_P->object) != LARGER);

  // Make sure both tree roots black.
  CGAL_multiset_assertion (_is_black (rootP));
  CGAL_multiset_assertion (_is_black (tree.rootP));

  // Splice max1_P (or min2_P) from its tree, but without deleting it.
  Node*   auxP = nullptr;

  if (max1_P != rootP)
  {
    // Splice max1_P from its current poisition in our tree.
    // We know it is has no right child, so we just have to connect its
    // left child with its parent.
    max1_P->parentP->rightP = max1_P->leftP;

    if (_is_valid (max1_P->leftP))
      max1_P->leftP->parentP = max1_P->parentP;

    // If max1_P is a black node, we have to fixup our tree.
    if (max1_P->color == Node::BLACK)
      _remove_fixup (max1_P->leftP, max1_P->parentP);

    auxP = max1_P;
  }
  else if (min2_P != tree.rootP)
  {
    // Splice min2_P from its current poisition in the other tree.
    // We know it is has no left child, so we just have to connect its
    // right child with its parent.
    if (min2_P->parentP != nullptr)
    {
      min2_P->parentP->leftP = min2_P->rightP;

      if (_is_valid (min2_P->rightP))
        min2_P->rightP->parentP = min2_P->parentP;

      // If min2_P is a black node, we have to fixup the other tree.
      if (min2_P->color == Node::BLACK)
        tree._remove_fixup (min2_P->rightP, min2_P->parentP);
    }

    auxP = min2_P;
  }
  else
  {
    // Both nodes are root nodes: Assign min2_P as the right child of the
    // root and color it red.
    rootP->rightP = min2_P;
    min2_P->parentP = rootP;
    min2_P->color = Node::RED;
    min2_P->leftP = nullptr;

    if (! _is_valid (min2_P->rightP))
    {
      // The other tree is of size 1 - set min2_P as the tree maximum.
      min2_P->rightP = &endNode;
      endNode.parentP = min2_P;
    }
    else
    {
      // The other tree is of size 2 - fixup from min2_P's right child and set
      // it as the tree maximum.
      _insert_fixup (min2_P->rightP);
      min2_P->rightP->rightP = &endNode;
      endNode.parentP = min2_P->rightP;
    }

    if (iSize > 0 && tree.iSize > 0)
      iSize += tree.iSize;
    else
      iSize = 0;

    // Clear the other tree (without actually deallocating the nodes).
    tree._shallow_clear();

    return;
  }

  // Mark that the maximal node in our tree is no longer the maximum.
  if (endNode.parentP != nullptr)
    endNode.parentP->rightP = nullptr;

  // Mark that the minimal node in the other tree is no longer the minimum.
  if (tree.beginNode.parentP != nullptr)
    tree.beginNode.parentP->leftP = nullptr;

  // Locate node1_P along the rightmost path in our tree and node2_P along the
  // leftmost path in the other tree, both having the same black-height.
  Node    *node1_P = rootP;
  Node    *node2_P = tree.rootP;
  size_t   iCurrBHeight;

  if (iBlackHeight <= tree.iBlackHeight)
  {
    // The other tree is taller (or both trees have the same height):
    // Go down the leftmost path of the other tree and locate node2_P.
    node2_P = tree.rootP;
    iCurrBHeight = tree.iBlackHeight;
    while (iCurrBHeight > iBlackHeight)
    {
      if (node2_P->color == Node::BLACK)
        iCurrBHeight--;
      node2_P = node2_P->leftP;
    }

    if (_is_red (node2_P))
      node2_P = node2_P->leftP;
    CGAL_multiset_assertion (_is_valid (node2_P));
  }
  else
  {
    // Our tree is taller:
    // Go down the rightmost path of our tree and locate node1_P.
    node1_P = rootP;
    iCurrBHeight = iBlackHeight;
    while (iCurrBHeight > tree.iBlackHeight)
    {
      if (node1_P->color == Node::BLACK)
        iCurrBHeight--;
      node1_P = node1_P->rightP;
    }

    if (_is_red (node1_P))
      node1_P = node1_P->rightP;
    CGAL_multiset_assertion (_is_valid (node2_P));
  }

  // Check which one of the tree roots have we reached.
  Node    *newRootP = nullptr;
  Node    *parentP;

  if (node1_P == rootP)
  {
    // We know that node2_P has the same number of black roots from it to
    // the minimal tree node as the number of black nodes from our tree root
    // to any of its leaves. We make rootP and node2_P siblings by moving
    // auxP to be their parent.
    parentP = node2_P->parentP;

    if (parentP == nullptr)
    {
      // Make auxP the root of the catenated tree.
      newRootP = auxP;
    }
    else
    {
      // The catenated tree will be rooted at the other tree's root.
      newRootP = tree.rootP;

      // Move auxP as the left child of parentP.
      parentP->leftP = auxP;
    }
  }
  else
  {
    // We know that node1_P has the same number of black roots from it to
    // the maximal tree node as the number of black nodes from the other tree
    // root to any of its leaves. We make tree.rootP and node1_P siblings by
    // moving auxP to be their parent.
    parentP = node1_P->parentP;

    CGAL_multiset_assertion (parentP != nullptr);

    // The catenated tree will be rooted at the current root of our tree.
    newRootP = rootP;

    // Move auxP as the right child of parentP.
    parentP->rightP = auxP;
  }

  // Move node1_P to be the left child of auxP, and node2_P to be its
  // right child. We also color this node red.
  auxP->parentP = parentP;
  auxP->color = Node::RED;
  auxP->leftP = node1_P;
  auxP->rightP = node2_P;

  node1_P->parentP = auxP;
  node2_P->parentP = auxP;

  // Set the catenated tree properties.
  if (rootP != newRootP)
  {
    // Take the black-height of the other tree.
    iBlackHeight = tree.iBlackHeight;
    rootP = newRootP;
  }

  if (iSize > 0 && tree.iSize > 0)
    iSize += tree.iSize;
  else
    iSize = 0;

  // Set the new maximal node in the tree (the minimal node remains unchanged).
  endNode.parentP = tree.endNode.parentP;
  endNode.parentP->rightP = &endNode;

  // Clear the other tree (without actually deallocating the nodes).
  tree._shallow_clear();

  // Perform a fixup of the tree invariants. This fixup will also take care
  // of updating the black-height of our catenated tree.
  _insert_fixup (auxP);

  return;
}

//---------------------------------------------------------
// Split the tree at a given position, such that it contains all objects
// in the range [begin, position) and all objects in the range
// [position, end) form a new output tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::split (iterator position,
                                                Self& tree)
{
  CGAL_multiset_precondition (tree.empty());

  // Check the extremal cases.
  if (position == begin())
  {
    // The tree should be copied to the output tree.
    tree._shallow_assign (*this);
    return;
  }
  else if (position == end())
  {
    // Nothing to do - all objects should remain in our tree.
    return;
  }

  // Prepare a vector describing the path from the given node to the tree
  // root (where SMALLER designates a left turn and LARGER designates a
  // right turn). Note that the length of this path (the depth of nodeP)
  // is at most twice the black-height of the tree.
  Node    *nodeP = position.nodeP;

  CGAL_multiset_precondition (_is_valid (nodeP));

  Node              *currP = nodeP;
  Comparison_result *path = new Comparison_result [2 * iBlackHeight];
  int                depth = 0;

  path[depth] = EQUAL;
  while (currP->parentP != nullptr)
  {
    depth++;
    if (currP == currP->parentP->leftP)
      path[depth] = SMALLER;
    else
      path[depth] = LARGER;

    currP = currP->parentP;
  }
  CGAL_multiset_assertion (currP == rootP);

  // Now go down the path and split the tree accordingly. We also keep
  // track of the black-height of the current node.
  size_t  iCurrBHeight = iBlackHeight;
  Self    leftTree;
  size_t  iLeftBHeight = 0;
  Node   *spineLeftP = nullptr;
  Node   *auxLeftP = nullptr;
  Self    rightTree;
  size_t  iRightBHeight = 0;
  Node   *spineRightP = nullptr;
  Node   *auxRightP = nullptr;
  Node   *childP = nullptr;
  Node   *nextP = nullptr;

  while (depth >= 0)
  {
    CGAL_multiset_assertion (_is_valid (currP));

    // If we encounter a black node, the black-height of both its left and
    // right subtrees is decremented.
    if (_is_black (currP))
      iCurrBHeight--;

    // Check in which direction we have to go.
    if (path[depth] != LARGER)
    {
      // We go left, so currP and its entire right sub-tree (T_r) should be
      // included in the right split tree.
      //
      //              (.) currP           .
      //              / \                 .
      //             /   \                .
      //           (.)    T_r             .
      //
      // This also covers that case where currP is the split node (nodeP).
      childP = currP->rightP;
      nextP = currP->leftP;

      if (_is_valid (childP) && rightTree.rootP == nullptr)
      {
        // Assing T_r to rightTree.
        rightTree.rootP = childP;
        rightTree.iBlackHeight = iCurrBHeight;

        // Make sure the root of rightTree is black.
        rightTree.rootP->parentP = nullptr;
        if (_is_red (rightTree.rootP))
        {
          rightTree.rootP->color = Node::BLACK;
          rightTree.iBlackHeight++;
        }

        // We store a black node along the leftmost spine of rightTree whose
        // black-hieght is exactly iRightBHeight.
        iRightBHeight = rightTree.iBlackHeight;
        spineRightP = rightTree.rootP;
      }
      else if (_is_valid (childP))
      {
        // Catenate T_r with the current rightTree.
        CGAL_multiset_assertion (_is_valid (spineRightP) &&
                                 _is_valid(auxRightP));

        // Make sure the root of T_r is black.
        size_t    iCurrRightBHeight = iCurrBHeight;

        if (_is_red (childP))
        {
          childP->color = Node::BLACK;
          iCurrRightBHeight++;
        }

        // Go down the leftmost path of rightTree until locating a black
        // node whose black height is exactly iCurrRightBHeight.
        CGAL_multiset_assertion (iRightBHeight >= iCurrRightBHeight);

        while (iRightBHeight > iCurrRightBHeight)
        {
          if (spineRightP->color == Node::BLACK)
            iRightBHeight--;
          spineRightP = spineRightP->leftP;
        }

        if (_is_red (spineRightP))
          spineRightP = spineRightP->leftP;
        CGAL_multiset_assertion (_is_valid (spineRightP));

        // Use the auxiliary node and make it the parent of T_r (which
        // becomes its left sub-tree) and spineRightP (which becomes its
        // right child). We color the auxiliary node red, as both its
        // children are black.
        auxRightP->parentP = spineRightP->parentP;
        auxRightP->color = Node::RED;
        auxRightP->leftP = childP;
        auxRightP->rightP = spineRightP;

        if (auxRightP->parentP != nullptr)
          auxRightP->parentP->leftP = auxRightP;
        else
          rightTree.rootP = auxRightP;

        childP->parentP = auxRightP;
        spineRightP->parentP = auxRightP;

        // Perform a fixup on the right tree.
        rightTree._insert_fixup (auxRightP);
        auxRightP = nullptr;

        // Note that childP is now located on the leftmost spine of
        // rightTree and its black-height is exactly iCurrRightBHeight.
        iRightBHeight = iCurrRightBHeight;
        spineRightP = childP;
      }

      // In case we have an auxiliary right node that has not been inserted
      // into the right tree, insert it now.
      if (auxRightP != nullptr)
      {
        if (rightTree.rootP != nullptr)
        {
          // The right tree is not empty. Traverse its leftmost spine to
          // locate the parent of auxRightP.
          while (_is_valid (spineRightP->leftP))
            spineRightP = spineRightP->leftP;

          auxRightP->parentP = spineRightP;
          auxRightP->color = Node::RED;
          auxRightP->rightP = nullptr;
          auxRightP->leftP = nullptr;

          spineRightP->leftP = auxRightP;

          // Perform a fixup on the right tree, following the insertion of
          // the auxiliary right node.
          rightTree._insert_fixup (auxRightP);
        }
        else
        {
          // auxRightP is the only node in the current right tree.
          rightTree.rootP = auxRightP;
          rightTree.iBlackHeight = 1;

          auxRightP->parentP = nullptr;
          auxRightP->color = Node::BLACK;
          auxRightP->rightP = nullptr;
          auxRightP->leftP = nullptr;
        }

        // Assign spineRightP to be the auxiliary node.
        spineRightP = auxRightP;
        iRightBHeight = (_is_black (spineRightP)) ? 1 : 0;
        auxRightP = nullptr;
      }

      // Mark currP as the auxiliary right node.
      auxRightP = currP;
    }

    if (path[depth] != SMALLER)
    {
      // We go right, so currP and its entire left sub-tree (T_l) should be
      // included in the left split tree.
      //
      //              (.) currP           .
      //              / \                 .
      //             /   \                .
      //          T_l    (.)              .
      //
      // An exception to this rule is when currP is the split node (nodeP),
      // so it should not be included in the left tree.
      childP = currP->leftP;
      nextP = currP->rightP;

      if (_is_valid (childP) && leftTree.rootP == nullptr)
      {
        // Assing T_l to leftTree.
        leftTree.rootP = childP;
        leftTree.iBlackHeight = iCurrBHeight;

        // Make sure the root of leftTree is black.
        leftTree.rootP->parentP = nullptr;
        if (_is_red (leftTree.rootP))
        {
          leftTree.rootP->color = Node::BLACK;
          leftTree.iBlackHeight++;
        }

        // We store a black node along the rightmost spine of leftTree whose
        // black-hieght is exactly iLeftBHeight.
        iLeftBHeight = leftTree.iBlackHeight;
        spineLeftP = leftTree.rootP;
      }
      else if (_is_valid (childP))
      {
        // Catenate T_l with the current leftTree.
        CGAL_multiset_assertion (_is_valid (spineLeftP) &&
                                 _is_valid(auxLeftP));

        // Make sure the root of T_l is black.
        size_t    iCurrLeftBHeight = iCurrBHeight;

        if (_is_red (childP))
        {
          childP->color = Node::BLACK;
          iCurrLeftBHeight++;
        }

        // Go down the rightmost path of leftTree until locating a black
        // node whose black height is exactly iCurrLeftBHeight.
        CGAL_multiset_assertion (iLeftBHeight >= iCurrLeftBHeight);

        while (iLeftBHeight > iCurrLeftBHeight)
        {
          if (spineLeftP->color == Node::BLACK)
            iLeftBHeight--;
          spineLeftP = spineLeftP->rightP;
        }

        if (_is_red (spineLeftP))
          spineLeftP = spineLeftP->rightP;
        CGAL_multiset_assertion (_is_valid (spineLeftP));

        // Use the auxiliary node and make it the parent of T_l (which
        // becomes its right sub-tree) and spineLeftP (which becomes its
        // left child). We color the auxiliary node red, as both its
        // children are black.
        auxLeftP->parentP = spineLeftP->parentP;
        auxLeftP->color = Node::RED;
        auxLeftP->leftP = spineLeftP;
        auxLeftP->rightP = childP;

        if (auxLeftP->parentP != nullptr)
          auxLeftP->parentP->rightP = auxLeftP;
        else
          leftTree.rootP = auxLeftP;

        childP->parentP = auxLeftP;
        spineLeftP->parentP = auxLeftP;

        // Perform a fixup on the left tree.
        leftTree._insert_fixup (auxLeftP);
        auxLeftP = nullptr;

        // Note that childP is now located on the rightmost spine of
        // leftTree and its black-height is exactly iCurrLeftBHeight.
        iLeftBHeight = iCurrLeftBHeight;
        spineLeftP = childP;
      }

      // In case we have an auxiliary left node that has not been inserted
      // into the left tree, insert it now.
      if (auxLeftP != nullptr)
      {
        if (leftTree.rootP != nullptr)
        {
          // The left tree is not empty. Traverse its rightmost spine to
          // locate the parent of auxLeftP.
          while (_is_valid (spineLeftP->rightP))
            spineLeftP = spineLeftP->rightP;

          auxLeftP->parentP = spineLeftP;
          auxLeftP->color = Node::RED;
          auxLeftP->rightP = nullptr;
          auxLeftP->leftP = nullptr;

          spineLeftP->rightP = auxLeftP;

          // Perform a fixup on the left tree, following the insertion of
          // the auxiliary left node.
          leftTree._insert_fixup (auxLeftP);
        }
        else
        {
          // auxLeftP is the only node in the left tree.
          leftTree.rootP = auxLeftP;
          leftTree.iBlackHeight = 1;

          auxLeftP->parentP = nullptr;
          auxLeftP->color = Node::BLACK;
          auxLeftP->rightP = nullptr;
          auxLeftP->leftP = nullptr;
        }

        // Assign spineLeftP to be the auxiliary node.
        spineLeftP = auxLeftP;
        iLeftBHeight = (_is_black (spineLeftP)) ? 1 : 0;
        auxLeftP = nullptr;
      }

      // Mark currP as the auxiliary right node.
      if (depth > 0)
        auxLeftP = currP;
    }

    // Proceed to the next step in the path.
    currP = nextP;
    depth--;
  }

  // It is now possible to free the path.
  delete[] path;

  CGAL_multiset_assertion (auxLeftP == nullptr && auxRightP == nodeP);

  // Fix the properties of the left tree: We know its minimal node is the
  // same as the current minimum.
  leftTree.beginNode.parentP = beginNode.parentP;
  leftTree.beginNode.parentP->leftP = &(leftTree.beginNode);

  // Traverse the rightmost path of the left tree to find the its maximum.
  CGAL_multiset_assertion (_is_valid (spineLeftP));

  while (_is_valid (spineLeftP->rightP))
    spineLeftP = spineLeftP->rightP;

  leftTree.endNode.parentP = spineLeftP;
  spineLeftP->rightP = &(leftTree.endNode);

  // Fix the properties of the right tree: We know its maximal node is the
  // same as the current maximum.
  rightTree.endNode.parentP = endNode.parentP;
  rightTree.endNode.parentP->rightP = &(rightTree.endNode);

  // We still have to insert the split node as the minimum node of the right
  // tree (we can traverse its leftmost path to find its parent).
  if (rightTree.rootP != nullptr)
  {
    while (_is_valid (spineRightP->leftP))
      spineRightP = spineRightP->leftP;

    nodeP->parentP = spineRightP;
    nodeP->color = Node::RED;
    nodeP->rightP = nullptr;
    nodeP->leftP = nullptr;

    spineRightP->leftP = nodeP;

    // Perform a fixup on the right tree, following the insertion of the
    // leftmost node.
    rightTree._insert_fixup (nodeP);
  }
  else
  {
    // nodeP is the only node in the right tree.
    rightTree.rootP = nodeP;
    rightTree.iBlackHeight = 1;

    nodeP->parentP = nullptr;
    nodeP->color = Node::BLACK;
    nodeP->rightP = nullptr;
    nodeP->leftP = nullptr;

    // In this case we also know the tree sizes:
    leftTree.iSize = iSize - 1;
    rightTree.iSize = 1;
  }

  // Set nodeP as the minimal node is the right tree.
  rightTree.beginNode.parentP = nodeP;
  nodeP->leftP = &(rightTree.beginNode);

  // Assign leftTree to (*this) and rightTree to the output tree.
  _shallow_assign (leftTree);

  tree.clear();
  tree._shallow_assign (rightTree);

  return;
}

//---------------------------------------------------------
// Move the contents of one tree to another without actually duplicating
// the nodes. This operation also clears the copied tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_shallow_assign (Self& tree)
{
  // Copy the assigned tree properties.
  rootP = tree.rootP;
  iSize = tree.iSize;
  iBlackHeight = tree.iBlackHeight;

  // Properly mark the minimal and maximal tree nodes.
  beginNode.parentP = tree.beginNode.parentP;

  if (beginNode.parentP != nullptr)
    beginNode.parentP->leftP = &beginNode;

  endNode.parentP = tree.endNode.parentP;

  if (endNode.parentP != nullptr)
    endNode.parentP->rightP = &endNode;

  // Clear the other tree (without actually deallocating the nodes).
  tree._shallow_clear();

  return;
}

//---------------------------------------------------------
// Clear the properties of the tree, without actually deallocating its nodes.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_shallow_clear ()
{
  rootP = nullptr;
  iSize = 0;
  iBlackHeight = 0;
  beginNode.parentP = nullptr;
  endNode.parentP = nullptr;

  return;
}

//---------------------------------------------------------
// Remove the given tree node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_remove_at (Node* nodeP)
{
  CGAL_multiset_precondition (_is_valid (nodeP));

  if (nodeP == rootP &&
      ! _is_valid (rootP->leftP) && ! _is_valid (rootP->rightP))
  {
    // In case of deleting the single object stored in the tree, free the root,
    // thus emptying the tree.
    _deallocate_node (rootP);

    rootP = nullptr;
    beginNode.parentP = nullptr;
    endNode.parentP = nullptr;
    iSize = 0;
    iBlackHeight = 0;

    return;
  }

  // Remove the given node from the tree.
  if (_is_valid (nodeP->leftP) && _is_valid (nodeP->rightP))
  {
    // If the node we want to remove has two children, find its successor,
    // which is the leftmost child in its right sub-tree and has at most
    // one child (it may have a right child).
    Node    *succP = _sub_minimum (nodeP->rightP);
    CGAL_multiset_assertion (_is_valid (succP));

    // Now physically swap nodeP and its successor. Notice this may temporarily
    // violate the tree properties, but we are going to remove nodeP anyway.
    // This way we have moved nodeP to a position were it is more convinient
    // to delete it.
    _swap (nodeP, succP);
  }

  // At this stage, the node we are going to remove has at most one child.
  Node        *childP = nullptr;

  if (_is_valid (nodeP->leftP))
  {
    CGAL_multiset_assertion (! _is_valid (nodeP->rightP));
    childP = nodeP->leftP;
  }
  else
  {
    childP = nodeP->rightP;
  }

  // Splice out the node to be removed, by linking its parent straight to the
  // removed node's single child.
  if (_is_valid (childP))
    childP->parentP = nodeP->parentP;

  if (nodeP->parentP == nullptr)
  {
    // If we are deleting the root, make the child the new tree node.
    rootP = childP;

    // If the deleted root is black, decrement the black height of the tree.
    if (nodeP->color == Node::BLACK)
      iBlackHeight--;
  }
  else
  {
    // Link the removed node parent to its child.
    if (nodeP == nodeP->parentP->leftP)
    {
      nodeP->parentP->leftP = childP;
    }
    else
    {
      nodeP->parentP->rightP = childP;
    }
  }

  // Fix-up the red-black properties that may have been damaged: If we have
  // just removed a black node, the black-height property is no longer valid.
  if (nodeP->color == Node::BLACK)
    _remove_fixup (childP, nodeP->parentP);

  // In case we delete the tree minimum of maximum, update the relevant
  // pointers.
  if (nodeP == beginNode.parentP)
  {
    beginNode.parentP = nodeP->successor();

    if (_is_valid (beginNode.parentP))
      beginNode.parentP->leftP = &beginNode;
    else
      beginNode.parentP = nullptr;
  }
  else if (nodeP == endNode.parentP)
  {
    endNode.parentP = nodeP->predecessor();

    if (_is_valid (endNode.parentP))
      endNode.parentP->rightP = &endNode;
    else
      endNode.parentP = nullptr;
  }

  // Delete the unnecessary node.
  _deallocate_node (nodeP);

  // Decrement the number of objects in the tree.
  if (iSize > 0)
    iSize--;

  return;
}

//---------------------------------------------------------
// Swap the location two nodes in the tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_swap (Node* node1_P,
                                                Node* node2_P)
{
  CGAL_multiset_assertion (_is_valid (node1_P));
  CGAL_multiset_assertion (_is_valid (node2_P));

  // Store the properties of the first node.
  typename Node::Node_color   color1 = node1_P->color;
  Node        *parent1_P = node1_P->parentP;
  Node        *right1_P = node1_P->rightP;
  Node        *left1_P = node1_P->leftP;

  // Copy the properties of the second node to the first node.
  node1_P->color = node2_P->color;

  if (node1_P != node2_P->parentP)
  {
    if (node2_P->parentP == nullptr)
    {
      rootP = node1_P;
    }
    else
    {
      if (node2_P->parentP->leftP == node2_P)
        node2_P->parentP->leftP = node1_P;
      else
        node2_P->parentP->rightP = node1_P;
    }

    node1_P->parentP = node2_P->parentP;
  }
  else
  {
    node1_P->parentP = node2_P;
  }

  if (node1_P != node2_P->rightP)
  {
    if (_is_valid (node2_P->rightP))
      node2_P->rightP->parentP = node1_P;

    node1_P->rightP = node2_P->rightP;
  }
  else
  {
    node1_P->rightP = node2_P;
  }

  if (node1_P != node2_P->leftP)
  {
    if (_is_valid (node2_P->leftP))
      node2_P->leftP->parentP = node1_P;

    node1_P->leftP = node2_P->leftP;
  }
  else
  {
    node1_P->leftP = node2_P;
  }

  // Copy the stored properties of the first node to the second node.
  node2_P->color = color1;

  if (node2_P != parent1_P)
  {
    if (parent1_P == nullptr)
    {
      rootP = node2_P;
    }
    else
    {
      if (parent1_P->leftP == node1_P)
        parent1_P->leftP = node2_P;
      else
        parent1_P->rightP = node2_P;
    }

    node2_P->parentP = parent1_P;
  }
  else
  {
    node2_P->parentP = node1_P;
  }

  if (node2_P != right1_P)
  {
    if (_is_valid (right1_P))
      right1_P->parentP = node2_P;

    node2_P->rightP = right1_P;
  }
  else
  {
    node2_P->rightP = node1_P;
  }

  if (node2_P != left1_P)
  {
    if (_is_valid (left1_P))
      left1_P->parentP = node2_P;

    node2_P->leftP = left1_P;
  }
  else
  {
    node2_P->leftP = node1_P;
  }

  // If one of the swapped nodes used to be the tree minimum, update
  // the properties of the fictitious before-the-begin node.
  if (beginNode.parentP == node1_P)
  {
    beginNode.parentP = node2_P;
    node2_P->leftP = &beginNode;
  }
  else if (beginNode.parentP == node2_P)
  {
    beginNode.parentP = node1_P;
    node1_P->leftP = &beginNode;
  }

  // If one of the swapped nodes used to be the tree maximum, update
  // the properties of the fictitious past-the-end node.
  if (endNode.parentP == node1_P)
  {
    endNode.parentP = node2_P;
    node2_P->rightP = &endNode;
  }
  else if (endNode.parentP == node2_P)
  {
    endNode.parentP = node1_P;
    node1_P->rightP = &endNode;
  }

  return;
}

//---------------------------------------------------------
// Swap the location two sibling nodes in the tree.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_swap_siblings (Node* node1_P,
                                                         Node* node2_P)
{
  CGAL_multiset_assertion (_is_valid (node1_P));
  CGAL_multiset_assertion (_is_valid (node2_P));

  // Store the properties of the first node.
  typename Node::Node_color   color1 = node1_P->color;
  Node        *right1_P = node1_P->rightP;
  Node        *left1_P = node1_P->leftP;

  // Copy the properties of the second node to the first node.
  node1_P->color = node2_P->color;

  node1_P->rightP = node2_P->rightP;
  if (_is_valid (node1_P->rightP))
    node1_P->rightP->parentP = node1_P;

  node1_P->leftP = node2_P->leftP;
  if (_is_valid (node1_P->leftP))
    node1_P->leftP->parentP = node1_P;

  // Copy the stored properties of the first node to the second node.
  node2_P->color = color1;

  node2_P->rightP = right1_P;
  if (_is_valid (node2_P->rightP))
    node2_P->rightP->parentP = node2_P;

  node2_P->leftP = left1_P;
  if (_is_valid (node2_P->leftP))
    node2_P->leftP->parentP = node2_P;

  // Swap the children of the common parent node.
  Node        *parent_P = node1_P->parentP;
  Node        *temp;

  CGAL_multiset_assertion (parent_P == node2_P->parentP);

  temp = parent_P->leftP;
  parent_P->leftP = parent_P->rightP;
  parent_P->rightP = temp;

  // If one of the swapped nodes used to be the tree minimum, update
  // the properties of the fictitious before-the-begin node.
  if (beginNode.parentP == node1_P)
  {
    beginNode.parentP = node2_P;
    node2_P->leftP = &beginNode;
  }
  else if (beginNode.parentP == node2_P)
  {
    beginNode.parentP = node1_P;
    node1_P->leftP = &beginNode;
  }

  // If one of the swapped nodes used to be the tree maximum, update
  // the properties of the fictitious past-the-end node.
  if (endNode.parentP == node1_P)
  {
    endNode.parentP = node2_P;
    node2_P->rightP = &endNode;
  }
  else if (endNode.parentP == node2_P)
  {
    endNode.parentP = node1_P;
    node1_P->rightP = &endNode;
  }

  return;
}

//---------------------------------------------------------
// Calculate the height of the subtree spanned by a given node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
size_t Multiset<Type, Compare, Allocator, UseCompactContainer>::_sub_height
    (const Node* nodeP) const
{
  CGAL_multiset_assertion (_is_valid (nodeP));

  // Recursively calculate the heights of the left and right sub-trees.
  size_t   iRightHeight = 0;
  size_t   iLeftHeight = 0;

  if (_is_valid (nodeP->rightP))
    iRightHeight = _sub_height (nodeP->rightP);

  if (_is_valid (nodeP->leftP))
    iLeftHeight = _sub_height (nodeP->leftP);

  // Return the maximal child sub-height + 1 (the current node).
  return ((iRightHeight > iLeftHeight) ? (iRightHeight + 1) :
                                         (iLeftHeight + 1));
}

//---------------------------------------------------------
// Calculate the height of the subtree spanned by a given node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
bool Multiset<Type, Compare, Allocator, UseCompactContainer>::_sub_is_valid
    (const Node* nodeP,
     size_t& sub_size,
     size_t& sub_bh) const
{
  // Make sure that the node is valid.
  if (! _is_valid (nodeP))
    return (false);

  // If the node is red, make sure that both it children are black (note that
  // nullptr nodes are also considered to be black).
  if (_is_red (nodeP) &&
      (! _is_black (nodeP->rightP) || ! _is_black (nodeP->leftP)))
  {
    return (false);
  }

  // Recursively calculate the black heights of the left and right sub-trees.
  size_t   iBlack = ((nodeP->color == Node::BLACK) ? 1 : 0);
  size_t   iLeftSize = 0;
  size_t   iLeftBHeight = 0;
  size_t   iRightSize = 0;
  size_t   iRightBHeight = 0;

  if (_is_valid (nodeP->leftP))
  {
    // Make sure that the object stored at nodeP is not smaller than the one
    // stored at its left child.
    if (comp_f (nodeP->object, nodeP->leftP->object) == SMALLER)
      return (false);

    // Recursively check that the left sub-tree is valid.
    if (! _sub_is_valid (nodeP->leftP, iLeftSize, iLeftBHeight))
      return (false);
  }

  if (_is_valid (nodeP->rightP))
  {
    // Make sure that the object stored at nodeP is not larger than the one
    // stored at its right child.
    if (comp_f (nodeP->object, nodeP->rightP->object) == LARGER)
      return (false);

    // Recursively check that the right sub-tree is valid.
    if (! _sub_is_valid (nodeP->rightP, iRightSize, iRightBHeight))
      return (false);
  }

  // Compute the size of the entire sub-tree.
  sub_size = iRightSize + iLeftSize + 1;

  // Make sure that the black heights of both sub-trees are equal.
  if (iRightBHeight != iLeftBHeight)
    return (false);
  sub_bh = iRightBHeight + iBlack;

  // If we reached here, the subtree is valid.
  return (true);
}

//---------------------------------------------------------
// Get the leftmost node in the sub-tree spanned by the given node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
typename Multiset<Type, Compare, Allocator, UseCompactContainer>::Node*
Multiset<Type, Compare, Allocator, UseCompactContainer>::_sub_minimum (Node* nodeP) const
{
  CGAL_multiset_assertion (_is_valid (nodeP));

  Node    *minP = nodeP;

  while (_is_valid (minP->leftP))
    minP = minP->leftP;
  return (minP);
}

//---------------------------------------------------------
// Get the rightmost node in the sub-tree spanned by the given node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
typename Multiset<Type, Compare, Allocator, UseCompactContainer>::Node*
Multiset<Type, Compare, Allocator, UseCompactContainer>::_sub_maximum (Node* nodeP) const
{
  CGAL_multiset_assertion (_is_valid (nodeP));

  Node    *maxP = nodeP;

  while (_is_valid (maxP->rightP))
    maxP = maxP->rightP;
  return (maxP);
}

//---------------------------------------------------------
// Left-rotate the sub-tree spanned by the given node:
//
//          |          RoateRight(y)            |
//          y         -------------->           x
//        /   \                               /   \       .
//       x     T3       RoatateLeft(x)       T1    y      .
//     /   \          <--------------            /   \    .
//    T1    T2                                  T2    T3
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_rotate_left (Node* xNodeP)
{
  // Get the right child of the node.
  Node        *yNodeP = xNodeP->rightP;

  CGAL_multiset_assertion (_is_valid (yNodeP));

  // Change its left subtree (T2) to x's right subtree.
  xNodeP->rightP = yNodeP->leftP;

  // Link T2 to its new parent x.
  if (_is_valid (yNodeP->leftP))
    yNodeP->leftP->parentP = xNodeP;

  // Assign x's parent to be y's parent.
  yNodeP->parentP = xNodeP->parentP;

  if (xNodeP->parentP == nullptr)
  {
    // Make y the new tree root.
    rootP = yNodeP;
  }
  else
  {
    // Assign a pointer to y from x's parent.
    if (xNodeP == xNodeP->parentP->leftP)
    {
      xNodeP->parentP->leftP = yNodeP;
    }
    else
    {
      xNodeP->parentP->rightP = yNodeP;
    }
  }

  // Assign x to be y's left child.
  yNodeP->leftP = xNodeP;
  xNodeP->parentP = yNodeP;

  return;
}

//---------------------------------------------------------
// Right-rotate the sub-tree spanned by the given node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_rotate_right (Node* yNodeP)
{
  // Get the left child of the node.
  Node        *xNodeP = yNodeP->leftP;

  CGAL_multiset_assertion (_is_valid (xNodeP));

  // Change its right subtree (T2) to y's left subtree.
  yNodeP->leftP = xNodeP->rightP;

  // Link T2 to its new parent y.
  if (_is_valid (xNodeP->rightP))
    xNodeP->rightP->parentP = yNodeP;

  // Assign y's parent to be x's parent.
  xNodeP->parentP = yNodeP->parentP;

  if (yNodeP->parentP == nullptr)
  {
    // Make x the new tree root.
    rootP = xNodeP;
  }
  else
  {
    // Assign a pointer to x from y's parent.
    if (yNodeP == yNodeP->parentP->leftP)
    {
      yNodeP->parentP->leftP = xNodeP;
    }
    else
    {
      yNodeP->parentP->rightP = xNodeP;
    }
  }

  // Assign y to be x's right child.
  xNodeP->rightP = yNodeP;
  yNodeP->parentP = xNodeP;

  return;
}

//---------------------------------------------------------
// Duplicate the entire sub-tree rooted at the given node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
typename Multiset<Type, Compare, Allocator, UseCompactContainer>::Node*
Multiset<Type, Compare, Allocator, UseCompactContainer>::_duplicate (const Node* nodeP)
{
  CGAL_multiset_assertion (_is_valid (nodeP));

  // Create a node of the same color, containing the same object.
  Node   *dupNodeP = _allocate_node(nodeP->object, nodeP->color);

  // Duplicate the children recursively.
  if (_is_valid (nodeP->rightP))
  {
    dupNodeP->rightP = _duplicate (nodeP->rightP);
    dupNodeP->rightP->parentP = dupNodeP;
  }

  if (_is_valid (nodeP->leftP))
  {
    dupNodeP->leftP = _duplicate (nodeP->leftP);
    dupNodeP->leftP->parentP = dupNodeP;
  }

  // Return the duplicated node.
  return (dupNodeP);
}

//---------------------------------------------------------
// Destroy the entire sub-tree rooted at the given node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_destroy (Node* nodeP)
{
  CGAL_multiset_assertion (_is_valid (nodeP));

  // Destroy the children recursively.
  if (_is_valid (nodeP->rightP))
    _destroy (nodeP->rightP);
  nodeP->rightP = nullptr;

  if (_is_valid (nodeP->leftP))
    _destroy (nodeP->leftP);
  nodeP->leftP = nullptr;

  // Free the subtree root node.
  _deallocate_node (nodeP);

  return;
}

//---------------------------------------------------------
// Fix-up the tree so it maintains the red-black properties after insertion.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_insert_fixup (Node* nodeP)
{
  CGAL_multiset_precondition (_is_red (nodeP));

  // Fix the red-black propreties: we may have inserted a red leaf as the
  // child of a red parent - so we have to fix the coloring of the parent
  // recursively.
  Node        *currP = nodeP;
  Node        *grandparentP;
  Node        *uncleP;

  while (currP != rootP && _is_red (currP->parentP))
  {
    // Get a pointer to the current node's grandparent (notice the root is
    // always black, so the red parent must have a parent).
    grandparentP = currP->parentP->parentP;
    CGAL_multiset_precondition (grandparentP != nullptr);

    if (currP->parentP == grandparentP->leftP)
    {
      // If the red parent is a left child, the uncle is the right child of
      // the grandparent.
      uncleP = grandparentP->rightP;

      if (_is_red (uncleP))
      {
        // If both parent and uncle are red, color them black and color the
        // grandparent red.
        // In case of a nullptr uncle, we treat it as a black node.
        currP->parentP->color = Node::BLACK;
        uncleP->color = Node::BLACK;
        grandparentP->color = Node::RED;

        // Move to the grandparent.
        currP = grandparentP;
      }
      else
      {
        // Make sure the current node is a left child. If not, left-rotate
        // the parent's sub-tree so the parent becomes the left child of the
        // current node (see _rotate_left).
        if (currP == currP->parentP->rightP)
        {
          currP = currP->parentP;
          _rotate_left (currP);
        }

        // Color the parent black and the grandparent red.
        currP->parentP->color = Node::BLACK;
        CGAL_multiset_assertion (grandparentP == currP->parentP->parentP);
        grandparentP->color = Node::RED;

        // Right-rotate the grandparent's sub-tree
        _rotate_right (grandparentP);
      }
    }
    else
    {
      // If the red parent is a right child, the uncle is the left child of
      // the grandparent.
      uncleP = grandparentP->leftP;

      if (_is_red (uncleP))
      {
        // If both parent and uncle are red, color them black and color the
        // grandparent red.
        // In case of a nullptr uncle, we treat it as a black node.
        currP->parentP->color = Node::BLACK;
        uncleP->color = Node::BLACK;
        grandparentP->color = Node::RED;

        // Move to the grandparent.
        currP = grandparentP;
      }
      else
      {
        // Make sure the current node is a right child. If not, right-rotate
        // the parent's sub-tree so the parent becomes the right child of the
        // current node.
        if (currP == currP->parentP->leftP)
        {
          currP = currP->parentP;
          _rotate_right (currP);
        }

        // Color the parent black and the grandparent red.
        currP->parentP->color = Node::BLACK;
        CGAL_multiset_assertion(grandparentP == currP->parentP->parentP);
        grandparentP->color = Node::RED;

        // Left-rotate the grandparent's sub-tree
        _rotate_left (grandparentP);
      }
    }
  }

  // Make sure that the root is black.
  if (_is_red (rootP))
  {
    // In case we color a red root black, we should increment the black
    // height of the tree.
    rootP->color = Node::BLACK;
    iBlackHeight++;
  }

  return;
}

//---------------------------------------------------------
// Fix-up the tree so it maintains the red-black properties after removal.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_remove_fixup (Node* nodeP,
                                                        Node* parentP)
{
  Node        *currP = nodeP;
  Node        *currParentP = parentP;
  Node        *siblingP;

  while (currP != rootP && _is_black (currP))
  {
    // Get a pointer to the current node's sibling (notice that the node's
    // parent must exist, since the node is not the rootP).
    if (currP == currParentP->leftP)
    {
      // If the current node is a left child, its sibling is the right
      // child of the parent.
      siblingP = currParentP->rightP;

      // Check the sibling's color. Notice that nullptr nodes are treated
      // as if they are colored black.
      if (_is_red (siblingP))
      {
        // In case the sibling is red, color it black and rotate.
        // Then color the parent red (and the grandparent is now black).
        siblingP->color = Node::BLACK;
        currParentP->color = Node::RED;
        _rotate_left (currParentP);
        siblingP = currParentP->rightP;
      }

      CGAL_multiset_assertion (_is_valid (siblingP));

      if (_is_black (siblingP->leftP) && _is_black (siblingP->rightP))
      {
        // If the sibling has two black children, color it red.
        siblingP->color = Node::RED;

        // The black-height of the entire sub-tree rooted at the parent is
        // now too small - fix it up recursively.
        currP = currParentP;
        currParentP = currParentP->parentP;

        // In case the current node is the tree root, we have just decreased
        // the black height of the entire tree.
        if (currP == rootP)
        {
          CGAL_multiset_assertion (currParentP == nullptr);
          iBlackHeight--;
        }
      }
      else
      {
        // In this case, at least one of the sibling's children is red.
        // It is therfore obvious that the sibling itself is black.
        if (_is_black (siblingP->rightP))
        {
          // The left child is red: Color it black, and color the sibling red.
          siblingP->leftP->color = Node::BLACK;
          siblingP->color = Node::RED;

          _rotate_right (siblingP);
          siblingP = currParentP->rightP;
        }

        // Color the parent black (it is now safe to color the sibling with
        // the same color the parent used to have) and rotate left around it.
        siblingP->color = currParentP->color;
        currParentP->color = Node::BLACK;
        if (_is_valid (siblingP->rightP))
          siblingP->rightP->color = Node::BLACK;
        _rotate_left (currParentP);

        // We set currP to be the root node in order to terminate the loop.
        currP = rootP;
      }
    }
    else
    {
      // If the current node is a right child, its sibling is the left
      // child of the parent.
      siblingP = currParentP->leftP;

      // Check the sibling's color. Notice that nullptr nodes are treated
      // as if they are colored black.
      if (_is_red (siblingP))
      {
        // In case the sibling is red, color it black and rotate.
        // Then color the parent red (and the grandparent is now black).
        siblingP->color = Node::BLACK;
        currParentP->color = Node::RED;
        _rotate_right (currParentP);

        siblingP = currParentP->leftP;
      }

      CGAL_multiset_assertion (_is_valid (siblingP));

      if (_is_black (siblingP->leftP) && _is_black (siblingP->rightP))
      {
        // If the sibling has two black children, color it red.
        siblingP->color = Node::RED;

        // The black-height of the entire sub-tree rooted at the parent is
        // now too small - fix it up recursively.
        currP = currParentP;
        currParentP = currParentP->parentP;

        // In case the current node is the tree root, we have just decreased
        // the black height of the entire tree.
        if (currP == rootP)
        {
          CGAL_multiset_assertion (currParentP == nullptr);
          iBlackHeight--;
        }
      }
      else
      {
        // In this case, at least one of the sibling's children is red.
        // It is therfore obvious that the sibling itself is black.
        if (_is_black (siblingP->leftP))
        {
          // The right child is red: Color it black, and color the sibling red.
          siblingP->rightP->color = Node::BLACK;
          siblingP->color = Node::RED;

          _rotate_left (siblingP);
          siblingP = currParentP->leftP;
        }

        // Color the parent black (it is now safe to color the sibling with
        // the same color the parent used to have) and rotate right around it.
        siblingP->color = currParentP->color;
        currParentP->color = Node::BLACK;
        if (_is_valid (siblingP->leftP))
          siblingP->leftP->color = Node::BLACK;
        _rotate_right (currParentP);

        // We set currP to be the root node in order to terminate the loop.
        currP = rootP;
      }
    }
  }

  // Make sure the current node is black.
  if (_is_red (currP))
  {
    currP->color = Node::BLACK;

    if (currP == rootP)
    {
      // In case we color a red root black, we should increment the black
      // height of the tree.
      iBlackHeight++;
    }
  }

  return;
}

//---------------------------------------------------------
// Allocate and initialize new tree node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
typename Multiset<Type, Compare, Allocator, UseCompactContainer>::Node*
Multiset<Type, Compare, Allocator, UseCompactContainer>::_allocate_node
        (const Type& object,
         typename Node::Node_color color)
{
  CGAL_multiset_assertion (color != Node::DUMMY_BEGIN &&
                           color != Node::DUMMY_END);

  Node* new_node = node_alloc.allocate(beginNode);
  new_node->init(object, color);
  return (new_node);
}

//---------------------------------------------------------
// De-allocate a tree node.
//
template <class Type, class Compare, typename Allocator, typename UseCompactContainer>
void Multiset<Type, Compare, Allocator, UseCompactContainer>::_deallocate_node (Node* nodeP)
{
  node_alloc.deallocate (nodeP);
}

} //namespace CGAL

#endif
