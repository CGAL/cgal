// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>,


#ifndef RED_BLACK_TREE_H
#define RED_BLACK_TREE_H


#include<CGAL/assertions.h>
#include<CGAL/enum.h>
#include<CGAL/memory.h>



CGAL_BEGIN_NAMESPACE

/*!
 * Container class for a red-black tree: A binary tree that satisfies the
 * following properties:
 * 1. Each node has a color, which is either red or black.
 * 2. Each red node has two red children (if a child is missing, it is 
 *    considered as a black node).
 * 3. The number of black nodes from every path from the tree root to a leaf
 *    is the same for all tree leaves (it is called the 'black depth' of the 
 *    tree).
 * Due to propeties 2-3, the depth of a red-black tree containing n nodes
 * is bounded by 2*log_2(n).
 *
 * The Red_black_tree template requires two template parmeters:
 * - The contained TYPE class represents the objects stored in the tree.
 *   It has to support the copy constructor and the assignment operator 
 *   (operator=).
 * - COMP is a functor used to define the order of objects of class TYPE:
 *   This class has to support an operator() that recieves two objects from
 *   the TYPE class and returns a negative, zero or a positive integer,
 *   depending on the comparison result.
 */

template <class TYPE, class COMP , typename Alloc = CGAL_ALLOCATOR(int)> 
class Red_black_tree
{

  struct Node;

  typedef typename Alloc::template rebind <Node> Node_alloc_rebind;
  typedef typename Node_alloc_rebind::other        Node_allocator;

  public:

  /*!
   * Representation of a node in a red-black tree.
   */
  struct Node
  {

    enum Color
    {
      Red,
      Black
    };

    typedef char Node_color;

    TYPE        object;             // The stored object.
    Node_color  color;              // The color of the node.
    Node        *parentP;           // Points on the parent node.
    Node        *rightP;            // Points on the right child of the node.
    Node        *leftP;             // Points on the left child of the node.


    Node() : 
      parentP(NULL),
      rightP(NULL),
      leftP(NULL)
    {}

    /*!
     * Constructor of a red-black tree node.
     * \param _object The object stored in the node.
     * \param _color The color of the node.
     */
    Node (const TYPE& _object,  Node_color _color) :
      object(_object),
      color(_color),
      parentP(NULL),
      rightP(NULL),
      leftP(NULL)
    {}

    void init(const TYPE& _object,  Node_color _color)
    {
      object  =  _object;
      color   =  _color;
     /* parentP =  NULL;
      rightP  =  NULL;
      leftP   =  NULL;*/
    }


    /*!
     * Recursive destructor for the entire sub-tree.
     */
    ~Node ()
    {}  

    void clear( Node_allocator& node_alloc)
    {
    // Delete the right sub-tree recursively. 
    if (rightP != NULL)
    {
      //delete rightP;
      rightP->clear(node_alloc);
      node_alloc.destroy(rightP); 
      node_alloc.deallocate(rightP, 1);
    }
    rightP = NULL;

    // Delete the left sub-tree recursively. 
    if (leftP != NULL)
    {
      //delete leftP;
      leftP->clear(node_alloc);
      node_alloc.destroy(leftP); 
      node_alloc.deallocate(leftP, 1);
    }
    leftP = NULL;
  }  

  };

  
  
  

public:

  /*!
   * Handles are used to identify objects in the tree.
   */
  class Handle
  {
    // Give the red-black tree class template access to the Handle's members.

    friend class Red_black_tree<TYPE,COMP,Alloc>;



  private:

    Node                 *nodeP;       // Points to a node in the tree.

    /*!
     * Private constructor.
     */
  

    Handle (Node *_nodeP) :
      nodeP(_nodeP)
    {}

  public:

    /*!
     * Deafult constructor.
     */
    Handle () :
      nodeP(NULL)
    {}

    /*!
     * Constructor of a NULL handle.
     */
    Handle (const void* p) :
      nodeP(reinterpret_cast<Node*>(const_cast<void*>(p)))
    {}

    /*!
     * Equality operators.
     */
    bool operator== (const void* p) const
    {
      return (nodeP == p);
    }

    bool operator!= (const void* p) const
    {
      return (nodeP != p);
    }


    Handle& operator++()
    {
      nodeP = Red_black_tree<TYPE,COMP,Alloc>::_successor(nodeP);
      return *this;
    }

    Handle operator++(int)
    {
      Handle temp = *this;
      ++*this;
      return temp;
    }

    Handle& operator--()
    {
      nodeP = Red_black_tree<TYPE,COMP,Alloc>::_predecessor(nodeP);
      return *this;
    }

    Handle operator--(int)
    {
      Handle temp = *this;
      --*this;
      return temp;
    }






    /*!
     * Get the object pointed by the handle.
     */
    const TYPE& operator* () const
    {
      CGAL_precondition(nodeP != NULL);

      return (nodeP->object);
    }

    TYPE& operator* () 
    {
      CGAL_precondition(nodeP != NULL);

      return (nodeP->object);
    }

  };

  friend class Handle;
  typedef Handle  iterator;
protected:

  // Data members:
  Node            *rootP;         // Pointer to the tree root.
  Node            *leftmostP;     // Points to the leftmost node (the minimum).
  Node            *rightmostP;   // Points to the righttmost node (the maximum).
  unsigned int    iSize;          // Number of objects stored in the tree.
  COMP            *compP;         // Used to compare the TYPE objects.
  bool            own_comp;       // Does the tree own its COMP object.
  Node_allocator  m_node_allocator; // allocator for the nodes
  Node            m_master_node;    //a default node object ,used by allocator

public:

  /*!
   * Default constructor. [takes O(1) operations]
   */
  Red_black_tree ();

  /*!
   * Constructor with a comparison object. [takes O(1) operations]
   * \param _compP A pointer to the comparison object to be used by the tree.
   */
  Red_black_tree (COMP* _compP);
  
  /*!
   * Copy constructor. [takes O(n) operations]
   * \param tree The copied tree.
   */
  Red_black_tree (const Red_black_tree<TYPE, COMP, Alloc>& tree);

  /*!
   * Destructor. [takes O(n) operations]
   */
  virtual ~Red_black_tree ();

  /*!
   * Assignment operator. [takes O(n) operations]
   * \param tree The copied tree.
   */
  const Red_black_tree<TYPE, COMP, Alloc>& 
  operator= (const Red_black_tree<TYPE, COMP, Alloc>& tree);

  /*!
   * Get the size of the tree. [takes O(1) operations]
   * \return The number of objects stored in the tree.
   */
  unsigned int size () const;

  /*!
   * Get the depth of the tree. [takes O(n) operations]
   * \return The length of the longest path from the root to a leaf node.
   */
  int depth () const;

  /*!
   * Check whether the tree contains an object. [takes O(log n) operations]
   * \param object The query object.
   * \return (true) if an equal object is found in the tree, otherwise (false).
   */
  bool contains (const TYPE& object) const;

  /*!
   * Get a handle to the given object. [takes O(log n) operations]
   * \param object The desired object.
   * \return A handle to the given object,
   *         or a NULL handle if no such object is found in the tree.
   */
  Handle get (const TYPE& object) const;

  /*!
   * Insert an object to the tree. [takes O(log n) operations]
   * \param object The object to be inserted.
   * \return A handle to the inserted object.
   */
  Handle insert (const TYPE& object);

  /*!
   * Insert an object to the tree, as the successor the given object.
   * [takes O(log n) operations at worst-case, but only O(1) amortized]
   * \param object The object to be inserted.
   * \param handle The handle after which the object should be inserted,
   *                or NULL to insert the object as the tree minimum.
   * \pre The handle is NULL, or is a valid handle created by the tree.
   *      The operation does not violate the tree properties.
   * \return A handle to the inserted object.
   */
  Handle insert_successor (const Handle& handle,
                           const TYPE& object);

  /*!
   * Insert an object to the tree, as the predecessor the given object.
   * [takes O(log n) operations at worst-case, but only O(1) amortized]
   * \param object The object to be inserted.
   * \param handle The handle before which the object should be inserted,
   *               or NULL to insert the object as the tree maximum.
   * \pre The handle is NULL, or is a valid handle created by the tree.
   *      The operation does not violate the tree properties.
   * \return A handle to the inserted object.
   */
  Handle insert_predecessor (const Handle& handle,
                             const TYPE& object);

  /*!
   * Remove an object from the tree. [takes O(log n) operations]
   * \param object The object to be removed.
   * \pre The object should be contained in the tree.
   */
  void remove (const TYPE& object);

  /*!
   * Remove the object pointed by the given handle. 
   * [takes O(log n) operations at worst-case, but only O(1) amortized]
   * \param handle A handle pointing the object to be removed from the tree.
   * \pre The handle must be a valid handle created by the tree.
   *      Note that the handles to removed objects become invalid.
   */
  void remove_at (const Handle& handle);

  /*!
   * Replace the object given by its handle by another object.
   * [takes O(1) operations]
   * \param handle A handle pointing the object to be replaced.
   * \param object The new object.
   * \pre The operation does not violate the tree properties.
   */
  void replace (const Handle& handle,
                const TYPE& object);
        
  /*!
   * Reset the binary tree. [takes O(n) operations]
   */
  void reset ();

  /*!
   * Get a handle to the tree minimum. [takes O(log n) operations]
   * \return A handle to the minimal object in the tree, 
   *         or a NULL handle if the tree is empty.
   */
  Handle minimum () const;

  /*!
   *  Get a handle to the tree maximum. [takes O(log n) operations]
   * \return A handle to the maximal object in the tree, 
   *         or a NULL handle if the tree is empty.
   */
  Handle maximum () const;

  /*!
   * Get a handle to the next object in the tree (according to the tree order).
   * [takes O(log n) operations at worst-case, but only O(1) amortized]
   * \param handle A handle to the current object.
   * \pre The handle must be a valid handle created by the tree.
   * \return The successor node, 
   *         or a NULL hanlde if we are at the tree maximum.
   */
  Handle successor (const Handle& handle) const;

  /*! 
   * Get a handle to the previous node in the tree (according to the tree 
   * order). [takes O(log n) operations at worst-case, but only O(1) amortized]
   * \param handle A handle to the current object.
   * \pre The handle must be a valid handle created by the tree.
   * \return The predecessor node, 
   *         or a NULL hanlde if we are at the tree minimum.
   */
  Handle predecessor (const Handle& handle) const;

  /*!
   * Get a handle to some node in the tree and check if its the minimum.
   * [takes o(logn) operations at worst-case , but only o(1) amortized]
   * \param handle A handle to some node.
   * \return true iff handle is the minimum.
   */
  bool is_minimum(const Handle& handle) const;


  /*!
   * Check validation of the tree. [takes o(n) operations]
   * \return true iff the tree is valid.
   *
   */
  bool is_valid() const;

  /*!
   * Get the first element whose key is not less than object.
   * [takes o(log(n)) operations.
   * \param object the object to be looking for lower bound.
   * \return The lower bound of object 
   *         or a NULL handle if there isn't one.
   */
  Handle lower_bound (const TYPE& object) const;

protected:

  /*!
   * Get a node that contains the given object.
   * \param object The desired object.
   * \return A node that contains a given object,
   *         or NULL if no such object is found in the tree.
   */
  Node* _get (const TYPE& object) const;

  /*!
   * Remove the object stored in the given tree node. 
   * \param nodeP The node storing the object to be removed from the tree.
   */
  void _remove_at (Node* nodeP);

  /*!
   * Get the next node in the tree (according to the tree order).
   * \param nodeP The current node.
   * \return The successor node, or NULL if nodeP is the tree maximum.
   */
  static Node* _successor (Node* nodeP) /*const*/;

  /*! 
   * Get the previous node in the tree (according to the tree order).
   * \param nodeP The current node.
   * \return The predecessor node, or NULL if nodeP is the tree minimum.
   */
  static Node* _predecessor (Node* nodeP) ;

  /*!
   * Calculate the depth of the sub-tree spanned by the given node.
   * \param nodeP The sub-tree root.
   * \return The sub-tree depth.
   */
  int _sub_depth (const Node* nodeP) const;

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
  Node* _duplicate (const Node* nodeP) const;

  /*!
   * Fix-up the red-black tree properties after an insertion operation.
   * \param nodeP The node that has just been inserted to the tree.
   */
  void _insert_fixup (Node* nodeP);

  /*!
   * Fix-up the red-black tree properties after a removal operation.
   * \param nodeP The child of the node that has just been removed from 
   *              the tree.
   */
  void _remove_fixup (Node* nodeP);

  
  Node* _lower_bound (const TYPE& object) const;

  Node* _allocate_node(const TYPE& object, typename Node::Node_color color);

};

//---------------------------------------------------------
// Default constructor.
//
template <class TYPE, class COMP, typename Alloc>
Red_black_tree<TYPE, COMP, Alloc>::Red_black_tree () :
  rootP(NULL),
  leftmostP(NULL),
  rightmostP(NULL),
  iSize(0),
  compP(NULL),
  own_comp(true)
{
  compP = new COMP;
}

//---------------------------------------------------------
// Constructor with a comparison object.
//
template <class TYPE, class COMP, typename Alloc>
Red_black_tree<TYPE, COMP, Alloc>::Red_black_tree (COMP* _compP) :
  rootP(NULL),
  leftmostP(NULL),
  rightmostP(NULL),
  iSize(0),
  compP(_compP),
  own_comp(false)
{
  CGAL_precondition(_compP != NULL);
}

//---------------------------------------------------------
// Copy constructor.
//
template <class TYPE, class COMP, typename Alloc>
Red_black_tree<TYPE, COMP, Alloc>::Red_black_tree 
        (const Red_black_tree<TYPE, COMP, Alloc>& tree) :
  rootP(NULL),
  leftmostP(NULL),
  rightmostP(NULL),
  iSize(tree.iSize),
  compP(NULL),
  own_comp(tree.own_comp)
{
  // If necessary, create a COMP object.
  if (own_comp)
    compP = new COMP;
  else
    compP = tree.compP;

  // Copy all the copied tree's nodes recursively.
  if (tree.rootP != NULL)
    rootP = _duplicate (tree.rootP);

  leftmostP = _sub_minimum (rootP);
  rightmostP = _sub_maximum (rootP);
}

//---------------------------------------------------------
// Destructor.
//
template <class TYPE, class COMP, typename Alloc>
Red_black_tree<TYPE, COMP, Alloc>::~Red_black_tree ()
{
  // If necessary, free the COMP object.
  if (own_comp)
    delete compP;
  compP = NULL;

  // Since the node destructor is recursive, it will delete the entire tree.
  if (rootP != NULL)
  {
    rootP->clear(m_node_allocator);
    m_node_allocator.destroy(rootP); 
    m_node_allocator.deallocate(rootP, 1);
  }
  rootP = NULL;
  leftmostP = NULL;
  rightmostP = NULL;
}

//---------------------------------------------------------
// Assignment operator.
//
template <class TYPE, class COMP, typename Alloc>
const Red_black_tree<TYPE, COMP, Alloc>& Red_black_tree<TYPE, COMP, Alloc>::operator=
        (const Red_black_tree<TYPE, COMP, Alloc>& tree)
{
  // Avoid self-assignment.
  if (this == &tree)
    return (*this);

  // Deal with the COMP object.
  if (own_comp)
    delete compP;

  own_comp = tree.own_comp;
  if (own_comp)
    compP = new COMP;
  else
    compP = tree.compP;

  // Remove all objects currently stored in the tree.
  reset();

  // Copy all the copied tree's nodes recursively.
  if (tree.rootP != NULL)
    rootP = _duplicate (tree.rootP);

  // Update the number of objects stored in the tree.
  iSize = tree.iSize;

  // Update the leftmost and rightmost pointers.
  leftmostP = _sub_minimum (rootP);
  rightmostP = _sub_maximum (rootP);

  return (*this);
}

//---------------------------------------------------------
// Get the number of objects stored in the tree.
//
template <class TYPE, class COMP, typename Alloc>
unsigned int Red_black_tree<TYPE, COMP, Alloc>::size () const
{
  return (iSize);
}

//---------------------------------------------------------
// Get the depth of the tree.
//
template <class TYPE, class COMP, typename Alloc>
int Red_black_tree<TYPE, COMP, Alloc>::depth () const
{
  if (rootP == NULL)
  {
    // Empty tree.
    return (0);
  }
  else
  {
    // Return the depth of the root's sub-tree (the entire tree).
    return (_sub_depth(rootP));
  }
}

//---------------------------------------------------------
// Check whether the tree contains the given object.
//
template <class TYPE, class COMP, typename Alloc>
bool Red_black_tree<TYPE, COMP, Alloc>::contains (const TYPE& object) const
{
    return (_get(object) != NULL);
}

//---------------------------------------------------------
// Return a pointer to the node containing the given object.
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP,Alloc>::Handle 
        Red_black_tree<TYPE, COMP, Alloc>::get (const TYPE& object) const
{
  return (Handle(_get(object)));
}

//---------------------------------------------------------
// Insert a new object to the tree.
//
template <class TYPE, class COMP, typename Alloc>
 typename Red_black_tree<TYPE, COMP, Alloc>::Handle
        Red_black_tree<TYPE, COMP, Alloc>::insert (const TYPE& object)
{
  if (rootP == NULL)
  {
    // In case the tree is empty, assign a new rootP.
    // Notice that the root is always black.
    //rootP = new Node (object, Node::Black);
    /*rootP = m_node_allocator.allocate(1);
    m_node_allocator.construct(rootP, Node());
    rootP -> init(object, Node::Black);*/
    rootP = _allocate_node(object, Node::Black);
    
    iSize = 1;

    // As the tree now contains a single node:
    leftmostP = rootP;
    rightmostP = rootP;

    return (Handle(rootP));
  }

  // Find a place for the new object, and insert it as a red leaf.    
  Node                      *currentP = rootP;
  //Node        *newNodeP = new Node (object, Node::Red);
  /*Node                      *newNodeP = m_node_allocator.allocate(1);
  m_node_allocator.construct(newNodeP, Node());
  newNodeP -> init(object, Node::Red);*/
  Node*  newNodeP = _allocate_node(object, Node::Red);

  Comparison_result         iCompResult;
  bool                      is_leftmost = true;
  bool                      is_rightmost = true;

  while (currentP != NULL)
  {
    // Compare the inserted object with the object stored in the current node.
    //if ((*compP)(object, currentP->object) > 0) // SMALLER(obj ,currentP->object)
    iCompResult = (*compP)(object, currentP->object);

    if(iCompResult == SMALLER)
    {
      is_rightmost = false;

      if (currentP->leftP == NULL)
      {
        // Insert the new leaf as the left child of the current node.
        currentP->leftP = newNodeP;
        newNodeP->parentP = currentP;
        currentP = NULL;            // In order to terminate the while loop.

        if (is_leftmost)
          leftmostP = newNodeP;
      }
      else
      {
        // Go to the left sub-tree.
        currentP = currentP->leftP;
      }
    }
    else if(iCompResult == LARGER)
    {
      is_leftmost = false;

      if (currentP->rightP == NULL)
      {
        // Insert the new leaf as the right child of the current node.
        currentP->rightP = newNodeP;
        newNodeP->parentP = currentP;
        currentP = NULL;            // In order to terminate the while loop.

        if (is_rightmost)
          rightmostP = newNodeP;
      }
      else
      {
        // Go to the right sub-tree.
        currentP = currentP->rightP;
      }
    }
    else // iCompResult == EQUAL
    {
      return (Handle(currentP));
    }
  }

  // Mark that a new node was added.
  iSize++;

  // Fix up the tree properties.
  _insert_fixup (newNodeP);  

  return (Handle(newNodeP));
}

//---------------------------------------------------------
// Insert a new object to the tree as the a successor of a given node.
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Handle
        Red_black_tree<TYPE, COMP, Alloc>::insert_successor (const Handle& handle,
                                                      const TYPE& object)
{
  Node  *nodeP = handle.nodeP;

  CGAL_assertion(nodeP == NULL || handle.treeP == this);

  if (rootP == NULL)
  {
    // In case the tree is empty, make sure that we did not recieve an invalid
    // handle.
    CGAL_precondition(nodeP == NULL);

    // Assign a new root node. Notice that the root is always black.
    //rootP = new Node (object, Node::Black);
    /*rootP = m_node_allocator.allocate(1);
    m_node_allocator.construct(rootP, Node());
    rootP -> init(object, Node::Black);*/
    rootP = _allocate_node(object, Node::Black);

    iSize = 1;

    // As the tree now contains a single node:
    leftmostP = rootP;
    rightmostP = rootP;

    return (Handle(rootP));
  }

  // Insert the new object as a red leaf, being the successor of nodeP.
  Node        *parentP;
  //Node        *newNodeP = new Node (object, Node::Red);
  /*Node *newNodeP = m_node_allocator.allocate(1);
  m_node_allocator.construct(newNodeP, Node());
  newNodeP -> init(object, Node::Red);*/
  Node *newNodeP = _allocate_node(object, Node::Red);

  if (nodeP == NULL)
  {
    // The new node should become the tree minimum: Place is as the left
    // child of the current minimal leaf.
    parentP = _sub_minimum(rootP);

    //CGAL_precondition((*compP)(object, parentP->object) >= 0);
    CGAL_precondition((*compP)(object, parentP->object) != LARGER);

    parentP->leftP = newNodeP;

    // As we inserted a new tree minimum:
    leftmostP = newNodeP;
  }
  else
  {
    // Make sure the insertion does not violate the tree order.
    CGAL_precondition_code(Node *_succP = _successor(nodeP));
    //CGAL_precondition((*compP)(nodeP->object, object) >= 0);
    CGAL_precondition((*compP)(nodeP->object, object) != LARGER);
    //CGAL_precondition(_succP == NULL || (*compP)(object, _succP->object) >= 0);
    CGAL_precondition(_succP == NULL || (*compP)(object, _succP->object) != LARGER);

    // In case given node has no right child, place the new node as its 
    // right child. Otherwise, place it at the leftmost position at the
    // sub-tree rooted at its right side.
    if (nodeP->rightP == NULL)
    {
      parentP = nodeP;
      parentP->rightP = newNodeP;
    }
    else
    {
      parentP = _sub_minimum (nodeP->rightP);
      parentP->leftP = newNodeP;
    }

    if (nodeP == rightmostP)
      // As we inserted a new tree maximum:
      rightmostP = newNodeP;
  }

  newNodeP->parentP = parentP;

  // Mark that a new node was added.
  iSize++;

  // Fix up the tree properties.
  _insert_fixup (newNodeP);  

  return (Handle (newNodeP));
}

//---------------------------------------------------------
// Insert a new object to the tree as the a predecessor of a given node.
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Handle
        Red_black_tree<TYPE, COMP, Alloc>::insert_predecessor (const Handle& handle,
                                                        const TYPE& object)
{
  Node  *nodeP = handle.nodeP;

  CGAL_assertion(nodeP == NULL || handle.treeP == this);

  if (rootP == NULL)
  {
    // In case the tree is empty, make sure that we did not recieve an invalid
    // handle.
    CGAL_precondition(nodeP == NULL);

    // Assign a new root node. Notice that the root is always black.
    //rootP = new Node (object, Node::Black);
    /*rootP = m_node_allocator.allocate(1);
    m_node_allocator.construct(rootP, Node());
    rootP -> init(object, Node::Black);*/
    rootP = _allocate_node(object, Node::Black);

    iSize = 1;

    // As the tree now contains a single node:
    leftmostP = rootP;
    rightmostP = rootP;

    return (Handle(rootP));
  }

  // Insert the new object as a red leaf, being the predecessor of nodeP.
  Node        *parentP;
  //Node        *newNodeP = new Node (object, Node::Red);
  /*Node *newNodeP = m_node_allocator.allocate(1);
  m_node_allocator.construct(newNodeP, Node());
  newNodeP -> init(object, Node::Red);*/
  Node *newNodeP = _allocate_node(object, Node::Red);

  if (nodeP == NULL)
  {
    // The new node should become the tree maximum: Place is as the right
    // child of the current maximal leaf.
    parentP = _sub_maximum(rootP);

   // CGAL_precondition((*compP)(object, parentP->object) <= 0);
    CGAL_precondition((*compP)(object, parentP->object) != SMALLER);

    parentP->rightP = newNodeP;

    // As we inserted a new tree maximum:
    rightmostP = newNodeP;
  }
  else
  {
    // Make sure the insertion does not violate the tree order.
    CGAL_precondition_code(Node *_predP = _predecessor(nodeP));
    //CGAL_precondition((*compP)(nodeP->object, object) <= 0);
    CGAL_precondition((*compP)(nodeP->object, object) != SMALLER);
    //CGAL_precondition(_predP == NULL || (*compP)(object, _predP->object) <= 0);
    CGAL_precondition(_predP == NULL || (*compP)(object, _predP->object) != SMALLER);

    // In case given node has no left child, place the new node as its 
    // left child. Otherwise, place it at the rightmost position at the
    // sub-tree rooted at its left side.
    if (nodeP->leftP == NULL)
    {
      parentP = nodeP;
      parentP->leftP = newNodeP;
    }
    else
    {
      parentP = _sub_maximum (nodeP->leftP);
      parentP->rightP = newNodeP;
    }

    if (nodeP == leftmostP)
      // As we inserted a new tree minimum:
      leftmostP = newNodeP;
  }

  newNodeP->parentP = parentP;

  // Mark that a new node was added.
  iSize++;

  // Fix up the tree properties.
  _insert_fixup (newNodeP);  

  return (Handle(newNodeP));
}

//---------------------------------------------------------
// Remove an object from the tree.
//
template <class TYPE, class COMP, typename Alloc>
void Red_black_tree<TYPE, COMP, Alloc>::remove (const TYPE& object)
{
  // Find the node containing the object that should be removed.
  Node        *nodeP = _get(object);
  
  CGAL_precondition (nodeP != NULL);
  
  // Remove the node.
  _remove_at (nodeP);
  return;
}

//---------------------------------------------------------
// Remove the object pointed by the given handle.
//
template <class TYPE, class COMP, typename Alloc>
void Red_black_tree<TYPE, COMP, Alloc>::remove_at (const Handle& handle)
{
  CGAL_precondition (handle.treeP == this);

  _remove_at (handle.nodeP);
  return;
}

//---------------------------------------------------------
// Replace the object pointed by the given handle.
//
template <class TYPE, class COMP, typename Alloc>
void Red_black_tree<TYPE, COMP, Alloc>::replace (const Handle& handle,
                                                 const TYPE& object)
{
  CGAL_precondition(handle.treeP == this);

  Node  *nodeP = handle.nodeP;
  CGAL_precondition(nodeP != NULL);

  // Make sure the replacement does not violate the tree order.
  CGAL_precondition_code(Node *_succP = _successor(nodeP));
  //CGAL_precondition(_succP == NULL || (*compP)(object, _succP->object) >= 0);
  CGAL_precondition(_succP == NULL || (*compP)(object, _succP->object) != LARGER);

  CGAL_precondition_code(Node *_predP = _predecessor(nodeP));
  //CGAL_precondition(_predP == NULL || (*compP)(object, _predP->object) <= 0);
  CGAL_precondition(_predP == NULL || (*compP)(object, _predP->object) != SMALLER);
  
  // Replace the object at nodeP.
  nodeP->object = object;

  return;
}
        
//---------------------------------------------------------
// Remove all objects from the tree.
//
template <class TYPE, class COMP, typename Alloc>
void Red_black_tree<TYPE, COMP, Alloc>::reset ()
{
  // Delete all the tree nodes recursively.
  if (rootP != NULL)
  {
   // delete rootP;
    rootP->clear(m_node_allocator);
    m_node_allocator.destroy(rootP); 
    m_node_allocator.deallocate(rootP, 1);
  }
  rootP = NULL;
  leftmostP = NULL;
  rightmostP = NULL;

  // Mark that there are no more objects in the tree.
  iSize = 0;

  return;
}

//---------------------------------------------------------
// Get a handle to the tree minimum.
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Handle 
        Red_black_tree<TYPE, COMP, Alloc>::minimum () const
{
  // Return the leftmost leaf in the tree (or NULL if the tree is empty).
  return (Handle(leftmostP));
}

//---------------------------------------------------------
// Get a handle the tree maximum.
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Handle
        Red_black_tree<TYPE, COMP, Alloc>::maximum () const
{
  // Return the rightmost leaf in the tree (or NULL if the tree is empty).
  return (Handle(rightmostP));
}

//---------------------------------------------------------
// Get a handle to the next node in the tree (according to the tree order).
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Handle
        Red_black_tree<TYPE, COMP, Alloc>::successor (const Handle& handle) const
{
  CGAL_precondition(handle.treeP == this);
  CGAL_precondition(handle.nodeP != NULL);

  // Get the successor.
  return (Handle(_successor(handle.nodeP)));
}

//---------------------------------------------------------
// Get a handle to the previous node in the tree (according to the tree order).
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Handle 
        Red_black_tree<TYPE, COMP, Alloc>::predecessor (const Handle& handle) const
{
  CGAL_precondition(handle.treeP == this);
  CGAL_precondition(handle.nodeP != NULL);

  // Get the predecessor.
  return (Handle(_predecessor(handle.nodeP)));
}

//---------------------------------------------------------
// Check if the tree is a valid one.
//
template <class TYPE, class COMP, typename Alloc>
bool Red_black_tree<TYPE, COMP, Alloc>::is_valid() const
{
  Node        *currentP = minimum().nodeP;
  Node        *prevP = currentP;
  
  while((currentP = _successor(currentP)) != NULL)
  {
    if(((*compP)(prevP->object , currentP->object)) == LARGER)
      return false;
    prevP = currentP;
  }
  return true;
}


//---------------------------------------------------------
// Find the lower bound of object
//

template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Handle  
        Red_black_tree<TYPE, COMP, Alloc>::lower_bound (const TYPE& object) const
{
  return (Handle(_lower_bound(object)));
}


//---------------------------------------------------------
// Check if handle is the minimum
//---------------------------------------------------------
template <class TYPE, class COMP, typename Alloc>
bool Red_black_tree<TYPE, COMP, Alloc>::is_minimum(const Handle& handle) const
{
  //TODO : add a precondition handle!=NULL
  if(_predecessor(handle.nodeP) == NULL)
    return true;
  return false;
}

//---------------------------------------------------------
// Return a pointer to the node containing the given object.
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Node* 
        Red_black_tree<TYPE, COMP, Alloc>::_get (const TYPE& object) const
{
  Node                      *currentP = rootP;
  Comparison_result         iCompResult;

  while (currentP != NULL)
  {
    //if ((iCompResult = (*compP)(object, currentP->object)) == 0)
    if ((iCompResult = (*compP)(object, currentP->object)) == EQUAL)
    {
      // In case of equality, we can return the current node.
      return (currentP);
    }
    else if (iCompResult == SMALLER)
    {
      // Go down to the left child.
      currentP = currentP->leftP;
    }
    else // iCompResult == LARGER
    {
      // Go down to the right child.
      currentP = currentP->rightP;
    }
  }

  // If we reached here, the object is not found in the tree.
  return (NULL);
}

//---------------------------------------------------------
// Return the pointer to the node of the lower bound of object
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Node*  
        Red_black_tree<TYPE, COMP, Alloc>::_lower_bound (const TYPE& object) const
{
  Node                      *currentP = rootP;
  Node                      *prevP = currentP;
  Comparison_result         iCompResult;

  if(rootP == NULL)
    return NULL;

  while (currentP != NULL)
  {
    if ((iCompResult = (*compP)(object, currentP->object)) == EQUAL)
    {
      // In case of equality, we can return the current node.
      return (currentP);
    }
    else if (iCompResult == SMALLER)
    {
      prevP = currentP;

      // Go down to the left child.
      currentP = currentP->leftP;
    }
    else // iCompResult == LARGER
    {
      prevP = currentP;

      // Go down to the right child.
      currentP = currentP->rightP;
    }
  }

  // If we reached here, the object is not found in the tree.
  if(iCompResult == SMALLER)
    return prevP;
  return (_successor(prevP));
}


//---------------------------------------------------------
// Remove the given tree node.
//
template <class TYPE, class COMP, typename Alloc>
void Red_black_tree<TYPE, COMP, Alloc>::_remove_at (Node* nodeP)
{
  CGAL_precondition (nodeP != NULL);

  // In case of deleting the single object stored in the tree, free the root,
  // thus emptying the tree.
  if (iSize == 1)
  {
    CGAL_assertion (nodeP == rootP);

    //delete rootP;
    rootP->clear(m_node_allocator);
    m_node_allocator.destroy(rootP); 
    m_node_allocator.deallocate(rootP, 1);
    rootP = NULL;
    leftmostP = NULL;
    rightmostP = NULL;

    iSize = 0;
    return;
  }

  // In case we delete the tree minimum of maximum, update the relevant pointers.
  if (nodeP == leftmostP)
    leftmostP = _successor (nodeP);
  else if (nodeP == rightmostP)
    rightmostP = _predecessor (nodeP);

  // Remove the given node from the tree.
  if (nodeP->leftP != NULL && nodeP->rightP != NULL)
  {
    // If the node we want to remove has two children, find its successor,
    // which is the leftmost child in its right sub-tree and has at most
    // one child (it may have a right child).
    Node    *succP = _sub_minimum (nodeP->rightP);

    // Now physically swap nodeP and its successor. Notice this may temporarily
    // violate the tree properties, but we are going to remove nodeP anyway.
    // This way we have moved nodeP to a position were it is more convinient
    // to delete it.
    bool            immediate_succ = (nodeP->rightP == succP);
    Node            *succ_parentP = succP->parentP;
    Node            *succ_leftP = succP->leftP;
    Node            *succ_rightP = succP->rightP;
    typename Node::Node_color succ_color = succP->color;

    succP->parentP = nodeP->parentP;
    succP->leftP = nodeP->leftP;
    succP->rightP = immediate_succ ? nodeP : nodeP->rightP;
    succP->color = nodeP->color;

    nodeP->parentP = immediate_succ ? succP : succ_parentP;
    nodeP->leftP = succ_leftP;
    nodeP->rightP = succ_rightP;
    nodeP->color = succ_color;

    CGAL_assertion(nodeP->parentP != NULL);
    if (! immediate_succ)
    { 
      if (succP == nodeP->parentP->leftP)
        nodeP->parentP->leftP = nodeP;
      else
        nodeP->parentP->rightP = nodeP;
    }

    if (nodeP->leftP != NULL)
      nodeP->leftP->parentP = nodeP;
    if (nodeP->rightP != NULL)
       nodeP->rightP->parentP = nodeP;

    if (succP->parentP != NULL)
    {
      if (nodeP == succP->parentP->leftP)
        succP->parentP->leftP = succP;
      else
        succP->parentP->rightP = succP;
    }
    else
    {
      rootP = succP;
    }

    if (succP->leftP != NULL)
      succP->leftP->parentP = succP;
    if (succP->rightP != NULL)
      succP->rightP->parentP = succP;
  }

  CGAL_assertion(nodeP->leftP == NULL || nodeP->rightP == NULL);

  // At this stage, the node we are going to remove has at most one child.
  Node        *childP = NULL;

  if (nodeP->leftP != NULL)
  {
    childP = nodeP->leftP;
  }
  else
  {
    childP = nodeP->rightP;
  }

  // Splice out the node to be removed, by linking its parent straight to the 
  // removed node's single child.
  if (childP != NULL)
    childP->parentP = nodeP->parentP;
    
  if (nodeP->parentP == NULL)
  {
    // If we are deleting the root, make the child the new tree node.
    rootP = childP;
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
  // just removed a black node, the black-depth property is no longer valid.
  if (nodeP->color == Node::Black && childP != NULL)
    _remove_fixup (childP);

  // Delete the un-necessary node (we nullify both its children because the 
  // node's destructor is recursive).
  nodeP->leftP = NULL;
  nodeP->rightP = NULL;
  //delete nodeP;
  nodeP->clear(m_node_allocator);
  m_node_allocator.destroy(nodeP); 
  m_node_allocator.deallocate(nodeP, 1);

  // Decrease the number of objects in the tree.
  iSize--;

  return;
}

//---------------------------------------------------------
// Get the next node in the tree (according to the tree order).
//
template <class TYPE, class COMP, typename Alloc>
 typename Red_black_tree<TYPE, COMP, Alloc>::Node* 
        Red_black_tree<TYPE, COMP, Alloc>::_successor (Node* nodeP) //const
{
  CGAL_precondition(nodeP != NULL);

  Node        *succP;

  if (nodeP->rightP != NULL)
  {
    // If there is a right child, the successor is the minimal object in 
    // the sub-tree spanned by this child.
    succP = nodeP->rightP;
    while (succP->leftP != NULL)
      succP = succP->leftP;
  }
  else
  {
    // Otherwise, go up the tree until reaching the parent from the left 
    // direction.
    const Node    *prevP = nodeP;

    succP = nodeP->parentP;
    while (succP != NULL && prevP == succP->rightP)
    {
      prevP = succP;
      succP = succP->parentP;
    }
  }

  return (succP);
}

//---------------------------------------------------------
// Get the previous node in the tree (according to the tree order).
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Node* 
        Red_black_tree<TYPE, COMP, Alloc>::_predecessor (Node* nodeP) 
{
  CGAL_precondition(nodeP != NULL);

  Node        *predP;

  if (nodeP->leftP != NULL)
  {
    // If there is a left child, the predecessor is the maximal object in 
    // the sub-tree spanned by this child.
    predP = nodeP->leftP;
    while (predP->rightP != NULL)
      predP = predP->rightP;
  }
  else
  {
    // Otherwise, go up the tree until reaching the parent from the right 
    // direction.
    const Node    *prevP = nodeP;

    predP = nodeP->parentP;
    while (predP != NULL && prevP == predP->leftP)
    {
      prevP = predP;
      predP = predP->parentP;
    }
  }

  return (predP);
}

//---------------------------------------------------------
// Calculate the depth of the subtree spanned by a given node.
//
template <class TYPE, class COMP, typename Alloc>
int Red_black_tree<TYPE, COMP, Alloc>::_sub_depth (const Node* nodeP) const
{
  // Recursively calculate the depth of the left and right sub-trees.
  int  iRightDepth = (nodeP->rightP == NULL) ? 0 : _sub_depth(nodeP->rightP);
  int  iLeftDepth = (nodeP->leftP == NULL) ? 0 : _sub_depth(nodeP->leftP);

  // Return the maximal child depth + 1 (the current node).
  return ((iRightDepth > iLeftDepth) ? (iRightDepth + 1) : (iLeftDepth + 1));
}

//---------------------------------------------------------
// Get the leftmost node in the sub-tree spanned by the given node.
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Node* 
        Red_black_tree<TYPE, COMP, Alloc>::_sub_minimum (Node* nodeP) const
{
  Node    *minP = nodeP;

  while (minP->leftP != NULL)
    minP = minP->leftP;
  return (minP);
}

//---------------------------------------------------------
// Get the rightmost node in the sub-tree spanned by the given node.
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Node* 
        Red_black_tree<TYPE, COMP, Alloc>::_sub_maximum (Node* nodeP) const
{
  Node    *maxP = nodeP;

  while (maxP->rightP != NULL)
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
template <class TYPE, class COMP, typename Alloc>
        void Red_black_tree<TYPE, COMP, Alloc>::_rotate_left (Node* xNodeP)
{
  // Get the right child of the node.
  Node        *yNodeP = xNodeP->rightP;

  // Change its left subtree (T2) to x's right subtree.
  xNodeP->rightP = yNodeP->leftP;

  // Link T2 to its new parent x.
  if (yNodeP->leftP != NULL)
    yNodeP->leftP->parentP = xNodeP;
    
  // Assign x's parent to be y's parent.
  yNodeP->parentP = xNodeP->parentP;

  if (xNodeP->parentP == NULL)
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
template <class TYPE, class COMP, typename Alloc>
        void Red_black_tree<TYPE, COMP, Alloc>::_rotate_right (Node* yNodeP)
{
  // Get the left child of the node.
  Node        *xNodeP = yNodeP->leftP;

  // Change its right subtree (T2) to y's left subtree.
  yNodeP->leftP = xNodeP->rightP;

  // Link T2 to its new parent y.
  if (xNodeP->rightP != NULL)
    xNodeP->rightP->parentP = yNodeP;
    
  // Assign y's parent to be x's parent.
  xNodeP->parentP = yNodeP->parentP;

  if (yNodeP->parentP == NULL)
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
// Return a pointer to a duplication of the given node.
//
template <class TYPE, class COMP, typename Alloc>
typename Red_black_tree<TYPE, COMP, Alloc>::Node* 
        Red_black_tree<TYPE, COMP, Alloc>::_duplicate (const Node* nodeP) const
{
  // Create a node of the same color, containing the same object.
  //Node*        dupNodeP = new Node (nodeP->object, nodeP->color);
  /*Node *dupNodeP = m_node_allocator.allocate(1);
  m_node_allocator.construct(dupNodeP, Node());
  dupNodeP -> init(nodeP->object, nodeP->color);*/
  Node *dupNodeP = _allocate_node(nodeP->object, nodeP->color);
    
  // Duplicate the children recursively.
  if (nodeP->rightP != NULL)
  {
    dupNodeP->rightP = _duplicate (nodeP->rightP);
    dupNodeP->rightP->parentP = dupNodeP;
  }
  else
  {
    dupNodeP->rightP = NULL;
  }

  if (nodeP->leftP != NULL)
  {
    dupNodeP->leftP = _duplicate (nodeP->leftP);
    dupNodeP->leftP->parentP = dupNodeP;
  }
  else
  {
    dupNodeP->leftP = NULL;
  }

  // Return the duplicated node.
  return (dupNodeP);
}

//---------------------------------------------------------
// Fix-up the tree so it maintains the red-black properties after insertion.
//
template <class TYPE, class COMP, typename Alloc>
void Red_black_tree<TYPE, COMP, Alloc>::_insert_fixup (Node* nodeP)
{
  CGAL_precondition(nodeP != NULL && nodeP->color == Node::Red);

  // Fix the red-black propreties: we may have inserted a red leaf as the 
  // child of a red parent - so we have to fix the coloring of the parent 
  // recursively.
  Node        *currP = nodeP;
  Node        *grandparentP;
  Node        *uncleP;

  while (currP != rootP && currP->parentP->color == Node::Red)
  {
    // Get a pointer to the current node's grandparent (notice the root is 
    // always black, so the red parent must have a parent).
    grandparentP = currP->parentP->parentP;
        
    if (currP->parentP == grandparentP->leftP)
    {
      // If the red parent is a left child, the uncle is the right child of 
      // the grandparent. 
      uncleP = grandparentP->rightP;

      if (uncleP != NULL && uncleP->color == Node::Red)
      {
        // If both parent and uncle are red, color them black and color the 
        // grandparent red.
        // In case of a NULL uncle, we treat it as a black node.
        currP->parentP->color = Node::Black;
        uncleP->color = Node::Black;
        grandparentP->color = Node::Red;

        // Move to the grandparent.
        currP = grandparentP;
      }
      else
      {
        // Make sure the current node is a right child. If not, left-rotate 
        // the parent's sub-tree so the parent becomes the right child of the 
        // current node (see _rotate_left).
        if (currP == currP->parentP->rightP)
        {
          currP = currP->parentP;
          _rotate_left (currP);
        }

        // Color the parent black and the grandparent red.
        currP->parentP->color = Node::Black;
        CGAL_assertion(grandparentP == currP->parentP->parentP);
        grandparentP->color = Node::Red;

        // Right-rotate the grandparent's sub-tree
        _rotate_right (grandparentP);
      }
    }
    else
    {
      // If the red parent is a right child, the uncle is the left child of 
      // the grandparent. 
      uncleP = grandparentP->leftP;

      if (uncleP != NULL && uncleP->color == Node::Red)
      {
        // If both parent and uncle are red, color them black and color the 
        // grandparent red.
        // In case of a NULL uncle, we treat it as a black node.
        currP->parentP->color = Node::Black;
        uncleP->color = Node::Black;
        grandparentP->color = Node::Red;

        // Move to the grandparent.
        currP = grandparentP;
      }
      else
      {
        // Make sure the current node is a left child. If not, right-rotate 
        // the parent's sub-tree so the parent becomes the left child of the 
        // current node.
        if (currP == currP->parentP->leftP)
        {
          currP = currP->parentP;
          _rotate_right (currP);
        }

        // Color the parent black and the grandparent red.
        currP->parentP->color = Node::Black;
        CGAL_assertion(grandparentP == currP->parentP->parentP);
        grandparentP->color = Node::Red;

        // Left-rotate the grandparent's sub-tree
        _rotate_left (grandparentP);
      }
    }
  }

  // Make sure that the root is black.
  rootP->color = Node::Black;

  return;
}

//---------------------------------------------------------
// Fix-up the tree so it maintains the red-black properties after removal.
//
template <class TYPE, class COMP, typename Alloc>
void Red_black_tree<TYPE, COMP, Alloc>::_remove_fixup (Node* nodeP)
{
  Node        *currP = nodeP;
  Node        *siblingP;

  while (currP != rootP && currP->color == Node::Black)
  {
    // Get a pointer to the current node's sibling (notice that the node's 
    // parent must exist, since the node is not the rootP).
    if (currP == currP->parentP->leftP)
    {
      // If the current node is a left child, its sibling is the right 
      // child of the parent.
      siblingP = currP->parentP->rightP;
      
      // Check the sibling's color. Notice that NULL nodes are treated
      // as if they are colored black.
      if (siblingP != NULL && siblingP->color == Node::Red)
      {
        // In case the sibling is red, color it black and rotate.
        // Then color the parent red (and the grandparent is now black).
        siblingP->color = Node::Black;
        currP->parentP->color = Node::Red;
        _rotate_left (currP->parentP);
        siblingP = currP->parentP->rightP;
      }
      
      if (siblingP != NULL &&
          (siblingP->leftP == NULL || 
           siblingP->leftP->color == Node::Black) && 
          (siblingP->rightP == NULL || 
           siblingP->rightP->color == Node::Black))
      {
        // If the sibling has two black children, color it red.
        siblingP->color = Node::Red;
        if (currP->parentP->color == Node::Red)
        {
          // If the parent is red, we can safely color it black and terminate
          // the fix-up process.
          currP->parentP->color = Node::Black;
          currP = rootP;          // In order to stop the while loop.
        }
        else
        {
          // The black depth of the entire sub-tree rooted at the parent is 
          // now too small - fix it up recursively.
          currP = currP->parentP;
        }
      }
      else
      {
        if (siblingP == NULL)
        {
          // Take special care of the case of a NULL sibling.
          if (currP->parentP->color == Node::Red)
          {
            currP->parentP->color = Node::Black;
            currP = rootP;          // In order to stop the while loop.
          }
          else
          {
            currP = currP->parentP;
          }
        }
        else
        {
          // In this case, at least one of the sibling's children is red. 
          // It is therfore obvious that the sibling itself is black.
          // RWRW: TO_CHECK: CGAL_assertion (siblingP->color == Node::Black); 
          
          if (siblingP->rightP != NULL &&
              siblingP->rightP->color == Node::Red)
          {
            // If the right child of the sibling is red, color it black and
            // rotate around the current parent.
            siblingP->rightP->color = Node::Black;
            _rotate_left (currP->parentP);
          }
          else
          {
            // If the left child of the sibling is red, rotate around the 
            // sibling, then rotate around the new sibling of our current
            // node.
            _rotate_right (siblingP);
            siblingP = currP->parentP->rightP;
            CGAL_assertion (siblingP != NULL);
            _rotate_left (siblingP);
          }
          
          // It is now safe to color the parent black and to terminate the 
          // fix-up process.
          if (currP->parentP->parentP != NULL)
            currP->parentP->parentP->color = currP->parentP->color;
          currP->parentP->color = Node::Black;
          currP = rootP;          // In order to stop the while loop.
        }
      }
    }
    else
    {
      // If the current node is a right child, its sibling is the left 
      // child of the parent.
      siblingP = currP->parentP->leftP;

      // Check the sibling's color. Notice that NULL nodes are treated
      // as if they are colored black.
      if (siblingP != NULL && siblingP->color == Node::Red)
      {
        // In case the sibling is red, color it black and rotate.
        // Then color the parent red (and the grandparent is now black).
        siblingP->color = Node::Black;
        currP->parentP->color = Node::Red;
        _rotate_right (currP->parentP);

        siblingP = currP->parentP->leftP;
      }

      if (siblingP != NULL &&
          (siblingP->leftP == NULL || 
           siblingP->leftP->color == Node::Black) && 
          (siblingP->rightP == NULL || 
           siblingP->rightP->color == Node::Black))
      {
        // If the sibling has two black children, color it red.
        siblingP->color = Node::Red;
        if (currP->parentP->color == Node::Red)
        {
          // If the parent is red, we can safely color it black and terminate
          // the fix-up process.
          currP->parentP->color = Node::Black;
          currP = rootP;          // In order to stop the while loop.
        }
        else
        {
          // The black depth of the entire sub-tree rooted at the parent is 
          // now too small - fix it up recursively.
          currP = currP->parentP;
        }
      }
      else
      {
        if (siblingP == NULL)
        {
          // Take special care of the case of a NULL sibling.
          if (currP->parentP->color == Node::Red)
          {
            currP->parentP->color = Node::Black;
            currP = rootP;          // In order to stop the while loop.
          }
          else
          {
            currP = currP->parentP;
          }
        }
        else
        {
          // In this case, at least one of the sibling's children is red. 
          // It is therfore obvious that the sibling itself is black.
          // RWRW: TO_CHECK: CGAL_assertion (siblingP->color == Node::Black);
          
          if (siblingP->leftP != NULL &&
              siblingP->leftP->color == Node::Red)
          {
            // If the left child of the sibling is red, color it black and
            // rotate around the current parent.
            siblingP->leftP->color = Node::Black;
            _rotate_right (currP->parentP);
          }
          else
          {
            // If the right child of the sibling is red, rotate around the 
            // sibling, then rotate around the new sibling of our current 
            // node.
            _rotate_left (siblingP);
            siblingP = currP->parentP->leftP;
            CGAL_assertion (siblingP != NULL);
            _rotate_right (siblingP);
          }

          // It is now safe to color the parent black and to terminate the 
          // fix-up process.
          if (currP->parentP->parentP != NULL)
            currP->parentP->parentP->color = currP->parentP->color;
          currP->parentP->color = Node::Black;
          currP = rootP;          // In order to stop the while loop.
        }
      }
    }
  }

  // The root can always be colored black.
  currP->color = Node::Black;

  return;
}


//------------------------------------------------
//
//
template <class TYPE, class COMP, typename Alloc>
 inline typename Red_black_tree<TYPE, COMP, Alloc>::Node* 
        Red_black_tree<TYPE, COMP, Alloc>::
        _allocate_node(const TYPE& object,typename Node::Node_color color)
{
  Node* new_node = m_node_allocator.allocate(1);
  m_node_allocator.construct(new_node, m_master_node);
  new_node->init(object, color);
  return new_node;
}

CGAL_END_NAMESPACE

#endif
