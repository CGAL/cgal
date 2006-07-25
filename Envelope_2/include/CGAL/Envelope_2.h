// ======================================================================
//
// file          : Envelope_2.h
// author(s)     : Ron Wein <wein@post.tau.ac.il> 
//
// ======================================================================
#ifndef CGAL_LU_ENVELOPE_2_H
#define CGAL_LU_ENVELOPE_2_H

#include "Allocator.h"
#include <list>

CGAL_BEGIN_NAMESPACE

template <class _Traits>
class Envelope_2
{
public:
  typedef _Traits                              Traits;
  typedef typename Traits::Point_2             Point_2;
  typedef typename Traits::X_monotone_curve_2  X_monotone_curve_2;
  typedef typename Traits::Curve_2             Curve_2;

  enum Envelope_type
  {
    LOWER,
    UPPER
  };

protected:

  /*!
   * A structure used to store an x-monotone curve along with a pointer to its
   * original curve.
   */
  struct M_curve_2
  {
    X_monotone_curve_2     xcv;          // The x-monotine curve.
    const Curve_2          *cvP;         // Pointer to its original curve.

    M_curve_2 (const X_monotone_curve_2& _xcv,
               const Curve_2 *_cvP) :
      xcv(_xcv),
      cvP(_cvP)
    {}
  };

protected:
  typedef std::list<X_monotone_curve_2>             X_curves_list;
  typedef typename X_curves_list::const_iterator    X_curves_iter;
  typedef std::list<M_curve_2>                      M_curves_list;
  typedef typename M_curves_list::const_iterator    M_curves_iter;
  typedef std::list<const M_curve_2 *>              M_curves_P_list;
  typedef typename M_curves_P_list::const_iterator  M_curves_P_iter;

  /*!
   * Representation of a vertex in the minimization (or maximization) diagram:
   * Holds the point it represents (the vertex of the lower/upper envelope) and
   * a list of all curves that pass at that point.
   */
  struct M_diagram_edge_1;

  struct M_diagram_vertex_1
  {
    Point_2          p;
    M_diagram_edge_1 *leftP;
    M_diagram_edge_1 *rightP;

    /*!
     * Constructor.
     */
    M_diagram_vertex_1 () :
      p(),
      leftP(NULL),
      rightP(NULL)
    {}
  };

  /*!
   * Representation of an edge in the minimization (or maximization) diagram:
   * Essentially it hold an interval and points to all curves that form the
   * lower/upper envelope above it (usually there's just one curve, unless the
   * input contains overlapping curves).
   */
  struct M_diagram_edge_1
  {
    const M_curve_2    *mcvP;
    M_diagram_vertex_1 *leftP;
    M_diagram_vertex_1 *rightP;

    /*!
     * Default constructor.
     */
    M_diagram_edge_1 () :
      mcvP(NULL),
      leftP(NULL),
      rightP(NULL)
    {}
  };

  /*!
   * Representation of a minimization (or maximization) diagram. 
   */
  struct M_diagram_1
  {
    Envelope_2<Traits> *envP;        // The lower envelope structure.
    M_diagram_vertex_1    *leftmostP;   // The leftmost vertex of the diagram.
    M_diagram_vertex_1    *rightmostP;  // The rightmost vertex of the diagram.

    /*!
     * Default constructor.
     */
    M_diagram_1 () :
      envP(NULL),
      leftmostP(NULL),
      rightmostP(NULL)
    {}

    /*!
     * Constructor.
     * \param _envP The envelope that owns the diagram.
     */
    M_diagram_1 (Envelope_2<Traits> *_envP) :
      envP(_envP),
      leftmostP(NULL),
      rightmostP(NULL)
    {}

    /*!
     * Destructor.
     */
    ~M_diagram_1 ()
    {
      reset();
    }

    /*!
     * Clear the diagram.
     */
    void reset ()
    {
      M_diagram_edge_1  *edgeP = NULL;

      while (leftmostP != NULL)
      {
        // Get a pointer to the next edge.
        edgeP = leftmostP->rightP;

        // Free the leftmost vertex.
        envP->vert_alloc.Free (leftmostP);

        // Update the leftmost vertex.
        if (edgeP != NULL)
        {
          leftmostP = edgeP->rightP;

          // Free the current edge.
          envP->edge_alloc.Free (edgeP);
        }
        else
        {
          leftmostP = NULL;
        }
      }
      leftmostP = NULL;
      rightmostP = NULL;

      return;
    }

    /*!
     * Append a vertex to the diagram: The new vertex that represents the 
     * given point as the new rightmost vertex of the diagram. The edge 
     * between the current rightmost vertex and the new one contains the same 
     * curves as the input edge.
     * \param p The point that the new vertex is associated with.
     * \param e The input edge, or NULL if the connecting edge should be empty.
     * \pre The diagram is not empty.
     */
    void append_vertex (const Point_2& p,
                        const M_diagram_edge_1* e)
    {
      CGAL_assertion(leftmostP != NULL && rightmostP != NULL);

      CGAL_assertion(envP->get_traits()->compare_xy_2_object() (rightmostP->p,
                                                     p) == SMALLER);

      // Create the new vertex and the new edge.
      M_diagram_vertex_1 *new_v = envP->vert_alloc.Allocate();
      M_diagram_edge_1   *new_e = envP->edge_alloc.Allocate();

      new_v->p = p;
      
      if (e != NULL)
        new_e->mcvP = e->mcvP;
      else
        new_e->mcvP = NULL;

      // Connect them to the right of the current rightmost vertex of outd.
      new_v->leftP = new_e;
      new_v->rightP = NULL;
      new_e->rightP = new_v;
      new_e->leftP = rightmostP;
      rightmostP->rightP = new_e;
      rightmostP = new_v;

      return;
    }    
  };
  friend struct M_diagram_1;

  // Data members:
  Traits         *traits;        // The traits object.
  bool           own_traits;     // Whether we own the traits object.
  Envelope_type  env_type;       // Either LOWER or UPPER.
  M_curves_list  m_curves_list;  // Stores all the x-monotone curves.
  M_diagram_1    diagram;        // The minimization (or maximization) diagram.

  Allocator<M_diagram_vertex_1>   vert_alloc;
  Allocator<M_diagram_edge_1>     edge_alloc;

  // Copy constructor and assignment operator - not supported.
  Envelope_2 (const Envelope_2& );
  const Envelope_2& operator= (const Envelope_2& );

public:

  /*!
   * Constructor.
   */
  Envelope_2 () :
    own_traits(true),
    env_type(LOWER),
    m_curves_list(),
    diagram()
  {
    diagram.envP = this;
    traits = new Traits;
  }

  /*!
   * Constructor with a traits object.
   * \param _traits The traits object.
   */
  Envelope_2 (const Traits& _traits) :
    own_traits(false),
    env_type(LOWER),
    m_curves_list(),
    diagram()
  {
    diagram.envP = this;
    traits = &_traits;
  }  

  /*!
   * Destructor.
   */
  virtual ~Envelope_2 ()
  {
    _reset();

    if (own_traits)
      delete traits;
  }

  /*!
   * Construct the lower (or upper) envelope to the given range of curves.
   * \param begin An iterator pointing at the beginning of the curves range. 
   * \param end A past-the-end iterator for the curves range.
   * \param type The envelope type (LOWER or UPPER).
   */
  template <class CurvesIterator>
  void insert (const CurvesIterator& begin,
               const CurvesIterator& end,
               const Envelope_type& type)
  {
    // Reset the current structures.
    _reset ();

    // Set the envelope type.
    env_type = type;

    // Prepare the list of x-monotone curves.
    CurvesIterator   it;
    for (it = begin; it != end; it++)
    {
      // Split the current curve to x-monotone sub-curves.

      std::list<CGAL::Object>    objects;
      traits->make_x_monotone_2_object()(*it, std::back_inserter(objects));

      std::list<CGAL::Object>::iterator obj_itr = objects.begin();

      for(; obj_itr != objects.end(); ++obj_itr)
      {
        X_monotone_curve_2 cv;
        if(CGAL::assign(cv, *obj_itr))
        {
          m_curves_list.push_back (M_curve_2(cv, &(*it)));
        }
      }
    }

    // Construct the lower/upper envelope.
    _construct_envelope ();
  }

  /*!
   * An iterator for the output minization (or maximization) diagram vertices.
   */
  class Vertex_iterator
  {
    friend class Envelope_2<Traits>;

  private:

    const M_diagram_1        *diagramP;
    const M_diagram_vertex_1 *vertexP;

    Vertex_iterator (const M_diagram_1& _diagram,
                     const M_diagram_vertex_1 *_vertexP) :
      diagramP(&_diagram),
      vertexP(_vertexP)
    {}

  public:

    /*!
     * Default constrcutor.
     */
    Vertex_iterator () :
      diagramP(NULL),
      vertexP(NULL)
    {}

    /*! 
     * Compare with another iterator.
     */
    bool operator== (const Vertex_iterator& vi) const
    {
      return (diagramP == vi.diagramP && vertexP == vi.vertexP);
    }

    bool operator!= (const Vertex_iterator& vi) const
    {
      return (diagramP != vi.diagramP || vertexP != vi.vertexP);
    }

    /*!
     * Increment the iterator.
     */
    void operator++ ()
    {
      CGAL_precondition(diagramP != NULL);

      if (vertexP == NULL)
        vertexP = diagramP->leftmostP;
      else if (vertexP->rightP != NULL)
        vertexP = vertexP->rightP->rightP;
      else
        vertexP = NULL;
    }

    void operator++ (int )
    {
      CGAL_precondition(diagramP != NULL);

      if (vertexP == NULL)
        vertexP = diagramP->leftmostP;
      else if (vertexP->rightP != NULL)
        vertexP = vertexP->rightP->rightP;
      else
        vertexP = NULL;
    }

    /*!
     * Decrement the iterator.
     */
    void operator-- ()
    {
      CGAL_precondition(diagramP != NULL);

      if (vertexP == NULL)
        vertexP = diagramP->rightmostP;
      else if (vertexP != NULL && vertexP->leftP != NULL)
        vertexP = vertexP->leftP->leftP;
      else
        vertexP = NULL;
    }

    void operator-- (int )
    {
      CGAL_precondition(diagramP != NULL);

      if (vertexP == NULL)
        vertexP = diagramP->rightmostP;
      else if (vertexP != NULL && vertexP->leftP != NULL)
        vertexP = vertexP->leftP->leftP;
      else
        vertexP = NULL;
    }

    /*!
     * Return the point associated with current vertex.
     */
    const Point_2& operator* () const
    {
      CGAL_precondition(vertexP != NULL);
      return (vertexP->p);
    }

    /*!
     * Return a curve to the right of the current vertex.
     */
    const Curve_2* right () const
    {
      CGAL_precondition(vertexP != NULL);

      if (vertexP->rightP == NULL || vertexP->rightP->mcvP == NULL)
        return (NULL);

      return (vertexP->rightP->mcvP->cvP);
    }

    /*!
     * Return a curve to the left of the current vertex.
     */
    const Curve_2* left () const
    {
      CGAL_precondition(vertexP != NULL);

      if (vertexP->leftP == NULL || vertexP->leftP->mcvP == NULL)
        return (NULL);

      return (vertexP->leftP->mcvP->cvP);
    }
  };

  friend class iterator;

  /*!
   * Get an iterator pointing at the leftmost vertex of the envelope.
   * \return An iterator to start the scan with.
   */
  Vertex_iterator begin () const
  {
    return (Vertex_iterator(diagram, diagram.leftmostP));
  }

  /*!
   * Get an iterator pointing at the rightmost vertex of the envelope.
   * \return An iterator to start the reverse scan with.
   */
  Vertex_iterator rbegin () const
  {
    return (Vertex_iterator(diagram, diagram.rightmostP));
  }

  /*!
   * Get an iterator indicating the last envelope vertex.
   * \return A past-the-end iterator.
   */
  Vertex_iterator end () const
  {
    return (Vertex_iterator(diagram, NULL));
  }

  /*!
   * Get an iterator indicating the last envelope vertex.
   * \return A past-the-end iterator for a reversed scan.
   */
  Vertex_iterator rend () const
  {
    return (Vertex_iterator(diagram, NULL));
  }

  /*!
   * Reset the envelope object.
   */
  void reset ()
  {
      _reset();
      return;
  }

  /*!
   * Get the traits object.
   * \return A pointer to the traits object.
   */
  Traits* get_traits () const
  {
    return (traits);
  }

protected:

  /*!
   * Reset the current structures.
   */
  void _reset ()
  {
    // Free the list of x-monotone curves.
    m_curves_list.clear();

    // Free the minimization diagram.
    diagram.reset();

    return;
  }

  /*!
   * Construct the lower/upper envelope of the x-montone curves currently
   * in the list.
   */
  void _construct_envelope ()
  {
    // Move all the vertical segments to a different list.
    M_curves_P_list reg_list;
    M_curves_P_list vert_list;
    M_curves_iter   m_iter;

    for (m_iter = m_curves_list.begin();
         m_iter != m_curves_list.end(); m_iter++)
    {
      if (traits->is_vertical_2_object()((*m_iter).xcv))
        vert_list.push_back(&(*m_iter));
      else
        reg_list.push_back(&(*m_iter));
    }

    // Construct the envelope for the non-vertical curves.
    _construct_envelope_non_vertical (reg_list,
                                      diagram);

    // Merge the vertical segments.
    if (vert_list.size() > 0)
      _merge_vertical_segments (vert_list,
                                diagram);

    return;
  }

  /*!
   * Construct the lower/upper envelope of the given list of non-vertical
   * curves.
   * \param curves_list The list of curves.
   * \param outd The output minimization (or maximization) diagram.
   */
  void _construct_envelope_non_vertical (const M_curves_P_list& curves_list,
                                         M_diagram_1& outd)
  {
    M_curves_P_iter    iter;

    if (curves_list.size() == 1)
    {
      // In case the list contains a single curve:
      iter = curves_list.begin();
      const M_curve_2 *mcvP = *iter;

      // Clear the output diagram.
      outd.reset();

      // Find the left and right end-points of the curve.
      Point_2           p_left;
      Point_2           p_right;

      p_left = traits->construct_min_vertex_2_object()(mcvP->xcv);
      p_right = traits->construct_max_vertex_2_object()(mcvP->xcv);

      // Create the two end vertices in the diagram.
      outd.leftmostP = vert_alloc.Allocate();
      outd.leftmostP->p = p_left;

      outd.rightmostP = vert_alloc.Allocate();
      outd.rightmostP->p = p_right;

      // Create the single diagram edge and link it to the vertices.
      M_diagram_edge_1 *edgeP = edge_alloc.Allocate();
      edgeP->mcvP = mcvP;
      edgeP->leftP = outd.leftmostP;
      edgeP->rightP = outd.rightmostP;

      outd.leftmostP->leftP = NULL;
      outd.leftmostP->rightP = edgeP;
      outd.rightmostP->leftP = edgeP;
      outd.rightmostP->rightP = NULL;
    }
    else
    {
      // In case of a list that contains at least 2 curves, break it to two 
      // lists, each containing one half of the curves.
      M_curves_P_list  list1, list2;
      M_diagram_1      d1(this), d2(this);

      iter = curves_list.begin();
      while (iter != curves_list.end())
      {
        list1.push_back(*iter);
        iter++;

        if (iter != curves_list.end())
        {
          list2.push_back(*iter);
          iter++;
        }
      }

      // Construct the diagrams (envelopes) for the two lists recursively 
      // and then merge the two diagrams to obtain the result.
      _construct_envelope_non_vertical (list1, d1);
      _construct_envelope_non_vertical (list2, d2);

      _merge_envelopes (d1, d2, outd);

      // Free the memory consumed by the two intermediate diagrams d1 and d2.
      d1.reset();
      d2.reset();
    }

    return;
  }

  /*
   * Merge two minimization (or maximization) diagrams.
   * \param d1 The first diagram, 
   *           representing the envelope of the curve set C1.
   * \param d1 The first diagram,
   *           representing the envelope of the curve set C1.
   * \param outd The merged diagram,
   *             representing the envelope of the union of C1 and C2.
   */
  void _merge_envelopes (const M_diagram_1& d1,
                         const M_diagram_1& d2,
                         M_diagram_1& outd)
  {
    std::cout<<"_merge_envelopes:\n";
    std::cout<<"d1.leftmostP: " << d1.leftmostP->p<<"\n";
    std::cout<<"d2.leftmostP: " << d2.leftmostP->p<<"\n";
    // First clear the output diagram.
    outd.reset();

    // Find which of the diagrams starts to the left, and just copy its
    // relevant portion to outd.
    const M_diagram_vertex_1 *v1 = d1.leftmostP;  // The current vertex of d1.
    const M_diagram_vertex_1 *v2 = d2.leftmostP;  // The current vertex of d2.
    const M_diagram_vertex_1 *u = NULL;           // The previous vertex.
    Comparison_result        res;
    bool                     same_x;

    // Construct the first vertex of the output diagram.
    res = _compare_vertices (v1, v2, same_x);

    outd.leftmostP = vert_alloc.Allocate();
    outd.leftmostP->leftP = NULL;
    outd.leftmostP->rightP = NULL;

    if (res == SMALLER)
    {
      outd.leftmostP->p = v1->p;
      u = v1;

      v1 = (v1->rightP != NULL) ? (v1->rightP->rightP) : NULL;
      if (same_x)
        v2 = (v2->rightP != NULL) ? (v2->rightP->rightP) : NULL;
    }
    else if (res == LARGER)
    {
      outd.leftmostP->p = v2->p;
      u = v2;

      v2 = (v2->rightP != NULL) ? (v2->rightP->rightP) : NULL;
      if (same_x)
        v1 = (v1->rightP != NULL) ? (v1->rightP->rightP) : NULL;
    }
    else // (res == EQUAL)
    {
      outd.leftmostP->p = v1->p;
      u = v1;

      v1 = (v1->rightP != NULL) ? (v1->rightP->rightP) : NULL;
      v2 = (v2->rightP != NULL) ? (v2->rightP->rightP) : NULL;
    }

    outd.rightmostP = outd.leftmostP;

    // Proceed and process the rest of the diagrams.
    while (v1 != NULL && v2 != NULL)
    {
      // Get pointers to two curves representing the lower (or upper) envelopes
      // of d1 and d2 in the current interval.
      const X_monotone_curve_2 *cv1P = _curve_to_left(v1);
      const X_monotone_curve_2 *cv2P = _curve_to_left(v2);

      // Find the next vertex to the right of (outd.rightmostP).
      res = _compare_vertices (v1, v2, same_x);

      if (cv1P != NULL && cv2P != NULL)
      {
        // Both curves are not NULL.
        _merge_two_intervals (u,
                              v1->leftP,
                              *cv1P, 
                              v2->leftP,
                              *cv2P,
                              (res == SMALLER) ? v1 : v2,
                              (res == SMALLER) ? 1 : 2,
                              outd);
      }
      else if (cv1P != NULL && cv2P == NULL)
      {
        // cv1P forms the current interval:
        _merge_single_interval (v1->leftP,
                                (res == SMALLER) ? v1 : v2,
                                (res != LARGER),
                                outd);
      }
      else if (cv1P == NULL && cv2P != NULL)
      {
        // cv1P forms the current interval:
        _merge_single_interval (v2->leftP,
                                (res == SMALLER) ? v1 : v2,
                                (res != SMALLER),
                                outd);
      }
      else
      {
        // An empty interval in both d1 and d2:
        if (res == SMALLER)
          outd.append_vertex (v1->p, NULL);
        else
          outd.append_vertex (v2->p, NULL);
      }

      // Proceed to the next vertex.
      if (res == SMALLER)
      {
        u = v1;
        v1 = (v1->rightP != NULL) ? (v1->rightP->rightP) : NULL;
        if (same_x)
          v2 = (v2->rightP != NULL) ? (v2->rightP->rightP) : NULL;
      }
      else if (res == LARGER)
      {
        u = v2;
        v2 = (v2->rightP != NULL) ? (v2->rightP->rightP) : NULL;
        if (same_x)
          v1 = (v1->rightP != NULL) ? (v1->rightP->rightP) : NULL;
      }
      else
      {
        u = v1;
        v1 = (v1->rightP != NULL) ? (v1->rightP->rightP) : NULL;
        v2 = (v2->rightP != NULL) ? (v2->rightP->rightP) : NULL;
      }

      // Make sure that u is not to the left of the rightmost vertex in the
      // output diagram.
      if (traits->compare_x_2_object() (u->p, outd.rightmostP->p) == SMALLER)
        u = outd.rightmostP;
    }

    // Append any diagram tails to the merged diagram.
    if (v1 != NULL)
    {
      while (v1 != NULL)
      {
        outd.append_vertex (v1->p, v1->leftP);
        v1 = (v1->rightP != NULL) ? (v1->rightP->rightP) : NULL;
      }
    }
    else if (v2 != NULL)
    {
      while (v2 != NULL)
      {
        outd.append_vertex (v2->p, v2->leftP);
        v2 = (v2->rightP != NULL) ? (v2->rightP->rightP) : NULL;
      }
    }

    return;
  }

  /*!
   * Compare two vertices.
   * \param v1 The first vertex.
   * \param v2 The second vertex.
   * \param same_x Output parameter: TRUE iff x(v1->p) = x(v2->p).
   * \return SMALLER if x(v1->p) < x(v2->p). Or, in case x(v1->p) = x(v2->p):
   *                 if we compute the lower envelope and y(v1->p) < y(v2->p),
   *                 if we compute the upper envelope and y(v1->p) > y(v2->p).
   *         LARGER if x(v1->p) > x(v2->p). Or, in case x(v1->p) = x(v2->p):
   *                if we compute the lower envelope and y(v1->p) > y(v2->p),
   *                if we compute the upper envelope and y(v1->p) < y(v2->p).
   *         EQUAL if v1->p = v2->0.
   */
  Comparison_result _compare_vertices (const M_diagram_vertex_1 *v1,
                                       const M_diagram_vertex_1 *v2,
                                       bool& same_x) const
  {
    Comparison_result   res = traits->compare_x_2_object() (v1->p, v2->p);

    if (res != EQUAL)
    {
      same_x = false;
      return (res);
    }
    else
    {
      same_x = true;
    }

    // In case x(v1->p) = x(v2->p):
    res = traits->compare_xy_2_object() (v1->p, v2->p);

    if ((env_type == LOWER && res == SMALLER) ||
        (env_type == UPPER && res == LARGER))
      return (SMALLER);
    else if ((env_type == LOWER && res == LARGER) ||
             (env_type == UPPER && res == SMALLER))
      return (LARGER);

    return (EQUAL);
  }

  /*!
   * Get a pointer to an x-monotone curve that lies to the left of the given 
   * vertex, or return NULL if there is no such curve.
   * \param v The vertex.
   * \return A pointer to a curve that lies to the left of v.
   */
  const X_monotone_curve_2 *_curve_to_left (const M_diagram_vertex_1 *v) const
  {
    if (v->leftP == NULL)
      // No edge to the left of v:
      return (NULL);
    else if (v->leftP->mcvP == NULL)
      // The edge to the left is empty (contains to curves)
      return (NULL);

    // Return a pointer to the first curve in the list.
    return (&(v->leftP->mcvP->xcv));
  }

  /*!
   * Deal with an interval which is non-empty in one of the merged diagram and
   * empty in the other.
   * \param e The non-empty edge.
   * \param v The next vertex (which is right of outd.rightmostP).
   * \param same_org Whether e and v originate from the same diagram.
   * \param outd The merged diagram.
   */
  void _merge_single_interval (const M_diagram_edge_1* e,
                               const M_diagram_vertex_1* v,
                               const bool& same_org,
                               M_diagram_1& outd) const
  {
    if (same_org)
    {
      // The non-empty edge ends at v, so we insert it to outd.
      outd.append_vertex (v->p, e);
      return;
    }

    // If v is not on e, we should insert it to the mered diagram only if it
    // is below (or above, in case of an upper envelope) the curves of e.
    const X_monotone_curve_2& cv = e->mcvP->xcv;
    Comparison_result         res = traits->compare_y_at_x_2_object() (v->p, cv);

    if ((env_type == LOWER && res == SMALLER) ||
        (env_type == UPPER && res == LARGER))
    {
      outd.append_vertex (v->p, e);
    }

    return;
  }

  /*!
   * Merge two non-empty intervals into the merged diagram.
   * \param u The previous vertex of the merged diagrams we have visited.
   * \param e1 The first non-empty edge.
   * \param cv1 A curve that belongs to the first edge in this interval.
   * \param e2 The second non-empty edge.
   * \param cv2 A curve that belongs to the second edge in this interval.
   * \param v The next vertex (which is right of outd.rightmostP).
   * \param org_v The origin of v: 1 if it is from e1, 2 if it is from e2.
   * \param outd The merged diagram.
   */
  void _merge_two_intervals (const M_diagram_vertex_1* u,
                             const M_diagram_edge_1* e1,
                             const X_monotone_curve_2& cv1,
                             const M_diagram_edge_1* e2,
                             const X_monotone_curve_2& cv2,
                             const M_diagram_vertex_1* v,
                             const int& org_v,
                             M_diagram_1& outd) const
  {
    // Compare the two curves to the right of (outd.rightmostP).
    Comparison_result res = traits->curves_compare_y_at_x (cv1, cv2,
                                                           u->p);

    
   /* if (org_v == 1)
    {
      res = traits->compare_y_at_x_2_object() (v->p, cv2);
    }
    else
    {
      res = traits->compare_y_at_x_2_object() (v->p, cv1);
      if(res == SMALLER)
        res = LARGER;
      else
        if(res == LARGER)
          res = SMALLER;
    }*/

    
    if (res == EQUAL)
      res = traits->compare_y_at_x_left_2_object() (cv1, cv2,
                                                    v->p);
    if (env_type == UPPER)
    {
      // Flip the result in case of an upper envelope.
      if (res == SMALLER)
        res = LARGER;
      else if (res == LARGER)
        res = SMALLER;
    }

    // Find the next intersection of the envelopes to the right of the current
    // rightmost point in the merged diagram.
    CGAL::Object       obj;
    X_monotone_curve_2 icv;
    Point_2            p1, p2;

    while (true)
    {
      std::list<CGAL::Object>    objects;
      traits->intersect_2_object()(cv1, cv2, std::back_inserter(objects));
      CGAL::Object obj;
      
      // Stop if no intersection has been found.
      if (objects.empty())
        break;

      for(std::list<CGAL::Object>::iterator itr = objects.begin();
          itr != objects.end();
          ++itr)
      {
        std::pair<Point_2, unsigned int>  x_pt; 
        if(CGAL::assign(x_pt, *itr))
        {
          if(traits->compare_xy_2_object()(x_pt.first, outd.rightmostP->p) == LARGER)
          {
            obj = make_object(x_pt.first);
            break;
          }
        }
        else
        {
          X_monotone_curve_2 cv;
          CGAL::assign(cv, *itr);
          CGAL_assertion(CGAL::assign(cv, *itr));
          Point_2 pt = traits->construct_min_vertex_2_object()(cv);
          if(traits->compare_xy_2_object()(pt, outd.rightmostP->p) == LARGER)
          {
            obj = *itr;
            break;
          }
        }
      }

      if(obj.is_empty())
        break;
     /* obj = traits->nearest_intersection_to_right (cv1, cv2,
                                                   outd.rightmostP->p);*/

      // Check for overlaps (when the returned object is a curve).
      if (CGAL::assign (icv, obj))
      {
        // Assign the leftmost endpoint in the overlapping curve to p1,
        // and the rightmost endpoint to p2.
        p1 = traits->construct_min_vertex_2_object() (icv);
        p2 = traits->construct_max_vertex_2_object() (icv);

        if (traits->compare_x_2_object() (p1, p2) != SMALLER)
        {
          Point_2 tmp = p1;
          p1 = p2;
          p2 = tmp;
        }

        // Disregard intersection points that are not to the left of v.
        if (traits->compare_x_2_object() (p1, v->p) != SMALLER)
          break;

        // Deal with overlaps:
        // If the overlap does not start at the rightmost vertex of the merged
        // diagram, insert p1 as a vertex.
        if (traits->compare_x_2_object() (p1, outd.rightmostP->p) == LARGER)
        {
          // RWRW: bug fix!
          if (res == EQUAL)
            res = SMALLER;

          if (res == SMALLER)
            outd.append_vertex (p1, e1);
          else if (res == LARGER)
            outd.append_vertex (p1, e2);
          else
            // This case should never occur:
            CGAL_assertion (res != EQUAL);
        }

        // Make the leftmost point of p2 and v->p a vertex. The edge to its
        // left contains all curves in e1 and in e2.
        bool  reached_v = (traits->compare_x_2_object() (p2, v->p) != SMALLER);

        if (reached_v)
          outd.append_vertex (v->p, e1);
        else
          outd.append_vertex (p2, e1);

        M_diagram_edge_1 *left_e = outd.rightmostP->leftP;
        left_e->mcvP = e2->mcvP;

        if (reached_v)
          return;

        // Compare the two curves to the right of the overlap.
        res = traits->compare_y_at_x_right_2_object() (cv1, cv2,
                                                   p2);

        if (env_type == UPPER)
        {
          // Flip the result in case of an upper envelope.
          if (res == SMALLER)
            res = LARGER;
          else if (res == LARGER)
            res = SMALLER;
        }
      }
      else if (CGAL::assign (p1, obj))
      {
        // Disregard intersection points that are not to the left of v.
        if (traits->compare_x_2_object() (p1, v->p) != SMALLER)
          break;

        // RWRW: bug fix!
        if (res == EQUAL)
            res = SMALLER;

        // Make p1 the new rightmost vertex of the merged diagram.
        if (res == SMALLER)
          outd.append_vertex (p1, e1);
        else if (res == LARGER)
          outd.append_vertex (p1, e2);
        else
          // This case should never occur (should be treated as an overlap):
          CGAL_assertion(res != EQUAL);

        // Compare the two curves to the right of the intersection point.
        res = traits->compare_y_at_x_right_2_object() (cv1, cv2,
                                                   p1);

        // This case should never occur (should be treated as an overlap):
        // RWRW: bug fix!
        if (res == EQUAL) res = SMALLER;
        // CGAL_assertion(res != EQUAL);

        if (env_type == UPPER)
        {
          // Flip the result in case of an upper envelope.
          if (res == SMALLER)
            res = LARGER;
          else if (res == LARGER)
            res = SMALLER;
        }
      }
      else
      {
        // This case should never occur:
        CGAL_assertion (false);
      }
    }

    // Check if v should also be inserted to the merged diagram.
    if (res == SMALLER)
    {
      // The final part of the interval is taken from e1.
      if (org_v == 1)
      {
        // In case v is also from e1, append it to the merged diagram.
        outd.append_vertex (v->p, e1);
      }
      else
      {
        // If v is from e2, check if it below (or above, in case of an upper
        // envelope) cv1 to insert it.
        res = traits->compare_y_at_x_2_object() (v->p, cv1);

        if (res == EQUAL ||
            (env_type == LOWER && res == SMALLER) ||
            (env_type == UPPER && res == LARGER))
        {
          outd.append_vertex (v->p, e1);
        }
      }
    }
    else if (res == LARGER)
    {
      // The final part of the interval is taken from e2.
      if (org_v == 2)
      {
        // In case v is also from e2, append it to the merged diagram.
        outd.append_vertex (v->p, e2);
      }
      else
      {
        // If v is from e1, check if it below (or above, in case of an upper
        // envelope) cv2 to insert it.
        res = traits->compare_y_at_x_2_object() (v->p, cv2);

        if (res == EQUAL ||
            (env_type == LOWER && res == SMALLER) ||
            (env_type == UPPER && res == LARGER))
        {
          outd.append_vertex (v->p, e2);
        }
      }
    }

    return;
  }

  /*!
  * A functor used to sort vertical segments by their x-coordinate.
  */
  class Vertical_strict_weak_ordering
  {
  private:    
    const Traits         *traits;

  public:
    Vertical_strict_weak_ordering (const Traits *_traits) :
      traits(_traits)
    {}

    bool operator() (const M_curve_2 *mcv1, const M_curve_2 *mcv2) const
    {
      return (traits->compare_x_2_object() (traits->construct_min_vertex_2_object()(mcv1->xcv),
              traits->construct_min_vertex_2_object()(mcv1->xcv)) == SMALLER);
    }
  };

  /*!
   * Merge the vertical segments into the lower/upper envelope given as a
   * minimization (or maximization) diagram.
   * \param vert_list The list of vertical segments.
   * \param outd The input minimization (or maximization) diagram.
   *             The function merges the vertical segments into this diagram.
   */
  void _merge_vertical_segments (M_curves_P_list& vert_list,
                                 M_diagram_1& outd)
  {
    // Sort the vertical segments by their increasing x-coordinate.
    Vertical_strict_weak_ordering vert_order(traits);

    // RWRW: Perform bubble-sort, as the following does not compile under Windows:
    // vert_list.sort (vert_order);
    typename M_curves_P_list::iterator  i_iter, j_iter;

    for (i_iter = vert_list.begin(); i_iter != vert_list.end(); i_iter++)
    {
      j_iter = i_iter;
      j_iter++;
      for (; j_iter != vert_list.end(); j_iter++)
      {
        if (! vert_order (*i_iter, *j_iter))
        {
          // Swap the pointers.
          const M_curve_2 *temp = *i_iter;
          *i_iter = *j_iter;
          *j_iter = temp;
        }
      }
    }

    // Go over all vertical segments that are to the left of the leftmost
    // vertex of the diagram.
    M_curves_P_iter    iter = vert_list.begin();
    M_diagram_vertex_1 *u = NULL;
    M_diagram_vertex_1 *v = outd.leftmostP;
    Comparison_result  res;
    Point_2            q;

    while (v != NULL)
    {
      while (iter != vert_list.end() &&
             traits->compare_x_2_object() (traits->construct_min_vertex_2_object()((*iter)->xcv),
                                v->p) == SMALLER)
      {
        // Get the lower (or the upper) point of the vertical segment.
        res = traits->compare_xy_2_object() (traits->construct_min_vertex_2_object()((*iter)->xcv),
          traits->construct_max_vertex_2_object()((*iter)->xcv));

        if ((env_type == LOWER && res == SMALLER) ||
          (env_type == UPPER && res == LARGER))
          q = traits->construct_min_vertex_2_object()((*iter)->xcv);
        else
          q = traits->construct_max_vertex_2_object()((*iter)->xcv);

        // Act according to the previous vertex u.
        if (u == NULL)
        {
          // The vertical segment is to the left of the leftmost diagram
          // vertex.
          M_diagram_vertex_1 *new_v1 = vert_alloc.Allocate();
          M_diagram_vertex_1 *new_v2 = vert_alloc.Allocate();
          M_diagram_edge_1   *new_e = edge_alloc.Allocate();
          M_diagram_edge_1   *empty_e = edge_alloc.Allocate();

          new_v1->p = q;
          new_v2->p = q;

          new_v1->leftP = NULL;
          new_v1->rightP = new_e;

          new_e->mcvP = *iter;
          new_e->leftP = new_v1;
          new_e->rightP = new_v2;

          new_v2->leftP = new_e;
          new_v2->rightP = empty_e;

          empty_e->mcvP = NULL;
          empty_e->leftP = new_v2;
          empty_e->rightP = outd.leftmostP;

          outd.leftmostP->leftP = empty_e;
          outd.leftmostP = new_v1;

          // Update the pointer to diagram vertex immediately to the left of v.
          u = new_v2;
        }
        else if (traits->compare_x_2_object() (q, u->p) == EQUAL)
        {
          // The vertical segment has the same x-coordinate as u.
          res = traits->compare_xy_2_object() (q, u->p);

          // Insert a new curve only if it is below (or above, in case of an
          // upper envelope) u->p.
          if ((env_type == LOWER && res == SMALLER) ||
            (env_type == UPPER && res == LARGER))
          {
            M_diagram_edge_1   *new_e = edge_alloc.Allocate();
            M_diagram_vertex_1 *new_v = vert_alloc.Allocate();

            new_v->p = q;

            new_e->mcvP = *iter;
            new_e->leftP = u;
            new_e->rightP = new_v;

            new_v->leftP = new_e;
            new_v->rightP = u->rightP;

            u->rightP->leftP = new_v;
            u->p = q;
            u->rightP = new_e;

            // Update the pointer to diagram vertex immediately to the 
            // left of v.
            u = new_v;
          }
        }
        else
        {
          // The vertical segment is placed in between u and v.
          bool                     add_q = false;
          const X_monotone_curve_2 *cvP = _curve_to_left(v);

          if (cvP == NULL)
          {
            // The edge between u and v is empty:
            add_q = true;
          }
          else
          {
            // Check whether q lies below (or above, in case of an upper
            // envelope) the curves of the edge to the left of u.
            res = traits->compare_y_at_x_2_object() (q, *cvP);

            add_q = (res == EQUAL) ||
                    (env_type == LOWER && res == SMALLER) ||
                    (env_type == UPPER && res == LARGER);
          }

          if (add_q)
          {
            // Cut the edge to the left of v and insert the vertical segment.
            M_diagram_vertex_1 *new_v1 = vert_alloc.Allocate();
            M_diagram_vertex_1 *new_v2 = vert_alloc.Allocate();
            M_diagram_edge_1   *new_e = edge_alloc.Allocate();
            M_diagram_edge_1   *dup_e = edge_alloc.Allocate();

            new_v1->p = q;
            new_v2->p = q;
            *dup_e = *(v->leftP);

            new_v1->leftP = v->leftP;
            new_v1->rightP = new_e;

            new_e->mcvP = *iter;
            new_e->leftP = new_v1;
            new_e->rightP = new_v2;

            new_v2->leftP = new_e;
            new_v2->rightP = dup_e;

            dup_e->leftP = new_v2;
            dup_e->rightP = v;

            v->leftP->rightP = new_v1;
            v->leftP = dup_e;

            // Update the pointer to diagram vertex immediately to the 
            // left of v.
            u = new_v2;
          }
        }

        // Move to the next vertical segment.
        iter++;
      }

      // Move to the next diagram vertex.
      u = v;
      if (v->rightP != NULL)
        v = v->rightP->rightP;
      else
        v = NULL;
    }

    // Deal with all segments located to the right of the diagram.
    while (iter != vert_list.end())
    {
      // Get the lower (or the upper) point of the vertical segment.
      res = traits->compare_xy_2_object() (traits->construct_min_vertex_2_object()((*iter)->xcv),
        traits->construct_max_vertex_2_object()((*iter)->xcv));

      if ((env_type == LOWER && res == SMALLER) ||
        (env_type == UPPER && res == LARGER))
        q = traits->construct_min_vertex_2_object()((*iter)->xcv);
      else
        q = traits->construct_max_vertex_2_object()((*iter)->xcv);

      // The vertical segment is to the right of the rightmost diagram vertex.
      M_diagram_vertex_1 *new_v1 = vert_alloc.Allocate();
      M_diagram_vertex_1 *new_v2 = vert_alloc.Allocate();
      M_diagram_edge_1   *new_e = edge_alloc.Allocate();
      M_diagram_edge_1   *empty_e = edge_alloc.Allocate();

      new_v1->p = q;
      new_v2->p = q;

      empty_e->mcvP = NULL;
      empty_e->leftP = outd.rightmostP;
      empty_e->rightP = new_v1;

      new_v1->leftP = empty_e;
      new_v1->rightP = new_e;

      new_e->mcvP = *iter;
      new_e->leftP = new_v1;
      new_e->rightP = new_v2;

      new_v2->leftP = new_e;
      new_v2->rightP = NULL;

      outd.rightmostP->rightP = empty_e;
      outd.rightmostP = new_v2;

      // Move to the next vertical segment.
      iter++;
    }

    return;
  }

};

CGAL_END_NAMESPACE

#endif

