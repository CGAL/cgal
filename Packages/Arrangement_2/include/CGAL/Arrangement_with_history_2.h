// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
#ifndef CGAL_ARRANGEMENT_WITH_HISTORY_2_H
#define CGAL_ARRANGEMENT_WITH_HISTORY_2_H

/*! \file
 * The header file is for the Arrangement_with_history_2<Traits,Dcel> class.
 */

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/In_place_list.h>

#include <set>

CGAL_BEGIN_NAMESPACE

template <class T_Traits, 
          class T_Dcel = Arr_default_dcel<Traits_> > 
class Arrangement_with_history_2 :
  public Arrangement_2<Arr_consolidated_curve_data_traits_2
    <T_Traits, typename T_Traits::Curve_2*>,
                       T_Dcel>
{
public:
  typedef T_Traits                                      Traits_2;
  typedef T_Dcel                                        Dcel;

  typedef typename Traits_2::Point_2                    Point_2;
  typedef typename Traits_2::Curve_2                    Curve_2;
  typedef typename Traits_2::X_monotone_curve_2         X_monotone_curve_2;

private:

  typedef Arr_consolidated_curve_data_traits_2<Traits_2,Curve_2*>
                                                        Data_traits;
  typedef Data_traits::Data_iterator                    Data_iterator;
  typedef Arrangement_2<Data_traits,Dcel>               Base;
  
  /*! Compare functor */
  struct Less_halfedge_handle
  {
    bool operator() (Halfedge_handle h1, Halfedge_handle h2) const
    {
      return (&(*h1) < &(*h2));
    }
  };

  /*! Retains the halfedges that correspond to a curve it inherits from
   */
  class Curve_halfedges : public Curve_2,
                          public In_place_list_base<Curve_halfedges>
  {
  private:
    typedef std::set<Halfedge_handle, Less_halfedge_handle>     Halfedges_set;

    Halfedges_set m_halfedges;
    
  public:
    typedef typename Halfedges_set::const_iterator              const_iterator;
    typedef typename Halfedges_set::iterator                    iterator;

    /*! Constructor */
    Curve_halfedges(const Curve_2 & curve) : Curve_2(curve) {}

    /*! Insert
     */
    iterator insert(Halfedge_handle he)
    {
      std::pair<iterator, bool> res = m_halfedges.insert(he);
      CGAL_assertion(res.second);
      return res.first;
    }
    
    /*! Erase
     */
    void erase(iterator pos)
    {
      m_halfedges.erase(pos);
    }
    
    /*! Erase
     */
    void erase(Halfedge_handle he)
    {
      size_t res = m_halfedges.erase(he);
      if (res == 0) res = m_halfedges.erase(he->twin());
      CGAL_assertion(res != 0);
    }

    /*! Begin */
    iterator begin()
    {
      return m_halfedges.begin();
    }

    /*! Begin */
    const_iterator begin() const
    {
      return m_halfedges.begin();
    }
    
    /*! End */
    iterator end()
    {
      return m_halfedges.end();
    }

    /*! End */
    const_iterator end() const
    {
      return m_halfedges.end();
    }
  };

  typedef CGAL_ALLOCATOR(Curve_halfedges)               Curves_alloc;
  typedef In_place_list<Curve_halfedges, false>         Curve_halfedges_list;

  /*! Observer */
  class Curve_halfedge_observer : public Arr_observer<Base>
  {
  private:
    /*! Register the halfedge in the curve halfedges
     */
     void _register_edge (Halfedge_handle e)
    {
      Curve_halfedges * curve_halfedges;
      Data_iterator di;
      for (di = e->curve().data_begin(); di != e->curve().data_end(); ++di)
      {
        curve_halfedges = static_cast<Curve_halfedges*>(*di);
        curve_halfedges->insert(e);
      }
    }

    /*! Unregister the halfedge in the curve halfedges
     */
     void _unregister_edge (Halfedge_handle e)
    {
      Curve_halfedges * curve_halfedges;
      Data_iterator di;
      for (di = e->curve().data_begin(); di != e->curve().data_end(); ++di)
      {
        curve_halfedges = static_cast<Curve_halfedges*>(*di);
        curve_halfedges->erase(e);
      }
    }

  public:
    /*!
     * Notification after the creation of a new edge.
     * \param e A handle to one of the twin halfedges that were created.
     */
    virtual void after_create_edge (Halfedge_handle e)
    {
      _register_edge(e);
    }

    /*!
     * Notification before the modification of an existing edge.
     * \param e A handle to one of the twin halfedges to be updated.
     * \param c The x-monotone curve to be associated with the edge.
     */
    virtual void before_modify_edge (Halfedge_handle e,
                                     const X_monotone_curve_2& /* c */)
    {
      _unregister_edge(e);
    }

    /*!
     * Notification after an edge was modified.
     * \param e A handle to one of the twin halfedges that were updated.
     */
    virtual void after_modify_edge (Halfedge_handle e)
    {
      _register_edge(e);
    }
    
    /*!
     * Notification before the splitting of an edge into two.
     * \param e A handle to one of the existing halfedges.
     * \param c1 The x-monotone curve to be associated with the first edge.
     * \param c2 The x-monotone curve to be associated with the second edge.
     */
    virtual void before_split_edge (Halfedge_handle e,
                                    const X_monotone_curve_2& /* c1 */,
                                    const X_monotone_curve_2& /* c2 */)
    {
      _unregister_edge(e);
    }

    /*!
     * Notification after an edge was split.
     * \param e1 A handle to one of the twin halfedges forming the first edge.
     * \param e2 A handle to one of the twin halfedges forming the second edge.
     */
    virtual void after_split_edge (Halfedge_handle e1, Halfedge_handle e2)
    {
      _register_edge(e1);
      _register_edge(e2);
    }

    /*!
     * Notification before the merging of two edges.
     * \param e1 A handle to one of the halfedges forming the first edge.
     * \param e2 A handle to one of the halfedges forming the second edge.
     * \param c The x-monotone curve to be associated with the merged edge.
     */
    virtual void before_merge_edge (Halfedge_handle e1, Halfedge_handle e2,
                                    const X_monotone_curve_2& /* c */)
    {
      _unregister_edge(e1);
      _unregister_edge(e2);
    }

    /*!
     * Notification after an edge was merged.
     * \param e A handle to one of the twin halfedges forming the merged edge.
     */
    virtual void after_merge_edge (Halfedge_handle e)
    {
      _register_edge(e);
    }

    /*!
     * Notification before the removal of an edge.
     * \param e A handle to one of the twin halfedges to be deleted.
     */
    virtual void before_remove_edge (Halfedge_handle e)
    {
      _unregister_edge(e);
    }
  };

  // Data members:
  Curves_alloc         m_curves_alloc;
  Curve_halfedges_list m_curves;
  Curve_halfedge_observer m_observer;
  
public:
  typedef Curve_halfedges_list::iterator                Curve_iterator;
  typedef Curve_iterator                                Curve_handle;
  typedef Curve_halfedges_list::const_iterator          Curve_const_iterator;
  typedef Curve_const_iterator                          Curve_const_handle;

  /*! Default constructor */
  Arrangement_with_history_2 ()
  {
    m_observer.attach(*this);
  }

  /*! Copy constructor. \todo not implemented yet */
  Arrangement_with_history_2 (const Self & arr);

  /*! Constructor given a traits object. */
  Arrangement_with_history_2 (Traits_2 * tr) :
    Base(static_cast<Data_traits*>(tr))
  {
    m_observer.attach(*this);
  }

  /*! Insert */
  Curve_handle insert(const Curve_2 & curve)
  {
    Curve_hafledges   *p_cv = m_curves_alloc.allocate (1);
    Curve_handle       ch;
    
    m_curves_alloc.construct (p_cv, curve);
    ch = m_curves.push_back (*p_cv);

    typename Data_traits::Curve_2 data_curve(curve, p_cv);
    insert(*this, data_curve);
    
    return (ch);
  }
  
  /*! Remove */
  size_t remove(Curve_handle & ch)
  {
  }
};

CGAL_END_NAMESPACE

#endif
