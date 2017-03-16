// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://balasmic@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_on_surface_2/include/CGAL/Arr_point_location/Trapezoidal_decomposition_2_impl.h $
// $Id: Arr_trapezoid_ric_pl_impl.h 56667 2010-06-09 07:37:13Z sloriot $
//
//
// Author(s)     : Michal Balas   <balasmic@post.tau.ac.il>
//                based on the methods implemented in Trapezoidal_decomposition.h by Oren Nechushtan

#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_FUNCTIONS_H
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_FUNCTIONS_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
* Member-function definitions for the Trapezoidal_decomposition_2<Traits>
* class.
*/

#include <algorithm>

#include <CGAL/assertions.h>
#include <CGAL/use.h>

namespace CGAL {

//-----------------------------------------------------------------------------
// Description:
//  splits the trapezoid with vertical line through v
//  assuming that he_bottom_ray_shoot & he_top_ray_shoot are in the
//  desired direction, such that v is their source
// Precondition:
//  The trapezoid is active and contains v in its closure
//
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Dag_node &
Trapezoidal_decomposition_2<Td_traits>
::split_trapezoid_by_vertex(Dag_node& split_node,
                            Vertex_const_handle v,
                            Halfedge_const_handle he)
{
  CGAL_precondition(!split_node.is_null());
  if (split_node.is_null())  return split_node;


  Td_map_item curr_item (split_node.get_data());
  CGAL_precondition(traits->is_active(curr_item));

  Dag_node left_node, right_node;

  if (traits->is_td_trapezoid(curr_item))
  {
    Td_active_trapezoid& tr (boost::get<Td_active_trapezoid>(curr_item));

    CGAL_warning(traits->is_in_closure(tr, traits->vtx_to_ce(v)));

    left_node.set_data(Td_active_trapezoid
                        (tr.left(), v, tr.bottom(), tr.top()));

    right_node.set_data(Td_active_trapezoid
                        (v, tr.right(), tr.bottom(), tr.top()));

    Td_active_trapezoid& left_tr  (boost::get<Td_active_trapezoid>(left_node.get_data()));
    Td_active_trapezoid& right_tr (boost::get<Td_active_trapezoid>(right_node.get_data()));

    CGAL_warning(traits->is_trapezoids_top_equal(left_tr,right_tr));
    CGAL_warning(traits->is_trapezoids_bottom_equal(left_tr,right_tr));
    CGAL_warning(left_tr.is_on_left_boundary() == tr.is_on_left_boundary());
    CGAL_warning(right_tr.is_on_right_boundary() == tr.is_on_right_boundary());

    left_tr.init_neighbours(tr.lb(),tr.lt(),right_node.get_data(),right_node.get_data());
    right_tr.init_neighbours(left_node.get_data(),left_node.get_data(),tr.rb(),tr.rt());
    if (!traits->is_empty_item(tr.lb()))
    {
      Td_active_trapezoid& lb(boost::get<Td_active_trapezoid>(tr.lb()));
      lb.set_rb(left_node.get_data());
    }
    if (!traits->is_empty_item(tr.lt()))
    {
      Td_active_trapezoid& lt(boost::get<Td_active_trapezoid>(tr.lt()));
      lt.set_rt(left_node.get_data());
    }
    if (!traits->is_empty_item(tr.rb()))
    {
      Td_active_trapezoid& rb(boost::get<Td_active_trapezoid>(tr.rb()));
      rb.set_lb(right_node.get_data());
    }
    if (!traits->is_empty_item(tr.rt()))
    {
      Td_active_trapezoid& rt(boost::get<Td_active_trapezoid>(tr.rt()));
      rt.set_lt(right_node.get_data());
    }

  }
  else //if the curr_item is an edge
  {
    Td_active_edge& e (boost::get<Td_active_edge>(curr_item));

    CGAL_warning(traits->is_in_closure(e, traits->vtx_to_ce(v)));

    left_node.set_data(Td_active_edge(e.halfedge()));

    right_node.set_data(Td_active_edge(e.halfedge()));

    Td_active_edge& left_e  (boost::get<Td_active_edge>(left_node.get_data()));
    Td_active_edge& right_e (boost::get<Td_active_edge>(right_node.get_data()));

    //CGAL_warning(left_e.is_on_left_boundary() == e.is_on_left_boundary());
    //CGAL_warning(right_e.is_on_right_boundary() == e.is_on_right_boundary());

    left_e.init_neighbours(boost::none);//left_e.init_neighbours(e.lb(),e.lt(),Td_map_item(),right_node.get_data());
    right_e.init_neighbours(e.next());//right_e.init_neighbours(left_node.get_data(),left_node.get_data(),e.rb(),e.rt());

  }

  // left and right are set to the point itself,
  // bottom and top are set to the ray shooting resulting curves at this
  // stage.
  //need to define the boundaries flag before creating the trapezoid
  Td_map_item vtx_item (build_vertex_map_item(v, he, &split_node));
  split_node.replace( vtx_item, left_node, right_node); //nodes depth are updated here
  update_largest_leaf_depth( (std::max)(left_node.depth(), right_node.depth()) );
  m_number_of_dag_nodes += 2;

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
  if (last_cv == curr_item)
    last_cv = vtx_item;
  if (prev_cv == curr_item)
    prev_cv = vtx_item;
#endif

  const Dag_node* left_ptr  = &split_node.left_child();
  const Dag_node* right_ptr = &split_node.right_child();

  boost::apply_visitor(set_dag_node_visitor((Dag_node*)left_ptr),left_ptr->get_data()); //(*left_ptr)->set_dag_node((Dag_node*)left_ptr);
  boost::apply_visitor(set_dag_node_visitor((Dag_node*)right_ptr),right_ptr->get_data()); //(*right_ptr)->set_dag_node((Dag_node*)right_ptr);

  //print_dag_addresses(*m_dag_root);

  return split_node;
}

template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Td_map_item
Trapezoidal_decomposition_2<Td_traits>
::build_vertex_map_item(Vertex_const_handle v,
                        Halfedge_const_handle he,
                        Dag_node* node)
{
  Curve_end ce(traits->vtx_to_ce(v));
  if ((traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce()) == ARR_INTERIOR)
      && (traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce()) == ARR_INTERIOR))
  {
    Td_active_vertex vtx (v, he, node);
    return vtx;
  }
  else
  {
    Td_active_fictitious_vertex vtx (v, he, node);
    return vtx;
  }
}


//-----------------------------------------------------------------------------
// Description:
//  the opposite operation for spliting the trapezoid with
//  vertical line through ce
// Precondition:
//  tr_node data is a td vertex and is active
template <class Td_traits>
void Trapezoidal_decomposition_2<Td_traits>
::undo_split_trapezoid_by_vertex(Dag_node& tr_node, const Curve_end& ce)
{
  CGAL_precondition(!tr_node.is_null());
  Td_map_item item = tr_node.get_data();
  CGAL_precondition(traits->is_active(item));
  CGAL_precondition(traits->is_td_vertex(item));

  //get the ds left child and right child nodes of
  //    tr_node (in the search structure)
  Dag_node tr_left_node  = tr_node.left_child();
  Dag_node tr_right_node = tr_node.right_child();

  //advances the ds nodes until ce is found
  search_using_dag(tr_left_node,  traits, ce, m_empty_he_handle);
  search_using_dag(tr_right_node, traits, ce, m_empty_he_handle);

  //make sure the trapezoids are active before merging them
  //CGAL_assertion(traits->is_active(tr_left_node.get_data()) && traits->is_active(tr_right_node.get_data()));

  // bool mrg_res =
  merge_if_possible(tr_left_node.get_data(), tr_right_node.get_data());

  CGAL_warning(!tr_left_node.is_inner_node());
  CGAL_warning(!tr_right_node.is_inner_node());


  deactivate_trapezoid( tr_right_node, &tr_left_node); //tr_right_node->remove(&tr_left_node);
  update_largest_leaf_depth( tr_left_node.depth()); //tr_left_node is not an inner node
  // mark root as deleted
  deactivate_vertex(tr_node); //tr_node->remove();
  //no need to update m_number_of_dag_nodes because the number of nodes did not change.
  // removed nodes were only marked as removed
}


template <class Td_traits>
void Trapezoidal_decomposition_2<Td_traits>
::deactivate_trapezoid (Dag_node& trpz_node, Dag_node* active_node) const
{
  CGAL_precondition(traits->is_active(trpz_node.get_data()));
  CGAL_precondition(traits->is_td_trapezoid(trpz_node.get_data()));
  if ( Td_active_trapezoid* trap = boost::get<Td_active_trapezoid>(&trpz_node.get_data()) )
    trap->clear_neighbors();
  trpz_node.set_data(Td_inactive_trapezoid());
  if (active_node)
    trpz_node.set_left_child(*active_node);
}

template <class Td_traits>
void Trapezoidal_decomposition_2<Td_traits>
::deactivate_vertex (Dag_node& vtx_node) const
{
  CGAL_precondition(traits->is_active(vtx_node.get_data()));
  CGAL_precondition(traits->is_td_vertex(vtx_node.get_data()));
  if (traits->is_fictitious_vertex(vtx_node.get_data()))
  {
    Td_active_fictitious_vertex& v(boost::get<Td_active_fictitious_vertex>(vtx_node.get_data()));
    vtx_node.set_data(Td_inactive_fictitious_vertex(v.vertex(), &vtx_node));
  }
  else
  {
    Td_active_vertex& v(boost::get<Td_active_vertex>(vtx_node.get_data()));
    vtx_node.set_data(Td_inactive_vertex(v.vertex(), &vtx_node));
  }
}

template <class Td_traits>
void Trapezoidal_decomposition_2<Td_traits>
::deactivate_edge (boost::shared_ptr<X_monotone_curve_2>& cv, Dag_node& edge_node) const
{
  CGAL_precondition(traits->is_active(edge_node.get_data()));
  CGAL_precondition(traits->is_td_edge(edge_node.get_data()));

  edge_node.set_data(Td_inactive_edge(cv, &edge_node));
}



//-----------------------------------------------------------------------------
// Description:
//  splits the trapezoid that corresponds to the root of the
//  trapezoidal tree with an input halfedge he
// Precondition:
//  The root trapezoid is active
//  The root trapezoid is devided by he or is equal to it and is vertical.
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Dag_node &
Trapezoidal_decomposition_2<Td_traits>
::split_trapezoid_by_halfedge(Dag_node& split_node,
                              Td_map_item& prev_e,
                              Td_map_item& prev_bottom_tr,
                              Td_map_item& prev_top_tr,
                              Halfedge_const_handle he)
{

  CGAL_warning(traits != NULL);
  CGAL_precondition(traits->is_active(split_node.get_data()));
  CGAL_precondition(traits->is_td_trapezoid(split_node.get_data()));

  Td_map_item curr_item (split_node.get_data());
  Td_active_trapezoid& split_tr = boost::get<Td_active_trapezoid>(curr_item);

  // sets left and right according to td_edge's source and target positions
  // sets bottom and top to Halfedge itself
  // no need to set the boundaries since its irrelevant for this trapezoid
  //   type (TD_EDGE)
  Td_active_edge sep(he, &split_node);


  //creates a one-way path for all the td edges
  //that represent the Halfedge.
  //rb() is used to retrieve the
  //next on path information

  Dag_node top_node(Td_active_trapezoid
                        (split_tr.left(), split_tr.right(), he, split_tr.top()));
                         //Td_active_trapezoid::TD_TRAPEZOID,
                         //split_tr.on_boundaries_flag() &
                         // (CGAL_TD_ON_LEFT_BOUNDARY  |
                         //  CGAL_TD_ON_RIGHT_BOUNDARY |
                         //  CGAL_TD_ON_TOP_BOUNDARY  )));
  Dag_node bottom_node(Td_active_trapezoid
                        (split_tr.left(),split_tr.right(), split_tr.bottom(), he));
                         //Td_active_trapezoid::TD_TRAPEZOID,
                         //split_tr.on_boundaries_flag() &
                         // (CGAL_TD_ON_LEFT_BOUNDARY  |
                         //  CGAL_TD_ON_RIGHT_BOUNDARY |
                         //  CGAL_TD_ON_BOTTOM_BOUNDARY )));

  Td_active_trapezoid& bottom = boost::get<Td_active_trapezoid>(bottom_node.get_data());
  Td_active_trapezoid& top    = boost::get<Td_active_trapezoid>(top_node.get_data());

  top.init_neighbours(prev_top_tr, split_tr.lt(), boost::none , split_tr.rt());
  bottom.init_neighbours(split_tr.lb(), prev_bottom_tr, split_tr.rb(), boost::none);

  if (!traits->is_empty_item(prev_bottom_tr))
  {
    Td_active_trapezoid& prev_btm (boost::get<Td_active_trapezoid>(prev_bottom_tr));
    prev_btm.set_rt(bottom_node.get_data());
  }
  if (!traits->is_empty_item(prev_top_tr))
  {
    Td_active_trapezoid& prev_top (boost::get<Td_active_trapezoid>(prev_top_tr));
    prev_top.set_rb(top_node.get_data());
  }
  if (!traits->is_empty_item(split_tr.lb()))
  {
    Td_active_trapezoid& lb (boost::get<Td_active_trapezoid>(split_tr.lb()));
    lb.set_rb(bottom_node.get_data());
  }
  if (!traits->is_empty_item(split_tr.lt()))
  {
    Td_active_trapezoid& lt (boost::get<Td_active_trapezoid>(split_tr.lt()));
    lt.set_rt(top_node.get_data());
  }
  if (!traits->is_empty_item(split_tr.rb()))
  {
    Td_active_trapezoid& rb (boost::get<Td_active_trapezoid>(split_tr.rb()));
    rb.set_lb(bottom_node.get_data());
  }
  if (!traits->is_empty_item(split_tr.rt()))
  {
    Td_active_trapezoid& rt (boost::get<Td_active_trapezoid>(split_tr.rt()));
    rt.set_lt(top_node.get_data());
  }
  split_node.replace(sep,bottom_node,top_node); //nodes depth are updated here
  update_largest_leaf_depth( (std::max)(bottom_node.depth(), top_node.depth()) );
  m_number_of_dag_nodes += 2; //two new nodes were added to the DAG

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
  if (last_cv == curr_item)
    last_cv = sep;
  if (prev_cv == curr_item)
    prev_cv = sep;
#endif

  const Dag_node* bottomPtr = &split_node.left_child();
  const Dag_node* topPtr    = &split_node.right_child();

  boost::apply_visitor(set_dag_node_visitor((Dag_node*)bottomPtr),
                       bottomPtr->get_data());
  boost::apply_visitor(set_dag_node_visitor((Dag_node*)topPtr),
                       topPtr->get_data());

  // Td_active_edge& new_e = boost::get<Td_active_edge>(split_node.get_data());
  if (!traits->is_empty_item(prev_e))
  {
    Td_active_edge& e ( boost::get<Td_active_edge>(prev_e));
    e.set_next(split_node.get_data());
  }
  //update these trapezoids pointers.
  // will be used for the next trapezoid that should be split
  //  by this Halfedge
  prev_bottom_tr = bottomPtr->get_data(); //(*bottomPtr).operator->();
  prev_top_tr    = topPtr->get_data(); //(*topPtr).operator->();
  prev_e         = split_node.get_data(); //tt.operator->();

  return split_node;
}


//-----------------------------------------------------------------------------
// Description:
//  update
//   tr.bottom()
//   vertical_ray_shoot downward from tr
//   tr.top()
//    vertical_ray_shoot upward from tr
//  update all the curves incident to the vertex that there's a new curve
//  starting from this vertex
//  this point must be an interior point and not a point on the boundaries,
//  since a point on the boundaries is related to one curve only
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Td_map_item &
Trapezoidal_decomposition_2<Td_traits>
::update_vtx_with_new_edge(Halfedge_const_handle he,
                           const Curve_end& ce,
                           Td_map_item& vtx_item,
                           const Locate_type&
                              CGAL_precondition_code(lt))
{
  CGAL_assertion(traits != NULL);
  CGAL_precondition(lt == POINT);
  CGAL_precondition(traits->is_active(vtx_item));

  //ee is interior
  CGAL_assertion((traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce())
                                                          == ARR_INTERIOR) &&
                 (traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce())
                                                          == ARR_INTERIOR) );
  //the underlying point of ce
  const Point& p = (ce.ce() == ARR_MIN_END) ?
                    traits->construct_min_vertex_2_object()(ce.cv()) :
                    traits->construct_max_vertex_2_object()(ce.cv()) ;

  //set cw to hold the halfedge whose source is p,
  // which is clockwise "smallest" starting from top (12 o'clock)
  Halfedge_const_handle cw_he (boost::apply_visitor(cw_he_visitor(), vtx_item));
  if (traits->compare_cw_around_point_2_object()
      (he->curve(), is_edge_to_right(he,p),
       cw_he->curve(),
       is_edge_to_right(cw_he,p), p) == SMALLER)
  {
    boost::apply_visitor(set_cw_he_visitor(he),vtx_item);//v_tr->set_top(he);
  }

  return vtx_item;
}

//-----------------------------------------------------------------------------
// Description:
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Td_map_item &
Trapezoidal_decomposition_2<Td_traits>
::insert_curve_at_vtx_using_dag (Halfedge_const_handle he,
                                   Vertex_const_handle v,
                                   Td_map_item& item,
                                   const Locate_type&
                                   CGAL_precondition_code(lt))
{
  CGAL_precondition(lt==TRAPEZOID || lt==UNBOUNDED_TRAPEZOID);

  Dag_node* node = boost::apply_visitor(dag_node_visitor(), item);

  CGAL_assertion(node != NULL);
  CGAL_assertion (he != m_empty_he_handle);


  //we need to use the halfedge whose source is v. //MICHAL: do we have to?
  if (he->source() == v)
 // if ((ce_pair.second == ARR_MIN_END && he->direction() == ARR_LEFT_TO_RIGHT) ||
 //     (ce_pair.second == ARR_MAX_END && he->direction() == ARR_RIGHT_TO_LEFT)  )
  {
    return *split_trapezoid_by_vertex(*node, v, he);
  }
  else
  {
    return *split_trapezoid_by_vertex(*node, v, he->twin());
  }
}

////-----------------------------------------------------------------------------
//// Description:
//template <class Td_traits>
//void Trapezoidal_decomposition_2<Td_traits>
//::set_trp_params_after_halfedge_update(Halfedge_const_handle old_he,
//                                       Halfedge_const_handle new_he,
//                                       Td_map_item& vtx_item)
//{
//  CGAL_precondition(old_he != m_empty_he_handle);
//  CGAL_precondition(traits->is_td_vertex(vtx_item));
//  CGAL_precondition(traits->is_active(vtx_item));
//
//  Halfedge_const_handle cw_he (boost::apply_visitor(cw_he_visitor(), vtx_item));
//
//  //make sure the cw_he is added in same direction as before
//  //  such that vtx_item is the source (done inside the set methods)
//  if (cw_he == old_he || cw_he->twin() == old_he)
//  {
//    boost::apply_visitor(set_cw_he_visitor(new_he), vtx_item);
//  }
//}

//-----------------------------------------------------------------------------
// Description:
template <class Td_traits>
void Trapezoidal_decomposition_2<Td_traits>
::update_vtx_cw_he_after_merge (const X_monotone_curve_2& old_cv,
                                Halfedge_const_handle new_he,
                                Td_map_item& vtx_item)
{
  CGAL_precondition(traits->is_td_vertex(vtx_item));
  CGAL_precondition(traits->is_active(vtx_item));

  Halfedge_const_handle cw_he (boost::apply_visitor(cw_he_visitor(), vtx_item));

  //make sure the cw_he is added in same direction as before
  //  such that v_tr is the source (done inside the set methods)
  if (traits->equal_2_object()(cw_he->curve(), old_cv))
  {
    boost::apply_visitor(set_cw_he_visitor(new_he), vtx_item);
  }
}

////MICHAL: currently not in use since split is implemented as removed and insert two
////         in case the split is done differenty - this method would have to be rewritten since it is obsolete
//template <class Td_traits>
//void Trapezoidal_decomposition_2<Td_traits>
//::set_trp_params_after_split_halfedge_update(Halfedge_const_handle new_he,
//                                            Td_map_item& vtx_item,
//                                            Halfedge_const_handle he1,
//                                            Halfedge_const_handle he2)
//{
//  CGAL_precondition(traits->is_td_vertex(vtx_item));
//  CGAL_precondition(traits->is_active(vtx_item));
//
//  Halfedge_const_handle top_he (boost::apply_visitor(top_he_visitor(), vtx_item));
//  Halfedge_const_handle bottom_he (boost::apply_visitor(bottom_he_visitor(), vtx_item));
//
//  //make sure the top & bottom are added in same direction as before
//  //  such that sep is the source (done inside the set methods)
//  if ((top_he == he1) || (top_he == he1->twin()) ||
//      (top_he == he2) || (top_he == he2->twin()) )
//  {
//     boost::apply_visitor(set_top_he_visitor(new_he), vtx_item); //v_tr.set_top(new_he);
//  }
//  if ((bottom_he == he1) || (bottom_he == he1->twin()) ||
//      (bottom_he == he2) || (bottom_he == he2->twin()) )
//  {
//    boost::apply_visitor(set_bottom_he_visitor(new_he), vtx_item); //v_tr.set_bottom(new_he);
//  }
//}

//-----------------------------------------------------------------------------
// Description:
//  update geometric boundary(top and bottom) for trapezoids
//  traveled along an iterator till end reached
//  precondition:
//  end==0 or end is on the path of the iterator
// postcondition:
//  end is pointer to the last trapezoid encountered,if any
template <class Td_traits>
void Trapezoidal_decomposition_2<Td_traits>
::update_map_items_after_merge (In_face_iterator& it,
                                Halfedge_const_handle old_he,
                                Halfedge_const_handle new_he,
                                Vertex_const_handle /* min_v */,
                                Vertex_const_handle /* max_v */,
                                Td_map_item& end)
{
  Td_map_item last_item = Td_map_item(0);

  while (!!it && !(*it == end))
  {
    Td_map_item curr_item = *it;
    CGAL_assertion(traits->is_active(curr_item));

    if (traits->is_td_edge(curr_item))
    {
      Td_active_edge& e(boost::get<Td_active_edge>(curr_item));
      if (e.halfedge() == old_he || e.halfedge() == old_he->twin())
        e.set_halfedge(new_he);
    }
    else if (traits->is_td_trapezoid(curr_item))
    {
      Td_active_trapezoid& tr(boost::get<Td_active_trapezoid>(curr_item));
      if (tr.bottom() == old_he || tr.bottom() == old_he->twin())
        tr.set_bottom(new_he);
      if (tr.top() == old_he || tr.top() == old_he->twin())
        tr.set_top(new_he);
    }
    else //if is_td_vertex
    {
      Halfedge_const_handle cw_he (boost::apply_visitor(cw_he_visitor(), curr_item));
      if (cw_he == old_he || cw_he == old_he->twin())
        boost::apply_visitor(set_cw_he_visitor(new_he), curr_item);
    }

    last_item = *it;
    ++it;
  }
  end = last_item;
}

//-----------------------------------------------------------------------------
// Description:
//  advances input Data structure using data structure,input point p and
//  possibly Halfedge p_he till
//  p is found(if p_he hadn't been given)
//  p_he is found(if p_he was given)
//  or
//  leaf node reached
// postcondition:
//  output is the closest active trapezoid to ce/p_he
// remark:
//  use this function with care!
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Locate_type
Trapezoidal_decomposition_2<Td_traits>
::search_using_dag (Dag_node& curr_node,
                    const Traits* traits,
                    const Point& p,
                    Halfedge_const_handle he,
                    Comparison_result up /*=EQUAL*/) const
{

  while(true)
  {
    //curr_node is the current pointer to node in the data structure
    //curr_item is the curent Td_map_item held in curr_node
    Td_map_item curr_item(curr_node.get_data());

    if (traits->is_td_vertex(curr_item))
    { // if the curr_item represents a vertex
      //bool is_fict_vtx = traits->is_fictitious_vertex(curr_item);
      //if ((is_fict_vtx  && is_end_point_left_low(p, *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)))) ||
      //    (!is_fict_vtx &&  is_end_point_left_low(p, boost::apply_visitor(point_for_vertex_visitor(), curr_item))) )
      if (is_end_point_left_low(p, curr_node))
      {
        curr_node = curr_node.left_child();
        continue;
      }
      //else if ((is_fict_vtx  && is_end_point_right_top(p, *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)))) ||
      //         (!is_fict_vtx && is_end_point_right_top(p, boost::apply_visitor(point_for_vertex_visitor(), curr_item))) )
      else if (is_end_point_right_top(p, curr_node))
      {
        curr_node = curr_node.right_child();
        continue;
      }
      //else if ((is_fict_vtx  && traits->equal_curve_end_2_object()(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), p)) ||
      //         (!is_fict_vtx && traits->equal_2_object()(boost::apply_visitor(point_for_vertex_visitor(), curr_item), p)) )
      else if (are_equal_end_points(p, curr_node))
      {
        if (he == m_empty_he_handle) //if he is the empty handle
        {
          if ( up == EQUAL )
          {      // point found!
            if (traits->is_active(curr_item))
              return POINT;
            curr_node = curr_node.left_child();
          }
          else if ( up == LARGER ) {          // vertical ray shut up
            curr_node = curr_node.right_child();
          }
          else /*if ( up == SMALLER ) */ {
            curr_node = curr_node.left_child();               // vertical ray shut down
          }
          continue;
        }
        else //if he was given
        {
          bool is_equal_to_he_min = traits->equal_curve_end_2_object()
                               (Curve_end(he,ARR_MIN_END), p);
          bool is_equal_to_he_max = traits->equal_curve_end_2_object()
                               (Curve_end(he,ARR_MAX_END), p);

          CGAL_assertion( is_equal_to_he_min || is_equal_to_he_max );
          CGAL_USE(is_equal_to_he_max);

          curr_node = is_equal_to_he_min ? curr_node.right_child() : curr_node.left_child();
          continue;
        }
      }
      else
      {
        CGAL_assertion(is_end_point_left_low(p,curr_node) ||
                       is_end_point_right_top(p,curr_node) ||
                       are_equal_end_points(p,curr_node));

                    //(is_fict_vtx &&
                    //  (is_end_point_left_low(p,*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item))) ||
                    //   is_end_point_right_top(p,*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item))) ||
                    //   traits->equal_curve_end_2_object()(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)),p))) ||
                    // (!is_fict_vtx &&
                    //  (is_end_point_left_low(p,boost::apply_visitor(point_for_vertex_visitor(), curr_item)) ||
                    //   is_end_point_right_top(p,boost::apply_visitor(point_for_vertex_visitor(), curr_item)) ||
                    //   traits->equal_2_object()(boost::apply_visitor(point_for_vertex_visitor(), curr_item),p))));

        return Locate_type();
      }
    }
    if (traits->is_td_edge(curr_item))
    { // if curr_item represents an edge,
      //   so top() is a real Halfedge with a curve() if curr_item is active
      //   or curr_item holds the curve if it is not active
      const X_monotone_curve_2& he_cv = *(boost::apply_visitor(cv_for_edge_visitor(), curr_item));

      Comparison_result cres = traits->compare_y_at_x_2_object()(p, he_cv);
      if (cres == SMALLER)
      {
        curr_node = curr_node.left_child();
        continue;
      }
      else if (cres == LARGER)
      {
        curr_node = curr_node.right_child();
        continue;
      }
      else
      {
        // p is on the curve itself
         CGAL_warning(
          (cres == EQUAL) &&
          (traits->compare_curve_end_x_2_object()
                        (p, Curve_end(he_cv,ARR_MAX_END)) != LARGER) &&
          (traits->compare_curve_end_x_2_object()
                        (p, Curve_end(he_cv,ARR_MIN_END)) != SMALLER));

        if (he == m_empty_he_handle) //if he was not given
        {
          // For a vertical curve, we always visit it after visiting
          // one of its endpoints.
          if ((up == EQUAL) || traits->is_vertical(curr_item))
          {
            if (traits->is_active(curr_item))
              return CURVE;
            curr_node = curr_node.left_child();
          }
          else if (up == LARGER)
          {
            curr_node = curr_node.right_child();
          }
          else
          { // if (up==SMALLER)
            curr_node = curr_node.left_child();
          }
          continue;
        }
        else //if he was given
        {
          //p is a parameter space interior point
          bool is_min_equal = traits->equal_curve_end_2_object()
                                      (Curve_end(he,ARR_MIN_END),
                                       Curve_end(he_cv,ARR_MIN_END));

          bool is_max_equal = traits->equal_curve_end_2_object()
                                      (Curve_end(he,ARR_MAX_END),
                                       Curve_end(he_cv,ARR_MAX_END));

          CGAL_warning (is_min_equal || is_max_equal);
          CGAL_USE(is_max_equal);

          Comparison_result res =
             is_min_equal ?
              traits->compare_cw_around_point_2_object()
                       (he_cv, is_curve_to_right(he_cv,p),
                        he->curve(), is_edge_to_right(he,p), p) :
              traits->compare_cw_around_point_2_object()
                       (he->curve(), is_edge_to_right(he,p),
                        he_cv, is_curve_to_right(he_cv,p), p ,false);

          switch(res)
          {
           case LARGER:
            curr_node = curr_node.right_child();
            break;
           case SMALLER:
            curr_node = curr_node.left_child();
            break;
           case EQUAL:
            switch(up)
            {
             case LARGER:
              curr_node = curr_node.right_child();
              break;
             case SMALLER:
              curr_node = curr_node.left_child();
              break;
             case EQUAL:
              if (traits->is_active(curr_item))
                return CURVE;
              curr_node = curr_node.left_child();
              break;
            }
            break;
          }
        }
      }
    }
    else
    {
      // if is_degenerate() == 0, meaning: curr_item is a real trapezoid
      if (traits->is_active(curr_item))
      {
        Td_active_trapezoid& tr = boost::get<Td_active_trapezoid>(curr_item);
        return tr.is_on_boundaries() ? UNBOUNDED_TRAPEZOID : TRAPEZOID;
      }
      curr_node = curr_node.left_child();
      continue;
    }
  }
}



////-----------------------------------------------------------------------------
//// Description:
////  advances input Data structure using data structure,input point p and
////  possibly Halfedge p_he till
////  p is found(if p_he hadn't been given)
////  p_he is found(if p_he was given)
////  or
////  leaf node reached
//// postcondition:
////  output is the closest active trapezoid to ce/p_he
//// remark:
////  use this function with care!
//template <class Td_traits>
//void
//Trapezoidal_decomposition_2<Td_traits>
//::search_and_print_using_dag (std::ostream& out,
//                              Dag_node& curr_node,
//                              const Traits* traits,
//                              const Point& p,
//                              Halfedge_const_handle he,
//                              Comparison_result up /*=EQUAL*/) const
//{
//  out << "QUERY: " << std::endl;
//  out << "x: " << CGAL::to_double(p.x())
//      << ", y: " << CGAL::to_double(p.y()) << std::endl;
//
//
//
//  while(true)
//  {
//    //curr_node is the current pointer to node in the data structure
//    //curr_item is the curent Td_map_item held in curr_node
//    Td_map_item curr_item(curr_node.get_data());
//
//    if (traits->is_td_vertex(curr_item))
//    {
//      // if the map item represents a vertex
//
//      out << " VERTEX : " ;
//      if (traits->is_active(curr_item))
//        out << " (active) ";
//      else
//        out << " (inactive) ";
//
//      // if the curr_item represents a vertex
//      bool is_fict_vtx = traits->is_fictitious_vertex(curr_item);
//      if (is_fict_vtx)
//      {
//        Curve_end vtx_ce(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)));
//        print_ce_data(vtx_ce.cv(), vtx_ce.ce(), out);
//      }
//      else
//      {
//        print_point_data(boost::apply_visitor(point_for_vertex_visitor(),curr_item), out);
//      }
//
//      if ((is_fict_vtx && is_end_point_left_low(p, *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)))) ||
//        (!is_fict_vtx &&  is_end_point_left_low(p, boost::apply_visitor(point_for_vertex_visitor(),curr_item))) )
//      {
//        out << " Going left " << std::endl;
//        curr_node = curr_node.left_child();
//        continue;
//      }
//      else if ((is_fict_vtx && is_end_point_right_top(p, *(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)))) ||
//               (!is_fict_vtx &&  is_end_point_right_top(p, boost::apply_visitor(point_for_vertex_visitor(),curr_item))) )
//      {
//        out << " Going right " << std::endl;
//        curr_node = curr_node.right_child();
//        continue;
//      }
//      else if ((is_fict_vtx && traits->equal_curve_end_2_object()(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), p)) ||
//               (!is_fict_vtx &&  traits->equal_2_object()(boost::apply_visitor(point_for_vertex_visitor(),curr_item), p)) )
//      {
//        out << " Equal to query " << std::endl;
//        if (he == m_empty_he_handle) //if he is the empty handle
//        {
//          if ( up == EQUAL )
//          {      // point found!
//            if (traits->is_active(curr_item))
//            {
//              out << " Found active point! " << std::endl;
//              return;
//            }
//            out << " (equal to inactive point) Going left " << std::endl;
//            curr_node = curr_node.left_child();
//          }
//          else if ( up == LARGER )
//          {          // vertical ray shut up
//            out << " Going right " << std::endl;
//            curr_node = curr_node.right_child();
//          }
//          else //if ( up == SMALLER )
//          {       // vertical ray shut down
//            out << " Going left " << std::endl;
//            curr_node = curr_node.left_child();
//          }
//          continue;
//        }
//        else //if he was given
//        {
//          //NOT GONNA HAPPEN
//        }
//      }
//      else
//      {
//        out << " Problem - comparing to point" << std::endl;
//        return;
//      }
//    }
//    if (traits->is_td_edge(curr_item))
//    {
//      // if curr_item represents an edge,
//      //   so top() is a real Halfedge with a curve() if curr_item is active
//      //   or curr_item holds the curve if it is not active
//      const X_monotone_curve_2& he_cv = *(boost::apply_visitor(cv_for_edge_visitor(), curr_item));
//
//      out << " EDGE : " ;
//      if (traits->is_active(curr_item))
//      {
//        out << " (active) ";
//        if (traits->is_vertical(curr_item))
//          out << " (vertical) ";
//      }
//      else
//        out << " (inactive) ";
//
//      print_cv_data(he_cv, out);
//
//      Comparison_result cres = traits->compare_y_at_x_2_object()(p, he_cv);
//      if (cres == SMALLER)
//      {
//        out << " Going left " << std::endl;
//        curr_node = curr_node.left_child();
//        continue;
//      }
//      else if (cres == LARGER)
//      {
//        out << " Going right " << std::endl;
//        curr_node = curr_node.right_child();
//        continue;
//      }
//      else
//      {
//        // p is on the CURVE (top_he = curr.top()) itself
//        out << " query is on the curve " << std::endl;
//
//        CGAL_assertion(
//          (cres == EQUAL) &&
//          (traits->compare_curve_end_x_2_object()
//                        (p, Curve_end(he_cv,ARR_MAX_END)) != LARGER) &&
//          (traits->compare_curve_end_x_2_object()
//                        (p, Curve_end(he_cv,ARR_MIN_END)) != SMALLER));
//        if (he == m_empty_he_handle) //if he was not given
//        {
//          // For a vertical curve, we always visit it after visiting
//          // one of its endpoints.
//          if ((up == EQUAL) || traits->is_vertical(curr_item))
//          {
//            if (traits->is_active(curr_item))
//            {
//              out << " On active curve " << std::endl;
//              return;
//            }
//            out << " (equal to inactive curve) Going left " << std::endl;
//            curr_node = curr_node.left_child();
//          }
//          else if (up == LARGER) {
//            out << " Going right " << std::endl;
//            curr_node = curr_node.right_child();
//          }
//          else { // if (up==SMALLER)
//            out << " Going left " << std::endl;
//            curr_node = curr_node.left_child();
//          }
//          continue;
//        }
//        else //if he was given
//        {
//          //NOT GONNA HAPPEN
//        }
//      }
//    }
//    else
//    {
//      // if item represents a trapezoid
//      out << " TRAPEZOID : " ;
//      if (traits->is_active(curr_item))
//      {
//        out << " (active) ";
//        Td_active_trapezoid tr =  boost::get<Td_active_trapezoid>(curr_item);
//        if (tr.is_on_boundaries())
//          out << " UNBOUNDED! ";
//        else
//          out << " BOUNDED! ";
//        return;
//      }
//      else
//        out << " (inactive) ";
//      out << " (on inactive trapezoid) Going left " << std::endl;
//      curr_node = curr_node.left_child();
//      continue;
//    }
//  }
//}

//-----------------------------------------------------------------------------
// Description:
//  advances input Data structure using data structure,input point ce and
//  possibly Halfedge p_he till
//  ce is found(if p_he hadn't been given)
//  p_he is found(if p_he was given)
//  or
//  leaf node reached
// postcondition:
//  output is the closest active trapezoid to ce/p_he
// remark:
//  use this function with care!
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Locate_type
Trapezoidal_decomposition_2<Td_traits>
::search_using_dag (Dag_node& curr_node,
                    const Traits* traits,
                    const Curve_end& ce,
                    Halfedge_const_handle he,
                    Comparison_result up /*=EQUAL*/) const
{
  if (he == m_empty_he_handle)
    return search_using_dag_with_cv (curr_node,traits,ce,NULL, up);
  else
    return search_using_dag_with_cv (curr_node,traits,ce,&he->curve(), up);
}

//-----------------------------------------------------------------------------
// Description:
//  advances input Data structure using data structure,input point ce and
//  possibly X_monotone_curve_2 p_cv till
//  ce is found(if p_cv hadn't been given)
//  p_cv is found(if p_cv was given)
//  or
//  leaf node reached
// postcondition:
//  output is the closest active trapezoid to ce/p_cv
// remark:
//  use this function with care!
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Locate_type
Trapezoidal_decomposition_2<Td_traits>
::search_using_dag_with_cv (Dag_node& curr_node,
                            const Traits* traits,
                            const Curve_end& ce,
                            const X_monotone_curve_2* p_cv,
                            Comparison_result up /*=EQUAL*/) const
{
  while(true)
  {
    //curr_node is the current pointer to node in the data structure
    //curr_item is the curent Td_map_item held in curr_node
    Td_map_item curr_item(curr_node.get_data());

    if (traits->is_td_vertex(curr_item))
    { // if the curr_item represents a vertex
      //bool is_fict_vtx = traits->is_fictitious_vertex(curr_item);
      //if ((is_fict_vtx  && is_end_point_right_top(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), ce)) ||
      //    (!is_fict_vtx && is_end_point_right_top(boost::apply_visitor(point_for_vertex_visitor(), curr_item), ce)) )
      if (is_end_point_left_low(ce, curr_node))
      {
        curr_node = curr_node.left_child();
        continue;
      }
      //else if ((is_fict_vtx  && is_end_point_left_low(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), ce)) ||
      //         (!is_fict_vtx && is_end_point_left_low(boost::apply_visitor(point_for_vertex_visitor(), curr_item), ce)) )
      else if (is_end_point_right_top(ce, curr_node))
      {
        curr_node = curr_node.right_child();
        continue;
      }
      //else if ((is_fict_vtx && traits->equal_curve_end_2_object()(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), ce)) ||
      //         (!is_fict_vtx &&  traits->equal_curve_end_2_object()(boost::apply_visitor(point_for_vertex_visitor(), curr_item), ce)) )
      else if (are_equal_end_points(ce, curr_node))
      {
        if (!p_cv) //if p_cv was not given
        {
          if ( up == EQUAL )
          {      // point found!
            if (traits->is_active(curr_item))
              return POINT;
            curr_node = curr_node.left_child();
          }
          else if ( up == LARGER )
          {          // vertical ray shut up
            curr_node = curr_node.right_child();
          }
          else
          { // if ( up == SMALLER )
            curr_node = curr_node.left_child();               // vertical ray shut down
          }
          continue;
        }
        else //if p_cv was given
        {
          bool is_equal_to_he_min = traits->equal_curve_end_2_object()
                               (Curve_end(*p_cv,ARR_MIN_END), ce);
          bool is_equal_to_he_max = traits->equal_curve_end_2_object()
                               (Curve_end(*p_cv,ARR_MAX_END), ce);

          CGAL_assertion( is_equal_to_he_min || is_equal_to_he_max );
          CGAL_USE(is_equal_to_he_max);

          curr_node = is_equal_to_he_min ? curr_node.right_child() : curr_node.left_child();

          continue;
        }
      }
      else
      {
        CGAL_assertion(is_end_point_left_low(ce, curr_node) ||
                       is_end_point_right_top(ce, curr_node) ||
                       are_equal_end_points(ce, curr_node));
          //(is_fict_vtx &&
          //            (is_end_point_left_low(ce,*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item))) ||
          //             is_end_point_right_top(ce,*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item))) ||
          //             traits->equal_curve_end_2_object()(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)),ce))) ||
          //           (!is_fict_vtx &&
          //            (is_end_point_left_low(ce,boost::apply_visitor(point_for_vertex_visitor(), curr_item)) ||
          //             is_end_point_right_top(ce,boost::apply_visitor(point_for_vertex_visitor(), curr_item)) ||
          //             traits->equal_curve_end_2_object()(boost::apply_visitor(point_for_vertex_visitor(), curr_item),ce))));
        return Locate_type();
      }
    }
    if (traits->is_td_edge(curr_item))
    { // if curr_item represents an edge,
      //   so top() is a real Halfedge with a curve() if curr_item is active
      //   or curr_item holds the curve if it is not active

      const X_monotone_curve_2& he_cv = *(boost::apply_visitor(cv_for_edge_visitor(), curr_item));
      Comparison_result cres = traits->compare_curve_end_y_at_x_2_object()(ce, he_cv);
      if (cres == SMALLER)
      {
        curr_node = curr_node.left_child();
        continue;
      }
      else if (cres == LARGER)
      {
        curr_node = curr_node.right_child();
        continue;
      }
      else
      {
        // ce is on the CURVE (he_cv = ce.cv()) itself
        CGAL_warning(
          (cres == EQUAL) &&
          (traits->compare_curve_end_x_2_object()
                        (ce, Curve_end(he_cv,ARR_MAX_END)) != LARGER) &&
          (traits->compare_curve_end_x_2_object()
                        (ce, Curve_end(he_cv,ARR_MIN_END)) != SMALLER));

        if (!p_cv) //if p_cv was not given
        {
          // For a vertical curve, we always visit it after visiting
          // one of its endpoints.
          if ((up == EQUAL) || traits->is_vertical(curr_item))
          {
            if (traits->is_active(curr_item))
              return CURVE;
            curr_node = curr_node.left_child();
          }
          else if (up == LARGER)
          {
            curr_node = curr_node.right_child();
          }
          else
          { // if (up==SMALLER)
            curr_node = curr_node.left_child();
          }
          continue;
        }
        else //if p_cv was given
        {
          Comparison_result res = EQUAL;

          if ((traits->parameter_space_in_x_2_object()
                      (ce.cv(), ce.ce()) == ARR_INTERIOR) &&
              (traits->parameter_space_in_y_2_object()
                      (ce.cv(), ce.ce()) == ARR_INTERIOR) )
          {
            //if ce is interior then there might be more than one curve
            // with ce as its endpoint
            bool is_min_equal = traits->equal_curve_end_2_object()
                                      (Curve_end(*p_cv,ARR_MIN_END),
                                       Curve_end(he_cv,ARR_MIN_END));

            bool is_max_equal = traits->equal_curve_end_2_object()
                                      (Curve_end(*p_cv,ARR_MAX_END),
                                       Curve_end(he_cv,ARR_MAX_END));

            CGAL_warning (is_min_equal || is_max_equal);
            CGAL_USE(is_max_equal);

            //the underlying point of ce
            const Point& p = (ce.ce() == ARR_MIN_END) ?
                           traits->construct_min_vertex_2_object()(ce.cv()) :
                           traits->construct_max_vertex_2_object()(ce.cv()) ;

            res = is_min_equal ?
                    traits->compare_cw_around_point_2_object()
                           (he_cv, is_curve_to_right(he_cv,p),
                            *p_cv, is_curve_to_right(*p_cv,p), p)  :
                    traits->compare_cw_around_point_2_object()
                           (*p_cv, is_curve_to_right(*p_cv,p),
                            he_cv, is_curve_to_right(he_cv,p), p ,false);
          }

          switch(res)
          {
           case LARGER:
            curr_node = curr_node.right_child();
          break;
           case SMALLER:
            curr_node = curr_node.left_child();
            break;
           case EQUAL:
            switch(up)
            {
             case LARGER:
              curr_node = curr_node.right_child();
              break;
             case SMALLER:
              curr_node = curr_node.left_child();
              break;
             case EQUAL:
              if (traits->is_active(curr_item))
                return CURVE;
              curr_node = curr_node.left_child();
              break;
            }
            break;
          }
        }
      }
    }
    else
    {
      // if is_degenerate() == 0, meaning: curr_item is a real trapezoid
      if (traits->is_active(curr_item))
      {
        Td_active_trapezoid& tr = boost::get<Td_active_trapezoid>(curr_item);
        return tr.is_on_boundaries() ? UNBOUNDED_TRAPEZOID : TRAPEZOID;
      }
      curr_node = curr_node.left_child();
      continue;
    }
  }
}






//-----------------------------------------------------------------------------
// Description:
//  advances input Data structure using data structure,input point ce and
//  possibly X_monotone_curve_2 p_cv till
//  ce is found(if p_cv hadn't been given)
//  p_cv is found(if p_cv was given)
//  or
//  leaf node reached
// postcondition:
//  output is the closest active trapezoid to ce/p_cv
// remark:
//  use this function with care!
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Locate_type
Trapezoidal_decomposition_2<Td_traits>
::search_using_dag_with_cv (Dag_node& curr_node,
                            const Traits* traits,
                            const Point& p,
                            const X_monotone_curve_2* p_cv,
                            Comparison_result up /*=EQUAL*/) const
{

  while(true)
  {
    //curr_node is the current pointer to node in the data structure
    //curr_item is the curent Td_map_item held in curr_node
    Td_map_item curr_item(curr_node.get_data());

    if (traits->is_td_vertex(curr_item))
    { // if the curr_item represents a vertex
      //bool is_fict_vtx = traits->is_fictitious_vertex(curr_item);
      //if ((is_fict_vtx  && is_end_point_right_top(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), p)) ||
      //    (!is_fict_vtx && is_end_point_right_top(boost::apply_visitor(point_for_vertex_visitor(), curr_item), p)) )
      if (is_end_point_left_low(p, curr_node))
      {
        curr_node = curr_node.left_child();
        continue;
      }
      //else if ((is_fict_vtx  && is_end_point_left_low(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), p)) ||
      //         (!is_fict_vtx && is_end_point_left_low(boost::apply_visitor(point_for_vertex_visitor(), curr_item), p)) )
      else if (is_end_point_right_top(p, curr_node))
      {
        curr_node = curr_node.right_child();
        continue;
      }
      //else if ((is_fict_vtx && traits->equal_curve_end_2_object()(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), p)) ||
      //         (!is_fict_vtx &&  traits->equal_2_object()(boost::apply_visitor(point_for_vertex_visitor(), curr_item), p)) )
      else if (are_equal_end_points(p, curr_node))
      {
        if (!p_cv) //if p_cv was not given
        {
          if ( up == EQUAL )
          {      // point found!
            if (traits->is_active(curr_item))
              return POINT;
            curr_node = curr_node.left_child();
          }
          else if ( up == LARGER )
          {          // vertical ray shut up
            curr_node = curr_node.right_child();
          }
          else
          { // if ( up == SMALLER )
            curr_node = curr_node.left_child();               // vertical ray shut down
          }
          continue;
        }
        else //if p_cv was given
        {
          bool is_equal_to_he_min = traits->equal_curve_end_2_object()
                               (Curve_end(*p_cv,ARR_MIN_END), p);
          bool is_equal_to_he_max = traits->equal_curve_end_2_object()
                               (Curve_end(*p_cv,ARR_MAX_END), p);

          CGAL_assertion( is_equal_to_he_min || is_equal_to_he_max );
          CGAL_USE(is_equal_to_he_max);

          curr_node = is_equal_to_he_min ? curr_node.right_child() : curr_node.left_child();

          continue;
        }
      }
      else
      {
        CGAL_assertion(is_end_point_left_low(p, curr_node) ||
                       is_end_point_right_top(p, curr_node) ||
                       are_equal_end_points(p, curr_node));
        //CGAL_assertion((is_fict_vtx &&
        //              (is_end_point_left_low(p,*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item))) ||
        //               is_end_point_right_top(p,*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item))) ||
        //               traits->equal_curve_end_2_object()(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)),p))) ||
        //             (!is_fict_vtx &&
        //              (is_end_point_left_low(p,boost::apply_visitor(point_for_vertex_visitor(), curr_item)) ||
        //               is_end_point_right_top(p,boost::apply_visitor(point_for_vertex_visitor(), curr_item)) ||
        //               traits->equal_2_object()(boost::apply_visitor(point_for_vertex_visitor(), curr_item),p))));
        return Locate_type();
      }
    }
    if (traits->is_td_edge(curr_item))
    { // if curr_item represents an edge,
      //   so top() is a real Halfedge with a curve() if curr_item is active
      //   or curr_item holds the curve if it is not active
      const X_monotone_curve_2& he_cv = *(boost::apply_visitor(cv_for_edge_visitor(), curr_item));
      Comparison_result cres = traits->compare_y_at_x_2_object()(p, he_cv);
      if (cres == SMALLER)
      {
        curr_node = curr_node.left_child();
        continue;
      }
      else if (cres == LARGER)
      {
        curr_node = curr_node.right_child();
        continue;
      }
      else
      {
        // p is on the CURVE itself
        CGAL_warning(
          (cres == EQUAL) &&
          (traits->compare_curve_end_x_2_object()
                        (p, Curve_end(he_cv,ARR_MAX_END)) != LARGER) &&
          (traits->compare_curve_end_x_2_object()
                        (p, Curve_end(he_cv,ARR_MIN_END)) != SMALLER));

        if (!p_cv) //if p_cv was not given
        {
          // For a vertical curve, we always visit it after visiting
          // one of its endpoints.
          if ((up == EQUAL) || traits->is_vertical(curr_item))
          {
            if (traits->is_active(curr_item))
              return CURVE;
            curr_node = curr_node.left_child();
          }
          else if (up == LARGER)
          {
            curr_node = curr_node.right_child();
          }
          else
          { // if (up==SMALLER)
            curr_node = curr_node.left_child();
          }
          continue;
        }
        else //if p_cv was given
        {
          Comparison_result res = EQUAL;

          //p is interior then there might be more than one curve
          // with p as its endpoint
          bool is_min_equal = traits->equal_curve_end_2_object()
                                    (Curve_end(*p_cv,ARR_MIN_END),
                                     Curve_end(he_cv,ARR_MIN_END));

          bool is_max_equal = traits->equal_curve_end_2_object()
                                    (Curve_end(*p_cv,ARR_MAX_END),
                                     Curve_end(he_cv,ARR_MAX_END));

          CGAL_warning (is_min_equal || is_max_equal);
          CGAL_USE(is_max_equal);

          res = is_min_equal ?
                  traits->compare_cw_around_point_2_object()
                         (he_cv, is_curve_to_right(he_cv,p),
                          *p_cv, is_curve_to_right(*p_cv,p), p)  :
                  traits->compare_cw_around_point_2_object()
                           (*p_cv, is_curve_to_right(*p_cv,p),
                            he_cv, is_curve_to_right(he_cv,p), p ,false);

          switch(res)
          {
           case LARGER:
            curr_node = curr_node.right_child();
          break;
           case SMALLER:
            curr_node = curr_node.left_child();
            break;
           case EQUAL:
            switch(up)
            {
             case LARGER:
              curr_node = curr_node.right_child();
              break;
             case SMALLER:
              curr_node = curr_node.left_child();
              break;
             case EQUAL:
              if (traits->is_active(curr_item))
                return CURVE;
              curr_node = curr_node.left_child();
              break;
            }
            break;
          }
        }
      }
    }
    else
    {
      // if is_degenerate() == 0, meaning: curr_item is a real trapezoid
      if (traits->is_active(curr_item))
      {
        Td_active_trapezoid& tr = boost::get<Td_active_trapezoid>(curr_item);
        return tr.is_on_boundaries() ? UNBOUNDED_TRAPEZOID : TRAPEZOID;
      }
      curr_node = curr_node.left_child();
      continue;
    }
  }
}





//-----------------------------------------------------------------------------
// Description:
//
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Dag_node
Trapezoidal_decomposition_2<Td_traits>
::container2dag (Nodes_map& ar, int left, int right,
                 int& num_of_new_nodes) const
{
  CGAL_warning(traits != NULL);

  if (right > left)
  {
    int d = (int)std::floor((double(right+left))/2);
    Dag_node& node(ar.find(d)->second);
    Td_map_item item(node.get_data());
    CGAL_assertion(traits->is_active(item));
    CGAL_assertion(traits->is_td_trapezoid(item));
    Dag_node& tr_node( ar.find(d)->second);
    Td_active_trapezoid& tr( boost::get<Td_active_trapezoid>(tr_node.get_data()));
    Vertex_const_handle v = tr.right();

    Curve_end ce(traits->vtx_to_ce(v));
    bool is_interior = traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce())
                      && traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce());

    Dag_node curr_node;
    Dag_node left_child  (container2dag(ar,left,d,num_of_new_nodes));
    Dag_node right_child (container2dag(ar,d+1,right,num_of_new_nodes));
    if (is_interior)
    {
      curr_node.replace(Td_map_item(Td_active_vertex(v,m_empty_he_handle)),left_child, right_child);
    }
    else
    {
      curr_node.replace(Td_map_item(Td_active_fictitious_vertex(v,m_empty_he_handle)),left_child, right_child);
    }

    num_of_new_nodes++;
    boost::apply_visitor(set_dag_node_visitor((Dag_node*)&(curr_node.left_child())),curr_node.left_child().get_data());
    boost::apply_visitor(set_dag_node_visitor((Dag_node*)&(curr_node.right_child())),curr_node.right_child().get_data());
    boost::apply_visitor(set_dag_node_visitor((Dag_node*)&curr_node),curr_node.get_data());
    //curr_node.left_child()->set_dag_node(&curr_node.left_child());
    //curr_node.right_child()->set_dag_node(&curr_node.right_child());
    //curr_node->set_dag_node(&curr_node);// fake temporary node
    deactivate_vertex(curr_node); //curr_node->remove(); // mark as deleted
    boost::apply_visitor(set_dag_node_visitor((Dag_node*)NULL),curr_node.get_data());//curr_node->set_dag_node(0);

    return curr_node;
  }
  else
  {
    return ar.find(left)->second;
  }
}

template <class Td_traits>
bool Trapezoidal_decomposition_2<Td_traits>
::is_last_edge(Halfedge_const_handle /* he */ , Td_map_item& vtx_item)
{
  CGAL_precondition(traits->is_td_vertex(vtx_item));
  CGAL_precondition(traits->is_active(vtx_item));

  Vertex_const_handle v (boost::apply_visitor(vertex_for_active_vertex_visitor(), vtx_item));

  typename Arrangement_on_surface_2::Halfedge_around_vertex_const_circulator first, second;
  first = second = v->incident_halfedges();
  ++second;

  if (second == first) //if he is the only halfedge around v -> return true
    return true;
  return false;
}

//-----------------------------------------------------------------------------
// Description:
//  if Halfedge or twin already inserted the latter is returned.
//  otherwise the left-low most edge-degenerate trapezoid that represents the
//  input Halfedge is returned
// Remark:
//  Given an edge-degenerate trapezoid representing a Halfedge,
//  all the other trapezoids representing the Halfedge can be extracted
//  via moving continously to the left and right neighbours.
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Td_map_item
Trapezoidal_decomposition_2<Td_traits>
::insert(Halfedge_const_handle he) //::insert_in_face_interior(Halfedge_const_handle he)
{
  //print_cv_data(he->curve());
  if (m_with_guarantees)
    update();

  // locate the input Halfedge end points in the Td_map_item Dag

  CGAL_assertion(traits != NULL);

  //get the two vertices of the halfedge
  Vertex_const_handle v1 = (he->direction() == ARR_LEFT_TO_RIGHT) ? he->source() : he->target();
  Vertex_const_handle v2 = (he->direction() == ARR_LEFT_TO_RIGHT) ? he->target() : he->source();

  //define the Curve end points (curve end = vertex)
  const Curve_end ce1(he, ARR_MIN_END); //MICHAL: to be removed?
  const Curve_end ce2(he, ARR_MAX_END); //MICHAL: to be removed?

  // make sure that the two endpoints  are not the same point
  CGAL_precondition(!traits->equal_curve_end_2_object()(ce1, ce2));

  Locate_type lt1,lt2;

  //should hold the trapezoids in which the edge endpoints should be located
  Td_map_item item1;
  Td_map_item item2;

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION

  locate_optimization(ce1,item1,lt1);

#else
  //location of the left endpoint of the edge we're inserting
  item1 = locate(ce1,lt1);

#endif


  //the inserted edge should not cut any existing edge
  if (lt1 == CURVE)
  {
    CGAL_precondition_msg(lt1 != CURVE, "Input is not planar as\
      one of the input point inside previously inserted Halfedge.");
    return Td_map_item(0);
  }

  //if the edge starts at vertex, we should not insert it into the DAG,
  //but we should update all the edges incident to the vertex.
  //else if this is a new vertex - insert a node to the DAG that will represent the new vertex.
  //the incident edges in this case is only the edge itself, and so it is a trivial operation.
  Td_map_item p1_item = (lt1 == POINT) ?
    update_vtx_with_new_edge(he,ce1,item1,lt1) :
    insert_curve_at_vtx_using_dag(he,v1,item1,lt1);


#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION

  locate_optimization(ce2,item2,lt2);
  locate_opt_empty();

#else
  // TODO(oren): locating the second endpoint. this is not necessary,
  // and time consuming.
  item2 = locate(ce2,lt2);

#endif

  if (lt2==CURVE)
  {
    CGAL_precondition_msg(lt2!=CURVE,"Input is not planar as\
      one of the input point inside previously inserted Halfedge.");
    return Td_map_item(0);
  }

  Td_map_item p2_item = (lt2 == POINT) ?
    update_vtx_with_new_edge(he,ce2,item2,lt2) :
    insert_curve_at_vtx_using_dag(he,v2,item2,lt2);

  // locate and insert end points of the input halfedge to the Td_map_item
  // Dag if needed
  Dag_node p1_node(*(boost::apply_visitor(dag_node_visitor(), p1_item)));
  //Dag_node p2_node(*(boost::apply_visitor(dag_node_visitor(), p2_item)));// Not in use

  // create the Td_map_item iterator for traveling along the Trapezoids that
  // intersect the input Halfedge, using left-low to right-high order
  In_face_iterator it = follow_curve(p1_node,he,LARGER);
  boost::optional<Td_active_trapezoid> curr_trp = boost::none;
  Td_map_item prev = p1_item;
  Td_map_item prev_bottom_tr = Td_map_item(0); //active_tr
  Td_map_item prev_top_tr = Td_map_item(0); //active_tr
  Td_map_item old_e = Td_map_item(0);
  Td_map_item old_top_tr = Td_map_item(0); //active_tr
  Td_map_item old_bottom_tr = Td_map_item(0); //active_tr

  Dag_node* node = 0;
  bool first_time = true;

  CGAL_assertion(!traits->is_empty_item(*it) && traits->is_td_trapezoid(*it));

  while(!!it) //this means as long as the iterator is valid
  {
    curr_trp = it.trp();
    CGAL_assertion(curr_trp != boost::none);
    prev_bottom_tr = (*curr_trp).lb();
    prev_top_tr = (*curr_trp).lt();

    // pass using it along cv
    it++;             //this is the logic of the iterator.
                      // the iterator goes to the next trapezoid right-high.
    node = (*curr_trp).dag_node();
    if(first_time)
    {
      Halfedge_const_handle top_he ((*curr_trp).top());
      if((top_he == he) || (top_he == he->twin())) {
        CGAL_warning((top_he != he) && (top_he != he->twin()));
        return Td_map_item(0);
      }
    }
    first_time = false;

    CGAL_assertion(node != NULL);
    split_trapezoid_by_halfedge (*node, old_e, old_bottom_tr, old_top_tr, he);



    if (node->is_inner_node())
    {
      // merge adjacent trapezoids on input Halfedge's bottom side if possible
      //   make sure we try to merge active trapezoids
      //node->left_child() is active
      Td_map_item new_btm_item = node->left_child().get_data();
      CGAL_assertion(traits->is_td_trapezoid(new_btm_item) && traits->is_active(new_btm_item));
      if(merge_if_possible(prev_bottom_tr, new_btm_item))
      {
        Dag_node* left_child_node = boost::apply_visitor(dag_node_visitor(),prev_bottom_tr);
        node->set_left_child(*left_child_node);
        old_bottom_tr = prev_bottom_tr;
        m_number_of_dag_nodes--; //update number of nodes in the DAG after merge
      }

      // merge adjacent trapezoids on input Halfedge's top side if possible
      //   make sure we try to merge active trapezoids
      //node->right_child() is active
      Td_map_item new_top_item = node->right_child().get_data();
      CGAL_assertion(traits->is_td_trapezoid(new_top_item) && traits->is_active(new_top_item));
      if(merge_if_possible(prev_top_tr, new_top_item))
      {
        Dag_node* right_child_node = boost::apply_visitor(dag_node_visitor(),
                                                          prev_top_tr);
        node->set_right_child(*right_child_node);
        old_top_tr = prev_top_tr;
        m_number_of_dag_nodes--; //update number of DAG nodes after merge
      }
      // update trapezoid's left/right neighbouring relations
      //MICHAL: if the assertion below fails then we need to check why
      CGAL_assertion(!traits->is_td_trapezoid(prev));
      if(traits->is_td_trapezoid(prev))
      {
        //MICHAL: if we reach here ->then this need to be uncommented
        //        I thought that prev is always a degenerate point
        //curr_trp->set_lb(prev);
        //curr_trp->set_lt(prev);
        //prev->set_rb(curr_trp);
        //prev->set_rt(curr_trp);
      }
    }
    else
    {

#ifdef CGAL_TD_DEBUG

      CGAL_assertion((*curr_trp).is_valid(traits));

#endif

      break;
    }

  }
  m_number_of_curves++;

#ifdef CGAL_TD_DEBUG
  write(std::cout,*m_dag_root,*traits) << std::endl;
  std::cout << "\nTD::insert() exited with data structure"
            << is_valid(*m_dag_root) << std::endl;
#endif

  return old_e;
}


//-----------------------------------------------------------------------------
// Description:
//
// Remark:
//  Assumes the map to be planar.
template <class Td_traits>
void Trapezoidal_decomposition_2<Td_traits>
::remove(Halfedge_const_handle he)
{
  //print_dag_addresses(*m_dag_root);

  if (m_with_guarantees)
    update();

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
  locate_opt_empty();
#endif

  CGAL_warning(traits != NULL);

  //calculating leftmost and rightmost curve ends of he
  const Curve_end leftmost(he,ARR_MIN_END);
  const Curve_end rightmost(he,ARR_MAX_END);

  //locating leftmost & rightmost curve ends
  Locate_type lt1,lt2;
  Td_map_item p1_item = locate(leftmost,lt1);
  Td_map_item p2_item = locate(rightmost,lt2);

  //both should be located on a point degenerate trapezoid
  CGAL_warning(lt1==POINT && lt2==POINT);

  if (lt1!=POINT || lt2!=POINT)
    return;

  CGAL_warning(boost::apply_visitor(dag_node_visitor(), p1_item) != NULL);
  CGAL_warning(boost::apply_visitor(dag_node_visitor(), p2_item) != NULL);

  //retrieve the Dag_nodes of the two point-degenerate trapezoid
  Dag_node& p1_node = *(boost::apply_visitor(dag_node_visitor(), p1_item));//*t1.dag_node();
  Dag_node& p2_node = *(boost::apply_visitor(dag_node_visitor(), p2_item));//*t2.dag_node();

  //calculate the immediate lower, central and upper neighbourhood of
  // the curve in the data structure
  //In_face_iterator btm_it(follow_curve(tt1,he,SMALLER));
  In_face_iterator btm_it(follow_curve(p1_node,he,SMALLER));
  In_face_iterator mid_it(follow_curve(p1_node,he,EQUAL));
  In_face_iterator top_it(follow_curve(p1_node,he,LARGER));

  bool inc_btm = true; //true if btm_it should be be incremented (not top_it)
  bool prev_inc_btm = false; //holds the inc_btm from previous iteration
  bool end_reached = false; //true if this is the last iteration

  //define the map of new DAG nodes of the new trapezoids
  Nodes_map new_array;
  int last_index[] = {0,0};
  int sz = 0;

  Td_active_trapezoid& btm_it_tr (btm_it.trp());
  Td_active_trapezoid& top_it_tr (top_it.trp());
  Vertex_const_handle left_v = btm_it_tr.left();
  Vertex_const_handle right_v;
  Td_map_item last_btm_tr_item = Td_map_item(0); //pointer to the last btm_it tr
  Td_map_item last_top_tr_item = Td_map_item(0); //pointer to the last top_it tr
  Td_map_item last_new_tr_item = Td_map_item(0); //last new trpz that was created
  Td_map_item old_tr_item      = Td_map_item(0); //old trpz on which the new is based


  CGAL_warning((top_it.trp()).left() == left_v);

  //-----------------------------------
  //1. remove adjacency at left end point
  //  first_cv_tr is the first trapezoid representing he (type TD_EDGE)
  Td_map_item first_edge_fragment_item (*mid_it);

  //-----------------------------------
  //2. update the map & the dag with new trapezoids which are merge of the
  //  trapezaoids above and below the removed halfedge.
  do {
    btm_it_tr = btm_it.trp();
    top_it_tr = top_it.trp();

    // decide which of btm_it,top_it to increment
    inc_btm = is_end_point_left_low(Curve_end(traits->vtx_to_ce(btm_it_tr.right())),
                                    Curve_end(traits->vtx_to_ce(top_it_tr.right())));
    // the current iterator that should be incremented
    In_face_iterator& curr_it =  inc_btm ? btm_it : top_it;
    Td_active_trapezoid& curr_it_tr (curr_it.trp());

    // reference to the last curr_it tr
    Td_map_item& last_tr_item (inc_btm ? last_btm_tr_item : last_top_tr_item);

    //set the new trpz right end
    right_v = curr_it_tr.right();

    // create a new trapezoid (the merge of top_it and btm_it)
    typename Nodes_map::value_type
      pair(sz,
           Dag_node(Td_active_trapezoid(left_v, right_v,
                                        btm_it_tr.bottom(),top_it_tr.top()),
                    (std::max)(btm_it_tr.dag_node()->depth(),
                              top_it_tr.dag_node()->depth())));
    new_array.insert(pair);
    //copy trapezoid data from btm and top trapezoids
    Dag_node& new_node = (new_array.find(sz))->second;
    ++sz;
    Td_active_trapezoid& tr(boost::get<Td_active_trapezoid>(new_node.get_data()));
    tr.set_dag_node(&new_node);
    tr.set_lb(btm_it_tr.lb());
    tr.set_lt(top_it_tr.lt());
    //new_node->set_dag_node(&new_node);
    //new_node->set_lb(btm_it->lb());
    //new_node->set_lt(top_it->lt());
    if (!traits->is_empty_item(last_new_tr_item))
    {
      if (traits->is_trpz_top_equal(last_new_tr_item, new_node.get_data()))
      {
        tr.set_lt(last_new_tr_item); //new_node->set_lt(last_new_tr);
      }

      if (traits->is_trpz_bottom_equal(last_new_tr_item, new_node.get_data()))
      {
        tr.set_lb(last_new_tr_item); //new_node->set_lb(last_new_tr);
      }
    }
    if (!traits->is_empty_item(tr.lb()))
    {
      Td_map_item lb_item = tr.lb();
      Td_active_trapezoid& lb_tr(boost::get<Td_active_trapezoid>(lb_item));
      lb_tr.set_rb(new_node.get_data());
    }
    if (!traits->is_empty_item(tr.lt()))
    {
      Td_map_item lt_item = tr.lt();
      Td_active_trapezoid& lt_tr(boost::get<Td_active_trapezoid>(lt_item));
      lt_tr.set_rt(new_node.get_data());
    }

    last_new_tr_item = new_node.get_data(); //new_node.operator->(); //get the last new trapezoid created

    //update arguments for next iteration:
    left_v = right_v; //new left is curr right
    last_btm_tr_item = *btm_it;
    last_top_tr_item = *top_it;

    old_tr_item = *curr_it;  //the old trpz on which the new is based
    curr_it++; //increment the iterator to the next trapezoid
    //curr_it_tr  = *(curr_it.trp());
    end_reached = !btm_it || !top_it;

    //copy neighbouring trapezoids in case top/btm are not the same for the old
    //  trapezoid and the next trapezoid after incrementing the old one
    if (!btm_it ||
        (inc_btm && !traits->is_trpz_bottom_equal(old_tr_item,*curr_it)))
    {
      Td_map_item rb (boost::apply_visitor(rb_visitor(),old_tr_item));
      if (!traits->is_empty_item(rb))
      {
        boost::apply_visitor(set_lb_visitor(last_new_tr_item),rb); //rb->set_lb(last_new_tr);
        boost::apply_visitor(set_rb_visitor(rb),last_new_tr_item); //last_new_tr->set_rb(rb);
      }
    }
    if (!top_it ||
        (!inc_btm && !traits->is_trpz_top_equal(old_tr_item,*curr_it)))
    {
      Td_map_item rt (boost::apply_visitor(rt_visitor(),old_tr_item));
      if (!traits->is_empty_item(rt))
      {
        boost::apply_visitor(set_lt_visitor(last_new_tr_item),rt); //rt->set_lt(last_new_tr);
        boost::apply_visitor(set_rt_visitor(rt),last_new_tr_item); //last_new_tr->set_rt(rt);
      }
    }

    //set the no longer relevant trapezoids as removed and add the new nodes
    //  as their replacement

    Dag_node* last_tr_node (boost::apply_visitor(dag_node_visitor(), last_tr_item));

    if (prev_inc_btm != inc_btm)
    {
      int num_of_new_nodes = 0;
      Dag_node tmp =
        container2dag (new_array, last_index[inc_btm ? 0 : 1],
                        sz-1, num_of_new_nodes);

      deactivate_trapezoid( *last_tr_node, &tmp); //last_tr->remove(&tmp);
      m_number_of_dag_nodes += num_of_new_nodes; //new vertex nodes (rooted at tmp) were added
      m_number_of_dag_nodes += 1; //new node (tmp) was added

      last_index[inc_btm ? 0 : 1] = sz;
      prev_inc_btm = inc_btm; //update for next iteration
      //tmp is the root of a sub graph.
      //The largest depth in this sub-graph may be the largest leaf depth
      update_largest_leaf_depth( tmp.max_depth() );

    }
    else
    {
      int num_of_new_nodes = 0;
      Dag_node tmp = container2dag (new_array, sz-1, sz-1, num_of_new_nodes);

      deactivate_trapezoid( *last_tr_node, &tmp); //last_tr->remove(&tmp);
      m_number_of_dag_nodes += 1; //new node (tmp) was added

      last_index[inc_btm ? 0 : 1] = sz;
      //tmp is a node with no children
      update_largest_leaf_depth( tmp.max_depth() );
    }

    // update the dag node pointer in the trapezoid
    const Dag_node* real = &last_tr_node->left_child();
    boost::apply_visitor(set_dag_node_visitor((Dag_node*)real),real->get_data()); //(*real)->set_dag_node((Dag_node*)real);
  }
  while(!end_reached);


  // get the iterator (btm_it or top_it) that holds the trapezoid that was
  //  not removed in the last iteration
  In_face_iterator& it = !prev_inc_btm ? btm_it : top_it;
  Td_active_trapezoid& tr (it.trp());
  // remove the last trapezoid to remove and set the last new trapezoid
  //  created as its replacement. update the relevant data
  Td_map_item rb = tr.rb();
  Td_map_item rt = tr.rt();

  int num_of_new_nodes = 0;
  Dag_node tmp =
    container2dag(new_array, last_index[!inc_btm ? 0 : 1],
                  static_cast<int>(new_array.size())-1, num_of_new_nodes);

  deactivate_trapezoid( *(tr.dag_node()), &tmp);  //it->remove(&tmp);
  //tmp is the root of a sub graph.
  //The largest depth in this sub-graph may be the largest leaf depth
  update_largest_leaf_depth( tmp.depth() );
  m_number_of_dag_nodes += num_of_new_nodes; //new node (tmp) was added

  const Dag_node* real = &tr.dag_node()->left_child();
  boost::apply_visitor(set_dag_node_visitor((Dag_node*)real),real->get_data()); //(*real)->set_dag_node((Dag_node*)real);

  if (!traits->is_empty_item(rb))
  {
    boost::apply_visitor(set_rb_visitor(rb),last_new_tr_item); //last_new_tr->set_rb(rb);
    boost::apply_visitor(set_lb_visitor(last_new_tr_item),rb); //rb->set_lb(last_new_tr);
  }
  if (!traits->is_empty_item(rt))
  {
    boost::apply_visitor(set_rt_visitor(rt),last_new_tr_item); //last_new_tr->set_rt(rt);
    boost::apply_visitor(set_lt_visitor(last_new_tr_item),rt); //rt->set_lt(last_new_tr);
  }

  //-----------------------------------
  //3. remove the trapezoids that represent the removed halfedge
  boost::shared_ptr<X_monotone_curve_2> removed_cv_ptr (new X_monotone_curve_2(he->curve()));
  Base_map_item_iterator last_edge_fragment_it = mid_it;//Base_trapezoid_iterator last_mid = mid_it;
  Dag_node* e_node = NULL;
  while(!!++mid_it)
  {
    e_node = boost::apply_visitor(dag_node_visitor(),*last_edge_fragment_it);
    deactivate_edge(removed_cv_ptr,*e_node); //last_mid->remove();
    last_edge_fragment_it = mid_it;
  }

  //-----------------------------------
  //4. remove adjacency at right end point

  //remove the final trapezoid representing the removed halfedge
  e_node = boost::apply_visitor(dag_node_visitor(),*last_edge_fragment_it);
  deactivate_edge(removed_cv_ptr,*e_node); //last_mid->remove();

  //-----------------------------------
  //5. if the halfedge vertices are now isolated, undo the split trapezoid
  //  by point(vtx) operation
  if (is_last_edge(he ,p1_item))
    undo_split_trapezoid_by_vertex (p1_node, leftmost);

  if (is_last_edge(he ,p2_item))
    undo_split_trapezoid_by_vertex (p2_node, rightmost);

  //-----------------------------------
  //6. reevaluating number of curves
  m_number_of_curves--;
}


//-----------------------------------------------------------------------------
// Description:
//
// preconditions:
//  p is not on an edge or a vertex.
template <class Td_traits>
typename Trapezoidal_decomposition_2<Td_traits>::Td_map_item&
Trapezoidal_decomposition_2<Td_traits>
::vertical_ray_shoot(const Point & p,Locate_type & lt,
                     const bool up_direction /*=true*/) const
{
#ifdef CGAL_TD_DEBUG
  CGAL_assertion(traits);
#endif

  // We replace the following locate with a direct call to
  // search_using_dag because we need to deal
  // with cases where the source of shoot is a point/curve.
  // reference t_p = locate(p,lt);

  Dag_node curr = *m_dag_root;

  lt = search_using_dag(curr, traits, p, m_empty_he_handle ,
                                  up_direction ?
                                  CGAL::LARGER : CGAL::SMALLER);

  Td_map_item& item(curr.get_data());

  if (traits->is_td_trapezoid(item))
  {
    CGAL_assertion(traits->is_active(item));
    /* using exact traits, it may happen that p is on the
       right side of the trapezoid directly under its
       right point(analogouly directly above its left point).
       with the trapezoid extending to the left.
       In this case vertical ray shoot upwards(downwards)
       doesn't returns c as output.

       Example.
       x---x
       p
       x------x
    */
    Td_active_trapezoid& tr (boost::get<Td_active_trapezoid>(item));

    if ((up_direction && !tr.is_on_right_boundary() &&
        (traits->compare_curve_end_x_2_object()
          (p,traits->vtx_to_ce(tr.right())) == EQUAL) &&
        (tr.is_on_left_boundary() ||
          !traits->equal_curve_end_2_object()(traits->vtx_to_ce(tr.left()),
                                              traits->vtx_to_ce(tr.right())))) ||
        (!up_direction && !tr.is_on_left_boundary() &&
        (traits->compare_curve_end_x_2_object()
          (p,traits->vtx_to_ce(tr.left())) == EQUAL) &&
        (tr.is_on_right_boundary() ||
          !traits->equal_curve_end_2_object()(traits->vtx_to_ce(tr.left()),
                                              traits->vtx_to_ce(tr.right())))))
    {
      // recalculate vertical ray shoot using locate on point
      return up_direction ?
          locate(traits->vtx_to_ce(tr.right()),lt) : locate(traits->vtx_to_ce(tr.left()),lt);
    }

    if (up_direction ? tr.is_on_top_boundary() : tr.is_on_bottom_boundary())
      lt = UNBOUNDED_TRAPEZOID;
    else
      lt = TRAPEZOID;

  }
  return item;
}

//MICHAL: commented due to inefficient depth update, remove and insert instead
//-----------------------------------------------------------------------------
//// Description:
////
//template <class Td_traits>
//void Trapezoidal_decomposition_2<Td_traits>
//::before_split_edge(const X_monotone_curve_2& cv,
//                   const X_monotone_curve_2& cv1,
//                   const X_monotone_curve_2& cv2)
//{
//#ifdef MICHAL_DEBUG
//  std::cout << "SPLITTING: --------------------------" << std::endl;
//  std::cout << "cv before split" << std::endl;
//  print_cv_data(cv);
//  std::cout << "cv1 before split" << std::endl;
//  print_cv_data(cv1);
//  std::cout << "cv2 before split" << std::endl;
//  print_cv_data(cv2);
//#endif
//
//#ifdef CGAL_TD_DEBUG
//  std::cout << "\nTD::before_split_edge(" << cv << "," << cv1 << "," << cv2
//            << ") called with "
//            << (is_valid(*m_dag_root) ? "valid" : "invalid")
//            << " data structure" <<  std::endl;
//  write(std::cout,*m_dag_root,*traits) << std::endl;
//#endif
//
//  if (m_with_guarantees) update();
//
//#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
//
//  locate_opt_empty();
//
//#endif
//
//#ifndef CGAL_TD_DEBUG
//
//  if (!traits)
//  {
//    CGAL_warning(traits != NULL);
//    return;
//  }
//  if (!traits->are_mergeable_2_object()(cv1,cv2))
//  {
//    CGAL_warning(traits->are_mergeable_2_object()(cv1,cv2));
//    return;
//  }
//
//#else
//
//  if (!traits->are_mergeable_2_object()(cv1,cv2))
//  {
//    std::cerr << "\ncv " << cv;
//    std::cerr << "\ncv1 " << cv1;
//    std::cerr << "\ncv1 " << cv2 << std::endl;
//  }
//  CGAL_precondition(traits != NULL);
//  CGAL_precondition(traits->are_mergeable_2_object()(cv1,cv2));
//
//#endif
//
//  // find the splitting point (curve end)
//  Curve_end split_ce =
//      traits->equal_curve_end_2_object()(Curve_end(cv1, ARR_MAX_END),
//                                        Curve_end(cv2, ARR_MIN_END) ) ?
//      Curve_end(cv1, ARR_MAX_END) : Curve_end(cv2, ARR_MAX_END);
//
//
//
//#ifndef CGAL_TD_DEBUG
//
//  CGAL_warning( is_end_point_left_low(Curve_end(cv,ARR_MIN_END),split_ce) );
//
//  CGAL_warning( is_end_point_right_top(Curve_end(cv,ARR_MAX_END),split_ce) );
//
//#else
//
//  CGAL_precondition( is_end_point_left_low(Curve_end(cv,ARR_MIN_END),split_ce) );
//
//  CGAL_precondition( is_end_point_right_top(Curve_end(cv,ARR_MAX_END),split_ce) );
//
//#endif
//
//  // find extremal points
//  const Curve_end leftmost = (traits->equal_curve_end_2_object()
//                                (Curve_end(cv1, ARR_MAX_END),
//                                 Curve_end(cv2, ARR_MIN_END) ))?
//                            Curve_end(cv1,ARR_MIN_END) :
//                            Curve_end(cv2,ARR_MIN_END) ;
//
//  const Curve_end rightmost = (traits->equal_curve_end_2_object()
//                                 (Curve_end(cv1, ARR_MAX_END),
//                                  Curve_end(cv2, ARR_MIN_END) ))?
//                            Curve_end(cv2,ARR_MAX_END) :
//                            Curve_end(cv1,ARR_MAX_END) ;
//
//  CGAL_assertion(traits->equal_curve_end_2_object()
//                              (Curve_end(cv,ARR_MIN_END), leftmost));
//  CGAL_assertion(traits->equal_curve_end_2_object()
//                              (Curve_end(cv,ARR_MAX_END), rightmost));
//
//  //locate the trapezoids of the extremal points
//  Locate_type lt1,lt2;
//
//  // representing trapezoids for extremal points
//  X_trapezoid& t1 = locate(leftmost,lt1);
//  X_trapezoid& t2 = locate(rightmost,lt2);
//
//#ifndef CGAL_TD_DEBUG
//
//  CGAL_warning(lt1==POINT && lt2==POINT);
//  CGAL_warning(t1.is_active() && t2.is_active());
//
//#else
//
//  CGAL_precondition(lt1==POINT && lt2==POINT);
//  CGAL_precondition(t1.is_active() && t2.is_active());
//  CGAL_warning(t1.dag_node() != NULL);
//  CGAL_warning(t2.dag_node() != NULL);
//
//#endif
//  m_before_split.m_cv_before_split = cv;
//  //set iterators for below curve, on curve & above curve
//  Dag_node& tt1 = *t1.dag_node();
//  m_before_split.m_p_btm_it = new In_face_iterator(follow_curve(tt1,m_before_split.m_cv_before_split,SMALLER));
//  m_before_split.m_p_mid_it = new In_face_iterator(follow_curve(tt1,m_before_split.m_cv_before_split,EQUAL));
//  m_before_split.m_p_top_it = new In_face_iterator(follow_curve(tt1,m_before_split.m_cv_before_split,LARGER));
//  //locate the splitting point in the trapezoidal map
//  //  should be found on a degenerate trapezoid representing a curve
//  Locate_type lt;
//  X_trapezoid& old_t = locate(split_ce,lt);
//
//#ifdef CGAL_TD_DEBUG
//
//  CGAL_assertion(lt==CURVE);
//  CGAL_precondition(old_t.is_active());
//  CGAL_warning(old_t.dag_node());
//
//#endif
//  //the DAG node of the curve trapezoid where the spiltting point is
//  Dag_node& old_split_node = *old_t.dag_node();
//
//  CGAL_assertion(traits->equal_curve_end_2_object()
//                  (traits->vtx_to_ce(old_t.left()),leftmost));
//
//  CGAL_assertion(traits->equal_curve_end_2_object()
//                  (traits->vtx_to_ce(old_t.right()),rightmost));
//
//
//
//  m_before_split.m_p_old_t = &old_t;
//  m_before_split.m_p_t1 = &t1;
//  m_before_split.m_p_t2 = &t2;
//
//}

//-----------------------------------------------------------------------------
// Description:
// Input:
//  1 whole curves
//  2 partial halfedge_handle-s
// precondition:
//  The two halfedges are valid
//  The first input curve is the union of the two halfedges.
//  The intersection of the latter is a point inside the
//  interior of the former.
//  The latter are ordered from left-down to right-up
// postcondition:
//  The first input curve is broken into two halfedges
//  corresponding to the input.




//MICHAL: commented due to inefficient depth update, remove and insert instead
//template <class Td_traits>
//void Trapezoidal_decomposition_2<Td_traits>
//::split_edge(const X_monotone_curve_2& cv,Halfedge_const_handle he1,
//               Halfedge_const_handle he2)
//{
//  //make sure both halfedges are valid
//  CGAL_precondition_code(Halfedge_const_handle invalid_he);
//  CGAL_precondition( he1 != invalid_he);
//  CGAL_precondition( he2 != invalid_he);
//
//#ifdef CGAL_TD_DEBUG
//  std::cout << "\nTD::split_edge(" << cv << "," << he1 << "," << he2
//            << ") called with " << (is_valid(*m_dag_root) ? "valid" : "invalid")
//            << " data structure" <<  std::endl;
//  write(std::cout,*m_dag_root,*traits) << std::endl;
//#endif
//
//
//  // find the splitting point (vertex & curve end)
//  Vertex_const_handle split_v =
//      traits->equal_curve_end_2_object()(Curve_end(he1, ARR_MAX_END),
//                                        Curve_end(he2, ARR_MIN_END) ) ?
//      ((he1->direction()== ARR_LEFT_TO_RIGHT) ?
//        he1->target() : he1->source())
//     :((he2->direction()== ARR_LEFT_TO_RIGHT) ?
//        he2->target() : he2->source());
//
//  Curve_end ce(traits->vtx_to_ce(split_v));
//
//  // find extremal points
//  const Curve_end leftmost = (traits->equal_curve_end_2_object()
//                                (Curve_end(he1, ARR_MAX_END),
//                                 Curve_end(he2, ARR_MIN_END) ))?
//                            Curve_end(he1,ARR_MIN_END) :
//                            Curve_end(he2,ARR_MIN_END) ;
//
//  const Curve_end rightmost = (traits->equal_curve_end_2_object()
//                                 (Curve_end(he1, ARR_MAX_END),
//                                  Curve_end(he2, ARR_MIN_END) ))?
//                            Curve_end(he2,ARR_MAX_END) :
//                            Curve_end(he1,ARR_MAX_END) ;
//
//
//  //MICHAL: new
//  X_trapezoid& t1 = *m_before_split.m_p_t1;
//  X_trapezoid& t2 = *m_before_split.m_p_t2;
//  X_trapezoid& old_t = *m_before_split.m_p_old_t;
//  In_face_iterator& bottom_it = *m_before_split.m_p_btm_it;
//  In_face_iterator& mid_it = *m_before_split.m_p_mid_it;
//  In_face_iterator& top_it = *m_before_split.m_p_top_it;
//  //MICHAL: new end
//
//  //the DAG node of the curve trapezoid where the spiltting point is
//  Dag_node& old_split_node = *old_t.dag_node();
//
//
//  // previous left and right sons of this DAG node
//  const Dag_node& old_left = old_split_node.left_child();
//  const Dag_node& old_right= old_split_node.right_child();
//
//  //define the left halfedge and the right halfedge, according
//  //  to the splitting point
//  //Halfedge_const_handle* p_left_he  = NULL;
//  //Halfedge_const_handle* p_right_he = NULL;
//  Halfedge_const_handle left_he  = he2;
//  Halfedge_const_handle right_he = he1;
//
//  if (traits->equal_curve_end_2_object() (Curve_end(he2,ARR_MIN_END), ce))
//  {
//    left_he  = he1; //p_left_he  = &he1;
//    right_he = he2; //p_right_he = &he2;
//  }
//  //else
//  //{
//  //  p_left_he  = &he2;
//  //  p_right_he = &he1;
//  //}
//
//  CGAL_assertion(old_t.is_active());
//
//  //updating the data structure:
//  //  the cv trpz where the splitting point is located will hold the point
//  //  its new left will hold the left cv (left and right children as before)
//  //  its new right will hold the right cv (left and right children asa before)
//
//  //defining the new left child node
//  const Dag_node& new_left_node =
//          Dag_node(X_trapezoid(old_t.left(), split_v,
//                               left_he, left_he, //*p_left_he, *p_left_he,
//                               X_trapezoid::TD_EDGE),
//                   old_left, old_right);
//
//  //defining the new right child node
//  const Dag_node& new_right_node =
//          Dag_node(X_trapezoid(split_v, old_t.right(),
//                               right_he, right_he, //*p_right_he, *p_right_he,
//                               X_trapezoid::TD_EDGE),
//                   old_left, old_right);
//
//  //defining the new parent node (which is the splitting point):
//
//  //split_v is the ARR_MAX_END of left_te, and ARR_MIN_END of right_te.
//  //need to send the halfedge whose source is split_v.
//  Halfedge_const_handle btm_he = left_he;
//  if (btm_he->direction() == ARR_LEFT_TO_RIGHT)
//    btm_he = btm_he->twin();
//  //Halfedge_const_handle* p_btm = p_left_he;
//  //if ((*p_btm)->direction() == ARR_LEFT_TO_RIGHT)
//  //  p_btm = &(*p_btm)->twin();
//
//  Halfedge_const_handle top_he = right_he;
//  if (top_he->direction() == ARR_RIGHT_TO_LEFT)
//    top_he = top_he->twin();
//  //Halfedge_const_handle* p_top = p_right_he;
//  //if ((*p_top)->direction() == ARR_RIGHT_TO_LEFT)
//  //  p_top = &(*p_top)->twin();
//
//
//  const Dag_node& new_pnt_node =
//          Dag_node(X_trapezoid(split_v, split_v, btm_he, top_he, //*p_btm, *p_top,
//                               X_trapezoid::TD_VERTEX),
//                   new_left_node, new_right_node);
//
//  X_trapezoid& new_left_t = *new_left_node;
//  X_trapezoid& new_right_t= *new_right_node;
//  X_trapezoid& new_t      = *new_pnt_node;
//
//  // locate trapezoid trees that correspond to the closest
//  //   trapezoids above and below ce
//  X_trapezoid* left_top_t   = top_it.operator->();
//  X_trapezoid* left_bottom_t= bottom_it.operator->();
//
//  while(is_end_point_left_low(traits->vtx_to_ce(left_top_t->right()),ce))
//    left_top_t = left_top_t->rb();
//
//  while(is_end_point_left_low(traits->vtx_to_ce(left_bottom_t->right()),ce))
//    left_bottom_t = left_bottom_t->rt();
//
//  Dag_node left_top    = *left_top_t->dag_node();
//  Dag_node left_bottom = *left_bottom_t->dag_node();
//
//  //replace the old curve cv with the new curves in the leftmost
//  //  and rightmost end points.
//  //the curve end ce belongs to cv interior
//  set_trp_params_after_split_halfedge_update (left_he , t1, he1, he2);
//  set_trp_params_after_split_halfedge_update (right_he, t2, he1, he2);
//
//  //set the point's lb() which is:
//  //     the first halfedge adjacent to the point emanating toward up
//  //     or right sweeped using a counter clockwise sweep
//  //     starting from up direction not including.
//  //set the point's rt() which is:
//  //     the first halfedge adjacent to the point emanating toward bottom
//  //     or left sweeped using a counter clockwise sweep
//  //     starting from bottom direction not including.
//  new_t.set_rt (&new_left_t);
//  new_t.set_lb (&new_right_t);
//
//  //set lb() and rt(0 of the two halfedges trapezoids
//  // rt() is the next clockwise degenerate_curve around
//  //      rightmost end_point (possibly himself)
//  // lb() is the next clockwise degenerate_curve around
//  //      leftmost end_point (possibly himself)
//  new_left_t.set_lb ((old_t.lb() != &old_t) ? old_t.lb() : &new_left_t);
//  new_left_t.set_rt (&new_right_t);
//  new_right_t.set_lb(&new_left_t);
//  new_right_t.set_rt((old_t.rt() != &old_t)? old_t.rt() : &new_right_t);
//
//  // update geometric boundary for trapezoids representing cv
//  // first update the trapezoids of the new left curve
//  X_trapezoid* prev = 0;
//  while(*mid_it != old_t)
//  {
//    mid_it->set_top(left_he);
//    mid_it->set_bottom(left_he);
//    mid_it->set_right(split_v);
//
//    //MICHAL: added this assertion just to be sure:
//    CGAL_assertion(mid_it->type() == X_trapezoid::TD_EDGE); //MICHAL: won't work, I removed the type() method
//
//    //MICHAL: who is prev?
//    prev = mid_it.operator->();
//    mid_it++;
//  }
//  if (prev)
//  {
//    prev->set_rb(&new_left_t); //MICHAL: I don't understand this
//  }
//  else // new_left_t is leftmost representative for he
//  {
//    //set_neighbours_after_split_halfedge_update (new_left_t, t1, he1, he2); //MICHAL: this method does nothing
//  }
//  if (t1.rt()==&old_t) t1.set_rt(&new_left_t);
//  if (t1.lb()==&old_t) t1.set_lb(&new_left_t);
//  mid_it++;
//  new_right_t.set_rb(mid_it.operator->()); //MICHAL: what does it do?
//
//  prev = 0;
//  while(!!mid_it)
//  {
//    mid_it->set_top(right_he);
//    mid_it->set_bottom(right_he);
//    mid_it->set_left(split_v);
//
//    //MICHAL: added this assertion just to be sure:
//    CGAL_assertion(mid_it->type() == X_trapezoid::TD_EDGE);
//
//    prev = mid_it.operator->();
//    mid_it++;
//  }
//  if (prev)
//  {
//    new_right_t.set_rb(old_t.rb()); //MICHAL: I don't understand this
//  }
//  else // new_right_t is rightmost representative for te
//  {
//    //set_neighbours_after_split_halfedge_update (new_right_t,t2,he1, he2,false); //MICHAL: this method does nothing
//  }
//  if (t2.rt()==&old_t) t2.set_rt(&new_right_t);
//  if (t2.lb()==&old_t) t2.set_lb(&new_right_t);
//
//  // update geometric boundary for trapezoids below te
//  while (*bottom_it != *left_bottom)
//  {
//
//#ifdef CGAL_TD_DEBUG
//
//    CGAL_assertion (traits->equal_2_object()(bottom_it->top()->curve(), cv));
//
//#endif
//    //MICHAL: added this assertion to see if it fails (if we reach an edge_end)
//    CGAL_assertion (bottom_it->type() == X_trapezoid::TD_TRAPEZOID);
//
//    bottom_it->set_top(left_he);
//    bottom_it++;
//  }
//
//#ifdef CGAL_TD_DEBUG
//
//  CGAL_assertion (*bottom_it==*left_bottom);
//
//#endif
//
//  Dag_node& bottom_tt = *bottom_it->dag_node();
//  bottom_it++;
//
//#ifdef CGAL_TD_DEBUG
//
//  CGAL_assertion(traits->is_in_closure (*bottom_tt, ce));
//
//#endif
//
//  split_trapezoid_by_vertex (bottom_tt, split_v, btm_he, top_he); //*p_btm, *p_top);
//  // set the splitting trapezoid to be the same one that splits the
//  // X_curve'like trapezoid
//  *bottom_tt = new_t;
//  // update top curves
//  bottom_tt.left_child()->set_top(left_he);
//  bottom_tt.right_child()->set_top(right_he);
//  // left and right are not neighbours.
//  bottom_tt.left_child()->set_rt(0);
//  bottom_tt.right_child()->set_lt(0);
//
//
//  while(!!bottom_it)
//  {
//
//#ifdef CGAL_TD_DEBUG
//
//    CGAL_assertion(traits->equal_2_object() (bottom_it->top()->curve(), cv));
//
//#endif
//    //MICHAL: added this assertion to see if it fails (if we reach an edge_end)
//    CGAL_assertion(bottom_it->type() == X_trapezoid::TD_TRAPEZOID);
//
//    bottom_it->set_top(right_he);
//    bottom_it++;
//  }
//
//  // update geometric boundary for trapezoids above cv
//  while (*top_it != *left_top)
//  {
//
//#ifdef CGAL_TD_DEBUG
//
//    CGAL_assertion(traits->equal_2_object() (top_it->bottom()->curve(), cv));
//
//#endif
//    //MICHAL: added this assertion to see if it fails (if we reach an edge_end)
//    CGAL_assertion(top_it->type() == X_trapezoid::TD_TRAPEZOID);
//
//    top_it->set_bottom(left_he);
//    top_it++;
//  }
//
//#ifdef CGAL_TD_DEBUG
//
//  CGAL_assertion(*top_it == *left_top);
//
//#endif
//
//  Dag_node &top_tt = *top_it->dag_node();
//  top_it++;
//
//#ifdef CGAL_TD_DEBUG
//
//  CGAL_assertion(traits->is_in_closure (*top_tt, ce));
//
//#endif
//
//  split_trapezoid_by_vertex (top_tt, split_v, btm_he, top_he);// left_he, right_he);
//  // set the splitting trapezoid to be the same one that splits the
//  // X_curve'like trapezoid
//  *top_tt = new_t;
//  // update bottom side
//  top_tt.left_child()->set_bottom(left_he);
//  top_tt.right_child()->set_bottom(right_he);
//  // left and right aren't neighbours
//  top_tt.left_child()->set_rb(0);
//  top_tt.right_child()->set_lb(0);
//
//  while(!!top_it)
//  {
//
//#ifndef CGAL_TD_DEBUG
//
//    CGAL_warning(traits->equal_2_object()(top_it->bottom()->curve(), cv));
//
//#else
//
//    if (!traits->equal_2_object()(top_it->bottom()->curve(), cv))
//      std::cout << "\ntop_it->bottom()->curve() "<< top_it->bottom()->curve()
//                << "\t cv= " << cv;
//    CGAL_assertion(traits->equal_2_object()(top_it->bottom()->curve() ,cv));
//
//#endif
//    //MICHAL: added this assertion to see if it fails (if we reach an edge_end)
//    CGAL_assertion(top_it->type() == X_trapezoid::TD_TRAPEZOID);
//
//    top_it->set_bottom(right_he);
//    top_it++;
//  }
//
//  // mark inactive trapezoids
//  //  depth of new_pnt_node is updated here (in the remove operation)
//  //   and also depth of the sub-DAG rooted at it
//  old_t.remove((Dag_node*)&new_pnt_node);
//  update_largest_leaf_depth( new_pnt_node.max_depth() ); //MICHAL: this is a recursive call for the sub-DAG --EXPENSIVE!
//  //adding nodes for the splitting-point and the two parts of the split curve
//  m_number_of_dag_nodes += 3;
//  old_t.set_curve_for_rem_he(m_before_split.m_cv_before_split); //MICHAL: added this so the trpz will hold the original curve before the split
//
//  const Dag_node* p_new       = &old_t.dag_node()->left_child();
//  const Dag_node* p_new_left  = &p_new->left_child();
//  const Dag_node* p_new_right = &p_new->right_child();
//  const Dag_node* p_old_left  = &p_new_left->left_child();
//  const Dag_node* p_old_right = &p_new_left->right_child();
//
//  (*p_new)->set_dag_node ((Dag_node*)p_new);
//  (*p_new_left)->set_dag_node ((Dag_node*)p_new_left);
//  (*p_new_right)->set_dag_node ((Dag_node*)p_new_right);
//  (*p_old_left)->set_dag_node ((Dag_node*)p_old_left);
//  (*p_old_right)->set_dag_node ((Dag_node*)p_old_right);
//
//#ifdef CGAL_TD_DEBUG
//
//  CGAL_assertion (old_split_node->is_valid(traits));
//  CGAL_assertion (new_pnt_node->is_valid(traits));
//  CGAL_assertion ((*p_new)->is_valid(traits));
//  CGAL_assertion ((*p_new_left)->is_valid(traits));
//  CGAL_assertion ((*p_new_right)->is_valid(traits));
//  CGAL_assertion ((*p_old_left)->is_valid(traits));
//  CGAL_assertion ((*p_old_right)->is_valid(traits));
//  CGAL_assertion (top_tt->is_valid(traits));
//  CGAL_assertion (bottom_tt->is_valid(traits));
//  CGAL_assertion (old_left->is_valid(traits));
//  CGAL_assertion (old_right->is_valid(traits));
//  CGAL_assertion (traits->is_degenerate_point(**p_new));
//  CGAL_assertion (traits->is_degenerate_curve(**p_new_left));
//  CGAL_assertion (traits->is_degenerate_curve(**p_new_right));
//  CGAL_assertion (traits->equal_curve_end_2_object()
//                  (Curve_end((*p_new_right)->bottom(), ARR_MIN_END),
//                   (*p_new)->right()) );
//  CGAL_assertion (traits->equal_curve_end_2_object()
//                  (Curve_end((*p_new_left)->top(), ARR_MAX_END),
//                   (*p_new)->left()) );
//#endif
//
//  // reevaluating number of curves
//  m_number_of_curves++;
//
//#ifdef CGAL_TD_DEBUG
//  std::cout << "\nTD::split_edge() exited with data structure"
//            << is_valid(*m_dag_root) << std::endl;
//  write(std::cout,*m_dag_root,*traits) << std::endl;
//#endif
//}
//



//-----------------------------------------------------------------------------
// Description:
//
template <class Td_traits>
void Trapezoidal_decomposition_2<Td_traits>
::merge_edge (Halfedge_const_handle he1,
              Halfedge_const_handle he2,
              const X_monotone_curve_2& cv)
{
  //make sure the halfedge is valid
  CGAL_precondition_code(Halfedge_const_handle invalid_he);
  CGAL_precondition( he1 != invalid_he);
  CGAL_precondition( he2 != invalid_he);

#ifdef CGAL_TD_DEBUG
  std::cout << "\nTD::merge_edge(" << he1->curve() << "," << he2->curve()
            << "," << cv
            << ") called with " << (is_valid(*m_dag_root) ? "valid" : "invalid")
            << " data structure" <<  std::endl;
  write(std::cout,*m_dag_root,*traits) << std::endl;
#endif

  if (m_with_guarantees)
    update();

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION

  locate_opt_empty();

#endif

  const X_monotone_curve_2& cv1 = he1->curve();
  const X_monotone_curve_2& cv2 = he2->curve();

  if (!traits)
  {
    CGAL_warning(traits != NULL);
    return;
  }
  if (!traits->are_mergeable_2_object() (cv1, cv2))
  {
    CGAL_warning(traits->are_mergeable_2_object() (cv1, cv2));
    return;
  }

  // Calculate the common/merged point (Curve_end) of cv1 and cv2.
  // There should be one!
  Curve_end ce = traits->equal_curve_end_2_object()
                    (Curve_end(cv1,ARR_MAX_END),Curve_end(cv2,ARR_MIN_END)) ?
                  Curve_end(cv1,ARR_MAX_END) :
    // [-- cv1 -->] p [-- cv2 -->] or [<-- cv2 --] p [<-- cv1 --]
                 traits->equal_curve_end_2_object()
                    (Curve_end(cv1,ARR_MIN_END),Curve_end(cv2,ARR_MAX_END)) ?
    // [<-- cv1 --] p [<-- cv2 --] or [-- cv2 -->] p [-- cv1 -->]
                  Curve_end(cv1,ARR_MIN_END) : //
                 traits->equal_curve_end_2_object()
                    (Curve_end(cv1,ARR_MIN_END),Curve_end(cv2,ARR_MIN_END)) ?
    // [<-- cv1 --] p [-- cv2 -->]
                  Curve_end(cv1,ARR_MIN_END) :
    // [-- cv1 -->] p [<-- cv2 --]
                  Curve_end(cv1,ARR_MAX_END);

  //find the halfedge that will contain the merged curve
  // [<-- cv1 --] p [<-- cv2 --] or [<-- cv1 --] p [-- cv2 -->]-> he1->twin()
  // [-- cv1 -->] p [-- cv2 -->] or [-- cv1 -->] p [<-- cv2 --]-> he1
  //  Notice the curve cv is not yet updated
  Halfedge_const_handle merged_he =
            (traits->equal_curve_end_2_object()
                (Curve_end(cv1, ARR_MIN_END), Curve_end(cv2, ARR_MAX_END)) ||
             traits->equal_curve_end_2_object()
                (Curve_end(cv1, ARR_MIN_END), Curve_end(cv2, ARR_MIN_END)) )
             ? he1->twin() :  he1;

#ifdef CGAL_TD_DEBUG
  // ce is interior to the union curve
  CGAL_precondition(
        is_end_point_left_low (Curve_end(cv, ARR_MIN_END), ce));
  CGAL_precondition(
        is_end_point_right_top (Curve_end(cv, ARR_MAX_END), ce));

#endif

  //get the leftmost & rightmost Curve_end-s

  Curve_end leftmost  (cv, ARR_MIN_END);
  Curve_end rightmost (cv, ARR_MAX_END);

  //locate the leftmost, rightmost and the merged point in the data structure
  Locate_type lt1,lt2,lt;
  Td_map_item leftp_item  = locate (leftmost, lt1);
  Td_map_item rightp_item = locate (rightmost, lt2);
  Td_map_item mrgp_item   = locate (ce, lt);


  //varifying that all trapezoids are not NULL and are of type POINT
  CGAL_warning (boost::apply_visitor(dag_node_visitor(),leftp_item) != NULL);
  CGAL_warning (boost::apply_visitor(dag_node_visitor(),rightp_item)!= NULL);
  CGAL_warning (boost::apply_visitor(dag_node_visitor(),mrgp_item)  != NULL);



  //define the left curve and the right curve, according
  //  to the common point (that is merged)
  const X_monotone_curve_2* p_left_cv  = &cv2;
  const X_monotone_curve_2* p_right_cv = &cv1;
  Halfedge_const_handle left_he = he2;
  Halfedge_const_handle right_he = he1;
  if (traits->equal_curve_end_2_object() (Curve_end (cv2, ARR_MIN_END), ce))
  {
    p_left_cv  = &cv1;
    p_right_cv = &cv2;
    left_he = he1;
    right_he = he2;
  }

#ifdef CGAL_TD_DEBUG

  //CGAL_assertion (traits->equal_curve_end_2_object()
  //                (leftp_item.left(), leftmost));
  //CGAL_assertion (traits->equal_curve_end_2_object()
  //                (rightp_item.left(), rightmost));

  CGAL_assertion (is_end_point_left_low(leftmost, ce));
  CGAL_assertion (is_end_point_left_low(ce, rightmost));
  //compare left cv min with leftmost
  CGAL_assertion (traits->equal_curve_end_2_object()
                  (Curve_end(*p_left_cv, ARR_MIN_END), leftmost));
  //compare left cv max with ce
  CGAL_assertion (traits->equal_curve_end_2_object()
                  (Curve_end(*p_left_cv, ARR_MAX_END), ce));
  //compare right cv min with ce
  CGAL_assertion (traits->equal_curve_end_2_object()
                  (Curve_end(*p_right_cv, ARR_MIN_END), ce));
  //compare right cv max with rightmost
  CGAL_assertion (traits->equal_curve_end_2_object()
                  (Curve_end(*p_right_cv, ARR_MAX_END), rightmost));

#endif

  //get the nodes of leftmost point and merge point
  Dag_node& leftp_node  = *(boost::apply_visitor(dag_node_visitor(),leftp_item));
  Dag_node& mrgp_node = *(boost::apply_visitor(dag_node_visitor(),mrgp_item));

  //set iterators for below left curve, on left curve & above left curve
  In_face_iterator
    btm_left_it (follow_curve (leftp_node, *p_left_cv, SMALLER)),
    mid_left_it (follow_curve (leftp_node, *p_left_cv, EQUAL)),
    top_left_it (follow_curve (leftp_node, *p_left_cv, LARGER));

  //set iterators for below right curve, on right curve & above right curve
  In_face_iterator
    btm_right_it (follow_curve (mrgp_node, *p_right_cv, SMALLER)),
    mid_right_it (follow_curve (mrgp_node, *p_right_cv, EQUAL)),
    top_right_it (follow_curve (mrgp_node, *p_right_cv, LARGER));


#ifdef CGAL_TD_DEBUG

  CGAL_assertion (btm_left_it.operator->());
  CGAL_assertion (mid_left_it.operator->());
  CGAL_assertion (top_left_it.operator->());
  CGAL_assertion (btm_right_it.operator->());
  CGAL_assertion (mid_right_it.operator->());
  CGAL_assertion (top_right_it.operator->());
  CGAL_assertion (traits->is_active(*btm_left_it.operator->()));
  CGAL_assertion (traits->is_active(*mid_left_it.operator->()));
  CGAL_assertion (traits->is_active(*top_left_it.operator->()));
  CGAL_assertion (traits->is_active(*btm_right_it.operator->()));
  CGAL_assertion (traits->is_active(*mid_right_it.operator->()));
  CGAL_assertion (traits->is_active(*top_right_it.operator->()));

#endif

  //replacing old curves with the new merged halfedge

  Td_map_item  left_cv_fraction_item = *mid_left_it;
  Td_map_item  on_cv_left            = Td_map_item(0);
  Td_map_item  on_cv_right           = *mid_right_it;
  Td_map_item  above_cv_left         = Td_map_item(0);
  Td_map_item  above_cv_right        = *top_right_it;
  Td_map_item  below_cv_left         = Td_map_item(0);
  Td_map_item  below_cv_right        = *btm_right_it;
  Td_map_item  right_cv_fraction_item= Td_map_item(0);
  Td_map_item  dummy1                = Td_map_item(0);
  Td_map_item  dummy2                = Td_map_item(0);

  Vertex_const_handle leftmost_v  =
    (left_he->direction() == ARR_LEFT_TO_RIGHT) ? left_he->source() :
                                                  left_he->target();
  Vertex_const_handle rightmost_v =
    (right_he->direction() == ARR_LEFT_TO_RIGHT) ? right_he->target() :
                                                   right_he->source();
  //replacing the given curve with a new Halfedge_handle along the trapezoids
  // starting at the iterator, until the end (last parameter) is reached.
  // updating the last param as the last updated trapzoid
  update_map_items_after_merge (mid_left_it, left_he, merged_he,
                                leftmost_v, rightmost_v, on_cv_left);
  update_map_items_after_merge (mid_right_it, right_he, merged_he,
                                leftmost_v, rightmost_v, right_cv_fraction_item);
  update_map_items_after_merge (top_left_it, left_he, merged_he,
                                leftmost_v, rightmost_v, above_cv_left);
  update_map_items_after_merge (top_right_it, right_he, merged_he,
                                leftmost_v, rightmost_v, dummy1);
  update_map_items_after_merge (btm_left_it, left_he, merged_he,
                                leftmost_v, rightmost_v, below_cv_left);
  update_map_items_after_merge (btm_right_it, right_he, merged_he,
                                leftmost_v, rightmost_v, dummy2);


  // merge trapezoids that were split by the upward and downward
  // vertical extensions from ce (the merged point)

  // make sure only active trapezoids are merged
  CGAL_assertion( traits->is_active(above_cv_left) && traits->is_active(above_cv_right) );
  CGAL_assertion( traits->is_active(below_cv_left) && traits->is_active(below_cv_right) );

  merge_if_possible (above_cv_left, above_cv_right);
  merge_if_possible (below_cv_left, below_cv_right);

  // mark older trapezoids as inactive - nodes depth are updated here
  Dag_node* above_cv_right_node (boost::apply_visitor(dag_node_visitor(),above_cv_right));
  Dag_node* above_cv_left_node  (boost::apply_visitor(dag_node_visitor(),above_cv_left));
  deactivate_trapezoid( *above_cv_right_node, above_cv_left_node); //above_cv_right->remove(above_cv_left->dag_node());
  Dag_node* below_cv_right_node (boost::apply_visitor(dag_node_visitor(),below_cv_right));
  Dag_node* below_cv_left_node  (boost::apply_visitor(dag_node_visitor(),below_cv_left));
  deactivate_trapezoid( *below_cv_right_node, below_cv_left_node); //below_cv_right->remove(below_cv_left->dag_node());
  update_largest_leaf_depth((std::max)(above_cv_left_node->depth(),
                                       below_cv_left_node->depth()));
  //no need to update m_number_of_dag_nodes because the number of nodes did not change.


#ifdef CGAL_TD_DEBUG

  CGAL_warning (!traits->is_empty_item(on_cv_left));
  CGAL_warning (!traits->is_empty_item(on_cv_right));

#endif

  // make the merged point's representative inactive
  deactivate_vertex(mrgp_node); //mrgp_node->remove();

  CGAL_assertion(traits->is_td_edge(on_cv_left) && traits->is_active(on_cv_left));

  Td_active_edge& e_left (boost::get<Td_active_edge>(on_cv_left));
  e_left.set_next(on_cv_right);

  CGAL_assertion(traits->is_td_edge(on_cv_right) && traits->is_active(on_cv_right));

  //replacing the curve in the end points' trapezoids themselves (updating top/ bottom)
  update_vtx_cw_he_after_merge (*p_left_cv, merged_he, leftp_item);
  update_vtx_cw_he_after_merge (*p_right_cv, merged_he, rightp_item); //MICHAL: maybe I should pass the he1 & he2?

#ifdef CGAL_TD_DEBUG

  CGAL_warning(!traits->is_empty_item(left_cv_fraction_item));
  CGAL_warning(!traits->is_empty_item(right_cv_fraction_item));

#endif

  // reevaluating number of curves
  m_number_of_curves--;

#ifdef CGAL_TD_DEBUG
  std::cout << "\nTD::merge_edge() exited with data structure"
            << is_valid(*m_dag_root) << std::endl;
  write(std::cout,*m_dag_root,*traits) << std::endl;
#endif

}


//-----------------------------------------------------------------------------
// Description:
//
template <class Td_traits>
unsigned long
Trapezoidal_decomposition_2<Td_traits>
::longest_query_path_length_rec(bool minus_inf, Dag_node& min_node,
                                bool plus_inf, Dag_node& max_node,
                                Dag_node& node)
{
  //if NULL
  if (node.is_null())
    return 0;

  //if not valid range or empty return 0
  if (!minus_inf && !plus_inf)
  {
    Td_map_item min_node_item(min_node.get_data());
    if (traits->is_fictitious_vertex(min_node_item))
    {
      const Curve_end min_ce(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),min_node_item)));

      Td_map_item max_node_item(max_node.get_data());
      if (traits->is_fictitious_vertex(max_node_item))
      {
        const Curve_end max_ce(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),max_node_item)));

        //min-fict, max-fict
        if (!is_end_point_left_low(min_ce, max_ce))
          return 0;
      }
      else
      {
        const Point& max_p(boost::apply_visitor(point_for_vertex_visitor(),max_node_item));

        //min-fict, max-pt
        //if smaller than the point represented by min_node
        if (!is_end_point_left_low(min_ce, max_p))
          return 0;
      }
    }
    else
    {
      const Point& min_p(boost::apply_visitor(point_for_vertex_visitor(),min_node_item));

      Td_map_item max_node_item(max_node.get_data());
      if (traits->is_fictitious_vertex(max_node_item))
      {
        const Curve_end max_ce(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),max_node_item)));

        //min-pt, max-fict
        if (!is_end_point_left_low(min_p, max_ce))
          return 0;
      }
      else
      {
        const Point& max_p(boost::apply_visitor(point_for_vertex_visitor(),max_node_item));

        //min-pt, max-pt
        //if smaller than the point represented by min_node
        if (!is_end_point_left_low(min_p, max_p))
          return 0;
      }
    }
  }


  //if node represents a trapezoid
  if (traits->is_td_trapezoid(node.get_data()) )
    return 1;

  //if node represents a curve
  if (!traits->is_td_vertex(node.get_data()) )
    return (1 + (std::max)(
                  longest_query_path_length_rec(minus_inf, min_node,
                                                plus_inf, max_node,
                                                node.left_child()) ,
                  longest_query_path_length_rec(minus_inf, min_node,
                                                plus_inf, max_node,
                                                node.right_child()) ));

  //if node represents a vertex
  Td_map_item curr_item(node.get_data());
  //bool is_fict_vtx = traits->is_fictitious_vertex(curr_item);

  //if node is smaller than min, use min in the next recursion min bound, otherwise use node
  Dag_node new_min_node = node;
  if (!minus_inf)
  {
    Td_map_item min_node_item(min_node.get_data());
    if (traits->is_fictitious_vertex(min_node_item))
    {
      const Curve_end min_ce(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),min_node_item)));
      //if smaller than the point represented by min_node
      //if ((is_fict_vtx  && is_end_point_left_low(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), min_ce)) ||
      //    (!is_fict_vtx && is_end_point_left_low(boost::apply_visitor(point_for_vertex_visitor(), curr_item), min_ce) ))
      if (is_end_point_right_top(min_ce, node))
      {
        new_min_node = min_node;
      }
    }
    else
    {
      const Point& min_p(boost::apply_visitor(point_for_vertex_visitor(),min_node_item));
      //if smaller than the point represented by min_node
      //if ((is_fict_vtx  && is_end_point_left_low(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), min_p)) ||
      //    (!is_fict_vtx && is_end_point_left_low(boost::apply_visitor(point_for_vertex_visitor(), curr_item), min_p) ))
      if (is_end_point_right_top(min_p, node))
      {
        new_min_node = min_node;
      }
    }
  }

  //if node is larger than max, use max in the next recursion max bound, otherwise use node
  Dag_node new_max_node = node;
  if (!plus_inf)
  {
    Td_map_item max_node_item(max_node.get_data());
    if (traits->is_fictitious_vertex(max_node_item))
    {
      const Curve_end max_ce(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),max_node_item)));
      //if larger than the point represented by max_node
      //if ((is_fict_vtx  && is_end_point_right_top(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), max_ce)) ||
      //    (!is_fict_vtx && is_end_point_right_top(boost::apply_visitor(point_for_vertex_visitor(), curr_item), max_ce) ))
      if (is_end_point_left_low(max_ce, node))
      {
        new_max_node = max_node;
      }
    }
    else
    {
      const Point& max_p(boost::apply_visitor(point_for_vertex_visitor(),max_node_item));
      //if smaller than the point represented by min_node
      //if ((is_fict_vtx  && is_end_point_right_top(*(boost::apply_visitor(curve_end_for_fict_vertex_visitor(),curr_item)), max_p)) ||
      //    (!is_fict_vtx && is_end_point_right_top(boost::apply_visitor(point_for_vertex_visitor(), curr_item), max_p) ))
      if (is_end_point_left_low(max_p, node))
      {
        new_max_node = max_node;
      }
    }
  }
  //o/w continue with updated parameters
  return (1 + (std::max)(
                  longest_query_path_length_rec(minus_inf, min_node,
                                                false, new_max_node,
                                                node.left_child()) ,
                  longest_query_path_length_rec(false, new_min_node,
                                                plus_inf, max_node,
                                                node.right_child()) ));

}



} //namespace CGAL

#endif
