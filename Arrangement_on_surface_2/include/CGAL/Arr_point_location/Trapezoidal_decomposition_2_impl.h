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
// $URL: svn+ssh://balasmic@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_on_surface_2/include/CGAL/Arr_point_location/Trapezoidal_decomposition_2_impl.h $
// $Id: Arr_trapezoid_ric_pl_impl.h 56667 2010-06-09 07:37:13Z sloriot $
// 
//
// Author(s)     : Michal Balas   <balasmic@post.tau.ac.il>
//                based on the methods implemented in Trapezoidal_decomposition.h by Oren Nechushtan

#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_FUNCTIONS_H
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_FUNCTIONS_H

/*! \file
* Member-function definitions for the Trapezoidal_decomposition_2<Traits>
* class.
*/

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
::split_trapezoid_by_vertex(Dag_node& tt,
                           Vertex_const_handle v,
                           Halfedge_const_handle he_bottom_ray_shoot,
                           Halfedge_const_handle he_top_ray_shoot)
{ 
#ifndef CGAL_TD_DEBUG
    
  CGAL_warning(!tt.is_null());
  if (tt.is_null())  return tt;
    
#else
    
  CGAL_precondition(!tt.is_null());
    
#endif
    
  X_trapezoid& curr = *tt;
  X_trapezoid* lb = curr.lb();
  X_trapezoid* lt = curr.lt();
  X_trapezoid* rb = curr.rb();
  X_trapezoid* rt = curr.rt();
  
#ifndef CGAL_TD_DEBUG
  CGAL_warning(curr.is_active());
  CGAL_warning(traits->is_in_closure(curr, v->curve_end()));
#else
  CGAL_precondition(curr.is_active());
  if (!traits->is_in_closure(curr,v))
  {
    std::cout << "\ncurr=";
    write(std::cout,curr,*traits) << "\tv=" << v;
  }
  CGAL_precondition(traits->is_in_closure(curr,v));
#endif
  
  // left and right are set to the point itself,
  // bottom and top are set to the ray shooting resulting curves at this
  // stage.
  //need to define the boundaries flag before creating the trapezoid
  X_trapezoid sep(v,v,he_bottom_ray_shoot,he_top_ray_shoot,
                  X_trapezoid::TD_VERTEX, 
                  build_boundaries_flag( v->curve_end()));

  Dag_node leftDS(X_trapezoid
                   (curr.left(), v, curr.bottom(), curr.top(),
                    X_trapezoid::TD_TRAPEZOID,
                    curr.on_boundaries_flag() &
                    (CGAL_TD_ON_LEFT_BOUNDARY |
                     CGAL_TD_ON_BOTTOM_BOUNDARY |
                     CGAL_TD_ON_TOP_BOUNDARY)));

  Dag_node rightDS(X_trapezoid
                   (v, curr.right(), curr.bottom(),curr.top(),
                    X_trapezoid::TD_TRAPEZOID, 
                    curr.on_boundaries_flag() &
                    (CGAL_TD_ON_RIGHT_BOUNDARY |
                     CGAL_TD_ON_BOTTOM_BOUNDARY |
                     CGAL_TD_ON_TOP_BOUNDARY)));
  
  
  X_trapezoid& left = *leftDS;
  X_trapezoid& right = *rightDS;
  
#ifndef CGAL_TD_DEBUG
  CGAL_warning(traits->is_trpz_top_equal(left,right));
  CGAL_warning(traits->is_trpz_bottom_equal(left,right));
  CGAL_warning(left.is_on_left_boundary() == curr.is_on_left_boundary());
  CGAL_warning(right.is_on_right_boundary() == curr.is_on_right_boundary());
#else
  CGAL_warning(traits->is_trpz_top_equal(left,right));
  CGAL_warning(traits->is_trpz_bottom_equal(left,right));
  CGAL_assertion(left.is_on_left_boundary() == curr.is_on_left_boundary());
  CGAL_assertion(right.is_on_right_boundary()== curr.is_on_right_boundary());
#endif
  
  if (!traits->is_degenerate_curve(curr))
  {
    left.init_neighbours(lb,lt,&right,&right);
    right.init_neighbours(&left,&left,rb,rt);
    if (lb) lb->set_rb(&left);
    if (lt) lt->set_rt(&left);
    if (rb) rb->set_lb(&right);
    if (rt) rt->set_lt(&right);
  }
  else //if the trapezoid curr is an edge
  {
    left.set_type(CGAL_TD_EDGE);
    left.set_bottom(he_bottom_ray_shoot);
    left.set_top(he_bottom_ray_shoot);
    right.set_type(CGAL_TD_EDGE);
    right.set_bottom(he_top_ray_shoot);
    right.set_top(he_top_ray_shoot);
    left.set_rt(&right);
    left.set_lb(lb);
    left.set_rb(0);
    right.set_lb(&left);
    right.set_rt(rt);
    right.set_rb(rb);
  }

  tt.replace( sep, leftDS, rightDS); //nodes depth are updated here
  update_largest_leaf_depth( std::max(leftDS.depth(), rightDS.depth()) );
  m_number_of_dag_nodes += 2;
  
  const Dag_node* leftPtr  = &tt.left_child();
  const Dag_node* rightPtr = &tt.right_child();

  (*leftPtr)->set_dag_node((Dag_node*)leftPtr);
  (*rightPtr)->set_dag_node((Dag_node*)rightPtr);
  
#ifdef CGAL_TD_DEBUG    
  CGAL_assertion( &left  == leftDS.operator->() );
  CGAL_assertion( &right == rightDS.operator->());
  CGAL_assertion( left   == *leftDS );
  CGAL_assertion( right  == *rightDS);
  CGAL_assertion( left   == *tt.left() );
  CGAL_assertion( right  == *tt.right());
 
  CGAL_assertion(**left.dag_node() == left);
  CGAL_assertion(**right.dag_node() == right);    
#endif
  
  return tt;
}


//-----------------------------------------------------------------------------
// Description:
//  the opposite operation for spliting the trapezoid with 
//  vertical line through ce 
// Precondition:
//  The root trapezoid is degenerate point (ce) and is active 
template <class Td_traits>
void Trapezoidal_decomposition_2<Td_traits>
::undo_split_trapezoid_by_vertex(Dag_node& tr_node, const Curve_end& ce)
{
    
#ifndef CGAL_TD_DEBUG
    
  if (tr_node.is_null()||
      !tr_node->is_active()||
      !traits->is_degenerate_point(*tr_node)||
      !traits->equal_curve_end_2_object()( tr_node->left()->curve_end(),ce))
  {
    CGAL_warning(!tr_node.is_null());
    CGAL_warning(tr_node->is_active());
    CGAL_warning(traits->is_degenerate_point(*tr_node));
    CGAL_warning(traits->equal_curve_end_2_object()
                  (tr_node->left()->curve_end(),ce));
    return;
  }
  
#else
  
  CGAL_precondition(!tr_node.is_null());
  CGAL_precondition(tr_node->is_active());
  CGAL_precondition(traits->is_degenerate_point(*tr_node));
  CGAL_precondition(traits->equal_curve_end_2_object()
                    (tr_node->left()->curve_end(),ce));
  
#endif
  //MICHAL: maybe use pointers to Dag_node instead? - why copy?
  //get the ds left child and right child nodes of 
  //    tr_node (in the search structure)
  Dag_node tr_left_node  = tr_node.left_child(); //MICHAL:is it ok to add &?
  Dag_node tr_right_node = tr_node.right_child(); //MICHAL:is it ok to add &?
  
  //advances the ds nodes until ce is found 
  search_using_dag(tr_left_node,  traits, ce, 0);
  search_using_dag(tr_right_node, traits, ce, 0);

  //make sure the trapezoids are active before merging them
  CGAL_assertion(tr_left_node->is_active() && tr_right_node->is_active());

  bool mrg_res = false;
#ifndef CGAL_TD_DEBUG
  mrg_res = merge_if_possible( &*tr_left_node, &*tr_right_node );
  
  CGAL_warning(!tr_left_node.is_inner_node());
  CGAL_warning(!tr_right_node.is_inner_node());
  CGAL_warning(tr_left_node->is_on_right_boundary() ==
               tr_right_node->is_on_right_boundary());

#else
  std::cout << "\nremove_split_trapezoid_by_point(){";
  std::cout << "\ntr_left_node=";
  write(std::cout,*tr_left_node,*traits);
  std::cout << "\ntr_right_node=";
  write(std::cout,*tr_right_node,*traits) << "}" << std::endl;

  mrg_res = merge_if_possible( &*tr_left_node, &*tr_right_node);

  std::cout << "\n->";
  write(std::cout,*tr_left_node,*traits) << std::endl;

  CGAL_postcondition(mrg_res);
  CGAL_assertion(!tr_left_node.is_inner_node());
  CGAL_assertion(!tr_right_node.is_inner_node());
  CGAL_assertion(tr_left_node->is_on_right_boundary() ==
                  tr_right_node->is_on_right_boundary());
  CGAL_assertion(**tr_left_node->dag_node() == *tr_left_node);
#endif
  
  tr_right_node->remove(&tr_left_node);
  update_largest_leaf_depth( tr_left_node.depth()); //tr_left_node is not an inner node
  // mark root as deleted
  tr_node->remove();
  //no need to update m_number_of_dag_nodes because the number of nodes did not change.
  // removed nodes were only marked as removed
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
::split_trapezoid_by_halfedge(Dag_node& tt, X_trapezoid*& prev,
                           X_trapezoid*& prev_bottom, X_trapezoid*& prev_top, 
                           Halfedge_const_handle he)
{
  
#ifndef CGAL_TD_DEBUG
  
  CGAL_warning(traits != NULL);
  
#else
  
  CGAL_assertion(traits != NULL);
  
#endif
  

  X_trapezoid& currt = *tt;
  
  CGAL_assertion( currt.is_active() );
  CGAL_assertion( currt.type() == X_trapezoid::TD_TRAPEZOID);

  // sets left and right according to td_edge's source and target positions
  // sets bottom and top to Halfedge itself
  // no need to set the boundaries since its irrelevant for this trapezoid 
  //   type (TD_EDGE)
  X_trapezoid sep(he->min_vertex(),
                  he->max_vertex(),
                  he, he, X_trapezoid::TD_EDGE);

  
  //creates a one-way path for all the edge-degenerate
  //trapezoids that represent the Halfedge.
  //rb() is used to retrieve the
  //next on path information

  Dag_node topBT(X_trapezoid
                        (currt.left(), currt.right(), he, currt.top(),
                         X_trapezoid::TD_TRAPEZOID,
                         currt.on_boundaries_flag() &
                          (CGAL_TD_ON_LEFT_BOUNDARY  |
                           CGAL_TD_ON_RIGHT_BOUNDARY |
                           CGAL_TD_ON_TOP_BOUNDARY  )));
  Dag_node bottomBT(X_trapezoid
                        (currt.left(),currt.right(), currt.bottom(), he,
                         X_trapezoid::TD_TRAPEZOID,
                         currt.on_boundaries_flag() &
                          (CGAL_TD_ON_LEFT_BOUNDARY  |
                           CGAL_TD_ON_RIGHT_BOUNDARY |
                           CGAL_TD_ON_BOTTOM_BOUNDARY )));

  X_trapezoid& bottom = *bottomBT;
  X_trapezoid& top    = *topBT;

  top.init_neighbours(prev_top, currt.lt(), 0, currt.rt());
  bottom.init_neighbours(currt.lb(), prev_bottom, currt.rb(), 0);

  if (prev_bottom) prev_bottom->set_rt(&bottom);
  if (prev_top)   prev_top->set_rb(&top);
  if (currt.lb()) currt.lb()->set_rb(&bottom);
  if (currt.lt()) currt.lt()->set_rt(&top);
  if (currt.rb()) currt.rb()->set_lb(&bottom);
  if (currt.rt()) currt.rt()->set_lt(&top);

  tt.replace(sep,bottomBT,topBT); //nodes depth are updated here
  update_largest_leaf_depth( std::max(bottomBT.depth(), topBT.depth()) );
  m_number_of_dag_nodes += 2; //two new nodes were added to the DAG
  
  const Dag_node* bottomPtr = &tt.left_child();
  const Dag_node* topPtr    = &tt.right_child();

  (*bottomPtr)->set_dag_node((Dag_node*)bottomPtr);
  (*topPtr)->set_dag_node((Dag_node*)topPtr);

  if (prev) prev->set_rb(tt.operator->());

  //update these trapezoids pointers.
  // will be used for the next trapezoid that should be split
  //  by this Halfedge
  prev_bottom = (*bottomPtr).operator->();
  prev_top    = (*topPtr).operator->();
  prev        = tt.operator->();

  return tt;
}


#if 0
//-----------------------------------------------------------------------------
// Description:
//  replace halfedge-vtx adjacency in the data structure with a new one
// precondition:
//  the halfedge represented by he_tr is top-right
//  relative to the vertex represented by sep
//  if and only if he_top_right=true
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::replace_curve_at_point_using_geometry(X_trapezoid& he_tr, 
                                       const X_trapezoid& sep,
                                       bool he_top_right=true)
{
  //MICHAL: I am not sure when is this method called - maybe should be removed?
  Curve_end ce( sep.left()->curve_end());
  
  Around_point_circulator circ(traits,ce,he_top_right ? sep.rt() : sep.lb());

  if (circ.operator->())
  {
    //if the curve-end ce is on the boundaries, there is only one edge 
    //  starting/ending at it so no need to go on
    //otherwise:  the curve end is interior , so other curves around it need to be checked

    //if ce is interior
    if ((traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce()) 
                                                    == ARR_INTERIOR) &&
        (traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce()) 
                                                      == ARR_INTERIOR)  )
    {
      //the underlying point of ce
      const Point& p = (ce.ce() == ARR_MIN_END) ?
                        traits->construct_min_vertex_2_object()(ce.cv()) :
                        traits->construct_max_vertex_2_object()(ce.cv()) ;
      
      //MICHAL: I think that the of should be removed 
      //and both should use the same while because circ->top() should be equal to circ->bottom()
      
      if (he_top_right)
      {
        while(traits->compare_cw_around_point_2_object ()
                 (circ->top()->curve(),
                  is_edge_to_right(circ->top(),p),
                  he_tr.top()->curve(), 
                  is_edge_to_right(he_tr.top(),p), p) != EQUAL)
        {
          circ++;
        }
      } 
      else
      {
        while(traits->compare_cw_around_point_2_object()
                (circ->bottom()->curve(), 
                 is_edge_to_right(circ->bottom(),p),
                 he_tr.top()->curve(), 
                 is_edge_to_right(he_tr.top(),p), p, false) != EQUAL)
        {
          circ++;
        }
      }
    }
    circ.replace(he_tr); 
  }
}

#endif //if 0

//-----------------------------------------------------------------------------
// Description:
//  replace halfedge-vertex adjacency in the data structure with a new one
// precondition:
//  the halfedge represented by he_tr is top-right
//  relative to the vertex represented by sep if and only if top=true
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::set_neighbours_after_merge_halfedge_update(X_trapezoid& he_tr, 
                                          const X_trapezoid& sep,
                                          const X_monotone_curve_2& mrg_cv,
                                          bool he_top_right /*=true*/)
{
  CGAL_precondition(sep.is_active());
  Curve_end ce( sep.left()->curve_end() );
  
  Around_point_circulator circ(traits,ce,he_top_right ? sep.rt() : sep.lb());
  if (circ.operator->())
  {
    //if the curve-end ce is on the boundaries, there is only one edge 
    //  starting/ending at it so no need to go on
    //otherwise:  the curve end is interior , so other curves around it need to be checked

    //if ce is interior
    if ((traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce()) 
                                                    == ARR_INTERIOR) &&
        (traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce()) 
                                                      == ARR_INTERIOR)  )
    {
      //the underlying point of ce
      const Point& p = (ce.ce() == ARR_MIN_END) ?
                        traits->construct_min_vertex_2_object()(ce.cv()) :
                        traits->construct_max_vertex_2_object()(ce.cv()) ;
      while (circ->top() != he_tr.top() && circ->top() != he_tr.top()->twin() &&
             traits->compare_cw_around_point_2_object()
               (circ->top()->curve(),is_edge_to_right(circ->top(),p),
                mrg_cv, is_curve_to_right(mrg_cv,p),p, he_top_right) != EQUAL )
      {
        circ++;
      }
    }
    circ.replace(he_tr); 
  }
}


//-----------------------------------------------------------------------------
// Description:
//  replace halfedge-vertex adjacency in the data structure with a new one
// precondition:
//  the halfedge represented by he_tr is top-right
//  relative to the vertex represented by sep if and only if top=true
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::set_neighbours_after_split_halfedge_update(X_trapezoid& he_tr, 
                                          const X_trapezoid& sep,
                                          Halfedge_const_handle he1, 
                                          Halfedge_const_handle he2, 
                                          bool he_top_right /*=true*/)
{
  CGAL_precondition(sep.is_active());
  Curve_end ce( sep.left()->curve_end());
  
  Around_point_circulator circ(traits, ce, he_top_right ? sep.rt() : sep.lb());
  if (circ.operator->())
  {
    //if the curve-end ce is on the boundaries, there is only one edge 
    //  starting/ending at it so no need to go on
    //otherwise:  the curve end is interior , so other curves around it need to be checked

    //if ce is interior
    if ((traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce()) 
                                                    == ARR_INTERIOR) &&
        (traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce()) 
                                                      == ARR_INTERIOR)  )
    {
      //the underlying point of ce
      const Point& p = (ce.ce() == ARR_MIN_END) ?
                        traits->construct_min_vertex_2_object()(ce.cv()) :
                        traits->construct_max_vertex_2_object()(ce.cv()) ;

      while((circ->top() != he1) && (circ->top() != he1->twin()) && //MICHAL: he comp
            (circ->top() != he2) && (circ->top() != he2->twin()) && //MICHAL: he comp
             traits->compare_cw_around_point_2_object ()
               (circ->top()->curve(),
                is_edge_to_right(circ->top(),p),
                he_tr.top()->curve(), 
                is_edge_to_right(he_tr.top(),p), p, he_top_right) != EQUAL)
      {
        circ++;
      }
    }
    circ.replace(he_tr); 
  }
}

//-----------------------------------------------------------------------------
// Description:
//  
// precondition:
//  
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::set_neighbours_after_halfedge_insertion (X_trapezoid& he_tr,
                                           X_trapezoid& v_tr)
{
  
#ifndef CGAL_TD_DEBUG
  
  CGAL_warning(traits != NULL);
  CGAL_warning(traits->is_degenerate_point(v_tr));
  CGAL_warning(v_tr.is_active());
  CGAL_warning(traits->is_degenerate_curve(he_tr));
  CGAL_warning(he_tr.is_active());
  CGAL_warning(
    traits->equal_curve_end_2_object()
      (v_tr.left()->curve_end(), he_tr.right()->curve_end()) ||
    traits->equal_curve_end_2_object()
      (v_tr.left()->curve_end(), he_tr.left()->curve_end())  );
  
#else
  
  CGAL_assertion(traits != NULL);
  CGAL_precondition(traits->is_degenerate_point(v_tr));
  CGAL_precondition(v_tr.is_active());
  CGAL_precondition(traits->is_degenerate_curve(he_tr));
  CGAL_precondition(he_tr.is_active());
  CGAL_precondition(
    traits->equal_curve_end_2_object()
      (v_tr.left()->curve_end(), he_tr.right()->curve_end()) ||
    traits->equal_curve_end_2_object()
      (v_tr.left()->curve_end(), he_tr.left()->curve_end())  );
  
#endif
  
  /* update (in this order)
     v_tr.lb()
       if no Halfedges adjacent to the point emanating toward up
       or right exist - returns null, otherwise return
       the first Halfedge sweeped using a counter clockwise sweep
       starting from up direction not including.
     v_tr.rt()
       if no Halfedges adjacent to the point emanating toward bottom
       or left exist returns null, otherwise return
       the first Halfedge sweeped using a counter clockwise sweep
       starting from bottom direction not including.
     he_tr.rt()
       next clockwise degenerate_curve around rightmost v_tr (possibly
       himself)
     he_tr.lb()
       next clockwise degenerate_curve around leftmost v_tr (possibly
       himself)
  */
  Halfedge_const_handle he = he_tr.top();
  const Curve_end ce( v_tr.left()->curve_end());
  X_trapezoid* rt = v_tr.rt();
  X_trapezoid* lb = v_tr.lb();

  if(traits->equal_curve_end_2_object()(ce, he_tr.left()->curve_end()))
  {  //if the end point value equals the he_tr curve left value

    if (!rt && !lb)
    { // empty circulator
      v_tr.set_rt(&he_tr);
      he_tr.set_lb(&he_tr);
    }
    else
    {
      //set circ[0] to first Td_Edge on a counter clockwise 
      //  sweep starting at te 
      Around_point_circulator circ(traits, ce, rt ? rt : lb);
      Around_point_circulator stopper = circ;
      // if !rt set circ to lb
      // otherwise advance as required

      //if ce is interior
      if ((traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce()) 
                                                      == ARR_INTERIOR) &&
          (traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce()) 
                                                      == ARR_INTERIOR)  )
      {
        //if ce is interior then there might be other curves starting/ending
        // at ce. therefore we need the comparison arround ce.
        // if ce is on the boundaries there is only one curve at ce.

        //the underlying point of ce
        const Point& p = (ce.ce() == ARR_MIN_END) ?
                          traits->construct_min_vertex_2_object()(ce.cv()) :
                          traits->construct_max_vertex_2_object()(ce.cv()) ;

#ifdef CGAL_TD_DEBUG
      
        Around_point_circulator first_circ(circ);
      
#endif
        while (traits->compare_cw_around_point_2_object ()
               (circ->top()->curve(), 
                is_edge_to_right(circ->top(), p),
                he->curve(), 
                is_edge_to_right(he,p), p)        == SMALLER)
        {
          circ++;
          if (circ == stopper)
            break;
          
#ifdef CGAL_TD_DEBUG
        
          CGAL_assertion(first_circ != circ);
          CGAL_assertion(circ->is_active()); 
#endif
        }
      
#ifdef CGAL_TD_DEBUG
      
        CGAL_assertion(traits->compare_cw_around_point_2_object()
                       (circ->top()->curve(), 
                        is_edge_to_right(circ->top(), p),
                        te->curve(), 
                        is_edge_to_right(te,p), p)   != EQUAL);
#endif
      }
      
      circ.insert(he_tr);
      // set v_tr.lb()
      // set v_tr.rt();
      if (lb)
      {
        Around_point_circulator lb_circ(traits, ce, lb);
        if (!rt) v_tr.set_rt(lb);
        if (lb_circ.operator->() == &he_tr) v_tr.set_lb(&he_tr);
      }
      else
      {
         //if ce is interior
        if ((traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce()) 
                                                      == ARR_INTERIOR) &&
            (traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce()) 
                                                      == ARR_INTERIOR)  )
        {
          //if ce is interior then there might be other curves starting/ending
          // at ce. therefore we need the comparison arround ce.
          // if ce is on the boundaries there is only one curve at ce.
          
          //the underlying point of ce
          const Point& p = (ce.ce() == ARR_MIN_END) ?
                          traits->construct_min_vertex_2_object()(ce.cv()) :
                          traits->construct_max_vertex_2_object()(ce.cv()) ;
          
          if (traits->compare_cw_around_point_2_object()
                (rt->top()->curve(), is_edge_to_right(rt->top(), p),
                 he->curve(), is_edge_to_right(he, p), p, false)  == SMALLER)
          {
            v_tr.set_rt(&he_tr);
          }
        }

      }
    }
  }
  else //if the end point value equals the he_tr curve right value
  {
    if (!rt && !lb)
    { // empty circulator
      v_tr.set_lb(&he_tr);
      he_tr.set_rt(&he_tr);
    }
    else
    {
      /* set circ[0] to first Halfedge on a counter clockwise 
         sweep starting at te */
      Around_point_circulator circ(traits,ce,lb ? lb : rt);
      Around_point_circulator stopper = circ;
      // if !lb set circ to rt
      // otherwise advance as required

      //if ce is interior
      if ((traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce()) 
                                                    == ARR_INTERIOR) &&
          (traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce()) 
                                                    == ARR_INTERIOR)  )
      {
        //if ce is interior then there might be other curves starting/ending
        // at ce. therefore we need the comparison arround ce.
        // if ce is on the boundaries there is only one curve at ce.
        
        //the underlying point of ce
        const Point& p = (ce.ce() == ARR_MIN_END) ?
                          traits->construct_min_vertex_2_object()(ce.cv()) :
                          traits->construct_max_vertex_2_object()(ce.cv()) ;
        
        while (traits->compare_cw_around_point_2_object()
               (circ->top()->curve(), is_edge_to_right(circ->top(),p),
                he->curve(), is_edge_to_right(he,p), p, false)  == SMALLER)
        {
          circ++;
          if (circ == stopper)
            break;
        }
      
#ifdef CGAL_TD_DEBUG
      
        CGAL_assertion(traits->compare_cw_around_point_2_object()
                  (circ->top()->curve(), is_edge_to_right(circ->top(),p),
                   he->curve(), is_edge_to_right(he,p), p, false) != EQUAL);
#endif

      } 
     
      circ.insert(he_tr);
      // set v_tr.lb()
      // set v_tr.rt();
      if (rt)
      { 
        Around_point_circulator rt_circ(traits,ce,rt);
        if (!lb) v_tr.set_lb(rt);
        if (rt_circ.operator->() == &he_tr) v_tr.set_rt(&he_tr);
      }
      else
      {
        //if ce is interior
        if ((traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce()) 
                                                    == ARR_INTERIOR) &&
            (traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce()) 
                                                    == ARR_INTERIOR)  )
        {
          //if ce is interior then there might be other curves starting/ending
          // at ce. therefore we need the comparison arround ce.
          // if ce is on the boundaries there is only one curve at ce.

          //the underlying point of ce
          const Point& p = (ce.ce() == ARR_MIN_END) ?
                          traits->construct_min_vertex_2_object()(ce.cv()) :
                          traits->construct_max_vertex_2_object()(ce.cv()) ;

          if(traits->compare_cw_around_point_2_object()
                       (lb->top()->curve(), is_edge_to_right(lb->top(),p),
                        he->curve(), is_edge_to_right(he,p), p) == SMALLER)
          {
            v_tr.set_lb(&he_tr);
          }
        } 
      }
    }
  }
}




//-----------------------------------------------------------------------------
// Description:
//  Update top(),bottom() for trapezoid
//  Update rt,lb
// remarks:
//  The point degenerate trapezoid representing a point (edge_end) ee holds as its top and
//  bottom curves
//  the output for a vertical ray shoot queries immidiately below the point
//  toward up and
//  immediately above the point toward down respectively.
//optimization:
//  Each degenerate X_curve trapezoid emanating from the point p holds a pointer
//  to the next
//  trapezoid in a clockwise sweep around ee(possibly to itself).
//  This pointer is stored in rt or lb depending on the trapezoid is top right
//  or bottom left of ee.
//  For the trapezoid representing ee, rt and lb hold the previous X_curve
//  degenerate trapezoid
//  in a clockwise sweep to the first top right and bottom left respectively.
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::remove_halfedge_at_vertex_using_geometry(const X_trapezoid& he_tr,
                                            X_trapezoid& v_tr)
{
  
#ifndef CGAL_TD_DEBUG
  
  CGAL_warning(traits != NULL);
  CGAL_warning(traits->is_degenerate_point(v_tr));
  CGAL_warning(v_tr.is_active());
  CGAL_warning(traits->is_degenerate_curve(he_tr));
  CGAL_warning(he_tr.is_active());
  CGAL_warning(
    traits->equal_curve_end_2_object()
      //MICHAL: is this better:( v_tr.left()->curve_end(), he_tr.top()->max_vertex()->curve_end()) ||
      (v_tr.left()->curve_end(), he_tr.right()->curve_end()) ||
    traits->equal_curve_end_2_object()
      //MICHAL: is this better:( v_tr.left()->curve_end(), he_tr.top()->min_vertex()->curve_end()) ||
       (v_tr.left()->curve_end(), he_tr.left()->curve_end()));
  
#else
  
  CGAL_assertion(traits != NULL);
  CGAL_precondition(traits->is_degenerate_point(v_tr));
  CGAL_precondition(v_tr.is_active());
  CGAL_precondition(traits->is_degenerate_curve(he_tr));
  CGAL_precondition(he_tr.is_active());
  CGAL_precondition(
    traits->equal_curve_end_2_object()
      (v_tr.left()->curve_end(), he_tr.right()->curve_end()) ||
    traits->equal_curve_end_2_object()
      (v_tr.left()->curve_end(), he_tr.left()->curve_end()));
  
  
#endif
  
  /* update (in this order)
     v_tr.lb()
       if no curves adjacent to the point eminating toward up
       or right exist returns null, otherwise return
       the first X_curve sweeped using a counter clockwise sweep
       starting from up direction not including.
     v_tr.rt()
       if no curves adjacent to the point eminating toward bottom
       or left exist returns null, otherwise return
       the first X_curve sweeped using a counter clockwise sweep
       starting from bottom direction not including.
     he_tr.rt()
       next clockwise degenerate_curve around rightmost v_tr (possibly
       himself)
     he_tr.lb()
       next clockwise degenerate_curve around leftmost v_tr (possibly
       himself)
  */
  Halfedge_const_handle he = he_tr.top();
  const Curve_end ce(v_tr.left()->curve_end()); 
  Around_point_circulator 
    prev_top(traits,ce,v_tr.rt()),
    prev_bottom(traits,ce,v_tr.lb());
  
  // update bottom
  if ((he == v_tr.bottom()) || (he->twin() == v_tr.bottom())) //MICHAL: he comp
  {
    Around_point_circulator bottom = (!!prev_bottom)? prev_bottom : prev_top;
    bottom++;
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_assertion(!!bottom);
    
#endif
    
    if (!bottom->is_on_bottom_boundary())
      v_tr.set_bottom(bottom->bottom());
    else
      v_tr.set_is_on_bottom_boundary(true);
  }

  // update top
  if ((he == v_tr.top()) || (he->twin() == v_tr.top())) //MICHAL: he comp
  {
    Around_point_circulator top = (!!prev_top)? prev_top : prev_bottom;
    top++;
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_assertion(!!top);
    
#endif
    
    if (!top->is_on_top_boundary())
      v_tr.set_top(top->top());
    else
      v_tr.set_is_on_top_boundary(true);
  }
  
  //update right top neighbour and left bottom neighbour
  bool b = (traits->compare_curve_end_xy_2_object()
                          (ce, he_tr.right()->curve_end()) == SMALLER);

  Around_point_circulator circ(traits, ce, b ? v_tr.rt(): v_tr.lb());
  
#ifdef CGAL_TD_DEBUG
  
  CGAL_precondition(!!circ);
  
#endif
  
  while(*circ != he_tr) circ++;
  X_trapezoid* removed = circ.operator->();
  circ.remove();
  if(!!circ)
  {
    X_trapezoid* effective_curr = circ[0];
    if (v_tr.rt() == removed)
      v_tr.set_rt(effective_curr);
    if (v_tr.lb() == removed)
      v_tr.set_lb(effective_curr);

    Around_point_circulator rt_circ(traits, ce,
                                    v_tr.rt());
    if (!!rt_circ)
    {
      rt_circ++;
      if (rt_circ.is_right_rotation())
        v_tr.set_rt(0);
    }

    Around_point_circulator lb_circ(traits, ce,v_tr.lb());
    if (!!lb_circ)
    {
      lb_circ++;
      if (!lb_circ.is_right_rotation())
        v_tr.set_lb(0);
    }
  }
  else
  {
    v_tr.set_rt(0);
    v_tr.set_lb(0);
  }
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
typename Trapezoidal_decomposition_2<Td_traits>::X_trapezoid &
Trapezoidal_decomposition_2<Td_traits>
::set_trp_params_after_halfedge_insertion (Halfedge_const_handle he,
                                            const Curve_end& ce,
                                            X_trapezoid*& v_tr,
                                            const Locate_type&
                                              CGAL_precondition_code(lt))
{
  CGAL_assertion(traits != NULL);
  CGAL_precondition(lt == POINT);
  //ee is interior
  CGAL_assertion((traits->parameter_space_in_x_2_object()(ce.cv(), ce.ce()) 
                                                          == ARR_INTERIOR) &&
                 (traits->parameter_space_in_y_2_object()(ce.cv(), ce.ce()) 
                                                          == ARR_INTERIOR) );
  //the underlying point of ce
  const Point& p = (ce.ce() == ARR_MIN_END) ?
                    traits->construct_min_vertex_2_object()(ce.cv()) :
                    traits->construct_max_vertex_2_object()(ce.cv()) ;
  
  //set top to hold the halfedge whose source is p, 
  // which is clockwise "smallest" starting from top (12 o'clock)
  if (traits->compare_cw_around_point_2_object()
      (he->curve(), is_edge_to_right(he,p),
       v_tr->top()->curve(), 
       is_edge_to_right(v_tr->top(),p), p) == SMALLER)
  {
    v_tr->set_top(he);
  }

  //set bottom to hold the halfedge whose source is p, 
  // which is clockwise "smallest" starting from bottom (6 o'clock)
  if (traits->compare_cw_around_point_2_object()
      (he->curve(), is_edge_to_right(he,p), 
       v_tr->bottom()->curve(),
       is_edge_to_right(v_tr->bottom(),p), p, false) == SMALLER)
  {
    v_tr->set_bottom(he);
  }

  return *v_tr;
}
  
//-----------------------------------------------------------------------------
// Description:
template <class Td_traits> 
typename Trapezoidal_decomposition_2<Td_traits>::X_trapezoid &
Trapezoidal_decomposition_2<Td_traits>
::insert_curve_at_point_using_dag (Halfedge_const_handle he,
                                   Vertex_const_handle v,
                                   X_trapezoid*& tr,
                                   const Locate_type&
                                   CGAL_precondition_code(lt))
{
  CGAL_precondition(lt==TRAPEZOID || lt==UNBOUNDED_TRAPEZOID);
  
  Dag_node *tt = tr->dag_node();
  
#ifdef CGAL_TD_DEBUG
  
  CGAL_assertion(tr->dag_node());
  
#endif
  

  //MICHAL: maybe insert into split_trapezoid_by_vertex
  CGAL_assertion_code(Halfedge_const_handle invalid_he);
  CGAL_assertion (he != invalid_he);
  
  Curve_end_pair ce_pair = v->curve_end();

  //MICHAL: 

  if (!traits->equal_2_object()
    (he->curve(), *(ce_pair.first))) {
    X_monotone_curve_2 cv1111 = he->curve();
    X_monotone_curve_2 cv2222 = *(ce_pair.first);
  }
  //CGAL_assertion(traits->equal_2_object()
  //                        (he->curve(), *(ce_pair.first)));
  
  //we need to send the halfedge whose source is v.
  if (he->source() == v) 
 // if ((ce_pair.second == ARR_MIN_END && he->direction() == ARR_LEFT_TO_RIGHT) ||
 //     (ce_pair.second == ARR_MAX_END && he->direction() == ARR_RIGHT_TO_LEFT)  )
  {
    return *split_trapezoid_by_vertex(*tt, v, he, he);
  }
  else
  {
    return *split_trapezoid_by_vertex(*tt, v, he->twin(), he->twin());
  }
}
        
//-----------------------------------------------------------------------------
// Description:
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::set_trp_params_after_halfedge_update(Halfedge_const_handle old_he,
                                       Halfedge_const_handle new_he,
                                       X_trapezoid& v_tr)
{
  CGAL_precondition_code(Halfedge_const_handle invalid_he);
  CGAL_precondition(old_he != invalid_he);
  CGAL_precondition(traits->is_degenerate_point(v_tr));
  CGAL_assertion(v_tr.is_active());

  //make sure the top & bottom are added in same direction as before
  //  such that v_tr is the source (done inside the set methods)
  if (v_tr.top() == old_he || v_tr.top()->twin() == old_he) //MICHAL: he comp
  {
    v_tr.set_top(new_he);
  }
  if (v_tr.bottom() == old_he || v_tr.bottom()->twin() == old_he) //MICHAL: he comp
  {
    v_tr.set_bottom(new_he);
  }
}

//-----------------------------------------------------------------------------
// Description:
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::set_trp_params_after_halfedge_update(const X_monotone_curve_2& old_cv,
                                       Halfedge_const_handle new_he,
                                       X_trapezoid& v_tr)
{
  

  CGAL_precondition(traits->is_degenerate_point(v_tr));
  CGAL_assertion(v_tr.is_active());

  //MICHAL: used to be:
  // if (!v_tr.is_on_top_boundary() && traits->equal_2_object()(v_tr.top(), old_cv))
  //but curve cannot be on top boundary
 
  //make sure the top & bottom are added in same direction as before
  //  such that v_tr is the source (done inside the set methods)
  if (traits->equal_2_object()(v_tr.top()->curve(), old_cv))
  {
    v_tr.set_top(new_he);
  }
  if (traits->equal_2_object()(v_tr.bottom()->curve(), old_cv)) 
  {
    v_tr.set_bottom(new_he);
  }
}

//-----------------------------------------------------------------------------
// Description:
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::set_trp_params_after_split_halfedge_update(Halfedge_const_handle new_he,
                                            X_trapezoid& v_tr,
                                            Halfedge_const_handle he1, 
                                            Halfedge_const_handle he2)
{
  
#ifdef CGAL_TD_DEBUG
  
  CGAL_precondition(traits->is_degenerate_point(sep));
  
#endif
  
  //make sure the top & bottom are added in same direction as before
  //  such that sep is the source (done inside the set methods)
  if ((v_tr.top() == he1) || (v_tr.top() == he1->twin()) || //MICHAL: he comp
      (v_tr.top() == he2) || (v_tr.top() == he2->twin()) )  //MICHAL: he comp
  {
    v_tr.set_top(new_he);
  }
  if ((v_tr.bottom() == he1) || (v_tr.bottom() == he1->twin()) || //MICHAL: he comp
      (v_tr.bottom() == he2) || (v_tr.bottom() == he2->twin()) )  //MICHAL: he comp
  {
    v_tr.set_bottom(new_he);
  }
}
 
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
::set_trp_params_after_halfedge_update (In_face_iterator& it,
                                        Halfedge_const_handle old_he,
                                        Halfedge_const_handle new_he,
                                        Vertex_const_handle min_v,
                                        Vertex_const_handle max_v,
                                        X_trapezoid*& end)
{
  X_trapezoid* last = 0;
 
  while (it.operator->() != end)
  {
    CGAL_assertion(it->is_active());

    if (!it->is_on_top_boundary() &&
        (it->top() == old_he || it->top() == old_he->twin())) //MICHAL: he comp
    {
      it->set_top(new_he);
    }
    if (!it->is_on_bottom_boundary() &&
        (it->bottom() == old_he || it->bottom() == old_he->twin())) //MICHAL: he comp
    {
      it->set_bottom(new_he);
    }
    //MICHAL: new code
    if (traits->is_degenerate_curve(*(it.operator->())))
    {
      it->set_left(min_v);
      it->set_right(max_v);
    }

    last = it.operator->();
    ++it;
  }
  end = last;
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
::search_using_dag (Dag_node& curr,
                    const Traits* traits,
                    const Point& p,
                    Halfedge_const_handle* p_he,
                    Comparison_result up /*=EQUAL*/) const
{
  //Halfedge_const_handle*  p_top; 
  Halfedge_const_handle top_he; 

#ifdef CGAL_TD_DEBUG
  X_trapezoid* old = NULL;
#endif
  
  while(true)
  {
#ifdef CGAL_TD_DEBUG
    // unbounded loop
    CGAL_assertion(curr.operator->() != old);
    old = curr.operator->();
#endif
    //curr is the current pointer to node in the data structure
    if (traits->is_degenerate_point(*curr))
    { // if the trapezoid (curr) represents a point
      
      bool pnt_exists = false;
      if (!curr->is_active() && !curr->is_on_boundaries())
        pnt_exists = true;
      
      // extract point (curve_end) from trapezoid
      const Curve_end left_ce(curr->is_active()? 
        curr->left()->curve_end() : curr->curve_end_for_rem_vtx());
      Point left_p;
      if (pnt_exists)
        left_p = curr->point_for_inner_rem_vtx();

      if ((!pnt_exists && is_end_point_left_low(p, left_ce)) ||
          (pnt_exists &&  is_end_point_left_low(p, left_p)) )
      {
        curr = curr.left_child();
        continue;
      }
      else if ((!pnt_exists && is_end_point_right_top(p, left_ce)) ||
               (pnt_exists &&  is_end_point_right_top(p, left_p)) )
      {
        curr = curr.right_child();
        continue;
      }
      else if ((!pnt_exists && traits->equal_curve_end_2_object()(left_ce, p)) ||
               (pnt_exists &&  traits->equal_2_object()(left_p, p)) )
      {
        if (!p_he) //if p_he was not given
        {
          if ( up == EQUAL ) 
          {      // point found!
            if (curr->is_active()) 
              return POINT;
            curr = curr.left_child();
          }
          else if ( up == LARGER ) {          // vertical ray shut up
            curr = curr.right_child();                      
          }
          else /*if ( up == SMALLER ) */ {
            curr = curr.left_child();               // vertical ray shut down
          }
          continue;
        }
        else //if p_he was given
        {
          bool is_equal_to_he_min = traits->equal_curve_end_2_object()
                               (Curve_end(*p_he,ARR_MIN_END), p);
          bool is_equal_to_he_max = traits->equal_curve_end_2_object()
                               (Curve_end(*p_he,ARR_MAX_END), p);
          
#ifndef CGAL_TD_DEBUG
          CGAL_warning( is_equal_to_he_min || is_equal_to_he_max );                 
#else
          CGAL_assertion( is_equal_to_he_min || is_equal_to_he_max );
#endif

          curr = is_equal_to_he_min ? curr.right_child() : curr.left_child();
          // (Oren 14/4/02) ??
                  
          continue;
        }
      }
      else
      {
              
        CGAL_assertion((!pnt_exists && 
                      (is_end_point_left_low(p,left_ce) ||
                       is_end_point_right_top(p,left_ce) ||
                       traits->equal_curve_end_2_object()(left_ce,p))) ||
                     (pnt_exists && 
                      (is_end_point_left_low(p,left_p) ||
                       is_end_point_right_top(p,left_p) ||
                       traits->equal_2_object()(left_p,p))));

        return Locate_type();
      }
      

/*

      // extract point (curve_end) from trapezoid
      p_left = &curr->left();
      const Curve_end left_ce((*p_left)->curve_end());

      if (is_end_point_left_low(p, left_ce))
      {
        curr = curr.left_child();
        continue;
      }
      else if (is_end_point_right_top(p, left_ce))
      {
        curr = curr.right_child();
        continue;
      }
      else if (traits->equal_curve_end_2_object()(left_ce, p))
      {
        if (!p_he) //if p_he was not given
        {
          if ( up == EQUAL ) 
          {      // point found!
            if (curr->is_active()) 
              return POINT;
            curr = curr.left_child();
          }
          else if ( up == LARGER ) {          // vertical ray shut up
            curr = curr.right_child();                      
          }
          else //if ( up == SMALLER )  
          {
            curr = curr.left_child();               // vertical ray shut down
          }
          continue;
        }
        else //if p_he was given
        {
          bool is_equal_to_he_min = traits->equal_curve_end_2_object()
                               (Curve_end(*p_he,ARR_MIN_END), p);
          bool is_equal_to_he_max = traits->equal_curve_end_2_object()
                               (Curve_end(*p_he,ARR_MAX_END), p);
          
#ifndef CGAL_TD_DEBUG
          CGAL_warning( is_equal_to_he_min || is_equal_to_he_max );                 
#else
          CGAL_assertion( is_equal_to_he_min || is_equal_to_he_max );
#endif

          curr = is_equal_to_he_min ? curr.right_child() : curr.left_child();
          // (Oren 14/4/02) ??
                  
          continue;
        }
      }
      else
      {
              
#ifndef CGAL_TD_DEBUG
        CGAL_warning(is_end_point_left_low(p,left_ce) ||
                     is_end_point_right_top(p,left_ce) ||
                     traits->equal_curve_end_2_object()(left_ce,p));
#else
        CGAL_assertion(is_end_point_left_low(p,left_ce) ||
                       is_end_point_right_top(p,left_ce) ||
                       traits->equal_curve_end_2_object()(left_ce,p));
#endif

        return Locate_type();
      }
      */
    }
    if (traits->is_degenerate_curve(*curr))
    { // if the trapezoid (curr) represents a curve, 
      //   so top() is a real Halfedge with a curve() if curr is active
      //   or curr holds the curve if curr is not active 
      const X_monotone_curve_2* p_he_cv = 
        (curr->is_active()) ? &curr->top()->curve() : &curr->curve_for_rem_he();

      Comparison_result cres = traits->compare_y_at_x_2_object()
                                                 (p, *p_he_cv);
      if (cres == SMALLER)
      {
        curr = curr.left_child();
        continue;
      }
      else if (cres == LARGER)
      {
        curr = curr.right_child();
        continue;
      }
      else
      {  
        // p is on the CURVE (top_he = curr.top()) itself 
#ifndef CGAL_TD_DEBUG      
         CGAL_warning(
          (cres == EQUAL) &&
          (traits->compare_curve_end_x_2_object()
                        (p, Curve_end(*p_he_cv,ARR_MAX_END)) != LARGER) &&
          (traits->compare_curve_end_x_2_object()
                        (p, Curve_end(*p_he_cv,ARR_MIN_END)) != SMALLER));
#else
              
       CGAL_assertion(
          (cres == EQUAL) &&
          (traits->compare_curve_end_x_2_object()
                        (p, Curve_end(*p_he_cv,ARR_MAX_END)) != LARGER) &&
          (traits->compare_curve_end_x_2_object()
                        (p, Curve_end(*p_he_cv,ARR_MIN_END)) != SMALLER));
#endif
        if (!p_he) //if p_he was not given
        {
          // For a vertical curve, we always visit it after visiting
          // one of its endpoints.
          if ((up == EQUAL) || traits->is_vertical(*curr)) {
            //std::cout << "EQUAL or VERTICAL" << std::endl;
            if (curr->is_active()) return CURVE; 
            curr = curr.left_child();
          }
          else if (up == LARGER) {
            curr = curr.right_child();
          }
          else { // if (up==SMALLER)
            curr = curr.left_child();
          }
          continue;
        }
        else //if p_he was given
        {
          //p is a parameter space interior point

          bool is_min_equal = traits->equal_curve_end_2_object()
                                      (Curve_end(*p_he,ARR_MIN_END), 
                                       Curve_end(*p_he_cv,ARR_MIN_END));

          bool is_max_equal = traits->equal_curve_end_2_object()
                                      (Curve_end(*p_he,ARR_MAX_END), 
                                       Curve_end(*p_he_cv,ARR_MAX_END));

#ifndef CGAL_TD_DEBUG          
          CGAL_warning (is_min_equal || is_max_equal);
#else
           if (!is_min_equal && !is_max_equal)
          {
            std::cerr << "\n*p_he_cv " << *p_he_cv;
            std::cerr << "\n*p_he "  << *p_he << std::endl;
            CGAL_assertion(is_min_equal || is_max_equal);
          }
#endif
          
           Comparison_result res = 
             is_min_equal ? 
              traits->compare_cw_around_point_2_object()
                       (*p_he_cv, is_curve_to_right(*p_he_cv,p),
                        (*p_he)->curve(), is_edge_to_right(*p_he,p), p) :
              traits->compare_cw_around_point_2_object()
                       ((*p_he)->curve(), is_edge_to_right(*p_he,p),
                        *p_he_cv, is_curve_to_right(*p_he_cv,p), p ,false);

          switch(res)
          {
           case LARGER:
            curr = curr.right_child();
            break;
           case SMALLER:
            curr = curr.left_child();
            break;
           case EQUAL:
            switch(up)
            {
             case LARGER:
              curr = curr.right_child();
              break;
             case SMALLER:
              curr = curr.left_child();
              break;
             case EQUAL:
              if (curr->is_active()) return CURVE;
              curr = curr.left_child();
              break;
                          
#ifdef CGAL_TD_DEBUG
             default:
              CGAL_assertion(up == LARGER || up == SMALLER || up == EQUAL);
              return Locate_type();
#endif
            }
            break;

#ifdef CGAL_TD_DEBUG
           default:
            CGAL_assertion(res == LARGER || res == SMALLER || res == EQUAL);
            return Locate_type();
#endif

          }
        }
      } 
    }
    else
    {
      // if is_degenerate() == 0, meaning: the trapezoid (curr)
      // is neither a point nor a curve , but a real trapezoid
      if (curr->is_active())
        return curr->is_on_boundaries() ? UNBOUNDED_TRAPEZOID : TRAPEZOID;
      curr = curr.left_child();
      continue;
    }
  }
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
void
Trapezoidal_decomposition_2<Td_traits>
::search_and_print_using_dag (std::ostream& out, Dag_node& curr,
                    const Traits* traits,
                    const Point& p,
                    Halfedge_const_handle* p_he,
                    Comparison_result up /*=EQUAL*/) const
{
  out << "QUERY: " << std::endl;
  out << "x: " << CGAL::to_double(p.x())
         << ", y: " << CGAL::to_double(p.y()) << std::endl;

  
  Halfedge_const_handle top_he; 

#ifdef CGAL_TD_DEBUG
  X_trapezoid* old = NULL;
#endif
  
  while(true)
  {
#ifdef CGAL_TD_DEBUG
    // unbounded loop
    CGAL_assertion(curr.operator->() != old);
    old = curr.operator->();
#endif
    //curr is the current pointer to node in the data structure
    if (traits->is_degenerate_point(*curr))
    { // if the trapezoid (curr) represents a point
      
      bool pnt_exists = false;
      if (!curr->is_active() && !curr->is_on_boundaries())
        pnt_exists = true;
      
      // extract point (curve_end) from trapezoid
      const Curve_end left_ce(curr->is_active()? 
        curr->left()->curve_end() : curr->curve_end_for_rem_vtx());
      Point left_p;
      if (pnt_exists)
        left_p = curr->point_for_inner_rem_vtx();
      
      out << " POINT : " ;
      if (curr->is_active())
        out << " (active) ";
      else
        out << " (inactive) ";
      print_ce_data(left_ce.cv(), left_ce.ce(), out);
      
      if ((!pnt_exists && is_end_point_left_low(p, left_ce)) ||
          (pnt_exists &&  is_end_point_left_low(p, left_p)) )
      {
        out << " Going left " << std::endl;
        curr = curr.left_child();
        continue;
      }
      else if ((!pnt_exists && is_end_point_right_top(p, left_ce)) ||
               (pnt_exists &&  is_end_point_right_top(p, left_p)) )
      {
        out << " Going right " << std::endl;
        curr = curr.right_child();
        continue;
      }
      else if ((!pnt_exists && traits->equal_curve_end_2_object()(left_ce, p)) ||
               (pnt_exists &&  traits->equal_2_object()(left_p, p)) )
      {
        out << " Equal to query " << std::endl;
        if (!p_he) //if p_he was not given
        {
          
          if ( up == EQUAL ) 
          {      // point found!
            if (curr->is_active()) 
            {
              out << " Found active point! " << std::endl;
              return;
            }
            out << " (equal to inactive point) Going left " << std::endl;
            curr = curr.left_child();
          }
          else if ( up == LARGER ) {          // vertical ray shut up
            out << " Going right " << std::endl;
            curr = curr.right_child();                      
          }
          else /*if ( up == SMALLER ) */ {
            out << " Going left " << std::endl;
            curr = curr.left_child();               // vertical ray shut down
          }
          continue;
        }
        else //if p_he was given
        {
          //NOT GONNA HAPPEN         
        }
      }
      else
      {
        out << " Problem - comparing to point" << std::endl;     
        return;
      }
      

    }
    if (traits->is_degenerate_curve(*curr))
    { // if the trapezoid (curr) represents a curve, 
      //   so top() is a real Halfedge with a curve() if curr is active
      //   or curr holds the curve if curr is not active 
      const X_monotone_curve_2* p_he_cv = 
        (curr->is_active()) ? &curr->top()->curve() : &curr->curve_for_rem_he();

      // if the trapezoid (curr) represents a curve, 
      //   so top() is a real Halfedge with a curve() if curr is active
      //   or curr holds the curve if curr is not active 
      out << " CURVE : " ;
      if (curr->is_active())
        out << " (active) ";
      else
        out << " (inactive) ";
      if (traits->is_vertical(*curr))
        out << " (vertical) ";

      print_cv_data(*p_he_cv, out);
      
      Comparison_result cres = traits->compare_y_at_x_2_object()
                                                 (p, *p_he_cv);
      if (cres == SMALLER)
      {
        out << " Going left " << std::endl;
        curr = curr.left_child();
        continue;
      }
      else if (cres == LARGER)
      {
        out << " Going right " << std::endl;
        curr = curr.right_child();
        continue;
      }
      else
      {  
        // p is on the CURVE (top_he = curr.top()) itself 
        out << " query is on the curve " << std::endl;
              
       CGAL_assertion(
          (cres == EQUAL) &&
          (traits->compare_curve_end_x_2_object()
                        (p, Curve_end(*p_he_cv,ARR_MAX_END)) != LARGER) &&
          (traits->compare_curve_end_x_2_object()
                        (p, Curve_end(*p_he_cv,ARR_MIN_END)) != SMALLER));
        if (!p_he) //if p_he was not given
        {
          // For a vertical curve, we always visit it after visiting
          // one of its endpoints.
          if ((up == EQUAL) || traits->is_vertical(*curr)) 
          {
            //std::cout << "EQUAL or VERTICAL" << std::endl;
            if (curr->is_active()) 
            {
              out << " On active curve or on vertical one  " << std::endl;            
              return; 
            }
            out << " (equal to inactive curve) Going left " << std::endl;
            curr = curr.left_child();
          }
          else if (up == LARGER) {
            out << " Going right " << std::endl;
            curr = curr.right_child();
          }
          else { // if (up==SMALLER)
            out << " Going left " << std::endl;
            curr = curr.left_child();
          }
          continue;
        }
        else //if p_he was given
        {
          //p is a parameter space interior point
          //NOT GONNA HAPPEN
        }
      } 
    }
    else
    {
      // if is_degenerate() == 0, meaning: the trapezoid (curr)
      // is neither a point nor a curve , but a real trapezoid
      out << " TRAPEZOID : " ;
      if (curr->is_active())
      {
        out << " (active) ";
        if (curr->is_on_boundaries())
          out << " UNBOUNDED! ";
        else
          out << " BOUNDED! ";
      }
      else
        out << " (inactive) ";
      if (curr->is_active())
        return;
      out << " (on inactive trapezoid) Going left " << std::endl;
      curr = curr.left_child();
      continue;
    }
  }

}
  
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
::search_using_dag (Dag_node& curr,
                    const Traits* traits,
                    const Curve_end& ce,
                    Halfedge_const_handle* p_he,
                    Comparison_result up /*=EQUAL*/) const
{
  if (p_he == NULL)
    return search_using_dag_with_cv (curr,traits,ce,NULL, up);
  else
  {
    CGAL_precondition_code(Halfedge_const_handle invalid_he);
    CGAL_precondition(*p_he != invalid_he);
    return search_using_dag_with_cv (curr,traits,ce,&(*p_he)->curve(), up);
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
::search_using_dag_with_cv (Dag_node& curr,
                            const Traits* traits,
                            const Curve_end& ce,
                            const X_monotone_curve_2* p_cv,
                            Comparison_result up /*=EQUAL*/) const
{

#ifdef MICHAL_DEBUG
  std::cout << "SEARCHING: ----------" << std::endl ;
  std::cout << "given ce: " << std::endl;
  print_ce_data(ce.cv(), ce.ce());
  
  if (p_cv)
  {
    std::cout << "p_cv is given : " << std::endl;
    print_cv_data(*p_cv);
  }
  std::cout << "search path:" << std::endl;

#endif

  //Halfedge_const_handle*  p_top; 
  Halfedge_const_handle  top_he; 
 // const X_monotone_curve_2* p_trp_cv;
  
#ifdef CGAL_TD_DEBUG
  X_trapezoid* old = NULL;
#endif
  
  while(true)
  {
#ifdef CGAL_TD_DEBUG
    // unbounded loop
    CGAL_assertion(curr.operator->() != old);
    old = curr.operator->();
#endif
    //curr is the current pointer to node in the data structure
    if (traits->is_degenerate_point(*curr))
    { // if the trapezoid (curr) represents a point
      
      bool pnt_exists = false;
      if (!curr->is_active() && !curr->is_on_boundaries())
        pnt_exists = true;
      
      // extract point (curve_end) from trapezoid
      const Curve_end left_ce(curr->is_active()? 
        curr->left()->curve_end() : curr->curve_end_for_rem_vtx());
      Point left_p;
      if (pnt_exists)
        left_p = curr->point_for_inner_rem_vtx();

#ifdef MICHAL_DEBUG
      std::cout << "point:" << std::endl;
      print_ce_data( left_ce.cv(), left_ce.ce() );
#endif

      
      if ((!pnt_exists && is_end_point_right_top(left_ce, ce)) ||
          (pnt_exists &&  is_end_point_right_top(left_p, ce)) )
      {
#ifdef MICHAL_DEBUG
        std::cout << "--go left" << std::endl;
#endif
        curr = curr.left_child();
        continue;
      }
      else if ((!pnt_exists && is_end_point_left_low(left_ce, ce)) ||
               (pnt_exists &&  is_end_point_left_low(left_p, ce)) )
      {
#ifdef MICHAL_DEBUG
        std::cout << "--go right" << std::endl;
#endif
        curr = curr.right_child();
        continue;
      }
      else if ((!pnt_exists && traits->equal_curve_end_2_object()(left_ce, ce)) ||
               (pnt_exists &&  traits->equal_curve_end_2_object()(left_p, ce)) )
      {
        if (!p_cv) //if p_cv was not given
        {
          if ( up == EQUAL ) 
          {      // point found!
            if (curr->is_active()) 
              return POINT;
            curr = curr.left_child();
#ifdef MICHAL_DEBUG
            std::cout << "--go left" << std::endl;
#endif          
     
          }
          else if ( up == LARGER ) {          // vertical ray shut up
            curr = curr.right_child();                      
#ifdef MICHAL_DEBUG
            std::cout << "--go right" << std::endl;
#endif 
          }
          else { // if ( up == SMALLER ) 
            curr = curr.left_child();               // vertical ray shut down
#ifdef MICHAL_DEBUG
            std::cout << "--go left" << std::endl;
#endif 

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

          curr = is_equal_to_he_min ? curr.right_child() : curr.left_child();
          // (Oren 14/4/02) ??
#ifdef MICHAL_DEBUG
        if (is_equal_to_he_min)
          std::cout << "--go right" << std::endl;
        else
          std::cout << "--go left" << std::endl;
#endif          
          
          continue;
        }
      }
      else
      {
              
        CGAL_assertion((!pnt_exists && 
                      (is_end_point_left_low(ce,left_ce) ||
                       is_end_point_right_top(ce,left_ce) ||
                       traits->equal_curve_end_2_object()(left_ce,ce))) ||
                     (pnt_exists && 
                      (is_end_point_right_top(left_p,left_ce) ||
                       is_end_point_left_low(left_p,left_ce) ||
                       traits->equal_curve_end_2_object()(left_p,ce))));

        return Locate_type();
      }
      
/*
      // extract point (edge_end) from trapezoid

      p_left = &curr->left();
      const Curve_end left_ce((*p_left)->curve_end());

      if (is_end_point_left_low(ce, left_ce))
      {
        curr = curr.left_child();
        continue;
      }
      else if (is_end_point_left_low(left_ce, ce))
      {
        curr = curr.right_child();
        continue;
      }
      else if (traits->equal_curve_end_2_object()
                        (left_ce, ce))
      {
        if (!p_cv) //if p_cv was not given
        {
          if ( up == EQUAL ) 
          {      // point found!
            if (curr->is_active()) 
              return POINT;
            curr = curr.left_child();
          }
          else if ( up == LARGER ) {          // vertical ray shut up
            curr = curr.right_child();                      
          }
          else //if ( up == SMALLER ) 
          {
            curr = curr.left_child();               // vertical ray shut down
          }
          continue;
        }
        else //if p_cv was given
        {
          bool is_equal_to_cv_min = traits->equal_curve_end_2_object()
                                          (Curve_end(*p_cv,ARR_MIN_END), ce);
          bool is_equal_to_cv_max = traits->equal_curve_end_2_object()
                                          (Curve_end(*p_cv,ARR_MAX_END), ce);
          
#ifndef CGAL_TD_DEBUG
          CGAL_warning( is_equal_to_cv_min || is_equal_to_cv_max );                 
#else
          CGAL_assertion( is_equal_to_cv_min || is_equal_to_cv_max );
#endif

          curr = is_equal_to_cv_min ? curr.right_child() : curr.left_child();
          // (Oren 14/4/02) ??
                  
          continue;
        }
      }
      else
      {
              
#ifndef CGAL_TD_DEBUG
        CGAL_warning(is_end_point_left_low(ce,left_ce) ||
                     is_end_point_left_low(left_ce,ce) ||
                     traits->equal_curve_end_2_object()(left_ce,ce));
#else
        CGAL_assertion(is_end_point_left_low(ce,left_ce) ||
                       is_end_point_left_low(left_ce,ce) ||
                       traits->equal_curve_end_2_object()(left_ce,ce));
#endif

        return Locate_type();
      }
      */
    }
    if (traits->is_degenerate_curve(*curr))
    { // if the trapezoid (curr) represents a curve
      const X_monotone_curve_2* p_he_cv = 
        (curr->is_active()) ? &curr->top()->curve() : &curr->curve_for_rem_he();
      
#ifdef MICHAL_DEBUG
      std::cout << "curve:" << std::endl;
      print_cv_data(*p_he_cv);
#endif

      Comparison_result cres = traits->compare_curve_end_y_at_x_2_object()
                                                   (ce, *p_he_cv);
      if (cres == SMALLER)
      {
        curr = curr.left_child();
#ifdef MICHAL_DEBUG
          std::cout << "--go left" << std::endl;
#endif          
     
        continue;
      }
      else if (cres == LARGER)
      {
        curr = curr.right_child();
#ifdef MICHAL_DEBUG
          std::cout << "--go right" << std::endl;
#endif          

        continue;
      }
      else
      {  
        // ce is on the CURVE (*p_he_cv = ce.cv()) itself 
#ifndef CGAL_TD_DEBUG      
         CGAL_warning(
          (cres == EQUAL) &&
          (traits->compare_curve_end_x_2_object()
                        (ce, Curve_end(*p_he_cv,ARR_MAX_END)) != LARGER) &&
          (traits->compare_curve_end_x_2_object()
                        (ce, Curve_end(*p_he_cv,ARR_MIN_END)) != SMALLER));
#else
              
       CGAL_assertion(
          (cres == EQUAL) &&
          (traits->compare_curve_end_x_2_object()
                        (ce, Curve_end(*p_he_cv,ARR_MAX_END)) != LARGER) &&
          (traits->compare_curve_end_x_2_object()
                        (ce, Curve_end(*p_he_cv,ARR_MIN_END)) != SMALLER));

#endif

        if (!p_cv) //if p_cv was not given
        {
          // For a vertical curve, we always visit it after visiting
          // one of its endpoints.
          if ((up == EQUAL) || traits->is_vertical(*curr)) {
            //std::cout << "EQUAL or VERTICAL" << std::endl;
            if (curr->is_active()) return CURVE;
            curr = curr.left_child(); 
#ifdef MICHAL_DEBUG
            std::cout << "--go left" << std::endl;
#endif          

          }
          else if (up == LARGER) {
            curr = curr.right_child();
#ifdef MICHAL_DEBUG
            std::cout << "--go right" << std::endl;
#endif 
          }
          else { // if (up==SMALLER) 
            curr = curr.left_child();
#ifdef MICHAL_DEBUG
            std::cout << "--go left" << std::endl;
#endif 
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
                                       Curve_end(*p_he_cv,ARR_MIN_END));

            bool is_max_equal = traits->equal_curve_end_2_object()
                                      (Curve_end(*p_cv,ARR_MAX_END), 
                                       Curve_end(*p_he_cv,ARR_MAX_END));

#ifndef CGAL_TD_DEBUG          
            CGAL_warning (is_min_equal || is_max_equal);
#else
            if (!is_min_equal && !is_max_equal)
            {
              std::cerr << "\n*p_he_cv " << *p_he_cv;
              std::cerr << "\n*p_cv "  << *p_cv << std::endl;
              CGAL_assertion(is_min_equal || is_max_equal);
            }
#endif

            //the underlying point of ce
            const Point& p = (ce.ce() == ARR_MIN_END) ?
                           traits->construct_min_vertex_2_object()(ce.cv()) :
                           traits->construct_max_vertex_2_object()(ce.cv()) ;
          
            res = is_min_equal ? 
                    traits->compare_cw_around_point_2_object()
                           (*p_he_cv, is_curve_to_right(*p_he_cv,p),
                            *p_cv, is_curve_to_right(*p_cv,p), p)  :
                    traits->compare_cw_around_point_2_object()
                           (*p_cv, is_curve_to_right(*p_cv,p),
                            *p_he_cv, is_curve_to_right(*p_he_cv,p), p ,false);
          }

          switch(res)
          {
           case LARGER:
            curr = curr.right_child();
#ifdef MICHAL_DEBUG
            std::cout << "--go right" << std::endl;
#endif 
            break;
           case SMALLER:
            curr = curr.left_child();
#ifdef MICHAL_DEBUG
            std::cout << "--go left" << std::endl;
#endif 
            break;
           case EQUAL:
            switch(up)
            {
             case LARGER:
              curr = curr.right_child();
#ifdef MICHAL_DEBUG
              std::cout << "--go right" << std::endl;
#endif 
              break;
             case SMALLER:
              curr = curr.left_child();
#ifdef MICHAL_DEBUG
              std::cout << "--go left" << std::endl;
#endif               
              break;
             case EQUAL:
              if (curr->is_active()) return CURVE;
              curr = curr.left_child();
#ifdef MICHAL_DEBUG
              std::cout << "--go left" << std::endl;
#endif 
              break;
                          
#ifdef CGAL_TD_DEBUG
             default:
              CGAL_assertion(up == LARGER || up == SMALLER || up == EQUAL);
              return Locate_type();
#endif
            }
            break;

#ifdef CGAL_TD_DEBUG
           default:
            CGAL_assertion(res == LARGER || res == SMALLER || res == EQUAL);
            return Locate_type();
#endif

          }
        }
      }

      /*
      top_he = curr->top(); //p_top = &curr->top();
      
      if (!curr->is_active())
      {
        curr = curr.left_child();
        continue;
      }
      Comparison_result cres = traits->compare_curve_end_y_at_x_2_object()
                                                   (ce, top_he->curve());
      if (cres == SMALLER)
      {
        curr = curr.left_child();
        continue;
      }
      else if (cres == LARGER)
      {
        curr = curr.right_child();
        continue;
      }
      else
      {  
        // ce is on the CURVE (top_he->curve() = ce.cv()) itself 
#ifndef CGAL_TD_DEBUG      
         CGAL_warning(
          (cres == EQUAL) &&
          (traits->compare_curve_end_x_2_object()
                        (ce, Curve_end(top_he,ARR_MAX_END)) != LARGER) &&
          (traits->compare_curve_end_x_2_object()
                        (ce, Curve_end(top_he,ARR_MIN_END)) != SMALLER));
#else
              
       CGAL_assertion(
          (cres == EQUAL) &&
          (traits->compare_curve_end_x_2_object()
                        (ce, Curve_end(top_he,ARR_MAX_END)) != LARGER) &&
          (traits->compare_curve_end_x_2_object()
                        (ce, Curve_end(top_he,ARR_MIN_END)) != SMALLER));
#endif
        if (!p_cv) //if p_cv was not given
        {
          // For a vertical curve, we always visit it after visiting
          // one of its endpoints.
          if ((up == EQUAL) || traits->is_vertical(*curr)) {
            //std::cout << "EQUAL or VERTICAL" << std::endl;
            if (curr->is_active()) return CURVE;
            curr = curr.left_child(); //curr should be active
          }
          else if (up == LARGER) {
            curr = curr.right_child();
          }
          else { // if (up==SMALLER) 
            curr = curr.left_child();
          }
          continue;
        }
        else //if p_cv was given
        {
          //   ce must be interior since there is already a 
          //   curve te which ce is one of its ends. and a curve end can be 
          //   an end point of more than one curves only if it's interior.

          CGAL_assertion(
                   (traits->parameter_space_in_x_2_object()
                                    (ce.cv(), ce.ce()) == ARR_INTERIOR) &&
                   (traits->parameter_space_in_y_2_object()
                                    (ce.cv(), ce.ce()) == ARR_INTERIOR) );


          bool is_min_equal = traits->equal_curve_end_2_object()
                                      (Curve_end(*p_cv,ARR_MIN_END), 
                                       Curve_end(top_he,ARR_MIN_END));

          bool is_max_equal = traits->equal_curve_end_2_object()
                                      (Curve_end(*p_cv,ARR_MAX_END), 
                                       Curve_end(top_he,ARR_MAX_END));

#ifndef CGAL_TD_DEBUG          
          CGAL_warning (is_min_equal || is_max_equal);
#else
           if (!is_min_equal && !is_max_equal)
          {
            std::cerr << "\ntop_he->curve() " << top_he->curve();
            std::cerr << "\np_cv "  << *p_cv << std::endl;
            CGAL_assertion(is_min_equal || is_max_equal);
          }
#endif
          //the underlying point of ce
          const Point& p = (ce.ce() == ARR_MIN_END) ?
                           traits->construct_min_vertex_2_object()(ce.cv()) :
                           traits->construct_max_vertex_2_object()(ce.cv()) ;
          
          Comparison_result res = 
            is_min_equal ? 
              traits->compare_cw_around_point_2_object()
                       (top_he->curve(), is_edge_to_right(top_he,p),
                        *p_cv, is_curve_to_right(*p_cv,p), p)  :
              traits->compare_cw_around_point_2_object()
                       ((*p_cv), is_curve_to_right(*p_cv,p),
                        top_he->curve(), is_edge_to_right(top_he,p), p ,false);
    
          switch(res)
          {
           case LARGER:
            curr = curr.right_child();
            break;
           case SMALLER:
            curr = curr.left_child();
            break;
           case EQUAL:
            switch(up)
            {
             case LARGER:
              curr = curr.right_child();
              break;
             case SMALLER:
              curr = curr.left_child();
              break;
             case EQUAL:
              if (curr->is_active()) return CURVE;
              curr = curr.left_child();
              break;
                          
#ifdef CGAL_TD_DEBUG
             default:
              CGAL_assertion(up == LARGER || up == SMALLER || up == EQUAL);
              return Locate_type();
#endif
            }
            break;

#ifdef CGAL_TD_DEBUG
           default:
            CGAL_assertion(res == LARGER || res == SMALLER || res == EQUAL);
            return Locate_type();
#endif

          }
        }
      }*/
    }
    else
    {
      // if is_degenerate() == 0, meaning: the trapezoid (curr)
      // is neither a point nor a curve , but a real trapezoid
      if (curr->is_active())
        return curr->is_on_boundaries() ? UNBOUNDED_TRAPEZOID : TRAPEZOID;
      curr = curr.left_child();
#ifdef MICHAL_DEBUG
      std::cout << "--go left" << std::endl;
#endif
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

#ifndef CGAL_TD_DEBUG
  
  CGAL_warning(traits != NULL);
  
#else
  
  CGAL_assertion(traits != NULL);
  
#endif
  
  if (right > left)
  {
    int d = (int)std::floor((double(right+left))/2);
    // Replacing operator [] of map with find to please MSVC 7
    Vertex_const_handle v = (ar.find(d)->second)->right();

    Dag_node curr(X_trapezoid
                    (&v,&v,0,0,CGAL_TD_VERTEX,false,false,false,false),
                  container2dag(ar,left,d,num_of_new_nodes),
                  container2dag(ar,d+1,right,num_of_new_nodes));
    num_of_new_nodes++;
    curr.left_child()->set_dag_node(&curr.left_child());
    curr.right_child()->set_dag_node(&curr.right_child());
    curr->set_dag_node(&curr);// fake temporary node
    curr->remove(); // mark as deleted
    curr->set_dag_node(0);

    return curr;
  }
  else
  {  // Replacing operator [] of map with find to please MSVC 7
    return ar.find(left)->second;
  }
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
typename Trapezoidal_decomposition_2<Td_traits>::X_trapezoid 
Trapezoidal_decomposition_2<Td_traits>
::insert(Halfedge_const_handle he) //::insert_in_face_interior(Halfedge_const_handle he)
{
#ifdef MICHAL_DEBUG
  std::cout << "INSERTING: --------------------------" << std::endl ;
   print_cv_data(he->curve());
#endif

  //Td_map_item mapitem = Td_halfedge();
#ifdef CGAL_TD_DEBUG
  std::cout << "\nTD::insert_in_face_interior(" << he << ") called with "
            << (is_valid(*m_dag_root) ? "valid" : "invalid") << " data structure"
            <<  std::endl;
  write(std::cout,*m_dag_root,*traits) << std::endl;
#endif

  if (m_needs_update) update();

  // locate the input Halfedge end points in the X_trapezoid Dag

  CGAL_assertion(traits != NULL);

  //get the two vertices of the halfedge
  Vertex_const_handle v1 = he->min_vertex();
  Vertex_const_handle v2 = he->max_vertex();

  //define the Curve end points (curve end = vertex)
  const Curve_end ce1(he, ARR_MIN_END); //MICHAL: to be removed?
  const Curve_end ce2(he, ARR_MAX_END); //MICHAL: to be removed?
  
  // make sure that the two endpoints  are not the same point
#ifndef CGAL_TD_DEBUG

  CGAL_warning(!traits->equal_curve_end_2_object()(ce1, ce2));
  
#else
  
  CGAL_precondition(!traits->equal_curve_end_2_object()(ce1, ce2));
  
#endif
  
  Locate_type lt1,lt2;

  //should hold the trapezoids in which the edge endpoints should be located
  X_trapezoid* tr1 = NULL; 
  X_trapezoid* tr2 = NULL;

#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
  
  locate_optimization(ce1,tr1,lt1);
  
#else
  //location of the left endpoint of the edge we're inserting
  tr1=&locate(ce1,lt1);
  
#endif
  
  //the inserted edge should not cut any existing edge
  if (lt1 == CURVE)
  {
    CGAL_precondition_msg(lt1!=CURVE,"Input is not planar as\
      one of the input point inside previously inserted Halfedge.");
    return X_trapezoid();
  }
  
  //if the edge starts at vertex, we should not insert it into the DAG, 
  //but we should update all the edges incident to the vertex. 
  //else if this is a new vertex - insert a node to the DAG that will represent the new vertex. 
  //the incident edges in this case is only the edge itself, and so it is a trivial operation.
  X_trapezoid& t_p1=
    (lt1==POINT) ?
    set_trp_params_after_halfedge_insertion(he,ce1,tr1,lt1) :
    insert_curve_at_point_using_dag(he,v1,tr1,lt1);
  
#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
  
  locate_optimization(ce2,tr2,lt2);
  locate_opt_empty();
  
#else
  // TODO(oren): locating the second endpoint. this is not necessary,
  // and time consuming. 
  tr2=&locate(ce2,lt2);
  
#endif
  
  if (lt2==CURVE)
  {
    CGAL_precondition_msg(lt2!=CURVE,"Input is not planar as\
      one of the input point inside previously inserted Halfedge.");
    return X_trapezoid();
  }
  
  X_trapezoid& t_p2= (lt2==POINT) ?
    set_trp_params_after_halfedge_insertion(he,ce2,tr2,lt2) :
    insert_curve_at_point_using_dag(he,v2,tr2,lt2);
  
  // locate and insert end points of the input halfedge to the X_trapezoid
  // Dag if needed
  Dag_node tt_p1(*t_p1.dag_node()); //MICHAL: is it ok to add & ?
  Dag_node tt_p2(*t_p2.dag_node()); //MICHAL: is it ok to add & ?
  
  // create the X_trapezoid iterator for traveling along the Trapezoids that
  // intersect the input Halfedge, using left-low to right-high order
  In_face_iterator it = follow_curve(tt_p1,he,LARGER);
  X_trapezoid* curr = 0;
  X_trapezoid* prev = &t_p1;
  X_trapezoid* prev_bottom = 0;
  X_trapezoid* prev_top = 0;
  X_trapezoid* old_output = it.operator->();
  X_trapezoid* old_top = 0;
  X_trapezoid* old_bottom = 0;
  
#ifndef CGAL_TD_DEBUG

  CGAL_warning(!traits->is_degenerate(*old_output));

#else
  
  CGAL_assertion(!traits->is_degenerate(*old_output));
  
#endif
  
  old_output = 0;
  Dag_node* tt = 0;
  bool first_time = true;

  while(!!it) //this means as long as the iterator is valid
  {
    curr = it.operator->();
    prev_bottom = curr->lb();
    prev_top = curr->lt();
    // pass using it along cv
    it++;             //this is the logic of the iterator.
                      // the iterator goes to the next trapezoid right-high.
    tt = curr->dag_node();
    if(first_time)
    {
      
#ifndef CGAL_TD_DEBUG
      
      if(!curr->is_on_top_boundary() && 
         ((curr->top() == he) || (curr->top() == he->twin()))) //MICHAL: he comp
      {
        CGAL_warning((curr->top() != he) && (curr->top() != he->twin()));//MICHAL: he comp
        return X_trapezoid();
      }
      
#else
      
      CGAL_precondition(curr->is_on_top_boundary() || 
                       ((curr->top() != he) && (curr->top() != he->twin())));//MICHAL: he comp
      
#endif
      
    }
    
    split_trapezoid_by_halfedge (*tt, old_output, old_bottom, old_top, he);
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_assertion(((**tt).top() == he) || ((**tt).top() == he->twin()));//MICHAL: he comp
    
#endif
    
    if(first_time)
    {
      set_neighbours_after_halfedge_insertion(*old_output,t_p1);
      first_time = false;
    }
    if (tt->is_inner_node())
    {
      // merge adjacent trapezoids on input Halfedge's bottom side if possible
      //   make sure we try to merge active trapezoids
      CGAL_assertion(!prev_bottom || prev_bottom->is_active()); //tt->left_child() is active
      if(merge_if_possible(prev_bottom, tt->left_child().operator->()))
      {
        tt->set_left_child(*(prev_bottom->dag_node()));
        old_bottom = prev_bottom;
        m_number_of_dag_nodes--; //update number of nodes in the DAG after merge
      }

      // merge adjacent trapezoids on input Halfedge's top side if possible
      //   make sure we try to merge active trapezoids
      CGAL_assertion(!prev_top || prev_top->is_active()); //tt->right_child() is active
      if(merge_if_possible(prev_top, tt->right_child().operator->()))
      {
        tt->set_right_child(*(prev_top->dag_node()));
        old_top = prev_top;
        m_number_of_dag_nodes--; //update number of nodes in the DAG after merge
      }
      
      // update trapezoid's left/right neighbouring relations
      if(!traits->is_degenerate(*prev) && !traits->is_degenerate(*curr))
      {
        curr->set_lb(prev);
        curr->set_lt(prev);
        prev->set_rb(curr);
        prev->set_rt(curr);
      }
    }
    else
    {
      
#ifdef CGAL_TD_DEBUG
      
      CGAL_assertion(curr->is_valid(traits));
      
#endif
      
      break;
    }
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_assertion(curr->is_valid(traits));
    
#endif
    
  }
  
#ifdef CGAL_TD_DEBUG
  
  CGAL_postcondition(traits->is_degenerate_curve(*old_output));
  CGAL_postcondition((old_output->top() == he) || 
                     (old_output->top() == he->twin())); //MICHAL: he comp
  
#endif
  
  set_neighbours_after_halfedge_insertion (*old_output, t_p2);
  m_number_of_curves++;
  
#ifdef CGAL_TD_DEBUG
  write(std::cout,*m_dag_root,*traits) << std::endl;
  std::cout << "\nTD::insert_in_face_interior() exited with data structure" 
            << is_valid(*m_dag_root) << std::endl;
#endif
  

  //print_dag_addresses(*m_dag_root);
  //std:: cout << "Largest leaf depth+1: " << (largest_leaf_depth() + 1) << std::endl;
  //std:: cout << "Number of DAG nodes: " << number_of_dag_nodes() << std::endl;
  //std::cout << "Longest query path: " << longest_query_path_length() << std::endl;

  return *old_output;
}


//-----------------------------------------------------------------------------
// Description:
// 
// Remark:
//  Assumes the map to be planar.
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::remove(Halfedge_const_handle he) //MICHAL: used to be: remove_in_face_interior
{
#ifdef CGAL_TD_DEBUG
  std::cout << "\nTD::remove_in_face_interior(" << he << ") called with "
            << (is_valid(*m_dag_root) ? "valid" : "invalid") << " data structure"
            <<  std::endl;
  write(std::cout,*m_dag_root,*traits) << std::endl;
#endif

  
  if (m_needs_update) update();
  
#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
  locate_opt_empty();
#endif
  
#ifndef CGAL_TD_DEBUG
  CGAL_warning(traits != NULL);
#else
  CGAL_assertion(traits != NULL);
#endif
  
  //calculating leftmost and rightmost curve ends of he
  const Curve_end leftmost(he,ARR_MIN_END);
  const Curve_end rightmost(he,ARR_MAX_END);

  //locating leftmost & rightmost curve ends
  Locate_type lt1,lt2;
  X_trapezoid& t1 = locate(leftmost,lt1);
  X_trapezoid& t2 = locate(rightmost,lt2);
  
  //both should be located on a point degenerate trapezoid
  CGAL_warning(lt1==POINT && lt2==POINT);

  if (!(lt1==POINT && lt2==POINT)) return;
  
#ifndef CGAL_TD_DEBUG
  
  CGAL_warning(t1.dag_node() != NULL);
  CGAL_warning(t2.dag_node() != NULL);
  
#endif
  //retrieve the Dag_nodes of the two point-degenerate trapezoid
  Dag_node& tt1 = *t1.dag_node();
  Dag_node& tt2 = *t2.dag_node();

  //calculate the immediate lower, central and upper neighbourhood of
  // the curve in the data structure
  In_face_iterator btm_it(follow_curve(tt1,he,SMALLER));
  In_face_iterator mid_it(follow_curve(tt1,he,EQUAL));
  In_face_iterator top_it(follow_curve(tt1,he,LARGER));
  
  bool inc_btm = true; //true if btm_it should be be incremented (not top_it)
  bool prev_inc_btm = false; //holds the inc_btm from previous iteration
  bool end_reached = false; //true if this is the last iteration
  
  //define the map of new DAG nodes of the new trapezoids
  Nodes_map new_array; 
  int last_index[] = {0,0};
  int sz = 0;

  Vertex_const_handle left_v = btm_it->left();
  Vertex_const_handle right_v;
  X_trapezoid* last_btm_tr      = NULL; //pointer to the last btm_it tr
  X_trapezoid* last_top_tr      = NULL; //pointer to the last top_it tr
  X_trapezoid* last_new_tr      = NULL; //last new trpz that was created
  X_trapezoid* old_tr           = NULL; //old trpz on which the new is based

  
#ifndef CGAL_TD_DEBUG
  CGAL_warning(traits->equal_curve_end_2_object()
              (top_it->left()->curve_end(), left_v->curve_end()));
#else
  CGAL_precondition(traits->equal_curve_end_2_object()
              (top_it->left()->curve_end(), left_v->curve_end()));
#endif

  //-----------------------------------
  //1. remove adjacency at left end point
  //  first_cv_tr is the first trapezoid representing he (type TD_EDGE)
  const X_trapezoid& first_cv_tr = *mid_it;  

#ifdef CGAL_TD_DEBUG
  CGAL_assertion((first_cv_tr.top() == he) || 
                 (first_cv_tr.top() == he->twin())); //MICHAL: he comp
  CGAL_assertion(traits->equal_curve_end_2_object()
                            (t1.left()->curve_end(), leftmost));
#endif
  
  remove_halfedge_at_vertex_using_geometry(first_cv_tr,t1);
  

  //-----------------------------------
  //2. update the map & the dag with new trapezoids which are merge of the
  //  trapezaoids above and below the removed halfedge.
  do {
    // decide which of btm_it,top_it to increment
    inc_btm = is_end_point_left_low(btm_it->right()->curve_end(),  
                                    top_it->right()->curve_end());
    // the current iterator that should be incremented
    In_face_iterator& curr_it =  inc_btm ? btm_it : top_it;
    // reference to the last curr_it tr
    X_trapezoid*& last_tr = inc_btm ? last_btm_tr : last_top_tr;
    
    //set the new trpz right end
    right_v = curr_it->right();
   
    CGAL_assertion( btm_it->is_active() && top_it->is_active() );

    // create a new trapezoid (the merge of top_it and btm_it)
    typename Nodes_map::value_type
      pair(sz,
           Dag_node(X_trapezoid(&left_v, &right_v, 
                                !btm_it->is_on_bottom_boundary() ?
                                  &btm_it->bottom() : 0,
                                !top_it->is_on_top_boundary() ?
                                  &top_it->top() : 0,
                                curr_it->type_flag(),
                                curr_it->is_on_left_boundary(),
                                curr_it->is_on_right_boundary(),
                                btm_it->is_on_bottom_boundary(),
                                top_it->is_on_top_boundary()) ,
                    std::max(btm_it->dag_node()->depth(),
                             top_it->dag_node()->depth())));
    new_array.insert(pair);
    //copy trapezoid data from btm and top trapezoids
    Dag_node& new_node = (new_array.find(sz))->second;
    ++sz;
    new_node->set_dag_node(&new_node);
    new_node->set_lb(btm_it->lb());
    new_node->set_lt(top_it->lt());
    if (last_new_tr)
    {
      if (traits->is_trpz_top_equal(*last_new_tr,*new_node))
      {
        new_node->set_lt(last_new_tr);
      }

      if (traits->is_trpz_bottom_equal(*last_new_tr,*new_node))
      {
        new_node->set_lb(last_new_tr);
      }
    }
    if (new_node->lb())
      new_node->lb()->set_rb(new_node.operator->());
    if (new_node->lt())
      new_node->lt()->set_rt(new_node.operator->());
    
    last_new_tr = new_node.operator->(); //get the last new trapezoid created

    //update arguments for next iteration:
    left_v = right_v; //new left is curr right
    last_btm_tr = btm_it.operator->(); 
    last_top_tr = top_it.operator->();
    
#ifdef CGAL_TD_DEBUG
    CGAL_warning(last_btm_tr);
    CGAL_warning(last_top_tr);
#endif
    
    old_tr = curr_it.operator->(); //the old trpz on which the new is based
    curr_it++; //increment the iterator to the next trapezoid
    end_reached = !btm_it || !top_it;
    
    //copy neighbouring trapezoids in case top/btm are not the same for the old 
    //  trapezoid and the next trapezoid after incrementing the old one
    if (!btm_it ||
        (inc_btm && !traits->is_trpz_bottom_equal(*old_tr,*curr_it)))
    {
      X_trapezoid* rb = old_tr->rb();
      if (rb) 
      {
        rb->set_lb(last_new_tr);
        last_new_tr->set_rb(rb);
      }
    }
    if (!top_it || 
        (!inc_btm && !traits->is_trpz_top_equal(*old_tr,*curr_it)))
    {
      X_trapezoid* rt = old_tr->rt();
      if (rt) 
      {
        rt->set_lt(last_new_tr);
        last_new_tr->set_rt(rt);
      }
    }

#ifdef CGAL_TD_DEBUG
    CGAL_assertion(last_new_tr->dag_node());
#endif

    //set the no longer relevant trapezoids as removed and add the new nodes
    //  as their replacement
    if (prev_inc_btm != inc_btm)
    {
      int num_of_new_nodes = 0;
      Dag_node tmp =
        container2dag (new_array, last_index[inc_btm ? 0 : 1], 
                        sz-1, num_of_new_nodes);
#ifdef CGAL_TD_DEBUG
      std::cout << "\nremove_in_face_interior allocated ";
      write(std::cout,tmp,*traits) << "\ninto ";
      write(std::cout,*last_tr,*traits,false) << std::endl;
#endif
      last_tr->remove(&tmp);
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

#ifdef CGAL_TD_DEBUG
      std::cout << "\nremove_in_face_interior allocated ";
      write(std::cout,tmp,*traits) << "\ninto ";
      write(std::cout,*last_tr,*traits,false) << std::endl;
#endif
      
      last_tr->remove(&tmp); 
      m_number_of_dag_nodes += 1; //new node (tmp) was added

      
      last_index[inc_btm ? 0 : 1] = sz;
      //tmp is a node with no children
      update_largest_leaf_depth( tmp.max_depth() ); 
      
    }

    // update the dag node pointer in the trapezoid
    const Dag_node* real = &last_tr->dag_node()->left_child();
    (*real)->set_dag_node((Dag_node*)real);
  }
  while(!end_reached);


  // get the iterator (btm_it or top_it) that holds the trapezoid that was 
  //  not removed in the last iteration
  In_face_iterator& it = !prev_inc_btm ? btm_it : top_it;
  
#ifdef CGAL_TD_DEBUG
  CGAL_warning(traits->equal_curve_end_2_object()(it->right(),rightmost));
#endif

  // remove the last trapezoid to remove and set the last new trapezoid
  //  created as its replacement. update the relevant data
  X_trapezoid* rb = it->rb();
  X_trapezoid* rt = it->rt();

  int num_of_new_nodes = 0;
  Dag_node tmp = 
    container2dag(new_array, last_index[!inc_btm ? 0 : 1], 
                  new_array.size()-1, num_of_new_nodes);

#ifdef CGAL_TD_DEBUG
  std::cout << "\nremove_in_face_interior allocated ";
  write(std::cout,tmp,*traits) << "\ninto ";
  write(std::cout,*it,*traits,false) << std::endl;
#endif

  it->remove(&tmp);
  //tmp is the root of a sub graph. 
  //The largest depth in this sub-graph may be the largest leaf depth
  update_largest_leaf_depth( tmp.depth() ); 
  m_number_of_dag_nodes += num_of_new_nodes; //new node (tmp) was added
  
  const Dag_node *real = &it->dag_node()->left_child();
  (*real)->set_dag_node((Dag_node*)real);
  
  if (rb) 
  {
    last_new_tr->set_rb(rb);
    rb->set_lb(last_new_tr);
  }
  if (rt) 
  {
    last_new_tr->set_rt(rt);
    rt->set_lt(last_new_tr);
  }
  
  //-----------------------------------
  //3. remove the trapezoids that represent the removed halfedge
  Base_trapezoid_iterator last_mid = mid_it;
  while(!!++mid_it)
  {
    
#ifdef CGAL_TD_DEBUG
    CGAL_warning(traits->is_degenerate_curve(*last_mid));
#endif
    
    last_mid->remove();
    last_mid=mid_it;
  }
  
  //-----------------------------------
  //4. remove adjacency at right end point
  
#ifdef CGAL_TD_DEBUG
  CGAL_assertion((last_mid->top() == he) || (last_mid->top() == he->twin())); //MICHAL: he comp
  CGAL_assertion(traits->equal_curve_end_2_object()(rightmost,t2.left()));
#endif
  
  remove_halfedge_at_vertex_using_geometry(*last_mid,t2);
  //remove the final trapezoid representing the removed halfedge
  last_mid->remove(); 
  
  //-----------------------------------
  //5. if the halfedge vertices are now isolated, undo the split trapezoid 
  //  by point(vtx) operation
  if (traits->is_isolated_point(t1)) 
    undo_split_trapezoid_by_vertex(tt1,leftmost);
  if (traits->is_isolated_point(t2)) 
    undo_split_trapezoid_by_vertex(tt2,rightmost);

  //-----------------------------------
  //6. reevaluating number of curves
  m_number_of_curves--;
  
#ifdef CGAL_TD_DEBUG
  std::cout << "\nTD::remove_in_face_interior() exited with data structure" 
            << is_valid(*m_dag_root) << std::endl;
  write(std::cout,*m_dag_root,*traits) << std::endl;
#endif

  //rebuild_if_necessary(Are_all_sides_oblivious_tag()); //MICHAL: should be commented!!!
 
  //print_dag_addresses(*m_dag_root);
  //std:: cout << "Largest leaf depth: " << largest_leaf_depth() << std::endl;
  //std:: cout << "Number of DAG nodes: " << number_of_dag_nodes() << std::endl;
 

}


//-----------------------------------------------------------------------------
// Description:
// 
// preconditions:
//  p is not on an edge or a vertex.
template <class Td_traits> 
typename Trapezoidal_decomposition_2<Td_traits>::X_trapezoid&
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
  
  Dag_node curr = *m_dag_root; //MICHAL: is it ok to add &?

#ifdef CGAL_TD_DEBUG
  CGAL_precondition(!!curr);
#endif
  
  lt = search_using_dag(curr, traits, p, NULL, 
                                  up_direction ?
                                  CGAL::LARGER : CGAL::SMALLER);
  X_trapezoid& t_p = *curr;
  
  
#ifdef CGAL_TD_DEBUG
  CGAL_warning(t_p.dag_node());
#endif
  X_trapezoid& tr = **t_p.dag_node();
  
  // tr should be non degenerate trapezoid
  CGAL_assertion(!traits->is_degenerate(tr));
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
  
  if ((up_direction && !tr.is_on_right_boundary() &&
       (traits->compare_x_2_object()(p,tr.right()) == EQUAL) && 
       (tr.is_on_left_boundary() ||
        !traits->equal_curve_end_2_object()(tr.left(),tr.right())))       ||
      (!up_direction && !tr.is_on_left_boundary() &&
       (traits->compare_x_2_object()(p,tr.left()) == EQUAL) && 
       (tr.is_on_right_boundary() ||
        !traits->equal_curve_end_2_object()(tr.left(),tr.right()))))
  {
    // recalculate vertical ray shoot using locate on point
    return up_direction ?
        locate(tr.right(),lt) : locate(tr.left(),lt);
  }
  
  if (up_direction ? tr.is_on_top_boundary() : tr.is_on_bottom_boundary())
  {
    lt = UNBOUNDED_TRAPEZOID;
  }
  else
  {
    Halfedge_handle he = up_direction ? tr.top() : tr.bottom();
    // Now we know that the trapezoid is bounded on the
    // direction of the shoot.
    lt = (traits->equal_curve_end_2_object()(p,Curve_end(he,ARR_MIN_END)) || 
         traits->equal_curve_end_2_object()(p,Curve_end(he,ARR_MAX_END))) ?  
       POINT : CURVE;
  }
  return tr;
}


//-----------------------------------------------------------------------------
// Description:
// 
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::before_split_edge(const X_monotone_curve_2& cv,
                   const X_monotone_curve_2& cv1, 
                   const X_monotone_curve_2& cv2)
{
#ifdef MICHAL_DEBUG
  std::cout << "SPLITTING: --------------------------" << std::endl;
  std::cout << "cv before split" << std::endl;
  print_cv_data(cv);
  std::cout << "cv1 before split" << std::endl;
  print_cv_data(cv1);
  std::cout << "cv2 before split" << std::endl;
  print_cv_data(cv2);
#endif

#ifdef CGAL_TD_DEBUG
  std::cout << "\nTD::before_split_edge(" << cv << "," << cv1 << "," << cv2 
            << ") called with " << (is_valid(*m_dag_root) ? "valid" : "invalid") 
            << " data structure" <<  std::endl;
  write(std::cout,*m_dag_root,*traits) << std::endl;
#endif
  
  if (m_needs_update) update();
  
#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
  
  locate_opt_empty();
  
#endif
    
#ifndef CGAL_TD_DEBUG
  
  if (!traits)
  {
    CGAL_warning(traits != NULL);
    return;
  }
  if (!traits->are_mergeable_2_object()(cv1,cv2))
  {
    CGAL_warning(traits->are_mergeable_2_object()(cv1,cv2));
    return;
  }
  
#else
  
  if (!traits->are_mergeable_2_object()(cv1,cv2))
  {
    std::cerr << "\ncv " << cv;
    std::cerr << "\ncv1 " << cv1;
    std::cerr << "\ncv1 " << cv2 << std::endl;
  }
  CGAL_precondition(traits != NULL);
  CGAL_precondition(traits->are_mergeable_2_object()(cv1,cv2));

#endif

  // find the splitting point (curve end)
  Curve_end split_ce = 
      traits->equal_curve_end_2_object()(Curve_end(cv1, ARR_MAX_END), 
                                        Curve_end(cv2, ARR_MIN_END) ) ?
      Curve_end(cv1, ARR_MAX_END) : Curve_end(cv2, ARR_MAX_END); 
  
  

#ifndef CGAL_TD_DEBUG
  
  CGAL_warning( is_end_point_left_low(Curve_end(cv,ARR_MIN_END),split_ce) );
  
  CGAL_warning( is_end_point_right_top(Curve_end(cv,ARR_MAX_END),split_ce) );
  
#else
  
  CGAL_precondition( is_end_point_left_low(Curve_end(cv,ARR_MIN_END),split_ce) );
  
  CGAL_precondition( is_end_point_right_top(Curve_end(cv,ARR_MAX_END),split_ce) );

#endif
  
  // find extremal points
  const Curve_end leftmost = (traits->equal_curve_end_2_object()
                                (Curve_end(cv1, ARR_MAX_END), 
                                 Curve_end(cv2, ARR_MIN_END) ))?
                            Curve_end(cv1,ARR_MIN_END) : 
                            Curve_end(cv2,ARR_MIN_END) ;

  const Curve_end rightmost = (traits->equal_curve_end_2_object()
                                 (Curve_end(cv1, ARR_MAX_END), 
                                  Curve_end(cv2, ARR_MIN_END) ))?
                            Curve_end(cv2,ARR_MAX_END) : 
                            Curve_end(cv1,ARR_MAX_END) ;

  CGAL_assertion(traits->equal_curve_end_2_object()
                              (Curve_end(cv,ARR_MIN_END), leftmost));
  CGAL_assertion(traits->equal_curve_end_2_object()
                              (Curve_end(cv,ARR_MAX_END), rightmost));

  //locate the trapezoids of the extremal points
  Locate_type lt1,lt2;
  
  // representing trapezoids for extremal points
  X_trapezoid& t1 = locate(leftmost,lt1);
  X_trapezoid& t2 = locate(rightmost,lt2);

#ifndef CGAL_TD_DEBUG
  
  CGAL_warning(lt1==POINT && lt2==POINT);
  CGAL_warning(t1.is_active() && t2.is_active());
  
#else
  
  CGAL_precondition(lt1==POINT && lt2==POINT);
  CGAL_precondition(t1.is_active() && t2.is_active());
  CGAL_warning(t1.dag_node() != NULL);
  CGAL_warning(t2.dag_node() != NULL);
  
#endif
  m_before_split.m_cv_before_split = cv;
  //set iterators for below curve, on curve & above curve
  Dag_node& tt1 = *t1.dag_node();
  m_before_split.m_p_btm_it = new In_face_iterator(follow_curve(tt1,m_before_split.m_cv_before_split,SMALLER));
  m_before_split.m_p_mid_it = new In_face_iterator(follow_curve(tt1,m_before_split.m_cv_before_split,EQUAL));
  m_before_split.m_p_top_it = new In_face_iterator(follow_curve(tt1,m_before_split.m_cv_before_split,LARGER));
  //locate the splitting point in the trapezoidal map
  //  should be found on a degenerate trapezoid representing a curve
  Locate_type lt;
  X_trapezoid& old_t = locate(split_ce,lt);
  
#ifdef CGAL_TD_DEBUG
  
  CGAL_assertion(lt==CURVE);
  CGAL_precondition(old_t.is_active());
  CGAL_warning(old_t.dag_node());
  
#endif
  //the DAG node of the curve trapezoid where the spiltting point is
  Dag_node& old_split_node = *old_t.dag_node();

  CGAL_assertion(traits->equal_curve_end_2_object()
                  (old_t.left()->curve_end(),leftmost)); 
  
  CGAL_assertion(traits->equal_curve_end_2_object()
                  (old_t.right()->curve_end(),rightmost)); 
  

  
  m_before_split.m_p_old_t = &old_t;
  m_before_split.m_p_t1 = &t1;
  m_before_split.m_p_t2 = &t2;
}

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
template <class Td_traits> 
void Trapezoidal_decomposition_2<Td_traits>
::split_edge(const X_monotone_curve_2& cv,Halfedge_const_handle he1, Halfedge_const_handle he2)
{     
  //make sure both halfedges are valid
  CGAL_precondition_code(Halfedge_const_handle invalid_he);
  CGAL_precondition( he1 != invalid_he);
  CGAL_precondition( he2 != invalid_he);

#ifdef CGAL_TD_DEBUG
  std::cout << "\nTD::split_edge(" << cv << "," << he1 << "," << he2 
            << ") called with " << (is_valid(*m_dag_root) ? "valid" : "invalid") 
            << " data structure" <<  std::endl;
  write(std::cout,*m_dag_root,*traits) << std::endl;
#endif
  
 
  // find the splitting point (vertex & curve end)
  Vertex_const_handle split_v = 
      traits->equal_curve_end_2_object()(Curve_end(he1, ARR_MAX_END), 
                                        Curve_end(he2, ARR_MIN_END) ) ?
      he1->max_vertex() : he2->max_vertex(); 
  
  Curve_end ce(split_v->curve_end());

  // find extremal points
  const Curve_end leftmost = (traits->equal_curve_end_2_object()
                                (Curve_end(he1, ARR_MAX_END), 
                                 Curve_end(he2, ARR_MIN_END) ))?
                            Curve_end(he1,ARR_MIN_END) : 
                            Curve_end(he2,ARR_MIN_END) ;

  const Curve_end rightmost = (traits->equal_curve_end_2_object()
                                 (Curve_end(he1, ARR_MAX_END), 
                                  Curve_end(he2, ARR_MIN_END) ))?
                            Curve_end(he2,ARR_MAX_END) : 
                            Curve_end(he1,ARR_MAX_END) ;


  //MICHAL: new
  X_trapezoid& t1 = *m_before_split.m_p_t1;
  X_trapezoid& t2 = *m_before_split.m_p_t2;
  X_trapezoid& old_t = *m_before_split.m_p_old_t;
  In_face_iterator& bottom_it = *m_before_split.m_p_btm_it;
  In_face_iterator& mid_it = *m_before_split.m_p_mid_it;
  In_face_iterator& top_it = *m_before_split.m_p_top_it;
  //MICHAL: new end

  //the DAG node of the curve trapezoid where the spiltting point is
  Dag_node& old_split_node = *old_t.dag_node();
  

  // previous left and right sons of this DAG node
  const Dag_node& old_left = old_split_node.left_child();
  const Dag_node& old_right= old_split_node.right_child();

  //define the left halfedge and the right halfedge, according
  //  to the splitting point
  //Halfedge_const_handle* p_left_he  = NULL;
  //Halfedge_const_handle* p_right_he = NULL;
  Halfedge_const_handle left_he  = he2;
  Halfedge_const_handle right_he = he1;
  
  if (traits->equal_curve_end_2_object() (Curve_end(he2,ARR_MIN_END), ce))
  {
    left_he  = he1; //p_left_he  = &he1;
    right_he = he2; //p_right_he = &he2;
  }
  //else
  //{
  //  p_left_he  = &he2;
  //  p_right_he = &he1;
  //}

  CGAL_assertion(old_t.is_active());

  //updating the data structure:
  //  the cv trpz where the splitting point is located will hold the point 
  //  its new left will hold the left cv (left and right children as before)
  //  its new right will hold the right cv (left and right children asa before)

  //defining the new left child node
  const Dag_node& new_left_node =
          Dag_node(X_trapezoid(old_t.left(), split_v,
                               left_he, left_he, //*p_left_he, *p_left_he,
                               X_trapezoid::TD_EDGE),
                   old_left, old_right);
  
  //defining the new right child node
  const Dag_node& new_right_node =
          Dag_node(X_trapezoid(split_v, old_t.right(),
                               right_he, right_he, //*p_right_he, *p_right_he,
                               X_trapezoid::TD_EDGE),
                   old_left, old_right);

  //defining the new parent node (which is the splitting point):

  //split_v is the ARR_MAX_END of left_te, and ARR_MIN_END of right_te.
  //need to send the halfedge whose source is split_v.
  Halfedge_const_handle btm_he = left_he;
  if (btm_he->direction() == ARR_LEFT_TO_RIGHT)
    btm_he = btm_he->twin();
  //Halfedge_const_handle* p_btm = p_left_he;
  //if ((*p_btm)->direction() == ARR_LEFT_TO_RIGHT)
  //  p_btm = &(*p_btm)->twin();

  Halfedge_const_handle top_he = right_he;
  if (top_he->direction() == ARR_RIGHT_TO_LEFT)
    top_he = top_he->twin();
  //Halfedge_const_handle* p_top = p_right_he;
  //if ((*p_top)->direction() == ARR_RIGHT_TO_LEFT)
  //  p_top = &(*p_top)->twin();

  
  const Dag_node& new_pnt_node =
          Dag_node(X_trapezoid(split_v, split_v, btm_he, top_he, //*p_btm, *p_top, 
                               X_trapezoid::TD_VERTEX), 
                   new_left_node, new_right_node);
  
  X_trapezoid& new_left_t = *new_left_node;
  X_trapezoid& new_right_t= *new_right_node;
  X_trapezoid& new_t      = *new_pnt_node;
  
  // locate trapezoid trees that correspond to the closest
  //   trapezoids above and below ce 
  X_trapezoid* left_top_t   = top_it.operator->();
  X_trapezoid* left_bottom_t= bottom_it.operator->();

  while(is_end_point_left_low(left_top_t->right()->curve_end(),ce))
    left_top_t = left_top_t->rb();

  while(is_end_point_left_low(left_bottom_t->right()->curve_end(),ce))
    left_bottom_t = left_bottom_t->rt();
  
  Dag_node left_top    = *left_top_t->dag_node();
  Dag_node left_bottom = *left_bottom_t->dag_node();

  //replace the old curve cv with the new curves in the leftmost 
  //  and rightmost end points.
  //the curve end ce belongs to cv interior
  set_trp_params_after_split_halfedge_update (left_he , t1, he1, he2);  
  set_trp_params_after_split_halfedge_update (right_he, t2, he1, he2);
 
  //set the point's lb() which is:
  //     the first halfedge adjacent to the point emanating toward up
  //     or right sweeped using a counter clockwise sweep
  //     starting from up direction not including.
  //set the point's rt() which is:   
  //     the first halfedge adjacent to the point emanating toward bottom
  //     or left sweeped using a counter clockwise sweep
  //     starting from bottom direction not including.
  new_t.set_rt (&new_left_t);
  new_t.set_lb (&new_right_t);

  //set lb() and rt(0 of the two halfedges trapezoids
  // rt() is the next clockwise degenerate_curve around 
  //      rightmost end_point (possibly himself)
  // lb() is the next clockwise degenerate_curve around 
  //      leftmost end_point (possibly himself)
  new_left_t.set_lb ((old_t.lb() != &old_t) ? old_t.lb() : &new_left_t);
  new_left_t.set_rt (&new_right_t);
  new_right_t.set_lb(&new_left_t);
  new_right_t.set_rt((old_t.rt() != &old_t)? old_t.rt() : &new_right_t);
  
  // update geometric boundary for trapezoids representing cv
  // first update the trapezoids of the new left curve
  X_trapezoid* prev = 0;
  while(*mid_it != old_t) 
  {
    mid_it->set_top(left_he);
    mid_it->set_bottom(left_he);
    mid_it->set_right(split_v); 

    //MICHAL: added this assertion just to be sure:
    CGAL_assertion(mid_it->type() == X_trapezoid::TD_EDGE);

    //MICHAL: who is prev?
    prev = mid_it.operator->();
    mid_it++;
  }
  if (prev)
  {
    prev->set_rb(&new_left_t); //MICHAL: I don't understand this
  }
  else // new_left_t is leftmost representative for he
  {
    set_neighbours_after_split_halfedge_update (new_left_t, t1, he1, he2);
  }
  if (t1.rt()==&old_t) t1.set_rt(&new_left_t);
  if (t1.lb()==&old_t) t1.set_lb(&new_left_t);
  mid_it++; 
  new_right_t.set_rb(mid_it.operator->()); //MICHAL: what does it do?
  
  prev = 0;
  while(!!mid_it) 
  {
    mid_it->set_top(right_he);
    mid_it->set_bottom(right_he);
    mid_it->set_left(split_v);

    //MICHAL: added this assertion just to be sure:
    CGAL_assertion(mid_it->type() == X_trapezoid::TD_EDGE);

    prev = mid_it.operator->();
    mid_it++;
  }
  if (prev)
  {
    new_right_t.set_rb(old_t.rb()); //MICHAL: I don't understand this
  }
  else // new_right_t is rightmost representative for te
  {
    set_neighbours_after_split_halfedge_update (new_right_t,t2,he1, he2,false);
  }
  if (t2.rt()==&old_t) t2.set_rt(&new_right_t);
  if (t2.lb()==&old_t) t2.set_lb(&new_right_t);
  
  // update geometric boundary for trapezoids below te
  while (*bottom_it != *left_bottom)
  {
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_assertion (traits->equal_2_object()(bottom_it->top()->curve(), cv));
    
#endif
    //MICHAL: added this assertion to see if it fails (if we reach an edge_end)
    CGAL_assertion (bottom_it->type() == X_trapezoid::TD_TRAPEZOID);
    
    bottom_it->set_top(left_he); 
    bottom_it++;
  }

#ifdef CGAL_TD_DEBUG
  
  CGAL_assertion (*bottom_it==*left_bottom);
  
#endif
  
  Dag_node& bottom_tt = *bottom_it->dag_node();
  bottom_it++;
  
#ifdef CGAL_TD_DEBUG
  
  CGAL_assertion(traits->is_in_closure (*bottom_tt, ce));
  
#endif
  
  split_trapezoid_by_vertex (bottom_tt, split_v, btm_he, top_he); //*p_btm, *p_top);
  // set the splitting trapezoid to be the same one that splits the 
  // X_curve'like trapezoid
  *bottom_tt = new_t; 
  // update top curves
  bottom_tt.left_child()->set_top(left_he);
  bottom_tt.right_child()->set_top(right_he);
  // left and right are not neighbours.
  bottom_tt.left_child()->set_rt(0);
  bottom_tt.right_child()->set_lt(0);
      

  while(!!bottom_it)
  {
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_assertion(traits->equal_2_object() (bottom_it->top()->curve(), cv));
    
#endif
    //MICHAL: added this assertion to see if it fails (if we reach an edge_end)
    CGAL_assertion(bottom_it->type() == X_trapezoid::TD_TRAPEZOID);
    
    bottom_it->set_top(right_he); 
    bottom_it++;
  }

  // update geometric boundary for trapezoids above cv
  while (*top_it != *left_top)
  {
    
#ifdef CGAL_TD_DEBUG
    
    CGAL_assertion(traits->equal_2_object() (top_it->bottom()->curve(), cv));
    
#endif
    //MICHAL: added this assertion to see if it fails (if we reach an edge_end)
    CGAL_assertion(top_it->type() == X_trapezoid::TD_TRAPEZOID);
    
    top_it->set_bottom(left_he); 
    top_it++;
  } 
  
#ifdef CGAL_TD_DEBUG
  
  CGAL_assertion(*top_it == *left_top);
  
#endif
  
  Dag_node &top_tt = *top_it->dag_node();
  top_it++;
  
#ifdef CGAL_TD_DEBUG
  
  CGAL_assertion(traits->is_in_closure (*top_tt, ce));
  
#endif

  split_trapezoid_by_vertex (top_tt, split_v, btm_he, top_he);// left_he, right_he);
  // set the splitting trapezoid to be the same one that splits the 
  // X_curve'like trapezoid
  *top_tt = new_t;
  // update bottom side
  top_tt.left_child()->set_bottom(left_he);
  top_tt.right_child()->set_bottom(right_he);
  // left and right aren't neighbours
  top_tt.left_child()->set_rb(0);
  top_tt.right_child()->set_lb(0);
  
  while(!!top_it)
  {
    
#ifndef CGAL_TD_DEBUG
    
    CGAL_warning(traits->equal_2_object()(top_it->bottom()->curve(), cv));
    
#else
    
    if (!traits->equal_2_object()(top_it->bottom()->curve(), cv))
      std::cout << "\ntop_it->bottom()->curve() "<< top_it->bottom()->curve()
                << "\t cv= " << cv;
    CGAL_assertion(traits->equal_2_object()(top_it->bottom()->curve() ,cv));
    
#endif
    //MICHAL: added this assertion to see if it fails (if we reach an edge_end)
    CGAL_assertion(top_it->type() == X_trapezoid::TD_TRAPEZOID);
    
    top_it->set_bottom(right_he); 
    top_it++;
  }
 
  // mark inactive trapezoids
  //  depth of new_pnt_node is updated here (in the remove operation)
  //   and also depth of the sub-DAG rooted at it
  old_t.remove((Dag_node*)&new_pnt_node); 
  update_largest_leaf_depth( new_pnt_node.max_depth() ); //MICHAL: this is a recursive call for the sub-DAG --EXPENSIVE!
  //adding nodes for the splitting-point and the two parts of the split curve
  m_number_of_dag_nodes += 3; 
  old_t.set_curve_for_rem_he(m_before_split.m_cv_before_split); //MICHAL: added this so the trpz will hold the original curve before the split
  
  const Dag_node* p_new       = &old_t.dag_node()->left_child();
  const Dag_node* p_new_left  = &p_new->left_child();
  const Dag_node* p_new_right = &p_new->right_child();
  const Dag_node* p_old_left  = &p_new_left->left_child();
  const Dag_node* p_old_right = &p_new_left->right_child();

  (*p_new)->set_dag_node ((Dag_node*)p_new);
  (*p_new_left)->set_dag_node ((Dag_node*)p_new_left);
  (*p_new_right)->set_dag_node ((Dag_node*)p_new_right);
  (*p_old_left)->set_dag_node ((Dag_node*)p_old_left);
  (*p_old_right)->set_dag_node ((Dag_node*)p_old_right);  

#ifdef CGAL_TD_DEBUG
  
  CGAL_assertion (old_split_node->is_valid(traits));
  CGAL_assertion (new_pnt_node->is_valid(traits));
  CGAL_assertion ((*p_new)->is_valid(traits));
  CGAL_assertion ((*p_new_left)->is_valid(traits));
  CGAL_assertion ((*p_new_right)->is_valid(traits));
  CGAL_assertion ((*p_old_left)->is_valid(traits));
  CGAL_assertion ((*p_old_right)->is_valid(traits));
  CGAL_assertion (top_tt->is_valid(traits));
  CGAL_assertion (bottom_tt->is_valid(traits));
  CGAL_assertion (old_left->is_valid(traits));
  CGAL_assertion (old_right->is_valid(traits));
  CGAL_assertion (traits->is_degenerate_point(**p_new));
  CGAL_assertion (traits->is_degenerate_curve(**p_new_left));
  CGAL_assertion (traits->is_degenerate_curve(**p_new_right));
  CGAL_assertion (traits->equal_curve_end_2_object()
                  (Curve_end((*p_new_right)->bottom(), ARR_MIN_END),
                   (*p_new)->right()) );
  CGAL_assertion (traits->equal_curve_end_2_object()
                  (Curve_end((*p_new_left)->top(), ARR_MAX_END),
                   (*p_new)->left()) );
#endif
      
  // reevaluating number of curves
  m_number_of_curves++;
  
#ifdef CGAL_TD_DEBUG
  std::cout << "\nTD::split_edge() exited with data structure" 
            << is_valid(*m_dag_root) << std::endl;
  write(std::cout,*m_dag_root,*traits) << std::endl;
#endif
  
}




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
  
  if (m_needs_update) update();
  
#ifndef CGAL_NO_TRAPEZOIDAL_DECOMPOSITION_2_OPTIMIZATION
  
  locate_opt_empty();
  
#endif
  
  const X_monotone_curve_2& cv1 = he1->curve();
  const X_monotone_curve_2& cv2 = he2->curve();
#ifndef CGAL_TD_DEBUG
  
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
  
#else

  if (!traits->are_mergeable_2_object() (cv1, cv2))
  {
    std::cerr << "\ncv " << cv;
    std::cerr << "\nhe1->curve() " << cv1;
    std::cerr << "\nhe2->curve() " << cv2 << std::endl;
  }    
  CGAL_assertion(traits != NULL);
  CGAL_precondition(traits->are_mergeable_2_object() (cv1, cv2));
  
#endif 

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
  X_trapezoid& lt_pt_t  = locate (leftmost, lt1);
  X_trapezoid& rt_pt_t  = locate (rightmost, lt2);
  X_trapezoid& mrg_pt_t = locate (ce, lt);
  

  //varifying that all trapezoids are not NULL and are of type POINT
#ifndef CGAL_TD_DEBUG
  
  CGAL_warning (lt_pt_t.dag_node()  != NULL);
  CGAL_warning (rt_pt_t.dag_node()  != NULL);
  CGAL_warning (mrg_pt_t.dag_node() != NULL);
  
#else
  
  CGAL_precondition (lt1==POINT && lt2==POINT && lt==POINT);
  CGAL_precondition (lt_pt_t.is_active() && rt_pt_t.is_active() 
                      && mrg_pt_t.is_active());
  CGAL_assertion (lt_pt_t.dag_node() != NULL);
  CGAL_assertion (rt_pt_t.dag_node() != NULL);
  CGAL_assertion (mrg_pt_t.dag_node() != NULL);
  
#endif
  
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
  
  CGAL_assertion (traits->equal_curve_end_2_object()
                  (lt_pt_t.left(), leftmost));
  CGAL_assertion (traits->equal_curve_end_2_object()
                  (rt_pt_t.left(), rightmost)); 
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
  Dag_node& lt_pt_node  = *lt_pt_t.dag_node();
  Dag_node& mrg_pt_node = *mrg_pt_t.dag_node();

  //set iterators for below left curve, on left curve & above left curve
  In_face_iterator
    btm_left_it (follow_curve (lt_pt_node, *p_left_cv, SMALLER)),
    mid_left_it (follow_curve (lt_pt_node, *p_left_cv, EQUAL)),
    top_left_it (follow_curve (lt_pt_node, *p_left_cv, LARGER));
  
  //set iterators for below right curve, on right curve & above right curve
  In_face_iterator
    btm_right_it (follow_curve (mrg_pt_node, *p_right_cv, SMALLER)),
    mid_right_it (follow_curve (mrg_pt_node, *p_right_cv, EQUAL)),
    top_right_it (follow_curve (mrg_pt_node, *p_right_cv, LARGER));
  
    
#ifdef CGAL_TD_DEBUG
  
  CGAL_assertion (btm_left_it.operator->());
  CGAL_assertion (mid_left_it.operator->());
  CGAL_assertion (top_left_it.operator->());
  CGAL_assertion (btm_right_it.operator->());
  CGAL_assertion (mid_right_it.operator->());
  CGAL_assertion (top_right_it.operator->());
  CGAL_assertion (btm_left_it->is_active());
  CGAL_assertion (mid_left_it->is_active());
  CGAL_assertion (top_left_it->is_active());
  CGAL_assertion (btm_right_it->is_active());
  CGAL_assertion (mid_right_it->is_active());
  CGAL_assertion (top_right_it->is_active());
  
#endif
  
  //replacing old curves with the new merged halfedge

  X_trapezoid*  leftmost_cv_t = mid_left_it.operator->();
  X_trapezoid*  mid_left      = 0;
  X_trapezoid*  mid_right     = mid_right_it.operator->();
  X_trapezoid*  top_left      = 0; //MICHAL: rename
  X_trapezoid*  top_right     = top_right_it.operator->(); //MICHAL: rename
  X_trapezoid*  btm_left      = 0; //MICHAL: rename
  X_trapezoid*  btm_right     = btm_right_it.operator->(); //MICHAL: rename
  X_trapezoid*  rightmost_cv_t= 0;
  X_trapezoid*  dummy         = 0; //MICHAL: rename
  X_trapezoid*  dummy2        = 0; //MICHAL: rename
  
  Vertex_const_handle leftmost_v  = leftmost_cv_t->left();
  Vertex_const_handle rightmost_v = mid_right->right();

  //replacing the given curve with a new Halfedge_handle along the trapezoids
  // starting at the iterator, until the end (last parameter) is reached. 
  // updating the last param as the last updated trapzoid
  set_trp_params_after_halfedge_update (mid_left_it, left_he, merged_he,
                                        leftmost_v, rightmost_v, mid_left);
  set_trp_params_after_halfedge_update (mid_right_it, right_he, merged_he,
                                      leftmost_v, rightmost_v, rightmost_cv_t);
  set_trp_params_after_halfedge_update (top_left_it, left_he, merged_he,
                                        leftmost_v, rightmost_v, top_left);
  set_trp_params_after_halfedge_update (top_right_it, right_he, merged_he,
                                        leftmost_v, rightmost_v, dummy);
  set_trp_params_after_halfedge_update (btm_left_it, left_he, merged_he,
                                        leftmost_v, rightmost_v, btm_left);
  set_trp_params_after_halfedge_update (btm_right_it, right_he, merged_he,
                                        leftmost_v, rightmost_v, dummy2);
  
  
  // merge trapezoids that were split by the upward and downward
  // vertical extensions from ce (the merged point)
  
  // make sure only active trapezoids are merged
  CGAL_assertion( top_left->is_active() && top_right->is_active() );
  CGAL_assertion( btm_left->is_active() && btm_right->is_active() );

#ifndef CGAL_TD_DEBUG
  
  merge_if_possible (top_left, top_right);
  merge_if_possible (btm_left, btm_right);

#else
  
  CGAL_warning (top_left);
  CGAL_warning (top_right);
  CGAL_warning (merge_if_possible (top_left,top_right));
  CGAL_warning (btm_left);
  CGAL_warning (btm_right);
  CGAL_warning (merge_if_possible (btm_left, btm_right));
  
  
#endif

  // mark older trapezoids as inactive - nodes depth are updated here
  top_right->remove(top_left->dag_node());
  btm_right->remove(btm_left->dag_node());
  update_largest_leaf_depth(std::max(top_left->dag_node()->depth(),
                                      btm_left->dag_node()->depth()));
  //no need to update m_number_of_dag_nodes because the number of nodes did not change.


#ifdef CGAL_TD_DEBUG
  
  CGAL_warning (mid_left);
  CGAL_warning (mid_right);
  CGAL_warning (tt->is_active());
  
#endif
  
  // make the merged point's representative inactive
  mrg_pt_node->remove();
  
  //MICHAL: added this assertion to see if it fails 
  CGAL_assertion(mid_left->type() == X_trapezoid::TD_EDGE);
  
  mid_left->set_rb(mid_right);
  mid_left->set_right(mid_right->right());
  
  //MICHAL: added this assertion to see if it fails 
  CGAL_assertion(mid_right->type() == X_trapezoid::TD_EDGE);
 
  mid_right->set_left(mid_left->left());
  mid_left->set_rt(0);
  mid_right->set_lb(0);

  //replacing the curve in the end points' trapezoids themselves (updating top/ bottom)
  set_trp_params_after_halfedge_update (*p_left_cv, merged_he, lt_pt_t);
  set_trp_params_after_halfedge_update (*p_right_cv, merged_he, rt_pt_t); //MICHAL: maybe I should pass the he1 & he2?
  
#ifdef CGAL_TD_DEBUG
  
  CGAL_warning(leftmost_cv_t  != NULL);
  CGAL_warning(rightmost_cv_t != NULL);
  
#endif

  //updating the connection between the edge trapezoids and the end points trapezoids
  set_neighbours_after_merge_halfedge_update(*leftmost_cv_t,lt_pt_t,cv,true);
  set_neighbours_after_merge_halfedge_update(*rightmost_cv_t,rt_pt_t,cv,false);

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
  //if node represents a curve or trapezoid
  if (!traits->is_degenerate_point(*node) )
    return (1 + std::max(
                  longest_query_path_length_rec(minus_inf, min_node,
                                                plus_inf, max_node,
                                                node.left_child()) ,
                  longest_query_path_length_rec(minus_inf, min_node,
                                                plus_inf, max_node,
                                                node.right_child()) ));
  //if this node represents a point
  //check if it is within param min & max
  
  bool pnt_exists = false;
  if (!node->is_active() && !node->is_on_boundaries())
    pnt_exists = true;
    
  // extract point (curve_end) from trapezoid
  const Curve_end vtx_ce(node->is_active()? 
      node->left()->curve_end() : node->curve_end_for_rem_vtx());
  Point vtx_p;
  if (pnt_exists)
    vtx_p = node->point_for_inner_rem_vtx();

  
  //check if not smaller than min
  if (!minus_inf)
  {
    // extract point (curve_end) from trapezoid
    const Curve_end min_ce(min_node->is_active()? 
       min_node->left()->curve_end() : min_node->curve_end_for_rem_vtx());

    //if smaller than the point represented by min_node 
    if ((!pnt_exists && is_end_point_left_low(vtx_ce, min_ce)) ||
        (pnt_exists &&  is_end_point_left_low(vtx_p, min_ce)) )
      return 0;
  }
  
  //check if not larger than max
  if (!plus_inf)
  {
    // extract point (curve_end) from trapezoid
    const Curve_end max_ce(max_node->is_active()? 
       max_node->left()->curve_end() : max_node->curve_end_for_rem_vtx());

    //if larger than the point represented by max_node 
    if ((!pnt_exists && is_end_point_right_top(vtx_ce, max_ce)) ||
        (pnt_exists &&  is_end_point_right_top(vtx_p, max_ce)) )
      return 0;
  }

  //o/w continue with updated parameters
  return (1 + std::max(
                  longest_query_path_length_rec(minus_inf, min_node,
                                                false, node,
                                                node.left_child()) ,
                  longest_query_path_length_rec(false, node,
                                                plus_inf, max_node,
                                                node.right_child()) ));
}



} //namespace CGAL

#endif
