// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, August 11
//
// file          : include/CGAL/R_Tree/R_tree_index_implementations.h
// package       : ExternalMemoryStructures (0.631)
// maintainer    : Philipp Kramer <kramer@inf.ethz.ch>
// chapter       : $CGAL_Chapter: Basic / External Data Structures $
// source        : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Gabriele Neyer<neyer@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Peter Widmayer <widmayer@inf.ethz.ch>)
//
// Split strategies for the R-tree.
// ======================================================================
#ifndef __R_TREE_index_implementations_H__
#define __R_TREE_index_implementations_H__
#ifdef CGAL_D_DOGN_Index 
  #define CGAL_DOGN_Index(cmd) cmd 
#else
  #define CGAL_DOGN_Index(cmd) 
#endif 


#include <CGAL/basic.h>
#include <vector>
#include <iostream>
#include <CGAL/R_tree.h>
#include <map> 

CGAL_BEGIN_NAMESPACE



// Guttman quadratic split
template <class Container, class I_F>
struct quadratic_split_leaf {
void operator()(Container *first, Container *last, 
		std::back_insert_iterator<std::vector<Container> > left,  
		std::back_insert_iterator<std::vector<Container> > right,
		 int minEntries,  int pagesize) {

  Container *seed1, *seed2, *I1, *I2, *ele, *Itmp;
  double max = -100000;
  double tmp;
  int M, left_nr, right_nr;;
  I_F i_f;
  M = last - first;
  if (M%2) M++;
  // Search the two points with maximal insertion Cost
  for (I1 = first; I1 != last; I1++)  {
    Itmp = I1;
    for (I2 = ++Itmp; I2 != last; I2++) {
      tmp = i_f.cost((*I1).key, (*I2).key);  
      if (tmp > max)  {
        max = tmp;
        seed1 = I1;
        seed2 = I2;
      }
    }
  }
  *left++ = *seed1;
  *right++ = *seed2;
  left_nr = 1;
  right_nr = 1;
   int to_left_size=i_f.written_size((*seed1).ele);
   int to_right_size=i_f.written_size((*seed2).ele);
   int elements= 0;
  std::list<bool> is_deleted;
  for (I1 = first; I1 != last; I1++)  {
    elements++;
    if (I1 != seed1 && I1 != seed2) 
      is_deleted.push_back(false);
    else
      is_deleted.push_back(true);
  }
  typename I_F::Key s1=(*seed1).key;
  typename I_F::Key s2=(*seed2).key;
  bool to_left=false;
  bool to_right=false;
   int act_size;
  double a,b;
  double Min=-10000;
  typename std::list<bool>::iterator del_ele, del_pos;
   int k;
  // insert the element with minimum inefficiency cost
  for (k=1; k< elements-1; k++){
    I1 = first;
    del_ele = is_deleted.begin();
    while(I1 != last && del_ele!=is_deleted.end())
      {
	//    for (I1 = first, del_ele = is_deleted.begin(); 
	//	 I1 != last, del_ele!=is_deleted.end(); I1++, del_ele++)  {
	if (I1 != seed1 && I1 != seed2 && !(*del_ele)) 
	  {
	    a = i_f.cost(s1,(*I1).key) ;
	    b = i_f.cost(s2, (*I1).key);
	    if(Min==-10000 || a<Min || b< Min)
	      {
		if(a==b){
		  if(left_nr<right_nr)
		    b=b+10.0;
		  else
		    a=a+10.0;
		}
		del_pos = del_ele;
		if (a< b){
		  to_left=true;
		  to_right=false;
		  Min = a;
		  ele = I1;
		}
		else 
		  {
		    to_left=false;
		    to_right=true;
		    Min = b;
		    ele = I1;
		  }
	      }
	  }
	I1++;
	del_ele++;
      }
    act_size=i_f.written_size((*ele).ele);
    
    // make sure that enough elements are in each container
    if((left_nr + elements - k - 2 < minEntries) || 
       (right_nr + elements - k - 2 < minEntries))
      {
	if(left_nr + elements - k - 2 < minEntries){
	  to_left_size +=act_size;
	  *left++=*ele;
	  *del_pos =true;
	  s1=i_f.unify(s1, (*ele).key);
	  left_nr++;
	}	
	else{
	  to_right_size+=act_size;
	  *right++=*ele;
	  *del_pos =true;
	  s2=i_f.unify(s2, (*ele).key);
	  right_nr++;
	}
      }
    else
      {
	if(to_left){
	  if(to_left_size + act_size <= pagesize)
	    {
	      to_left_size+=act_size;
	      *left++=*ele;
	      *del_pos =true;
	      s1=i_f.unify(s1, (*ele).key);
	      left_nr++;
	    }
	  else
	    {
	      to_right_size+=act_size;
	      *right++=*ele;
	      *del_pos =true;
	      s2=i_f.unify(s2, (*ele).key);
	      right_nr++;
	    }
	}
	
	if(to_right){
	  if(to_right_size + act_size <= pagesize)
	    {
	      to_right_size+=act_size;
	      *right++=*ele;
	      *del_pos =true;
	      s2=i_f.unify(s2, (*ele).key);
	      right_nr++;
	    }
	  else
	    {
	      to_left_size+=act_size;
	      *left++=*ele;
	      *del_pos =true;
	      s1=i_f.unify(s1, (*ele).key);
	      left_nr++;
	    }
	}
      }
    to_right=to_left=false;
    Min=-10000;
  }
}
};



// same split method as above - for inner vertices of the tree
// Guttman quadratic split
template <class Container, class I_F>
struct quadratic_split_node {
void operator()(Container *first, Container *last, 
		Container *left, Container *right,
		int minEntries, int maxEntries) {
  Container *seed1, *seed2, *I1, *I2, *ele, *Itmp;
  double max = -100000;
  double tmp;
  int left_nr, right_nr;;
  I_F i_f;
  // Search the two points with maximal insertion Cost
  for (I1 = first; I1 != last; I1++)  {
    Itmp = I1;
    for (I2 = ++Itmp; I2 != last; I2++) {
      tmp = i_f.cost((*I1).key, (*I2).key);  
      if (tmp > max)  {
        max = tmp;
        seed1 = I1;
        seed2 = I2;
      }
    }
  }
  *left++ = *seed1;
  *right++ = *seed2; 
  left_nr = 1;
  right_nr = 1;
  int elements= 0;
  std::list<bool> is_deleted;
  for (I1 = first; I1 != last; I1++)  {
    elements++;
    if (I1 != seed1 && I1 != seed2) 
      is_deleted.push_back(false);
    else
      is_deleted.push_back(true);
  }
  typename I_F::Key s1=(*seed1).key;
  typename I_F::Key s2=(*seed2).key;
  bool to_left=false;
  bool to_right=false;
  double a,b;
  double Min=-10000;
  typename std::list<bool>::iterator del_ele, del_pos;
  int k;
  for (k=1; k< elements-1; k++){
    I1 = first; 
    del_ele = is_deleted.begin();
    while(I1 != last && del_ele!=is_deleted.end())
      {
	//    for (I1 = first, del_ele = is_deleted.begin();
	//    I1 != last, del_ele!=is_deleted.end(); ++I1, ++del_ele )  {
	if (I1 != seed1 && I1 != seed2 && *del_ele==false) 
	  {
	    a = i_f.cost(s1,(*I1).key) ;
	    b = i_f.cost(s2, (*I1).key);
	    if(Min==-10000 || a<Min || b< Min){
	      if(a==b)
		{
		  if(left_nr<right_nr)
		    b=b+10.0;
		  else
		    a=a+10.0;
		}
	      del_pos = del_ele;
	      if (a< b)
		{
		  to_left=true;
		  to_right=false;
		  Min = a;
		  ele = I1;
		}
	      else 
		{
		  to_left=false;
		  to_right=true;
		  Min = b;
		  ele = I1;
		}
	    }
	  }
	I1++; 
	del_ele++; 
      }
    // make sure that enough elements are in each container
    if((left_nr + elements - k - 2 < minEntries) || 
       (right_nr + elements - k - 2 < minEntries))
      {
	if(left_nr + elements - k - 2 < minEntries){
	  *left++=*ele;
	  *del_pos =true;
	  s1=i_f.unify(s1, (*ele).key);
	  left_nr++;
	}	
	else{
	  *right++=*ele;
	  *del_pos =true;
	  s2=i_f.unify(s2, (*ele).key);
	  right_nr++;
	}
      }
    else
      {
	if(to_left){
	  if(left_nr+1<=maxEntries)
	    {
	      *left++=*ele;
	      *del_pos =true;
	      s1=i_f.unify(s1, (*ele).key);
	      left_nr++;
	    }
	  else{
	    *right++=*ele;
	    *del_pos =true;
	    s2=i_f.unify(s2, (*ele).key);
	    right_nr++;
	  }
	}
	if(to_right){	  
	  if(right_nr+1<=maxEntries)
	    {
	      *right++=*ele;
	      *del_pos =true;
	      s2=i_f.unify(s2, (*ele).key);
	      right_nr++;
	    }
	}
      }	
    to_right=to_left=false;
    Min=-10000;
  }
}
};



// Beckmann Kriegel... Split -- data entries
template <class Container, class I_F, class Sort_axis>
struct star_split_leaf {
  void operator()(Container *first, Container *last, 
		  std::back_insert_iterator<std::vector<Container> > left,  
		  std::back_insert_iterator<std::vector<Container> > right,
		   int IO_min_cap,  int pagesize) {
    I_F i_f;
    //    sort_axis_data_2_dim<Container> sorting_routine;
    Sort_axis sorting_routine;
    typedef typename I_F::Key Key;
    Key first_group, second_group, best_key_first, best_key_second;
    double max = -100000;
    double tmp,best_dist;
    int best_axis, split_axis=0;
    Container *it_cont;
     int i,j, minEntries, maxEntries,absEntries=0;
    int abs_size=0;
    CGAL_DOGN_Index(
	       std::cerr << "DATA ENTRY SPLIT\n";
	       std::cerr << "Elements that have to be split:";
	       )
      for(it_cont=first;it_cont!=last;it_cont++)
      {
	absEntries++;
	CGAL_DOGN_Index(
		   i_f.dump((*it_cont).key);
		   std::cerr << "size:" <<  (*it_cont).ele.size()
		             << std::endl;
		   )
	abs_size+= (*it_cont).ele.size();
      }
    CGAL_DOGN_Index(
	       std::cerr << "pagesize is " << pagesize << "abs size is "
	                 << abs_size;
	       std::cerr << " end\n";
	       )
    while(sorting_routine(split_axis, first, last))
      {
	CGAL_DOGN_Index(
		   std::cerr << "The new sorting sequence is:\n";
		   for(it_cont=first;it_cont!=last;it_cont++)
		   {
		     i_f.dump((*it_cont).key);
		   }
		   std::cerr << std::endl;
		   )
	//compute minEntries, maxEntries
	minEntries=absEntries; maxEntries=0;
	 int maxSize1=0, maxSize2=0;
	it_cont= first;
	while(maxSize1<pagesize && it_cont!=last)
	  {
	    maxSize1+=(*it_cont).ele.size();
	    maxEntries++;
	    it_cont++;
	  }
	if(maxSize1>pagesize) maxEntries--;
	it_cont=last;
	it_cont--;
	while(maxSize2<pagesize && it_cont!=first)
	  {
	    maxSize2+=(*it_cont).ele.size();
	    minEntries--;
	    it_cont--;
	  }
	if(maxSize2>pagesize) minEntries++;
	if((minEntries < IO_min_cap)&& (absEntries>=2*minEntries))
	  minEntries=IO_min_cap;
	CGAL_DOGN_Index(
		   std::cerr << "MinEntries, MaxEntries, AbsEntries are : "
		             << minEntries;
		   std::cerr << " " << maxEntries << " " << absEntries
		             << std::endl;
		   )
	//compute the values for each solution
	//for(i=minEntries;i<=maxEntries-minEntries;i++)
	for(i=minEntries;i<=maxEntries;i++)
	  {
	    j=0;
	    for(it_cont=first;it_cont!=last;it_cont++)
	      {
		if(j<i)
		  {
		    first_group=i_f.unify(first_group,(*it_cont).key);
		    j++;
		  }
		else
		  second_group=i_f.unify(second_group,(*it_cont).key);
	      }
	    tmp=i_f.cost(i_f.intersection(first_group,second_group));
	    if (tmp < max || max==-100000)
	      {
		max=tmp;
		best_dist = i;
		best_axis=split_axis;
		best_key_first=first_group;
		best_key_second=second_group;
	      }
	    else
	      if(tmp==max)
		{
		  if(i_f.cost(best_key_first)+i_f.cost(best_key_second) > 
		     i_f.cost(first_group)+i_f.cost(second_group))
		    {
		      best_dist = i;
		      best_axis=split_axis;
		      best_key_first=first_group;
		      best_key_second=second_group;
		    }
		}
	  }
	split_axis++;
      }
    //distribute the entries into two groups
    sorting_routine(best_axis, first, last);
    it_cont=first;
    CGAL_DOGN_Index(
	       std::cerr << "output in left list\n";
	       )
    for(i=0;i<best_dist;i++)
      {
	CGAL_DOGN_Index(
		   i_f.dump((*it_cont).key);
		   )
	*left++=*it_cont++;
      }
    CGAL_DOGN_Index(
	       std::cerr << "output in right list\n";
	       )
    while(it_cont!=last)
      {
	CGAL_DOGN_Index(
		   i_f.dump((*it_cont).key);
		   )
	*right++=*it_cont++;
      }
    CGAL_DOGN_Index(
	       std::cerr << "end\n";
	       )

  }
};


// same split method as above - for inner vertices of the tree
// Kriegel Beckmann...  split
template <class Container, class I_F, class Sort_axis>
struct star_split_node {
  void operator()(Container *first, Container *last, 
		  Container *left, Container *right,
		  int minEntries, int maxEntries) {
    //    sort_axis_key_2_dim<Container> sorting_routine;
    Sort_axis sorting_routine;
    I_F i_f;
    Container *it_cont;
    typedef typename I_F::Key Key;
    Key first_group, second_group, best_key_first, best_key_second;
    double max = -100000;
    double best_dist,tmp;
    int i,j,best_axis,split_axis=0;
    while(sorting_routine(split_axis, first, last))
      {
	CGAL_DOGN_Index(
		   std::cerr << "The new INNER sorting sequence is:\n";
		   for(it_cont=first;it_cont!=last;it_cont++)
		   {
		     i_f.dump((*it_cont).key);
		   }
		   std::cerr << std::endl;
		   )
	//compute the values for each solution
	//	for(i=minEntries;i<=maxEntries-minEntries;i++)
	for(i=minEntries;i<=maxEntries;i++)
	  {
	    j=0;
	    for(it_cont=first;it_cont!=last;it_cont++)
	      {
		if(j<i)
		  {
		    first_group=i_f.unify(first_group,(*it_cont).key);
		    j++;
		  }
		else
		  second_group=i_f.unify(second_group,(*it_cont).key);
	      }
	    tmp=i_f.cost(i_f.intersection(first_group,second_group));
	    if (tmp < max || max==-100000)
	      {
		max=tmp;
		best_dist = i;
		best_axis=split_axis;
		best_key_first=first_group;
		best_key_second=second_group;
	      }
	    else
	      if(tmp==max)
		{
		  if(i_f.cost(best_key_first)+i_f.cost(best_key_second) > 
		     i_f.cost(first_group)+i_f.cost(second_group))
		    {
		      best_dist = i;
		      best_axis=split_axis;
		      best_key_first=first_group;
		      best_key_second=second_group;
		    }
		}
	  }
	split_axis++;
      }
    //distribute the entries into two groups
    sorting_routine(best_axis, first, last);
    it_cont=first;
    for(i=0;i<best_dist;i++)
      *left++=*it_cont++;
    while(it_cont!=last)
      *right++=*it_cont++;
  }

};




//Guttmann choose subtree strategy
template <class Container, class I_F>
struct choose_subtree{
  typedef typename I_F::Key Key;
  
  Container* operator()(Key& k, Container *first, Container *last,  int level){
    I_F traits;
    std::list<Container *> ties;
    std::list<Container *>::iterator ti;
    bool ties_true=false;
    double cost = 0, best_cost = 111111, min;
    Container *start_pos = first ,*best_pos;
    while(start_pos!=last && (*start_pos).deleted)
      start_pos++;
    best_cost=traits.cost(k,(*start_pos).key);
    best_pos=start_pos;
    start_pos++;
    while(start_pos!=last)
      {
	if(!(*start_pos).deleted)
	  {
	    cost = traits.cost(k,(*start_pos).key);
	    if(cost < best_cost)
	      {
		best_pos=start_pos;
		best_cost = cost;
		ties_true=false;
	      }
	    else
	      if(best_cost==cost){
		if(ties.back()!=best_pos)
		  { //empty ties first!
		    while(!ties.empty())
		      ties.pop_front();
		    ties.push_back(best_pos);
		  }
		ties.push_back(start_pos);
		ties_true=true;
		best_cost = cost;
		best_pos=start_pos;
		CGAL_DOGN_Index(
			   std::cerr << "Best pos changed to ";
			   traits.dump((*best_pos).key,0);
			   std::cerr << std::endl;
			   )
	      }
	  }
	start_pos++;
      }
      CGAL_DOGN_Index(
                std::cerr << "vor ties true ";
		)
    if(ties_true)
      {
	ti=ties.begin();
	Container * do_it=*(ti);
	min = traits.cost((*do_it).key);
	best_pos=do_it;
	ti++;
	while(ti!=ties.end())
	  {
	    do_it=*(ti);
	    cost=traits.cost((*do_it).key);
	    if(cost<min)
	      {
		min=cost;
		best_pos=do_it;
	      }
	    ti++;
	  }
      }
      CGAL_DOGN_Index(
                std::cerr << "vor return ";
		)
    return best_pos;
  }
};


// den normalen min area choose subtree strategy dazutun
//Beckmann-Kriegel... choose subtree strategy
template <class Container, class I_F>
struct star_choose_subtree{
  typedef typename I_F::Key Key;
  
  Container* operator()(Key& k, Container *first, Container *last,  int level){
    I_F traits;
    std::list<Container *> ties;
    std::list<Container *>::iterator ti;
    bool ties_true=false;
    double cost = 0, best_cost = 111111, min, overlap;
    double best_overlap=-10001;
    Container *i, *start_pos = first ,*best_pos, *act_pos;
    Key act_union, key_inter;
    while(start_pos!=last && (*start_pos).deleted)
      start_pos++;
    best_pos=act_pos=start_pos;
    if(level>1) //child pointers do not point to leaves
      {	
	best_cost=traits.cost(k,(*start_pos).key);
	start_pos++;
	while(start_pos!=last)
	  {
	    if(!(*start_pos).deleted)
	      {
		cost = traits.cost(k,(*start_pos).key);
		if(cost < best_cost)
		  {
		    best_pos=start_pos;
		    best_cost = cost;
		    ties_true=false;
		  }
		else
		  if(best_cost==cost){
		    if(ties.back()!=best_pos)
		      { //empty ties first!
			while(!ties.empty())
			  ties.pop_front();
			ties.push_back(best_pos);
		      }
		    ties.push_back(start_pos);
		    ties_true=true;
		    best_cost = cost;
		    best_pos=start_pos;
		    CGAL_DOGN_Index(
			       std::cerr << "Best pos changed to ";
			       traits.dump((*best_pos).key,0);
			       std::cerr << std::endl;
			       )
		  }
	      }
	    start_pos++;
	  }
	  CGAL_DOGN_Index(
	            std::cerr << "vor ties true";
		    )
	if(ties_true)
	  {
	    ti=ties.begin();
	    Container * do_it=*(ti);
	    min = traits.cost((*do_it).key);
	    best_pos=do_it;
	    ti++;
	    while(ti!=ties.end())
	      {
		do_it=*(ti);
		cost=traits.cost((*do_it).key);
		if(cost<min)
		  {
		    min=cost;
		    best_pos=do_it;
		  }
		ti++;
	      }
	  }
      }
    else
      {
      CGAL_DOGN_Index(
                std::cerr << "im else ";
		)
	for(;start_pos!=last;start_pos++)
	  {
	    if(!(*start_pos).deleted)
	      {
		act_pos=start_pos;
		act_union = traits.unify(k,(*start_pos).key);
		overlap=0;
		i=first;
		for ( ; i!=last;i++) 	    {
		  if((i!=act_pos) && (!(*i).deleted))
		    {
		      key_inter = traits.intersection((*i).key, act_union);
		      overlap += traits.cost(key_inter);
		    }
		}
		CGAL_DOGN_Index(
			   std::cerr << "The overlap : " << overlap
			             << std::endl;
			   )
		if(best_overlap > overlap || best_overlap<-10000){
		  best_overlap = overlap;
		  best_pos=start_pos;
		  ties_true=false;
		  CGAL_DOGN_Index( 
			     std::cerr << "Best pos changed to ";
			     traits.dump((*start_pos).key,0); 
			     std::cerr << std::endl;
			     )
		}
		else
		  if(best_overlap==overlap){
		    if(ties.back()!=best_pos)
		      { //empty ties first!
			while(!ties.empty())
			  ties.pop_front();
			ties.push_back(best_pos);
		      }
		    ties.push_back(start_pos);
		    ties_true=true;
		    best_overlap = overlap;
		    best_pos=start_pos;
		    CGAL_DOGN_Index(
			       std::cerr << "Best pos changed to ";
			       traits.dump((*best_pos).key,0);
			       std::cerr << std::endl;
			       )
		  }
	      }
	  }
	CGAL_DOGN_Index(
	           std::cerr << "vor ties true 2  ";
		   )
	if(ties_true) //choose the position with minimum area
	  {
	    ti=ties.begin();
	    Container * do_it=*(ti);
	    min = traits.cost((*do_it).key,k);
	    best_pos=do_it;
	    ti++;
	    while(ti!=ties.end())
	      {
		do_it=*(ti);
		cost=traits.cost((*do_it).key,k);
		if(cost<min)
		  {
		    min=cost;
		    best_pos=do_it;
		  }
		ti++;
	      }
	  }
      }
      CGAL_DOGN_Index(
                std::cerr << "vor return cccc ";
		)
    return best_pos;
  }
  //*********End Beckmann Kriegel ... **********************

};


template <class Container, class I_F>
struct reinsertion_node {
  typedef typename I_F::Key Key;

  void operator()(Container *first, Container *last, 
			std::back_insert_iterator<std::vector<Container>
			     > new_elements,  
			std::back_insert_iterator<std::vector<Container>
			     > to_reinsert,
			 int minEntries,  int maxEntries) 
  {
    I_F traits;
    std::multimap<double,Container, std::less<double> > to_split;
    std::multimap<double,Container, std::less<double> >::iterator splitIt;
    double dist;
    std::pair<double,Container> p;
    Container *act; 
    Container element;
    Key bound=first->key;
    for(act=first;act!=last;act++)
      bound=traits.unify(bound, act->key);
    for(act=first;act!=last;act++)
      {
	dist=traits.center_dist((*act).key,bound);
	p.first=dist; 
	p.second=*act;
	to_split.insert(p);
      }
    //put the first 1-100/reinsertion_part% of the entries in YY, 
    //the rest should be reinserted
    //we perform close reinsert since it outperformed far reinsert
     int max;
     int reinsertion_part = 30;
    if(0<reinsertion_part < 1){
      double x=(1-(100/reinsertion_part));
      max = ( int)(maxEntries*x);  
    }
    else
      max = ( int)(maxEntries*(0.7)); 
    if(max<minEntries)
      max=minEntries;
     int count=0;
    for(splitIt=to_split.begin();splitIt!=to_split.end();splitIt++)
      {
	CGAL_DOGN_Index(
		   traits.dump((*splitIt).second.key);
		   )
	if(count < max)
	  {
	    (*new_elements)++=(*splitIt).second;
	  }
	else
	  {
	    (*to_reinsert)++=(*splitIt).second;
	  }
	count++;
      }
  }
};

template <class Container, class I_F>
struct reinsertion_leaf {
  typedef typename I_F::Key Key;
  void operator()(Container *first, Container *last, 
		  std::back_insert_iterator<std::vector<Container>
		       > new_elements,  
		  std::back_insert_iterator<std::vector<Container>
		       > to_reinsert,
		   int IO_min_cap,  int pagesize) 
  {
    I_F traits;
    std::multimap<double,Container, std::less<double> > to_split;
    Container *act; 
    Container element;
    double dist; 
    std::pair<double,Container> p;
    Key bound=first->key;
    for(act=first;act!=last;act++)
      bound=traits.unify(bound, act->key);
    for(act=first;act!=last;act++)
      {
	dist=traits.center_dist((*act).key,bound);
	p.first=dist; p.second=*act;
	to_split.insert(p);
      }
    std::multimap<double,Container, std::less<double> >::iterator splitIt;
    //put the first (1-100/reinsertion_part)% of the entries in YY, 
    //the rest should be reinserted
    //we perform close reinsert since it outperformed far reinsert
    double max;
     int reinsertion_part = 30;
    if(0<reinsertion_part < 1){
      double x=(1-(100/reinsertion_part));
      max = x*(pagesize);
    }
    else
      max = 0.7*pagesize;
    int offset=0;
    bool full=false;
    for(splitIt=to_split.begin();splitIt!=to_split.end();splitIt++)
      {
	CGAL_DOGN_Index(
		   traits.dump((*splitIt).second.key);
		   std::cerr << "size: " << (*splitIt).second.size()
		             << std::endl;
		   std::cerr << "offset: " << offset << " max: " << max
		             << std::endl;
		   )
	if(offset + (*splitIt).second.size() <max && !full)
	  {
	    (*new_elements)++=(*splitIt).second;
	    offset += (*splitIt).second.size();
	  }
	else
	  {
	    full=true;
	    (*to_reinsert)++=(*splitIt).second;
	  }
      }
  }
};


template <class Container, class I_F>
struct dummy_reinsertion_node {
  typedef typename I_F::Key Key;

  void operator()(Container *first, Container *last, 
			std::back_insert_iterator<std::vector<Container>
			     > new_elements,  
			std::back_insert_iterator<std::vector<Container>
			     > to_reinsert,
			 int minEntries,  int maxEntries) 
  {

  }
};

template <class Container, class I_F>
struct dummy_reinsertion_leaf {
  typedef typename I_F::Key Key;
  void operator()(Container *first, Container *last, 
		  std::back_insert_iterator<std::vector<Container>
		       > new_elements,  
		  std::back_insert_iterator<std::vector<Container>
		       > to_reinsert,
		   int minEntries,  int pagesize) 
  {
  }
};









CGAL_END_NAMESPACE

#endif




