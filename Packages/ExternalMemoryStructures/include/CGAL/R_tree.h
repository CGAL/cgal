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
// file          : include/CGAL/R_Tree/R_tree.h
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
// Implementation of the R_Tree. Needs file IO_tree_traits and all files that 
// provide the interfaces of file R_Tree_interface.h
// ======================================================================
#ifndef __R_TREE_H__
#define __R_TREE_H__


//*********************************************************************
// Debug output handling
//*********************************************************************
#ifdef CGAL_D_DOGN_Data 
  #define CGAL_DOGN_Data(cmd) cmd 
#else
  #define CGAL_DOGN_Data(cmd) 
#endif 

#ifdef CGAL_D_DOGN_Control 
  #define CGAL_DOGN_Control(cmd) cmd 
#else
  #define CGAL_DOGN_Control(cmd) 
#endif 

#ifdef CGAL_D_DOGN_ControlB
  #define CGAL_DOGN_ControlB(cmd) cmd 
#else
  #define CGAL_DOGN_ControlB(cmd) 
#endif 

#ifdef CGAL_D_DOGN_Insert
  #define CGAL_DOGN_Insert(cmd) cmd 
#else
  #define CGAL_DOGN_Insert(cmd) 
#endif 

#ifdef CGAL_D_DOGN_Storage
  #define CGAL_DOGN_Storage(cmd) cmd 
#else
  #define CGAL_DOGN_Storage(cmd) 
#endif 

#ifdef CGAL_D_DOGN_Search
  #define CGAL_DOGN_Search(cmd) cmd 
#else
  #define CGAL_DOGN_Search(cmd) 
#endif 

#ifdef CGAL_D_DOGN_Delete
  #define CGAL_DOGN_Delete(cmd) cmd 
#else
  #define CGAL_DOGN_Delete(cmd) 
#endif 

#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <list>
#include <map>
#include <vector>
#include <unistd.h>

CGAL_BEGIN_NAMESPACE

// this information suffices to retrieve a node from disk
// it is used especially for the iterator
struct Node_info{
   int level;
   int pnext;
};

   void print_node_infos(std::list<Node_info*>& node_stack)
  {
    std::list<Node_info*>::iterator it;
    std::cerr << "Node_info_list &&&&&&&&&&&\n";
    if(!node_stack.empty())
    {
	for(it=node_stack.begin();it!=node_stack.end();it++)
	  {
	    std::cerr << "level" << (*it)->level << std::endl;
	    std::cerr << "pnext" << (*it)->pnext << std::endl;
	    std::cerr << std::endl;
	  }
     }
  }


// Each entry of an inner node is of this type. 
// pnext points to its offset in the file (not the physical offset).
template <class Traits>
class R_tree_value {
public:
  typedef typename Traits::Key Key;
  typedef R_tree_value *iterator;
  Key key;
  long pnext;
  bool deleted;
  Traits traits;
  R_tree_value(){deleted=true;}
  ~R_tree_value(){}
  R_tree_value(const R_tree_value& t)
  {
    pnext=t.pnext;
    key=t.key;
    deleted=t.deleted;
  }
  R_tree_value  & operator=(const R_tree_value &t)  
  {
    pnext=t.pnext;
    key=t.key;
    deleted=t.deleted;
    return *this;
  }
  int size()
  {
    return (int)(sizeof(long) +sizeof(bool) + traits.written_size_key(key));
  }
  void read(char **s)
  {
    int slong=(int)sizeof(long);
    int sbool=(int)sizeof(bool);
    char *from_long=new char[slong];
    char *from_bool=new char[sbool];
    int i,r=0;
    traits.read_key(s,key);
    r=traits.written_size_key(key);
    for (i=0; i<slong; i++)
       from_long[i] = (*s)[r+i];
    pnext=*((long *)from_long);
    r+=slong;
    for (i=0; i<sbool; i++)
       from_bool[i] = (*s)[r+i];
    deleted=*((bool *)from_bool);
    delete[] from_long;
    delete[] from_bool;
    CGAL_DOGN_Data(
	      std::cerr << "R-tree-value:";
	      traits.dump(key);
	      std::cerr << "--pnext:" << pnext << "--" << deleted << std::endl;
	      )
    
  }

  void write(char **s)
  {
    int slong=(int)sizeof(long);
    int sbool=(int)sizeof(bool);
    char *from_long=new char[slong];
    char *from_bool=new char[sbool];
    int i,r=0;
    traits.write_key(s,key); 
    r=traits.written_size_key(key);
    memcpy(from_long,(char *)(&pnext),sizeof(long));
    for (i=0; i<slong; i++)
      (*s)[i+r] = from_long[i];
    r += slong;
    memcpy(from_bool,(char *)(&deleted),sizeof(bool));
    for (i=0; i<sbool; i++)
      (*s)[i+r] = from_bool[i];
    delete[] from_bool;
    delete[] from_long;
    CGAL_DOGN_Data(
	      //testing
	      Key test;
	      std::cerr << "R-tree-value:";
	      traits.read_key(s,test);
	      traits.dump(test);
	      std::cerr << "=====test keys\n";
	      )
  }

};

// Each data entry is of this type. 
template <class Traits>
class R_tree_leaf_data {
public:
  typedef typename Traits::Data Data;
  typedef typename Traits::Key  Key;
  typedef R_tree_leaf_data *iterator;
  //  bool deleted;
  Data ele;
  Key key;
  Traits traits;
  size_t size(){
    return traits.written_size(ele) + traits.written_size_key(key);
    //    return sizeof(*this);
  }
  R_tree_leaf_data(){}
  ~R_tree_leaf_data(){}
  R_tree_leaf_data(const R_tree_leaf_data& t)
  {
    ele=t.ele;
    key=t.key;
  }
  R_tree_leaf_data  & operator=(const R_tree_leaf_data &t)  
  {
    ele=t.ele;
    key=t.key;
    return *this;
  }
  void read(char** s)
  {
    int r=0;
    traits.read(s,ele);
    r=traits.written_size(ele);
    traits.read_key(&((*s)+r),key);
  }

  void write(char **s)
  {
    int r=0;
    traits.write(s,data); 
    r=traits.written_size(ele);
    traits.write_key(&((*s)+r),key); 
  }
};

//abstract base class defining the base type of a node and of a leaf
template<class v_type, class Data, class Key, class Traits,
  const  int IO_min_cap_nodes,const  int IO_max_cap_nodes, 
  const  int IO_min_cap_leaves, const  int IO_page_size>
class abstract_node_entry
{
public:
   int eles;
  long pthis;		       
   int level;
  abstract_node_entry(){}
  virtual ~abstract_node_entry(){}
  virtual v_type & begin()=0;
  virtual v_type & end()=0;
  virtual bool overfilled( int offset)=0;
  virtual bool underfilled()=0;
  virtual void clear()=0;
  virtual void read(char **s)=0;
  virtual void write (char **s) =0;
  virtual size_t size() =0;
};


// inner node can store Max_cap v_type elements
template<class v_type,class Data, class Key, class Traits, 
  const  int IO_min_cap_nodes, const  int IO_max_cap_nodes, 
  const  int IO_min_cap_leaves, const  int IO_page_size>
class node_entry : public abstract_node_entry<v_type,Data,Key,Traits,
	IO_min_cap_nodes,IO_max_cap_nodes,IO_min_cap_leaves,IO_page_size>
{
public:
  v_type ele[IO_max_cap_nodes];
  
  node_entry(int the_level=1){
    level=the_level;
    pthis=-1;
    //    is_init=true;
    int i;
    for (i=0; i < IO_max_cap_nodes; ++i) 
      {
	ele[i].deleted = true;
      }
    eles = 0;
  }
		     
  void clear(){
    eles = 0;
     int i;
    for (i=0; i < IO_max_cap_nodes; ++i) {
      ele[i].deleted = true;
    }
  }
    
  virtual v_type & begin() { 
    return ele[0];
  }

  virtual v_type & end() { return ele[eles]; }
    
  virtual size_t size() {
    return 2*sizeof(int) + sizeof(long) + IO_max_cap_nodes*ele[0].size();
   }
    
  bool overfilled( int number)
  {
    return (eles <= IO_max_cap_nodes ?false : true);
  }
    
  bool underfilled()
  {
    return (eles < IO_min_cap_nodes ?true : false);
  }
    
    
  virtual void read(char **s)
  {
    int slong=(int)sizeof(long);
    int sint=(int)sizeof(int);
    char *from_int=new char[sint];
    char *from_long=new char[slong];
    int i,r=0;
    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i];
    r += sint;
    eles=*((int *)from_int);

    for (i=0; i<slong; i++)
      from_long[i]=(*s)[i+r];
    r += slong;
    pthis=*((long *)from_long);

    for (i=0; i<sint; i++)
      from_int[i]=(*s)[i+r];
    r += sint;
    level=*((int *)from_int);
    char **k;
    k = new char *;
    int ele_size=ele[0].size();
    for (i=0; i<IO_max_cap_nodes; i++)
      {
	*k=((*s)+r);
	ele[i].read(k);
	r+=ele_size;
      }
    delete[] from_int;
    delete k;
    delete[] from_long;
  }
  
  virtual void write (char **s) 
  { 
    int slong=(int)sizeof(long);
    int sint=(int)sizeof(int);
    char *from_int=new char[sint];
    char *from_long=new char[slong];
    int i,r=0;
    memcpy(from_int,(char *)(&eles),sint);
    for (i=0; i<sint; i++)
      (*s)[i] = from_int[i];
    r += sint;

    memcpy(from_long,(char *)(&pthis),slong);
    for (i=0; i<slong; i++)
      (*s)[i+r] = from_long[i];
    r += slong;

    memcpy(from_int,(char *)(&level),sint);
    for (i=0; i<sint; i++)
      (*s)[i+r] = from_int[i];
    r += sint;

    int ele_size=ele[0].size();
    char **k;
    k = new char *;
    for (i=0; i<IO_max_cap_nodes; i++){
      *k=((*s)+r);
      ele[i].write(k);
      r+=ele_size;
    }
    delete k;
    delete[] from_int;
    delete[] from_long;
  }
};
  
//leaf node can contain data in its field raw_data.
template<class v_type, class Data, class Key, class Traits,
  const  int IO_min_cap_nodes, const  int IO_max_cap_nodes, 
  const  int IO_min_cap_leaves, const  int IO_page_size>
class data_entry : public abstract_node_entry<v_type,Data,Key,Traits,
    IO_min_cap_nodes, IO_max_cap_nodes,IO_min_cap_leaves, IO_page_size>
{
public:
  v_type fake;
  int leaf_offset;
  char raw_data[IO_page_size];
  data_entry(){
    pthis=-1;
    level=0;
    leaf_offset = 0;
    eles = 0;
  }
  
  void clear(){
    eles = 0;   
    leaf_offset = 0;
  }
  virtual size_t size() {
    //    return 3*sizeof(int) + sizeof(long) + IO_page_size*sizeof(char)+2;
    return 3*sizeof(int) + sizeof(long) + IO_page_size*sizeof(char);
  }
  virtual v_type & begin(){
    return fake; }
  virtual v_type & end(){
    return fake;}
  
  void delete_key( int offset,  int num)
  {
    char *tmp1=raw_data + offset;
    char *tmp2=raw_data +offset +num;
    memcpy(tmp1, tmp2, leaf_offset -num -offset);
    leaf_offset = leaf_offset - num;
    eles--;
  }
  bool overfilled( int leafnumber)
  {
    return (available_space(leafnumber)>0 ?false : true);
  }
  
  bool underfilled()
  {
    if (eles<IO_min_cap_leaves)
      return true;
    else
      return false;
  }
  
   int available_space ( int offset) {
    return(IO_page_size - offset);
  }
  
  virtual void read(char **s){ 
    int slong=(int)sizeof(long);
    int sint=(int)sizeof(int);
    char *from_int=new char[sint];
    char *from_long=new char[slong];
    int i,r=0;
    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i];
    r += sint;
    eles=*((int *)from_int);

    for (i=0; i<slong; i++)
      from_long[i]=(*s)[i+r];
    r += slong;
    pthis=*((long *)from_long);

    for (i=0; i<sint; i++)
      from_int[i]=(*s)[i+r];
    r += sint;
    level=*((int *)from_int);

    for (i=0; i<sint; i++)
      from_int[i]=(*s)[i+r];
    r += sint;
    leaf_offset=*((int *)from_int);
    
    for (i=0; i<IO_page_size; i++)
      raw_data[i]=(*s)[i+r];
    r += IO_page_size;
    delete[] from_int;
    delete[] from_long;
    
    CGAL_DOGN_Data(
	      char **ff;
	      ff=new char *;
	      Data element;
	      Key key;
	      Traits traits;
	       int offset=0;
	      std::cerr << "ausgabe von read :" << eles << std::endl;
	      std::cerr << *s << std::endl;
	      for(i=0;i<eles; i++)
	      {
		*ff=raw_data +offset;
		traits.read(ff,element);
		key = traits.build(element);
		traits.dump(key);
		offset += traits.written_size(element);
		std::cerr <<"new offset:" << offset << std::endl;
	      }
	      std::cerr << std::endl;
	      delete ff;
	      )
  }
  
  virtual void write (char **s) {
    int slong=(int)sizeof(long);
    int sint=(int)sizeof(int);
    char *from_int=new char[sint];
    char *from_long=new char[slong];
    int i,r=0;
    memcpy(from_int,(char *)(&eles),sint);
    for (i=0; i<sint; i++)
      (*s)[i] = from_int[i];
    r += sint;
    memcpy(from_long,(char *)(&pthis),slong);
    for (i=0; i<slong; i++)
      (*s)[i+r] = from_long[i];
    r += slong;
    memcpy(from_int,(char *)(&level),sint);
     for (i=0; i<sint; i++)
       (*s)[i+r] = from_int[i];
     r += sint;
     memcpy(from_int,(char *)(&leaf_offset),sint);
     for (i=0; i<sint; i++)
       (*s)[i+r] = from_int[i];
     r += sint;

     for (i=0; i<IO_page_size; i++)
       (*s)[i+r] = raw_data[i];
     r += IO_page_size;
     delete[] from_int;
     delete[] from_long;

     CGAL_DOGN_Data(
	      char **ff;
	      ff=new char *;
	      Data element;
	      Key key;
	      Traits traits;
	      int offset=0;
	      std::cerr << "ausgabe von write :";
	      for(i=0;i<eles; i++)
	      {
		*ff=raw_data + offset;
		traits.read(ff,element);
		key = traits.build(element);
		traits.dump(key);
		offset += traits.written_size(element);
		std::cerr <<"new offset:" << offset << std::endl;
	      }
	      std::cerr << std::endl;
	      data_entry d;
	      d.read(s);
	      offset=0;
	      std::cerr << "ausgabe von read :" << d.eles << std::endl;
	      std::cerr << *s << std::endl;
	      for(i=0;i<d.eles; i++)
	      {
		*ff=d.raw_data + offset;
		traits.read(ff, element);
		key = traits.build(element);
		traits.dump(key);
		offset += traits.written_size(element);
		std::cerr <<"new offset:" << offset << std::endl;
	      }
	      std::cerr << std::endl;
	      delete ff;
	      )
  }
};
  


// Range Tree node: Implementation of a node of the R_Tree.
// We differ between an inner node containing only keys of data
// elements and leaf nodes which contain the data.
// Data_traits: R_Tree_interface
// DATA: Data
// V: R_Tree_node_type<Key>
// D: R_Tree_data_type<Key>
// S1: split function for inner vertices
// S2: split function for leaf vertices
// page_size: size of the buffer for data. 
// E.g. disk block size - sizeof(extra data).
// long name conflict: exchange R_tree_node R_Tree_node
// reinserted when an overflow occurs. 
template <class Data_traits, class R_tree_index,
  class IO_tree_traits_nodes, class IO_tree_traits_leaves, 
   int IO_min_cap_nodes,  
   int IO_max_cap_nodes,  int IO_min_cap_leaves, 
   int IO_page_size>
class R_tree_node {
public:
  typedef R_tree_node Node;
  typedef IO_tree_traits_nodes IO_intf_nodes;
  typedef typename Data_traits::Key Key;
  typedef typename Data_traits::Data Data;
  typedef Data_traits Traits;
  typedef R_tree_value<Traits> v_type;
  typedef R_tree_leaf_data<Traits> d_type;
  
  Traits traits;
  
  //definition of the inner node type depending on the template parameters
  typedef node_entry<v_type, Data, Key, Traits, IO_min_cap_nodes, 
             IO_max_cap_nodes, IO_min_cap_leaves, IO_page_size> Node_entry;
  //definition of the leaf node type
  typedef data_entry<v_type, Data, Key, Traits, IO_min_cap_nodes, 
             IO_max_cap_nodes,IO_min_cap_leaves, IO_page_size> Data_entry;
  
  typedef IO_tree_traits_nodes Tree;
  typedef IO_tree_traits_leaves Leaf_data;
  typedef size_t size_type;
  typedef v_type* iterator;
  typedef iterator const_iterator;
  typedef d_type* d_iterator;
  long data_size;
  long node_size;
  char **data_space;
  char **node_space;
  char **ff;
  
protected:
  Tree *tree;
  Leaf_data *leaf_data;
  typename R_tree_index::Split_node ps;
  typename R_tree_index::Split_leaf ps_leaf;
  typename R_tree_index::Choose_subtree choose_subtree;
  typename R_tree_index::Reinsertion_leaf reinsertion_leaf;
  typename R_tree_index::Reinsertion_node reinsertion_node;

public:
  // XX is always used to store the current node
  abstract_node_entry<v_type,Data,Key,Traits,IO_min_cap_nodes, 
                      IO_max_cap_nodes, IO_min_cap_leaves, IO_page_size> *XX;

  R_tree_node() {
    Data_entry d;
    Node_entry n;
    data_size=d.size();
    node_size=n.size();

    CGAL_DOGN_Data(std::cerr << "Size of data: "
                             << data_size << " size of node: " 
	                     << node_size << std::endl;) 
    data_space=new char *;
    *data_space= new char[data_size];
    node_space=new char *;
    *node_space= new char[node_size];
    ff=new char *;
  }
  
  // A new node is constructed. Either a leaf node or an inner node.
  // If level=0 then a leaf is constructed
  R_tree_node( int level)
  {
    Data_entry d;
    Node_entry n;
    data_size=d.size();
    node_size=n.size();
    CGAL_DOGN_Data(std::cerr << "Size of data: " << data_size << 
              " size of node: " << node_size << std::endl;) 
    data_space=new char *;
    *data_space= new char[data_size];
    node_space=new char *;
    *node_space= new char[node_size];
    ff=new char *;
    if(level==0)
      {
	XX = new Data_entry();
      }
    else
      {
	XX = new Node_entry(level);
      }
  }
  
  // deletion of the node
  virtual ~R_tree_node() 
    { 
      if(dynamic_cast<Data_entry*>(XX))
	delete dynamic_cast<Data_entry*>(XX);
      else
	if(dynamic_cast<Node_entry*>(XX))
	  delete dynamic_cast<Node_entry*>(XX);

      delete[] *data_space;
      delete data_space;
      delete[] *node_space;
      delete node_space;
      delete ff;
      CGAL_DOGN_ControlB(std::cerr << "Delete r_tree_node called\n";)  
    }
  
   void clear() {
    XX->clear();
  }
  
  iterator begin() { 
    return &XX->begin(); 
  }
  iterator end() { 
    return &XX->end();}
  
  
   void open(Tree &t, Leaf_data &l) { tree = &t; leaf_data = &l; }


  // The data of the node (either leaf or not) is read.
  // Therefore, the data at position num in the file is retrieved
  // by a database call.
   bool get( int level, long num) {
    if (level==0)
      {
  	if(dynamic_cast<Data_entry*>(XX)){
	    if ((*leaf_data).get(num,data_space))
	      (*(dynamic_cast<Data_entry*>(XX))).read(data_space);
	    else{
	      CGAL_DOGN_Control(std::cerr << "Could not get leaf data number: " 
			   << num << std::endl;)
	      return false;
	    }
	}
	else
	  {
	    delete dynamic_cast<Node_entry*>(XX);
	    XX = new Data_entry();
	    CGAL_DOGN_Data(std::cerr << "Size of data in get: "<< data_size 
		      << std::endl;)
	    if ((*leaf_data).get(num,data_space))
		(*(dynamic_cast<Data_entry*>(XX))).read(data_space);
	    else
	      {
		CGAL_DOGN_Control(
			     std::cerr << "Could not get leaf data" << num 
			     << std::endl;
			     )
		return false;
	      }
	  }
	CGAL_DOGN_Data(
		  std::cerr << "GET a data entry \n";
		  dump(*(dynamic_cast<Data_entry*>(XX)),0,0);
		  std::cerr << "GET a data entry -- end \n";
		  )
      }
    else{
      if(dynamic_cast<Node_entry*>(XX)){
	  if((*tree).get(num, node_space))
	    (*(dynamic_cast<Node_entry*>(XX))).read(node_space);
	  else
	    {
	      CGAL_DOGN_Control(
			   std::cerr << "Could not get tree data: number" 
			   << num << std::endl;
			   )
	      return false;
	    }
      }
      else
	{   
	  delete dynamic_cast<Data_entry*>(XX);
	  XX = new Node_entry();
	  if((*tree).get(num,node_space))
	      (*(dynamic_cast<Node_entry*>(XX))).read(node_space);
	  else
	    {
	      CGAL_DOGN_Control(
			   std::cerr << "Could not get tree data: number"  
			   << num << std::endl;
			   )
	      return false;
	    }
	}
    }
    return true;
  }

  // The data of the node (either leaf or not) is written.
  // Therefore, the data at position num in the file is written
  // by a database call.
   void put(long num = -1) 
  {
    if(num>-1)
      XX->pthis=num;
    if (XX->level==0)
      { 
	if(dynamic_cast<Data_entry*>(XX))
	  {
	    if (XX->pthis > -1)
	      {
		(*(dynamic_cast<Data_entry*>(XX))).write(data_space);
		(*leaf_data).update(XX->pthis,data_space);
	      }
	    else 
	      {
		//set the number where the data is positioned
		(*(dynamic_cast<Data_entry*>(XX))).pthis = 
		  (*leaf_data).get_pos();
		CGAL_DOGN_Storage(
			     std::cerr << "New pos: " << XX->pthis 
			     << std::endl;
			     )
		  (*(dynamic_cast<Data_entry*>(XX))).write(data_space);
		  (*leaf_data).insert(XX->pthis,data_space);
	      }
	  }
	CGAL_DOGN_Control(
	else
	  std::cerr << "Error wrong Node_entry\n";
	)
      }
    else
      {
	if(dynamic_cast<Node_entry*>(XX))
	  { 
	    if (XX->pthis > -1)
	      {
		(*(dynamic_cast<Node_entry*>(XX))).write(node_space);
		(*tree).update(XX->pthis, node_space);
		(*(dynamic_cast<Node_entry*>(XX))).read(node_space); //test
	      }
	    else 
	      {
		XX->pthis = (*tree).get_pos();
		CGAL_DOGN_Storage(
			     std::cerr << "New pos: " << XX->pthis 
			     << std::endl;
			     )
		(*(dynamic_cast<Node_entry*>(XX))).write(node_space);
		(*tree).insert(XX->pthis,node_space);
		(*(dynamic_cast<Node_entry*>(XX))).read(node_space); //test
	      }
	  }
        CGAL_DOGN_Control(
	else
	  std::cerr << "Error wrong Node_entry\n";)
      }
  }

   // YY is written to the database
  void put_data(Data_entry& YY,long num = -1)
  {
    if(num>-1)
      YY.pthis=num;
    if (YY.pthis > -1)
      {
	YY.write(data_space);
	(*leaf_data).update(YY.pthis,data_space);
      }
    else 
      {
	//set the number where the data is positioned
	YY.pthis = (*leaf_data).get_pos();
	CGAL_DOGN_ControlB(std::cerr << "New pos: " << YY.pthis << std::endl;)
	YY.write(data_space);
	(*leaf_data).insert(YY.pthis,data_space);
      }
  }

   // YY is written to the database
   void put_node(Node_entry& YY, long num=-1)
  {
    if(num>-1)
      YY.pthis=num;

    if (YY.pthis > -1)
      {
	YY.write(node_space);
	(*tree).update(YY.pthis, node_space);
      }
    else 
      {
	YY.pthis = (*tree).get_pos();
	CGAL_DOGN_ControlB(std::cerr << "New pos: " << YY.pthis << std::endl;)
	YY.write(node_space);
	(*tree).insert(YY.pthis,node_space);
      }
  }
	

  
  // the data on the associated file at offset num is marked deleted.
   void erase(long num = -1) {
    if (XX->level==0)
      {
	if (num > -1) XX->pthis = num;
	if (XX->pthis > -1) (*leaf_data).erase(XX->pthis);
      }
    else
      {
	if (num > -1) XX->pthis = num;
	if (XX->pthis > -1) (*tree).erase(XX->pthis);
      }
  }

  //*****************Beckmann Kriegel ******************
  // A data element is inserted into the tree. The best position in 
  // respect to the minimum cost is determined and then the element is
  // either there inserted if the available space is sufficiently large
  // or the node is split.
   void insert_star(Data& d, v_type& v, std::vector<v_type> &children, 
		   std::vector<v_type> &root_children,
		   std::vector<bool>& Reinsert,  int MaxLevel,
		   std::pair< int,long> *LevelPthis)
  {
    CGAL_DOGN_Insert(
		std::cerr << "Insert_star\n";
		if(dynamic_cast<Data_entry*>(XX))
		dump((* dynamic_cast<Data_entry*>(XX)),0,0);
		if(dynamic_cast<Node_entry*>(XX))
		dump_this_node((* dynamic_cast<Node_entry*>(XX)),0,0);
		)
      if (XX->level==0) // leaf insertion
      {
	if(dynamic_cast<Data_entry*>(XX))
	  insert_star(d, children,root_children, Reinsert, MaxLevel, 
                      LevelPthis);
        CGAL_DOGN_Control(
	else
	  std::cerr << "Error dynamic cast \n";)
      }
    else 
      {
	 int pos=0;
	typename std::vector<v_type>::iterator i;
	long XXpthis=XX->pthis;
	 int XXlevel=XX->level;
	CGAL_DOGN_Insert(
		    std::cerr << "pthis: " << XX->pthis << std::endl;
		    )
	if(dynamic_cast<Node_entry*>(XX))
	  {
	    Node_entry YY=*(dynamic_cast<Node_entry*>(XX));
	    //choose the subtree in which to position v
	    iterator pos_it = choose_subtree(v.key, &YY.ele[0],
                                     &YY.ele[IO_max_cap_nodes],YY.level);
	    while (pos_it!=&YY.ele[pos]&&pos<YY.eles) pos++;
	    CGAL_DOGN_Insert(
			std::cerr << "The choosed position is:" 
			<< pos << std::endl;
			)
	  }
	CGAL_DOGN_Control(
	else
	  std::cerr << "Wrong dynamic cast 1\n";
	)
	std::pair< int,long> p(pos,XXpthis);
	LevelPthis[XXlevel]=p;
	(*(dynamic_cast<Node_entry*>(XX))).ele[pos].key= traits.unify(
                      (*(dynamic_cast<Node_entry*>(XX))).ele[pos].key, v.key);
	put();
	// decent the tree, call the insertion routine for the chosen subtree
	get(XX->level-1,(*(dynamic_cast<Node_entry*>(XX))).ele[pos].pnext);
	insert_star(d, v, children, root_children, Reinsert, MaxLevel, 
		    LevelPthis);
	get(XXlevel,XXpthis); // get the actual XX back
	if (!children.empty()) { //the lower tree was split, insert children
	  (*(dynamic_cast<Node_entry*>(XX))).ele[pos]= (*children.begin());
	  children.erase(children.begin());
	  if(dynamic_cast<Node_entry*>(XX)){
	    for (i = children.begin(); i != children.end(); ++i)
	      // the  children are simply put into the vertex
	      insert_star((*i), *(dynamic_cast<Node_entry*>(XX)));
	  }
	  CGAL_DOGN_Control(
	  else
	    std::cerr << "Wrong dynamic cast 2\n";
	  )
	  children.erase(children.begin(), children.end());
	}
	put(XXpthis); //update the vertex
	split_node_star(children,root_children,Reinsert,MaxLevel,LevelPthis);
      }
  }
  //*****************end Beckmann Kriegel ******************

  // special insertion which is only called if the place of insertion
  // was determined before and when there is enough place in the vertex.
   bool insert_special_star( int i, v_type &v, Node_entry& YY)  {
      YY.ele[i] = v;
      YY.eles++;
      put();
      return true;
  }

   // a node is reinserted on the same level
   void  forced_reinsert_node(v_type& v, std::vector<v_type> &children,
			     std::vector<v_type> &root_children, 
			      int  OnLevel, std::vector<bool>& Reinsert,
			      int MaxLevel,
			     std::pair< int,long> *LevelPthis)
  {
    CGAL_DOGN_Insert(
		std::cerr << "FORCED REINSERT NODE\n";
		dump_this_node((* dynamic_cast<Node_entry*>(XX)),0,0);
		)
     int act_level=XX->level;
    long act_pthis=XX->pthis;
    if (XX->level==OnLevel) //the right level to insert the node
      {
	if(dynamic_cast<Node_entry*>(XX))
	insert_star(v,(* dynamic_cast<Node_entry*>(XX)));  
	put();
	split_node_star(children,root_children,Reinsert,MaxLevel, LevelPthis);
      }
    else //a subtree in which to position the entry has to be chosen
      {
	 int pos=0;
	std::vector<v_type>::iterator i;
	if(dynamic_cast<Node_entry*>(XX))
	  {
	    Node_entry YY=*(dynamic_cast<Node_entry*>(XX));
	    iterator pos_it = choose_subtree(v.key, &YY.ele[0],
					     &YY.ele[IO_max_cap_nodes],
					     YY.level);
	    while (pos_it!=&YY.ele[pos]&&pos<YY.eles) pos++;
	    CGAL_DOGN_Insert(std::cerr << "The choosed position is:" << pos 
			<< std::endl;)
	  }
	CGAL_DOGN_Control(
	else
	  std::cerr << "Wrong dynamic cast 1\n";
	)
	(*(dynamic_cast<Node_entry*>(XX))).ele[pos].key= 
	  traits.unify((*(dynamic_cast<Node_entry*>(XX))).ele[pos].key, v.key);
	put();
	get(XX->level-1,(*(dynamic_cast<Node_entry*>(XX))).ele[pos].pnext);
	if(dynamic_cast<Node_entry*>(XX))
	  forced_reinsert_node(v, children, 
			       root_children,
			       OnLevel,Reinsert,MaxLevel,LevelPthis);   
        CGAL_DOGN_Control(
	else
	  std::cerr << "Wrong dynamic cast 3\n";
	)
	get(act_level,act_pthis);
	// may be the reinsertion produces new splits
	CGAL_DOGN_Insert(
		    std::cerr << "Reinsertion produces another split\n ";
		    )
	if (!children.empty()) 
	  { 
	    (*(dynamic_cast<Node_entry*>(XX))).ele[pos]= (*children.begin());
	    children.erase(children.begin());
	    if(dynamic_cast<Node_entry*>(XX)){
	      for (i = children.begin(); i != children.end(); ++i)
		// the  children are simply put into the vertex
		insert_star((*i), *(dynamic_cast<Node_entry*>(XX)));
	    }
	    CGAL_DOGN_Control(
	    else
	      std::cerr << "Wrong dynamic cast 4\n";)
	    children.erase(children.begin(), children.end());
	  }
	put();
	split_node_star(children, root_children ,Reinsert, MaxLevel,
			LevelPthis);
      }
  }
  
   // a node is reinserted. This procedure is called first. The root node
   // is retrieved and then the node is either inserted or the right subtree
   // is searched by calling forced_reinsert
   void  root_forced_reinsert(v_type& v,std::vector<v_type>& root_children, 
			      int  OnLevel, std::vector<bool>& Reinsert,
			      int MaxLevel,
			     std::pair< int,long> *LevelPthis)
  {
    CGAL_DOGN_Insert(
		std::cerr << "ROOT FORCED REINSERT NODE\n";
		dump_this_node((* dynamic_cast<Node_entry*>(XX)),0,0);
		std::cerr << "Rootpthis is:" << LevelPthis[MaxLevel].second 
                          << std::endl;
		)
    std::vector<v_type> act_children;
    int pos;
    get(MaxLevel,LevelPthis[MaxLevel].second);
    if (XX->level==OnLevel)
      {
	CGAL_DOGN_Insert(
		    std::cerr << "Insert in root level" << std::endl;
		    )
	insert_star(v,(* dynamic_cast<Node_entry*>(XX)));  
	put();
	split_node_star(act_children,root_children,Reinsert,MaxLevel, 
			LevelPthis);
	if(!act_children.empty()) {
	  root_children.insert(root_children.begin(),*act_children.begin());
	  act_children.erase(act_children.begin());
	}
	while(!act_children.empty()){
	  root_children.push_back(*act_children.begin());
	  act_children.erase(act_children.begin());
	}

      }
    else
      {
	pos=0;
	if(dynamic_cast<Node_entry*>(XX))    
	  {
	    Node_entry YY=*(dynamic_cast<Node_entry*>(XX));
	    iterator pos_it = choose_subtree(v.key, &YY.ele[0],
					     &YY.ele[IO_max_cap_nodes],
					     YY.level);
	    while (pos_it!=&YY.ele[pos]&&pos<YY.eles) pos++;
	    CGAL_DOGN_Insert(
			std::cerr << "The choosed position is:" << pos 
			<< std::endl;
			)
	  }
	CGAL_DOGN_Control(
	else
	  std::cerr << "Wrong dynamic cast 5\n";)
	(*(dynamic_cast<Node_entry*>(XX))).ele[pos].key= 
	  traits.unify((*(dynamic_cast<Node_entry*>(XX))).ele[pos].key, v.key);
	put();
	get(XX->level-1,(*(dynamic_cast<Node_entry*>(XX))).ele[pos].pnext); 
	forced_reinsert_node(v, act_children,root_children, OnLevel,
			     Reinsert,MaxLevel,LevelPthis);
	if(!act_children.empty()) {
	  root_children.insert(root_children.begin(),*act_children.begin());
	  act_children.erase(act_children.begin());
	}
	while(!act_children.empty()){
	  root_children.push_back(*act_children.begin());
	  act_children.erase(act_children.begin());
	}
	if(!root_children.empty()) {
	   CGAL_DOGN_Insert(
		       std::cerr << "Rootpthis is:" 
		       << LevelPthis[MaxLevel].second << std::endl;
		       )
	  get(MaxLevel,LevelPthis[MaxLevel].second);
	  (*(dynamic_cast<Node_entry*>(XX))).ele[pos]=(*root_children.begin());
	  root_children.erase(root_children.begin());
	  while(XX->eles < IO_max_cap_nodes && !root_children.empty())
	    {
	      (*(dynamic_cast<Node_entry*>(XX))).ele[XX->eles]=
		*(root_children.begin());
	      root_children.erase(root_children.begin());
	      XX->eles++;
	    }
	  //no more childs fit into the root node. Give them back to the root
	  put();
	}
      }
  }



//the elements are sorted according to the reinsertion order. Then the first
//70% are left in the node. The rest has to be reinserted
   void compute_to_split_assign(Node_entry &YY, Key &knew,
		 std::back_insert_iterator<std::vector<v_type> > to_reinsert)
  {
    std::vector<v_type> new_elements;
    std::vector<v_type>::iterator It;
    int count=0;
    reinsertion_node((v_type *) &YY.ele[0], (v_type *) &YY.ele[YY.eles], 
    		     std::back_inserter(new_elements), 
    		     to_reinsert,
    		     IO_min_cap_nodes, IO_max_cap_nodes );    

    YY.clear();
    knew=(*new_elements.begin()).key;
    for(It=new_elements.begin();It!=new_elements.end();It++)
      {
	CGAL_DOGN_Insert(
		    traits.dump((*It).key);
		    )
	YY.ele[count]=(*It);
	YY.eles++; 
	knew = traits.unify(knew, (*It).key);
	count++;
      }
  }

   void reinsert_node( std::vector<v_type>& root_children,
				     int OnLevel, 
				    std::vector<bool>& Reinsert,
				     int MaxLevel,
				    std::pair< int,long> *LevelPthis) {

    std::vector<v_type> to_reinsert;
    std::vector<v_type>::iterator reIt; 
    Key knew;
    CGAL_DOGN_Insert(
		std::cerr << "REINSERT NODE\n";
		dump_this_node((* dynamic_cast<Node_entry*>(XX)),0,0);
    		)
    //for all entries compute the distance between the conters of their 
    //rectangles and the center of the bounding rectangle of it.
    //access to the elements is in increasing order
     int act_level=XX->level;
    long act_pthis=XX->pthis;
    compute_to_split_assign(*(dynamic_cast<Node_entry*>(XX)),knew,
			    std::back_inserter(to_reinsert));
    put(); 
    adapt_parent_keys(LevelPthis,knew,act_level+1,MaxLevel);
    for(reIt=to_reinsert.begin();reIt!=to_reinsert.end();reIt++)
      {
	CGAL_DOGN_Insert(
		    traits.dump((*reIt).key);
		    )
	root_forced_reinsert((*reIt), root_children, OnLevel, 
			     Reinsert,MaxLevel,LevelPthis);
	get(act_level,act_pthis);
      }
  }



  void reinsert_data_split(Data &d, Data_entry &YY, Key &knew, 
			   std::vector<d_type> &to_reinsert)
  {
    CGAL_DOGN_Insert(
		std::cerr << "Reinsert data split\n";
		dump((* dynamic_cast<Data_entry*>(XX)),0,0);
		)
    Data element;
    int i,offset=0;
    std::vector<d_type> new_elements;
    d_type *V;
    V= new d_type[YY.eles+2];
    d_type   val;
    val.key = traits.build(d);
    val.ele = d;
    V[0]=val;
    for(i=0;i<YY.eles; i++)
      {
	*ff=YY.raw_data + offset;
	traits.read(ff,element);
	V[i+1].key = traits.build(element);
	V[i+1].ele = element;
	offset += traits.written_size(element);
      }
    reinsertion_leaf(&V[0], &V[YY.eles+1], std::back_inserter(new_elements), 
		    std::back_inserter(to_reinsert), IO_min_cap_leaves,
		    IO_page_size); 
    delete[] V;
    d_type *dt;  
    YY.leaf_offset=0;
    YY.eles=0;
    knew =new_elements.begin()->key;
    for(dt=new_elements.begin();dt!=new_elements.end();dt++)
      {
	CGAL_DOGN_Insert(
		    traits.dump((*dt).key);
		    std::cerr << "size: " << traits.written_size((*dt).ele) 
		    << std::endl;
		    )

	if(YY.leaf_offset + traits.written_size((*dt).ele) < 
	   IO_page_size)
	  {
	    *ff=YY.raw_data + YY.leaf_offset;
	    traits.write(ff,(*dt).ele);
	    YY.leaf_offset += traits.written_size((*dt).ele);
	    YY.eles++;
	    knew = traits.unify(knew, dt->key);
	    CGAL_DOGN_Insert(
			std::cerr << "\n will stay in leaf: ";
			traits.dump(dt->key);
			)
	  }
	else
	  {
	    CGAL_DOGN_Insert(
	                std::cerr << "reinsertion_leaf has distributed
			              to many data elements in the old
				      leaf -- the rest will be inserted
				      into the other reinsertion vector";
			)
	    to_reinsert.push_back(*dt);
	  }
      }  
    CGAL_DOGN_Insert(
		std::cerr << "the bounding box of the leaf is now:";
		traits.dump(knew);
		)
  }


   void  reinsert_data(Data &d, std::vector<v_type>& root_children,
		       int OnLevel, std::vector<bool>& Reinsert, 
		       int MaxLevel, std::pair< int,long> *LevelPthis) 
  {
    std::vector<d_type> to_reinsert;
    std::vector<d_type>::iterator reIt;
    CGAL_DOGN_Insert(
		std::cerr << "REINSERT DATA\n";
		dump((* dynamic_cast<Data_entry*>(XX)),0,0);
		)
    int act_level=XX->level;
    long act_pthis=XX->pthis;
    Key knew;
    //determine the elements that are to be reinserted
    reinsert_data_split(d, *(dynamic_cast<Data_entry*>(XX)),knew,to_reinsert);
    put();
     int pos;
     int count=0;
    std::pair< int,long> *IndexPthis= new std::pair< int,long>[MaxLevel+1] ;
    CGAL_DOGN_Insert(
		std::cerr << "The new key of the leaf is: " << std::endl;
		traits.dump(knew);
		std::cerr << "Call adapt parent keys\n";
		)
    adapt_parent_keys(LevelPthis,knew,act_level+1,MaxLevel);
    CGAL_DOGN_Insert(
	     std::cerr << "END of adapt parent keys, now reinsert elements\n";
		)
    for(reIt=to_reinsert.begin();reIt!=to_reinsert.end();reIt++)
      {
	CGAL_DOGN_Insert(
		    std::cerr << "\n We now reinsert: ";
		    traits.dump((*reIt).key);
		    std::cerr << "Rootpthis is:" << 
		                 LevelPthis[MaxLevel].second << std::endl;
		    )
	std::vector<v_type> act_children;
	get(MaxLevel,LevelPthis[MaxLevel].second);
	v_type v;
	v.deleted = false;
	v.key= (*reIt).key;
	//use now ordinary insert since it has not to be 
	//inserted on a special level
	pos=0;
	if(dynamic_cast<Node_entry*>(XX))
	  {
	    Node_entry YY=*(dynamic_cast<Node_entry*>(XX));
	    iterator pos_it = choose_subtree(v.key, &YY.ele[0],
					     &YY.ele[IO_max_cap_nodes],
					     YY.level);
	    while (pos_it!=&YY.ele[pos]&&pos<YY.eles) pos++;
	    CGAL_DOGN_Insert(
			std::cerr << "The choosed position is:" << 
			pos << std::endl;
			)
	  }
	CGAL_DOGN_Control(
		    else
		    std::cerr << "Wrong dynamic cast 7\n";
		    )
	(*(dynamic_cast<Node_entry*>(XX))).ele[pos].key= 
	  traits.unify((*(dynamic_cast<Node_entry*>(XX))).ele[pos].key, v.key);
	put();
	get(XX->level-1,(*(dynamic_cast<Node_entry*>(XX))).ele[pos].pnext);
	insert_star((*reIt).ele, v, act_children,root_children,Reinsert,
		    MaxLevel,IndexPthis);
	if(!act_children.empty()) {
	  root_children.insert(root_children.begin(),*act_children.begin());
	  act_children.erase(act_children.begin());}
	while(!act_children.empty()){
	  root_children.push_back(*act_children.begin());
	  act_children.erase(act_children.begin());
	}
	if(!root_children.empty()) {
	  CGAL_DOGN_Insert(
		      std::cerr << "Rootpthis is:" 
		      << LevelPthis[MaxLevel].second << std::endl;
		      )
	  get(MaxLevel,LevelPthis[MaxLevel].second);
	  if(dynamic_cast<Node_entry*>(XX))
	    {
	      (*(dynamic_cast<Node_entry*>(XX))).ele[pos]=
		                           *(root_children.begin());
	      root_children.erase(root_children.begin());
	      while(XX->eles < IO_max_cap_nodes && !root_children.empty())
		{
		  (*(dynamic_cast<Node_entry*>(XX))).ele[XX->eles]=
		    *(root_children.begin());
		  root_children.erase(root_children.begin());
		  XX->eles++;
		}
	      //no more childs fit into the root node. 
              //Give them back to the root
	    }
	  put();
	}
	count++;
      }
    get(act_level,act_pthis);
    delete [] IndexPthis;
  }



  // a inner node is split.
   void split_node_star(std::vector<v_type> &children,
		       std::vector<v_type> &root_children, 
		       std::vector<bool>& Reinsert,
		       int MaxLevel,
		       std::pair< int,long> *LevelPthis) {
     int act_level=XX->level; 
     if (XX->eles >= IO_max_cap_nodes) 
       {
	 //it makes no sense to reinsert in the root level
	 CGAL_DOGN_Insert(
		     std::cerr << "(((((((((REINSERT)))))))))))" 
		     << act_level << std::endl;
		     )
	   if(!Reinsert[act_level] || act_level==MaxLevel)
	     {
	       split_star(children,*(dynamic_cast<Node_entry*>(XX)));
	     }
	   else
	     {
	       CGAL_DOGN_Insert(
			   std::cerr << "(((((((((REINSERT-REALLY)))))))))))" 
			   << act_level << std::endl;
			   )
		 Reinsert[act_level]=false;
		 reinsert_node(root_children, act_level,Reinsert,MaxLevel,
			       LevelPthis);
	     }
       }
  }
  //end R_star

  void direct_insert(Data &d, Data_entry& YY)  {
    *ff=YY.raw_data + YY.leaf_offset;
    traits.write(ff,d);
    YY.leaf_offset += traits.written_size(d);
    YY.eles++;
    CGAL_DOGN_ControlB(
		  std::cerr << "direct insert: number of elements" 
		  << YY.eles << std::endl;
		  )
 }


  // the element is inserted into a leaf node.
  // If the available space is not large enough, the
  // data of the leaf node is split according to the split
  // strategy
   bool must_be_reinserted(Data &d, 
	std::vector<v_type>& children, std::vector<v_type>& root_children,
	Data_entry& YY,std::vector<bool>& Reinsert, int MaxLevel,
        std::pair< int,long> *LevelPthis)  {

    CGAL_DOGN_Insert(
		std::cerr << "Must Be Reinserted\n";
		dump((* dynamic_cast<Data_entry*>(XX)),0,0);
		)
    if (( int)YY.available_space(YY.leaf_offset)< (int)traits.written_size(d))
      {
	CGAL_DOGN_Insert(
	            std::cerr << "available_space :" << 
		        YY.available_space(YY.leaf_offset) << std::endl;
		    std::cerr << "size of data:" << traits.written_size(d) 
		    << std::endl;
		    std::cerr << "pagesize:" << IO_page_size << std::endl;
		    std::cerr << "(((((((((REINSERT))))))))))) 0 " 
		    << std::endl;
		    )
	if(!Reinsert[0] || MaxLevel==0 )    //no reinsertion but a split
	  {
	    std::vector<Data>  V;
	    std::vector<d_type> left, right;
	    d_type *old;
	    old= new d_type[YY.eles+2];
	    Data_entry *nleft=new Data_entry();
	    Data_entry *nright=new Data_entry();
	    Data          element;
	    d_type   val;
	    v_type vleft, vright;
	    val.key = traits.build(d);
	    val.ele = d;
	    old[0]=val;
	     int i, offset =0;
	    for(i=0;i<YY.eles; i++)
	      {
		*ff=YY.raw_data + offset;
		traits.read(ff, element);
		V.push_back(element);
		old[i+1].key = traits.build(element);
		old[i+1].ele = element;
		offset += traits.written_size(element);
	      }
	    
	    ps_leaf(&old[0], &old[YY.eles+1],  std::back_inserter(left), 
		    std::back_inserter(right), IO_min_cap_leaves,
		    IO_page_size);
	    delete[] old;
	    
	    d_type *dt;
	    Key kleft =left.begin()->key;
	    Key kright =right.begin()->key;
	    for(dt=left.begin();dt!=left.end();dt++)
	      {
		direct_insert(dt->ele,*nleft);
		kleft = traits.unify(kleft, dt->key);
	      }
	    for(dt=right.begin();dt!=right.end();dt++)
	      {
		direct_insert(dt->ele, *nright);
		kright = traits.unify(kright, dt->key);
	      }
	    
	    CGAL_DOGN_Insert(
			std::cerr << "SPLIT the left and right keys are: " 
			<< std::endl;
			traits.dump(kleft);
			traits.dump(kright);
			std::cerr << std::endl;
			)
	    
	    vleft.key = kleft;  
	    put_data(*nleft, YY.pthis);
	    vleft.pnext = YY.pthis;
	    vleft.deleted = false;
	    
	    vright.key = kright;
	    put_data(*nright,-1);
	    vright.pnext = nright->pthis;
	    vright.deleted = false;

	    children.erase(children.begin(), children.end());
	    children.push_back(vleft);
	    children.push_back(vright);
	    get(0, YY.pthis);
	    delete nleft;
	    delete nright;
	    return false;
	  }
	else{ // do reinsertion - return true
	  CGAL_DOGN_Insert( 
		      std::cerr << "(((((((((REINSERT))))))))))) of 0 level" 
		      << std::endl;
		      )
	    return true;
	}
      }
    else
      {
	*ff=YY.raw_data + YY.leaf_offset;
	traits.write(ff,d);
	YY.leaf_offset += traits.written_size(d);
	YY.eles++;
	put();
	return false;
      }
    return true;
  }
  //******************* end R_star

  // the element is inserted into a leaf node.
  // If the available space is not large enough, the
  // data of the leaf node is split according to the split
  // strategy
   bool insert_star(Data &d, 
		   std::vector<v_type>& children, 
		   std::vector<v_type>& root_children,
		   std::vector<bool>& Reinsert,  int MaxLevel,
		   std::pair< int,long> *LevelPthis)  {
     int act_level=XX->level; 
    long act_pthis=XX->pthis;
    if( must_be_reinserted(d,children, root_children,
			   *(dynamic_cast<Data_entry*>(XX)),Reinsert, 
			   MaxLevel, LevelPthis))
      {
	CGAL_DOGN_Insert(
		    std::cerr << "(((((((((((REINSERT REALLY)))))))))" 
		    << act_level << std::endl;
		    )
	get(act_level,act_pthis);
	Reinsert[0]=false;
	reinsert_data(d,root_children,0,Reinsert, MaxLevel,LevelPthis);
      }
    return true;
  }

  
  // the data is inserted into an inner node. The node is split
  // if necessary
   bool insert_star(v_type &v, Node_entry& YY)  {
    if ((YY.eles >= 0) && (YY.eles < IO_max_cap_nodes)) {
      YY.ele[YY.eles] = v;
      YY.eles++;
      put();
      return true;
    }
    return false;
  }
protected:
  
 // a node YY is split into two.
  void split_star(std::vector<v_type> &children, Node_entry &YY) {
    v_type vleft, vright;
    long act_pthis=YY.pthis;
     int act_level=YY.level;
    Node_entry *XXLeft=new Node_entry(YY.level);
    Node_entry *XXRight=new Node_entry(YY.level);
    ps(&YY.ele[0], &YY.ele[YY.eles], &XXLeft->ele[0], &XXRight->ele[0] , 
       IO_min_cap_nodes, IO_max_cap_nodes );
    
    vleft.key = build(*XXLeft);
    XXLeft->pthis=act_pthis;
    CGAL_DOGN_Insert(
		dump_this_node((*dynamic_cast<Node_entry*>(XX)),0,0);
		dump_this_node(*XXLeft,0,0);
		)
    put_node(*XXLeft, act_pthis);
    vleft.pnext = act_pthis;
    vleft.deleted = false;
    
    vright.key = build(*XXRight);

    CGAL_DOGN_Insert(
		std::cerr << "SPLIT the left and right keys are: " 
		<< std::endl;
		traits.dump(vleft.key);
		traits.dump(vright.key);
		std::cerr << std::endl;
		)
    put_node(*XXRight,-1);
    vright.pnext = XXRight->pthis;
    vright.deleted = false;

    CGAL_DOGN_Insert(
		dump_this_node(*XXRight,0,0);
		dump_this_node(*XXLeft,0,0);
		get(XXRight->level,XXRight->pthis);
		dump_this_node((*dynamic_cast<Node_entry*>(XX)),0,0);
		get(XXLeft->level,XXLeft->pthis);
		dump_this_node((*dynamic_cast<Node_entry*>(XX)),0,0);
		)
    children.erase(children.begin(), children.end());
    children.push_back(vleft);
    children.push_back(vright);
    get(act_level,act_pthis);
    delete XXLeft;
    delete XXRight;
  }
  
public:
  // the key of a node YY is build. It is the union of all keys
   Key build(Node_entry &YY) {
    Key k;
     int i=0;
    while( i<IO_max_cap_nodes && YY.ele[i].deleted)
      i++;
    if(i<IO_max_cap_nodes)
      {
	YY.eles=1;
	k = YY.ele[i].key;
	i++;
	for ( ;i<IO_max_cap_nodes; i++) {
	  if (!YY.ele[i].deleted) {
	    YY.eles++;
	    k=traits.unify(k,YY.ele[i].key);
	  }
	}
      }
    return k;
  }

  // the key of an data entry is build
   Key build(Data_entry &YY) 
  {     
    Data element;
    Key key;
     int i, offset=0;
    if(YY.eles >0)
      {
	*ff=YY.raw_data + offset;
	traits.read(ff, element);
	key=traits.build(element);
	offset += traits.written_size(element);
      }
    for(i=1;i<YY.eles; i++)
      {
	*ff=YY.raw_data + offset;
	traits.read(ff,element);
	key = traits.unify(key,traits.build(element));
	offset += traits.written_size(element);
      }
    return key;
  }
  
  //the data entry is searched for a key a
   bool find_key_data(Key& a, Data_entry& YY)
  {
    Data element;
    Key key;
     int i, offset=0;
    for(i=0;i<YY.eles; i++)
      {
	*ff=YY.raw_data + offset;
	traits.read(ff,element);
	key = traits.build(element);
	if(traits.equal_key(key,a))
	  return true;
	offset += traits.written_size(element);
      }
    return false;
  }

  bool get_next_pos(Key& a, int& pos, long& pnext, Node_entry& YY)
  {
    bool found=false;
    while(!found && pos<YY.eles)
      {
	if (traits.include(YY.ele[pos].key,a)) //compare before
	  {
	    found=true;
	    pnext = YY.ele[pos].pnext;
	  }
	  pos++;
      }
    return found;
  }

  //the data entry is searched for a key a
   bool find_key_node(Key& a)
  {
     int act_level=XX->level;
    long act_pthis=XX->pthis;
     int pos=0;
    long pnext=0;    
    while(get_next_pos(a, pos,pnext,*(dynamic_cast<Node_entry*>(XX)) ))
      {
	get(act_level-1,pnext);
	if(find_key(a))
	  return true;
	get(act_level,act_pthis);    
      }
    return false;
  }

  // an entry is searched for a key a
   bool find_key(Key& a) {
    if(dynamic_cast<Data_entry*>(XX))
      return find_key_data(a, *(dynamic_cast<Data_entry*>(XX)));
    if(dynamic_cast<Node_entry*>(XX))
      return find_key_node(a);
  }
  
  //the number of times key a is in YY is returned.
   long search_data(Key& a, Data_entry& YY)
  {
    long k = 0;
    Data element;
    Key key;
     int i, offset=0;
    for(i=0;i<YY.eles; i++)
      {
	*ff=YY.raw_data + offset;
	traits.read(ff,element);
	key = traits.build(element);
	if(traits.equal_key(key,a))
	  k++;
	offset += traits.written_size(element);
      }
    return k;
  }
  
  //the number of times key a is in YY is returned.
   long search_node(Key& a)
  {
     int act_level=XX->level;
    long act_pthis=XX->pthis;
    long k=0;
    long pnext;
     int pos=0;
    while(get_next_pos(a, pos,pnext,*(dynamic_cast<Node_entry*>(XX)) ))
      {
	get(act_level-1,pnext);
	k+= search(a);
	get(act_level,act_pthis);
      }
    return k;
  }

  //the number of times key a is in a node is returned.
   long search(Key& a) { 
    if(dynamic_cast<Data_entry*>(XX))
      return search_data(a, *(dynamic_cast<Data_entry*>(XX)));
    if(dynamic_cast<Node_entry*>(XX))
      return search_node(a);
  }
  
  
  
  //the number of times key a is in included in XX is returned.
   bool find_key_include(Key& a) {
    if(dynamic_cast<Data_entry*>(XX))
      return find_key_include_data(a, *(dynamic_cast<Data_entry*>(XX)));
    if(dynamic_cast<Node_entry*>(XX))
      return find_key_include_node(a);
    return false;
  }
  
  // true is returned if a is included by an element
   bool find_key_include_data(Key& a, Data_entry& YY)
  {
    Data element;
    Key key;
     int i, offset=0;
    for(i=0;i<YY.eles; i++)
      {
	*ff=YY.raw_data + offset;
	traits.read(ff,element);
        key = traits.build(element);
	//returns true if key includes a
	if(traits.include(key,a))
	  return true;
	offset += traits.written_size(element);
      }
    return false;
  }
  
  bool get_next_pos_include(Key& a, int& pos,  int& pnext, Node_entry& YY)
  {
    bool found=false;
    while(!found && pos<YY.eles)
      {
	if (traits.include(YY.ele[pos].key,a)) 
	  {
	    found=true;
	    pnext = YY.ele[pos].pnext;
	  }
	  pos++;
      }
    return found;
  }

  //the data entry is searched for a key a
   bool find_key_include_node(Key& a)
  {
     int act_level=XX->level;
    long act_pthis=XX->pthis;
     int pos=0, pnext;
    while(get_next_pos_include(a,pos,pnext,*(dynamic_cast<Node_entry*>(XX)) ))
      {
	get(act_level-1,pnext);
	if(find_key_include(a))
	  return true;
	get(act_level,act_pthis);    
      }
    return false;
  }


  // true is returned if a is intersected by an element
   bool find_key_intersect_data(Key& a, Data_entry& YY)
  {
    Data element;
    Key key;
     int i, offset=0;

    for(i=0;i<YY.eles; i++)
      {
	*ff=YY.raw_data + offset;
	traits.read(ff,element);
	key = traits.build(element);
	if(traits.intersect(key,a))
	  return true;
	offset += traits.written_size(element);
      }
    return false;
  }
  
  bool get_next_pos_inter(Key& a,  int& pos,  int& pnext, Node_entry& YY)
  {
    bool found=false;
    while(!found && pos<YY.eles)
      {
	if (traits.intersect(YY.ele[pos].key,a)) 
	  {
	    found=true;
	    pnext = YY.ele[pos].pnext;
	  }
	  pos++;
      }
    return found;
  }

  // true is returned if a is intersected by an element
   bool find_key_intersect_node(Key& a)
  {
     int pos=0,pnext;
     int act_level=XX->level;
    long act_pthis=XX->pthis;

     while(get_next_pos_inter(a, pos,pnext,*(dynamic_cast<Node_entry*>(XX)) ))
       {
	 get(act_level-1,pnext);
	 if(find_key_intersect(a))
	   return true;
	 get(act_level,act_pthis);
       }
     return false;
  }
  
    // true is returned if a is intersected by an element
   bool find_key_intersect(Key& a) {
    if(dynamic_cast<Data_entry*>(XX))
      return find_key_intersect_data(a, *(dynamic_cast<Data_entry*>(XX)));
    if(dynamic_cast<Node_entry*>(XX))
      return find_key_intersect_node(a);
    return false;
  }

  //********************INTERSECTION ITERATOR **************************

   // the node infos are colleced in a stack. Therefore the iteration 
   // process can resume from the last node info on the stack
   void collect_node_infos(Node_entry& YY,
			  std::list<Node_info*>& node_stack, Key k){
     int i;
    Node_info *n;
    for (i = 0; i < YY.eles; i++) 
      {
	if (traits.intersect(YY.ele[i].key,k)) 
	  {
	    n=new Node_info;
	    n->level=YY.level-1; 
	    n->pnext=YY.ele[i].pnext;
	    node_stack.push_front(n);
	    CGAL_DOGN_Search(
			std::cerr << "element put on stack with key: ";
			traits.dump(YY.ele[i].key,4);
			std::cerr << std::endl;
			)
	  }
      }
  }

   //same as above. But for the end iterator everything is reversed
   void reverse_collect_node_infos(Node_entry& YY,
				  std::list<Node_info*>& node_stack, Key k){
     int i;
    Node_info *n;
    for (i = 0; i<YY.eles; i++) 
      {
	if (traits.intersect(YY.ele[i].key,k)) 
	  {
	    n=new Node_info;
	    n->level=YY.level-1; //???
	    n->pnext=YY.ele[i].pnext;
	    node_stack.push_back(n);
	    CGAL_DOGN_Search(
			std::cerr << "element put on stack with key: ";
			traits.dump(YY.ele[i].key,4);
			std::cerr << "level" << n->level << std::endl;
			std::cerr << "pnext" << n->pnext << std::endl;
			std::cerr << std::endl;
			)
	 }
      }
  }


   void test_actual_data(Data_entry &actual_data)
  {
     int element_offset=0;
     int next_element=0;
    Data element;
    Key key;
    while (next_element < actual_data.eles)
      {
	*ff=actual_data.raw_data + element_offset;
	traits.read(ff,element);
	key = traits.build(element);
	traits.dump(key);
	next_element++;
	element_offset+=traits.written_size(element);
      }
  }

   bool found_init_iterator_data(Data_entry& YY, 
				 std::list<Node_info*>& node_stack,
				 Data_entry &actual_data,  int& next_element, 
				 Key k, Data &next_data,  int& element_offset)
  {
    Data element;
    Key key;
    element_offset=0;
    next_element=0;
    bool found=false;
    while (next_element < YY.eles && !found)
      {
	*ff=YY.raw_data + element_offset;
	traits.read(ff,element);
	key = traits.build(element);
	CGAL_DOGN_Search(
		    std::cerr << "data init 1\n";
		    traits.dump(key,4);
		    )
	if(traits.intersect(key,k))
	{
	  next_data=element;
	  found = true;
	  for( int i=0;i<IO_page_size;i++)
	    actual_data.raw_data[i]=YY.raw_data[i];
	  actual_data.eles=YY.eles;
	  CGAL_DOGN_Search(
	  	  std::cerr<< "-------------Test actual data:\n";
	  	  test_actual_data(actual_data);
	  	  std::cerr << "------------end\n";
		  )
	}
	element_offset += traits.written_size(element);
	next_element++;
      }
    return found;
  }

   bool init_iterator_data(std::list<Node_info*>& node_stack,
			  Data_entry &actual_data,  int& next_element, Key k, 
			  Data &next_data,  int& element_offset)
  {
    if(!found_init_iterator_data((*dynamic_cast<Data_entry*>(XX)),
				 node_stack,actual_data, next_element, k,
				 next_data,element_offset))
      {
	Node_info *ni;
	if(!node_stack.empty())
	  {
	    ni=(*node_stack.begin());
	    node_stack.pop_front();
	    CGAL_DOGN_Search(
			std::cerr << "data init 5 -- node stack not empty\n";
			std::cerr << "00000000pop from stack \n";
			std::cerr << "stack length " << node_stack.size()
			<< std::endl;
			std::cerr << "level:" << ni->level << std::endl;
			std::cerr << "pnext:" << ni->pnext << std::endl;
			)
	    get(ni->level,ni->pnext);
	    delete ni;
	    return init_iterator(node_stack,actual_data,next_element,
				 k,next_data,element_offset);
	  }
	return false;
      }
    return true;
  }


   bool found_reverse_init_iterator_data(Data_entry& YY, 
					std::list<Node_info*>& node_stack,
					Data_entry &actual_data,
					 int& next_element, Key k, 
					Data &next_data,  int& element_offset)
  {
    Data element;
    Key key;
    element_offset=0;
    next_element=0;
    bool found=false;
    while (next_element < YY.eles)
      {
	*ff=YY.raw_data + element_offset;
	traits.read(ff,element);
	key = traits.build(element);
	if(traits.intersect(key,k))
	{
	  next_data=element;
	  found = true;
	  for( int i=0;i<IO_page_size;i++)
	    actual_data.raw_data[i]=YY.raw_data[i];
	  actual_data.eles=YY.eles;
	}
	element_offset += traits.written_size(element);
	next_element++;
      }
    return found;
  }

   bool reverse_init_iterator_data(std::list<Node_info*>& node_stack,
				  Data_entry &actual_data, 
				   int& next_element, Key k, 
				  Data &next_data,  int& element_offset)
  {  
    Node_info *ni;
    if(!found_reverse_init_iterator_data((*dynamic_cast<Data_entry*>(XX)),
					 node_stack,actual_data, next_element,
					 k, next_data,element_offset))
      {
	if(!node_stack.empty())
	  {
	    ni=*(node_stack.begin());
	    node_stack.pop_front();
	    get(ni->level,ni->pnext);
	    delete ni;
	    return reverse_init_iterator(node_stack,actual_data,next_element,
					 k,next_data,element_offset);
	  }
	return false;
      }
    return true;
  }


   bool init_iterator_node(std::list<Node_info*>& node_stack,
			  Data_entry &actual_data, int& next_element, Key k, 
			  Data &next_data,  int& element_offset)
  {
    Node_info *ni;
    collect_node_infos((*dynamic_cast<Node_entry*>(XX)) ,node_stack,k);
    if(!node_stack.empty())
      {
	ni=*(node_stack.begin());
	node_stack.pop_front();
	CGAL_DOGN_Search(
		    std::cerr << "init_iterator node 2\n";
		    std::cerr << "0000000000pop from stack\n";
		    std::cerr << "stack length " << node_stack.size()
		    << std::endl;
		    std::cerr << "level:" << ni->level << std::endl;
		    std::cerr << "pnext:" << ni->pnext << std::endl;
		    std::cerr << "init_iterator node 3\n";
		    )
	get(ni->level,ni->pnext); 
	delete ni;
	return init_iterator(node_stack,actual_data,next_element,
			       k,next_data,element_offset);
      }
    return false;
  }

   bool reverse_init_iterator_node(std::list<Node_info*>& node_stack,
			     Data_entry &actual_data, int& next_element, 
			     Key k, Data &next_data,  int& element_offset)
  {
    Node_info *ni;
    
    reverse_collect_node_infos((*dynamic_cast<Node_entry*>(XX)),node_stack,k);
    if(!node_stack.empty())
      {
	ni=*(node_stack.begin());
	node_stack.pop_front();
	get(ni->level,ni->pnext);
	delete ni;
	return reverse_init_iterator(node_stack,actual_data,next_element,
				     k,next_data,element_offset);
      }
    return false;
  }


   bool init_iterator(std::list<Node_info*>& node_stack,
		     Data_entry &actual_data, int& next_element, 
		     Key k, Data &next_data,  int& element_offset)
  {
    if(dynamic_cast<Data_entry*>(XX))
      return init_iterator_data(node_stack,
			  actual_data,next_element,k,next_data,element_offset);
    if(dynamic_cast<Node_entry*>(XX))
      return init_iterator_node(node_stack,
			  actual_data,next_element,k,next_data,element_offset);
    return false;
  }

   bool reverse_init_iterator(std::list<Node_info*>& node_stack,
			     Data_entry& actual_data, int& next_element, 
			     Key k, Data& next_data,  int& element_offset)
  {
    if(dynamic_cast<Data_entry*>(XX))
      return reverse_init_iterator_data(node_stack,
				   actual_data,next_element,k,next_data,
				   element_offset);
     if(dynamic_cast<Node_entry*>(XX))
      return reverse_init_iterator_node(node_stack,
				   actual_data,next_element,k,next_data,
				   element_offset);
    return false;
  }


 // find the next element of the iteration
 bool next_element(std::list<Node_info*>& node_stack,
		    Data_entry& actual_data, int& next_element, 
		    Key k, Data& next_data,  int& element_offset)     
  {
    CGAL_DOGN_Search(
		std::cerr << "next element 1\n";
		std::cerr << next_element << std::endl;
		std::cerr << element_offset << std::endl;
		std::cerr<< "-------------Test actual data:\n";
		test_actual_data(actual_data);
		std::cerr << "------------end\n";
		)
    Node_info *ni;
    Data element;
    Key key;
    bool found=false;
    while (next_element < actual_data.eles && !found)
      {
	*ff=actual_data.raw_data + element_offset;
	traits.read(ff,element);
	key = traits.build(element);
	CGAL_DOGN_Search(
		    std::cerr << "next element 2\n";
		    traits.dump(key,4);
		    )
	if(traits.intersect(key,k))
	{
	  next_data=element;
	  found = true;
	}
	element_offset += traits.written_size(element);
	next_element++;
      }
    if(!found)
      {
	CGAL_DOGN_Search(  
		    std::cerr << "next element 4\n";
		    )
	if(!node_stack.empty())
	  {
	    ni=*(node_stack.begin());
	     
	    int zz=ni->level;
	    CGAL_DOGN_Search(
			 int z=ni->level;
			std::cerr << "next element 5\n";
			std::cerr << "stack length " << node_stack.size()
			<< std::endl;
			std::cerr << "0000000000pop from stack\n";
			std::cerr << "level:" << ni->level << std::endl;
			std::cerr << "pnext:" << ni->pnext << std::endl;
			std::cerr << "level:z " << z << std::endl;
			std::cerr << "level:" << ni->level << std::endl;
			std::cerr << "level:" << ni->level << std::endl;
			z=ni->level;
			)
	    get(zz,ni->pnext);
	    node_stack.pop_front();
	    delete ni;
	    return init_iterator(node_stack,actual_data,next_element,
				 k,next_data,element_offset);
	  }
	return false;
      }
    return true;
  }

  //********************ENCLOSING ITERATOR **************************
 //everything for the enclosing iterator

   void collect_encl_node_infos(Node_entry& YY,
			       std::list<Node_info*>& node_stack, Key k){
    iterator i;
    Node_info *n;
    for (i = &YY.ele[0]; i !=&YY.ele[YY.eles];i++) 
      {
	if (traits.include((*i).key,k)) 
	  {
	    n=new Node_info;
	    n->level=YY.level-1; // or like that???
	    n->pnext=(*i).pnext;
	    node_stack.push_front(n);
	    CGAL_DOGN_Search(
			std::cerr << "element put on stack with key: ";
			traits.dump((*i).key,4);
			std::cerr << std::endl;
			)
	  }
      }
  }

   void reverse_collect_encl_node_infos(Node_entry& YY,
				  std::list<Node_info*>& node_stack, Key k){
    iterator i;
    Node_info *n;	    
    for (i = &YY.ele[0]; i !=&YY.ele[YY.eles];i++) 
      {
	if (traits.intersect((*i).key,k)) 
	  {
	    n=new Node_info;
	    n->level=YY.level-1; //or like that ???
	    n->pnext=(*i).pnext;
	    node_stack.push_back(n);
	    CGAL_DOGN_Search(
			std::cerr << "element put on stack with key: ";
			traits.dump((*i).key,4);
			std::cerr << std::endl;
			)
	  }
      }
  }

   bool found_init_encl_iterator_data(Data_entry& YY, std::list<Node_info*>& 
				      node_stack, Data_entry &actual_data, 
				      int& next_element, Key k, 
			              Data &next_data,  int& element_offset)
  {
    Data element;
    Key key;
    element_offset=0;
    next_element=0;
    bool found=false;
    while (next_element < YY.eles && !found)
      {
	*ff=YY.raw_data + element_offset;
	traits.read(ff,element);
	key = traits.build(element);
	CGAL_DOGN_Search(
		    traits.dump(key,4);
		    )
	if(traits.include(key,k))
	{
	  next_data=element;
	  found = true;
	  for( int i=0;i<IO_page_size;i++)
	    actual_data.raw_data[i]=YY.raw_data[i];
	  actual_data.eles=YY.eles;
	}
	element_offset += traits.written_size(element);
	next_element++;
      }
    return found;
  }


   bool init_encl_iterator_data(std::list<Node_info*>& node_stack,
			       Data_entry &actual_data, int& next_element, 
			       Key k, Data &next_data,  int& element_offset)
  {
    Node_info *ni;
    if(!found_init_encl_iterator_data((*dynamic_cast<Data_entry*>(XX)),
				       node_stack,actual_data, next_element, k,
				      next_data,element_offset))
      {
	if(!node_stack.empty())
	  {
	    ni=(*node_stack.begin());
	    node_stack.pop_front();
	    get(ni->level,ni->pnext);
	    delete ni;
	    return init_encl_iterator(node_stack,actual_data,next_element,
				      k,next_data,element_offset);
	  }
	return false;
      }
    return true;
  }



   bool found_reverse_init_encl_iterator_data(Data_entry& YY, 
					     std::list<Node_info*>& node_stack,
					     Data_entry &actual_data,
					      int& next_element, Key k, 
					     Data &next_data, 
					      int& element_offset)
    {
      Data element;
      Key key;
      element_offset=0;
      next_element=0;
      bool found=false;
      while (next_element < YY.eles)
	{
	  *ff=YY.raw_data + element_offset;
	  traits.read(ff,element);
	  key = traits.build(element);
	  if(traits.include(key,k))
	    {
	      next_data=element;
	      found = true;
	      for( int i=0;i<IO_page_size;i++)
		actual_data.raw_data[i]=YY.raw_data[i];
	      actual_data.eles=YY.eles;
	    }
	  element_offset += traits.written_size(element);
	  next_element++;
	}
      return found;
    }

   bool reverse_init_encl_iterator_data(std::list<Node_info*>& node_stack,
				       Data_entry &actual_data,
				        int& next_element, Key k, 
				       Data &next_data,  int& element_offset)
    {
      Node_info *ni;
      if(!found_reverse_init_encl_iterator_data(
                      (*dynamic_cast<Data_entry*>(XX)),
		      node_stack,actual_data, next_element, k,
		      next_data,element_offset))
	{
	  if(!node_stack.empty())
	    {
	      ni=*(node_stack.begin());
	      node_stack.pop_front();
	      get(ni->level,ni->pnext);
	      delete ni;
	      return reverse_init_encl_iterator(node_stack,actual_data,
						next_element,
						k,next_data,element_offset);
	    }
	  return false;
	}
      return true;
    }


   bool init_encl_iterator_node(std::list<Node_info*>& node_stack,
			       Data_entry &actual_data,
			       int& next_element, Key k, 
			       Data &next_data,  int& element_offset)
  {
    Node_info *ni;
    //    std::cerr << "init_iterator node 1\n";
    collect_encl_node_infos((*dynamic_cast<Node_entry*>(XX)),node_stack,k);
    if(!node_stack.empty())
      {
	ni=*(node_stack.begin());
	node_stack.pop_front();
	get(ni->level,ni->pnext);
	delete ni;
	return init_encl_iterator(node_stack,actual_data,next_element,
				  k,next_data,element_offset);
      }
    return false;
  }


   bool reverse_init_encl_iterator_node(std::list<Node_info*>& node_stack,
				       Data_entry &actual_data,
				        int& next_element, Key k, 
				       Data &next_data,  int& element_offset)
  {
    Node_info *ni;
    
    reverse_collect_node_infos((*dynamic_cast<Node_entry*>(XX)),node_stack,k);
    if(!node_stack.empty())
      {
	ni=*(node_stack.begin());
	node_stack.pop_front();
	get(ni->level,ni->pnext);
	delete ni;
	return reverse_init_encl_iterator(node_stack,actual_data,next_element,
					  k,next_data,element_offset);
      }
    return false;
  }


   bool init_encl_iterator(std::list<Node_info*>& node_stack,
			  Data_entry &actual_data, int& next_element, 
			  Key k, Data &next_data,  int& element_offset)
  {
    if(dynamic_cast<Data_entry*>(XX))
      return init_encl_iterator_data(node_stack,
				actual_data,next_element,k,next_data,
				element_offset);
    if(dynamic_cast<Node_entry*>(XX))
      return init_encl_iterator_node(node_stack,
				actual_data,next_element,k,next_data,
				element_offset);
    return false;
  }

   bool reverse_init_encl_iterator(std::list<Node_info*>& node_stack,
				  Data_entry& actual_data, int& next_element, 
				  Key k, Data& next_data,  int& element_offset)
  {
    if(dynamic_cast<Data_entry*>(XX))
      return reverse_init_encl_iterator_data(node_stack,
					actual_data,next_element,k,next_data,
					element_offset);
    if(dynamic_cast<Node_entry*>(XX))
      return reverse_init_encl_iterator_node(node_stack,
					actual_data,next_element,k,next_data,
					element_offset);
    return false;
  }

   bool next_encl_element(std::list<Node_info*>& node_stack,
			 Data_entry& actual_data, int& next_element, 
			 Key k, Data& next_data,  int& element_offset)     
  {
    Node_info *ni;
    Data element;
    Key key;
    bool found=false;
    while (next_element < actual_data.eles && !found)
      {
	*ff=actual_data.raw_data + element_offset;
	traits.read(ff,element);
	key = traits.build(element);
	CGAL_DOGN_Search(
		    traits.dump(key,4);
		    )
	if(traits.include(key,k))
	{
	  next_data=element;
	  found = true;
	}
	element_offset += traits.written_size(element);
	next_element++;
      }
    if(!found)
      {
	if(!node_stack.empty())
	  {
	    ni=*(node_stack.begin());
	     int  zz=ni->level;
	    get(zz,ni->pnext);
	    node_stack.pop_front();
	    delete ni;
	    return init_encl_iterator(node_stack,actual_data,next_element,
					k,next_data,element_offset);
	  }
	return false;
      }
    return true;
  }


  //********************COMPARE ITERATOR **************************


   void collect_compare_node_infos(Node_entry& YY,
			       std::list<Node_info*>& node_stack, Key k){
    iterator i;
    Node_info *n;
    for (i = &YY.ele[0]; i !=&YY.ele[YY.eles];i++) 
      {
	if (traits.compare((*i).key,k)) 
	  {
	    n=new Node_info;
	    n->level=YY.level-1; // or like that???
	    n->pnext=(*i).pnext;
	    node_stack.push_front(n);
	    CGAL_DOGN_Search(
			std::cerr << "element put on stack with key: ";
			traits.dump((*i).key,4);
			std::cerr << std::endl;
			)
	  }
      }
  }

   void reverse_collect_compare_node_infos(Node_entry& YY,
				  std::list<Node_info*>& node_stack, Key k){
    iterator i;
    Node_info *n;
    for (i = &YY.ele[0]; i !=&YY.ele[YY.eles];i++) 
      {
	if (traits.compare((*i).key,k)) 
	  {
	    n=new Node_info;
	    n->level=YY.level-1; //or like that ???
	    n->pnext=(*i).pnext;
	    node_stack.push_back(n);
	    CGAL_DOGN_Search(
			std::cerr << "element put on stack with key: ";
			traits.dump((*i).key,4);
			std::cerr << std::endl;
			)
	  }
      }
  }

   bool found_init_compare_iterator_data(Data_entry& YY, 
			       std::list<Node_info*>& node_stack,
			       Data_entry &actual_data, int& next_element, 
			       Key k, Data &next_data,  int& element_offset)
  {
    Data element;
    Key key;
    element_offset=0;
    next_element=0;
    bool found=false;
    while (next_element < YY.eles && !found)
      {
	*ff=YY.raw_data + element_offset;
	traits.read(ff,element);
	key = traits.build(element);
	CGAL_DOGN_Search(
		    traits.dump(key,4);
		    )
	if(traits.compare(key,k))
	{
	  next_data=element;
	  found = true;
	  for( int i=0;i<IO_page_size;i++)
	    actual_data.raw_data[i]=YY.raw_data[i];
	  actual_data.eles=YY.eles;
	}
	element_offset += traits.written_size(element);
	next_element++;
      }
    return found;
  }


   bool init_compare_iterator_data(std::list<Node_info*>& node_stack,
			       Data_entry &actual_data, int& next_element, 
			       Key k, Data &next_data,  int& element_offset)
  {
    Node_info *ni;
    if(!found_init_compare_iterator_data((*dynamic_cast<Data_entry*>(XX)),
					 node_stack,actual_data, 
					 next_element, k,
				         next_data,element_offset))
      {
	if(!node_stack.empty())
	  {
	    ni=(*node_stack.begin());
	    node_stack.pop_front();
	    get(ni->level,ni->pnext);
	    delete ni;
	    return init_compare_iterator(node_stack,actual_data,next_element,
				      k,next_data,element_offset);
	  }
	return false;
      }
    return true;
  }



   bool found_reverse_init_compare_iterator_data(Data_entry& YY, 
					     std::list<Node_info*>& node_stack,
					     Data_entry &actual_data,
					      int& next_element, Key k, 
					     Data &next_data, 
					      int& element_offset)
    {
      Data element;
      Key key;
      element_offset=0;
      next_element=0;
      bool found=false;
      while (next_element < YY.eles)
	{
	  *ff=YY.raw_data + element_offset;
	  traits.read(ff,element);
	  key = traits.build(element);
	  if(traits.compare(key,k))
	    {
	      next_data=element;
	      found = true;
	      for( int i=0;i<IO_page_size;i++)
		actual_data.raw_data[i]=YY.raw_data[i];
	      actual_data.eles=YY.eles;
	    }
	  element_offset += traits.written_size(element);
	  next_element++;
	}
      return found;
    }

   bool reverse_init_compare_iterator_data(std::list<Node_info*>& node_stack,
				       Data_entry &actual_data,
				        int& next_element, Key k, 
				       Data &next_data,  int& element_offset)
    {
      Node_info *ni;
      if(!found_reverse_init_compare_iterator_data(
                               (*dynamic_cast<Data_entry*>(XX)),
				node_stack,actual_data, next_element, k,
				next_data,element_offset))
	{
	  if(!node_stack.empty())
	    {
	      ni=*(node_stack.begin());
	      node_stack.pop_front();
	      get(ni->level,ni->pnext);
	      delete ni;
	      return reverse_init_compare_iterator(
	                        node_stack,actual_data,next_element,
				k,next_data,element_offset);
	    }
	  return false;
	}
      return true;
    }


   bool init_compare_iterator_node(std::list<Node_info*>& node_stack,
			       Data_entry &actual_data,
			       int& next_element, Key k, 
			       Data &next_data,  int& element_offset)
  {
    Node_info *ni;
    //    std::cerr << "init_iterator node 1\n";
    collect_compare_node_infos((*dynamic_cast<Node_entry*>(XX)),node_stack,k);
    if(!node_stack.empty())
      {
	ni=*(node_stack.begin());
	node_stack.pop_front();
	get(ni->level,ni->pnext);
	delete ni;
	return init_compare_iterator(node_stack,actual_data,next_element,
				       k,next_data,element_offset);
      }
    return false;
  }


   bool reverse_init_compare_iterator_node(std::list<Node_info*>& node_stack,
				       Data_entry &actual_data,
				        int& next_element, Key k, 
				       Data &next_data,  int& element_offset)
  {
    Node_info *ni;
    
    reverse_collect_node_infos((*dynamic_cast<Node_entry*>(XX)),node_stack,k);
    if(!node_stack.empty())
      {
	ni=*(node_stack.begin());
	node_stack.pop_front();
	get(ni->level,ni->pnext);
	delete ni;
	return reverse_init_compare_iterator(node_stack,
                                             actual_data,next_element,
					     k,next_data,element_offset);
      }
    return false;
  }


   bool init_compare_iterator(std::list<Node_info*>& node_stack,
			  Data_entry &actual_data, int& next_element, 
			  Key k, Data &next_data,  int& element_offset)
  {
    if(dynamic_cast<Data_entry*>(XX))
      return init_compare_iterator_data(node_stack,
				actual_data,next_element,k,next_data,
				element_offset);
    if(dynamic_cast<Node_entry*>(XX))
      return init_compare_iterator_node(node_stack,
				actual_data,next_element,k,next_data,
				element_offset);
    return false;
  }

   bool reverse_init_compare_iterator(std::list<Node_info*>& node_stack,
				  Data_entry& actual_data, int& next_element, 
				  Key k, Data& next_data,  int& element_offset)
  {
    if(dynamic_cast<Data_entry*>(XX))
      return reverse_init_compare_iterator_data(node_stack,
					actual_data,next_element,k,next_data,
					element_offset);
    if(dynamic_cast<Node_entry*>(XX))
      return reverse_init_compare_iterator_node(node_stack,
					actual_data,next_element,k,next_data,
					element_offset);
    return false;
  }

   bool next_compare_element(std::list<Node_info*>& node_stack,
			 Data_entry& actual_data, int& next_element, 
			 Key k, Data& next_data,  int& element_offset)     
  {
    Node_info *ni;
    Data element;
    Key key;
    bool found=false;
    while (next_element < actual_data.eles && !found)
      {
	*ff=actual_data.raw_data + element_offset;
	traits.read(ff,element);
	key = traits.build(element);
	CGAL_DOGN_Search(
		    traits.dump(key,4);
		    )
	if(traits.compare(key,k))
	{
	  next_data=element;
	  found = true;
	}
	element_offset += traits.written_size(element);
	next_element++;
      }
    if(!found)
      {
	if(!node_stack.empty())
	  {
	    ni=*(node_stack.begin());
	     int  zz=ni->level;
	    get(zz,ni->pnext);
	    node_stack.pop_front();
	    delete ni;
	    return init_compare_iterator(node_stack,actual_data,next_element,
					k,next_data,element_offset);
	  }
	return false;
      }
    return true;
  }


  //*************************SEARCH ***************************




  // the elements are returned that a intersects.
   std::back_insert_iterator<std::list<Data> > search_intersect(Key& a, 
		    Data_entry& YY,
		    std::back_insert_iterator<std::list<Data> >& result)
  {
    Data element;
    Key key;
     int i, offset=0;
    for(i=0;i<YY.eles; i++)
      {
	*ff=YY.raw_data + offset;
	traits.read(ff,element);
	key = traits.build(element);
	if(traits.intersect(key,a))
	  (*result++)=element;
	offset += traits.written_size(element);
      }
    return result;
  }

   Key bounding_box(Data_entry& YY)
  {
    Key key,k;
    Data element;
     int i, offset=0;
     *ff=YY.raw_data + offset;
    traits.read(ff,element);
    key = traits.build(element);
    for(i=0;i<YY.eles; i++)
    {
      *ff=YY.raw_data + offset;
      traits.read(ff,element);
      key = traits.build(element);
      k=traits.unify(k,key);
      offset += traits.written_size(element);
      CGAL_DOGN_ControlB(
		    traits.dump(k,4);
		    )
    }
    return k;
  }

   Key bounding_box(Node_entry& YY)
  {
     int i;
    Key k;
    k=YY.ele[0].key;
    for (i = 1; i < YY.eles; i++) 
    {
       k=traits.unify(k,YY.ele[i].key);
       CGAL_DOGN_ControlB(
		     traits.dump(k,4);
		     )
    }
    return k;
  }


   Key bounding_box()
  {
    if(dynamic_cast<Data_entry*>(XX))
      return bounding_box(*(dynamic_cast<Data_entry*>(XX)));
    if(dynamic_cast<Node_entry*>(XX))
      return bounding_box(*(dynamic_cast<Node_entry*>(XX)));
    Key a; //fake
    return a;
  }
    

  //the key at pos i was deleted. Another key is copied to that place
   void reorder_keys(Node_entry& YY,  int i)
  {
     int ki=YY.eles;
    while(ki>i && YY.ele[ki].deleted)
      ki--;
    if(ki>i)
      {
	YY.ele[i]=YY.ele[ki];
	YY.ele[ki].deleted=true;
      }
  }

   //after a reinsertion call the parent keys have to be adapted
   void adapt_parent_keys(std::pair< int,long>* LevelPthis, Key k,  int level,
			  int MaxLevel)
  {
    CGAL_DOGN_Insert(
		std::cerr << "ADAPT PARENT KEYS\n";
		if(dynamic_cast<Data_entry*>(XX))
		dump((* dynamic_cast<Data_entry*>(XX)),0,0);
		if(dynamic_cast<Node_entry*>(XX))
		dump_this_node((* dynamic_cast<Node_entry*>(XX)),0,0);
		)
    if(level<=MaxLevel)
      {
	get(level,LevelPthis[level].second);
	(*(dynamic_cast<Node_entry*>(XX))).ele[LevelPthis[level].first].key=k;
	put();
	k = build(*(dynamic_cast<Node_entry*>(XX)));
	adapt_parent_keys(LevelPthis,k,level+1,MaxLevel);
      }
  }

   void adapt_parent_keys(std::pair< int,long>* LevelPthis,  int level,
			  int MaxLevel,std::vector<v_type>& root_children,
			 std::vector<bool>& Reinsert)
  {
    CGAL_DOGN_Insert(
		std::cerr << "ADAPT PARENT KEYS\n";
		if(dynamic_cast<Data_entry*>(XX))
		dump((* dynamic_cast<Data_entry*>(XX)),0,0);
		if(dynamic_cast<Node_entry*>(XX))
		dump_this_node((* dynamic_cast<Node_entry*>(XX)),0,0);
		)
    if(level<=MaxLevel)
      {
	get(level,LevelPthis[level].second);
	(*(dynamic_cast<Node_entry*>(XX))).ele[LevelPthis[level].first].deleted=
	                                                                   true;
	XX->eles--;
	reorder_keys(*(dynamic_cast<Node_entry*>(XX)),LevelPthis[level].first);
	put();
	CGAL_DOGN_Insert(	
		    dump_this_node((* dynamic_cast<Node_entry*>(XX)),0,0);
		    )
	if(XX->eles < IO_min_cap_nodes && level < MaxLevel) 
	  //must be deleted
	  {
	    adapt_parent_keys(LevelPthis,level+1,MaxLevel,
	                      root_children,Reinsert);
	    // get the child node again and reinsert all element on level
	    get(level,LevelPthis[level].second);
	     int j;
	    for(j=0;j<XX->eles;j++)
	      {    
		CGAL_DOGN_Insert(
			    dump_this_node((* dynamic_cast<Node_entry*>(XX)),
			                   0,0);
			    )
		v_type v=((*dynamic_cast<Node_entry*>(XX))).ele[j];
		root_forced_reinsert(v,root_children, level,Reinsert,
				     MaxLevel,LevelPthis);
		get(level,LevelPthis[level].second);	
	      }
	    erase(LevelPthis[level].second);
	  }
	else
	  {
	    Key k;
    	    k = build(*(dynamic_cast<Node_entry*>(XX)));
	    CGAL_DOGN_Insert(
			traits.dump(k);    
			dump_this_node((* dynamic_cast<Node_entry*>(XX)),0,0);
			)
	    adapt_parent_keys(LevelPthis,k,level+1,MaxLevel);
	  }
      }
  }

  // the elements are deleted that have key a.
   bool XXdelete_key_data(Key& a,
			 Data& result,
			 std::pair< int,long> *LevelPthis,
			 std::vector<bool>& Reinsert,  int MaxLevel,
			 std::vector<v_type>& children)
  {
    bool found=false;
    CGAL_DOGN_Delete(
		std::cerr << "XXdelete_key_data\n";
		if(dynamic_cast<Data_entry*>(XX))
		dump((* dynamic_cast<Data_entry*>(XX)),0,0);
		if(dynamic_cast<Node_entry*>(XX))
		dump_this_node((* dynamic_cast<Node_entry*>(XX)),0,0);
		)
    Data element;
    Key key;
    bool touched = false;
     int i, num_eles = XX->eles;
     int offset=0, act_level=XX->level;
    long act_pthis=XX->pthis;
    std::list<Data> to_del;
    i=0;
    while(i<num_eles && !found)
      {
	*ff=((*(dynamic_cast<Data_entry*>(XX))).raw_data + offset);
	traits.read(ff,element);
	key = traits.build(element);
	if(traits.equal_key(key,a))
	  {
	    touched = true;
	    found = true;
	    result=element;
	    (*(dynamic_cast<Data_entry*>(XX))).delete_key(offset, 
					traits.written_size(element));
	    put();
	    CGAL_DOGN_Delete(
			dump((* dynamic_cast<Data_entry*>(XX)),0,0);
			)
	  }
	else
	  offset += traits.written_size(element);
	i++;
      }
    if(XX->underfilled() && MaxLevel > 0) //otherwise a node can be underfilled
      {
	// the vertex is eliminated, the parent keys are adapted
	adapt_parent_keys( LevelPthis,XX->level+1, MaxLevel, 
			   children, Reinsert); 
	get(act_level,act_pthis);	
	CGAL_DOGN_Delete(
		    dump((* dynamic_cast<Data_entry*>(XX)),0,0);
		    )
	Data element;
	std::vector<v_type> act_children;
	v_type v;
	 int j,offset=0;
	num_eles=XX->eles;
	for(j=0;j<num_eles;j++)
	  { 
	    *ff=(*(dynamic_cast<Data_entry*>(XX))).raw_data + offset;
	    traits.read(ff,element);
	    offset += traits.written_size(element); 
	    get(MaxLevel,LevelPthis[MaxLevel].second);
	    v.deleted = false;
	    v.key= traits.build(element);
	    //use now ordenary insert since it has not to be 
	    //inserted on a special level
	    insert_star(element, v, act_children,children,Reinsert, 
			MaxLevel,LevelPthis);
	    if(!act_children.empty()) {
	      children.insert(children.begin(),*act_children.begin());
	      act_children.erase(act_children.begin());
	    }
	    while(!act_children.empty()){
	      children.push_back(*act_children.begin());
	      act_children.erase(act_children.begin());
	    }
	    if(!children.empty()) {
	      CGAL_DOGN_Delete(
			  std::cerr << "Rootpthis is:" << 
			  LevelPthis[MaxLevel].second << std::endl;
			  )
	      get(MaxLevel,LevelPthis[MaxLevel].second);
	      if(dynamic_cast<Node_entry*>(XX))
		{
		  int pos=LevelPthis[MaxLevel].first; //the changed index
		  (*(dynamic_cast<Node_entry*>(XX))).ele[pos]=*
		                                             (children.begin());
		  children.erase(children.begin());
		  while(XX->eles < IO_max_cap_nodes && !children.empty())
		    {
		      (*(dynamic_cast<Node_entry*>(XX))).ele[XX->eles]=
			*(children.begin());
		      children.erase(children.begin());
		      XX->eles++;
		    }
		  //no more childs fit into the root node. 
		  //Give them back to the root
		}
	      put();
	    }
	    get(act_level,act_pthis);	
	  }
	erase();
      }
    else
      if(touched)
	{
	  Key k= build(*(dynamic_cast<Data_entry*>(XX)));
	  adapt_parent_keys(LevelPthis,k,XX->level+1,MaxLevel); 
	  // the vertex is eliminated, the parent keys are adapted
	  get(act_level,act_pthis);
	  put();
	  offset = 0;
	  CGAL_DOGN_Delete(
		      std::cerr << "Key touched. New elements are: ";
		      dump(*(dynamic_cast<Data_entry*>(XX)),1,0);
		      std::cerr << "\n\n";
		      )
	}
    return found;
  }
  

   bool delete_key(Key& a, 
	   Data & result, 
	   std::pair< int,long> *LevelPthis,
	   std::vector<bool>& Reinsert,  int MaxLevel, 
	   std::vector<v_type>& children)
  {
    bool found=false;
    CGAL_DOGN_Delete(
		std::cerr << "Delete Key\n";
		if(dynamic_cast<Data_entry*>(XX))
		dump((* dynamic_cast<Data_entry*>(XX)),0,0);
		if(dynamic_cast<Node_entry*>(XX))
		dump_this_node((* dynamic_cast<Node_entry*>(XX)),0,0);
		)
     int i, act_level=XX->level;
    long act_pthis=XX->pthis;
    long child_pthis;
    for (i = 0; i< XX->eles; i++) 
      {
	if(!(*(dynamic_cast<Node_entry*>(XX))).ele[i].deleted)
	  {
	   // true if key includes a
	   if (traits.include((*(dynamic_cast<Node_entry*>(XX))).ele[i].key,a))
	      {
		std::pair< int,long> p(i,act_pthis);
		LevelPthis[act_level]=p;
		child_pthis=(*(dynamic_cast<Node_entry*>(XX))).ele[i].pnext;
		get(act_level-1,child_pthis);
		if(dynamic_cast<Node_entry*>(XX))
		  found = delete_key(a, result, LevelPthis ,Reinsert, 
				     MaxLevel, children);
		if(dynamic_cast<Data_entry*>(XX))
		  found = XXdelete_key_data(a, result, LevelPthis, Reinsert, 
					    MaxLevel, children);
		if (!get(act_level,act_pthis))
		  {
		    CGAL_DOGN_Control(
	std::cerr << "should not happen - node was deleted, everything ok\n";
				 )
		      return false;
		  }
		if(found)
		  return true;
	      }
	  }
      }
    return found;
  }

  
  // the nodes YY and ZZ are tried to put together.
   void add_or_split(std::vector<v_type>& children, Node_entry& YY, 
		    Node_entry& ZZ, std::vector<bool>& Reinsert)
  {
    int i; 
    long act_pthis=YY.pthis;
    int act_level=YY.level;
    if(YY.eles + ZZ.eles < IO_max_cap_nodes)
      {
	for(i=0;i<IO_max_cap_nodes ;i++)
	  {
	    if(!ZZ.ele[i].deleted)
	      insert_star(ZZ.ele[i], YY);
	  }
      }
    else
      {
	v_type vleft, vright; 
	Node_entry *XXLeft=new Node_entry(YY.level);
	Node_entry *XXRight=new Node_entry(YY.level);
	v_type keys[2*IO_max_cap_nodes];   //[YY.eles+ZZ.eles];
	 int j,k=0;
	for (j=0; j < IO_max_cap_nodes; ++j) 
	{
	  if (!YY.ele[j].deleted)
	    keys[k++]=YY.ele[j];
	}
	for (j=0; j < IO_max_cap_nodes; ++j) 
	{
	  if (!ZZ.ele[j].deleted)
	    keys[k++]=ZZ.ele[j];
	}
	ps(keys, &keys[YY.eles+ZZ.eles], &XXLeft->ele[0], &XXRight->ele[0],
	   IO_min_cap_nodes, IO_max_cap_nodes);   
	vleft.key = build(*XXLeft);
	XXLeft->pthis=act_pthis; 
	put_node(*XXLeft, act_pthis);
	vleft.pnext = act_pthis;
	vleft.deleted = false;

	CGAL_DOGN_Insert(
		    std::cerr << "SPLIT the left and right keys are: " 
		    << std::endl;
		    traits.dump(vleft.key);
		    traits.dump(vright.key);
		    std::cerr << std::endl;
		    )
	put_node(*XXRight,-1);
	vright.pnext = XXRight->pthis;
	vright.deleted = false;
	
	children.erase(children.begin(), children.end());
	children.push_back(vleft);
	children.push_back(vright);

	CGAL_DOGN_Insert(
		    dump_this_node(*XXRight,0,0);
		    dump_this_node(*XXLeft,0,0);
		    get(XXRight->level,XXRight->pthis);
		    dump_this_node((*dynamic_cast<Node_entry*>(XX)),0,0);
		    get(XXLeft->level,XXLeft->pthis);
		    dump_this_node((*dynamic_cast<Node_entry*>(XX)),0,0);
		    )
	get(act_level,act_pthis);
      }
  }
 
  // the data nodes YY and ZZ are tried to put together
   void add_or_split(std::vector<v_type>& children, Data_entry& YY, 
		     Data_entry& ZZ, 
		     std::vector<bool>& Reinsert)
  {
    if (YY.available_space(YY.leaf_offset)<ZZ.leaf_offset)
      {
	CGAL_DOGN_Insert(
		    std::cerr << "available_space :" << 
		    YY.available_space(YY.leaf_offset) << std::endl;
		    std::cerr << "size of data: ZZ style" << ZZ.leaf_offset 
		    << std::endl;
		    std::cerr << "pagesize:" << IO_page_size << std::endl;
		    )
	std::vector<Data>  V;
	std::vector<d_type> left, right;
	d_type *old;
	old= new d_type[YY.eles+YY.eles +2];
	Data_entry *nleft=new Data_entry();
	Data_entry *nright=new Data_entry();
	Data          element;
	v_type vleft, vright;
	 int i,j, offset =0;
	for(i=0;i<YY.eles; i++)
	  {
	    *ff=YY.raw_data + offset;
	    traits.read(ff,element);
	    V.push_back(element);
	    old[i].key = traits.build(element);
	    old[i].ele = element;
	    offset += traits.written_size(element);
	  }
	offset=0;
	for(j=0;j<ZZ.eles; i++)
	  {
	    *ff=ZZ.raw_data + offset;
	    traits.read(ff,element);
	    V.push_back(element);
	    old[j+i].key = traits.build(element);
	    old[j+i].ele = element;
	    offset += traits.written_size(element);
	  }
	ps_leaf(&old[0], &old[j+i], std::back_inserter(left), 
		std::back_inserter(right), IO_min_cap_leaves, 
		IO_page_size);
	delete[] old;
	d_type *dt;
	Key kleft =left.begin()->key;
	Key kright =right.begin()->key;
	for(dt=left.begin();dt!=left.end();dt++)
	  {
	    direct_insert(dt->ele,*nleft);
	    kleft = traits.unify(kleft, dt->key);
	  }
	for(dt=right.begin();dt!=right.end();dt++)
	  {
	    direct_insert(dt->ele, *nright);
	    kright = traits.unify(kright, dt->key);
	  }
	CGAL_DOGN_Insert(
		    std::cerr << "SPLIT the left and right keys are: " 
		    << std::endl;
		    traits.dump(kleft);
		    traits.dump(kright);
		    std::cerr << std::endl;
		    )

	vleft.key = kleft;
	put_data(*nleft, YY.pthis);
	vleft.pnext = YY.pthis;
	vleft.deleted = false;
	
	vright.key = kright;
	put_data(*nright,-1);
	vright.pnext = nright->pthis;
	vright.deleted = false;
	
	children.erase(children.begin(), children.end());
	children.push_back(vleft);
	children.push_back(vright);
	delete nleft;
	delete nright;
	get(0, YY.pthis);
      }
    else
      {
	Data  element;
	 int i, offset =0;
	for(i=0;i<ZZ.eles; i++)
	  {
	    *ff=ZZ.raw_data + offset;
	    traits.read(ff,element);
	    *ff=YY.raw_data + YY.leaf_offset;
	    traits.write(ff,element);
	    YY.leaf_offset += traits.written_size(element);
	    YY.eles++;
	    offset += traits.written_size(element);
	  }
	put();
      }
  }


  // tree information is given out
   void dump(Data_entry & YY,  int levelX,  int indent)
  {
    std::cerr << "********************************************\n";
    std::cerr << " dump of data\n";
    Data element;
    Key key;
     int i, offset=0;
    for(i=0;i<YY.eles; i++)
      {
	*ff=YY.raw_data + offset;
	traits.read(ff,element);
	offset += traits.written_size(element);
       	key = traits.build(element);
	for( int j=0;j<indent; j++)
	  std::cerr << "  ";
        traits.dump(key);
      }
    std::cerr << std::endl;
    std::cerr << "********************************************\n";
  }

   void dump_this_node(Node_entry& YY,  int levelX,  int indent){
     int j;
    iterator i;
    std::cerr << "============================================\n";
    std::cerr << " dump of node\n";
    for (i = &YY.ele[0];i !=&YY.ele[YY.eles]; i++)  {
      for(j=0;j<indent; j++)
	std::cerr << "  ";
      std::cerr << "pnext " << (*i).pnext;
      std::cerr << ((*i).deleted ? "*" : "") << "\t";
      for(j=0;j<indent; j++)
	std::cerr << "  ";
      traits.dump((*i).key,0);
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
    std::cerr << "============================================\n";
  }

  
   void dump_node( int levelX,  int indent){
     int j;
     int act_level=XX->level; long act_pthis=XX->pthis;
    dump_this_node(*(dynamic_cast<Node_entry*>(XX)), levelX, indent);
    std::cerr << std::endl;
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
    std::cerr << "Dump of children: \n";
     int l=0;
    while (l<XX->eles)
      {
	std::cerr << "parent has number: " << act_pthis << std::endl;
	for(j=0;j<indent; j++)
	  std::cerr << "  ";
	traits.dump((*(dynamic_cast<Node_entry*>(XX))).ele[l].key);
	std::cerr << "going down: " << act_pthis << std::endl;
	std::cerr << "level: " << act_level-1 << " pnext: " 
	 << (*(dynamic_cast<Node_entry*>(XX))).ele[l].pnext << std::endl;
	if(get(act_level-1,(*(dynamic_cast<Node_entry*>(XX))).ele[l].pnext))
	  dump(levelX,indent+1);
	get(act_level,act_pthis);
	std::cerr << "going up: " << act_pthis << std::endl;
	for(j=0;j<indent; j++)
	  std::cerr << "  ";
	traits.dump((*(dynamic_cast<Node_entry*>(XX))).ele[l].key);
	std::cerr << std::endl;
	l++;
      }
    std::cerr << std::endl;
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  }
  
  
   void dump( int levelX,  int indent) {
    if(dynamic_cast<Data_entry*>(XX))
      dump(*(dynamic_cast<Data_entry*>(XX)), levelX, indent);
    if(dynamic_cast<Node_entry*>(XX))
      dump_node(levelX, indent);
  }
  
  
   void dump_vertex(Data_entry& YY, bool ok,  int levelX,  int speach)
  {
  }

  void dump_vertex_simple(Node_entry & YY)
  {
    int i;
    std::cerr << std::endl;
    std::cerr << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
    std::cerr << "Output of a node\n";
    for(i=0;i<YY.eles;i++)
      traits.dump(YY.ele[i].key);
    std::cerr << std::endl;
    std::cerr << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
  }

  
   void dump_vertex_node(bool ok,  int levelX,  int speach)
  {
    Key elem;
     int act_level=XX->level;
    long act_pthis=XX->pthis;
    int i, pnext;
    std::cerr << std::endl;
    std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++\n";
    std::cerr << "Output of a node\n";
    for(i=0;i<XX->eles;i++)
      {
	pnext=(*(dynamic_cast<Node_entry*>(XX))).ele[i].pnext;
	std::cerr << "Adress is " << pnext << std::endl;
	if (pnext >= 0)
	  {
	    get(act_level-1,pnext);
	    if(dynamic_cast<Data_entry*>(XX))
	      elem = build(*(dynamic_cast<Data_entry*>(XX)));
	    if(dynamic_cast<Node_entry*>(XX))
	      elem = build( *(dynamic_cast<Node_entry*>(XX)));
	    get(act_level,act_pthis);
	    if (!(traits.equal_key(elem,
		 (*(dynamic_cast<Node_entry*>(XX))).ele[i].key))) 
	      {
		if (speach > 0) 
		  {
		    std::cerr << "!!! Error !!! Invalid parent key ";
		    traits.dump((*(dynamic_cast<Node_entry*>(XX))).ele[i].key);
		    std::cerr << " instead of ";
		    traits.dump(elem);
		  }
		std::cerr << "not ok since invalid parent key" << std::endl;
		ok = false;
		(*(dynamic_cast<Node_entry*>(XX))).ele[i].key = elem;
		put();
	      }
	  } 
	else 
	  {
	    if (speach > 0) 
	      {
		std::cerr << "!!! Error !!! Invalid tree pointer ";
		std::cerr << "(past end of tree structure) in node ";
	      }
	    std::cerr << "not ok since invalid tree pointer" << std::endl;
	    ok = false;
	  }
      }
     int c,j;
    for (c = 0, j = 0; j<IO_max_cap_nodes; ++j)  
      {
	if (!(*(dynamic_cast<Node_entry*>(XX))).ele[j].deleted)
	  c++;
      }
    if (c != XX->eles) 
      {
	if (speach > 0) 
	  {
	    std::cerr << "!!! Error !!! Wrong number of children in node ";
	    std::cerr << ", " << c << " instead of " << XX->eles << std::endl;
	  }
	std::cerr << "not ok since wrong number of children " << std::endl;
	ok = false;
      }
    if (!ok) put();
    for (j=0;j<IO_max_cap_nodes; j++)  
      {
	if(!(*(dynamic_cast<Node_entry*>(XX))).ele[j].deleted)
	get(act_level-1, (*(dynamic_cast<Node_entry*>(XX))).ele[j].pnext);
	isvalid(levelX, speach);
      }
    get(act_level,act_pthis);
    if (speach > 1)
      std::cerr << "\tnode number " << act_pthis << (ok ? " ok" : " not ok") 
		<< std::endl;
    std::cerr << std::endl;
    std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++\n";
  }

  
  // the tree structure is tested.
   bool isvalid( int levelX = 0,  int speach = 0) 
  {
    bool ok = true;
    Key elem;
    
    if (speach > 1)
      std::cerr << "Testing\tnode number " << XX->pthis << std::endl;
    if(!XX->level>0)
      {   
	if(XX->overfilled(XX->eles)){
	  if (speach > 0) {
	    std::cerr << "!!! Error !!! Too many children in node ";
	    std::cerr << ", " << XX->eles << " instead of " 
		      << IO_max_cap_nodes << std::endl;
	  }
	  std::cerr << "not ok since overfilled" << std::endl;
	  ok = false;
	}
	else
	  if (XX->eles > IO_max_cap_nodes) 
	    {
	      if (speach > 0) 
		{
		  std::cerr << "!!! Error !!! Too many children in node ";
		  std::cerr << ", " << XX->eles << " instead of " 
			    << IO_max_cap_nodes << std::endl;
		}
	      std::cerr << "not ok since eles > IO_max_cap" << std::endl;
	      ok = false;
	    }
	if(XX->underfilled())
	  {
	    if (speach > 0) 
	      {
		std::cerr << "!!! Warning !!! Too few children in node ";
		std::cerr << ", " << XX->eles << " instead of " 
			  << IO_min_cap_nodes << std::endl;
		std::cerr << "Tree needs reorganization" << std::endl;
	      }
	    std::cerr << "not ok in underfilled";
	    ok = false;
	  }
      }
    if(dynamic_cast<Data_entry*>(XX))
      dump_vertex(*(dynamic_cast<Data_entry*>(XX)), ok, levelX, speach);
    if(dynamic_cast<Node_entry*>(XX))
      dump_vertex_node( ok, levelX, speach);
    return ok;
  }
  
};

//this structure contains the necessary information of the tree
// that is its heights, if it is empty and the address of the root
class header_type {
public:
  int tree_height;
  bool empty;
  long rootpthis;
  //  size_t size() { return sizeof(int) + sizeof(bool) + sizeof(long)+3;}
  size_t size() { return sizeof(int) + sizeof(bool) + sizeof(long);}

  void read(char **s){
    int sint=(int)sizeof(int);
    int slong=(int)sizeof(long);
    int sbool=(int)sizeof(bool);
    char *from_int=new char[sint];
    char *from_bool=new char[sbool];
    char *from_long=new char[slong];
    int i,r=0;
    for (i=0; i<sint; i++)
      from_int[i] = (*s)[i];
    r += sint;
    tree_height=*((int *)from_int);

    for (i=0; i<slong; i++)
      from_long[i] = (*s)[i+r];
    rootpthis=*((long *)from_long);
    r += slong;

    for (i=0; i<sbool; i++)
      from_bool[i] = (*s)[i+r];
    empty=*((bool *)from_bool);

    delete[] from_bool;
    delete[] from_int;
    delete[] from_long;
  }
  void write (char **s) {
    int sint=(int)sizeof(int);
    int slong=(int)sizeof(long);
    int sbool=(int)sizeof(bool);

    char *from_int=new char[sint];
    char *from_bool=new char[sbool];
    char *from_long=new char[slong];
    int i,r=0;
    memcpy(from_int,(char *)(&tree_height),sint);
    for (i=0; i<sint; i++)
      (*s)[i] = from_int[i];
    r += sint;
    memcpy(from_long,(char *)(&rootpthis),slong);
     for (i=0; i<slong; i++)
       (*s)[i+r] = from_long[i];
     r += slong;

    memcpy(from_bool,(char *)(&empty),sbool);
     for (i=0; i<sbool; i++)
       (*s)[i+r] = from_bool[i];
     r += sbool;

     delete[] from_bool;
     delete[] from_int;
     delete[] from_long;
  }
};


  
// Implementation of the R_Tree. 
template <class Traits, class R_tree_index, class R_tree_storage>
class R_tree {
public:
  const static bool reinsertions = R_tree_index::reinsertions;
  const static  int IO_min_cap_nodes=R_tree_storage::IO_min_cap_nodes;
  const static  int IO_max_cap_nodes=R_tree_storage::IO_max_cap_nodes;
  const static  int IO_min_cap_leaves=
    R_tree_storage::IO_min_cap_leaves;
  const static  int IO_page_size=R_tree_storage::IO_page_size;
  const static bool headerextern=R_tree_storage::headerextern;
  typedef typename Traits::Data Data;
  typedef typename Traits::Key Key;
  typedef typename R_tree_storage::IO_tree_traits_nodes IO_intf_nodes;
  typedef typename R_tree_storage::IO_tree_traits_leaves IO_intf_leaves;
  typedef typename R_tree_storage::IO_tree_traits_nodes IO_tree_traits_nodes;
  typedef typename R_tree_storage::IO_tree_traits_leaves IO_tree_traits_leaves;
  typedef R_tree_node<Traits,R_tree_index,
    IO_tree_traits_nodes, IO_tree_traits_leaves, IO_min_cap_nodes, 
    IO_max_cap_nodes, IO_min_cap_leaves,
    IO_page_size> Node;
  typedef typename Node::Tree Tree;
  typedef typename Node::Leaf_data Leaf_data;
  typedef typename Node::v_type v_type;
  typedef typename Node::Node_entry Node_entry;
  typedef typename Node::Data_entry Data_entry;
  typedef typename Node::size_type size_type;
  typedef R_tree<Traits, R_tree_index,R_tree_storage> RT;
  Traits traits;

protected:

  Leaf_data leaf_data;
  Tree tree;
  std::fstream head;
  header_type header;
   int header_number;
  char fhead[80], ftree[80], fdata[80], fleaf_data[80];
  Node *root;
   int speach;
   int rootlevel;
  long rootpthis;
  char **header_string;
  int header_size;

public:
  std::vector<bool> Reinsert;

  R_tree() {root=0;}

  // initialization of the R_Tree from the files h,t,leaf_data_file.
  // Either these files are empty or they were created in a previous run.
  // Otherwise no guarantee can be given.
  R_tree(char *h, char *t, char *leaf_data_file) {
    //allocate memory for the header
    header_string = new (char *);
    *header_string= new char[header.size()];
    header_size=header.size();
    open(h, t, leaf_data_file);
    if(header.empty){
      root = new Node(0); //set level of root to 0
      root->open(tree, leaf_data); 
      header.rootpthis=rootpthis=root->XX->pthis=leaf_data.get_pos();
      CGAL_DOGN_ControlB(
		    std::cerr << "Rootpthis " << rootpthis;
		    )
      (*(dynamic_cast<Data_entry*>(root->XX))).write((root->data_space));
      leaf_data.insert(root->XX->pthis,(root->data_space));
      rootlevel=0;
    }
    else
      {
	if(header.tree_height==0){
	  rootlevel=0;
	  root = new Node(0);
	}
	else 
	  {
	    root = new Node(header.tree_height);
	    rootlevel = header.tree_height;
	    rootpthis=header.rootpthis;
	  }
      }
    root->open(tree, leaf_data); 
    root->get(rootlevel,rootpthis);
    int i;
    if(reinsertions)
    {
      for(i=0;i<=header.tree_height;i++)
	Reinsert.push_back(true);
	}
    else
      for(i=0;i<=header.tree_height;i++)
	Reinsert.push_back(false);
  }
  virtual ~R_tree() 
    {
      if(header.empty)
	{
	  root->get(rootlevel,rootpthis);
	  root->erase();
	}
      if(headerextern)
	{
	  header.write(header_string);
	  head.seekp(0);
	  head.write(*header_string,header_size);
	  head.close();
	}
      tree.close();
      leaf_data.close();
      delete[] *header_string;
      delete header_string;
      if(root!=0)
        delete root; 
      CGAL_DOGN_ControlB(
		    std::cerr << "R_tree destructor called\n";
		    )
    }
  
protected:
  // the files are opend and the tree initialized
  void open(char *h, char *t, char *leaf_data_file) {
    Node_entry n_entry;
    Data_entry d_entry;
    if(headerextern){
      //fake the SGI CC - compiler: it does not want to open empty 
      //files for reading 
      std::fstream fake_open;
      fake_open.open(h,std::fstream::out|std::fstream::app);
      CGAL_DOGN_Control(
		   if(!fake_open.is_open())
		   std::cerr << "fake open is not open";
		   else
		   std::cerr << "header fake open is open";
		   )
      fake_open.close();
      head.open(h,(std::fstream::in|std::fstream::out|std::fstream::ate));
      CGAL_DOGN_Control(
	        //SGI CC compiler does not like to open empty files for input
		   if(!head.is_open()) 
		   {
	    std::cerr << "Header file did not open! Data will not be saved.\n";
		   }
		   )
    }
    tree.open(n_entry.size(),t, 
             (std::fstream::binary|std::fstream::in|std::fstream::out));
    leaf_data.open(d_entry.size(),leaf_data_file, 
             (std::fstream::binary|std::fstream::in|std::fstream::out));
    strcpy(fhead, h);
    strcpy(ftree, t);
    strcpy(fleaf_data, leaf_data_file);
    
    if(headerextern){
      head.seekp(0, std::fstream::end);
      int teller;
      teller = head.tellp();
      if (0>=teller || leaf_data.number_of_elements()==0) { // no entries
	CGAL_DOGN_ControlB(	
		      std::cerr<<"no entries for header\n";
		      )
	header.empty=true;
	header.rootpthis = -1;
	header.tree_height = 0;
	header.write(header_string);
	head.write(*header_string,header_size);
      }
      else{
	CGAL_DOGN_ControlB(
		      std::cerr<<"entries for header\n";
		      )
	head.seekg(0);
	head.get(*header_string,header_size);
	header.read(header_string);
      }
    }
    else{
     header.empty=true;
     header.tree_height = 0;
     header.rootpthis =-1;
    }
    CGAL_DOGN_ControlB(
		  std::cerr << "Header information: \n";
		  std::cerr << "empty: " << header.empty << " treeheight " 
		  << header.tree_height;
		  std::cerr << "Rootpthis: " << header.rootpthis << std::endl;
		  )
    empty = header.empty;
    tree_height = header.tree_height;
  }

  // the offset of rootXX in the file is returned.
   long get_pnext(Node_entry& rootXX)
  {    
    return rootXX.ele[0].pnext;
  }
  
public:
  bool empty;
   int tree_height;


  // insertion of an element into the tree
   void insert(Data& elem) {
    std::vector<v_type> children, root_children;
    std::pair< int,long> *IndexPthis= new std::pair< int,long>[rootlevel+1];
    typename std::vector<v_type>::iterator i;
    v_type v;
    
    v.deleted = false;
    v.key= traits.build(elem);
    root->get(rootlevel,rootpthis);
    root->insert_star(elem, v, children,root_children, Reinsert,
		      rootlevel,IndexPthis);
    if(header.empty==true)
    {
      header.empty=false;      
      if(headerextern){
	header.write(header_string);
	head.seekp(0);
	head.write(*header_string,header_size);
      }
    }

    root->get(rootlevel,rootpthis);
    CGAL_DOGN_Insert(
	   std::cerr << "Output of the root element right after insert" 
	             << std::endl;
	   if(rootlevel>0)
	   {
	     if(dynamic_cast<Node_entry*>(root->XX))
	      root->dump_vertex_simple(*(dynamic_cast<Node_entry*>(root->XX)));
	   }
	   std::cerr << "end\n";
		)
    while(!root_children.empty()){
      children.push_back(*root_children.begin());
      root_children.erase(root_children.begin());
    }
    if (!children.empty()) {
      if(reinsertions)
	Reinsert.push_back(true);
      else
	Reinsert.push_back(false);
      header.tree_height++;
      rootlevel++;
      if(dynamic_cast<Node_entry*>(root->XX)){
	root->XX->pthis = tree.get_pos();
	(*(dynamic_cast<Node_entry*>(root->XX))).write((root->node_space));
	tree.insert(root->XX->pthis, (root->node_space));
	(*children.begin()).pnext = root->XX->pthis;
	root->XX->pthis = rootpthis;
	root->clear();
	root->XX->level=rootlevel;
      }
      else
      {
        delete root;
	root = new Node(1); 
	root->open(tree, leaf_data);
	header.rootpthis=rootpthis=root->XX->pthis = tree.get_pos();
	(*(dynamic_cast<Node_entry*>(root->XX))).write((root->node_space));
	tree.insert(root->XX->pthis,(root->node_space));
	CGAL_DOGN_Insert(
		    std::cerr << "root pthis moved to : " << root->XX->pthis 
		    << std::endl;
		    )
      }
      if(dynamic_cast<Node_entry*>(root->XX))
	{
	  for (i = children.begin(); i != children.end(); i++){
	    root->insert_star((*i), *(dynamic_cast<Node_entry*>(root->XX) ));
	    CGAL_DOGN_Insert(
			std::cerr << "number " << i->pnext << std::endl;
			)
	  }
	  root->put();
       	}

      if(headerextern){
	header.write(header_string);
	head.seekp(0);
	head.write(*header_string,header_size);
      }
    }
    CGAL_DOGN_Insert(
	    dump(0,0);
	    std::cerr << "Output of the root element right after insert" 
	              << std::endl;
	    if(rootlevel>0)
	    {
	     if(dynamic_cast<Node_entry*>(root->XX))
	      root->dump_vertex_simple(*(dynamic_cast<Node_entry*>(root->XX)));
	    }
	    std::cerr << "end\n";
	        )
    delete[] IndexPthis;
  }


  // returns true if a is contained in an key.
   bool find_key_include(Key& a) {
      root->get(rootlevel,rootpthis);
    return root->find_key_include(a);
  }
  
  // returns true if a intersects another element.
    bool find_key_intersect(Key& a) {
      root->get(rootlevel,rootpthis);
    return root->find_key_intersect(a);
  }
  

  // all elements with key a are deleted and returned
   bool delete_key(Key& a, Data& result) 
  { 
    bool found=false;
    std::pair< int,long> *IndexPthis= new  std::pair< int,long>[rootlevel+1];
    std::vector<v_type> children;
    std::vector<v_type>::iterator i;
    root->get(rootlevel,rootpthis);
    if(dynamic_cast<Data_entry*>(root->XX))
      found = root->XXdelete_key_data(a,result,IndexPthis,Reinsert,rootlevel, 
				      children);
    else
      if(dynamic_cast<Node_entry*>(root->XX))
	found = root->delete_key(a,result,IndexPthis,Reinsert,rootlevel,
				 children);
    root->get(rootlevel,rootpthis);
    if (!children.empty()) {
      header.empty=false;
      header.tree_height++;
      rootlevel++;
      // copy root at another place, the place is given to the first child.
      if(dynamic_cast<Node_entry*>(root->XX)){
	root->XX->pthis=tree.get_pos();
	(*(dynamic_cast<Node_entry*>(root->XX))).write((root->node_space));
	tree.insert(root->XX->pthis, (root->node_space));
	(*children.begin()).pnext = root->XX->pthis;
	root->XX->pthis = rootpthis;
	root->clear();
	root->XX->level=rootlevel;
      }
      else
      {
        delete root;
	root = new Node(1); //give root level 1 now 
	root->open(tree, leaf_data);
	root->XX->pthis=tree.get_pos();
	(*(dynamic_cast<Node_entry*>(root->XX))).write((root->node_space));
	tree.insert(root->XX->pthis, (root->node_space));
	CGAL_DOGN_Delete(
		    std::cerr << "number: " << root->XX->pthis << std::endl;
		    )
      }
      if(dynamic_cast<Node_entry*>(root->XX))
	{
	  for (i = children.begin(); i != children.end(); i++){
	    root->insert_star((*i), *(dynamic_cast<Node_entry*>(root->XX) ));
	    CGAL_DOGN_Delete(
			std::cerr << "number " << i->pnext << std::endl;
			)
	  }
	  root->put();
       	}
      CGAL_DOGN_Delete(
		  dump(0,0);
		  )
      if(headerextern){
	header.write(header_string);
	head.seekp(0);
	head.write(*header_string,header_size);
      }
    }
    else{
      if(root->XX->eles<1 && dynamic_cast<Data_entry*>(root->XX))
	{   
	  root->erase();
	  header.empty=true;
	  delete root;
	  root = new Node(0);
	  rootlevel=0;
	  root->open(tree, leaf_data);
	  root->XX->pthis=-1;
	  root->put(-1);
	  header.rootpthis=rootpthis=root->XX->pthis;
	  Reinsert.erase(Reinsert.end());
	  if(reinsertions)
	    Reinsert.push_back(true);
	  else
	    Reinsert.push_back(false);
	  if(headerextern){
	    header.write(header_string);
	    head.seekp(0);
	    head.write(*header_string,header_size);
	  }
	}
      while(root->XX->eles<=1 && dynamic_cast<Node_entry*>(root->XX) )
	{
	  header.tree_height--;
	  rootlevel--;
	  Reinsert.erase(Reinsert.end());
	  Node *n= new Node(root->XX->level-1);
	  n->open(tree, leaf_data);
	  long pnext = get_pnext(*(dynamic_cast<Node_entry*>(root->XX)));
	  n->get(root->XX->level-1, pnext);
	  n->erase();
	  n->XX->pthis=-1;
	  root->erase();
	  delete root;
	  n->put(-1);
	  header.rootpthis=rootpthis=n->XX->pthis;
	  root = n;
	  if(headerextern){
	    header.write(header_string);
	    head.seekp(0);
	    head.write(*header_string,header_size);
	  }
	}
    }
    CGAL_DOGN_Delete(
		std::cerr << " $$$$$$$$AFTER DELETE KEY :";
		traits.dump(a);
		dump(0,0);
		)
    delete [] IndexPthis;
    return found;
  }



  //the level of the root is returned
  int get_rootlevel()
  {
    return rootlevel;
  }
  
  bool get_reinsertion_flag(int the_level)
  {
    if (the_level <= rootlevel)
      return Reinsert[the_level];
    return false;
  }

  //if value=true => the next time a node in level the_level
  //is split, the vertices are reinserted on the same level
  //if value=false => reinsertion at level the_level is prohibited
  bool set_reinsertion_level( int the_level,bool value)
  {
    if(the_level<=rootlevel && reinsertions)
      {
	Reinsert[the_level]=value;
	return true;
      }
    return false;
  }


  // information about the tree is given out
   void dump( int levelX=0,  int indent=0) {
     root->get(rootlevel,rootpthis);
    
    std::cerr << "Dump of R-tree" << std::endl;
    std::cerr << "_________________________" << std::endl;
    std::cerr << "Tree Height " << header.tree_height << std::endl;
    std::cerr << "Number of nodes in tree\t\t" << tree.number_of_elements() 
	      << std::endl;
    std::cerr << "Number of leaves in tree\t\t" 
	      << leaf_data.number_of_elements() << std::endl;
    std::cerr << "Associated files are" << std::endl;
    std::cerr << "\tHeader\t\t\t" << fhead << std::endl;
    std::cerr << "\tTree\t\t\t" << ftree << std::endl;
    std::cerr << "\tData\t\t\t" << fleaf_data << std::endl;
    if(!header.empty)
      root->dump(levelX,0);
  }
  
  // return true if the tree structure is valid
  bool isvalid( int levelX = 0) {
      root->get(rootlevel,rootpthis);
    return root->isvalid(levelX, speach);
  }

  // class iterator is used to iterate through the data elements
class iterator;
friend class iterator;
  public:

  iterator begin()
  { 
    root->get(rootlevel,rootpthis);
    Key k=root->bounding_box();
    return iterator( k, this, 1,0); 
  }
  iterator end()
  {   
    root->get(rootlevel,rootpthis);
    Key k=root->bounding_box();
    return iterator( k, this, 0,0);
  }

  iterator begin(Key start)
  { 
    root->get(rootlevel,rootpthis);
    return iterator(start, this, 1,0); 
  }
  iterator end(Key stop)
  { 
    root->get(rootlevel,rootpthis);
    return iterator(stop, this, 0,0);
  }
  iterator begin_intersect(Key start)
  { 
    root->get(rootlevel,rootpthis);
    return iterator(start, this, 1,0); 
  }
  iterator end_intersect(Key stop)
  { 
    root->get(rootlevel,rootpthis);
    return iterator(stop, this, 0,0);
  }


  iterator begin_enclose(Key start)
  { 
    root->get(rootlevel,rootpthis);
    return iterator(start, this, 1,1); 
  }
  iterator end_enclose(Key stop)
  { 
    root->get(rootlevel,rootpthis);
    return iterator(stop, this, 0,1);
  }

  iterator begin_compare(Key start)
  { 
    root->get(rootlevel,rootpthis);
    return iterator(start, this, 1,2); 
  }

  iterator end_compare(Key stop)
  { 
    root->get(rootlevel,rootpthis);
    return iterator(stop, this, 0,2);
  }


  class iterator{
  protected:
    std::list<Node_info *> node_stack;
    std::list<Node_info *> end_node_stack;

    // this type is 0 for intersection, 1 for enclosing, 2 for compare
     int iterator_compare_type;
     int next_element, end_next_element;
     int offset, end_offset;
    Data iterator_element;
    Key win;
    Data_entry actual_data, end_actual_data;
    RT *rtree;
    bool valid;
    bool past_the_end;

  public:

    iterator():valid(false){ next_element=0;}
    ~iterator(){
      Node_info *ni;
      CGAL_DOGN_Search(  
		  print_node_infos(node_stack);
		  print_node_infos(end_node_stack);
		  )
      while(!node_stack.empty()){
	ni=(*node_stack.begin());
	node_stack.erase(node_stack.begin());
	delete ni;
      }
      while(!end_node_stack.empty()){
	ni=(*end_node_stack.begin());
	end_node_stack.erase(end_node_stack.begin());
	delete ni;
      }
    }


    iterator(Key v_win, RT *r_tree,  bool is_start, int intersection_type)
            : win(v_win) 
    { 
      //set the type of the compare function of the iterator.
      iterator_compare_type=intersection_type;
      Node_info *ni;
      r_tree->root->get(r_tree->rootlevel,r_tree->rootpthis);
      if(is_start)
	{
	  past_the_end = false; // mark iterator as non past the end iterator
	  while(!node_stack.empty()){
	    ni=(*node_stack.begin());
	    node_stack.erase(node_stack.begin());
	    delete ni;
	  }
	}
      else
	{
	  past_the_end = true; // mark iterator as past the end iterator
	  while(!end_node_stack.empty()){
	    ni=(*end_node_stack.begin());
	    end_node_stack.erase(end_node_stack.begin());
	    delete ni;
	  }
	}
      rtree=r_tree;
      if(intersection_type==0)
	{
	  CGAL_DOGN_Search(
		      std::cerr << "INTERSECTION ITERATOR\n";
		      )
	  if(is_start)
	    {
	      CGAL_DOGN_Search(
			  std::cerr << "iterator 2\n";
			  )
	      if (r_tree->traits.intersect(v_win,r_tree->root->bounding_box()))
		{
		  CGAL_DOGN_Search(   
			      std::cerr << "iterator 3\n";
			      )
		  if (r_tree->root->init_iterator(node_stack,actual_data,
						next_element,
						win,iterator_element,offset))
		    valid = true;
		  else
		    valid = false;
		  r_tree->root->get(r_tree->rootlevel,r_tree->rootpthis);
		  CGAL_DOGN_Search(
			      std::cerr << 
			      std::cerr<< "-------------Test actual data:\n";
			      rtree->root->test_actual_data(actual_data);
			      std::cerr << "------------end\n";
			      )
		    r_tree->root->get(r_tree->rootlevel,r_tree->rootpthis);
		}
	      else
		valid = false;
	    }
	  else
	    {
	      CGAL_DOGN_Search(
			  std::cerr << "iterator 5\n";
			  )
	      if (r_tree->traits.intersect(v_win,r_tree->root->bounding_box()))
		{
		  CGAL_DOGN_Search(
			      std::cerr << "iterator 6\n";
			      )
		  if(r_tree->root->reverse_init_iterator(end_node_stack,
						       end_actual_data,
						       end_next_element,win,
						       iterator_element,
						       end_offset))
		    valid = true;
		  else
		    valid = false;
		  r_tree->root->get(r_tree->rootlevel,r_tree->rootpthis);
		}
	      else
		valid = false;
	    }
	}
      else
	{
	  if(intersection_type == 1)
	    {
	      CGAL_DOGN_Search(
			  std::cerr << "ENCLOSING ITERATOR\n";
			  )
	      if(is_start)
		{
		  CGAL_DOGN_Search(
			      std::cerr << "iterator 2\n";
			      )
		 if (r_tree->traits.include(r_tree->root->bounding_box(),v_win))
		    {
		      if (r_tree->root->init_encl_iterator(node_stack,
							   actual_data,
							   next_element,
							   win,
							   iterator_element,
							   offset))
			valid = true;
		      else
			valid = false;
		      CGAL_DOGN_Search(
				std::cerr << 
				std::cerr<< "-------------Test actual data:\n";
				rtree->root->test_actual_data(actual_data);
				std::cerr << "------------end\n";
				  )
		    }
		  else
		    valid = false;
		  r_tree->root->get(r_tree->rootlevel,r_tree->rootpthis);
		}
	      else
		{
		 if (r_tree->traits.include(r_tree->root->bounding_box(),v_win))
		  {
		    if(r_tree->root->reverse_init_encl_iterator(end_node_stack,
							      end_actual_data,
							      end_next_element,
							      win,
							      iterator_element,
							      end_offset))
			valid = true;
		      else
			valid = false;
		      r_tree->root->get(r_tree->rootlevel,r_tree->rootpthis);
		    }
		  else
		    valid = false;
		}
	    }
	  else
	    {
	      if(intersection_type == 2)
		{
		  CGAL_DOGN_Search(
			      std::cerr << "COMPARE ITERATOR\n";
			      )
		  if(is_start)
		    {
		      CGAL_DOGN_Search(
				  std::cerr << "iterator 2\n";
				  )
		      if (r_tree->traits.compare(r_tree->root->bounding_box(),
						 v_win))
			{
			  if (r_tree->root->init_compare_iterator(node_stack,
							      actual_data,
							      next_element,
							      win,
							      iterator_element,
							      offset))
			    valid = true;
			  else
			    valid = false;
			}
		      else
			valid = false;
		      r_tree->root->get(r_tree->rootlevel,r_tree->rootpthis);
		    }
		  else
		    {
		      if (r_tree->traits.compare(r_tree->root->bounding_box(),
						 v_win))
			{
			  if(r_tree->root->reverse_init_compare_iterator(
						     end_node_stack,
						     end_actual_data,
						     end_next_element,win,
						     iterator_element,
						     end_offset))
			    valid = true;
			  else
			    valid = false;
			  r_tree->root->get(r_tree->rootlevel,
					    r_tree->rootpthis);
			}
		      else
			valid = false;
		    }
		}
	    }
	}
      CGAL_DOGN_Search(
		  print_node_infos(node_stack);
		  )
    }

    bool operator==(const iterator& x) const 
      { 
	if(valid && x.valid){
	  if(past_the_end == x.past_the_end)
	    {
	      rtree->root->get(rtree->rootlevel,rtree->rootpthis);
	      if(rtree->traits.equal_key(rtree->traits.build(
	         x.iterator_element),rtree->traits.build(iterator_element)))
		return true;
	    }
	  return false;
	}
	if(!valid && !x.valid){
	  if(rtree->traits.equal_key(win,x.win))
	    return true;
	  return false;
	}
	if(!valid){
	  if(!past_the_end && x.past_the_end)
	    if(rtree->traits.equal_key(win,x.win))
	      return true;
	  return false;
	}
	else
	  if(past_the_end && !x.past_the_end)
	    if(rtree->traits.equal_key(win,x.win))
	      return true;
	return false;
      }

    bool operator!=(const iterator& x) const 
    { 
      return !(*this==x);
    }

    iterator& operator=(const iterator& x)
      {
	rtree->root->get(rtree->rootlevel,rtree->rootpthis);
	if (&x==this)
	  return *this;
	iterator_compare_type=x.iterator_compare_type;
	Node_info *ni;
	while(!node_stack.empty()){
	  ni=(*node_stack.begin());
	  node_stack.erase(node_stack.begin());
	  delete ni;
	}
	std::list<Node_info *>::const_iterator node_it= x.node_stack.begin();
	for(;node_it!=x.node_stack.end();node_it++){
	  ni=new Node_info;
	  ni->level=(*node_it)->level;
	  ni->pnext=(*node_it)->pnext;
	  node_stack.push_back(ni);
	}
	while(!end_node_stack.empty()){
	  ni=(*end_node_stack.begin());
	  end_node_stack.erase(end_node_stack.begin());
	  delete ni;
	}
	node_it= x.end_node_stack.begin();
	for(;node_it!=x.end_node_stack.end();node_it++){
	  ni=new Node_info;
	  ni->level=(*node_it)->level;
	  ni->pnext=(*node_it)->pnext;
	  end_node_stack.push_back(ni);
	}
	next_element=x.next_element;
	offset=x.offset;
	iterator_element=x.iterator_element;
	win=x.win;
	for(int i=0;i<IO_page_size;i++)
	  actual_data.raw_data[i]=x.actual_data.raw_data[i];
	actual_data.eles=x.actual_data.eles;
	rtree=x.rtree;
	valid = x.valid;
	past_the_end = x.past_the_end;
	return *this;
    }
      
    Data operator*() const 
    { 	  
      return iterator_element;
    }

    bool is_valid()
    {
      return valid;
    }

    iterator& operator++(){
      if(past_the_end)
	{
	  valid=false;
	  return *this;
	}
      rtree->root->get(rtree->rootlevel,rtree->rootpthis);
      CGAL_DOGN_Search(
		  std::cerr << "++++++operator\n";
		  std::cerr<< "-------------Test actual data:\n";
		  rtree->root->test_actual_data(actual_data);
		  std::cerr << "------------end\n";
		  print_node_infos(node_stack);
		  )
      rtree->root->get(rtree->rootlevel,rtree->rootpthis);
      if(iterator_compare_type==0) //intersection iterator
	{
	  CGAL_DOGN_Search(
		      std::cerr << "INTERSECTION ITERATOR\n";
		      )
	  if(rtree->root->next_element(node_stack, actual_data, next_element, 
				       win, iterator_element, offset))
	    valid=true;
	  else
	    valid=false;
	}
      else
	{
	  if(iterator_compare_type==1) //enclosing iterator
	    {
	      CGAL_DOGN_Search(
			  std::cerr << "ENCLOSE ITERATOR\n";
			  )
	      if(rtree->root->next_encl_element(node_stack, actual_data, 
						next_element, 
						win, iterator_element, offset))
		valid=true;
	      else
		valid=false;
	    }
	  else
	    {
	      if(iterator_compare_type==2) //compare iterator
		{
		  CGAL_DOGN_Search(
			      std::cerr << "COMPARE ITERATOR\n";
			      )
	       if(rtree->root->next_compare_element(node_stack, actual_data,
						    next_element, 
						    win, iterator_element, 
						    offset))
		    valid=true;
		  else
		    valid=false;
		}
	    }
	}
      CGAL_DOGN_Search(
		  std::cerr << "nach operator++\n";
		  std::cerr<< "-------------Test actual data:\n";
		  rtree->root->test_actual_data(actual_data);
		  std::cerr << "------------end\n";
		  )
      return *this;
    }

    //  iterator operator++(int) { 
    //  iterator tmp = *this;
    //  ++*this;
    //  return tmp;
    // }
  };
};


CGAL_END_NAMESPACE

#endif

