// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ---------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.1-I-17$
// release_date  : $CGAL_Date: 1999/09/11 $
//
// file          : include/CGAL/R_Tree/examples/ExternalMemoryStructures/r_star_star_tree_internal.C
// chapter       : $CGAL_Chapter: Basic / External Data Structures $
// package       : $CGAL_Package: External Data Structures$
// source        : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Gabriele Neyer<neyer@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Peter Widmayer <widmayer@inf.ethz.ch>)
//
// Example instanciation of the R_Tree.h 
// ============================================================================
#include <CGAL/R_tree.h>
#include <CGAL/R_tree_key.h>
#include <CGAL/R_star_tree_index.h>
#include <CGAL/R_tree_traits_implementation.h>
#include <CGAL/R_tree_internal_storage.h>

//Definition of the data type
struct Data{
public:
  typedef CGAL::R_tree_key_2 Key;
  Key key;
  size_t size(void) const {
    return sizeof(*this);
  }
  void read(char ** s) {
    key.read(s);
  }
  void write(char ** s) {
    key.write(s);
  }    
  void dump(int level =0){
    key.dump();
  }
};

//definition of the R_tree_traits - depending on the Data type
typedef CGAL::R_tree_traits<Data> TTraits; 

typedef Data::Key Key;

/* definition of the R_Tree that contains Data elements, uses 
   Star Tree index structure and stores the elements in 
   internal memory */
typedef CGAL::R_tree<TTraits, CGAL::R_star_tree_index<TTraits>,
  CGAL::R_tree_internal_storage> R_Tree_Inst;



int main() {
  TTraits traits;
  Data  elem; 
  int k; 
  Key key= Key(0,4,4,8);
  std::vector<Data > source;

  /* creation of R_tree associated to the files: */
  R_Tree_Inst r_star_tree("__star_tree.head","__star_tree.dat", "__star_leaf_data.dat");

  //create some arbitrary bounding boxes
  elem.key=key;
  source.push_back(elem);

  key.xmin=2; key.ymin=2; key.xmax=9; key.ymax=9;
  elem.key=key;
  source.push_back(elem);

  key.xmin=5; key.ymin=6; key.xmax=5; key.ymax=6;
  elem.key=key;
  source.push_back(elem);

  key.xmin=0;key.ymin=0;key.xmax=2;key.ymax=2;
  elem.key=key;
  source.push_back(elem);
  
  key.xmin=12; key.ymin=3; key.xmax=14; key.ymax=9;
  elem.key=key;
  source.push_back(elem);

  key.xmin=6; key.ymin=2; key.xmax=7; key.ymax=9;
  elem.key=key;
  source.push_back(elem);

  key.xmin=3; key.ymin=8; key.xmax=9; key.ymax=12;
  elem.key=key;
  source.push_back(elem);

  key.xmin=0; key.ymin=1; key.xmax=5; key.ymax=4;
  elem.key=key;
  source.push_back(elem);

  key.xmin=2; key.ymin=3; key.xmax=3; key.ymax=7;
  elem.key=key;
  source.push_back(elem);

  key.xmin=1; key.ymin=5; key.xmax=3; key.ymax=8;
  elem.key=key;
  source.push_back(elem);

  key.xmin=4; key.ymin=9; key.xmax=10; key.ymax=11;
  elem.key=key;
  source.push_back(elem);

  key.xmin=6; key.ymin=5; key.xmax=7; key.ymax=8;
  elem.key=key;
  source.push_back(elem);

  key.xmin=3; key.ymin=3; key.xmax=8; key.ymax=8;
  elem.key=key;
  source.push_back(elem);

  key.xmin=4; key.ymin=4; key.xmax=9; key.ymax=9;
  elem.key=key;
  source.push_back(elem);

  /* Insertion of elements */
  for (k=0;k<14;++k) {
    r_star_tree.insert(source[k]);
  }
  r_star_tree.dump();


  /* Iteration through all elements of the tree */
  std::cerr<< "\n Iteration through all elements of the tree\n";
  R_Tree_Inst::iterator it_begin=r_star_tree.begin();
  R_Tree_Inst::iterator it_end=r_star_tree.end();
  while(it_begin != it_end){
    std::cerr << std::endl;
    (*it_begin).dump();
    ++it_begin;
  }
  std::cerr << std::endl;
  (*it_begin).dump();
  std::cerr<< "\n End of iteration through all elements of the tree\n";

  /* Iteration through all elements of the tree that have non empty
     intersection with source[0].key=(0,4,4,8) */
  std::cerr<< "\n Iteration through all elements of the tree\n";
  std::cerr<< "that have non empty intersection with source[0].key=(0,4,4,8)\n";
  it_begin=r_star_tree.begin(source[0].key);
  it_end=r_star_tree.end(source[0].key);
  while(it_begin != it_end){
    std::cerr << std::endl;
    (*it_begin).dump();
    ++it_begin;
  }
  std::cerr<< "\n End of iteration through the query elements of the tree\n";


  /* Iteration through all elements of the tree that ENCLOSE
     source[2].key=(5,6,5,6) */
  std::cerr<< "\n Iteration through all elements of the tree\n";
  std::cerr<< "that enclose source[2].key=(5,6,5,6)\n";
  it_begin=r_star_tree.begin_enclose(source[2].key);
  it_end=r_star_tree.end_enclose(source[2].key);
  while(it_begin != it_end){
    std::cerr << std::endl;
    (*it_begin).dump();
    ++it_begin;
  }
  std::cerr<< "\n End of iteration through the query elements of the tree\n";

  /* Iteration through all elements of the tree that COMPARE
     source[5].key=(6,2,7,9) */
  std::cerr<< "\n Iteration through all elements of the tree\n";
  std::cerr<< "that compare source[5].key=(6,2,7,9)\n";
  it_begin=r_star_tree.begin_compare(source[5].key);
  it_end=r_star_tree.end_compare(source[5].key);
  while(it_begin != it_end){
    std::cerr << std::endl;
    (*it_begin).dump();
    ++it_begin;
  }
  std::cerr<< "\n End of iteration through the query elements of the tree\n";

  std::cerr << "\n Check for elements that intersect source[1].key=(2,9,2,9)\n";
  if(!r_star_tree.find_key_intersect(source[1].key))
    {
      std::cerr << "\n no key intersection of "; 
      traits.dump(traits.build(source[1]));
    }
  else
    std::cerr << "\n key intersection = true";

  std::cerr << "\n Check for elements that intersects source[1].key=(2,9,2,9)\n";
  if(!r_star_tree.find_key_include(source[1].key))
    {
      std::cerr << "\n no key include of"; 
      traits.dump(traits.build(source[1]));
    }
  else
    std::cerr << "\n key include = true";

  Data data_del;
  std::cerr << "\n Deletion of all data with key source[0].key=(0,4,4,8)\n";
  while(r_star_tree.delete_key(source[0].key,data_del))
    traits.dump(data_del.key);
  std::cerr << "\n Deletion of all data with key source[1].key=(2,9,2,9)\n";
  while(r_star_tree.delete_key(source[1].key,data_del))
    traits.dump(data_del.key);
  std::cerr << "\n Deletion of all data with key source[2].key=(5,6,5,6)\n";
  while(r_star_tree.delete_key(source[2].key,data_del))
    traits.dump(data_del.key);

  std::cerr << "\n Check for elements that intersect source[1].key=(2,9,2,9)\n";
  if(!r_star_tree.find_key_intersect(source[1].key))
    {
      std::cerr << "\n no key intersect of"; 
      traits.dump(traits.build(source[1]));
    }
  else
    std::cerr << "\n key intersect = true";
}
