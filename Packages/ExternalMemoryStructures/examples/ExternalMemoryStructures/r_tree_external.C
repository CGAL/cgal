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
// file          : include/CGAL/R_Tree/examples/ExternalMemoryStructures/r_star_tree_external.C
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
#include <CGAL/R_tree_index.h>
#include <CGAL/R_tree_traits_implementation.h>
#include <CGAL/R_tree_external_storage.h>
#define NUMBER 16


#if defined(__sgi) && !defined(__GNUC__)
#error "the fstream.h version for SGI CC compilers is buggy"
#else


//Definition of the data type
struct Data{
public:
  typedef CGAL::R_tree_key_2 Key;
  Key key;
  size_t size(void) const {
    return key.size();
  }
  void read(char ** s) {
    key.read(s);
    //    *this = *(Data *)s;
  }
  void write(char ** s) {
    key.write(s);
    //    memcpy(s,this,sizeof(*this));
  }    
  void dump(int level =0){
    key.dump();
  }
};

//definition of the R_tree_traits - depending on the Data type
typedef CGAL::R_tree_traits<Data> TTraits; 

typedef Data::Key Key;

/* definition of the R_Tree that contains Data elements, uses 
   Guttmann index structure and stores the elements in 
   external memory */
typedef CGAL::R_tree<TTraits, CGAL::R_tree_index<TTraits>,
  CGAL::R_tree_external_storage> R_Tree_Inst;



int main() {
  /* creation of R_tree associated to the files: */
  R_Tree_Inst r_tree("__tree.head","__tree.dat", "__leaf_data.dat");

  TTraits traits;
  Data  elem; 
  int k; 
  Key key= Key(0,2,0,2);
  std::vector<Data > source;

  //create the squares
  for (k=0;k<NUMBER;++k) {
    elem.key=key;
    source.push_back(elem);
    key.xmin++; key.ymin++; key.xmax++; key.ymax++;
  }
  std::cerr << "Main 2\n";
  /* Output of the restored tree */
  std::cerr << "\n Output of the restored tree\n";
  r_tree.dump();
  /* Insertion of elements */
  for (k=0;k<NUMBER;++k) {
    r_tree.insert(source[k]);
  }
  r_tree.dump();


  /* Iteration through all elements of the tree */
  std::cerr<< "\n Iteration through all elements of the tree\n";
  R_Tree_Inst::iterator it_begin=r_tree.begin();
  R_Tree_Inst::iterator it_end=r_tree.end();
  while(it_begin != it_end){
    std::cerr << std::endl;
    (*it_begin).dump();
    ++it_begin;
  }
  std::cerr << std::endl;
  (*it_begin).dump();
  std::cerr<< "\n End of iteration through all elements of the tree\n";

  /* Iteration through all elements of the tree that have non empty
     intersection with source[0].key=(0,2,0,2) */
  std::cerr<< "\n Iteration through all elements of the tree\n";
  std::cerr<< "that have non empty intersection with source[0].key=(0,2,0,2)\n";
  it_begin=r_tree.begin(source[0].key);
  it_end=r_tree.end(source[0].key);
  while(it_begin != it_end){
    std::cerr << std::endl;
    (*it_begin).dump();
    ++it_begin;
  }
  std::cerr<< "\n End of iteration through the query elements of the tree\n";


  /* Iteration through all elements of the tree that ENCLOSE
     source[2].key=(2,4,2,4) */
  std::cerr<< "\n Iteration through all elements of the tree\n";
  std::cerr<< "that enclose source[2].key=(2,4,2,4)\n";
  it_begin=r_tree.begin_enclose(source[2].key);
  it_end=r_tree.end_enclose(source[2].key);
  while(it_begin != it_end){
    std::cerr << std::endl;
    (*it_begin).dump();
    ++it_begin;
  }
  std::cerr<< "\n End of iteration through the query elements of the tree\n";

  /* Iteration through all elements of the tree that COMPARE
     source[4].key=(4,8,4,8) */
  std::cerr<< "\n Iteration through all elements of the tree\n";
  std::cerr<< "that compare source[4].key=(4,8,4,8)\n";
  it_begin=r_tree.begin_compare(source[4].key);
  it_end=r_tree.end_compare(source[4].key);
  while(it_begin != it_end){
    std::cerr << std::endl;
    (*it_begin).dump();
    ++it_begin;
  }
  std::cerr<< "\n End of iteration through the query elements of the tree\n";

  std::cerr << "\n Check for elements that intersects source[1].key=(1,2,1,2)\n";
  if(!r_tree.find_key_intersect(source[1].key))
    {
      std::cerr << "\n no key intersection of "; 
      traits.dump(traits.build(source[1]));
    }
  else
    std::cerr << "\n key intersection = true";

  std::cerr << "\n Check for elements that intersects source[1].key=(1,2,1,2)\n";
  if(!r_tree.find_key_include(source[1].key))
    {
      std::cerr << "\n no key include of"; 
      traits.dump(traits.build(source[1]));
    }
  else
    std::cerr << "\n key include = true";

  Data data_del;
  std::cerr << "\n Deletion of all data with key source[0].key=(0,2,0,2)\n";
  while(r_tree.delete_key(source[0].key,data_del))
    traits.dump(data_del.key);
  std::cerr << "\n Deletion of all data with key source[1].key=(1,3,1,3)\n";
  while(r_tree.delete_key(source[1].key,data_del))
    traits.dump(data_del.key);
  std::cerr << "\n Deletion of all data with key source[2].key=(2,4,2,4)\n";
  while(r_tree.delete_key(source[2].key,data_del))
    traits.dump(data_del.key);

  std::cerr << "\n Check for elements that intersect source[1].key=(1,2,1,2)\n";
  if(!r_tree.find_key_intersect(source[1].key))
    {
      std::cerr << "\n no key intersect of"; 
      traits.dump(traits.build(source[1]));
    }
  else
    std::cerr << "\n key intersect = true";
  r_tree.dump();
}
#endif
