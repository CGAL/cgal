// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $Revision:
// release_date  : $Date:
//
// file          : $RCSfile: Data_structure_using_octree_3.h
// package       : Data_structure_for_queries
// maintainer    : 
// revision      : 
//
// author(s)     : Marie Samozino
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================
#ifndef DATA_STRUCTURE_USING_OCTREE
#define DATA_STRUCTURE_USING_OCTREE

#include <string>
#include <CGAL/basic.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <utility>
#include <map> 
#include <set>
#include <list>

#include "Constrained_Element.h"
#include "Constrained_Vertex.h"
#include "Constrained_Edge.h"
#include "Constrained_Facet.h"

#include "CGAL/squared_distance_3.h"

#include <CGAL/iterator.h>



//Corresponding to the  "octree neighboor presence table" in the
//thesis of P.Bhattacharya.
//The following function shows wether for a given node its neighbor
//in a particular direction is in the same octant or not.A false
//entrie indicates the presence of neighbor of the given node in the
//same octant as the node is and true otherwise.

//dans les tableau indice 1 = indice du cube , indice 2 = direction
//Pierre utilise des char au lieu de int --> pourquoi ??
static const unsigned int PSI[8][6] = {{4,4,2,2,1,1},
				       {5,5,3,3,0,0},
				       {6,6,0,0,3,3},
				       {7,7,1,1,2,2},
				       {0,0,6,6,5,5},
				       {1,1,7,7,4,4},
				       {2,2,4,4,7,7},
				       {3,3,5,5,6,6}};

// F B U D L R
static const bool NU[8][6] = {{false, true,  false, true,  false, true },
			      {false, true,  false, true,  true,  false},
			      {false, true,  true,  false, false, true },
			      {false, true,  true,  false, true,  false},
			      {true,  false, false, true,  false, true },
			      {true,  false, false, true,  true,  false},
			      {true,  false, true,  false, false, true },
			      {true,  false, true,  false, true,  false}};


CGAL_BEGIN_NAMESPACE

template <class Gt>
class Data_structure_using_octree_3
{
  /* typedef's */
public: 
  typedef Gt                                   Geom_traits;
  typedef typename Gt::FT                      Nb;
  typedef typename Gt::Triangle_3              Triangle;  
  typedef typename Gt::Point_3                 Point;  
  typedef typename Gt::Plane_3                 Plane; 
  typedef typename Gt::Vector_3                Vector;
  typedef typename Gt::Iso_cuboid_3            Iso_Cuboid;
  typedef typename Gt::Segment_3               Segment;
  typedef typename Gt::Line_3                  Line;
  typedef typename Gt::Object_3                Object;
  typedef typename Gt::Ray_3                   Ray;

  //CAUTION : can i use the kernel bbox ?
  typedef CGAL::Bbox_3                         Bbox;

  // Constraints
  typedef Constrained_Element<Gt>              C_elt;
  typedef Constrained_Vertex<Gt>               C_vertex;
  typedef Constrained_Edge<Gt>                 C_edge;
  typedef Constrained_Facet<Gt>                C_facet;

  typedef std::list<C_elt*>                    C_elt_list;
  typedef typename C_elt_list :: iterator      C_elt_iterator;

  typedef std::map<int,C_elt*>                 Constraint_map;
  typedef std::list<int>                       Ind_constraint_list;


  //Name coresponds to the index node. it's a number includ between 0
  //and 7. it's the path from the root to the node. 

  typedef std::list<int>             Name;
  typedef struct                     Node;
  typedef Node*                      Octree;

  // data structure for the octree
  struct Node { 
    Iso_Cuboid                 boite;
    Name                       nom;
    Octree                     fils[8];
    Octree                     pere;
    Ind_constraint_list        i_c_l;
  };

  //a direction
  typedef enum Direction {FRONT  = 0,
        BACK   = 1,
        UP     = 2,
        DOWN   = 3,
        LEFT   = 4,
        RIGHT  = 5};

  //typedef enum Direction {FRONT, BACK, UP, DOWN, LEFT, RIGHT};


public:
  Octree                        arbre;
  Constraint_map                c_m;

private :
  Geom_traits                   gt;

  // for the octree
  double                        niv_max;

  //to keep in memory the visited edge
  std::list<Segment>                         ls;

  //f = number of facet, s = number of constraints
  int count,f,s,nbp;

  double tl,tri;

  //to know which constraint we want to take into account.
  bool    use_vertex, use_edge,  use_facet;

  Bbox    space_bbox;

  Nb PREC;

public:
  // le constructeur
  Data_structure_using_octree_3()
  {
    gt = Geom_traits();
    count=0; 
    f = 0;
    s = 0;
    nbp = 0;
    tl = 0;tri = 0;
    use_vertex = false;
    use_edge   = false;
    use_facet   = true;
    space_bbox = Bbox();
    PREC = 0.0000000000001;
  }

  Data_structure_using_octree_3(bool v, bool e, bool f )
  {
    gt = Geom_traits();
    count=0; 
    f = 0;
    s = 0;nbp = 0;
    tl = 0;tri = 0;
    use_vertex = v;
    use_edge   = e;
    use_facet   = f; 
    PREC = 0.0000000000001;
  }

  void setUse_vertex  (bool b){ use_vertex = b; }
  void setUse_edge    (bool b){ use_edge   = b; }
  void setUse_facet    (bool b){ use_facet   = b; }

public :

  //initialisation
  template <class Point_output_iterator>
  void input(std::ifstream& ifs,
       Point_output_iterator out = Emptyset_iterator());

  void     input      (  C_elt_iterator it, C_elt_iterator end);

  //Queries
  Nb       lfs        ( Point p); 
  Nb       lfs2       ( Point p);        
  Object   intersect  ( Point p, Point q);
  Object   intersect  ( Ray r);    



private :
  void     create_data_structure       (  );
  void     add_constrained_vertex      ( Point p );
  void     add_constrained_edge        ( Point p, Point q);
  void     add_constrained_facet       ( Point p, Point q, Point r);

  void     construct_octree             ( Node* o,
    Ind_constraint_list icl); 

  Constraint_map     intersection      ( Iso_Cuboid boite, Constraint_map* cm); 


  bool     is_leaf                    ( Node* o);
  bool     visited                    ( Segment s);

  //octree traversal function
  std::list<int>     ind_neighbor_son     ( Direction d);
  std::list<Node*>   go_to_neighbor       ( Node* o, Direction d);
  std::list<Node*>   find_face_neighbor   ( Node* noeud,Direction d);
  Node*              go_to_node           ( Octree o, Name p);

  //localisation d'un point dans l'octree
  Node*             locate                ( Octree o,Point p); 


  //Queries
  Nb                lfs_cont              ( Point p, 
                                            Node* noeud,
                                            std::map<int,Nb>* distance,
                                            std::list<int>* ordre,
                                            std::map<Name,bool>*
                                            visited_cube,
                                            Nb min);
  Nb                lfs_simp              ( Point p, 
                                            Node* noeud,
                                            std::map<int,Nb>* distance,
                                            std::list<int>* ordre,
                                            std::map<Name,bool>*
                                            visited_cube,
                                            Nb min);
  int               sci_aux               ( Point p,
                                            Segment s,
                                            Node* noeud, 
                                            std::map<Name,bool>*
                                            visited_cube);      

  //to compute distance
  bool              squared_distance_to_constraint_in_node( Point p,
                                                            Node* noeud,
                                                            std::map<int,Nb>* md,
                                                            std::list<int>* od);

  //to check or find intersection 
  Nb                squared_distance2          ( Point p,Triangle t);

  Nb                squared_distance_to_side  ( Point p,
                                                Iso_Cuboid c,
                                                Direction d);

  bool              does_intersect        ( Segment seg,
                                            Iso_Cuboid c,
                                            Direction d);
  bool              does_intersect        ( Ray ray,
                                            Iso_Cuboid c,
                                            Direction d);

  Point             intersection          ( Segment seg,
                                            Iso_Cuboid c,
                                            Direction d);
  bool              do_intersect2         ( Iso_Cuboid c,
                                            Segment seg);
  Point             intersection          ( Ray ray,
                                            Iso_Cuboid c,
                                            Direction d);
};          

//----------------------------------------------------------------------------

//add a vertex to the constraints list
template <class Gt>
void
Data_structure_using_octree_3 <Gt> :: add_constrained_vertex(Point p )
{
  if (nbp == 0)
    //it's the first point
    space_bbox = p.bbox();
  else
    space_bbox = space_bbox + p.bbox();

  nbp ++;
  if (use_vertex)
  {
    //create a new constraint
    std::pair<int,C_elt*>    c_p;
    C_vertex*                cv = new C_vertex (p );
    c_p.first  = s;
    c_p.second = cv;
    c_m.insert(c_p);
    s++;
  }
}

//----------------------------------------------------------------------------
//return true if s is in ls... make sure that we will a
//constrained_edge is not duplicate
template <class Gt>
bool
Data_structure_using_octree_3 <Gt> ::visited( Segment  s) 
{

  typename Geom_traits::Equal_3 equal_3 =  gt.equal_3_object();

  if (ls.empty())
    return false;
  else {

    typename std::list<Segment>::iterator 
      lsit  = ls.begin(),
      lsend = ls.end();

    while ( lsit != lsend ){
      if (equal(s,(*lsit)))
        return true;

      lsit++; 
    }

    return false;
  }
}

//----------------------------------------------------------------------------
// add a constrained edge to the constraints list
template <class Gt>
void
Data_structure_using_octree_3 <Gt> :: add_constrained_edge(Point p, 
                                                           Point q)
{
  if (use_edge) {
    typename Geom_traits::Construct_segment_3 
      construct_segment =  gt.construct_segment_3_object();
    Segment seg = construct_segment(p,q);

    if (!(visited(seg))){ 
      ls.push_front(seg);

      std::pair<int,C_elt*>    c_p;
      C_edge*                  ce = new C_edge(seg);
      c_p.first  = s;
      c_p.second = ce;
      c_m.insert(c_p);

      s++;

    }
  }

}
//----------------------------------------------------------------------------

// add a constrained_facet to the constraints list
template <class Gt>
void
Data_structure_using_octree_3 <Gt> :: 
add_constrained_facet(Point p, Point q, Point r){

  if (use_facet) {
    std::pair<int,C_elt*>    c_p;
    C_facet* cf = new C_facet (p,q,r);
    c_p.first  = s;
    c_p.second = cf;
    c_m.insert(c_p);


    f++;
    s++;
  }
}
//----------------------------------------------------------------------------
//create the data structure... initialisation + construct the data
//structure via Construct_octree
template <class Gt>
void
Data_structure_using_octree_3 <Gt> :: create_data_structure()

{

  Point 
    p = Point(space_bbox.xmin(),space_bbox.ymin(),space_bbox.zmin()),
    q = Point(space_bbox.xmax(),space_bbox.ymax(),space_bbox.zmax());

  typename Geom_traits::Construct_iso_cuboid_3
    construct_iso_cuboid =  gt.construct_iso_cuboid_3_object();
  typename Geom_traits::Compute_volume_3  compute_volume =
    gt.compute_volume_3_object();

  Iso_Cuboid boite = construct_iso_cuboid(p,q);


  double 
    vol = CGAL::to_double(compute_volume(boite));

  //CAUTION : this criteria have no sens... must find an other one.
  // we want a depth max such that the volume of the cuboid isn't
  // smaller than 10^-3
  niv_max =std::log(vol)/log((double)8) + 3*
    std::log((double)10)/log((double)8);
  int niv = 5;
  // the deph of the tree mustn't be greater than 5 otherwise the
  // tree will be too big ie the time of contruction too long  
  if ( niv_max > niv   )
    niv_max = niv;

  arbre = new Node();
  arbre->nom = Name();
  arbre->boite = boite;
  arbre->i_c_l = std::list<int>::list();
  arbre->pere = NULL;

  Ind_constraint_list icl;

  typename Constraint_map ::iterator
    cmit  = c_m.begin(), 
    cmend = c_m.end();

  while ( cmit != cmend )
  { 
    icl.push_front((*cmit).first);
    ++cmit;
  }

  construct_octree(arbre, icl);
}

//----------------------------------------------------------------------------- 
// construc the octree
template <class Gt>
void 
Data_structure_using_octree_3<Gt> :: construct_octree ( Node* o,
                                                       Ind_constraint_list icl)
{ 

  typename Geom_traits:: Construct_iso_cuboid_3 construct_iso_cuboid =
    gt.construct_iso_cuboid_3_object();
  typename Geom_traits:: Construct_midpoint_3 construct_midpoint =
    gt.construct_midpoint_3_object();
  typename Geom_traits::Construct_centroid_3  construct_centroid =
    gt.construct_centroid_3_object();
  typename Geom_traits::Construct_vertex_3  construct_vertex =
    gt.construct_vertex_3_object();

  typename Geom_traits::Compute_volume_3  compute_volume =
    gt.compute_volume_3_object();

  //threshold has to leave which the cell is not subdivised
  unsigned s = 25; 


  // go through all constraints of the node
  Ind_constraint_list ::iterator
    iclit = icl.begin(), 
    iclend = icl.end();

  //std::cout<<"creation de la liste de contraintes"<<std::endl;

  while ( iclit != iclend )
  {
    if ( c_m[(*iclit)]->does_intersect( o->boite ))
    {
      o->i_c_l.push_front((*iclit));
    }
    ++iclit;
  }

  //std::cout<<"creation de la liste des fils"<<std::endl;
  if ((o->i_c_l.size()>=s) && (o->nom.size()<=niv_max))
  {
    //create the 8 children of the node
    for (int l=0;l<8;l++)
    {
      o->fils[l] = new Node();
      o->fils[l]->nom = std::list<int>::list(o->nom);
      o->fils[l]->nom.push_back(l);
      o->fils[l]->i_c_l = std::list<int>::list(); 
      o->fils[l]->pere = o;
    }

    Point s[8];
    for (int i=0; i<8;i++) {
      s[i] = construct_vertex(o->boite,i);
    }

    Point m01, m03, m05, m27, m67, m47, m2347, m0123, m5647, m0354,
      m0156,m1267, centre;

    m01    = construct_midpoint( s[0],s[1]);
    m03    = construct_midpoint( s[0],s[3]);
    m05    = construct_midpoint( s[0],s[5]);
    m27    = construct_midpoint( s[2],s[7]);
    m67    = construct_midpoint( s[6],s[7]);
    m47    = construct_midpoint( s[4],s[7]);
    m0123  = construct_centroid( s[0],s[1],s[2],s[3]);
    m2347  = construct_centroid( s[4],s[7],s[2],s[3]); 
    m5647  = construct_centroid( s[5],s[6],s[4],s[7]); 
    m0354  = construct_centroid( s[4],s[0],s[5],s[3]); 
    m0156  = construct_centroid( s[0],s[1],s[5],s[6]); 
    m1267  = construct_centroid( s[1],s[2],s[6],s[7]);
    centre = construct_midpoint( m0156,m2347);

    o->fils[0]->boite = construct_iso_cuboid(m05,m5647);
    o->fils[1]->boite = construct_iso_cuboid(m0156,m67);
    o->fils[2]->boite = construct_iso_cuboid(s[0],centre);
    o->fils[3]->boite = construct_iso_cuboid(m01,m1267);
    o->fils[4]->boite = construct_iso_cuboid(m0354,m47);
    o->fils[5]->boite = construct_iso_cuboid(centre,s[7]);
    o->fils[6]->boite = construct_iso_cuboid(m03,m2347);
    o->fils[7]->boite = construct_iso_cuboid(m0123,m27);

    for(int i = 0; i <8; i++){
      //std::cerr << "|\r ";
      construct_octree(o->fils[i], o->i_c_l);
    }

  }
  else{
    for (int l=0;l<8;l++){
      o->fils[l] = new Node();
      o->fils[l] = NULL;
    }
  }

  /*std::cout<<"Noeud cree : ";
  for( Name::iterator it = o->nom.begin(); it != o->nom.end(); it++)
    std::cout<<*it<<"   ";
  std::cout<<std::endl;
  std::cout<<"\tnombre de contraintes : "<<o->i_c_l.size()<<std::endl;*/
}
//----------------------------------------------------------------------------
template <class Gt>
void 
Data_structure_using_octree_3<Gt> ::  
input( C_elt_iterator it, C_elt_iterator end ){

  while (it != end) {
    //compute the bbox
    if (nbp == 0)
      //it's the first point
      space_bbox = (*it)->getBBox();
    else
      space_bbox = space_bbox + (*it)->getBBox();
    nbp ++;

    std::pair<int,C_elt*>    c_p;
    c_p.first  = s;
    c_p.second = (*it);
    c_m.insert(c_p);
    s++;
    it++;
  }

  create_data_structure();



}

//----------------------------------------------------------------------------
//read a oFF file and construct the corresponding data structure
template <class Gt>
template <class Point_output_iterator>
void 
Data_structure_using_octree_3<Gt> ::  
input(std::ifstream& ifs,
      Point_output_iterator out){
  int nb_points,nb_faces, ent2;
  char entete;

  // lecture de off
  entete = ifs.get();
  entete = ifs.get();
  entete = ifs.get();


  // lecture des 3 nombres
  ifs >> nb_points;
  ifs >> nb_faces;
  ifs >> ent2;

  int i1,i2,n,s = 0;
  double x,y,z,f1;

  std::cout<<"nombre de points : "<<nb_points<<std::endl;
  std::cout<<"nombre de faces : "<<nb_faces<<std::endl;
  std::cout<<"ent2 : "<<ent2<<std::endl;

  std::vector<Point> vvect(nb_points);

  //std::list<C_elt*> lc_input;

  // points
  while (s < nb_points){
    ifs >> x; ifs >> y; ifs >> z;
    vvect[s] = Point (x,y,z);
    *out++ = vvect[s];
    //idea : instead of compute the bbox in constrained_vertex why
    //don't you compute it here ?
    add_constrained_vertex(vvect[s]);
    s++;
  }

  std::cout<<"number of vertices : "<<s<<std::endl;

  //linking
  int nb_points_of_face;
  int i = 0,k;
  ifs >> nb_points_of_face;

  Point *p = new Point[nb_points_of_face];

  while (!ifs.eof()) { 
    std::set<int> lind;
    lind.clear();
    ifs >> i1;
    lind.insert(i1);

    p[0] = vvect[i1];

    n = 0;
    k = 1;
    for (i = 2 ; i <= nb_points_of_face ; i++) {
      ifs >> i2;
      // pour eviter d'avoir deux fois le meme sommet ds la face
      if (lind.find(i2) == lind.end()){
        p[k]=vvect[i2];
        n++;
        lind.insert(i2); 
      }
      k++;
    }
    if (n >= 2) {// ruse pour traiter uniquement le cas de faces triangulaires
      add_constrained_facet(p[0],p[1],p[2]);
      //C_facet* cf = new C_facet (p[0],p[1],p[2]);
      //lc_input.push_front(cf);

    }
    ifs >> f1;

    while ((!ifs.eof()) && (f1 < 2)) ifs >> f1;

    nb_points_of_face = (int) f1;
  }

  std::cout<<"number of facets : "<<f<<std::endl;

  //creation de la structure
  create_data_structure();
  std::cout<<"data structure created"<<std::endl;

}


//-----------------------------------------------------------------------------
// find the constraints in cm which intersect the iso_cuboide boite
template <class Gt>
typename Data_structure_using_octree_3<Gt> :: Constraint_map
Data_structure_using_octree_3<Gt> :: intersection (Iso_Cuboid      boite, 
                                                   Constraint_map* cm)
{
  Constraint_map result = *( new  Constraint_map);
  typename Constraint_map::iterator
    cmit = cm->begin(), 
    cmend = cm->end();

  while ( cmit != cmend ){
    if ((*cmit)->does_intersect(boite) ){
      result.push_back((*cmit));
    }
    ++cmit;
  }
  return result;

} 


//-----------------------------------------------------------------------------
//auxiliary function : return the children indices of a block such
//that this sub-block have a side includ in the block's side in
//direction D
template <class Gt>
std :: list<int>
Data_structure_using_octree_3 <Gt> :: ind_neighbor_son (Direction d)
{

  std::list<int> ind ;
  switch (d)
  {
  case BACK :
    ind.push_front(0) ;
    ind.push_front(1) ;
    ind.push_front(2) ;
    ind.push_front(3) ;
    break;
  case FRONT  :
    ind.push_front(4) ;
    ind.push_front(5) ;
    ind.push_front(6) ;
    ind.push_front(7) ;
    break;
  case DOWN : 
    ind.push_front(0) ;
    ind.push_front(1) ;
    ind.push_front(4) ;
    ind.push_front(5) ;
    break;
  case UP   : 
    ind.push_front(2) ;
    ind.push_front(3) ;
    ind.push_front(6) ;
    ind.push_front(7) ;
    break;
  case RIGHT  : 
    ind.push_front(0) ;
    ind.push_front(2) ;
    ind.push_front(4) ;
    ind.push_front(6) ;
    break;
  case LEFT : 
    ind.push_front(1) ;
    ind.push_front(3) ;
    ind.push_front(5) ;
    ind.push_front(7) ;
    break;
  default : //it mustn't be
    CGAL_kernel_assertion_msg(false,"error .... ");
    return ind;
  }

  return ind;
}

//-----------------------------------------------------------------------------
//we found the neighbor node "voisin" of a given node, we want to
//know all the neighbor leafs of this node. So we go down the tree,
//choosing in voisin each child which its corresponding block have a
//side includ in side of voisin's block in direction D.
template <class Gt>
std :: list<typename Data_structure_using_octree_3 <Gt> :: Node*>
Data_structure_using_octree_3 <Gt> :: go_to_neighbor ( Node* voisin,
                                                      Direction d)
{

  //CGAL :: Timer t0;
  //t0.start();

  //std::cerr<<"entree dans gotoneighbor"<<endl;

  std::list<Node*> result, aux;

  if (is_leaf(voisin)){ 
    result.push_front(voisin);
    return result;
  }
  else {
    //do the same with each child.
    std::list<int> ind_fils = ind_neighbor_son (d);
    std::list<int> :: iterator 
      ifit  = ind_fils.begin(),
      ifend = ind_fils.end();

    while (ifit != ifend){
      if (voisin->fils[(*ifit)] != NULL ) {
        if (is_leaf(voisin->fils[(*ifit)])){
          result.push_front(voisin->fils[(*ifit)]);     
        }
        else{
          aux = go_to_neighbor(voisin->fils[(*ifit)],d);
          result.merge (aux);
        }
      }


      ifit++;
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
//return the neighbor leaf of a given node noeud in the diection D.
template <class Gt>
std :: list<typename Data_structure_using_octree_3 <Gt> :: Node*>
Data_structure_using_octree_3 <Gt> :: 
find_face_neighbor (Node* noeud, Direction d)
{

  //curant neighbor node
  Node* voisin = &*noeud;

  //n --> for going down the tree(step 2) keep in memory the path...
  //p --> index of the node
  Name  n, p = noeud->nom;
  n.clear();

  Name :: iterator 
    pbeg = --p.begin(),
    pit  = --p.end();


  while (pit != pbeg)
  {
    if ( NU[(*pit)][d]) //( nu((*pit),d))
    {//compute n for step 2
      n.push_front(PSI[(*pit)][d]);  //psi((*pit),d));

      //go up to the father
      if (voisin->pere != NULL )
        voisin =&*( voisin->pere);
      else 
      {// it mustn't be
        CGAL_kernel_assertion_msg(false,"error : voisin->pere = NULL");
        break;
      }

      --pit;
    }
    else
    {
      //the neighbor is a present in the same octant of voisin
      //voisin->pere is the first common ancestor  of the node and
      //its neighbor

      if ( voisin->pere != NULL )
        voisin = &*( voisin->pere);
      else {
        CGAL_kernel_assertion_msg(false,"error : voisin->pere = NULL");
        //break;
      }

      //we compute the neighbor node in the octant
      n.push_front(PSI[(*pit)][d]);   //psi((*pit),d)); 

      // go down the tree
      Name::iterator
        nit = n.begin(),
        nend = n.end();

      while (nit != nend)
      {
        if ( !(is_leaf(voisin)))
          voisin =&*( voisin->fils[(*nit)]);
        else
          break;
        ++nit;
      }

      pit = pbeg;
    }
  }

  //we found the neighbor node of noeud in the direction D. Now,  we
  //want to know all the neighbor leafs of this node.
  std::list< Node* > result, aux;

  if ( voisin->nom.size() !=0) { 
    if (is_leaf(voisin)) {
      result.push_front(voisin);
      return result;
    }
    else {
      std::list<int> ind_fils = ind_neighbor_son (d);
      std::list<int> :: iterator 
        ifit  = ind_fils.begin(),
        ifend = ind_fils.end();

      while (ifit != ifend)
      {
        if (voisin->fils[(*ifit)] != NULL ) {

          if (is_leaf(voisin->fils[(*ifit)])){
            result.push_front(voisin->fils[(*ifit)]);
          }
          else {
            aux = go_to_neighbor(voisin->fils[(*ifit)],d);
            result.merge (aux);
          }
        }
        ifit++;
      }
    }
  }

  return result;
}
//-----------------------------------------------------------------------------
//return the node correspomding to the name p
template <class Gt>
typename Data_structure_using_octree_3 <Gt> :: Node*
Data_structure_using_octree_3 <Gt> :: go_to_node (Octree o, Name p) {

  if (p.empty()){
    return o;
  }
  else {
    int pn = p.front();
    p.pop_front ();
    if (o->fils[pn] != NULL )
      return  go_to_node (o->fils[pn],p);
    else {
      return o;
    }
  }
} 

//-----------------------------------------------------------------------------
//return true iff o is a leaf
template <class Gt>
bool
Data_structure_using_octree_3 <Gt> :: is_leaf(Node* o) 
{
  for (int i = 0; i<8 ; i++ ) 
  {
    if (o->fils[i] != NULL)
      return false;
  }
  return true;
}

//-----------------------------------------------------------------------------
// locate a point p in the octree
template <class Gt>
typename Data_structure_using_octree_3 <Gt> :: Node*
Data_structure_using_octree_3 <Gt> :: locate(Octree o, Point p)
{
  typename Geom_traits:: Has_on_bounded_side_3  has_on_bounded_side
    =gt.has_on_bounded_side_3_object();
  typename Geom_traits:: Has_on_boundary_3  has_on_boundary
    =gt.has_on_boundary_3_object();

  /*std::cout<<"localisation dans : ";
  for( Name::iterator it = o->nom.begin(); it != o->nom.end(); it++)
    std::cout<<*it<<"   ";
  std::cout<<std::endl;*/
  //std::cout<<o->boite<<std::endl;

  if ( (has_on_bounded_side(o->boite,p))||
    (has_on_boundary(o->boite ,p))) 
  {
      if (is_leaf(o)) {
        return o;
      }
      else {
        //look at the child whose block contain p
        for (int i=0;i<8;i++){
          if ( (has_on_bounded_side(o->fils[i]->boite ,p)) 
            ||  (has_on_boundary(o->fils[i]->boite ,p))){
              return locate(o->fils[i] ,p);
            }
        }
      }
    }
  //else
  //  std::cout<<"is not in bounded side..."<<std::endl;

    //std::cout<<"the point p lie outside the bounding box"<<std::endl;
    return NULL;
}



//-----------------------------------------------------------------------------
template <class Gt>
typename Data_structure_using_octree_3<Gt> :: Nb
Data_structure_using_octree_3<Gt> :: squared_distance2(Point p,Triangle t)
{

  typename Geom_traits::Construct_vertex_3  construct_vertex =
    gt.construct_vertex_3_object(); 

  typename Geom_traits::Construct_segment_3 
    construct_segment =  gt.construct_segment_3_object();

  typename Geom_traits:: Compute_squared_distance_3
    compute_squared_distance = gt.compute_squared_distance_3_object();

  Nb res,restemp;
  Point 
    a0 = construct_vertex(t,0),
    a1 = construct_vertex(t,1),
    a2 = construct_vertex(t,2);

  Vector v0,v1,v2;
  v0 = a2 - a1;
  v1 = a0 - a2; 
  v2 = a1 - a0; 

  //normal of the plane a0a1a2
  Vector N = cross_product(v0,v1);

  //let p' the projected point of p on the plane a0a1a2

  bool b0=0,b1=0,b2=0;
  //test if p' is anticlockwise of v0
  if (cross_product(v0,p - a0)*N >0)
    b0 = 1;
  //test if p' is anticlockwise of v1
  if (cross_product(v1,p - a1)*N >0)
    b0 = 1;
  //test if p' is anticlockwise of v2
  if (cross_product(v2,p - a2)*N >0)
    b0 = 1;

  if (b0*b1*b2) {
    res = (N*(p - a1)*
      (N*(p - a1)/(N*N)));
  }
  else{
    Segment
      e1 = construct_segment (a0,a1),
      e2 = construct_segment (a0,a2),
      e3 = construct_segment (a1,a2);

    res = squared_distance(p,e1);
    restemp = squared_distance(p,e2);

    if (res > restemp) res = restemp;

    restemp = squared_distance(p,e3);
    if (res > restemp) res = restemp; 
  }
  return res;

}
//----------------------------------------------------------------------------
//compute the distance between a point and the side of c in direction d
template <class Gt>
typename Data_structure_using_octree_3<Gt> :: Nb
Data_structure_using_octree_3<Gt> ::squared_distance_to_side( Point p,
                                                             Iso_Cuboid c, 
                                                             Direction d)
{
  typename Geom_traits::Construct_plane_3  construct_plane =
    gt.construct_plane_3_object();
  typename Geom_traits::Construct_vertex_3  construct_vertex =
    gt.construct_vertex_3_object();

  Point q,r,s,t;
  switch (d){
case FRONT :
  t = construct_vertex(c,0);
  q = construct_vertex(c,1);
  r = construct_vertex(c,6);
  s = construct_vertex(c,5);
  break;
case BACK  :
  t = construct_vertex(c,2);
  q = construct_vertex(c,3);
  r = construct_vertex(c,4);
  s = construct_vertex(c,7);
  break;
case UP    : 
  t = construct_vertex(c,6);
  q = construct_vertex(c,5);
  r = construct_vertex(c,4);
  s = construct_vertex(c,7);
  break;
case DOWN  : 
  t = construct_vertex(c,0);
  q = construct_vertex(c,1);
  r = construct_vertex(c,2);
  s = construct_vertex(c,3);
  break;
case LEFT  : 
  t = construct_vertex(c,0);
  q = construct_vertex(c,3);
  r = construct_vertex(c,4);
  s = construct_vertex(c,5);
  break;
case RIGHT : 
  t = construct_vertex(c,1);
  q = construct_vertex(c,2);
  r = construct_vertex(c,7);
  s = construct_vertex(c,6);
  break;
default : //it mustn't be
  CGAL_kernel_assertion_msg(false,"error .... ");
  return -4;
  }

  typename Geom_traits::Construct_triangle_3
    construct_triangle =  gt.construct_triangle_3_object();
  Triangle
    t1 = construct_triangle(t,q,s),
    t2 = construct_triangle(q,r,s);

  Nb 
    d1 = squared_distance2(p,t1),
    d2 = squared_distance2(p,t2);

  if (d1<d2) return d1;
  else return d2;

}
//-----------------------------------------------------------------------------
//compute the distance between a point p and all constraints whose
//intersect the node noeu.
//input  : a point p and a node noeud
//output : the map md such that md[i]= the distance between p and the
//ith constraint
//         the list od, well-ordered such that for i<j, md[od[i]] is
//         nearest p than md[od[j]]
template <class Gt>
bool 
Data_structure_using_octree_3<Gt> :: 
squared_distance_to_constraint_in_node(Point p,
                                       Node* noeud,
                                       std::map<int,Nb>* md,
                                       std::list<int>* od)
{  

  //CGAL :: Timer t0;
  //t0.start();

  Nb                   dist;
  std::pair<int,Nb>   d_aux;

  if ( noeud != NULL ){
    //compute the distance between p and all the constraints of noeud.
    Ind_constraint_list :: iterator
      iclit  = noeud->i_c_l.begin(),
      iclend = noeud->i_c_l.end();

    while (iclit != iclend){
      if ( md->find((*iclit)) == md->end() ) {
        // Compute the distance between p and the constraint (*iclit)
        dist = c_m[(*iclit)]->sqared_distance_to_sheet(p);
        d_aux.first  = (*iclit);
        d_aux.second = dist;
        md->insert(d_aux);

        if (od->empty()) // it's the first constraint
          od->push_front((*iclit));
        else { // find where we have to put the constraint in od
          std::list<int> :: iterator 
            it  = od->begin(),
            end = od->end();

          while ((it != end) && ( (*md)[(*it)]<dist))
            it++;

          od->insert(it,(*iclit));
        }
      }      
      iclit++;
    }
    return true;
  }
  return false;

}


//-----------------------------------------------------------------------------
//lfs auxiliary function... compute the lfs with shewchuk definition
template <class Gt>
typename Data_structure_using_octree_3<Gt> :: Nb
Data_structure_using_octree_3<Gt> :: 
lfs_cont( Point p, 
         Node* noeud,
         std::map<int,Nb>* distance,
         std::list<int>* ordre,
         std::map<Name,bool>* visited_cube,
         typename Data_structure_using_octree_3<Gt> :: Nb min)
{

  std :: list<Node*>   voisins;

  (*visited_cube)[noeud->nom] = true;

  if (noeud->i_c_l.size()!= 0){
    //compute the distance between p and all constraints of the node
    //whose containt it
    squared_distance_to_constraint_in_node(p,noeud,distance,ordre);

    //we go through ordre until we find 2 separates constraints 
    //nearest p
    std::list<int> :: iterator 
      oit  = ordre->begin(),
      oend = ordre->end();

    bool encours = true;
    int i =0;
    while ((oit != oend) && (encours) )
    {
      std::list<int> :: iterator oit2  = ordre->begin();
      int j = 0;

      while ( (oit2 != oit) && (encours) )
      {
        std::cout<<" couple "<<i<<", "<<j<<std::endl;

        if ( !( c_m[(*oit2)]->does_intersect( c_m[(*oit)])) )
        {
          Nb d = (*distance)[(*oit2)];

          std::cout<<" d = "<<d<<std::endl;

          if ( d>PREC ) {
            //min = d;
            min = (*distance)[(*oit)];
            encours = false;
          }
        }
        else
          std::cout<<"(c_m[(*oit2)]->does_intersect( c_m[(*oit)]))"<<std::endl;


        oit2++;
        j++;
      }
      oit++;
      i++;
      std::cout<<"\tmin = "<<min<<std::endl;
    }
  }

  //visiting neighbor in each direction
  Direction tdir[6] = {FRONT, BACK, UP, DOWN, LEFT, RIGHT};
  Nb l;

  for (int i=0; i<6;i++)  {
    Nb dc =  squared_distance_to_side (p ,noeud->boite,tdir[i]);

    if ( (min == -1) || ( dc < min ) ){
      voisins = find_face_neighbor(noeud,tdir[i]);
      typename std::list<Node*> :: iterator
        vit  = voisins.begin(),
        vend = voisins.end();

      while (vit != vend){

        if (!((*visited_cube)[(*vit)->nom])){
          l=lfs_cont(p,(*vit),distance,ordre,visited_cube,min);
          if ( (min == -1) || (l < min))
            min = l;
        }
        vit++;
      }
    }
  }
  return min;
}
//-----------------------------------------------------------------------------
//lfs auxiliary function... compute the lfs with the definition 2
template <class Gt>
typename Data_structure_using_octree_3<Gt> ::Nb   
Data_structure_using_octree_3<Gt> :: lfs_simp( Point p, 
                                              Node* noeud,
                                              std::map<int,Nb>* distance,
                                              std::list<int>* ordre,
                                              std::map<Name,bool>*
                                              visited_cube,
                                              Nb min)
{
  std :: list<Node*>   voisins;

  (*visited_cube)[noeud->nom] = true;

  if (noeud->i_c_l.size()!= 0){
    //compute the distance between p and all constraints of the node
    //whose containt it
    squared_distance_to_constraint_in_node(p,noeud,distance,ordre);

    //we go through ordre until we find a constraint such that its
    //distance to p is greater than PREC 
    //ie find the nearest constraint (of noeud) whose not containt p.
    std::list<int> :: iterator 
      oit  = ordre->begin(),
      oend = ordre->end();

    bool encours = true;
    while ((oit != oend) && (encours) ){
      Nb d = (*distance)[(*oit)];
      if (d>PREC) { 
        min = d;
        encours = false;
      }
      oit++;
    }

  } 

  //visiting neighbor in each direction
  Direction tdir[6] = {FRONT, BACK, UP, DOWN, LEFT, RIGHT};
  Nb l;

  //std::cout<<"visite des voisins"<<std::endl;
  for (int i=0; i<6;i++)  {
    Nb dc =  squared_distance_to_side (p ,noeud->boite,tdir[i]);

    if ( (min == -1) || ( dc < min ) ){

      voisins = find_face_neighbor(noeud,tdir[i]);

      typename std::list<Node*> :: iterator
        vit  = voisins.begin(),
        vend = voisins.end();

      while (vit != vend){

        if (!((*visited_cube)[(*vit)->nom])) {
          l=lfs_simp(p,(*vit),distance,ordre,visited_cube,min);

          if ( (min == -1) || (l < min))
            min = l;
        }
        vit++;
      }
    }
  }
  return min;
}
//-----------------------------------------------------------------------------
//Compute the lfs with shewchuk definition
template <class Gt>
typename Data_structure_using_octree_3<Gt> :: Nb
Data_structure_using_octree_3<Gt> :: lfs(Point p){


  //   lfs_time.start();

  std::list<int>       ordre ;
  std::map<int,Nb>     distance;
  std::map<Name,bool> visited_cube;

  //locate the point p in the octree

  Node* noeud = locate(arbre,p);

  //   CGAL::Timer t0;
  //   std::cout<<count++<<" - point traite : "<<p<<endl;

  //   t0.start();
  Nb l=-3;
  if (noeud != NULL){
    l= lfs_cont(p,noeud,&distance,&ordre,&visited_cube,-1 );
  }
  else
    std::cout<<"Noeud == NULL"<<std::endl;

  //   t0.stop();

  //std::cout<<"lfs calculee en "<<t0.time()<<endl;
  

  //   tl = tl +t0.time();

  //   lfs_time.stop();

  return l;
}
//-----------------------------------------------------------------------------
//Compute the lfs with the definition 2
template <class Gt>
typename Data_structure_using_octree_3<Gt> :: Nb
Data_structure_using_octree_3<Gt> :: lfs2(Point p)
{


  //   lfs_time.start();
  //std::cout<<"Point : "<<p<<std::endl;

  std::list<int>       ordre ;
  std::map<int,Nb>     distance;
  std::map<Name,bool> visited_cube;

  //locate the point p in the octree
  Node* noeud = locate(arbre,p);

  //   CGAL::Timer t0;
  //   std::cout<<count++<<" - point traite : "<<p<<endl;
  //   t0.start();

  Nb l=-3;
  if (noeud != NULL)
  {
    l= lfs_simp(p,noeud,&distance,&ordre,&visited_cube,-1 );
  }

  //   t0.stop();

  //std::cout<<"lfs calculee en "<<t0.time()<<endl;
  //std::cout<<"lfs de ce point = "<<l<<endl<<endl;
  //std::cout<<tl<<endl;

  //   tl = tl +t0.time();

  //   lfs_time.stop();

  return l;
}
//----------------------------------------------------------------------------
//Return true iff seg intersect the side of c in the direction d
template <class Gt>
bool
Data_structure_using_octree_3<Gt> :: does_intersect( Segment seg,
                                                    Iso_Cuboid c,
                                                    Direction d)
{

  typename Geom_traits::Construct_vertex_3  construct_vertex =
    gt.construct_vertex_3_object();
  typename Geom_traits:: Is_degenerate_3  is_degenerate
    =gt.is_degenerate_3_object();
  typename Geom_traits::Construct_triangle_3
    construct_triangle =  gt.construct_triangle_3_object();

  if (is_degenerate(seg))
    return false;
  if (is_degenerate(c))
    return false;

  //q, r, s and t are the vertex of the side of c in the direction d
  Point q,r,s,t;
  switch (d)
  {
  case FRONT :
    t = construct_vertex(c,0);
    q = construct_vertex(c,1);
    r = construct_vertex(c,6);
    s = construct_vertex(c,5);
    break;
  case BACK  :
    t = construct_vertex(c,3);
    q = construct_vertex(c,2);
    r = construct_vertex(c,7);
    s = construct_vertex(c,4);
    break;
  case UP    : 
    t = construct_vertex(c,6);
    q = construct_vertex(c,5);
    r = construct_vertex(c,4);
    s = construct_vertex(c,7);
    break;
  case DOWN  : 
    t = construct_vertex(c,1);
    q = construct_vertex(c,0);
    r = construct_vertex(c,3);
    s = construct_vertex(c,2);
    break;
  case LEFT  : 
    t = construct_vertex(c,0);
    q = construct_vertex(c,3);
    r = construct_vertex(c,4);
    s = construct_vertex(c,5);
    break;
  case RIGHT : 
    t = construct_vertex(c,1);
    q = construct_vertex(c,2);
    r = construct_vertex(c,7);
    s = construct_vertex(c,6);
    break;
  default : //ca ne doit pas arriver !!
    CGAL_kernel_assertion_msg(false,"error .... ");
    return -4;
  }

  Triangle
    t1 = construct_triangle(t,q,r),
    t2 = construct_triangle(t,r,s);  

  return ( (do_intersect (t1,seg)) 
    || (do_intersect (t2,seg)) );

}
//----------------------------------------------------------------------------
//Return true iff ray intersect the side of c in the direction d
template <class Gt>
bool
Data_structure_using_octree_3<Gt> :: does_intersect( Ray ray,
                                                    Iso_Cuboid c,
                                                    Direction d)
{

  typename Geom_traits::Construct_plane_3  construct_plane =
    gt.construct_plane_3_object();
  typename Geom_traits::Construct_vertex_3  construct_vertex =
    gt.construct_vertex_3_object();
  typename Geom_traits:: Is_degenerate_3  is_degenerate
    =gt.is_degenerate_3_object();


  if (is_degenerate(ray))
    return false;
  if (is_degenerate(c))
    return false;


  Point q,r,s,t;
  switch (d)
  {
  case FRONT :
    t = construct_vertex(c,0);
    q = construct_vertex(c,1);
    r = construct_vertex(c,6);
    s = construct_vertex(c,5);
    break;
  case BACK  :
    t = construct_vertex(c,2);
    q = construct_vertex(c,3);
    r = construct_vertex(c,4);
    s = construct_vertex(c,7);
    break;
  case UP    : 
    t = construct_vertex(c,6);
    q = construct_vertex(c,5);
    r = construct_vertex(c,4);
    s = construct_vertex(c,7);
    break;
  case DOWN  : 
    t = construct_vertex(c,0);
    q = construct_vertex(c,1);
    r = construct_vertex(c,2);
    s = construct_vertex(c,3);
    break;
  case LEFT  : 
    t = construct_vertex(c,0);
    q = construct_vertex(c,3);
    r = construct_vertex(c,4);
    s = construct_vertex(c,5);
    break;
  case RIGHT : 
    t = construct_vertex(c,1);
    q = construct_vertex(c,2);
    r = construct_vertex(c,7);
    s = construct_vertex(c,6);
    break;
  default : //ca ne doit pas arriver !!
    CGAL_kernel_assertion_msg(false,"error .... ");
    return -4;
  }

  typename Geom_traits::Construct_triangle_3 
    construct_triangle =  gt.construct_triangle_3_object();
  Triangle  
    t1 = construct_triangle(q,s,r),
    t2 = construct_triangle(q,s,t);  

  return ( (do_intersect(t1,ray)) || 
    (do_intersect(t2,ray)) );
}

//----------------------------------------------------------------------------
//Construct the intersection between seg and the side of c in the
//direction d... 
//PRECONDITION:  seg intersect the side of c in the direction d
template <class Gt>
typename Data_structure_using_octree_3<Gt> :: Point
Data_structure_using_octree_3<Gt> :: intersection ( Segment seg,
                                                   Iso_Cuboid c,
                                                   Direction d)
{

  typename Geom_traits::Construct_plane_3  construct_plane =
    gt.construct_plane_3_object();
  typename Geom_traits::Construct_vertex_3  construct_vertex =
    gt.construct_vertex_3_object();
  typename Geom_traits::Intersect_3  intersect =
    gt.intersect_3_object();
  typename Geom_traits:: Assign_3 assign =
    gt.assign_3_object();

  typename Geom_traits:: Is_degenerate_3  is_degenerate
    =gt.is_degenerate_3_object();


  CGAL_kernel_exactness_precondition( ! c.is_degenerate() ) ;
  CGAL_kernel_exactness_precondition( ! seg.is_degenerate() ) ;

  Point q,r,s,t;
  switch (d){
case FRONT :
  t = construct_vertex(c,0); 
  q = construct_vertex(c,1);
  r = construct_vertex(c,6);
  s = construct_vertex(c,5);
  break;
case BACK  :
  t = construct_vertex(c,2);
  q = construct_vertex(c,3);
  r = construct_vertex(c,4);
  s = construct_vertex(c,7);
  break;
case UP    : 
  t = construct_vertex(c,6);
  q = construct_vertex(c,5);
  r = construct_vertex(c,4);
  s = construct_vertex(c,7);
  break;
case DOWN  :
  t = construct_vertex(c,0);
  q = construct_vertex(c,1);
  r = construct_vertex(c,2);
  s = construct_vertex(c,3);
  break;
case LEFT  : 
  t = construct_vertex(c,0);
  q = construct_vertex(c,3);
  r = construct_vertex(c,4);
  s = construct_vertex(c,5);
  break;
case RIGHT :
  t = construct_vertex(c,1);
  q = construct_vertex(c,2);
  r = construct_vertex(c,7);
  s = construct_vertex(c,6);
  break;
default : //ca ne doit pas arriver !!
  CGAL_kernel_assertion_msg(false,"error .... ");
  //return NULL;
  }


  Plane plan = construct_plane(s,q,r);
  typename Gt::Object_3 o = intersect(plan,seg);

  Point p ;
  assign(p,o);

  return p;

}

//----------------------------------------------------------------------------
//return true iff seg intersect the iso_cuboid c
template <class Gt>
bool
Data_structure_using_octree_3<Gt> :: do_intersect2(Iso_Cuboid c, Segment seg)  
{
  typename Geom_traits::Construct_vertex_3 
    construct_vertex =  gt.construct_vertex_3_object();
  typename Geom_traits:: Has_on_bounded_side_3  has_on_bounded_side
    =gt.has_on_bounded_side_3_object();
  typename Geom_traits::Construct_triangle_3
    construct_triangle =  gt.construct_triangle_3_object();
  typename Geom_traits:: Has_on_boundary_3  has_on_boundary
    =gt.has_on_boundary_3_object();

  typename Geom_traits:: Is_degenerate_3  is_degenerate
    =gt.is_degenerate_3_object();


  if (is_degenerate(seg))
    return false;
  if (is_degenerate(c))
    return false;

  Point 
    p = construct_vertex(seg,0),
    q = construct_vertex(seg,1);

  //test if one of segment vertex  lie on the bounded side or on the
  // boundary of c
  if (has_on_bounded_side(c,p) ||  has_on_boundary(c,q) ) {
    return true;
  }
  if (has_on_bounded_side(c,q) ||  has_on_boundary(c,q) ){   
    return true;
  }

  //neither of segment vertex lie on the bounded side or on the
  // boundary of c ---> construct triangle whose make the sides of c
  // and test if segment intersect these triangles.
  Point sc[8];
  for (int i=0; i<8;i++){
    sc[i] = construct_vertex(c,i);
    //std::cout << i<<" - "<<sc[i]<<endl; 
  }

  Triangle tf;

  tf = construct_triangle(sc[1],sc[3],sc[0]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[1],sc[3],sc[2]);
  if (do_intersect(tf,seg)) return true;

  tf = construct_triangle(sc[1],sc[6],sc[0]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[0],sc[6],sc[5]);
  if (do_intersect(tf,seg)) return true;

  tf = construct_triangle(sc[1],sc[7],sc[6]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[1],sc[7],sc[2]);
  if (do_intersect(tf,seg)) return true;

  tf = construct_triangle(sc[2],sc[4],sc[3]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[2],sc[4],sc[7]);
  if (do_intersect(tf,seg)) return true;     

  tf = construct_triangle(sc[3],sc[5],sc[0]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[3],sc[5],sc[4]);
  if (do_intersect(tf,seg)) return true; 

  tf = construct_triangle(sc[5],sc[4],sc[6]);  
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[6],sc[4],sc[7]);
  if (do_intersect(tf,seg)) return true;

  return false;
}
//----------------------------------------------------------------------------
//Construct the intersection between ray and the side of c in the
//direction d... 
//PRECONDITION:  seg intersect the side of c in the direction d
template <class Gt>
typename Data_structure_using_octree_3<Gt> :: Point
Data_structure_using_octree_3<Gt> :: intersection (Ray ray,
                                                   Iso_Cuboid c,
                                                   Direction d)
{

  typename Geom_traits::Construct_plane_3  construct_plane =
    gt.construct_plane_3_object();
  typename Geom_traits::Construct_vertex_3  construct_vertex =
    gt.construct_vertex_3_object();
  typename Geom_traits::Intersect_3  intersect =
    gt.intersect_3_object();
  typename Geom_traits:: Assign_3 assign =
    gt.assign_3_object();

  typename Geom_traits:: Is_degenerate_3  is_degenerate
    =gt.is_degenerate_3_object();


  CGAL_kernel_exactness_precondition( ! c.is_degenerate() ) ;
  CGAL_kernel_exactness_precondition( ! ray.is_degenerate() ) ;

  Point q,r,s,t;
  switch (d)
  {
  case FRONT :
    //t = construct_vertex(c,0); 
    q = construct_vertex(c,1);
    r = construct_vertex(c,6);
    s = construct_vertex(c,5);
    break;
  case BACK  :
    //t = construct_vertex(c,2);
    q = construct_vertex(c,3);
    r = construct_vertex(c,4);
    s = construct_vertex(c,7);
    break;
  case UP    : 
    //t = construct_vertex(c,6);
    q = construct_vertex(c,5);
    r = construct_vertex(c,4);
    s = construct_vertex(c,7);
    break;
  case DOWN  :
    //t = construct_vertex(c,0);
    q = construct_vertex(c,1);
    r = construct_vertex(c,2);
    s = construct_vertex(c,3);
    break;
  case LEFT  : 
    //t = construct_vertex(c,0);
    q = construct_vertex(c,3);
    r = construct_vertex(c,4);
    s = construct_vertex(c,5);
    break;
  case RIGHT :
    //t = construct_vertex(c,1);
    q = construct_vertex(c,2);
    r = construct_vertex(c,7);
    s = construct_vertex(c,6);
    break;
  default : //ca ne doit pas arriver !!
    CGAL_kernel_assertion_msg(false,"error .... ");
    //return NULL;
  }


  Plane plan = construct_plane(q,r,s);
  typename Gt::Object_3 o = intersect(plan,ray);

  Point p ;
  assign(p,o);

  return p;

}


//----------------------------------------------------------------------------
//Check if there is a constraint whose intersect the segment pq... if
//this intersection exist then compute it.
template <class Gt>
typename Data_structure_using_octree_3<Gt> :: Object
Data_structure_using_octree_3<Gt> :: intersect( Point p, 
                                               Point q)
{
  //   ri_time.start();

  typename Geom_traits::Construct_segment_3  construct_segment =
    gt.construct_segment_3_object();

  typename Geom_traits:: Is_degenerate_3  is_degenerate
    =gt.is_degenerate_3_object();

  std::map<Name,bool> visited_cube;

  //locate the point p in the octree
  Node* noeud = locate(arbre,p);
  int ind_fc = -1;

  Segment s = construct_segment(p,q);

  if (is_degenerate(s))
    return Object();


  if (noeud != NULL) {
    ind_fc = sci_aux(p,s,noeud,&visited_cube);
  }
  else{
    //p doesn't lie in the bounding box so locate q in the octree
    noeud = locate(arbre,q);
    if (noeud != NULL){
      ind_fc = sci_aux(q,s,noeud,&visited_cube); 
    }
    else{
      //p and q don't lie in the bounding box so chek if the segment
      //pq intersect a side of the bounding box, if an intersection
      //exist, compute it and locate it
      Direction tdir[6]  = {FRONT, BACK, UP, DOWN, LEFT, RIGHT};

      int i= 0;
      bool find = false;

      while ((i<6) && (!(find))){

        if (does_intersect(s,arbre->boite,tdir[i])) {
          //we found the side of the bbox s intersect.
          find = true;

          //locate the intersection point
          Point r = intersection(s,arbre->boite,tdir[i]); 
          noeud = locate(arbre,r);
          if (noeud != NULL) {
            ind_fc = sci_aux(r,s,noeud,&visited_cube);
          }
        }
        i++;
      }
    }
  }

  // ri_time.stop();

  if (ind_fc == -1) {
    //we didn't find a constraint intersected by s
    //std::cout<<"false"<<endl<<endl;
    return Object();
  }
  else {
    //std::cout<<"true"<<endl<<endl;
    return  c_m[ind_fc]->intersection(s);
  }


}

//----------------------------------------------------------------------------
// auxiliary function of intersect
template <class Gt>
int
Data_structure_using_octree_3<Gt> :: sci_aux (Point p,
                                              Segment s,
                                              Node* noeud, 
                                              std::map<Name,bool>*
                                              visited_cube)
{

  int ind_min = -1;
  std :: list<Node*>   voisins;

  //mark the cube --> we don't want to viste it again
  (*visited_cube)[noeud->nom] = true;


  //affichage
  // Name :: iterator 
  //     nomit = noeud->nom.begin(),
  //     nomend = noeud->nom.end();

  //   while (nomit != nomend)
  //     {
  //       std::cerr<<(*nomit);
  //       nomit++;
  //     }
  //   std::cerr<<endl;
  //std::cerr<<"boite correspondante "<<noeud->boite<<endl;
  //std::cerr<<"nombre de contraintes du noeud : "<<noeud->i_c_l.size()<<endl;
  //fin affichage

  Ind_constraint_list :: iterator
    iclit  = noeud->i_c_l.begin(),
    iclend = noeud->i_c_l.end();

  bool find = false;  

  while ((iclit != iclend) && (!(find)) )
  {
    if (c_m[(*iclit)]->does_intersect(s)) {
      ind_min = (*iclit); 
      //stop when the first constaint is found
      find = true;
    }   
    iclit++;   
  }

  if (!(find)) {
    //we didn't find constraint whose intersect the segment, so we
    //look at the neighbor
    Direction tdir[6]  = {FRONT, BACK,  UP, DOWN, LEFT,  RIGHT};

    int i = 0;

    while ((i<6) && (!(find))) {
      if (does_intersect(s,noeud->boite,tdir[i])) {

        voisins = find_face_neighbor(noeud,tdir[i]);

        typename std::list<Node*> :: iterator
          vit  = voisins.begin(),
          vend = voisins.end();

        while ((vit != vend) && (!(find))){
          if (!((*visited_cube)[(*vit)->nom])) {

            if (do_intersect2((*vit)->boite,s)) { 
              ind_min =  sci_aux(p,s,(*vit),visited_cube );
              if (ind_min != -1 )
                find = true;
            }  
          } 
          vit++;
        }
      }
      i++;
    }
  } 

  return ind_min;

}
//----------------------------------------------------------------------------
//Check if there is a constraint whose intersect ray... if
//this intersection exist then compute it.
template <class Gt>
typename Data_structure_using_octree_3<Gt> :: Object
Data_structure_using_octree_3<Gt> :: intersect(Ray r)
{

  // ri_time.start();

  typename Geom_traits::Construct_segment_3  construct_segment =
    gt.construct_segment_3_object();
  typename Geom_traits::Construct_point_on_3  construct_point_on =
    gt.construct_point_on_3_object();

  std::map<Name,bool> visited_cube;

  std::set<Point, typename Geom_traits::Compare_xyz_3> lp;

  //find the side of the bbox is intersected by ray
  Direction tdir[6]  = {FRONT, BACK, UP, DOWN, LEFT, RIGHT};
  for (int i=0;i<6;i++){
    if (does_intersect(r,arbre->boite,tdir[i])) 
    {
      Point p = intersection(r,arbre->boite,tdir[i]);
      lp.insert(p);
      //std::cout<<"une intersection trouvee "<<p<<std::endl;
    }
  }

  Point origin = construct_point_on(r, 0);

  Point p,q;

  if (lp.size()==2) {
    p = *(lp.begin());
    q = *(--lp.end());
    //std::cout<<"deux intersections "<<p<<" et "<<q<<std::endl;
    return intersect(p,q);
  }
  else if (lp.size()==1) {
    p = *(lp.begin());
    //std::cout<<"une intersection "<< p <<std::endl;
    return intersect(p,origin);
  }
  else
  {  
    //std::cout<<"pas d'intersection"<<std::endl;
    return Object(); 
  }
}



CGAL_END_NAMESPACE
#endif
