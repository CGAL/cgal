#ifndef CGAL_TRIANGULATION_DS_FACE_2_H
#define CGAL_TRIANGULATION_DS_FACE_2_H

#include <CGAL/Triangulation_short_names_2.h>

template <class Vb, class Fb >
class  CGAL_Triangulation_ds_vertex_2 ;



template < class Vb, class Fb >
class  CGAL_Triangulation_ds_face_2
  : public Fb
{
public:
  //typedef typename Fb::Triangle Triangle;
  typedef CGAL_Triangulation_ds_vertex_2<Vb,Fb> Vertex;
  typedef CGAL_Triangulation_ds_face_2<Vb,Fb> Face;

  // creators
  CGAL_Triangulation_ds_face_2()
    : Fb()
  {}
    
  CGAL_Triangulation_ds_face_2(Vertex* v0, Vertex* v1, Vertex* v2)
    :  Fb(v0,v1,v2)
  {}
    
  CGAL_Triangulation_ds_face_2(Vertex* v0, Vertex* v1, Vertex* v2,
				Face* n0, Face* n1, Face* n2)
    :  Fb(v0,v1,v2,n0,n1,n2)
  {}

   CGAL_Triangulation_ds_face_2( const Face * f)
    :  Fb()
  {
    set_vertices(f->vertex(0), f->vertex(1), f->vertex(2));
    set_neighbors(f->neighbor(0), f->neighbor(1), f->neighbor(2));
  }

  //setting
  inline 
  void set_vertex(int i, Vertex* v)
  {
    Fb::set_vertex(i,v);
  }
    
    
  inline 
   void set_neighbor(int i, Face* n)
  {
    Fb::set_neighbor(i,n);
  }

  //  void set_vertices() inherited
      
  inline 
  void set_vertices(Vertex* v0,
		    Vertex* v1,
		    Vertex* v2)
  {
    Fb::set_vertices(v0,v1,v2);
   }
    
  //   void set_neighbors() inherited
     
  inline
  void set_neighbors(Face* n0,
		     Face* n1,
		     Face* n2)
  {
    Fb::set_neighbors(n0,n1,n2);
  }

  //Vertex Access Member Functions
  Vertex* vertex(int i) const
  {
    return( (Vertex*) (Fb::vertex(i)));
  } 

 inline 
 bool has_vertex(const Vertex* v) const
  {
    return (Fb::has_vertex(v));
  }
    
    
  inline 
  bool has_vertex(const Vertex* v, int& i) const
  {
    return (Fb::has_vertex(v,i));
  }
    
  inline 
  int index(const Vertex* v) const
  {
    return(Fb::vertex_index(v));
  }

  // Neighbors Access Functions
  inline 
  Face* neighbor(int i) const
  {
    return ((Face*) Fb::neighbor(i));
  }
    
  inline 
  bool has_neighbor(const Face* n) const
  {
    return (Fb::has_neighbor(n));
  }
    
    
  inline 
  bool has_neighbor(const Face* n, int& i) const
  {
    return (Fb::has_neighbor(n,i));
  }
    
    
  inline 
  int index(const Face* n) const
  {
    return(Fb::face_index(n));
  }
    
  //Miscelleanous
  int dimension()
  {
    if (vertex(2) != NULL) {return 2;}
    else return( vertex(1) != NULL ? 1 : 0);
  }



  //Additionnal Operations

  //the following function has been moved to the tds class
//   void insert_in_face(Vertex*& v)
//   {
//     CGAL_triangulation_precondition(v != NULL);
//     Vertex* v0 = vertex(0);
//     Vertex* v2 = vertex(2);
//     Vertex* v1 = vertex(1);
//     
//     Face* f1 = neighbor(1);
//     Face* f2 = neighbor(2);
//     int i1,i2 ;
//     if (f1 != NULL) {i1= cw(f1->index(vertex(cw(1))));}
//     if (f2 != NULL) {i2 = cw(f2->index(vertex(cw(2))));}
//     
//     Face* n1ptr = new Face(v0, v, v2,
// 			   this, f1, NULL);
//     
//     Face* n2ptr = new Face(v0, v1, v,
// 			   this, NULL, f2);
// 
//     n1ptr->set_neighbor(2, n2ptr);
//     n2ptr->set_neighbor(1, n1ptr);
//     if (f1 != NULL) {f1->set_neighbor(i1,n1ptr);}
//     if (f2 != NULL) {f2->set_neighbor(i2,n2ptr);}
// 
//     set_vertex(0, v);
//     set_neighbor(1, n1ptr);
//     set_neighbor(2, n2ptr);
// 
//     if( v0->face() == this  ) {  v0->set_face(n2ptr); }
//     v->set_face(this);
//    }
    
//the following function has been moved to the tds class
//   void insert_on_edge(const Vertex* v, int i)
//   {
//     
//       CGAL_triangulation_precondition(v != NULL); 
//       Face* n = neighbor(i);
// 
//       // The following seems natural, but it may fail if the faces
//       // this and n are neighbors on two edges (1-dim triangulation,
//       // with infinite faces
//       // int in = n->index(this);
// 
//       int in;
//       CGAL_triangulation_assertion( n->has_vertex(vertex(cw(i),in)));
//       in = cw(in);
//       insert_in_face(v);
//       n->flip(in); 
//   }
// 

  //no longer necessary
//   bool insert_outside(const Vertex* v, int i)
//   {
//     CGAL_triangulation_precondition(v != NULL);
//         if(neighbor(i) != NULL) {
//             cerr << "Insert_outside impossible as neighbor face already exists" << endl;
//             return false;
//         }
//         Face* f = new Face(v, vertex(cw(i)), vertex(ccw(i)),
//                            this, NULL, NULL);
//         set_neighbor(i, f);
//         v->set_face(f);
//         Vertex* w = vertex(ccw(i));
//         if(w != NULL){
//             w->set_face(f);
//         }
//     
//         return true;
//     }
    

  // this function has been moved to the Tds class
//   bool remove(Vertex* v)
//   {
//     CGAL_triangulation_precondition(v != NULL);
//     CGAL_triangulation_precondition( has_vertex(v));
//      int i = index(v);
//     
//         Face* left, *right, *ll, *rr;
//     
//         left = neighbor(cw(i));
//         right = neighbor(ccw(i));
//     
// //         if(left == NULL || right == NULL) {
// //             if(left == NULL && right == NULL) {
// //                 Face* n = neighbor(i);
// //                 if(n != NULL) {
// //                     int ni = n->index(this);
// //                     n->set_neighbor(ni, NULL);
// //                     Vertex* q = vertex(cw(i));
// //                     if(q != NULL){
// //                         q->set_face(n);
// //                     }
// //                 } else {
// //                     cerr << "removal of boundary vertex failed as face has"
// //                          << "no neighbors" << endl;
// //                 }
// //                 handle().Delete();
// //                 v.Delete();
// //                 return true;
// //             } else {
// //                 cerr << "removal of boundary vertex with degree != 2 failed";
// //                 cerr << endl;
// //                 return false;
// //             }
// //         }
// //     
//         if(v->degree() != 3){
//             cerr << "removal of internal vertex with degree != 3 failed";
//             cerr << endl;
//             return false;
//         }
//     
//         int li = left->index(this);
// 	int ri = right->index(this);
//         Vertex* q = left->vertex(li);
// 	CGAL_triangulation_assertion( left->vertex(li) == right->vertex(ri));
//     
//         ll = left->neighbor(cw(li));
//         if(ll != NULL) {
//             int lli = ll->index(left);
//             ll->set_neighbor(lli, this);
//         } 
//         set_neighbor(cw(i), ll);
// 	if (vertex(ccw(i))->face() == left) vertex(ccw(i))->set_face(this);    
//         
//     
//         
//         rr = right->neighbor(ccw(ri));
//         if(rr != NULL) {
//             int rri = rr->index(right);
//             rr->set_neighbor(rri, this);
//         } 
//         set_neighbor(ccw(i), rr);
// 	if (vertex(cw(i))->face() == right) vertex(cw(i))->set_face(this);  
//         
// 	set_vertex(i, q);
// 	if (q->face() == right || q->face() == left) {
// 	   q->set_face(this);
// 	}
// 
// 	delete right;
//     	delete left;
//         
//         delete v;
//         return true;   
//     }
//     

//   void flip(int i)
//   {
//     Face* n  = neighbor(i);
//     
//     Vertex* v_cw = vertex(cw(i));
//     Vertex*  v_ccw = vertex(ccw(i));
// 
//     // we should not attempt to flip two faces which are adjacent on two edges
//     // This configuration happens in 1-dim triangulation
// 
//     int ni;
//     CGAL_triangulation_assertion( n->has_vertex(v_cw,ni));
//     ni = cw(ni);
//     CGAL_triangulation_assertion( vertex(i) != n->vertex(ni));
//     CGAL_triangulation_assertion( this == n->neighbor(ni) );
// 
//     // Old stuff
//     // The following seems natural, but it may fail if the faces
//     // this and n are neighbors on two edges (1-dim triangulation,
//     // with infinite faces
//     // int ni = n->index(this);
//     //
//     //  int ni = cw(n->index(v_cw));
//     //  CGAL_triangulation_assertion( this == n->neighbor(ni) );
//     
//     // bl == bottom left, tr == top right
//     Face* tr  = neighbor(ccw(i));
//     Face* bl = n->neighbor(ccw(ni));
//     int bli, tri;
// 
//     // Old stuff which seems natural
//     // but makes problem if f and tr or n and bl are incident through two edges
//     // bli = bl->index(n);
//     // tri = tr->index(f);
//     bli = 3 - ( bl->index(n->vertex(ni)) + bl->index(n->vertex(cw(ni))) );
//     tri = 3 - (tr->index(vertex(i)) + tr->index(vertex(cw(i))));
//     
//     set_vertex(cw(i), n->vertex(ni));
//     n->set_vertex(cw(ni), vertex(i));
//     
//     // update the neighborhood relations
//     set_neighbor(i, bl);
//     bl->set_neighbor(bli, this);
//            
//     set_neighbor(ccw(i), n);
//     n->set_neighbor(ccw(ni), this);
//     
//     n->set_neighbor(ni, tr);
//     tr->set_neighbor(tri,n);    
//     
//     if(v_cw->face() == this) {
//       v_cw->set_face(n);
//     }
//     
//     if(v_ccw->face() == n) {
//       v_ccw->set_face(this);
//     }
//   }

   bool is_valid(bool verbose = false, int level = 0) const
  {
    bool result = Fb::is_valid();
    for(int i = 0; i < 3; i++) {
      Face* n = neighbor(i);
            
      // The following seems natural, but it may fail if the faces
      // this and n are neighbors on two edges (1-dim triangulation,
      // with infinite faces
      // int ni = n->index(this);

      //  int ni = cw(n->index(vertex(cw(i))));
      // CGAL_triangulation_assertion( this == n->neighbor(ni) );
      // result = result && (vertex(cw(i)) == n->vertex(ccw(ni)));
      // result = result && (vertex(ccw(i)) == n->vertex(cw(ni)));

      int in;
      if (! n->has_vertex(vertex(cw(i)),in )) return false;
      in = cw(in); 
      result = result && ( this == n->neighbor(in) );
      result = result && (vertex(ccw(i)) == n->vertex(cw(in)));

    }
    return result;
  }
   

};

#endif CGAL_TRIANGULATION_DS_FACE_2_H
