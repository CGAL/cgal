// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
//
// Author(s)     : Kapelushnik Lior <liorkape@post.tau.ac.il>

/*! \file
 * spherical arrangements of none intersecting arcs of great circles on a sphere
 */

#ifndef CGAL_SPHERICAL_MAP_H_
#define CGAL_SPHERICAL_MAP_H_

#include "CGAL/Cubical_gaussian_map_3.h"

#include <CGAL/Sphere_dcel.h>
#include <CGAL/Sphere_traits.h>

#include <CGAL/basic.h>

#include <list>

// bug fix include
#include <set>
#include <CGAL/Origin.h>

CGAL_BEGIN_NAMESPACE

/*
 this class represents a spherical map

 SphericalMapDcel - the data structures representing spherical map elements
 Sphere_traits - a traits class that contains a definition for spherical arc as curve
*/
template <class SphericalMapDcel, class Sphere_traits>
class Spherical_map: public SphericalMapDcel {
  
public:
   // spherical map types decleration
   typedef typename Sphere_traits::X_monotone_curve_2  X_monotone_curve_2;
   typedef typename Sphere_traits::Curve_2    Curve_2;
   typedef typename Sphere_traits::Direction_3    Direction_3;
   typedef typename SphericalMapDcel::Vertex_iterator   Vertex_iterator;
   typedef typename SphericalMapDcel::Vertex_handle   Vertex_handle;
   typedef typename SphericalMapDcel::Vertex     Vertex;
   typedef typename SphericalMapDcel::Halfedge_iterator  Halfedge_iterator;
   typedef typename SphericalMapDcel::Halfedge_handle   Halfedge_handle;
   typedef typename SphericalMapDcel::Halfedge    Halfedge;
   typedef typename SphericalMapDcel::Ccb_halfedge_circulator  Ccb_halfedge_circulator;
   typedef typename SphericalMapDcel::Face_iterator  Face_iterator;
   typedef typename SphericalMapDcel::Face_handle   Face_handle;
   typedef typename SphericalMapDcel::Face    Face;
   typedef typename SphericalMapDcel::Holes_iterator  Holes_iterator;
   typedef typename SphericalMapDcel::FaceList::size_type Size;
   typedef typename SphericalMapDcel::CGM    CGM;
   typedef typename Sphere_traits::Point_3    Point_3;
   
   /*
  empty constructor
  initializes map state with no inserted elements, also initializes the
  cubical gaussian map
   */
   Spherical_map(): SphericalMapDcel(), m_cgm(), m_dirty(true),
  m_vertices(), m_halfedges(), m_faces() {  
    m_cgm.init_arrangements();
   }

  /*
   copy constructor,
   copy the internal cubical representation
  */
  Spherical_map(const Spherical_map &smap):
    m_cgm(smap.m_cgm), m_dirty(true), m_vertices(), m_halfedges(), m_faces() {}

  
   /*
    destructor
  clears the pointers from internal elements to spherical elements
   */  
   ~Spherical_map() {
     clearData();  
   }

   /*
    assignment operator  
   */   
   Spherical_map &operator=(const Spherical_map &smap) {     
     clearData();
     m_cgm = smap;
     m_dirty = true;     
     return *this;
   }
  
  /*
   get the internal cubical gaussian map that represents the spherical map

   return value - the cubical gaussian map representing the spherical arrangement
  */
  const CGM *get_cgm() const {   
    return &m_cgm;
  }

  /*
   get the internal cubical gaussian map that represents the spherical map

   return value - the cubical gaussian map representing the spherical arrangement
   non const return value
  */
  CGM *get_cgm() {   
    return &m_cgm;
  }

  /*
   insert a general great circle arc that does not intersect other spherical curves

   cv - the general great circle arc
  */
  void insert(const X_monotone_curve_2 &cv) {
    // get arc's endpoints and plane normal direction
    Direction_3 stDir = cv.getStDir(),
    enDir = cv.getEnDir();
    Direction_3 normDir = cv.getPlaneDir();
    // check arc's angle state (less then eauals or larger then 180 degrees
    if (Direction_3(CGAL::cross_product(stDir.vector(), enDir.vector())) == normDir) {
      // if the cross product of start and end vectors has the direction of the
      // plane normal, the angle is less than 180 degrees (counter clock wise)
      insertLess180(stDir, enDir);
    } else if (stDir == -enDir) {
      // if endpoints in opposite direction the arc has exactly 180 degrees
      Direction_3 midDir; // a middle direction
      midDir = Direction_3(CGAL::cross_product(normDir.vector(),
        stDir.vector()));
      // insert curve diveded to two arcs with 90 degrees
      insertLess180(stDir, midDir);
      insertLess180(midDir, enDir);
    } else  {
      // the arc has more than 180 degree, divide to 3 smaller arcs
      insertLess180(stDir, -enDir);
      insertLess180(-enDir, -stDir);
      insertLess180(-stDir, enDir);
    }
  }

  /*
     insert a circle that does not intersect other spherical curves

   cirPlaneDir - a normal to the plane of the circle
  */
  void insertCircle(Direction_3 &cirPlaneDir) {
    // find a direction on the plane
    // the plane has form ax+by+cz=0 so a point is find by setting two
    // axis distance to 1 and substituting in the formula, x=1,y=1,z=(-a-b)/c
    Direction_3 aDir;
    if (cirPlaneDir.dz() != 0) {
      aDir =
        Direction_3(1,1,(-cirPlaneDir.dx()-cirPlaneDir.dy())/cirPlaneDir.dz());
    } else if (cirPlaneDir.dy() != 0) {
      aDir =
        Direction_3(1,(-cirPlaneDir.dx()-cirPlaneDir.dz())/cirPlaneDir.dy(), 1);
    } else if (cirPlaneDir.dx() != 0) {
      aDir =
        Direction_3((-cirPlaneDir.dy()-cirPlaneDir.dz())/cirPlaneDir.dx(), 1, 1);

    }
    // find another direction on the plane
    Direction_3 aDir2(CGAL::cross_product(cirPlaneDir.vector(), aDir.vector()));
    // add curves between found direction, both long and short curve
    X_monotone_curve_2 cv1(aDir, aDir2, true);
    X_monotone_curve_2 cv2(aDir, aDir2, false);
    insert(cv1);
    insert(cv2);
  }

   // bug definitions
   // insert with bug bypass patch
   // insert a curve of less than 180 degrees
  /*
   insert a curve with less than 180 degrees between endpoints

   stDir - the arc source endpoint direction
   enDir - the arc target endpoint direction
  */
  void insertLess180(Direction_3 stDir, Direction_3 enDir) {
    Projected_normal prSource, prTarget;
    // create cubical projected normals from the directions
    prSource.compute_projection(stDir.vector());
    prTarget.compute_projection(enDir.vector());

    ProjIterator pit, pit2; // iterators for projected normals set
    // search for source projected normal in set of inserted projected normals
    pit = m_projNorms.find(prSource);

    bool foundSrc = false, 
    foundTrg = false; // so far the projected normals were not found in the set

    if (pit == m_projNorms.end()) {
      // projected normal was not inserted earlier
    } else {
      // projected normal already inserted
      foundSrc = true;
      prSource = *pit; // take the already inserted projected normal
    }
    // search for target projected normal in set of inserted projected normals
    pit2 = m_projNorms.find(prTarget);
    if (pit2 == m_projNorms.end()) {
      // projected normal was not inserted earlier
    } else {
      // projected normal already inserted
      foundTrg = true;
      prTarget = *pit2; // take the already inserted projected normal
    }
    // insert the arc to the cubical gaussian map
    m_cgm.insert(prSource, prTarget, false); 
    // update cubical map vertices with extended vertices data
    markProjNorm(prSource, stDir);
    markProjNorm(prTarget, enDir);
    m_dirty = true; // internal structure no longer valid
    // add projected normal to set of projected normals if not already there
    if (!foundTrg) {
      m_projNorms.insert(prTarget);
    }
    if (!foundSrc) {
      m_projNorms.insert(prSource);
    }
  }
   

     // not working with bug, insert without a projected normals cache
   void insert2(const X_monotone_curve_2 &cv) {
     Direction_3 stDir = cv.getStDir(),
       enDir = cv.getEnDir();
     Projected_normal prSource, prTarget;
     prSource.compute_projection(stDir.vector());
     prTarget.compute_projection(enDir.vector());
     m_cgm.insert(prSource, prTarget, false);
     markProjNorm(prSource, stDir);
     markProjNorm(prTarget, enDir);
     m_dirty = true;
   }

  /*
   get an iterator to the first spherical vertex

   return value - an iterator to the first spherical vertex
  */
  Vertex_iterator vertices_begin() {
    update();
    return m_vertices.begin();
  }

  /*
   get an iterator to the past the end spherical vertex

   return value - an iterator to the past the end spherical vertex
  */
  Vertex_iterator vertices_end() {
    update();
    return m_vertices.end();
  }

  /*
   get an iterator to the first spherical halfedge

   return value - an iterator to the first spherical halfedge
  */
  Halfedge_iterator halfedges_begin() {
    update();
    return m_halfedges.begin();
  }

  /*
   get an iterator to the past the end spherical halfedge

   return value - an iterator to the past the end spherical halfedge
  */
  Halfedge_iterator halfedges_end() {
    update();
    return m_halfedges.end();
  }

  /*
   get an iterator to the first spherical face

   return value - an iterator to the first spherical face
  */
  Face_iterator faces_begin() {
    update();
    return m_faces.begin();
  }

  /*
   get an iterator to the past the end spherical face

   return value - an iterator to the past the end spherical face
  */
  Face_iterator faces_end() {
    update();
    return m_faces.end();
  }
  
  /*
   get the number of faces in the spherical arrangement

   return value - the number of faces in the spherical map
  */
  Size number_of_faces() {
    update();
    return m_faces.size();
  }

  /*
   get the number of halfedges in the spherical arrangement

   return value - the number of halfedges in the spherical map
  */
  Size number_of_halfedges() {
    update();
    return m_halfedges.size();
  }

  /*
   get the number of vertices in the spherical arrangement

   return value - the number of vertices in the spherical map
  */
  Size number_of_vertices() {
    update();
    return m_vertices.size();
  }

   /*
     update the internal representation 
   according to the cubical gaussian map representation
    */
  void update() {
    if (!m_dirty) { // already updated, data valid, no need to update again
      return;
    }
    // clear data from previous update
    clearData();
    // loop over all cubical maps and update structure
    unsigned int i;
    CGM_vertex_iterator vit;    
    CGM_halfedge_iterator heit;    
    CGM_face_iterator fit;    
    for (i=0;i<6;i++) { // loop over cubical faces arrangements
      // get a cubical face arrangement
      CGM_planar_map &curMap = m_cgm.get_arrangement(i); 
      // update vertices
      for (vit=curMap.vertices_begin(); vit!=curMap.vertices_end(); ++vit) {
        //loop over all map vertices
        if (vit->is_rep()) { // representative vertex, wrap and add
          m_vertices.push_back(Vertex(vit, &m_cgm));
          // update cubical vertex to point to spherical vertex
          Vertex_iterator *newVerIt = new Vertex_iterator;
          *newVerIt = m_vertices.end();
          --(*newVerIt);
          vit->setSphereVertex(newVerIt);
        }
      }
      // update halfedges
      for (heit=curMap.halfedges_begin(); heit!=curMap.halfedges_end(); ++heit) {
        // loop over all map halfedges
        heit->setMark(false);
        // check if halfedge is real and has a real target
        if ((!(heit->face())->is_unbounded()) &&
          (heit->get_is_real()) && (heit->target())->getReal()) {
          // a spherical halfedge is represented by a real cubical halfedge that
          // has a real target so if an arc span over more then one halfedges
          // on the cube, only one of these halfedges which has a real target
          // will represent the spherical halfedge
          m_halfedges.push_back(Halfedge(heit, &m_cgm));
          // update cubical halfedge to point to spherical halfedge
          Halfedge_iterator *newHEIt = new Halfedge_iterator;
          *newHEIt = m_halfedges.end();
          --(*newHEIt);
          heit->setSphereEdge(newHEIt);
        }
      }
      // update faces
      for (fit=curMap.faces_begin(); fit != curMap.faces_end(); ++fit) {
        // loop over all map faces
        if (!(fit->is_unbounded()) && !(fit->getMark())) {
          // an unmarked face, add spherical face and mark
          typename SphericalMapDcel::IntFaceList curFaces;
          // mark internal faces without deletion of mark
          curFaces.clear();
          // get a list of internal faces connected to the face
          SphericalMapDcel::connectedFacesRec(m_cgm, curFaces, fit);
          // add spherical face to list of faces
          m_faces.push_back(Face(fit, &m_cgm));
          // update cubical faces to point to the spherical face
          Face_iterator *newFaceIt = new Face_iterator;
          *newFaceIt = m_faces.end();
          --(*newFaceIt);
          typename SphericalMapDcel::IntFaceList::iterator flit;
          for (flit=curFaces.begin(); flit!=curFaces.end(); ++flit) {
            (*flit)->setSphereFace(newFaceIt);
          }
        }
      }

      // unmark the cgm faces of current cubical face map
      for (fit=curMap.faces_begin(); fit != curMap.faces_end(); ++fit) {
        fit->setMark(false);
      }
    }

    if (m_halfedges.size() == 0) {
      // one face of entire cube without outer ccb
      m_faces.begin()->m_doesOuterExist = false;
    }
       
    m_dirty = false; // the data structure is currently valid
  }

  /*
   print the content of the spherical map

   os - output stream
   return value - the output stream after writing the spherical map state to it
  */
  std::ostream &print(std::ostream &os) {
    update(); // make sure spherical data is valid
    os << "#------------------------------ Printing spherical arrangement" <<
      std::endl;
    os << "# ------------------------------------------------------------------" << 
      std::endl;       
    os << "#Printing number of vertices halfedges and faces in spherical arrangement" << std::endl;
    os << number_of_vertices() << " " <<
      number_of_halfedges() << " " <<
      number_of_faces() << std::endl;
    // print vertices information
    os << "# " << number_of_vertices() << " vertices" << std::endl;
    os << "# ------------------------------------------" << std::endl;
    Vertex_iterator vit;
    // print direction of each vertex
    for (vit = vertices_begin(); vit != vertices_end(); ++vit) {
      os << vit->direction() << std::endl;
    }
    // print halfedges information
    os << "# " << number_of_halfedges() << " halfedges" << std::endl;
    os << "# ------------------------------------------" << std::endl;
    Halfedge_iterator hit;
    // for each halfedge, print it's source direction and curve
    for (hit = halfedges_begin(); hit != halfedges_end(); ++hit) {
      os << hit->source()->direction() << " " <<
        hit->curve() << std::endl;
    }
    // print faces information
    os << "# " << number_of_faces() << " faces" << std::endl;
    os << "# ------------------------------------------" << std::endl;
    Face_iterator fit;
    // for each face, print it's outer ccb and holes
    for (fit = faces_begin(); fit != faces_end(); ++fit) {
      fit->update(); // make sure face has valid data
      os << "# writing face" << std::endl;
      os << "# ------------------------------------------" << std::endl;
      Ccb_halfedge_circulator chc;
      Ccb_halfedge_circulator circEn;
      unsigned int numEdges;
      if (!fit->does_outer_ccb_exist()) {
        // no outer ccb, map is empty
        os << "# UNBOUNDED" << std::endl;
        os << "# number halfedges on outer boundary" <<
          std::endl << 0 << std::endl;
      } else {
        // arrangement not empty, has outer ccb
        os << "# outer ccb" << std::endl;
        os << "# number halfedges on outer boundary" << std::endl;
           
        chc = fit->outer_ccb();
        circEn = chc;
        // count number of halfedges on outer ccb
        numEdges = 0;
        do {
          ++numEdges;
          ++chc;
        } while (chc != circEn);
        os << numEdges << std::endl;;
        // write the curves of outer ccb
        do {
          os << chc->curve() << std::endl;
          ++chc;
        } while (chc != circEn);
      }
      // write holes
      os << "# number of holes" << std::endl;
      os << fit->m_holes.size() << std::endl;
      Holes_iterator holeIt;
      for (holeIt = fit->holes_begin(); holeIt != fit->holes_end(); ++holeIt) {
        // for each hole
        os << "# inner ccb" << std::endl;
        os << "# number halfedges on inner boundary" << std::endl;
        chc = *holeIt;
        circEn = chc;
        // count number of halfedges on hole
        numEdges = 0;
        do {
          ++numEdges;
          ++chc;
        } while (chc != circEn);
        os << numEdges << std::endl;;
        // write curves of hole
        do {
          os << chc->curve() << std::endl;
          ++chc;
        } while (chc != circEn);
      }
      os << "# finished writing face" << std::endl;
      os << "# ------------------------------------------" << std::endl;       
    }
    os << "#------------------------------ End of spherical arrangement" << 
      std::endl;
    os << "# ------------------------------------------------------------------" << std::endl;       

    return os;
  }
 
private:
  // internal types
  typedef typename CGM::Arrangement      CGM_planar_map;
  typedef typename CGM::Projected_normal    Projected_normal;
  typedef typename SphericalMapDcel::VertexList   VertexList;
  typedef typename CGM::Arr_vertex_iterator         CGM_vertex_iterator;
  typedef typename SphericalMapDcel::HalfedgeList   HalfedgeList;
  typedef typename CGM::Arr_halfedge_iterator         CGM_halfedge_iterator;
  typedef typename SphericalMapDcel::FaceList     FaceList;
  typedef typename CGM::Arr_face_iterator         CGM_face_iterator;
  typedef typename Sphere_traits::Vector_3      Vector_3;

  /*
    clear internal lists and pointers to map elements from cubical elements
  */
  void clearData(void) {
    Vertex_iterator clVerIt;
    Halfedge_iterator clHEIt;
    Face_iterator clFaceIt;
    // clear all pointers from internal cubical vertices to spherical vertices
    for (clVerIt = m_vertices.begin(); clVerIt != m_vertices.end(); ++clVerIt) {
      delete (Vertex_iterator *)(clVerIt->sphereVertexPointer());
    }
    // clear all pointers from internal cubical halfedges to spherical halfedges
    for (clHEIt = m_halfedges.begin(); clHEIt != m_halfedges.end(); ++clHEIt) {
      delete (Halfedge_handle *)(clHEIt->sphereHalfedgePointer());
    }
    // clear all pointers from internal cubical faces to spherical faces
    for (clFaceIt = m_faces.begin(); clFaceIt != m_faces.end(); ++clFaceIt) {
      delete (Face_handle *)(clFaceIt->sphereFacePointer());
    }
    // clear the lists of spherical elements
    m_vertices.clear();
    m_halfedges.clear();
    m_faces.clear();
  }


  // bug fix definitions
  typedef typename Sphere_traits::Point_3    Point_3;
  /*
   a function object designed to compare two projected normals,
   two projected normals has the same relation as the directions
   of their projected normals
  */
  struct CompProjs {
    /*
     the function object operator, compares two projected normals

     p1, p2 - the two projected normals to compare
     return value - true if the direction of p1 < direction of p2
    */
    bool operator() (const Projected_normal &p1, const Projected_normal &p2) const
    {
      Point_3 pt1 = p1.get_projected_normal();
      Point_3 pt2 = p2.get_projected_normal();
      return (pt1<pt2);
    }
  };

   // a type for set of projection normals that are compared by their direction
   typedef std::set<Projected_normal, CompProjs>    ProjSet;
   typedef typename ProjSet::iterator          ProjIterator;

   ProjSet m_projNorms; // set of inserted projected normals

   CGM m_cgm; // the cgm representation of the map
   bool m_dirty; // are the internal representaions up to date
   VertexList m_vertices; // list of sphere vertices
   HalfedgeList m_halfedges; // list of sphere halfedges
   FaceList m_faces; // list of sphere faces


  /*
     update cgm vertices of projected normal, make sure one of the
     vertices is a representative, set them as real and set direction

   projNorm - the projected normal to update
   dir - the direction inserted for the projection normal
  */
  void markProjNorm(Projected_normal & projNorm, Direction_3 & dir) {
    unsigned int i;
    unsigned int oneVer; // will hold an index to a vertex of projected normal
    bool isRep = false; // is there a representative vertex
    for (i=0;i<3;i++) {
      if (projNorm.is_vertex_set(i)) { // if found a vertex of projected normal
        (projNorm.get_vertex(i))->setReal(); // mark vertex as real vertex
        (projNorm.get_vertex(i))->set_direction(dir); // set vertex direction
        if ((projNorm.get_vertex(i))->is_rep()) {
          // if vertex is representative, update the corresponding flag
          isRep = true;
        }
        oneVer = i; // i is a number of a vertex that is set
      }
    }
    if (!isRep) { // mark a new representative if none already exists
      projNorm.get_vertex(oneVer)->set_rep();
    }
  }
};

/*
 output operator for a spherical map

 os - the output stream
 sphere - the sphere to output
*/
template <class SphericalMapDcel, class Sphere_traits>
std::ostream & operator<<(std::ostream &os,
        Spherical_map<SphericalMapDcel, Sphere_traits>
        &sphere) 
{
  // print using the spherical print method
  return sphere.print(os);
}

CGAL_END_NAMESPACE

#endif
