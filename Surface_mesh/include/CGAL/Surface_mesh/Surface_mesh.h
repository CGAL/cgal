//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public License
// as published by the Free Software Foundation, version 2.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//=============================================================================


#ifndef CGAL_SURFACE_MESH_H
#define CGAL_SURFACE_MESH_H

#include <iterator>
#include <algorithm>
#include <utility>
#include <iostream>
#include <cstddef>

#include <boost/array.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>

#include <CGAL/circulator.h>
#include <CGAL/assertions.h>
#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/Surface_mesh/IO.h>
#include <CGAL/Surface_mesh/Properties.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {
  /// \ingroup PkgSurface_mesh
  /// This class is a data structure that can be used as halfedge data structure or polyhedral
  /// surface. It is an alternative to the classes `HalfedgeDS` and `Polyhedron_3`
  /// defined in the packages  \ref PkgHDSSummary and \ref PkgPolyhedronSummary. 
  /// The main difference is that it is indexed based and not pointer based,
  /// and that the mechanism for adding information to vertices, halfedges,
  /// and faces is much simpler and done at runtime and not at compile time.
  /// @tparam P The type of the "point" property of a vertex. There is no requirement on `P`,
  ///         besides being default constructible and assignable. 
  ///         In typical use cases it will be a 2D or 3D point type.


template <typename P>
class Surface_mesh
{

    typedef Surface_mesh<P> Self;

    template<typename>
    class Handle_iterator;
public:
    /// \name Basic Types
    ///
    ///@{

    /// The point type.
    typedef P Point;

    /// The type used to represent an index.
    typedef int size_type;

    ///@}

    /// \name Basic Elements
    ///
    ///@{

    /// Base class for vertex, halfedge, edge, and face descriptor. 
    ///
    /// \attention Note that `Descriptor` is not a model of the concept `Handle`,
    /// because it cannot be dereferenced.
    /// \sa `Vertex_descriptor`, `Halfedge_descriptor`, `Edge_descriptor`, `Face_descriptor`.
    template<typename T>
    class Descriptor
    {
    public:
        /// Constructor. %Default construction creates an invalid descriptor.
        /// We write -1, which is <a href="http://en.cppreference.com/w/cpp/concept/numeric_limits">
        /// <tt>std::numeric_limits<size_type>::max()</tt></a>
        /// as `size_type` is an unsigned type. 
        explicit Descriptor(size_type _idx=-1) : idx_(_idx) {}

        /// Get the underlying index of this descriptor
        size_type idx() const { return idx_; }

        /// reset descriptor to be invalid (index=-1)
        void reset() { idx_=-1; }

        /// return whether the descriptor is valid, i.e., the index is not equal to -1.
        bool is_valid() const { return idx_ != -1; }

        /// are two descriptors equal?
        bool operator==(const T& _rhs) const {
            return idx_ == _rhs.idx_;
        }

        /// are two descriptors different?
        bool operator!=(const T& _rhs) const {
            return idx_ != _rhs.idx_;
        }

        /// Comparison by index.
        bool operator<(const T& _rhs) const {
            return idx_ < _rhs.idx_;
        }

        /// increments the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// increment.
        Descriptor& operator++() { ++idx_; return *this; }
        /// decrements the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// decrement.
        Descriptor& operator--() { --idx_; return *this; }

        /// increments the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// increment.
        Descriptor operator++(int) { Descriptor tmp(*this); ++idx_; return tmp; }
        /// decrements the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// decrement.
        Descriptor operator--(int) { Descriptor tmp(*this); --idx_; return tmp; }
    private:
        size_type idx_;
    };


    /// This class represents a vertex.
    /// \sa `Halfedge_descriptor`, `Edge_descriptor`, `Face_descriptor`
    class Vertex_descriptor : public Descriptor<Vertex_descriptor>
    {
    public:
        /// %Default constructor (with invalid index).
        explicit Vertex_descriptor(size_type _idx=-1) : Descriptor<Vertex_descriptor>(_idx) {}

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Vertex_descriptor const& v)
        {
            return (os << 'v' << v.idx());
        }
    };

    /// This class represents a halfedge.
    /// \sa `Vertex_descriptor`, `Edge_descriptor`, `Face_descriptor`
    class Halfedge_descriptor : public Descriptor<Halfedge_descriptor>
    {
    public:
        /// %Default constructor (with invalid index)
        explicit Halfedge_descriptor(size_type _idx=-1) : Descriptor<Halfedge_descriptor>(_idx) {}

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Halfedge_descriptor const& h)
        {
            return (os << 'h' << h.idx());
        }
    };

    /// This class represents a face
    /// \sa `Vertex_descriptor`, `Halfedge_descriptor`, `Edge_descriptor`
    class Face_descriptor : public Descriptor<Face_descriptor>
    {
    public:
        /// %Default constructor (with invalid index)
        explicit Face_descriptor(size_type _idx=-1) : Descriptor<Face_descriptor>(_idx) {}

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Face_descriptor const& f)
        {
            return (os << 'f' << f.idx());
        }
    };

    /// This class represents an edge.
    ///
    /// Although it does not derive
    /// from `Descriptor` it is conceptually the same. \sa `Vertex_descriptor`,
    /// `Halfedge_descriptor`, `Face_descriptor`
    class Edge_descriptor
    {
    public:
        /// %Default constructor (with invalid index).
        Edge_descriptor(size_type idx=-1) : halfedge_(idx * 2) { }

        /// constructs an `Edge_descriptor` from a halfedge.
        Edge_descriptor(Halfedge_descriptor he) : halfedge_(he) { }
        /// @cond CGAL_DOCUMENT_INTERNALS
        /// returns the internal halfedge.
        Halfedge_descriptor halfedge() const { return halfedge_; }

        /// returns the underlying index of this descriptor.
        size_type idx() const { return halfedge_.idx() / 2; }

        /// resets descriptor to be invalid (index=-1)
        void reset() { halfedge_.reset(); }

        /// returns whether the descriptor is valid, i.e., the index is not equal to -1.
        bool is_valid() const { return halfedge_.is_valid(); }

        /// Are two descriptors equal?
        bool operator==(const Edge_descriptor& other) const { return this->idx() == other.idx(); }

        /// Are two descriptors different?
        bool operator!=(const Edge_descriptor& other) const { return this->idx() != other.idx(); }

        /// compares by index.
        bool operator<(const Edge_descriptor& other) const { return this->idx() < other.idx();}

        /// decrements the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// decrement.
        Edge_descriptor& operator--() { halfedge_ = Halfedge_descriptor(halfedge_.idx() - 2); return *this; }

        /// increments the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// increment.
        Edge_descriptor& operator++() { halfedge_ = Halfedge_descriptor(halfedge_.idx() + 2); return *this; }

        /// decrements internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// decrement.
        Edge_descriptor operator--(int) { Edge_descriptor tmp(*this); halfedge_ = Halfedge_descriptor(halfedge_.idx() - 2); return tmp; }

        /// increments internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// increment.
        Edge_descriptor operator++(int) { Edge_descriptor tmp(*this); halfedge_ = Halfedge_descriptor(halfedge_.idx() + 2); return tmp; }

        /// @endcond 

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Edge_descriptor const& e)
        {
            return (os << 'e' << e.idx() << " on " << e.halfedge());
        }
    private:
        Halfedge_descriptor halfedge_;
    };


 
    ///@}
private: //-------------------------------------------------- connectivity types

    /// This type stores the vertex connectivity
    /// \sa `Halfedge_connectivity`, `Face_connectivity`
    struct Vertex_connectivity
    {
        /// an incoming halfedge per vertex (it will be a border halfedge for border vertices)
        Halfedge_descriptor  halfedge_;
    };


    /// This type stores the halfedge connectivity
    /// \sa `Vertex_connectivity`, `Face_connectivity`
    struct Halfedge_connectivity
    {
        /// face incident to halfedge
        Face_descriptor      face_;
        /// vertex the halfedge points to
        Vertex_descriptor    vertex_;
        /// next halfedge within a face (or along a border)
        Halfedge_descriptor  next_halfedge_;
        /// previous halfedge within a face (or along a border)
        Halfedge_descriptor  prev_halfedge_;
    };


    /// This type stores the face connectivity
    /// \sa `Vertex_connectivity`, `Halfedge_connectivity`
    struct Face_connectivity
    {
        /// a halfedge that is part of the face
        Halfedge_descriptor  halfedge_;
    };

private: //------------------------------------------------------ iterator types
    template<typename Descriptor_>
    class Descriptor_iterator
      : public boost::iterator_facade< Descriptor_iterator<Descriptor_>,
                                       Descriptor_,
                                       std::bidirectional_iterator_tag
                                       >
    {
        typedef boost::iterator_facade< Descriptor_iterator<Descriptor_>,
                                        Descriptor_,
                                        std::bidirectional_iterator_tag
                                        > Facade;
    public:
        Descriptor_iterator() : hnd_(), mesh_(NULL) {}
        Descriptor_iterator(const Descriptor_& h, const Surface_mesh* m)
          : hnd_(h), mesh_(m) {
            if (mesh_ && mesh_->has_garbage())
              while (mesh_->has_valid_index(hnd_) && mesh_->is_removed(hnd_)) ++hnd_;
        }
    private:
        friend class boost::iterator_core_access;
        void increment()
        {
            ++hnd_;
            CGAL_assertion(mesh_ != NULL);
            while (mesh_->has_garbage() && mesh_->has_valid_index(hnd_) && mesh_->is_removed(hnd_)) ++hnd_;
        }

        void decrement()
        {
            --hnd_;
            CGAL_assertion(mesh_ != NULL);
            while (mesh_->has_garbage() && mesh_->has_valid_index(hnd_) && mesh_->is_removed(hnd_)) --hnd_;
        }

        bool equal(const Descriptor_iterator& other) const
        {
            return this->hnd_ == other.hnd_;
        }

        Descriptor_& dereference() const { return const_cast<Descriptor_&>(hnd_); }

        Descriptor_ hnd_;
        const Surface_mesh* mesh_;
    };
public:
    /// \name Iterator and Range Types
    ///
    ///@{

    /// \brief This class iterates linearly over all vertices. 
    ///
    /// A model of `BidirectionalIterator` with value type `Vertex_descriptor`.
    /// \sa `vertices()`
    /// \sa `Halfedge_iterator`, `Edge_iterator`, `Face_iterator`
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Vertex_iterator;
#else
    typedef Descriptor_iterator<Vertex_descriptor> Vertex_iterator;
#endif 

    /// \brief The range over all vertex descriptors.
    ///
    /// A model of <a href="http://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a>.
    /// \sa `vertices()`
    /// \sa `Halfedge_range`, `Edge_range`, `Face_range`
    typedef std::pair<Vertex_iterator, Vertex_iterator> Vertex_range;

    /// \brief This class iterates linearly over all halfedges.
    ///
    /// A model of `BidirectionalIterator` with value type `Halfedge_descriptor`.
    /// \sa `halfedges()`
    /// \sa `Vertex_iterator`, `Edge_iterator`, `Face_iterator`
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Halfedge_iterator;
#else
    typedef Descriptor_iterator<Halfedge_descriptor> Halfedge_iterator;
#endif

    /// \brief The range over all halfedge descriptors.
    ///
    /// A model of <a href="http://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a>.
    /// \sa `halfedges()`
    /// \sa `Vertex_range`, `Edge_range`, `Face_range`
    typedef std::pair<Halfedge_iterator, Halfedge_iterator> Halfedge_range;

    /// \brief This class iterates linearly over all edges.
    ///
    /// A model of `BidirectionalIterator` with value type `Halfedge_descriptor`.
    /// \sa `edges()`
    /// \sa `Vertex_iterator`, `Halfedge_iterator`, `Face_iterator`
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Edge_iterator;
#else
    typedef Descriptor_iterator<Edge_descriptor> Edge_iterator;
#endif

    /// \brief The range over all edge descriptors.
    ///
    /// A model of <a href="http://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a>.
    /// \sa `edges()`
    /// \sa `Halfedge_range`, `Vertex_range`, `Face_range`
    typedef std::pair<Edge_iterator, Edge_iterator> Edge_range;

    /// \brief This class iterates linearly over all faces.
    ///
    /// A model of `BidirectionalIterator` with value type `Face_descriptor`.
    /// \sa `faces()`
    /// \sa `Vertex_iterator`, `Halfedge_iterator`, `Edge_iterator`
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Face_iterator;
#else
    typedef Descriptor_iterator<Face_descriptor> Face_iterator;
#endif

    /// \brief The range over all face descriptors.
    ///
    /// A model of <a href="http://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a>.
    /// \sa `faces()`
    /// \sa `Vertex_range`, `Halfedge_range`, `Edge_range`
    typedef std::pair<Face_iterator, Face_iterator> Face_range;


  typedef CGAL::Vertex_around_target_iterator<Surface_mesh> Vertex_around_target_iterator;
  typedef std::pair<Vertex_around_target_iterator,Vertex_around_target_iterator> Vertex_around_target_range;

  typedef CGAL::Halfedge_around_target_iterator<Surface_mesh>  Halfedge_around_target_iterator;
  typedef std::pair<Halfedge_around_target_iterator,Halfedge_around_target_iterator> Halfedge_around_target_range;

  typedef CGAL::Face_around_target_iterator<Surface_mesh>  Face_around_target_iterator;
  typedef std::pair<Face_around_target_iterator,Face_around_target_iterator> Face_around_target_range;

  typedef CGAL::Vertex_around_face_iterator<Surface_mesh>  Vertex_around_face_iterator;
  typedef std::pair<Vertex_around_face_iterator,Vertex_around_face_iterator> Vertex_around_face_range;

  typedef CGAL::Halfedge_around_face_iterator<Surface_mesh>  Halfedge_around_face_iterator;
  typedef std::pair<Halfedge_around_face_iterator,Halfedge_around_face_iterator> Halfedge_around_face_range;

  typedef CGAL::Face_around_face_iterator<Surface_mesh>  Face_around_face_iterator;
  typedef std::pair<Face_around_face_iterator,Face_around_face_iterator> Face_around_face_range;

   
    /// @cond CGAL_BEGIN_END
    /// Start iterator for vertices.
    Vertex_iterator vertices_begin() const
    {
        return Vertex_iterator(Vertex_descriptor(0), this);
    }

    /// End iterator for vertices.
    Vertex_iterator vertices_end() const
    {
        return Vertex_iterator(Vertex_descriptor(num_vertices()), this);
    }
    /// @endcond

    /// returns the iterator range of the vertices of the mesh.
    Vertex_range vertices() const {
      return std::make_pair(vertices_begin(), vertices_end());
    }

    /// @cond CGAL_BEGIN_END
    /// Start iterator for halfedges.
    Halfedge_iterator halfedges_begin() const
    {
        return Halfedge_iterator(Halfedge_descriptor(0), this);
    }

    /// End iterator for halfedges.
    Halfedge_iterator halfedges_end() const
    {
        return Halfedge_iterator(Halfedge_descriptor(num_halfedges()), this);
    }
    /// @endcond

    /// returns the iterator range of the halfedges of the mesh.
    Halfedge_range halfedges() const {
      return std::make_pair(halfedges_begin(), halfedges_end());
    }


    /// @cond CGAL_BEGIN_END
    /// Start iterator for edges.
    Edge_iterator edges_begin() const
    {
        return Edge_iterator(Edge_descriptor(0), this);
    }

    /// End iterator for edges.
    Edge_iterator edges_end() const
    {
        return Edge_iterator(Edge_descriptor(num_edges()), this);
    }
    /// @endcond

    /// returns the iterator range of the edges of the mesh.
    Edge_range edges() const
    {
        return std::make_pair(edges_begin(), edges_end());
    }


    /// @cond CGAL_BEGIN_END
    /// Start iterator for faces.
    Face_iterator faces_begin() const
    {
        return Face_iterator(Face_descriptor(0), this);
    }

    /// End iterator for faces.
    Face_iterator faces_end() const
    {
        return Face_iterator(Face_descriptor(num_faces()), this);
    }
    /// @endcond

    /// returns the iterator range of the faces of the mesh.
    Face_range faces() const {
      return std::make_pair(faces_begin(), faces_end());
    }

    /// returns the iterator range for vertices around vertex `target(h)`, starting at `source(h)`.
    Vertex_around_target_range vertices_around_target(Halfedge_descriptor h) const
    {
      return CGAL::vertices_around_target(h,*this);
    }

    /// returns the iterator range for incoming halfedges around vertex `target(h)`, starting at `h`.
    Halfedge_around_target_range halfedges_around_target(Halfedge_descriptor h) const
    {
      return CGAL::halfedges_around_target(h,*this);
    }

    /// returns the iterator range for faces around vertex `target(h)`, starting at `face(h)`.
    Face_around_target_range faces_around_target(Halfedge_descriptor h) const
    {
      return CGAL::faces_around_target(h,*this);
    }

    /// returns the iterator range for vertices around face `face(h)`, starting at `target(h)`.
    Vertex_around_face_range vertices_around_face(Halfedge_descriptor h) const
     {
       return CGAL::vertices_around_face(h,*this);
     }

    /// returns the iterator range for halfedges around face `face(h)`, starting at `h`.
    Halfedge_around_face_range halfedges_around_face(Halfedge_descriptor h) const
    {
      return CGAL::halfedges_around_face(h,*this);
    }

    /// returns the iterator range for halfedges around face `face(h)`, starting at `h`.
    Face_around_face_range faces_around_face(Halfedge_descriptor h) const
    {
       return CGAL::faces_around_face(h,*this);
    }

    ///@}


public: 
    /// \name Circulator Types
    ///
    /// The following circulators enable to iterate through the elements around a face or vertex.
    /// As explained in the \ref SurfaceMeshOrientation "User Manual", we can speak of a  
    /// *clockwise* or *counterclockwise*
    /// traversal, by looking at the surface from the right side.  
    ///@{

    /// \brief This class circulates clockwise through all 
    /// one-ring neighbors of a vertex. 
    ///  A model of `BidirectionalCirculator` with value type `Vertex_descriptor`.
    /// \sa `Halfedge_around_target_circulator`, `Face_around_target_circulator`

  typedef CGAL::Vertex_around_target_circulator<Surface_mesh> Vertex_around_target_circulator;



    /// \brief This class circulates clockwise through all incident faces of a vertex.
    ///  A model of `BidirectionalCirculator` with value type `Face_descriptor`.
    /// \sa `Vertex_around_target_circulator`, `Halfedge_around_target_circulator`

  typedef CGAL::Face_around_target_circulator<Surface_mesh> Face_around_target_circulator;


    /// \brief This class circulates clockwise through all halfedges around a vertex that have this vertex as target.
    ///  A model of `BidirectionalCirculator` with value type `Halfedge_descriptor`.
    /// \sa `Vertex_around_target_circulator`, `Halfedge_around_target_circulator`

  typedef CGAL::Halfedge_around_target_circulator<Surface_mesh> Halfedge_around_target_circulator;


    /// \brief This class circulates clockwise through all halfedges around a vertex that have this vertex as source.
    ///  A model of `BidirectionalCirculator` with value type `Halfedge_descriptor`.
    /// \sa `Vertex_around_target_circulator`, `Halfedge_around_target_circulator`

  typedef CGAL::Halfedge_around_source_circulator<Surface_mesh> Halfedge_around_source_circulator;

    /// \brief This class circulates counterclockwise through all vertices around a face.
    ///  A model of `BidirectionalCirculator` with value type `Vertex_descriptor`.

  typedef  CGAL::Vertex_around_face_circulator<Surface_mesh> Vertex_around_face_circulator;


    /// \brief This class circulates counterclockwise through all halfedges around a face.
    ///  A model of `BidirectionalCirculator` with value type `Halfedge_descriptor`.

  typedef  CGAL::Halfedge_around_face_circulator<Surface_mesh> Halfedge_around_face_circulator;

   /// \brief This class circulates counterclockwise through all faces around a face.
   ///  A model of `BidirectionalCirculator` with value type `Face_descriptor`.
   ///  Note that face descriptor is the same after `operator++`, if the faces share 
   ///  several halfedges.

  typedef  CGAL::Face_around_face_circulator<Surface_mesh> Face_around_face_circulator;
  /// @}

  /// @cond CGAL_DOCUMENT_INTERNALS
  // typedefs which make it easier to write the partial specialisation of boost::graph_traits

  typedef Vertex_descriptor   vertex_descriptor;
  typedef P                   vertex_property_type;
  typedef Halfedge_descriptor halfedge_descriptor;
  typedef Edge_descriptor     edge_descriptor;
  typedef Face_descriptor     face_descriptor;

  typedef Vertex_iterator     vertex_iterator;
  typedef Halfedge_iterator   halfedge_iterator;
  typedef Edge_iterator       edge_iterator;
  typedef Face_iterator      face_iterator;
  typedef CGAL::Out_edge_iterator<Self>     out_edge_iterator;

  typedef boost::undirected_tag             directed_category;
  typedef boost::disallow_parallel_edge_tag edge_parallel_category; 

  struct traversal_category : public virtual boost::bidirectional_graph_tag,
                              public virtual boost::vertex_list_graph_tag,
                              public virtual boost::edge_list_graph_tag
  {};

  typedef size_type vertices_size_type;
  typedef size_type halfedges_size_type;
  typedef size_type edges_size_type;
  typedef size_type faces_size_type;
  typedef size_type degree_size_type;

 /// @endcond
public:

    /// \name Construction, Destruction, Assignment
    ///
    ///  Copy constructors as well as assignment do also copy simplices marked as removed.
    ///@{

    /// %Default constructor.
    Surface_mesh();

    /// Copy constructor: copies `rhs` to `*this`. Performs a deep copy of all properties.
    Surface_mesh(const Surface_mesh& rhs) { *this = rhs; }

    /// assigns `rhs` to `*this`. Performs a deep copy of all properties.
    Surface_mesh& operator=(const Surface_mesh& rhs);

    /// assigns `rhs` to `*this`. Does not copy custom properties.
    Surface_mesh& assign(const Surface_mesh& rhs);

    ///@}

public:

    /// \name Adding Vertices, Edges, and Faces
    ///@{

   /// adds a new vertex, and resizes vertex properties if necessary.
    Vertex_descriptor add_vertex()
    {
      if(vertices_freelist_ != -1){
        size_type idx = vertices_freelist_;
        vertices_freelist_ = vconn_[Vertex_descriptor(vertices_freelist_)].halfedge_.idx();
        --removed_vertices_;
        vremoved_[Vertex_descriptor(idx)] = false;
        return Vertex_descriptor(idx);
      } else {
        vprops_.push_back();
        return Vertex_descriptor(num_vertices()-1);
      }
    }

    /// adds a new vertex, resizes vertex properties if necessary,
    /// and sets the point property to `p`.
    Vertex_descriptor add_vertex(const Point& p) 
    {
        Vertex_descriptor v = add_vertex();
        vpoint_[v] = p;
        return v;
    }



public:

    /// adds a new edge, and resizes edge and halfedge properties if necessary.
    Halfedge_descriptor add_edge()
    {
      Halfedge_descriptor h0, h1;
      if(edges_freelist_ != -1){
        size_type idx = edges_freelist_;
        edges_freelist_ = hconn_[Halfedge_descriptor(edges_freelist_)].next_halfedge_.idx();
        --removed_edges_;
        eremoved_[Edge_descriptor(Halfedge_descriptor(idx))] = false;
        return Halfedge_descriptor(idx);
      } else {
        eprops_.push_back();
        hprops_.push_back();
        hprops_.push_back();

        return Halfedge_descriptor(num_halfedges()-2);
      }
    }

    /// adds two opposite halfedges, and resizes edge and halfedge properties if necessary.
    /// \returns the halfedge with `v1` as target
    Halfedge_descriptor add_edge(Vertex_descriptor v0, Vertex_descriptor v1)
    {
        CGAL_assertion(v0 != v1);
        Halfedge_descriptor h = add_edge();

        set_target(h, v1);
        set_target(opposite(h), v0);

        return h;
    }

    /// adds a new face, and resizes face properties  if necessary.
    Face_descriptor add_face()
    {
      if(faces_freelist_ != -1){
        size_type idx = faces_freelist_;
        faces_freelist_ = fconn_[Face_descriptor(faces_freelist_)].halfedge_.idx();
        --removed_faces_;
        fremoved_[Face_descriptor(idx)] = false;
        return Face_descriptor(idx);
      } else {
        fprops_.push_back();
        return Face_descriptor(num_faces()-1);
      }
    }

    /// adds a new face with vertices from a random access container.
    template <typename RandomAccessContainer>
    Face_descriptor add_face(const RandomAccessContainer& vertices);

    /// \todo Offer a variadic version
    /// adds a new triangle connecting vertices `v0`, `v1`, `v2`
    Face_descriptor add_face(Vertex_descriptor v0, Vertex_descriptor v1, Vertex_descriptor v2)
    {
        boost::array<Vertex_descriptor, 3> 
            v = {{v0, v1, v2}};
        return add_face(v);
    }

  /// adds a new quad connecting vertices `v0`, `v1`, `v2`, `v3`.
    Face_descriptor add_face(Vertex_descriptor v0, Vertex_descriptor v1, Vertex_descriptor v2, Vertex_descriptor v3)
    {
        boost::array<Vertex_descriptor, 4> 
            v = {{v0, v1, v2, v3}};
        return add_face(v);
    }
    ///@}


 
    /// \name Low-Level Removal Functions 
    ///
    /// Although the elements are only marked as removed
    /// their connectivity and properties should not be used.
    ///
    /// \warning Functions in this group do not adjust any of
    /// connected elements and usually leave the surface mesh in an
    /// invalid state.
    /// 
    ///
    /// @{

    /// removes vertex `v` from the halfedge data structure without
    /// adjusting anything.
    void remove_vertex(Vertex_descriptor v)
    {
        if(!vremoved_) vremoved_ = property_map<Vertex_descriptor, bool>("v:removed", false);
        vremoved_[v] = true; ++removed_vertices_; garbage_ = true;
        vconn_[v].halfedge_ = Halfedge_descriptor(vertices_freelist_);
        vertices_freelist_ = v.idx();
    }

    /// removes the two halfedges corresponding to `e` from the halfedge data structure without
    /// adjusting anything.
    void remove_edge(Edge_descriptor e)
    {
        eremoved_[e] = true; ++removed_edges_; garbage_ = true;
        hconn_[Halfedge_descriptor(e.idx() << 1)].next_halfedge_ = Halfedge_descriptor(edges_freelist_ );
        edges_freelist_ = (e.idx() << 1);
    }

    /// removes  face `f` from the halfedge data structure without
    /// adjusting anything.

    void remove_face(Face_descriptor f)
    {
        if(!fremoved_) fremoved_ = property_map<Face_descriptor, bool>("f:removed", false);
        fremoved_[f] = true; ++removed_faces_; garbage_ = true;
        fconn_[f].halfedge_ = Halfedge_descriptor(faces_freelist_);
        faces_freelist_ = f.idx();
    }


    ///@}


    /// \name Memory Management
    ///
    /// Functions to check the number of used elements, the amount of space
    /// allocated for elements, to clear the structure, and to perform
    /// garbage collection.
    ///@{

 
    /// returns the number of used and removed vertices in the mesh.
    size_type num_vertices() const { return (size_type) vprops_.size(); }

    /// returns the number of used and removed halfedges in the mesh.
    size_type num_halfedges() const { return (size_type) hprops_.size(); }

    /// returns the number of used and removed edges in the mesh.
    size_type num_edges() const { return (size_type) eprops_.size(); }

    /// returns the number of used and removed faces in the mesh.
    size_type num_faces() const { return (size_type) fprops_.size(); }


    /// returns the number of removed vertices in the mesh.
    size_type num_removed_vertices() const { return removed_vertices_; }

    /// returns the number of removed halfedges in the mesh.
    size_type num_removed_halfedges() const { return 2*removed_edges_; }

    /// returns the number of removed edges in the mesh.
    size_type num_removed_edges() const { return removed_edges_; }

    /// returns the number of removed faces in the mesh.
    size_type num_removed_faces() const { return removed_faces_; }


    /// returns `true` iff the mesh is empty, i.e., has no used vertices.
    bool empty() const { return num_vertices() == num_removed_vertices(); }

    /// removes all vertices, edges and faces. Collects garbage and clears all properties.
    void clear();

    /// @cond CGAL_DOCUMENT_INTERNALS
    /// removes unused memory from vectors. This shrinks the storage
    /// of all properties to the minimal required size.
    /// \attention Invalidates all existing references to properties.
    void shrink_to_fit()
    {
        vprops_.shrink_to_fit();
        hprops_.shrink_to_fit();
        eprops_.shrink_to_fit();
        fprops_.shrink_to_fit();
    }
    /// @endcond

    /// reserves space for vertices, faces, edges and their currently
    /// associated properties.
    void reserve(size_type nvertices,
                 size_type nedges,
                 size_type nfaces )
    {
        vprops_.reserve(nvertices);
        hprops_.reserve(2*nedges);
        eprops_.reserve(nedges);
        fprops_.reserve(nfaces);
    }

    /// checks if any vertices, faces, or edges are marked as removed.
    /// \sa collect_garbage
    bool has_garbage() const { return garbage_; }

    /// really removes vertices, edges, and faces which were marked removed.
    /// \sa `has_garbage()`
    /// \todo Add a version which creates a property for the previous position
    /// of an element. We really would need a mapping old->new, but we no longer
    /// have arrays of the right size for old properties. How to solve that elegantly?
  void collect_garbage();

    /// returns whether vertex `v` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Vertex_descriptor v) const
    {
        return vremoved_[v];
    }
    /// returns whether halfedge `h` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Halfedge_descriptor h) const
    {
        return eremoved_[edge(h)];
    }
    /// returns whether edge `e` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Edge_descriptor e) const
    {
        return eremoved_[e];
    }
    /// returns whether face `f` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Face_descriptor f) const
    {
        return fremoved_[f];
    }

    ///@}

    /// @cond CGAL_DOCUMENT_INTERNALS
    ///
    /// \name Simple Validity Checks
    ///
    /// Functions in this group check if the index is valid, that is between
    /// `0` and the currently allocated maximum amount of the
    /// elements. They do not check if an element is marked as removed.
    ///@{

    /// returns whether the index of vertex `v` is valid, that is within the current array bounds.
    bool has_valid_index(Vertex_descriptor v) const
    {
        return (0 <= v.idx()) && (v.idx() < (int)num_vertices());
    }
    /// returns whether the index of halfedge `h` is valid, that is within the current array bounds.
    bool has_valid_index(Halfedge_descriptor h) const
    {
        return (0 <= h.idx()) && (h.idx() < (int)num_halfedges());
    }
    /// returns whether the index of edge `e` is valid, that is within the current array bounds.
    bool has_valid_index(Edge_descriptor e) const
    {
        return (0 <= e.idx()) && (e.idx() < (int)num_edges());
    }
    /// returns whether the index of face `f` is valid, that is within the current array bounds.
    bool has_valid_index(Face_descriptor f) const
    {
        return (0 <= f.idx()) && (f.idx() < (int)num_faces());
    }

    /// @}
    /// @endcond
    
    /// \name Validity Checks
    ///
    /// Functions in this group perform checks for structural
    /// consistency of the object. They are expensive and should only
    /// be used in debug configurations.

    ///@{

    /// perform an expensive validity check on the data structure and
    /// print found errors to `std::cerr` when `verbose == true`.
  bool is_valid(bool verbose = true) const
    {
        bool valid = true;
        for(Halfedge_iterator it = halfedges_begin(); it != halfedges_end(); ++it) { 
            valid = valid && next(*it).is_valid();
            valid = valid && opposite(*it).is_valid();
            if(!valid) {
                std::cerr << "Integrity of halfedge " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (opposite(*it) != *it);
            valid = valid && (opposite(opposite(*it)) == *it);
            if(!valid) {
                std::cerr << "Integrity of opposite halfedge of " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (next(prev(*it)) == *it);
            if(!valid) {
                std::cerr << "Integrity of previous halfedge of " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (prev(next(*it)) == *it);
            if(!valid) {
                std::cerr << "Integrity of next halfedge of " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && target(*it).is_valid();
            if(!valid) {
                std::cerr << "Integrity of vertex of halfedge " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (target(*it) == target(opposite(next(*it))));
            if(!valid) {
                std::cerr << "Halfedge vertex of next opposite is not the same for " << *it << "."  << std::endl;
                break;
            }
        }

        for(Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it) { 
            if(halfedge(*it).is_valid()) {
                // not an isolated vertex
                valid = valid && (target(halfedge(*it)) == *it);
                if(!valid) {
                    std::cerr << "Halfedge of " << *it << " is not an incoming halfedge." << std::endl;
                    break;
                }
            }
        }


        return valid;
    }

    /// performs a validity check on a single vertex.
    bool is_valid(Vertex_descriptor v) const {
        Halfedge_descriptor h = vconn_[v].halfedge_;
        if(h!= null_halfedge() && (!has_valid_index(h) || is_removed(h))) {
            std::cerr << "Vertex connectivity halfedge error in " << v.idx()
                      << " with " << h.idx() << std::endl;
            return false;
        }
        return true;
    }

    /// performs a validity check on a single halfedge.
    bool is_valid(Halfedge_descriptor h) const {
        Face_descriptor f = hconn_[h].face_;
        Vertex_descriptor v = hconn_[h].vertex_;
        Halfedge_descriptor hn = hconn_[h].next_halfedge_;
        Halfedge_descriptor hp = hconn_[h].prev_halfedge_;

        bool valid = true;
        // don't validate the face if this is a border halfedge
        if(!is_border(h)) {
            if(!has_valid_index(f) || is_removed(f)) {
                std::cerr << "Halfedge connectivity Face "
                          << (!has_valid_index(f) ? "invalid" : "removed")
                          << " in " << h.idx() << std::endl;
                valid = false;
            }
        }

        if(!has_valid_index(v) || is_removed(v)) {
            std::cerr << "Halfedge connectivity Vertex "
                      << (!has_valid_index(v) ? "invalid" : "removed")
                      << " in " << h.idx() << std::endl;
            valid = false;
        }

        if(!has_valid_index(hn) || is_removed(hn)) {
            std::cerr << "Halfedge connectivity hnext "
                      << (!has_valid_index(hn) ? "invalid" : "removed")
                      << " in " << h.idx() << std::endl;
            valid = false;
        }
        if(!has_valid_index(hp) || is_removed(hp)) {
            std::cerr << "Halfedge connectivity hprev "
                      << (!has_valid_index(hp) ? "invalid" : "removed")
                      << " in " << h.idx() << std::endl;
            valid = false;
        }
        return valid;
    }

    /// performs a validity check on a single face.
    bool is_valid(Face_descriptor f) const {
        Halfedge_descriptor h = fconn_[f].halfedge_;
        if(!has_valid_index(h) || is_removed(h)) {
            std::cerr << "Face connectivity halfedge error in " << f.idx()
                      << " with " << h.idx() << std::endl;
            return false;
        }
        return true;
    }

    ///@}



    /// \name Low-Level Connectivity
    ///@{

    /// returns the number of incident halfedges of vertex `v`.
    size_type degree(Vertex_descriptor v) const;

    /// returns the number of incident halfedges of face `f`.
    size_type degree(Face_descriptor f) const;

    /// returns an incoming halfedge of vertex `v`.
    /// If `v` is a border vertex this will be a border halfedge.
    /// \invariant `target(halfedge(v)) == v`
    Halfedge_descriptor halfedge(Vertex_descriptor v) const
    {
        return vconn_[v].halfedge_;
    }

    /// sets the incoming halfedge of vertex `v` to `h`.
    void set_halfedge(Vertex_descriptor v, Halfedge_descriptor h)
    {
        vconn_[v].halfedge_ = h;
    }

    /// returns whether `v` is isolated, i.e., not incident to any face.
    /// \todo Should it be "not incident to a halfedge"?
    bool is_isolated(Vertex_descriptor v) const
    {
        return !halfedge(v).is_valid();
    }

    /// returns the vertex the halfedge `h` emanates from.
    Vertex_descriptor source(Halfedge_descriptor h) const
    {
        return target(opposite(h));
    }

    /// returns the vertex the halfedge `h` points to.
    Vertex_descriptor target(Halfedge_descriptor h) const
    {
        return hconn_[h].vertex_;
    }

    /// sets the vertex the halfedge `h` points to to `v`.
    void set_target(Halfedge_descriptor h, Vertex_descriptor v)
    {
        hconn_[h].vertex_ = v;
    }

    /// returns the face incident to halfedge `h`.
    Face_descriptor face(Halfedge_descriptor h) const
    {
        return hconn_[h].face_;
    }

    /// sets the incident face to halfedge `h` to `f`.
    void set_face(Halfedge_descriptor h, Face_descriptor f)
    {
        hconn_[h].face_ = f;
    }

    /// returns the next halfedge within the incident face.
    Halfedge_descriptor next(Halfedge_descriptor h) const
    {
        return hconn_[h].next_halfedge_;
    }

    /// @cond CGAL_DOCUMENT_INTERNALS
    // sets the next halfedge of `h` within the face to `nh`.
    void set_next_only(Halfedge_descriptor h, Halfedge_descriptor nh)
    {
      hconn_[h].next_halfedge_ = nh;
    }

    // sets previous halfedge of `h` to `nh`.
    void set_prev_only(Halfedge_descriptor h, Halfedge_descriptor nh)
    {
      if(h != null_halfedge()){
        hconn_[h].prev_halfedge_ = nh;
      }
    }
    /// @endcond  

    /// sets the next halfedge of `h` within the face to `nh` and
    /// the previous halfedge of `nh` to `h`.
    void set_next(Halfedge_descriptor h, Halfedge_descriptor nh)
    {
      set_next_only(h, nh);
      set_prev_only(nh, h);
    }

    /// returns the previous halfedge within the incident face.
    Halfedge_descriptor prev(Halfedge_descriptor h) const
    {
        return hconn_[h].prev_halfedge_;
    }

    /// returns the opposite halfedge of `h`.
    Halfedge_descriptor opposite(Halfedge_descriptor h) const
    {
        return Halfedge_descriptor((h.idx() & 1) ? h.idx()-1 : h.idx()+1);
    }

    /// returns `opposite(next(h))`, that is the next halfedge \ref SurfaceMeshOrientation 
    /// "clockwise" around the target vertex of `h`. 
    Halfedge_descriptor next_around_target(Halfedge_descriptor h) const
    {
        return prev(opposite(h));
    }

    /// returns `opposite(next(h))`, that is the next halfedge \ref SurfaceMeshOrientation
    /// "counterclockwise" around the target vertex of `h`. 
    Halfedge_descriptor prev_around_target(Halfedge_descriptor h) const
    {
        return opposite(next(h));
    }

    /// returns the edge that contains halfedge `h` as one of its two halfedges.
    Edge_descriptor edge(Halfedge_descriptor h) const
    {
        return Edge_descriptor(h);
    }

    /// returns the halfedge corresponding to the edge `e`.
    Halfedge_descriptor halfedge(Edge_descriptor e) const
    {
        return Halfedge_descriptor(e.halfedge());
    }

    /// returns the i'th halfedge of edge `e`, for `i=0` or `1`.
    Halfedge_descriptor halfedge(Edge_descriptor e, unsigned int i) const
    {
        CGAL_assertion(i<=1);
        return Halfedge_descriptor((e.idx() << 1) + i);
    }

    /// returns the i'th vertex of edge `e`, for `i=0` or `1`.
    Vertex_descriptor vertex(Edge_descriptor e, unsigned int i) const
    {
        CGAL_assertion(i<=1);
        return target(halfedge(e, i));
    }

    /// returns a halfedge of face `f`.
    Halfedge_descriptor halfedge(Face_descriptor f) const
    {
        return fconn_[f].halfedge_;
    }

    /// sets the halfedge of face `f` to `h`.
    void set_halfedge(Face_descriptor f, Halfedge_descriptor h)
    {
        fconn_[f].halfedge_ = h;
    }

    /// finds a halfedge between two vertices. Returns a default constructed
    /// `Halfedge_descriptor`, if  `source` and  `target` are not connected.
    Halfedge_descriptor halfedge(Vertex_descriptor source, Vertex_descriptor target) const;

    ///@}


    /// \name Borders
    ///@{

    /// returns whether `v` is a border vertex.
    bool is_border(Vertex_descriptor v) const
    {
        Halfedge_descriptor h(halfedge(v));
        return (!(h.is_valid() && face(h).is_valid()));
    }

    /// returns whether `h` is a border halfege, i.e., if its face equals `Face_descriptor()`.
    bool is_border(Halfedge_descriptor h) const
    {
        return !face(h).is_valid();
    }


    /// returns whether `e` is a border edge, i.e., if any 
    /// of its two halfedges is a border halfedge.
    bool is_border(Edge_descriptor e) const
    {
      return is_border(e.halfedge()) || is_border(opposite(e.halfedge()));
    }


    /// returns whether `f` is a border face, i.e., if the opposite 
    /// of any of its adjacent halfedges is a border halfedge.
    bool is_border(Face_descriptor f) const
    {
        Halfedge_descriptor h  = halfedge(f);
        Halfedge_descriptor hh = h;
        do
        {
            if (is_border(opposite(h)))
                return true;
            h = next(h);
        }
        while (h != hh);
        return false;
    }
  
  /// restores the constant time border property for vertex `v`.
  void fix_border(Vertex_descriptor v)
  {
    if(halfedge(v) == null_halfedge())
      return;
    Halfedge_around_target_circulator hatc(halfedge(v),*this), done(hatc);
    do {
      if(is_border(*hatc)){
        set_halfedge(v,*hatc);
        return;
      }
    }while(++hatc != done);
  }

  /// restores the constant time border property for all vertices 
  /// around the face associated to `h`.
  void fix_border(Halfedge_descriptor h)
  {
    if(is_border(h)){
      Halfedge_around_face_circulator hafc(h,*this),done(hafc);
      do {
        set_halfedge(target(*hafc),*hafc);
      }while(++hafc != done);
    } else {
       Vertex_around_face_circulator vafc(h,*this),done(vafc);
      do {
        fix_border(*vafc);
      }while(++vafc != done);
    }
  }

  /// restores the constant time border property for all vertices
  /// of the surface mesh.
  void fix_border()
  {
    BOOST_FOREACH(Vertex_descriptor vd, vertices()){
      fix_border(vd);
    }
  }


    ///@}


private: //--------------------------------------------------- property handling

    /// @cond BROKEN_DOC
    typedef boost::fusion::map<
      boost::fusion::pair< typename Surface_mesh::Vertex_descriptor, Property_container<Vertex_descriptor> (Surface_mesh::*)>,
      boost::fusion::pair< typename Surface_mesh::Halfedge_descriptor, Property_container<Halfedge_descriptor> (Surface_mesh::*)>,
      boost::fusion::pair< typename Surface_mesh::Edge_descriptor, Property_container<Edge_descriptor> (Surface_mesh::*)>,
      boost::fusion::pair< typename Surface_mesh::Face_descriptor, Property_container<Face_descriptor> (Surface_mesh::*)> >
        map_type;

    map_type pmap_;
    /// @endcond

    public:
 

 /*! \name Property Handling

    */
    ///@{


    /// adds a property map named `name` with value type `T` and default `t`
    /// for descriptor type `Desc`. Returns an invalid property map if a property
    /// map named `name` already exists.

  
    template<class Desc, class T>
    Property_map<Desc, T>
    add_property_map(const std::string& name, const T t=T()) {
        return (this->*boost::fusion::at_key<Desc>(pmap_)).template add<T>(name, t);
    }

 
    /// returns a property map named `name` with value type `T` for descriptor type `Desc`. 
    /// Returns an invalid property map, if no property map of
    /// matching type named `name` exists.
    template <class Desc, class T>
    Property_map<Desc, T> get_property_map(const std::string& name) const
    {
        return (this->*boost::fusion::at_key<Desc>(pmap_)).template get<T>(name);
    }



  /// @cond CGAL_DOCUMENT_INTERNALS
    /// retrieves a property with key `Desc` and value type `T` with 
    /// `name`. If no such property exists one is added with the default
    /// value `t`.
    template <class Desc, class T>
    Property_map<Desc, T> property_map(const std::string& name, const T t=T())
    {
        return (this->*boost::fusion::at_key<Desc>(pmap_)).template get_or_add<T>(name, t);
    }
  /// @endcond

    /// removes property map `p`. The memory allocated for that property map is
    /// freed.
    template<class Desc, class T>
    void remove_property_map(Property_map<Desc, T>& p)
    {
        (this->*boost::fusion::at_key<Desc>(pmap_)).remove(p);
    }

    /// @cond CGAL_DOCUMENT_INTERNALS
    /// returns the std::type_info of the value type of the
    /// property identified by `name`.  `typeid(void)` if `name`
    /// does not identify any property.
    ///
    /// @tparam Desc The key type of the property. 

    template<class Desc>
    const std::type_info& get_property_type(const std::string& name)
    {
        return (this->*boost::fusion::at_key<Desc>(pmap_)).get_type(name);
    }
  /// @endcond

    /// returns a list of all strings that describe properties with the key type `Desc`.
    /// @tparam Desc The key type of the properties.
    template<class Desc>
    std::vector<std::string> properties() const
    {
        return (this->*boost::fusion::at_key<Desc>(pmap_)).properties();
    }

    /// returns the property for "v:point".
    Property_map<Vertex_descriptor, Point>
    points() const { return vpoint_; }

    /// returns the point associated to vertex `v`.
    const Point&
    point(Vertex_descriptor v) const { return vpoint_[v]; }

    /// returns the point associated to vertex `v`.
    Point&
    point(Vertex_descriptor v) { return vpoint_[v]; }

    /// @cond CGAL_DOCUMENT_INTERNALS
    /// prints property statistics to the stream `out`. The output is human-readable but
    /// not machine-friendly.  
    ///
    void property_stats(std::ostream& out = std::cout) const;
    /// @endcond
    ///@}


 /// \name Null Elements
    ///@{

  /// returns `Vertex_descriptor(-1)`.
  static Vertex_descriptor null_vertex()
  {
    return vertex_descriptor(-1);
  }

  /// returns `Edge_descriptor(-1)`.
  static Edge_descriptor null_edge()
  {
    return edge_descriptor(-1);
  }
  /// returns `Halfedge_descriptor(-1)`.
  static Halfedge_descriptor null_halfedge()
  {
    return halfedge_descriptor(-1);
  }
  /// returns `Face_descriptor(-1)`.
  static Face_descriptor null_face()
  {
    return face_descriptor(-1);
  }
  /// @}

  
private: //--------------------------------------------------- helper functions


    /// make sure that the incoming halfedge of vertex v is a border halfedge
    /// if `v` is a border vertex.
    void adjust_incoming_halfedge(Vertex_descriptor v);

private: //------------------------------------------------------- private data
    Property_container<Vertex_descriptor> vprops_;
    Property_container<Halfedge_descriptor> hprops_;
    Property_container<Edge_descriptor> eprops_;
    Property_container<Face_descriptor> fprops_;

    Property_map<Vertex_descriptor, Vertex_connectivity>      vconn_;
    Property_map<Halfedge_descriptor, Halfedge_connectivity>  hconn_;
    Property_map<Face_descriptor, Face_connectivity>          fconn_;

    Property_map<Vertex_descriptor, bool>  vremoved_;
    Property_map<Edge_descriptor, bool>    eremoved_;
    Property_map<Face_descriptor, bool>    fremoved_;

    Property_map<Vertex_descriptor, Point>   vpoint_;

    size_type removed_vertices_;
    size_type removed_edges_;
    size_type removed_faces_;

    size_type vertices_freelist_;
    size_type edges_freelist_;
    size_type faces_freelist_;
    bool garbage_;
};

  /*! \addtogroup PkgSurface_mesh
   *
   * @{
   */

  /// \relates Surface_mesh
  /// Inserts the surface mesh in an output stream in Ascii OFF format.
  template <typename P>
  std::ostream& operator<<(std::ostream& os, const Surface_mesh<P>& sm)
  {
    return os;
  }
  /// \relates Surface_mesh
  /// Extracts the surface mesh from an input stream in Ascii OFF format.
  template <typename P>
  std::istream& operator>>(std::istream& is, Surface_mesh<P>& sm)
  {
    return is;
  }
 /*! @} */

template <typename P>
Surface_mesh<P>::
Surface_mesh()
  : pmap_(boost::fusion::make_pair< Surface_mesh::Vertex_descriptor >(&Surface_mesh::vprops_)
          , boost::fusion::make_pair< Surface_mesh::Halfedge_descriptor >(&Surface_mesh::hprops_)
          , boost::fusion::make_pair< Surface_mesh::Edge_descriptor >(&Surface_mesh::eprops_)
          , boost::fusion::make_pair< Surface_mesh::Face_descriptor >(&Surface_mesh::fprops_))
{
    // allocate standard properties
    // same list is used in operator=() and assign()
    vconn_    = add_property_map<Vertex_descriptor, Vertex_connectivity>("v:connectivity");
    hconn_    = add_property_map<Halfedge_descriptor, Halfedge_connectivity>("h:connectivity");
    fconn_    = add_property_map<Face_descriptor, Face_connectivity>("f:connectivity");
    vpoint_   = add_property_map<Vertex_descriptor, Point>("v:point");
    vremoved_ = add_property_map<Vertex_descriptor, bool>("v:removed", false);
    eremoved_ = add_property_map<Edge_descriptor, bool>("e:removed", false);
    fremoved_ = add_property_map<Face_descriptor, bool>("f:removed", false);

    removed_vertices_ = removed_edges_ = removed_faces_ = 0;
    vertices_freelist_ = edges_freelist_ = faces_freelist_ = -1;
    garbage_ = false;
}


//-----------------------------------------------------------------------------
template <typename P>
Surface_mesh<P>&
Surface_mesh<P>::
operator=(const Surface_mesh<P>& rhs)
{
    if (this != &rhs)
    {
        pmap_ = rhs.pmap_;

        // deep copy of property containers
        vprops_ = rhs.vprops_;
        hprops_ = rhs.hprops_;
        eprops_ = rhs.eprops_;
        fprops_ = rhs.fprops_;

        // property handles contain pointers, have to be reassigned
        vconn_    = property_map<Vertex_descriptor, Vertex_connectivity>("v:connectivity");
        hconn_    = property_map<Halfedge_descriptor, Halfedge_connectivity>("h:connectivity");
        fconn_    = property_map<Face_descriptor, Face_connectivity>("f:connectivity");
        vremoved_ = property_map<Vertex_descriptor, bool>("v:removed");
        eremoved_ = property_map<Edge_descriptor, bool>("e:removed");
        fremoved_ = property_map<Face_descriptor, bool>("f:removed");
        vpoint_   = property_map<Vertex_descriptor, P>("v:point");

        // how many elements are removed?
        removed_vertices_  = rhs.removed_vertices_;
        removed_edges_     = rhs.removed_edges_;
        removed_faces_     = rhs.removed_faces_;
        vertices_freelist_ = rhs.vertices_freelist_;
        edges_freelist_    = rhs.edges_freelist_;
        faces_freelist_    = rhs.faces_freelist_;
        garbage_           = rhs.garbage_;
    }

    return *this;
}


//-----------------------------------------------------------------------------
template <typename P>
Surface_mesh<P>&
Surface_mesh<P>::
assign(const Surface_mesh<P>& rhs)
{
    if (this != &rhs)
    {
        // clear properties
        vprops_.clear();
        hprops_.clear();
        eprops_.clear();
        fprops_.clear();

        // allocate standard properties
        vconn_    = add_property_map<Vertex_descriptor, Vertex_connectivity>("v:connectivity");
        hconn_    = add_property_map<Halfedge_descriptor, Halfedge_connectivity>("h:connectivity");
        fconn_    = add_property_map<Face_descriptor, Face_connectivity>("f:connectivity");
        vpoint_   = add_property_map<Vertex_descriptor, P>("v:point");
        vremoved_ = add_property_map<Vertex_descriptor, bool>("v:removed", false);
        eremoved_ = add_property_map<Edge_descriptor, bool>("e:removed", false);
        fremoved_ = add_property_map<Face_descriptor, bool>("f:removed", false);

        // copy properties from other mesh
        vconn_.array()     = rhs.vconn_.array();
        hconn_.array()     = rhs.hconn_.array();
        fconn_.array()     = rhs.fconn_.array();
        vpoint_.array()    = rhs.vpoint_.array();
        vremoved_.array()  = rhs.vremoved_.array();
        eremoved_.array()  = rhs.eremoved_.array();
        fremoved_.array()  = rhs.fremoved_.array();

        // resize (needed by property containers)
        vprops_.resize(rhs.num_vertices());
        hprops_.resize(rhs.num_halfedges());
        eprops_.resize(rhs.num_edges());
        fprops_.resize(rhs.num_faces());

        // how many elements are removed?
        removed_vertices_  = rhs.removed_vertices_;
        removed_edges_     = rhs.removed_edges_;
        removed_faces_     = rhs.removed_faces_;
        vertices_freelist_ = rhs.vertices_freelist_;
        edges_freelist_    = rhs.edges_freelist_;
        faces_freelist_    = rhs.faces_freelist_;
        garbage_           = rhs.garbage_;
    }

    return *this;
}

//-----------------------------------------------------------------------------
template <typename P>
void
Surface_mesh<P>::
clear()
{
    vprops_.resize(0);
    hprops_.resize(0);
    eprops_.resize(0);
    fprops_.resize(0);

    vprops_.shrink_to_fit();
    hprops_.shrink_to_fit();
    eprops_.shrink_to_fit();
    fprops_.shrink_to_fit();

    removed_vertices_ = removed_edges_ = removed_faces_ = 0;
    vertices_freelist_ = edges_freelist_ = faces_freelist_ = -1;
    garbage_ = false;
}

//-----------------------------------------------------------------------------
/// @cond CGAL_DOCUMENT_INTERNALS
template <typename P>
void
Surface_mesh<P>::
property_stats(std::ostream& out) const
{
    std::vector<std::string> props;

    out << "vertex properties:\n";
    props = properties<Vertex_descriptor>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;

    out << "halfedge properties:\n";
    props = properties<Halfedge_descriptor>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;

    out << "edge properties:\n";
    props = properties<Edge_descriptor>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;

    out << "face properties:\n";
    props = properties<Face_descriptor>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;
}
/// @endcond

//-----------------------------------------------------------------------------
template <typename P>
typename Surface_mesh<P>::Halfedge_descriptor
Surface_mesh<P>::
halfedge(Vertex_descriptor source, Vertex_descriptor target) const
{
    CGAL_assertion(has_valid_index(source) && has_valid_index(target));

    Halfedge_descriptor h  = halfedge(target);
    const Halfedge_descriptor hh = h;

    if (h.is_valid())
    {
        do
        {
            if (this->source(h) == source)
              return h;
            h = next_around_target(h);
        }
        while (h != hh);
    }

    return Halfedge_descriptor();
}


//-----------------------------------------------------------------------------
template <typename P>
void
Surface_mesh<P>::
adjust_incoming_halfedge(Vertex_descriptor v)
{
    Halfedge_descriptor h  = halfedge(v);
    Halfedge_descriptor hh = h;

    if (h.is_valid())
    {
        if (target(h) != v)
        {
            // wrong target, flip
            h = opposite(h);
            hh = h;
            set_halfedge(v, h);
        }

        do
        {
            if (is_border(h))
            {
                set_halfedge(v, h);
                return;
            }
            h = next_around_target(h);
        }
        while (h != hh);
    }
}

//-----------------------------------------------------------------------------

 /// @cond CGAL_DOCUMENT_INTERNALS

template <typename P>
template <typename RandomAccessContainer>
typename Surface_mesh<P>::Face_descriptor
Surface_mesh<P>::add_face(const RandomAccessContainer& vertices)
{
    Vertex_descriptor                   v;
    unsigned int             i, ii, n((int)vertices.size()), id;
    std::vector<Halfedge_descriptor>    halfedges(n);
    std::vector<bool>        is_new(n), needs_adjust(n, false);
    Halfedge_descriptor                 inner_next, inner_prev,
    outer_next, outer_prev,
    border_next, border_prev,
    patch_start, patch_end;

    // cache for set_next and vertex' set_halfedge
    typedef std::pair<Halfedge_descriptor, Halfedge_descriptor>  NextCacheEntry;
    typedef std::vector<NextCacheEntry>    NextCache;

    NextCache    next_cache;
    next_cache.reserve(3*n);


    // don't allow degenerated faces
    CGAL_assertion (n > 2);

    // test for topological errors
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        if ( !is_border(vertices[i]) )
        {
            std::cerr << "Surface_meshT::add_face: complex vertex " << vertices[i] << std::endl;
            return Face_descriptor();
        }

        halfedges[i] = halfedge(vertices[i], vertices[ii]);
        is_new[i]    = !halfedges[i].is_valid();
        if (!is_new[i] && !is_border(halfedges[i]))
        {
            std::cerr << "Surface_meshT::add_face: complex edge " << halfedges[i] << std::endl;
            std::cerr << std::boolalpha << is_border(halfedges[i]) << std::endl;
            std::cerr << target(halfedges[i]) << std::endl;

            return Face_descriptor();
        }
    }

    // re-link patches if necessary
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        if (!is_new[i] && !is_new[ii])
        {
            inner_prev = halfedges[i];
           inner_next = halfedges[ii];

            if (next(inner_prev) != inner_next)
            {
                // here comes the ugly part... we have to relink a whole patch

                // search a free gap
                // free gap will be between border_prev and border_next
                outer_prev = opposite(inner_next);
                outer_next = opposite(inner_prev);
                border_prev = outer_prev;
                do
                    border_prev = opposite(next(border_prev));
                while (!is_border(border_prev) || border_prev==inner_prev);
                border_next = next(border_prev);
                CGAL_assertion(is_border(border_prev));
                CGAL_assertion(is_border(border_next));


                // ok ?
                if (border_next == inner_next)
                {
                    std::cerr << "Surface_meshT::add_face: patch re-linking failed\n";
                    return Face_descriptor();
                }

                // other halfedges' descriptors
                patch_start = next(inner_prev);
                patch_end   = prev(inner_next);

                // relink
                next_cache.push_back(NextCacheEntry(border_prev, patch_start));
                next_cache.push_back(NextCacheEntry(patch_end, border_next));
                next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
            }
        }
    }
    // create missing edges
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
      if (is_new[i]){
            halfedges[i] = add_edge(vertices[i], vertices[ii]);
            assert(source(halfedges[i]) == vertices[i]);
      }
    // create the face
    Face_descriptor f = add_face();
    set_halfedge(f, halfedges[n-1]);

    // setup halfedges
    for (i=0, ii=1; i<n; ++i, ++ii, ii%=n)
    {
        v          = vertices[ii];
        inner_prev = halfedges[i];
        inner_next = halfedges[ii];

        id = 0;
        if (is_new[i])  id |= 1;
        if (is_new[ii]) id |= 2;

        if (id)
        {
            outer_prev = opposite(inner_next);
            outer_next = opposite(inner_prev);

            // set outer links
            switch (id)
            {
                case 1: // prev is new, next is old

                    border_prev = prev(inner_next);
                    next_cache.push_back(NextCacheEntry(border_prev, outer_next));
                    set_halfedge(v, border_prev);
                    break;

                case 2: // next is new, prev is old

                    border_next = next(inner_prev);
                    next_cache.push_back(NextCacheEntry(outer_prev, border_next));
                    set_halfedge(v, outer_prev);
                    break;

                case 3: // both are new

                    if (!halfedge(v).is_valid())
                    {
                        set_halfedge(v, outer_prev);
                        next_cache.push_back(NextCacheEntry(outer_prev, outer_next));
                    }
                    else
                    {
                        border_prev = halfedge(v);
                        border_next = next(border_prev);
                        next_cache.push_back(NextCacheEntry(border_prev, outer_next));
                        next_cache.push_back(NextCacheEntry(outer_prev, border_next));
                    }
                    break;
            }

            // set inner link
            next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
        }
        else {
            needs_adjust[ii] = true; // (halfedge(v) == inner_next); the code is over adjusting
        }

        // set face descriptor
        set_face(halfedges[i], f);
    }



    // process next halfedge cache
    typename NextCache::const_iterator ncIt(next_cache.begin()), ncEnd(next_cache.end());
    for (; ncIt != ncEnd; ++ncIt)
        set_next(ncIt->first, ncIt->second);



    // adjust vertices' halfedge descriptor
    for (i=0; i<n; ++i)
      if (true) //(needs_adjust[i])
            adjust_incoming_halfedge(vertices[i]);

    return f;
}

  /// @endcond

//-----------------------------------------------------------------------------
template <typename P>
typename Surface_mesh<P>::size_type
Surface_mesh<P>::
degree(Vertex_descriptor v) const
{
    size_type count(0);

    Vertex_around_target_circulator vvit(halfedge(v), *this);
    Vertex_around_target_circulator vvend = vvit;
    if(vvit) do
    {
        ++count;
    } while (++vvit != vvend);

    return count;
}


//-----------------------------------------------------------------------------
template <typename P>
typename Surface_mesh<P>::size_type
Surface_mesh<P>::
degree(Face_descriptor f) const
{
    size_type count(0);

    Vertex_around_face_circulator fvit(halfedge(f),*this);
    Vertex_around_face_circulator fvend = fvit;
    if(fvit) do {
        ++count;
    } while (++fvit != fvend);

    return count;
}

template <typename P>
void
Surface_mesh<P>::
collect_garbage()
{
    int  i, i0, i1,
    nV(num_vertices()),
    nE(num_edges()),
    nH(num_halfedges()),
    nF(num_faces());

    Vertex_descriptor    v;
    Halfedge_descriptor  h;
    Face_descriptor      f;


    // setup descriptor mapping
    Property_map<Vertex_descriptor, Vertex_descriptor>      vmap = add_property_map<Vertex_descriptor, Vertex_descriptor>("v:garbage-collection");
    Property_map<Halfedge_descriptor, Halfedge_descriptor>  hmap = add_property_map<Halfedge_descriptor, Halfedge_descriptor>("h:garbage-collection");
    Property_map<Face_descriptor, Face_descriptor>          fmap = add_property_map<Face_descriptor, Face_descriptor>("f:garbage-collection");
    for (i=0; i<nV; ++i)
        vmap[Vertex_descriptor(i)] = Vertex_descriptor(i);
    for (i=0; i<nH; ++i)
        hmap[Halfedge_descriptor(i)] = Halfedge_descriptor(i);
    for (i=0; i<nF; ++i)
        fmap[Face_descriptor(i)] = Face_descriptor(i);



    // really remove vertices
    if (nV > 0)
    {
        i0=0;  i1=nV-1;

        while (1)
        {
            // find first removed and last un-removed
            while (!vremoved_[Vertex_descriptor(i0)] && i0 < i1)  ++i0;
            while ( vremoved_[Vertex_descriptor(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            vprops_.swap(i0, i1);
        };

        // remember new size
        nV = vremoved_[Vertex_descriptor(i0)] ? i0 : i0+1;
    }

    // really remove edges
    if (nE > 0)
    {
        i0=0;  i1=nE-1;

        while (1)
        {
            // find first removed and last un-removed
            while (!eremoved_[Edge_descriptor(i0)] && i0 < i1) ++i0;
            while ( eremoved_[Edge_descriptor(i1)] && i0 < i1) --i1;
            if (i0 >= i1) break;

            // swap
            eprops_.swap(i0, i1);
            hprops_.swap(2*i0,   2*i1);
            hprops_.swap(2*i0+1, 2*i1+1);
        };

        // remember new size
        nE = eremoved_[Edge_descriptor(i0)] ? i0 : i0+1;
        nH = 2*nE;
    }


    // really remove faces
    if (nF > 0)
    {
        i0=0;  i1=nF-1;

        while (1)
        {
            // find 1st removed and last un-removed
            while (!fremoved_[Face_descriptor(i0)] && i0 < i1)  ++i0;
            while ( fremoved_[Face_descriptor(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            fprops_.swap(i0, i1);
        };

        // remember new size
        nF = fremoved_[Face_descriptor(i0)] ? i0 : i0+1;
    }


    // update vertex connectivity
    for (i=0; i<nV; ++i)
    {
        v = Vertex_descriptor(i);
        if (!is_isolated(v))
            set_halfedge(v, hmap[halfedge(v)]);
    }


    // update halfedge connectivity
    for (i=0; i<nH; ++i)
    {
        h = Halfedge_descriptor(i);
        set_target(h, vmap[target(h)]);
        set_next(h, hmap[next(h)]);
        if (!is_border(h))
            set_face(h, fmap[face(h)]);
    }


    // update descriptors of faces
    for (i=0; i<nF; ++i)
    {
        f = Face_descriptor(i);
        set_halfedge(f, hmap[halfedge(f)]);
    }

    // remove descriptor maps
    remove_property_map<Vertex_descriptor>(vmap);
    remove_property_map<Halfedge_descriptor>(hmap);
    remove_property_map<Face_descriptor>(fmap);

    // finally resize arrays
    vprops_.resize(nV); vprops_.shrink_to_fit();
    hprops_.resize(nH); hprops_.shrink_to_fit();
    eprops_.resize(nE); eprops_.shrink_to_fit();
    fprops_.resize(nF); fprops_.shrink_to_fit();

    removed_vertices_ = removed_edges_ = removed_faces_ = 0;
    vertices_freelist_ = edges_freelist_ = faces_freelist_ = -1;
    garbage_ = false;
}

} // CGAL

#endif /* CGAL_SURFACE_MESH_H */

