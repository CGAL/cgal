// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
// Copyright (C) 2014 GeometryFactory
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_SURFACE_MESH_H
#define CGAL_SURFACE_MESH_H

#include <CGAL/license/Surface_mesh.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh/IO.h>
#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/Surface_mesh/Properties.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/circulator.h>
#include <CGAL/Handle_hash_function.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/property_map.h>

#include <boost/cstdint.hpp>
#include <boost/array.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

namespace CGAL {

#ifndef DOXYGEN_RUNNING
    /// Base class for vertex, halfedge, edge, and face index.
    ///
    /// \attention Note that `Index` is not a model of the concept `Handle`,
    /// because it cannot be dereferenced.
    /// \sa `Vertex_index`, `Halfedge_index`, `Edge_index`, `Face_index`.
    template<typename T>
    class SM_Index
    {
    public:
    typedef boost::uint32_t size_type;
        /// Constructor. %Default construction creates an invalid index.
        /// We write -1, which is <a href="https://en.cppreference.com/w/cpp/types/numeric_limits">
        /// <tt>(std::numeric_limits<size_type>::max)()</tt></a>
        /// as `size_type` is an unsigned type.
        explicit SM_Index(size_type _idx=(std::numeric_limits<size_type>::max)()) : idx_(_idx) {}

        /// Get the underlying index of this index
        operator size_type() const { return idx_; }

        /// reset index to be invalid (index=(std::numeric_limits<size_type>::max)())
        void reset() { idx_=(std::numeric_limits<size_type>::max)(); }

        /// return whether the index is valid, i.e., the index is not equal to `%std::numeric_limits<size_type>::max()`.
        bool is_valid() const {
          size_type inf = (std::numeric_limits<size_type>::max)();
          return idx_ != inf;
        }

        // Compatibility with OpenMesh handle
        size_type idx() const {
          return idx_;
        }
        // For convenience
        size_type id() const {
          return idx_;
        }

        /// increments the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// increment.
        SM_Index& operator++() { ++idx_; return *this; }
        /// decrements the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// decrement.
        SM_Index& operator--() { --idx_; return *this; }

        /// increments the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// increment.
        SM_Index operator++(int) { SM_Index tmp(*this); ++idx_; return tmp; }
        /// decrements the internal index. This operation does not
        /// guarantee that the index is valid or undeleted after the
        /// decrement.
        SM_Index operator--(int) { SM_Index tmp(*this); --idx_; return tmp; }

        SM_Index operator+=(std::ptrdiff_t n) { idx_ = size_type(std::ptrdiff_t(idx_) + n); return *this; }

    protected:
        size_type idx_;
    };

  template <class T>
  std::size_t hash_value(const SM_Index<T>&  i)
  {
    std::size_t ret = i;
    return ret;
  }


    // Implementation for Surface_mesh::Vertex_index
    class SM_Vertex_index
 : public SM_Index<SM_Vertex_index>
    {
    public:

        SM_Vertex_index() : SM_Index<SM_Vertex_index>((std::numeric_limits<size_type>::max)()) {}

        explicit SM_Vertex_index(size_type _idx) : SM_Index<SM_Vertex_index>(_idx) {}

        template<class T> bool operator==(const T&) const = delete;
        template<class T> bool operator!=(const T&) const = delete;
        template<class T> bool operator<(const T&) const = delete;

        /// are two indices equal?
        bool operator==(const SM_Vertex_index& _rhs) const {
            return this->idx_ == _rhs.idx_;
        }

        /// are two indices different?
        bool operator!=(const SM_Vertex_index& _rhs) const {
            return this->idx_ != _rhs.idx_;
        }

        /// Comparison by index.
        bool operator<(const SM_Vertex_index& _rhs) const {
          return this->idx_ < _rhs.idx_;
        }


        friend std::ostream& operator<<(std::ostream& os, SM_Vertex_index const& v)
        {
          return (os << 'v' << (size_type)v );
        }
    };

    // Implementation of Surface_mesh::Halfedge_index
    class SM_Halfedge_index
      : public SM_Index<SM_Halfedge_index>
    {
    public:

        // Workaround for a bug in g++4.4 in ADL for function next:
        // we provide the types needed for std::iterator_traits<Surface_mesh::halfedge_index>,
        // although this descriptor is not an iterator.
        typedef void iterator_category;
        typedef void value_type;
        typedef void difference_type;
        typedef void pointer;
        typedef void reference;

        SM_Halfedge_index() : SM_Index<SM_Halfedge_index>((std::numeric_limits<size_type>::max)()) {}

        explicit SM_Halfedge_index(size_type _idx) : SM_Index<SM_Halfedge_index>(_idx) {}

        template<class T> bool operator==(const T&) const = delete;
        template<class T> bool operator!=(const T&) const = delete;
        template<class T> bool operator<(const T&) const = delete;

        /// are two indices equal?
        bool operator==(const SM_Halfedge_index& _rhs) const {
            return this->idx_ == _rhs.idx_;
        }

        /// are two indices different?
        bool operator!=(const SM_Halfedge_index& _rhs) const {
            return this->idx_ != _rhs.idx_;
        }

        /// Comparison by index.
        bool operator<(const SM_Halfedge_index& _rhs) const {
          return this->idx_ < _rhs.idx_;
        }

        friend std::ostream& operator<<(std::ostream& os, SM_Halfedge_index const& h)
        {
          return (os << 'h' << (size_type)h );
        }
    };

    /// Implementation of Surfae_mesh::Face_index
    class SM_Face_index
      : public SM_Index<SM_Face_index>
    {
    public:

        SM_Face_index() : SM_Index<SM_Face_index>((std::numeric_limits<size_type>::max)()) {}

        explicit SM_Face_index(size_type _idx) : SM_Index<SM_Face_index>(_idx) {}

        template<class T> bool operator==(const T&) const = delete;
        template<class T> bool operator!=(const T&) const = delete;
        template<class T> bool operator<(const T&) const = delete;

        /// are two indices equal?
        bool operator==(const SM_Face_index& _rhs) const {
            return this->idx_ == _rhs.idx_;
        }

        /// are two indices different?
        bool operator!=(const SM_Face_index& _rhs) const {
            return this->idx_ != _rhs.idx_;
        }

        /// Comparison by index.
        bool operator<(const SM_Face_index& _rhs) const {
          return this->idx_ < _rhs.idx_;
        }

        friend std::ostream& operator<<(std::ostream& os, SM_Face_index const& f)
        {
          return (os << 'f' << (size_type)f );
        }
    };

    /// Implementation of Surface_mesh::Edge_index
    class SM_Edge_index
    {
    public:
        typedef boost::uint32_t size_type;

        SM_Edge_index() : halfedge_((std::numeric_limits<size_type>::max)()) { }

        explicit SM_Edge_index(size_type idx) : halfedge_(idx * 2) { }

        explicit SM_Edge_index(SM_Halfedge_index he) : halfedge_(he) { }

        // returns the internal halfedge.
        SM_Halfedge_index halfedge() const { return halfedge_; }

        // returns the underlying index of this index.
        operator size_type() const { return (size_type)halfedge_ / 2; }

        // compatibility with OpenMesh handles
        size_type idx() const { return (size_type)halfedge_ / 2; }

        // resets index to be invalid (index=(std::numeric_limits<size_type>::max)())
        void reset() { halfedge_.reset(); }

        // returns whether the index is valid, i.e., the index is not equal to (std::numeric_limits<size_type>::max)().
        bool is_valid() const { return halfedge_.is_valid(); }


        template<class T> bool operator==(const T&) const = delete;
        template<class T> bool operator!=(const T&) const = delete;
        template<class T> bool operator<(const T&) const = delete;

        // Are two indices equal?
        bool operator==(const SM_Edge_index& other) const { return (size_type)(*this) == (size_type)other; }

        // Are two indices different?
        bool operator!=(const SM_Edge_index& other) const { return (size_type)(*this) != (size_type)other; }

        // compares by index.
        bool operator<(const SM_Edge_index& other) const { return (size_type)(*this) < (size_type)other;}

        // decrements the internal index. This operation does not
        // guarantee that the index is valid or undeleted after the
        // decrement.
        SM_Edge_index& operator--() { halfedge_ = SM_Halfedge_index((size_type)halfedge_ - 2); return *this; }

        // increments the internal index. This operation does not
        // guarantee that the index is valid or undeleted after the
        // increment.
        SM_Edge_index& operator++() { halfedge_ = SM_Halfedge_index((size_type)halfedge_ + 2); return *this; }

        // decrements internal index. This operation does not
        // guarantee that the index is valid or undeleted after the
        // decrement.
        SM_Edge_index operator--(int) { SM_Edge_index tmp(*this); halfedge_ = SM_Halfedge_index((size_type)halfedge_ - 2); return tmp; }

        // increments internal index. This operation does not
        // guarantee that the index is valid or undeleted after the
        // increment.
        SM_Edge_index operator++(int) { SM_Edge_index tmp(*this); halfedge_ = SM_Halfedge_index((size_type)halfedge_ + 2); return tmp; }

        SM_Edge_index operator+=(std::ptrdiff_t n) { halfedge_ = SM_Halfedge_index(size_type(std::ptrdiff_t(halfedge_) + 2*n)); return *this; }

        // prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, SM_Edge_index const& e)
        {
          return (os << 'e' << (size_type)e << " on " << e.halfedge());
        }

        friend  std::size_t hash_value(const SM_Edge_index&  i)
        {
          return i;
        }

    private:
        SM_Halfedge_index halfedge_;
    };
#endif

  /// \ingroup PkgSurface_mesh
  /// This class is a data structure that can be used as halfedge data structure or polyhedral
  /// surface. It is an alternative to the classes `HalfedgeDS` and `Polyhedron_3`
  /// defined in the packages  \ref PkgHalfedgeDS and \ref PkgPolyhedron.
  /// The main difference is that it is indexed based and not pointer based,
  /// and that the mechanism for adding information to vertices, halfedges, edges,
  /// and faces is much simpler and done at runtime and not at compile time.
  /// When elements are removed, they are only marked as removed, and a garbage
  /// collection function must be called to really remove them.
  /// @tparam P The type of the \em point property of a vertex. There is no requirement on `P`,
  ///         besides being default constructible and assignable.
  ///         In typical use cases it will be a 2D or 3D point type.
  /// \cgalModels `MutableFaceGraph` and `FaceListGraph`
  ///
  /// \sa \ref PkgBGLConcepts "Graph Concepts"

template <typename P>
class Surface_mesh
{

    typedef Surface_mesh<P> Self;

    template<typename>
    class Handle_iterator;
public:

#ifndef DOXYGEN_RUNNING
  template <class I, class T>
  struct Property_map : Properties::Property_map_base<I, T, Property_map<I, T> >
  {
    typedef Properties::Property_map_base<I, T, Property_map<I, T> > Base;
    typedef typename Base::reference reference;
    Property_map() = default;
    Property_map(const Base& pm): Base(pm) {}
  };

  template <typename Key, typename T>
  struct Get_property_map {
    typedef Property_map<Key, T> type;
  };
#endif // DOXYGEN_RUNNING

    /// \name Basic Types
    ///
    ///@{

    /// The point type.
    typedef P Point;

    /// The type used to represent an index.
    typedef boost::uint32_t size_type;

    ///@}

    /// \name Basic Elements
    ///
    ///@{


#ifdef DOXYGEN_RUNNING

    /// This class represents a vertex.
    /// \cgalModels `Index`
    /// \cgalModels `LessThanComparable`
    /// \cgalModels `Hashable`
    /// \sa `Halfedge_index`, `Edge_index`, `Face_index`
    class Vertex_index
    {
    public:
        /// %Default constructor.
        Vertex_index(){}

        Vertex_index(size_type _idx){}

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Vertex_index const& v)
        {}
    };
#else
  typedef SM_Vertex_index Vertex_index;
#endif

#ifdef DOXYGEN_RUNNING

    /// This class represents a halfedge.
    /// \cgalModels `Index`
    /// \cgalModels `LessThanComparable`
    /// \cgalModels `Hashable`
    /// \sa `Vertex_index`, `Edge_index`, `Face_index`
    class Halfedge_index
    {
    public:
        /// %Default constructor
        Halfedge_index(){}

        Halfedge_index(size_type _idx){}

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Halfedge_index const& h)
        {
        }

    };
#else
  typedef SM_Halfedge_index Halfedge_index;
#endif

#ifdef DOXYGEN_RUNNING
    /// This class represents a face
    /// \cgalModels `Index`
    /// \cgalModels `LessThanComparable`
    /// \cgalModels `Hashable`
    /// \sa `Vertex_index`, `Halfedge_index`, `Edge_index`
    class Face_index
    {
    public:
        /// %Default constructor
        Face_index(){}

        Face_index(size_type _idx){}

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Face_index const& f)
        {}
    };
#else
  typedef SM_Face_index Face_index;
#endif

#ifdef DOXYGEN_RUNNING
    /// This class represents an edge.
    /// \cgalModels `Index`
    /// \cgalModels `LessThanComparable`
    /// \cgalModels `Hashable`
    /// \sa `Vertex_index`, `Halfedge_index`, `Face_index`
    class Edge_index
    {
    public:
        /// %Default constructor
        Edge_index(){}

        Edge_index(size_type idx){}

        /// constructs an `Edge_index` from a halfedge.
        Edge_index(Halfedge_index he){}

        /// prints the index and a short identification string to an ostream.
        friend std::ostream& operator<<(std::ostream& os, typename Surface_mesh::Edge_index const& e)
        {}
    };
#else
  typedef SM_Edge_index Edge_index;
#endif


    ///@}
#ifndef CGAL_TEST_SURFACE_MESH
private: //-------------------------------------------------- connectivity types
#endif

    /// This type stores the vertex connectivity
    /// \sa `Halfedge_connectivity`, `Face_connectivity`
    struct Vertex_connectivity
    {
        /// an incoming halfedge per vertex (it will be a border halfedge for border vertices)
        Halfedge_index  halfedge_;
    };


    /// This type stores the halfedge connectivity
    /// \sa `Vertex_connectivity`, `Face_connectivity`
    struct Halfedge_connectivity
    {
        /// face incident to halfedge
        Face_index      face_;
        /// vertex the halfedge points to
        Vertex_index    vertex_;
        /// next halfedge within a face (or along a border)
        Halfedge_index  next_halfedge_;
        /// previous halfedge within a face (or along a border)
        Halfedge_index  prev_halfedge_;
    };


    /// This type stores the face connectivity
    /// \sa `Vertex_connectivity`, `Halfedge_connectivity`
    struct Face_connectivity
    {
        /// a halfedge that is part of the face
        Halfedge_index  halfedge_;
    };

private: //------------------------------------------------------ iterator types
    template<typename Index_>
    class Index_iterator
      : public boost::iterator_facade< Index_iterator<Index_>,
                                       Index_,
                                       std::random_access_iterator_tag,
                                       Index_
                                       >
    {
        typedef boost::iterator_facade< Index_iterator<Index_>,
                                        Index_,
                                        std::random_access_iterator_tag,
                                        Index_> Facade;
    public:
        Index_iterator() : hnd_(), mesh_(nullptr) {}
        Index_iterator(const Index_& h, const Surface_mesh* m)
          : hnd_(h), mesh_(m) {
          if (mesh_ && mesh_->has_garbage()){
              while (mesh_->has_valid_index(hnd_) && mesh_->is_removed(hnd_)) ++hnd_;
          }
        }
    private:
        friend class boost::iterator_core_access;
        void increment()
        {
            ++hnd_;
            CGAL_assertion(mesh_ != nullptr);

            if(mesh_->has_garbage())
              while ( mesh_->has_valid_index(hnd_) && mesh_->is_removed(hnd_)) ++hnd_;
        }

        void decrement()
        {
            --hnd_;
            CGAL_assertion(mesh_ != nullptr);
            if(mesh_->has_garbage())
               while ( mesh_->has_valid_index(hnd_) && mesh_->is_removed(hnd_)) --hnd_;
        }

        void advance(std::ptrdiff_t n)
        {
            CGAL_assertion(mesh_ != nullptr);

            if (mesh_->has_garbage())
            {
              if (n > 0)
                for (std::ptrdiff_t i = 0; i < n; ++ i)
                  increment();
              else
                for (std::ptrdiff_t i = 0; i < -n; ++ i)
                  decrement();
            }
            else
              hnd_ += n;
        }

        std::ptrdiff_t distance_to(const Index_iterator& other) const
        {
            if (mesh_->has_garbage())
            {
              bool forward = (other.hnd_ > hnd_);

              std::ptrdiff_t out = 0;
              Index_iterator it = *this;
              while (!it.equal(other))
              {
                if (forward)
                {
                  ++ it;
                  ++ out;
                }
                else
                {
                  -- it;
                  -- out;
                }
              }
              return out;
            }

            // else
            return std::ptrdiff_t(other.hnd_) - std::ptrdiff_t(this->hnd_);
        }

        bool equal(const Index_iterator& other) const
        {
            return this->hnd_ == other.hnd_;
        }

        Index_ dereference() const { return hnd_; }

        Index_ hnd_;
        const Surface_mesh* mesh_;

    };
public:
    /// \name Range Types
    ///
    /// Each range `R` in this section has a nested type `R::iterator`,
    /// is convertible to `std::pair<R::iterator,R::iterator>`, so that one can use `boost::tie()`,
    /// and can be used with `BOOST_FOREACH()`, as well as with the C++11 range based for-loop.

    ///@{

#ifndef DOXYGEN_RUNNING
    typedef Index_iterator<Vertex_index> Vertex_iterator;
#endif

    /// \brief The range over all vertex indices.
    ///
    /// A model of <a href="https://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Vertex_index`.
    /// \sa `vertices()`
    /// \sa `Halfedge_range`, `Edge_range`, `Face_range`
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type Vertex_range;
#else
    typedef Iterator_range<Vertex_iterator> Vertex_range;
#endif

#ifndef DOXYGEN_RUNNING
    typedef Index_iterator<Halfedge_index> Halfedge_iterator;
#endif

    /// \brief The range over all halfedge indices.
    ///
    /// A model of <a href="https://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Halfedge_index`.
    /// \sa `halfedges()`
    /// \sa `Vertex_range`, `Edge_range`, `Face_range`
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type Halfedge_range;
#else
    typedef Iterator_range<Halfedge_iterator> Halfedge_range;
#endif

#ifndef DOXYGEN_RUNNING
    typedef Index_iterator<Edge_index> Edge_iterator;
#endif

    /// \brief The range over all edge indices.
    ///
    /// A model of <a href="https://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Edge_index`.
    /// \sa `edges()`
    /// \sa `Halfedge_range`, `Vertex_range`, `Face_range`
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type Edge_range;
#else
    typedef Iterator_range<Edge_iterator> Edge_range;
#endif


#ifndef DOXYGEN_RUNNING
    typedef Index_iterator<Face_index> Face_iterator;
#endif
    /// \brief The range over all face indices.
    ///
    /// A model of <a href="https://www.boost.org/libs/range/doc/html/range/concepts/bidirectional_range.html">BidirectionalRange</a> with value type `Face_index`.
    /// \sa `faces()`
    /// \sa `Vertex_range`, `Halfedge_range`, `Edge_range`
 #ifdef DOXYGEN_RUNNING
    typedef unspecified_type Face_range;
#else
   typedef Iterator_range<Face_iterator> Face_range;
#endif

#ifndef DOXYGEN_RUNNING

  typedef CGAL::Vertex_around_target_iterator<Surface_mesh> Vertex_around_target_iterator;
  typedef Iterator_range<Vertex_around_target_iterator> Vertex_around_target_range;

  typedef CGAL::Halfedge_around_target_iterator<Surface_mesh>  Halfedge_around_target_iterator;
  typedef Iterator_range<Halfedge_around_target_iterator> Halfedge_around_target_range;

  typedef CGAL::Face_around_target_iterator<Surface_mesh>  Face_around_target_iterator;
  typedef Iterator_range<Face_around_target_iterator> Face_around_target_range;

  typedef CGAL::Vertex_around_face_iterator<Surface_mesh>  Vertex_around_face_iterator;
  typedef Iterator_range<Vertex_around_face_iterator> Vertex_around_face_range;

  typedef CGAL::Halfedge_around_face_iterator<Surface_mesh>  Halfedge_around_face_iterator;
  typedef Iterator_range<Halfedge_around_face_iterator> Halfedge_around_face_range;

  typedef CGAL::Face_around_face_iterator<Surface_mesh>  Face_around_face_iterator;
  typedef Iterator_range<Face_around_face_iterator> Face_around_face_range;
#endif

    /// @cond CGAL_BEGIN_END
    /// Start iterator for vertices.
    Vertex_iterator vertices_begin() const
    {
        return Vertex_iterator(Vertex_index(0), this);
    }

    /// End iterator for vertices.
    Vertex_iterator vertices_end() const
    {
        return Vertex_iterator(Vertex_index(num_vertices()), this);
    }
    /// @endcond


    /// returns the iterator range of the vertices of the mesh.
    Vertex_range vertices() const {
      return make_range(vertices_begin(), vertices_end());
    }


    /// @cond CGAL_BEGIN_END
    /// Start iterator for halfedges.
    Halfedge_iterator halfedges_begin() const
    {
        return Halfedge_iterator(Halfedge_index(0), this);
    }

    /// End iterator for halfedges.
    Halfedge_iterator halfedges_end() const
    {
        return Halfedge_iterator(Halfedge_index(num_halfedges()), this);
    }
    /// @endcond


    /// returns the iterator range of the halfedges of the mesh.
    Halfedge_range halfedges() const {
      return make_range(halfedges_begin(), halfedges_end());
    }


    /// @cond CGAL_BEGIN_END
    /// Start iterator for edges.
    Edge_iterator edges_begin() const
    {
        return Edge_iterator(Edge_index(0), this);
    }

    /// End iterator for edges.
    Edge_iterator edges_end() const
    {
        return Edge_iterator(Edge_index(num_edges()), this);
    }
    /// @endcond


    /// returns the iterator range of the edges of the mesh.
    Edge_range edges() const
    {
        return make_range(edges_begin(), edges_end());
    }


    /// @cond CGAL_BEGIN_END
    /// Start iterator for faces.
    Face_iterator faces_begin() const
    {
        return Face_iterator(Face_index(0), this);
    }

    /// End iterator for faces.
    Face_iterator faces_end() const
    {
        return Face_iterator(Face_index(num_faces()), this);
    }
    /// @endcond

    /// returns the iterator range of the faces of the mesh.
    Face_range faces() const {
      return make_range(faces_begin(), faces_end());
    }

#ifndef DOXYGEN_RUNNING
    /// returns the iterator range for vertices around vertex `target(h)`, starting at `source(h)`.
    Vertex_around_target_range vertices_around_target(Halfedge_index h) const
    {
      return CGAL::vertices_around_target(h,*this);
    }

    /// returns the iterator range for incoming halfedges around vertex `target(h)`, starting at `h`.
    Halfedge_around_target_range halfedges_around_target(Halfedge_index h) const
    {
      return CGAL::halfedges_around_target(h,*this);
    }

    /// returns the iterator range for faces around vertex `target(h)`, starting at `face(h)`.
    Face_around_target_range faces_around_target(Halfedge_index h) const
    {
      return CGAL::faces_around_target(h,*this);
    }

    /// returns the iterator range for vertices around face `face(h)`, starting at `target(h)`.
    Vertex_around_face_range vertices_around_face(Halfedge_index h) const
     {
       return CGAL::vertices_around_face(h,*this);
     }

    /// returns the iterator range for halfedges around face `face(h)`, starting at `h`.
    Halfedge_around_face_range halfedges_around_face(Halfedge_index h) const
    {
      return CGAL::halfedges_around_face(h,*this);
    }

    /// returns the iterator range for halfedges around face `face(h)`, starting at `h`.
    Face_around_face_range faces_around_face(Halfedge_index h) const
    {
       return CGAL::faces_around_face(h,*this);
    }

#endif

    ///@}


public:

#ifndef DOXYGEN_RUNNING
    /// \name Circulator Types
    ///
    /// The following circulators enable to iterate through the elements around a face or vertex.
    /// As explained in the \ref SurfaceMeshOrientation "User Manual", we can speak of a
    /// *clockwise* or *counterclockwise*
    /// traversal, by looking at the surface from the right side.
    ///@{

    /// \brief This class circulates clockwise through all
    /// one-ring neighbors of a vertex.
    ///  A model of `BidirectionalCirculator` with value type `Vertex_index`.
    /// \sa `Halfedge_around_target_circulator`, `Face_around_target_circulator`

  typedef CGAL::Vertex_around_target_circulator<Surface_mesh> Vertex_around_target_circulator;



    /// \brief This class circulates clockwise through all incident faces of a vertex.
    ///  A model of `BidirectionalCirculator` with value type `Face_index`.
    /// \sa `Vertex_around_target_circulator`, `Halfedge_around_target_circulator`

  typedef CGAL::Face_around_target_circulator<Surface_mesh> Face_around_target_circulator;


    /// \brief This class circulates clockwise through all halfedges around a vertex that have this vertex as target.
    ///  A model of `BidirectionalCirculator` with value type `Halfedge_index`.
    /// \sa `Vertex_around_target_circulator`, `Halfedge_around_target_circulator`

  typedef CGAL::Halfedge_around_target_circulator<Surface_mesh> Halfedge_around_target_circulator;


    /// \brief This class circulates clockwise through all halfedges around a vertex that have this vertex as source.
    ///  A model of `BidirectionalCirculator` with value type `Halfedge_index`.
    /// \sa `Vertex_around_target_circulator`, `Halfedge_around_target_circulator`

  typedef CGAL::Halfedge_around_source_circulator<Surface_mesh> Halfedge_around_source_circulator;

    /// \brief This class circulates counterclockwise through all vertices around a face.
    ///  A model of `BidirectionalCirculator` with value type `Vertex_index`.

  typedef  CGAL::Vertex_around_face_circulator<Surface_mesh> Vertex_around_face_circulator;


    /// \brief This class circulates counterclockwise through all halfedges around a face.
    ///  A model of `BidirectionalCirculator` with value type `Halfedge_index`.

  typedef  CGAL::Halfedge_around_face_circulator<Surface_mesh> Halfedge_around_face_circulator;

   /// \brief This class circulates counterclockwise through all faces around a face.
   ///  A model of `BidirectionalCirculator` with value type `Face_index`.
   ///  Note that the face index is the same after `operator++`, if the neighboring faces share
   ///  several halfedges.

  typedef  CGAL::Face_around_face_circulator<Surface_mesh> Face_around_face_circulator;
  /// @}
#endif

  /// @cond CGAL_DOCUMENT_INTERNALS
  // typedefs which make it easier to write the partial specialisation of boost::graph_traits

  typedef Vertex_index   vertex_index;
  typedef P                   vertex_property_type;
  typedef Halfedge_index halfedge_index;
  typedef Edge_index     edge_index;
  typedef Face_index     face_index;

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

    /// Move constructor.
    Surface_mesh(Surface_mesh&& sm)
      : vprops_(std::move(sm.vprops_))
      , hprops_(std::move(sm.hprops_))
      , eprops_(std::move(sm.eprops_))
      , fprops_(std::move(sm.fprops_))
      , vconn_(std::move(sm.vconn_))
      , hconn_(std::move(sm.hconn_))
      , fconn_(std::move(sm.fconn_))
      , vremoved_(std::move(sm.vremoved_))
      , eremoved_(std::move(sm.eremoved_))
      , fremoved_(std::move(sm.fremoved_))
      , vpoint_(std::move(sm.vpoint_))
      , removed_vertices_(std::exchange(sm.removed_vertices_, 0))
      , removed_edges_(std::exchange(sm.removed_edges_, 0))
      , removed_faces_(std::exchange(sm.removed_faces_, 0))
      , vertices_freelist_(std::exchange(sm.vertices_freelist_,(std::numeric_limits<size_type>::max)()))
      , edges_freelist_(std::exchange(sm.edges_freelist_,(std::numeric_limits<size_type>::max)()))
      , faces_freelist_(std::exchange(sm.faces_freelist_,(std::numeric_limits<size_type>::max)()))
      , garbage_(std::exchange(sm.garbage_, false))
      , recycle_(std::exchange(sm.recycle_, true))
      , anonymous_property_(std::exchange(sm.anonymous_property_, 0))
    {}

    /// assigns `rhs` to `*this`. Performs a deep copy of all properties.
    Surface_mesh& operator=(const Surface_mesh& rhs);

    /// move assignment
    Surface_mesh& operator=(Surface_mesh&& sm)
    {
      vprops_ = std::move(sm.vprops_);
      hprops_ = std::move(sm.hprops_);
      eprops_ = std::move(sm.eprops_);
      fprops_ = std::move(sm.fprops_);
      vconn_ = std::move(sm.vconn_);
      hconn_ = std::move(sm.hconn_);
      fconn_ = std::move(sm.fconn_);
      vremoved_ = std::move(sm.vremoved_);
      eremoved_ = std::move(sm.eremoved_);
      fremoved_ = std::move(sm.fremoved_);
      vpoint_ = std::move(sm.vpoint_);
      removed_vertices_ = std::exchange(sm.removed_vertices_, 0);
      removed_edges_ = std::exchange(sm.removed_edges_, 0);
      removed_faces_ = std::exchange(sm.removed_faces_, 0);
      vertices_freelist_ = std::exchange(sm.vertices_freelist_, (std::numeric_limits<size_type>::max)());
      edges_freelist_ = std::exchange(sm.edges_freelist_,(std::numeric_limits<size_type>::max)());
      faces_freelist_ = std::exchange(sm.faces_freelist_,(std::numeric_limits<size_type>::max)());
      garbage_ = std::exchange(sm.garbage_, false);
      recycle_ = std::exchange(sm.recycle_, true);
      anonymous_property_ = std::exchange(sm.anonymous_property_, 0);
      return *this;
    }

    /// assigns `rhs` to `*this`. Does not copy custom properties.
    Surface_mesh& assign(const Surface_mesh& rhs);

    ///@}

public:

    /// \name Adding Vertices, Edges, and Faces
    ///@{

   /// adds a new vertex, and resizes vertex properties if necessary.
    Vertex_index add_vertex()
    {
      size_type inf = (std::numeric_limits<size_type>::max)();
      if(recycle_ && (vertices_freelist_ != inf)){
        size_type idx = vertices_freelist_;
        vertices_freelist_ = (size_type)vconn_[Vertex_index(vertices_freelist_)].halfedge_;
        --removed_vertices_;
        vremoved_[Vertex_index(idx)] = false;
        vprops_.reset(Vertex_index(idx));
        return Vertex_index(idx);
      } else {
        vprops_.push_back();
        return Vertex_index(num_vertices()-1);
      }
    }

    /// adds a new vertex, resizes vertex properties if necessary,
    /// and sets the \em point property to `p`.
    /// \note Several vertices may have the same point property.
    Vertex_index add_vertex(const Point& p)
    {
        Vertex_index v = add_vertex();
        vpoint_[v] = p;
        return v;
    }



public:

    /// adds a new edge, and resizes edge and halfedge properties if necessary.
    Halfedge_index add_edge()
    {
      Halfedge_index h0, h1;
      size_type inf = (std::numeric_limits<size_type>::max)();
      if(recycle_ && (edges_freelist_ != inf)){
        size_type idx = edges_freelist_;
        edges_freelist_ = (size_type)hconn_[Halfedge_index(edges_freelist_)].next_halfedge_;
        --removed_edges_;
        eremoved_[Edge_index(Halfedge_index(idx))] = false;
        hprops_.reset(Halfedge_index(idx));
        hprops_.reset(opposite(Halfedge_index(idx)));
        eprops_.reset(Edge_index(Halfedge_index(idx)));
        return Halfedge_index(idx);
      } else {
        eprops_.push_back();
        hprops_.push_back();
        hprops_.push_back();

        return Halfedge_index(num_halfedges()-2);
      }
    }

    /// adds two opposite halfedges, and resizes edge and halfedge properties if necessary.
    /// Sets the targets of the halfedge to the given vertices, but does not modify the halfedge
    /// associated to the vertices.
    /// \note The function does not check whether there is already an edge between the vertices.
    /// \returns the halfedge with `v1` as target

    Halfedge_index add_edge(Vertex_index v0, Vertex_index v1)
    {
        CGAL_assertion(v0 != v1);
        Halfedge_index h = add_edge();

        set_target(h, v1);
        set_target(opposite(h), v0);

        return h;
    }

    /// adds a new face, and resizes face properties if necessary.
    Face_index add_face()
    {
      size_type inf = (std::numeric_limits<size_type>::max)();
      if(recycle_ && (faces_freelist_ != inf)){
        size_type idx = faces_freelist_;
        faces_freelist_ = (size_type)fconn_[Face_index(faces_freelist_)].halfedge_;
        --removed_faces_;
        fprops_.reset(Face_index(idx));
        fremoved_[Face_index(idx)] = false;
        return Face_index(idx);
      } else {
        fprops_.push_back();
        return Face_index(num_faces()-1);
      }
    }

    /// if possible, adds a new face with vertices from a range with value type `Vertex_index`.
    /// The function adds halfedges between successive vertices if they are not yet indicent to halfedges,
    /// or updates the connectivity of halfedges already in place.
    /// Resizes halfedge, edge, and face properties if necessary.
    /// \returns the face index of the added face, or `Surface_mesh::null_face()` if the face could not be added.
    template <typename Range>
    Face_index add_face(const Range& vertices);


    /// adds a new triangle connecting vertices `v0`, `v1`, `v2`.
    /// \returns the face index of the added face, or `Surface_mesh::null_face()` if the face could not be added.
    Face_index add_face(Vertex_index v0, Vertex_index v1, Vertex_index v2)
    {
        boost::array<Vertex_index, 3>
            v = {{v0, v1, v2}};
        return add_face(v);
    }

    /// adds a new quad connecting vertices `v0`, `v1`, `v2`, `v3`.
    /// \returns the face index of the added face, or `Surface_mesh::null_face()` if the face could not be added.
    Face_index add_face(Vertex_index v0, Vertex_index v1, Vertex_index v2, Vertex_index v3)
    {
        boost::array<Vertex_index, 4>
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
    void remove_vertex(Vertex_index v)
    {
        vremoved_[v] = true; ++removed_vertices_; garbage_ = true;
        vconn_[v].halfedge_ = Halfedge_index(vertices_freelist_);
        vertices_freelist_ = (size_type)v;
    }

    /// removes the two halfedges corresponding to `e` from the halfedge data structure without
    /// adjusting anything.
    void remove_edge(Edge_index e)
    {
        eremoved_[e] = true; ++removed_edges_; garbage_ = true;
        hconn_[Halfedge_index((size_type)e << 1)].next_halfedge_ = Halfedge_index(edges_freelist_ );
        edges_freelist_ = ((size_type)e << 1);
    }

    /// removes  face `f` from the halfedge data structure without
    /// adjusting anything.

    void remove_face(Face_index f)
    {
        fremoved_[f] = true; ++removed_faces_; garbage_ = true;
        fconn_[f].halfedge_ = Halfedge_index(faces_freelist_);
        faces_freelist_ = (size_type)f;
    }


    ///@}


    /// \name Memory Management
    ///
    /// Functions to check the number of elements, the amount of space
    /// allocated for elements, and to clear the structure.
    ///@{

  /// returns the number of vertices in the mesh.
  size_type number_of_vertices() const
  {
    return num_vertices() - number_of_removed_vertices();
  }

  /// returns the number of halfedges in the mesh.
  size_type number_of_halfedges() const
  {
    return num_halfedges() - number_of_removed_halfedges();
  }

  /// returns the number of edges in the mesh.
  size_type number_of_edges() const
  {
    return num_edges() - number_of_removed_edges();
  }

  /// returns the number of faces in the mesh.
  size_type number_of_faces() const
  {
    return num_faces() - number_of_removed_faces();
  }

    /// returns `true` iff the mesh is empty, i.e., has no vertices, halfedges and faces.
    bool is_empty() const
  {
    return ( num_vertices() == number_of_removed_vertices()
             && num_halfedges() == number_of_removed_halfedges()
             && num_faces() == number_of_removed_faces());
  }

  /// removes all vertices, halfedge, edges and faces. Collects garbage but keeps all property maps.
  void clear_without_removing_property_maps();

    /// removes all vertices, halfedge, edges and faces. Collects garbage and removes all property maps added by a call to `add_property_map()` for all simplex types.
    ///
    /// After calling this method, the object is the same as a newly constructed object. The additional property maps are also removed and must thus be re-added if needed.
    void clear();


    /// reserves space for vertices, halfedges, edges, faces, and their currently
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

      void resize(size_type nvertices,
                 size_type nedges,
                 size_type nfaces )
    {
        vprops_.resize(nvertices);
        hprops_.resize(2*nedges);
        eprops_.resize(nedges);
        fprops_.resize(nfaces);
    }

  /// copies the simplices from `other`, and copies values of
  /// properties that already exist under the same name in `*this`.
  /// In case `*this` has a property that does not exist in `other`
  /// the copied simplices get the default value of the property.
  bool join(const Surface_mesh& other)
  {
    // increase capacity
    const size_type nv = num_vertices(), nh = num_halfedges(), nf = num_faces();
    resize(num_vertices()+  other.num_vertices(),
            num_edges()+  other.num_edges(),
            num_faces()+  other.num_faces());

    // append properties in the free space created by resize
    vprops_.transfer(other.vprops_);
    hprops_.transfer(other.hprops_);
    fprops_.transfer(other.fprops_);
    eprops_.transfer(other.eprops_);

    // translate halfedge index in vertex -> halfedge
    for(size_type i = nv; i < nv+other.num_vertices(); i++){
      Vertex_index vi(i);
      if(vconn_[vi].halfedge_ != null_halfedge()){
        vconn_[vi].halfedge_ = Halfedge_index(size_type(vconn_[vi].halfedge_)+nh);
      }
    }
    // translate halfedge index in face -> halfedge
    for(size_type i = nf; i < nf+other.num_faces(); i++){
      Face_index fi(i);
      if(fconn_[fi].halfedge_ != null_halfedge()){
        fconn_[fi].halfedge_ = Halfedge_index(size_type(fconn_[fi].halfedge_)+nh);
      }
    }
    // translate indices in halfedge -> face, halfedge -> target, halfedge -> prev, and halfedge -> next
    for(size_type i = nh; i < nh+other.num_halfedges(); i++){
      Halfedge_index hi(i);
      if(hconn_[hi].face_ != null_face()){
        hconn_[hi].face_ = Face_index(size_type(hconn_[hi].face_)+nf);
      }
      if( hconn_[hi].vertex_ != null_vertex()){
        hconn_[hi].vertex_ = Vertex_index(size_type(hconn_[hi].vertex_)+nv);
      }
      if(hconn_[hi].next_halfedge_ != null_halfedge()){
        hconn_[hi].next_halfedge_ = Halfedge_index(size_type(hconn_[hi].next_halfedge_)+nh);
      }
      if(hconn_[hi].prev_halfedge_ != null_halfedge()){
        hconn_[hi].prev_halfedge_ = Halfedge_index(size_type(hconn_[hi].prev_halfedge_)+nh);
      }
    }
    size_type inf_value = (std::numeric_limits<size_type>::max)();

    // merge vertex free list
    if(other.vertices_freelist_ != inf_value){
      Vertex_index vi(nv+other.vertices_freelist_);
      Halfedge_index inf((std::numeric_limits<size_type>::max)());
      // correct the indices in the linked list of free vertices copied (due to vconn_ translation)
      while(vconn_[vi].halfedge_ != inf){
        Vertex_index corrected_vi = Vertex_index(size_type(vconn_[vi].halfedge_)+nv-nh);
        vconn_[vi].halfedge_ = Halfedge_index(corrected_vi);
        vi = corrected_vi;
      }
      // append the vertex free linked list of `this` to the copy of `other`
      vconn_[vi].halfedge_ = Halfedge_index(vertices_freelist_);
      // update the begin of the vertex free linked list
      vertices_freelist_ = nv + other.vertices_freelist_;
    }
    // merge face free list
    if(other.faces_freelist_ != inf_value){
      Face_index fi(nf+other.faces_freelist_);
      Halfedge_index inf((std::numeric_limits<size_type>::max)());
      // correct the indices in the linked list of free faces copied (due to fconn_ translation)
      while(fconn_[fi].halfedge_ != inf){
        Face_index corrected_fi = Face_index(size_type(fconn_[fi].halfedge_)+nf-nh);
        fconn_[fi].halfedge_ = Halfedge_index(corrected_fi);
        fi = corrected_fi;
      }
      // append the face free linked list of `this` to the copy of `other`
      fconn_[fi].halfedge_ = Halfedge_index(faces_freelist_);
      // update the begin of the face free linked list
      faces_freelist_ = nf + other.faces_freelist_;
    }
    // merge edge free list
    if(other.edges_freelist_ != inf_value){
      Halfedge_index hi(nh+other.edges_freelist_);
      Halfedge_index inf((std::numeric_limits<size_type>::max)());
      while(hconn_[hi].next_halfedge_ != inf){
        hi = hconn_[hi].next_halfedge_;
      }
      // append the halfedge free linked list of `this` to the copy of `other`
      hconn_[hi].next_halfedge_ = Halfedge_index(edges_freelist_);
      // update the begin of the halfedge free linked list
      edges_freelist_ = nh + other.edges_freelist_;
    }
    // update garbage infos
    garbage_ = garbage_ || other.garbage_;
    removed_vertices_ += other.removed_vertices_;
    removed_edges_ += other.removed_edges_;
    removed_faces_ += other.removed_faces_;
    return true;
  }

    ///@}


    /// \name Garbage Collection
    ///
    /// While removing elements only marks them as removed
    /// garbage collection really removes them.
    /// The API in this section allows to check whether
    /// an element is removed, to get the number of
    /// removed elements, and to collect garbage.
    /// The number of elements together with the number of  removed elements is
    /// an upperbound on the index, and is needed
    /// by algorithms that temporarily store a
    /// property in a vector of the appropriate size.
    /// Note however that by garbage collecting elements get new indices.
    /// In case you store indices in an auxiliary data structure
    /// or in a property these indices are potentially no longer
    /// referring to the right elements.
    /// When adding elements, by default elements that are marked as removed
    /// are recycled.

    ///@{
#ifndef DOXYGEN_RUNNING
   /// returns the number of used and removed vertices in the mesh.
    size_type num_vertices() const { return (size_type) vprops_.size(); }

    /// returns the number of used and removed halfedges in the mesh.
    size_type num_halfedges() const { return (size_type) hprops_.size(); }

    /// returns the number of used and removed edges in the mesh.
    size_type num_edges() const { return (size_type) eprops_.size(); }

    /// returns the number of used and removed faces in the mesh.
    size_type num_faces() const { return (size_type) fprops_.size(); }

#endif

    /// returns the number of vertices in the mesh which are marked removed.
    size_type number_of_removed_vertices() const { return removed_vertices_; }

    /// returns the number of halfedges in the mesh which are marked removed.
    size_type number_of_removed_halfedges() const { return 2*removed_edges_; }

    /// returns the number of edges in the mesh which are marked removed.
    size_type number_of_removed_edges() const { return removed_edges_; }

    /// returns the number offaces in the mesh which are marked removed.
    size_type number_of_removed_faces() const { return removed_faces_; }



    /// returns whether vertex `v` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Vertex_index v) const
    {
        return vremoved_[v];
    }
    /// returns whether halfedge `h` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Halfedge_index h) const
    {
        return eremoved_[edge(h)];
    }
    /// returns whether edge `e` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Edge_index e) const
    {
        return eremoved_[e];
    }
    /// returns whether face `f` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Face_index f) const
    {
        return fremoved_[f];
    }

    /// checks if any vertices, halfedges, edges, or faces are marked as removed.
    /// \sa collect_garbage
    bool has_garbage() const { return garbage_; }

    /// really removes vertices, halfedges, edges, and faces which are marked removed.
    /// \sa `has_garbage()`
    /// \attention By garbage collecting elements get new indices.
    /// In case you store indices in an auxiliary data structure
    /// or in a property these indices are potentially no longer
    /// referring to the right elements.
    void collect_garbage();

    //undocumented convenience function that allows to get old-index->new-index information
    template <typename Visitor>
    void collect_garbage(Visitor& visitor);

    /// controls the recycling or not of simplices previously marked as removed
    /// upon addition of new elements.
    /// When set to `true` (default value), new elements are first picked in the garbage (if any)
    /// while if set to `false` only new elements are created.
    void set_recycle_garbage(bool b);

    /// Getter
    bool does_recycle_garbage() const;

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
    bool has_valid_index(Vertex_index v) const
    {
      return ((size_type)v < num_vertices());
    }

    /// returns whether the index of halfedge `h` is valid, that is within the current array bounds.
    bool has_valid_index(Halfedge_index h) const
    {
      return ((size_type)h < num_halfedges());
    }
    /// returns whether the index of edge `e` is valid, that is within the current array bounds.
    bool has_valid_index(Edge_index e) const
    {
      return ((size_type)e < num_edges());
    }
    /// returns whether the index of face `f` is valid, that is within the current array bounds.
    bool has_valid_index(Face_index f) const
    {
        return ((size_type)f < num_faces());
    }

    /// @}
    /// @endcond

    /// \name Validity Checks
    ///
    /// Functions in this group perform checks for structural
    /// consistency of a complete surface mesh, or an individual element.
    /// They are expensive and should only be used in debug configurations.

    ///@{

    /// perform an expensive validity check on the data structure and
    /// print found errors to `std::cerr` when `verbose == true`.
    bool is_valid(bool verbose = false) const
    {
        bool valid = true;
        size_type vcount = 0, hcount = 0, fcount = 0;
        for(Halfedge_iterator it = halfedges_begin(); it != halfedges_end(); ++it) {
            ++hcount;
            valid = valid && next(*it).is_valid();
            valid = valid && opposite(*it).is_valid();
            if(!valid) {
                if (verbose)
                  std::cerr << "Integrity of halfedge " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (opposite(*it) != *it);
            valid = valid && (opposite(opposite(*it)) == *it);
            if(!valid) {
              if (verbose)
                std::cerr << "Integrity of opposite halfedge of " << *it << " corrupted."  << std::endl;
              break;
            }

            valid = valid && (next(prev(*it)) == *it);
            if(!valid) {
                if (verbose)
                  std::cerr << "Integrity of previous halfedge of " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (prev(next(*it)) == *it);
            if(!valid) {
                if (verbose)
                  std::cerr << "Integrity of next halfedge of " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && target(*it).is_valid();
            if(!valid) {
                if (verbose)
                  std::cerr << "Integrity of vertex of halfedge " << *it << " corrupted."  << std::endl;
                break;
            }

            valid = valid && (target(*it) == target(opposite(next(*it))));
            if(!valid) {
                if (verbose)
                  std::cerr << "Halfedge vertex of next opposite is not the same for " << *it << "."  << std::endl;
                break;
            }
        }

        for(Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it) {
          ++vcount;
            if(halfedge(*it).is_valid()) {
                // not an isolated vertex
                valid = valid && (target(halfedge(*it)) == *it);
                if(!valid) {
                    if (verbose)
                      std::cerr << "Halfedge of " << *it << " is not an incoming halfedge." << std::endl;
                    break;
                }
            }
        }
        for(Face_iterator it = faces_begin(); it != faces_end(); ++it) {
          ++fcount;
        }

        valid = valid && (vcount == number_of_vertices());
        if(!valid && verbose){
          std::cerr << "#vertices: iterated: " << vcount << " vs number_of_vertices(): " << number_of_vertices()<< std::endl;
        }

        valid = valid && (hcount == number_of_halfedges());
        if(!valid && verbose){
          std::cerr << "#halfedges: iterated: " << hcount << " vs number_of_halfedges(): " << number_of_halfedges()<< std::endl;
        }

        valid = valid && (fcount == number_of_faces());
        if(!valid && verbose){
          std::cerr << "#faces: iterated: " << fcount << " vs number_of_faces(): " << number_of_faces()<< std::endl;
        }

        size_type inf = (std::numeric_limits<size_type>::max)();
        size_type vfl = vertices_freelist_;
        size_type rv = 0;
        while(vfl != inf){
          vfl = (size_type)vconn_[Vertex_index(vfl)].halfedge_;
          rv++;
        }
        valid = valid && ( rv == removed_vertices_ );


        size_type efl = edges_freelist_;
        size_type re = 0;
        while(efl != inf){
          efl = (size_type)hconn_[Halfedge_index(efl)].next_halfedge_;
          re++;
        }
        valid = valid && ( re == removed_edges_ );

        size_type ffl = faces_freelist_;
        size_type rf = 0;
        while(ffl != inf){
          ffl = (size_type)fconn_[Face_index(ffl)].halfedge_;
          rf++;
        }
        valid = valid && ( rf == removed_faces_ );

        return valid;
    }

    /// performs a validity check on a single vertex.
    bool is_valid(Vertex_index v,
                  bool verbose = false) const
    {
        Verbose_ostream verr(verbose);

        if(!has_valid_index(v))
        {
          verr << "Vertex has invalid index: " << (size_type)v << std::endl;
          return false;
        }

        Halfedge_index h = vconn_[v].halfedge_;
        if(h != null_halfedge() && (!has_valid_index(h) || is_removed(h))) {
          verr << "Vertex connectivity halfedge error: Vertex " << (size_type)v
               << " with " << (size_type)h << std::endl;
          return false;
        }
        return true;
    }

    /// performs a validity check on a single halfedge.
    bool is_valid(Halfedge_index h,
                  bool verbose = false) const
    {
        Verbose_ostream verr(verbose);

        if(!has_valid_index(h))
        {
          verr << "Halfedge has invalid index: " << (size_type)h << std::endl;
          return false;
        }

        Face_index f = hconn_[h].face_;
        Vertex_index v = hconn_[h].vertex_;
        Halfedge_index hn = hconn_[h].next_halfedge_;
        Halfedge_index hp = hconn_[h].prev_halfedge_;

        bool valid = true;
        // don't validate the face if this is a border halfedge
        if(!is_border(h)) {
            if(!has_valid_index(f) || is_removed(f)) {
                verr << "Halfedge connectivity error: Face "
                     << (!has_valid_index(f) ? "invalid" : "removed")
                     << " in " << (size_type)h << std::endl;
                valid = false;
            }
        }

        if(!has_valid_index(v) || is_removed(v)) {
            verr << "Halfedge connectivity error: Vertex "
                 << (!has_valid_index(v) ? "invalid" : "removed")
                 << " in " << (size_type)h << std::endl;
            valid = false;
        }

        if(!has_valid_index(hn) || is_removed(hn)) {
            verr << "Halfedge connectivity error: hnext "
                 << (!has_valid_index(hn) ? "invalid" : "removed")
                 << " in " << (size_type)h << std::endl;
            valid = false;
        }
        if(!has_valid_index(hp) || is_removed(hp)) {
            verr << "Halfedge connectivity error: hprev "
                 << (!has_valid_index(hp) ? "invalid" : "removed")
                 << " in " << (size_type)h << std::endl;
            valid = false;
        }
        return valid;
    }


    /// performs a validity check on a single edge.
    bool is_valid(Edge_index e,
                  bool verbose = false) const
    {
      Verbose_ostream verr(verbose);

      if(!has_valid_index(e))
      {
        verr << "Edge has invalid index: " << (size_type)e << std::endl;
        return false;
      }

      Halfedge_index h = halfedge(e);
      return is_valid(h, verbose) && is_valid(opposite(h), verbose);
    }


    /// performs a validity check on a single face.
    bool is_valid(Face_index f,
                  bool verbose = false) const
    {
        Verbose_ostream verr(verbose);

        if(!has_valid_index(f))
        {
          verr << "Face has invalid index: " << (size_type)f << std::endl;
          return false;
        }

        Halfedge_index h = fconn_[f].halfedge_;
        if(!has_valid_index(h) || is_removed(h)) {
          verr << "Face connectivity halfedge error: Face " << (size_type)f
               << " with " << (size_type)h << std::endl;
          return false;
        }
        return true;
    }

    ///@}



    /// \name Low-Level Connectivity
    ///@{

    /// returns the vertex the halfedge `h` points to.
    Vertex_index target(Halfedge_index h) const
    {
        return hconn_[h].vertex_;
    }

    /// sets the vertex the halfedge `h` points to to `v`.
    void set_target(Halfedge_index h, Vertex_index v)
    {
        hconn_[h].vertex_ = v;
    }

    /// returns the face incident to halfedge `h`.
    Face_index face(Halfedge_index h) const
    {
        return hconn_[h].face_;
    }

    /// sets the incident face to halfedge `h` to `f`.
    void set_face(Halfedge_index h, Face_index f)
    {
        hconn_[h].face_ = f;
    }

    /// returns the next halfedge within the incident face.
    Halfedge_index next(Halfedge_index h) const
    {
        return hconn_[h].next_halfedge_;
    }

    /// returns the previous halfedge within the incident face.
    Halfedge_index prev(Halfedge_index h) const
    {
        return hconn_[h].prev_halfedge_;
    }

    /// @cond CGAL_DOCUMENT_INTERNALS
    // sets the next halfedge of `h` within the face to `nh`.
    void set_next_only(Halfedge_index h, Halfedge_index nh)
    {
      hconn_[h].next_halfedge_ = nh;
    }

    // sets previous halfedge of `h` to `nh`.
    void set_prev_only(Halfedge_index h, Halfedge_index nh)
    {
      if(h != null_halfedge()){
        hconn_[h].prev_halfedge_ = nh;
      }
    }
    /// @endcond

    /// sets the next halfedge of `h` within the face to `nh` and
    /// the previous halfedge of `nh` to `h`.
    void set_next(Halfedge_index h, Halfedge_index nh)
    {
      set_next_only(h, nh);
      set_prev_only(nh, h);
    }

    /// returns an incoming halfedge of vertex `v`.
    /// If `v` is a border vertex this will be a border halfedge.
    /// \invariant `target(halfedge(v)) == v`
    Halfedge_index halfedge(Vertex_index v) const
    {
        return vconn_[v].halfedge_;
    }

    /// sets the incoming halfedge of vertex `v` to `h`.
    void set_halfedge(Vertex_index v, Halfedge_index h)
    {
        vconn_[v].halfedge_ = h;
    }


    /// returns a halfedge of face `f`.
    Halfedge_index halfedge(Face_index f) const
    {
        return fconn_[f].halfedge_;
    }

    /// sets the halfedge of face `f` to `h`.
    void set_halfedge(Face_index f, Halfedge_index h)
    {
        fconn_[f].halfedge_ = h;
    }

    /// returns the opposite halfedge of `h`. Note that there is no function `set_opposite()`.
    Halfedge_index opposite(Halfedge_index h) const
    {
        return Halfedge_index(((size_type)h & 1) ? (size_type)h-1 : (size_type)h+1);
    }

    ///@}

    /// \name Low-Level Connectivity Convenience Functions
    ///@{

    /// returns the vertex the halfedge `h` emanates from.
    Vertex_index source(Halfedge_index h) const
    {
        return target(opposite(h));
    }

    /// returns `opposite(next(h))`, that is the next halfedge \ref SurfaceMeshOrientation
    /// "clockwise" around the target vertex of `h`.
    Halfedge_index next_around_target(Halfedge_index h) const
    {
        return opposite(next(h));
    }

    /// returns `prev(opposite(h))`, that is the previous halfedge \ref SurfaceMeshOrientation
    /// "clockwise" around the target vertex of `h`.
    Halfedge_index prev_around_target(Halfedge_index h) const
    {
        return prev(opposite(h));
    }

    /// returns `next(opposite(h))`, that is the next halfedge \ref SurfaceMeshOrientation
    /// "clockwise" around the source vertex of `h`.
    Halfedge_index next_around_source(Halfedge_index h) const
    {
        return next(opposite(h));
    }

    /// returns `opposite(prev(h))`, that is the previous halfedge \ref SurfaceMeshOrientation
    /// "clockwise" around the source vertex of `h`.
    Halfedge_index prev_around_source(Halfedge_index h) const
    {
        return opposite(prev(h));
    }

    /// returns the i'th vertex of edge `e`, for `i=0` or `1`.
    Vertex_index vertex(Edge_index e, unsigned int i) const
    {
        CGAL_assertion(i<=1);
        return target(halfedge(e, i));
    }

    /// finds a halfedge between two vertices. Returns a default constructed
    /// `Halfedge_index`, if  `source` and  `target` are not connected.
    Halfedge_index halfedge(Vertex_index source, Vertex_index target) const;

    ///@}


    /// \name Switching between Halfedges and Edges
    ///@{

    /// returns the edge that contains halfedge `h` as one of its two halfedges.
    Edge_index edge(Halfedge_index h) const
    {
        return Edge_index(h);
    }

    /// returns the halfedge corresponding to the edge `e`.
    Halfedge_index halfedge(Edge_index e) const
    {
        return Halfedge_index(e.halfedge());
    }

    /// returns the i'th halfedge of edge `e`, for `i=0` or `1`.
    Halfedge_index halfedge(Edge_index e, unsigned int i) const
    {
        CGAL_assertion(i<=1);
        return Halfedge_index(((size_type)e << 1) + i);
    }

    ///@}


    /// \name Degree Functions
    ///@{

    /// returns the number of incident halfedges of vertex `v`.
    size_type degree(Vertex_index v) const;

    /// returns the number of incident halfedges of face `f`.
    size_type degree(Face_index f) const;

    ///@}



    /// \name Borders
    ///
    ///  A halfedge, or edge is on the border of a surface mesh
    /// if it is incident to a `null_face()`.  A vertex is on a border
    /// if it is isolated or incident to a border halfedge. While for a halfedge and
    /// edge this is a constant time operation, for a vertex it means
    /// to look at all incident halfedges.  If algorithms operating on a
    /// surface mesh maintain that the halfedge associated to a border vertex is
    /// a border halfedge, this is a constant time operation too.
    /// This section provides functions to check if an element is on a
    /// border and to change the halfedge associated to a border vertex.
    ///@{

    /// returns whether `v` is a border vertex.
    /// \cgalAdvancedBegin
    /// With the default value for
    /// `check_all_incident_halfedges` the function iteratates over the incident halfedges.
    /// With `check_all_incident_halfedges == false` the function returns `true`, if the incident
    /// halfedge associated to vertex `v` is a border halfedge, or if the vertex is isolated.
    /// \cgalAdvancedEnd
    /// \attention If the data contained in the `Surface_mesh` is not a 2-manifold, then
    /// this operation is not guaranteed to return the right result.
  bool is_border(Vertex_index v, bool check_all_incident_halfedges = true) const
    {
        Halfedge_index h(halfedge(v));
        if (h == null_halfedge()){
          return true;
        }
        if(check_all_incident_halfedges){
          Halfedge_around_target_circulator hatc(h,*this), done(hatc);
          do {
            if(is_border(*hatc)){
              return true;
            }
          }while(++hatc != done);
          return false;
        }
        return is_border(h);
    }

    /// returns whether `h` is a border halfege, that is if its incident face is `sm.null_face()`.
    bool is_border(Halfedge_index h) const
    {
        return !face(h).is_valid();
    }


    /// returns whether `e` is a border edge, i.e., if any
    /// of its two halfedges is a border halfedge.
    bool is_border(Edge_index e) const
    {
      return is_border(e.halfedge()) || is_border(opposite(e.halfedge()));
    }

  /// iterates over the incident halfedges and sets the incident halfedge
  /// associated to vertex `v` to a border halfedge and returns `true` if it exists.
  bool set_vertex_halfedge_to_border_halfedge(Vertex_index v)
  {
    if(halfedge(v) == null_halfedge()){
      return false;
    }
    Halfedge_around_target_circulator hatc(halfedge(v),*this), done(hatc);
    do {
      if(is_border(*hatc)){
        set_halfedge(v,*hatc);
        return true;
      }
    }while(++hatc != done);
    return false;
  }

  /// applies `set_vertex_halfedge_to_border_halfedge(Vertex_index)` on all vertices
  /// around the face associated to `h`.
  void set_vertex_halfedge_to_border_halfedge(Halfedge_index h)
  {
    if(is_border(h)){
      Halfedge_around_face_circulator hafc(h,*this),done(hafc);
      do {
        set_halfedge(target(*hafc),*hafc);
      }while(++hafc != done);
    } else {
       Vertex_around_face_circulator vafc(h,*this),done(vafc);
      do {
        set_vertex_halfedge_to_border_halfedge(*vafc);
      }while(++vafc != done);
    }
  }

  /// applies `set_vertex_halfedge_to_border_halfedge(Vertex_index)` on all vertices
  /// of the surface mesh.
  void set_vertex_halfedge_to_border_halfedge()
  {
    for(Halfedge_index h : halfedges()){
      if(is_border(h)){
          set_halfedge(target(h),h);
        }
    }
  }


    /// returns whether `v` is isolated, i.e., incident to `Surface_mesh::null_halfedge()`.
    bool is_isolated(Vertex_index v) const
    {
        return !halfedge(v).is_valid();
    }

    ///@}


private: //--------------------------------------------------- property handling

  // Property_selector maps an index type to a property_container, the
  // dummy is necessary to make it a partial specialization (full
  // specializations are only allowed at namespace scope).
  template<typename, bool = true>
  struct Property_selector {};

  template<bool dummy>
  struct Property_selector<typename CGAL::Surface_mesh<P>::Vertex_index, dummy> {
    CGAL::Surface_mesh<P>* m_;
    Property_selector(CGAL::Surface_mesh<P>* m) : m_(m) {}
    Properties::Property_container<Self,
                                   typename CGAL::Surface_mesh<P>::Vertex_index>&
    operator()() { return m_->vprops_; }
    void resize_property_array() { m_->vprops_.resize_property_array(3); }
  };
  template<bool dummy>
  struct Property_selector<typename CGAL::Surface_mesh<P>::Halfedge_index, dummy> {
    CGAL::Surface_mesh<P>* m_;
    Property_selector(CGAL::Surface_mesh<P>* m) : m_(m) {}
    Properties::Property_container<Self,
                                   typename CGAL::Surface_mesh<P>::Halfedge_index>&
    operator()() { return m_->hprops_; }
    void resize_property_array() { m_->hprops_.resize_property_array(1); }
  };
  template<bool dummy>
  struct Property_selector<typename CGAL::Surface_mesh<P>::Edge_index, dummy> {
    CGAL::Surface_mesh<P>* m_;
    Property_selector(CGAL::Surface_mesh<P>* m) : m_(m) {}
    Properties::Property_container<Self,
                                   typename CGAL::Surface_mesh<P>::Edge_index>&
    operator()() { return m_->eprops_; }
    void resize_property_array() { m_->eprops_.resize_property_array(1); }
  };
  template<bool dummy>
  struct Property_selector<typename CGAL::Surface_mesh<P>::Face_index, dummy> {
    CGAL::Surface_mesh<P>* m_;
    Property_selector(CGAL::Surface_mesh<P>* m) : m_(m) {}
    Properties::Property_container<Self,
                                   typename CGAL::Surface_mesh<P>::Face_index>&
    operator()() { return m_->fprops_; }
    void resize_property_array() { m_->fprops_.resize_property_array(2); }
  };

    public:


 /*! \name Property Handling

 A `Properties::Property_map<I,T>` allows to associate properties of type `T` to a vertex, halfdge, edge, or face index type I.
 Properties can be added, and looked up with a string, and they can be removed at runtime.
 The \em point property of type `P` is associated to the string "v:point".

    */
    ///@{

  /// Model of `LvaluePropertyMap` with `I` as key type and `T` as value type, where `I`
  /// is either a vertex, halfedge, edge, or face index type.
#ifdef DOXYGEN_RUNNING
  template <class I, class T>
  using Property_map = unspecified_type;

#else


#endif

    /// adds a property map named `name` with value type `T` and default `t`
    /// for index type `I`. Returns the property map together with a Boolean
    /// that is `true` if a new map was created. In case it already exists
    /// the existing map together with `false` is returned.


    template<class I, class T>
    std::pair<Property_map<I, T>, bool>
    add_property_map(std::string name=std::string(), const T t=T()) {
      if(name.empty()){
        std::ostringstream oss;
        oss << "anonymous-property-" << anonymous_property_++;
        name = std::string(oss.str());
      }
      return Property_selector<I>(this)().template add<T>(name, t);
    }

    /// returns a property map named `name` with key type `I` and value type `T`,
    /// and a Boolean that is `true` if the property exists.
    /// In case it does not exist the Boolean is `false` and the behavior of
    /// the property map is undefined.
    template <class I, class T>
    std::pair<Property_map<I, T>,bool> property_map(const std::string& name) const
    {
      return Property_selector<I>(const_cast<Surface_mesh*>(this))().template get<T>(name);
    }


    /// removes property map `p`. The memory allocated for that property map is
    /// freed.
    template<class I, class T>
    void remove_property_map(Property_map<I, T>& p)
    {
      (Property_selector<I>(this)()).template remove<T>(p);
    }


    /// removes all property maps for index type `I` added by a call to `add_property_map<I>()`.
    /// The memory allocated for those property maps is freed.
    template<class I>
    void remove_property_maps()
    {
        Property_selector<I>(this).resize_property_array();
    }

    /// removes all property maps for all index types added by a call to `add_property_map()`.
    /// The memory allocated for those property maps is freed.
    void remove_all_property_maps()
    {
        remove_property_maps<Vertex_index>();
        remove_property_maps<Face_index>();
        remove_property_maps<Edge_index>();
        remove_property_maps<Halfedge_index>();
    }

    /// @cond CGAL_DOCUMENT_INTERNALS
    /// returns the std::type_info of the value type of the
    /// property identified by `name`.  `typeid(void)` if `name`
    /// does not identify any property.
    ///
    /// @tparam I The key type of the property.

    template<class I>
    const std::type_info& property_type(const std::string& name)
    {
      return Property_selector<I>(this)().get_type(name);
    }
    /// @endcond

    /// returns a vector with all strings that describe properties with the key type `I`.
    /// @tparam I The key type of the properties.
    template<class I>
    std::vector<std::string> properties() const
    {
      return Property_selector<I>(const_cast<Self*>(this))().properties();
    }

    /// returns the property for the string "v:point".
    Property_map<Vertex_index, Point>
    points() const { return vpoint_; }

    Property_map<Vertex_index, Point>&
    points() { return vpoint_; }

    /// returns the point associated to vertex `v`.
    const Point&
    point(Vertex_index v) const { return vpoint_[v]; }

    /// returns the point associated to vertex `v`.
    Point&
    point(Vertex_index v) { return vpoint_[v]; }

    /// @cond CGAL_DOCUMENT_INTERNALS
    /// prints property statistics to the stream `out`. The output is human-readable but
    /// not machine-friendly.
    ///
    void property_stats(std::ostream& out = std::cout) const;
    /// @endcond
    ///@}


 /// \name Null Elements
    ///@{

  /// returns `Vertex_index(std::numeric_limits<size_type>::%max())`.
  static Vertex_index null_vertex()
  {
    return vertex_index((std::numeric_limits<size_type>::max)());
  }

  /// returns `Edge_index(std::numeric_limits<size_type>::%max())`.
  static Edge_index null_edge()
  {
    return edge_index((std::numeric_limits<size_type>::max)());
  }
  /// returns `Halfedge_index(std::numeric_limits<size_type>::%max())`.
  static Halfedge_index null_halfedge()
  {
    return halfedge_index((std::numeric_limits<size_type>::max)());
  }
  /// returns `Face_index(std::numeric_limits<size_type>::%max())`.
  static Face_index null_face()
  {
    return face_index((std::numeric_limits<size_type>::max)());
  }
  /// @}

#if defined(CGAL_SURFACE_MESH_TEST_SUITE)
  Vertex_index vertex_freelist() const
  {
    return Vertex_index(vertices_freelist_);
  }

  Face_index face_freelist() const
  {
    return Face_index(faces_freelist_);
  }

  Edge_index edge_freelist() const
  {
    return Edge_index(edges_freelist_>>1);
  }
#endif

private: //--------------------------------------------------- helper functions


    /// make sure that the incoming halfedge of vertex v is a border halfedge
    /// if `v` is a border vertex.
    void adjust_incoming_halfedge(Vertex_index v);

private: //------------------------------------------------------- private data
    Properties::Property_container<Self, Vertex_index> vprops_;
    Properties::Property_container<Self, Halfedge_index> hprops_;
    Properties::Property_container<Self, Edge_index> eprops_;
    Properties::Property_container<Self, Face_index> fprops_;

    Property_map<Vertex_index, Vertex_connectivity>      vconn_;
    Property_map<Halfedge_index, Halfedge_connectivity>  hconn_;
    Property_map<Face_index, Face_connectivity>          fconn_;

    Property_map<Vertex_index, bool>  vremoved_;
    Property_map<Edge_index, bool>    eremoved_;
    Property_map<Face_index, bool>    fremoved_;

    Property_map<Vertex_index, Point>   vpoint_;

    size_type removed_vertices_;
    size_type removed_edges_;
    size_type removed_faces_;

    size_type vertices_freelist_;
    size_type edges_freelist_;
    size_type faces_freelist_;
    bool garbage_;
    bool recycle_;

    size_type anonymous_property_;
};

  /*! \addtogroup PkgSurface_mesh
   *
   * @{
   */

  /// \relates Surface_mesh
  /// Inserts `other` into `sm`.
  /// Shifts the indices of vertices of `other` by `sm.number_of_vertices() + sm.number_of_removed_vertices()`
  /// and analogously for halfedges, edges, and faces.
  /// Copies entries of all property maps which have the same name in `sm` and `other`.
  /// that is, property maps which are only in `other` are ignored.
  /// Also copies elements which are marked as removed, and concatenates the freelists of `sm` and `other`.

  template <typename P>
  Surface_mesh<P>& operator+=(Surface_mesh<P>& sm, const Surface_mesh<P>& other)
  {
    sm.join(other);
    return sm;
  }

  /// \relates Surface_mesh
  ///
  /// This operator calls `write_OFF(std::ostream& os, const CGAL::Surface_mesh& sm)`.
   template <typename P>
  std::ostream& operator<<(std::ostream& os, const Surface_mesh<P>& sm)
  {
    IO::write_OFF(os, sm);
    return os;
  }

  /// \relates Surface_mesh
  /// Extracts the surface mesh from an input stream in OFF
  /// and appends it to the surface mesh `sm`.
  ///
  /// This operator calls `read_OFF(std::istream& is, CGAL::Surface_mesh& sm)`.
  template <typename P>
  std::istream& operator>>(std::istream& is, Surface_mesh<P>& sm)
  {
    IO::read_OFF(is, sm);
    return is;
  }

 /*! @} */

template <typename P>
Surface_mesh<P>::
Surface_mesh()
{
    // allocate standard properties
    // same list is used in operator=() and assign()
    vconn_    = add_property_map<Vertex_index, Vertex_connectivity>("v:connectivity").first;
    hconn_    = add_property_map<Halfedge_index, Halfedge_connectivity>("h:connectivity").first;
    fconn_    = add_property_map<Face_index, Face_connectivity>("f:connectivity").first;
    vpoint_   = add_property_map<Vertex_index, Point>("v:point").first;
    vremoved_ = add_property_map<Vertex_index, bool>("v:removed", false).first;
    eremoved_ = add_property_map<Edge_index, bool>("e:removed", false).first;
    fremoved_ = add_property_map<Face_index, bool>("f:removed", false).first;

    removed_vertices_ = removed_edges_ = removed_faces_ = 0;
    vertices_freelist_ = edges_freelist_ = faces_freelist_ = (std::numeric_limits<size_type>::max)();
    garbage_ = false;
    recycle_ = true;
    anonymous_property_ = 0;
}


//-----------------------------------------------------------------------------
template <typename P>
Surface_mesh<P>&
Surface_mesh<P>::
operator=(const Surface_mesh<P>& rhs)
{
    if (this != &rhs)
    {
        // deep copy of property containers
        vprops_ = rhs.vprops_;
        hprops_ = rhs.hprops_;
        eprops_ = rhs.eprops_;
        fprops_ = rhs.fprops_;

        // property handles contain pointers, have to be reassigned
        vconn_    = property_map<Vertex_index, Vertex_connectivity>("v:connectivity").first;
        hconn_    = property_map<Halfedge_index, Halfedge_connectivity>("h:connectivity").first;
        fconn_    = property_map<Face_index, Face_connectivity>("f:connectivity").first;
        vremoved_ = property_map<Vertex_index, bool>("v:removed").first;
        eremoved_ = property_map<Edge_index, bool>("e:removed").first;
        fremoved_ = property_map<Face_index, bool>("f:removed").first;
        vpoint_   = property_map<Vertex_index, P>("v:point").first;

        // how many elements are removed?
        removed_vertices_  = rhs.removed_vertices_;
        removed_edges_     = rhs.removed_edges_;
        removed_faces_     = rhs.removed_faces_;
        vertices_freelist_ = rhs.vertices_freelist_;
        edges_freelist_    = rhs.edges_freelist_;
        faces_freelist_    = rhs.faces_freelist_;
        garbage_           = rhs.garbage_;
        recycle_           = rhs.recycle_;
        anonymous_property_ = rhs.anonymous_property_;
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
        vconn_    = add_property_map<Vertex_index, Vertex_connectivity>("v:connectivity").first;
        hconn_    = add_property_map<Halfedge_index, Halfedge_connectivity>("h:connectivity").first;
        fconn_    = add_property_map<Face_index, Face_connectivity>("f:connectivity").first;
        vpoint_   = add_property_map<Vertex_index, P>("v:point").first;
        vremoved_ = add_property_map<Vertex_index, bool>("v:removed", false).first;
        eremoved_ = add_property_map<Edge_index, bool>("e:removed", false).first;
        fremoved_ = add_property_map<Face_index, bool>("f:removed", false).first;

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
        recycle_           = rhs.recycle_;
        anonymous_property_ = rhs.anonymous_property_;
    }

    return *this;
}

//-----------------------------------------------------------------------------
template <typename P>
void
Surface_mesh<P>::
clear()
{
  clear_without_removing_property_maps();
  remove_all_property_maps();
}

template <typename P>
void
Surface_mesh<P>::
clear_without_removing_property_maps()
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
  vertices_freelist_ = edges_freelist_ = faces_freelist_ = (std::numeric_limits<size_type>::max)();
  garbage_ = false;
  recycle_ = true;
  anonymous_property_ = 0;
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
    props = properties<Vertex_index>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;

    out << "halfedge properties:\n";
    props = properties<Halfedge_index>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;

    out << "edge properties:\n";
    props = properties<Edge_index>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;

    out << "face properties:\n";
    props = properties<Face_index>();
    for (unsigned int i=0; i<props.size(); ++i)
        out << "\t" << props[i] << std::endl;
}
/// @endcond

//-----------------------------------------------------------------------------
template <typename P>
typename Surface_mesh<P>::Halfedge_index
Surface_mesh<P>::
halfedge(Vertex_index source, Vertex_index target) const
{
    CGAL_assertion(has_valid_index(source) && has_valid_index(target));

    Halfedge_index h  = halfedge(target);
    const Halfedge_index hh = h;

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

    return Halfedge_index();
}


//-----------------------------------------------------------------------------
template <typename P>
void
Surface_mesh<P>::
adjust_incoming_halfedge(Vertex_index v)
{
    Halfedge_index h  = halfedge(v);
    Halfedge_index hh = h;

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
template <typename Range>
typename Surface_mesh<P>::Face_index
Surface_mesh<P>::add_face(const Range& r)
{
  return CGAL::Euler::add_face(r, *this);
}

  /// @endcond

//-----------------------------------------------------------------------------
template <typename P>
typename Surface_mesh<P>::size_type
Surface_mesh<P>::
degree(Vertex_index v) const
{
    Halfedge_index h = halfedge(v);

    if(h == null_halfedge()){
      return 0;
    }
    size_type count(0);
    Halfedge_index done = h;
    do {
      ++count;
      h = opposite(next(h));
    }while(h != done);

    return count;
}


//-----------------------------------------------------------------------------
template <typename P>
typename Surface_mesh<P>::size_type
Surface_mesh<P>::
degree(Face_index f) const
{
    size_type count(0);
    if(halfedge(f) == null_halfedge()){
      return 0;
    }
    Vertex_around_face_circulator fvit(halfedge(f),*this);
    Vertex_around_face_circulator fvend = fvit;
    if(fvit) do {
        ++count;
    } while (++fvit != fvend);

    return count;
}

template <typename P> template< typename Visitor>
void
Surface_mesh<P>::
collect_garbage(Visitor &visitor)
{
    if (!has_garbage())
    {
      return;
    }

    int  i, i0, i1,
    nV(num_vertices()),
    nE(num_edges()),
    nH(num_halfedges()),
    nF(num_faces());

    Vertex_index    v;
    Halfedge_index  h;
    Face_index      f;


    // setup index mapping%
    Property_map<Vertex_index, Vertex_index>      vmap = add_property_map<Vertex_index, Vertex_index>("v:garbage-collection").first;
    Property_map<Halfedge_index, Halfedge_index>  hmap = add_property_map<Halfedge_index, Halfedge_index>("h:garbage-collection").first;
    Property_map<Face_index, Face_index>          fmap = add_property_map<Face_index, Face_index>("f:garbage-collection").first;
    for (i=0; i<nV; ++i)
        vmap[Vertex_index(i)] = Vertex_index(i);
    for (i=0; i<nH; ++i)
        hmap[Halfedge_index(i)] = Halfedge_index(i);
    for (i=0; i<nF; ++i)
        fmap[Face_index(i)] = Face_index(i);



    // really remove vertices
    if (nV > 0)
    {
        i0=0;  i1=nV-1;

        while (1)
        {
            // find first removed and last un-removed
            while (!vremoved_[Vertex_index(i0)] && i0 < i1)  ++i0;
            while ( vremoved_[Vertex_index(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            vprops_.swap(i0, i1);
        };

        // remember new size
        nV = vremoved_[Vertex_index(i0)] ? i0 : i0+1;
    }

    // really remove edges
    if (nE > 0)
    {
        i0=0;  i1=nE-1;

        while (1)
        {
            // find first removed and last un-removed
            while (!eremoved_[Edge_index(i0)] && i0 < i1) ++i0;
            while ( eremoved_[Edge_index(i1)] && i0 < i1) --i1;
            if (i0 >= i1) break;

            // swap
            eprops_.swap(i0, i1);
            hprops_.swap(2*i0,   2*i1);
            hprops_.swap(2*i0+1, 2*i1+1);
        };

        // remember new size
        nE = eremoved_[Edge_index(i0)] ? i0 : i0+1;
        nH = 2*nE;
    }


    // really remove faces
    if (nF > 0)
    {
        i0=0;  i1=nF-1;

        while (1)
        {
            // find 1st removed and last un-removed
            while (!fremoved_[Face_index(i0)] && i0 < i1)  ++i0;
            while ( fremoved_[Face_index(i1)] && i0 < i1)  --i1;
            if (i0 >= i1) break;

            // swap
            fprops_.swap(i0, i1);
        };

        // remember new size
        nF = fremoved_[Face_index(i0)] ? i0 : i0+1;
    }


    // update vertex connectivity
    for (i=0; i<nV; ++i)
    {
        v = Vertex_index(i);
        if (!is_isolated(v))
            set_halfedge(v, hmap[halfedge(v)]);
    }


    // update halfedge connectivity
    for (i=0; i<nH; ++i)
    {
        h = Halfedge_index(i);
        set_target(h, vmap[target(h)]);
        set_next(h, hmap[next(h)]);
        if (!is_border(h))
            set_face(h, fmap[face(h)]);
    }


    // update indices of faces
    for (i=0; i<nF; ++i)
    {
        f = Face_index(i);
        set_halfedge(f, hmap[halfedge(f)]);
    }

    //apply visitor before invalidating the maps
    visitor(vmap, hmap, fmap);
    // remove index maps
    remove_property_map<Vertex_index>(vmap);
    remove_property_map<Halfedge_index>(hmap);
    remove_property_map<Face_index>(fmap);

    // finally resize arrays
    vprops_.resize(nV); vprops_.shrink_to_fit();
    hprops_.resize(nH); hprops_.shrink_to_fit();
    eprops_.resize(nE); eprops_.shrink_to_fit();
    fprops_.resize(nF); fprops_.shrink_to_fit();

    removed_vertices_ = removed_edges_ = removed_faces_ = 0;
    vertices_freelist_ = edges_freelist_ = faces_freelist_ = -1;
    garbage_ = false;
}

#ifndef DOXYGEN_RUNNING
namespace collect_garbage_internal {
struct Dummy_visitor{
  template<typename A, typename B, typename C>
  void operator()(const A&, const B&, const C&)
  {}
};

}
#endif

template <typename P>
void
Surface_mesh<P>::
collect_garbage()
{
  collect_garbage_internal::Dummy_visitor visitor;
  collect_garbage(visitor);
}


template <typename P>
void
Surface_mesh<P>::
set_recycle_garbage(bool b)
{
  recycle_ = b;
}


template <typename P>
bool
Surface_mesh<P>::
does_recycle_garbage() const
{
  return recycle_;
}


namespace internal{
  namespace handle {
    template <>
    struct Hash_functor<SM_Vertex_index>{
      std::size_t
      operator()(const SM_Vertex_index i)
      {
        return i;
      }
    };

    template <>
    struct Hash_functor<SM_Halfedge_index>{
      std::size_t
      operator()(const SM_Halfedge_index i)
      {
        return i;
      }
    };

    template <>
    struct Hash_functor<SM_Edge_index>{
      std::size_t
      operator()(const SM_Edge_index i)
      {
        return i;
      }
    };

    template <>
    struct Hash_functor<SM_Face_index>{
      std::size_t
      operator()(const SM_Face_index i)
      {
        return i;
      }
    };
  }
}

} // namespace CGAL

#ifndef DOXYGEN_RUNNING

namespace std {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash
#endif

#ifndef CGAL_CFG_NO_STD_HASH

  template <>
  struct hash<CGAL::SM_Halfedge_index >
    : public CGAL::cpp98::unary_function<CGAL::SM_Halfedge_index, std::size_t> {

    std::size_t operator()(const CGAL::SM_Halfedge_index& i) const
    {
      return i;
    }
  };

  template <>
  struct hash<CGAL::SM_Vertex_index >
    : public CGAL::cpp98::unary_function<CGAL::SM_Vertex_index, std::size_t>  {

    std::size_t operator()(const CGAL::SM_Vertex_index& i) const
    {
      return i;
    }
  };

  template <>
  struct hash<CGAL::SM_Face_index >
    : public CGAL::cpp98::unary_function<CGAL::SM_Face_index, std::size_t> {

    std::size_t operator()(const CGAL::SM_Face_index& i) const
    {
      return i;
    }
  };

  template <>
  struct hash<CGAL::SM_Edge_index >
    : public CGAL::cpp98::unary_function<CGAL::SM_Edge_index, std::size_t>  {

    std::size_t operator()(const CGAL::SM_Edge_index& i) const
    {
      return i;
    }
  };
#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace std

namespace boost {
  template <>
  struct hash<CGAL::SM_Vertex_index > {
    std::size_t operator()(const CGAL::SM_Vertex_index& i) const
    {
      return i;
    }
  };

} // namespace boost

#endif // DOXYGEN_RUNNING

#include <CGAL/enable_warnings.h>

#endif /* CGAL_SURFACE_MESH_H */
