//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
//
// SPDX-License-Identifier: GPL-2.0-only
//
//=============================================================================


#ifndef SURFACE_MESH_H
#define SURFACE_MESH_H


//== INCLUDES =================================================================


#include "Vector.h"
#include "types.h"
#include "properties.h"


//== CLASS DEFINITION =========================================================


/// \addtogroup surface_mesh surface_mesh
/// @{

/// A halfedge data structure for polygonal meshes.
class Surface_mesh
{

public: //------------------------------------------------------ topology types


    /// Base class for all topology types (internally it is basically an index)
    /// \sa Vertex, Halfedge, Edge, Face
    class Base_descriptor
    {
    public:

        /// constructor
        explicit Base_descriptor(int _idx=-1) : idx_(_idx) {}

        /// Get the underlying index of this handle
        int idx() const { return idx_; }

        /// reset handle to be invalid (index=-1)
        void reset() { idx_=-1; }

        /// return whether the handle is valid, i.e., the index is not equal to -1.
        bool is_valid() const { return idx_ != -1; }

        /// are two handles equal?
        bool operator==(const Base_descriptor& _rhs) const {
            return idx_ == _rhs.idx_;
        }

        /// are two handles different?
        bool operator!=(const Base_descriptor& _rhs) const {
            return idx_ != _rhs.idx_;
        }

        /// compare operator useful for sorting handles
        bool operator<(const Base_descriptor& _rhs) const {
            return idx_ < _rhs.idx_;
        }

    private:
        friend class Vertex_iterator;
        friend class Halfedge_iterator;
        friend class Edge_iterator;
        friend class Face_iterator;
        friend class Surface_mesh;
        int idx_;
    };


    /// this type represents a vertex (internally it is basically an index)
    ///  \sa Halfedge, Edge, Face
    struct Vertex : public Base_descriptor
    {
        /// default constructor (with invalid index)
        explicit Vertex(int _idx=-1) : Base_descriptor(_idx) {}
        std::ostream& operator<<(std::ostream& os) const { return os << 'v' << idx(); }
    };


    /// this type represents a halfedge (internally it is basically an index)
    /// \sa Vertex, Edge, Face
    struct Halfedge : public Base_descriptor
    {
        /// default constructor (with invalid index)
        explicit Halfedge(int _idx=-1) : Base_descriptor(_idx) {}
    };


    /// this type represents an edge (internally it is basically an index)
    /// \sa Vertex, Halfedge, Face
    struct Edge : public Base_descriptor
    {
        /// default constructor (with invalid index)
        explicit Edge(int _idx=-1) : Base_descriptor(_idx) {}
    };


    /// this type represents a face (internally it is basically an index)
    /// \sa Vertex, Halfedge, Edge
    struct Face : public Base_descriptor
    {
        /// default constructor (with invalid index)
        explicit Face(int _idx=-1) : Base_descriptor(_idx) {}
    };




public: //-------------------------------------------------- connectivity types

    /// This type stores the vertex connectivity
    /// \sa Halfedge_connectivity, Face_connectivity
    struct Vertex_connectivity
    {
        /// an outgoing halfedge per vertex (it will be a boundary halfedge for boundary vertices)
        Halfedge  halfedge_;
    };


    /// This type stores the halfedge connectivity
    /// \sa Vertex_connectivity, Face_connectivity
    struct Halfedge_connectivity
    {
        /// face incident to halfedge
        Face      face_;
        /// vertex the halfedge points to
        Vertex    vertex_;
        /// next halfedge within a face (or along a boundary)
        Halfedge  next_halfedge_;
        /// previous halfedge within a face (or along a boundary)
        Halfedge  prev_halfedge_;
    };


    /// This type stores the face connectivity
    /// \sa Vertex_connectivity, Halfedge_connectivity
    struct Face_connectivity
    {
        /// a halfedge that is part of the face
        Halfedge  halfedge_;
    };




public: //------------------------------------------------------ property types

    /// Vertex property of type T
    /// \sa Halfedge_property, Edge_property, Face_property
    template <class T> class Vertex_property : public Property<T>
    {
    public:

        /// default constructor
        explicit Vertex_property() {}
        explicit Vertex_property(Property<T> p) : Property<T>(p) {}

        /// access the data stored for vertex \c v
        typename Property<T>::reference operator[](Vertex v)
        {
            return Property<T>::operator[](v.idx());
        }

        /// access the data stored for vertex \c v
        typename Property<T>::const_reference operator[](Vertex v) const
        {
            return Property<T>::operator[](v.idx());
        }
    };


    /// Halfedge property of type T
    /// \sa Vertex_property, Edge_property, Face_property
    template <class T> class Halfedge_property : public Property<T>
    {
    public:

        /// default constructor
        explicit Halfedge_property() {}
        explicit Halfedge_property(Property<T> p) : Property<T>(p) {}

        /// access the data stored for halfedge \c h
        typename Property<T>::reference operator[](Halfedge h)
        {
            return Property<T>::operator[](h.idx());
        }

        /// access the data stored for halfedge \c h
        typename Property<T>::const_reference operator[](Halfedge h) const
        {
            return Property<T>::operator[](h.idx());
        }
    };


    /// Edge property of type T
    /// \sa Vertex_property, Halfedge_property, Face_property
    template <class T> class Edge_property : public Property<T>
    {
    public:

        /// default constructor
        explicit Edge_property() {}
        explicit Edge_property(Property<T> p) : Property<T>(p) {}

        /// access the data stored for edge \c e
        typename Property<T>::reference operator[](Edge e)
        {
            return Property<T>::operator[](e.idx());
        }

        /// access the data stored for edge \c e
        typename Property<T>::const_reference operator[](Edge e) const
        {
            return Property<T>::operator[](e.idx());
        }
    };


    /// Face property of type T
    /// \sa Vertex_property, Halfedge_property, Edge_property
    template <class T> class Face_property : public Property<T>
    {
    public:

        /// default constructor
        explicit Face_property() {}
        explicit Face_property(Property<T> p) : Property<T>(p) {}

        /// access the data stored for face \c f
        typename Property<T>::reference operator[](Face f)
        {
            return Property<T>::operator[](f.idx());
        }

        /// access the data stored for face \c f
        typename Property<T>::const_reference operator[](Face f) const
        {
            return Property<T>::operator[](f.idx());
        }
    };




public: //------------------------------------------------------ iterator types

    /// this class iterates linearly over all vertices
    /// \sa vertices_begin(), vertices_end()
    /// \sa Halfedge_iterator, Edge_iterator, Face_iterator
    class Vertex_iterator
    {
    public:

        /// Default constructor
        Vertex_iterator(Vertex v=Vertex(), const Surface_mesh* m=NULL) : hnd_(v), mesh_(m)
        {
            if (mesh_ && mesh_->garbage()) while (mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
        }

        /// get the vertex the iterator refers to
        Vertex operator*()  const { return  hnd_; }

        /// are two iterators equal?
        bool operator==(const Vertex_iterator& rhs) const
        {
            return (hnd_==rhs.hnd_);
        }

        /// are two iterators different?
        bool operator!=(const Vertex_iterator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment iterator
        Vertex_iterator& operator++()
        {
            ++hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
            return *this;
        }

        /// pre-decrement iterator
        Vertex_iterator& operator--()
        {
            --hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) --hnd_.idx_;
            return *this;
        }

    private:
        Vertex  hnd_;
        const Surface_mesh* mesh_;
    };


    /// this class iterates linearly over all halfedges
    /// \sa halfedges_begin(), halfedges_end()
    /// \sa Vertex_iterator, Edge_iterator, Face_iterator
    class Halfedge_iterator
    {
    public:

        /// Default constructor
        Halfedge_iterator(Halfedge h=Halfedge(), const Surface_mesh* m=NULL) : hnd_(h), mesh_(m)
        {
            if (mesh_ && mesh_->garbage()) while (mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
        }

        /// get the halfedge the iterator refers to
        Halfedge operator*()  const { return  hnd_; }

        /// are two iterators equal?
        bool operator==(const Halfedge_iterator& rhs) const
        {
            return (hnd_==rhs.hnd_);
        }

        /// are two iterators different?
        bool operator!=(const Halfedge_iterator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment iterator
        Halfedge_iterator& operator++()
        {
            ++hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
            return *this;
        }

        /// pre-decrement iterator
        Halfedge_iterator& operator--()
        {
            --hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) --hnd_.idx_;
            return *this;
        }

    private:
        Halfedge  hnd_;
        const Surface_mesh* mesh_;
    };


    /// this class iterates linearly over all edges
    /// \sa edges_begin(), edges_end()
    /// \sa Vertex_iterator, Halfedge_iterator, Face_iterator
    class Edge_iterator
    {
    public:

        /// Default constructor
        Edge_iterator(Edge e=Edge(), const Surface_mesh* m=NULL) : hnd_(e), mesh_(m)
        {
            if (mesh_ && mesh_->garbage()) while (mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
        }

        /// get the edge the iterator refers to
        Edge operator*()  const { return  hnd_; }

        /// are two iterators equal?
        bool operator==(const Edge_iterator& rhs) const
        {
            return (hnd_==rhs.hnd_);
        }

        /// are two iterators different?
        bool operator!=(const Edge_iterator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment iterator
        Edge_iterator& operator++()
        {
            ++hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
            return *this;
        }

        /// pre-decrement iterator
        Edge_iterator& operator--()
        {
            --hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) --hnd_.idx_;
            return *this;
        }

    private:
        Edge  hnd_;
        const Surface_mesh* mesh_;
    };


    /// this class iterates linearly over all faces
    /// \sa faces_begin(), faces_end()
    /// \sa Vertex_iterator, Halfedge_iterator, Edge_iterator
    class Face_iterator
    {
    public:

        /// Default constructor
        Face_iterator(Face f=Face(), const Surface_mesh* m=NULL) : hnd_(f), mesh_(m)
        {
            if (mesh_ && mesh_->garbage()) while (mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
        }

        /// get the face the iterator refers to
        Face operator*()  const { return  hnd_; }

        /// are two iterators equal?
        bool operator==(const Face_iterator& rhs) const
        {
            return (hnd_==rhs.hnd_);
        }

        /// are two iterators different?
        bool operator!=(const Face_iterator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment iterator
        Face_iterator& operator++()
        {
            ++hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) ++hnd_.idx_;
            return *this;
        }

        /// pre-decrement iterator
        Face_iterator& operator--()
        {
            --hnd_.idx_;
            assert(mesh_);
            while (mesh_->garbage() && mesh_->is_valid(hnd_) && mesh_->is_deleted(hnd_)) --hnd_.idx_;
            return *this;
        }

    private:
        Face  hnd_;
        const Surface_mesh* mesh_;
    };




public: //---------------------------------------------------- circulator types

    /// this class circulates through all one-ring neighbors of a vertex
    /// \sa Halfedge_around_vertex_circulator, Face_around_vertex_circulator
    class Vertex_around_vertex_circulator
    {
    public:

        /// default constructor
        Vertex_around_vertex_circulator(const Surface_mesh* m=NULL, Vertex v=Vertex())
        : mesh_(m)
        {
            if (mesh_) halfedge_ = mesh_->halfedge(v);
        }

        /// are two circulators equal?
        bool operator==(const Vertex_around_vertex_circulator& rhs) const
        {
            assert(mesh_);
            return ((mesh_==rhs.mesh_) && (halfedge_==rhs.halfedge_));
        }

        /// are two circulators different?
        bool operator!=(const Vertex_around_vertex_circulator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment (rotate couter-clockwise)
        Vertex_around_vertex_circulator& operator++()
        {
            assert(mesh_);
            halfedge_ = mesh_->ccw_rotated_halfedge(halfedge_);
            return *this;
        }

        /// pre-decrement (rotate clockwise)
        Vertex_around_vertex_circulator& operator--()
        {
            assert(mesh_);
            halfedge_ = mesh_->cw_rotated_halfedge(halfedge_);
            return *this;
        }

        /// get the vertex the circulator refers to
        Vertex operator*()  const
        {
            assert(mesh_);
            return mesh_->to_vertex(halfedge_);
        }

        /// cast to bool: true if vertex is not isolated
        operator bool() const { return halfedge_.is_valid(); }

    private:
        const Surface_mesh*  mesh_;
        Halfedge         halfedge_;
    };


    /// this class circulates through all outgoing halfedges of a vertex
    /// \sa Vertex_around_vertex_circulator, Face_around_vertex_circulator
    class Halfedge_around_vertex_circulator
    {
    public:

        /// default constructor
        Halfedge_around_vertex_circulator(const Surface_mesh* m=NULL, Vertex v=Vertex())
        : mesh_(m)
        {
            if (mesh_) halfedge_ = mesh_->halfedge(v);
        }

        /// are two circulators equal?
        bool operator==(const Halfedge_around_vertex_circulator& rhs) const
        {
            assert(mesh_);
            return ((mesh_==rhs.mesh_) && (halfedge_==rhs.halfedge_));
        }

        /// are two circulators different?
        bool operator!=(const Halfedge_around_vertex_circulator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment (rotate couter-clockwise)
        Halfedge_around_vertex_circulator& operator++()
        {
            assert(mesh_);
            halfedge_ = mesh_->ccw_rotated_halfedge(halfedge_);
            return *this;
        }

        /// pre-decrement (rotate clockwise)
        Halfedge_around_vertex_circulator& operator--()
        {
            assert(mesh_);
            halfedge_ = mesh_->cw_rotated_halfedge(halfedge_);
            return *this;
        }

        /// get the halfedge the circulator refers to
        Halfedge operator*() const { return halfedge_; }

        /// cast to bool: true if vertex is not isolated
        operator bool() const { return halfedge_.is_valid(); }

    private:
        const Surface_mesh*  mesh_;
        Halfedge         halfedge_;
    };


    /// this class circulates through all incident faces of a vertex
    /// \sa Vertex_around_vertex_circulator, Halfedge_around_vertex_circulator
    class Face_around_vertex_circulator
    {
    public:

        /// construct with mesh and vertex (vertex should not be isolated!)
        Face_around_vertex_circulator(const Surface_mesh* m=NULL, Vertex v=Vertex())
        : mesh_(m)
        {
            if (mesh_)
            {
                halfedge_ = mesh_->halfedge(v);
                if (halfedge_.is_valid() && mesh_->is_boundary(halfedge_))
                    operator++();
            }
        }

        /// are two circulators equal?
        bool operator==(const Face_around_vertex_circulator& rhs) const
        {
            assert(mesh_);
            return ((mesh_==rhs.mesh_) && (halfedge_==rhs.halfedge_));
        }

        /// are two circulators different?
        bool operator!=(const Face_around_vertex_circulator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment (rotates counter-clockwise)
        Face_around_vertex_circulator& operator++()
        {
            assert(mesh_ && halfedge_.is_valid());
            do {
                halfedge_ = mesh_->ccw_rotated_halfedge(halfedge_);
            } while (mesh_->is_boundary(halfedge_));
            return *this;
        }

        /// pre-decrement (rotate clockwise)
        Face_around_vertex_circulator& operator--()
        {
            assert(mesh_ && halfedge_.is_valid());
            do
                halfedge_ = mesh_->cw_rotated_halfedge(halfedge_);
            while (mesh_->is_boundary(halfedge_));
            return *this;
        }

        /// get the face the circulator refers to
        Face operator*() const
        {
            assert(mesh_ && halfedge_.is_valid());
            return mesh_->face(halfedge_);
        }

        /// cast to bool: true if vertex is not isolated
        operator bool() const { return halfedge_.is_valid(); }

    private:
        const Surface_mesh*  mesh_;
        Halfedge         halfedge_;
    };


    /// this class circulates through the vertices of a face
    /// \sa Halfedge_around_face_circulator
    class Vertex_around_face_circulator
    {
    public:

        /// default constructor
        Vertex_around_face_circulator(const Surface_mesh* m=NULL, Face f=Face())
        : mesh_(m)
        {
            if (mesh_) halfedge_ = mesh_->halfedge(f);
        }

        /// are two circulators equal?
        bool operator==(const Vertex_around_face_circulator& rhs) const
        {
            assert(mesh_);
            return ((mesh_==rhs.mesh_) && (halfedge_==rhs.halfedge_));
        }

        /// are two circulators different?
        bool operator!=(const Vertex_around_face_circulator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment (rotates counter-clockwise)
        Vertex_around_face_circulator& operator++()
        {
            assert(mesh_ && halfedge_.is_valid());
            halfedge_ = mesh_->next_halfedge(halfedge_);
            return *this;
        }

        /// pre-decrement (rotates clockwise)
        Vertex_around_face_circulator& operator--()
        {
            assert(mesh_ && halfedge_.is_valid());
            halfedge_ = mesh_->prev_halfedge(halfedge_);
            return *this;
        }

        /// get the vertex the circulator refers to
        Vertex operator*() const
        {
            assert(mesh_ && halfedge_.is_valid());
            return mesh_->to_vertex(halfedge_);
        }

    private:
        const Surface_mesh*  mesh_;
        Halfedge         halfedge_;
    };


    /// this class circulates through all halfedges of a face
    /// \sa Vertex_around_face_circulator
    class Halfedge_around_face_circulator
    {
    public:

        /// default constructor
        Halfedge_around_face_circulator(const Surface_mesh* m=NULL, Face f=Face())
        : mesh_(m)
        {
            if (mesh_) halfedge_ = mesh_->halfedge(f);
        }

        /// are two circulators equal?
        bool operator==(const Halfedge_around_face_circulator& rhs) const
        {
            assert(mesh_);
            return ((mesh_==rhs.mesh_) && (halfedge_==rhs.halfedge_));
        }

        /// are two circulators different?
        bool operator!=(const Halfedge_around_face_circulator& rhs) const
        {
            return !operator==(rhs);
        }

        /// pre-increment (rotates counter-clockwise)
        Halfedge_around_face_circulator& operator++()
        {
            assert(mesh_ && halfedge_.is_valid());
            halfedge_ = mesh_->next_halfedge(halfedge_);
            return *this;
        }

        /// pre-decrement (rotates clockwise)
        Halfedge_around_face_circulator& operator--()
        {
            assert(mesh_ && halfedge_.is_valid());
            halfedge_ = mesh_->prev_halfedge(halfedge_);
            return *this;
        }

        /// get the halfedge the circulator refers to
        Halfedge operator*() const { return halfedge_; }

    private:
        const Surface_mesh*  mesh_;
        Halfedge         halfedge_;
    };




public: //-------------------------------------------- constructor / destructor

    /// \name Construct, destruct, assignment
    //@{

    /// default constructor
    Surface_mesh();

    // destructor
    ~Surface_mesh();

    /// copy constructor: copies \c rhs to \c *this. performs a deep copy of all properties.
    Surface_mesh(const Surface_mesh& rhs) { operator=(rhs); }

    /// assign \c rhs to \c *this. performs a deep copy of all properties.
    Surface_mesh& operator=(const Surface_mesh& rhs);

    /// assign \c rhs to \c *this. does not copy custom properties.
    Surface_mesh& assign(const Surface_mesh& rhs);

    //@}




public: //------------------------------------------------------------- file IO

    /// \name File IO
    //@{

    /// read mesh from file \c filename. file extension determines file type.
    /// \sa write(const std::string& filename)
    bool read(const std::string& filename);

    /// write mesh to file \c filename. file extensions determines file type.
    /// \sa read(const std::string& filename)
    bool write(const std::string& filename) const;

    //@}




public: //----------------------------------------------- add new vertex / face

    /// \name Add new elements by hand
    //@{

    /// add a new vertex with position \c p
    Vertex add_vertex(const Point& p);

    /// add a new face with vertex list \c vertices
    /// \sa add_triangle, add_quad
    Face add_face(const std::vector<Vertex>& vertices);

    /// add a new triangle connecting vertices \c v1, \c v2, \c v3
    /// \sa add_face, add_quad
    Face add_triangle(Vertex v1, Vertex v2, Vertex v3);

    /// add a new quad connecting vertices \c v1, \c v2, \c v3, \c v4
    /// \sa add_triangle, add_face
    Face add_quad(Vertex v1, Vertex v2, Vertex v3, Vertex v4);

    //@}




public: //--------------------------------------------------- memory management

    /// \name Memory Management
    //@{

    /// returns number of (deleted and valid) vertices in the mesh
    unsigned int vertices_size() const { return (unsigned int) vprops_.size(); }
    /// returns number of (deleted and valid)halfedge in the mesh
    unsigned int halfedges_size() const { return (unsigned int) hprops_.size(); }
    /// returns number of (deleted and valid)edges in the mesh
    unsigned int edges_size() const { return (unsigned int) eprops_.size(); }
    /// returns number of (deleted and valid)faces in the mesh
    unsigned int faces_size() const { return (unsigned int) fprops_.size(); }


    /// returns number of vertices in the mesh
    unsigned int n_vertices() const { return vertices_size() - deleted_vertices_; }
    /// returns number of halfedge in the mesh
    unsigned int n_halfedges() const { return halfedges_size() - 2*deleted_edges_; }
    /// returns number of edges in the mesh
    unsigned int n_edges() const { return edges_size() - deleted_edges_; }
    /// returns number of faces in the mesh
    unsigned int n_faces() const { return faces_size() - deleted_faces_; }


    /// returns true iff the mesh is empty, i.e., has no vertices
    unsigned int empty() const { return n_vertices() == 0; }


    /// clear mesh: remove all vertices, edges, faces
    void clear();

    /// remove unused memory from vectors
    void free_memory();

    /// reserve memory (mainly used in file readers)
    void reserve(unsigned int nvertices,
                 unsigned int nedges,
                 unsigned int nfaces );


    /// remove deleted vertices/edges/faces
    void garbage_collection();


    /// returns whether vertex \c v is deleted
    /// \sa garbage_collection()
    bool is_deleted(Vertex v) const
    {
        return vdeleted_[v];
    }
    /// returns whether halfedge \c h is deleted
    /// \sa garbage_collection()
    bool is_deleted(Halfedge h) const
    {
        return edeleted_[edge(h)];
    }
    /// returns whether edge \c e is deleted
    /// \sa garbage_collection()
    bool is_deleted(Edge e) const
    {
        return edeleted_[e];
    }
    /// returns whether face \c f is deleted
    /// \sa garbage_collection()
    bool is_deleted(Face f) const
    {
        return fdeleted_[f];
    }


    /// return whether vertex \c v is valid, i.e. the index is stores it within the array bounds.
    bool is_valid(Vertex v) const
    {
        return (0 <= v.idx()) && (v.idx() < (int)vertices_size());
    }
    /// return whether halfedge \c h is valid, i.e. the index is stores it within the array bounds.
    bool is_valid(Halfedge h) const
    {
        return (0 <= h.idx()) && (h.idx() < (int)halfedges_size());
    }
    /// return whether edge \c e is valid, i.e. the index is stores it within the array bounds.
    bool is_valid(Edge e) const
    {
        return (0 <= e.idx()) && (e.idx() < (int)edges_size());
    }
    /// return whether face \c f is valid, i.e. the index is stores it within the array bounds.
    bool is_valid(Face f) const
    {
        return (0 <= f.idx()) && (f.idx() < (int)faces_size());
    }

    //@}




public: //---------------------------------------------- low-level connectivity

    /// \name Low-level connectivity
    //@{

    /// returns an outgoing halfedge of vertex \c v.
    /// if \c v is a boundary vertex this will be a boundary halfedge.
    Halfedge halfedge(Vertex v) const
    {
        return vconn_[v].halfedge_;
    }

    /// set the outgoing halfedge of vertex \c v to \c h
    void set_halfedge(Vertex v, Halfedge h)
    {
        vconn_[v].halfedge_ = h;
    }

    /// returns whether \c v is a boundary vertex
    bool is_boundary(Vertex v) const
    {
        Halfedge h(halfedge(v));
        return (!(h.is_valid() && face(h).is_valid()));
    }

    /// returns whether \c v is isolated, i.e., not incident to any face
    bool is_isolated(Vertex v) const
    {
        return !halfedge(v).is_valid();
    }


    /// returns the vertex the halfedge \c h points to
    Vertex to_vertex(Halfedge h) const
    {
        return hconn_[h].vertex_;
    }

    /// returns the vertex the halfedge \c h emanates from
    Vertex from_vertex(Halfedge h) const
    {
        return to_vertex(opposite_halfedge(h));
    }

    /// sets the vertex the halfedge \c h points to to \c v
    void set_vertex(Halfedge h, Vertex v)
    {
        hconn_[h].vertex_ = v;
    }

    /// returns the face incident to halfedge \c h
    Face face(Halfedge h) const
    {
        return hconn_[h].face_;
    }

    /// sets the incident face to halfedge \c h to \c f
    void set_face(Halfedge h, Face f)
    {
        hconn_[h].face_ = f;
    }

    /// returns the next halfedge within the incident face
    Halfedge next_halfedge(Halfedge h) const
    {
        return hconn_[h].next_halfedge_;
    }

    /// sets the next halfedge of \c h within the face to \c nh
    void set_next_halfedge(Halfedge h, Halfedge nh)
    {
        hconn_[h].next_halfedge_ = nh;
        hconn_[nh].prev_halfedge_ = h;
    }

    /// returns the previous halfedge within the incident face
    Halfedge prev_halfedge(Halfedge h) const
    {
        return hconn_[h].prev_halfedge_;
    }

    /// returns the opposite halfedge of \c h
    Halfedge opposite_halfedge(Halfedge h) const
    {
        return Halfedge((h.idx() & 1) ? h.idx()-1 : h.idx()+1);
    }

    /// returns the halfedge that is rotated counter-clockwise around the
    /// start vertex of \c h. it is the opposite halfedge of the previous halfedge of \c h.
    Halfedge ccw_rotated_halfedge(Halfedge h) const
    {
        return opposite_halfedge(prev_halfedge(h));
    }

    /// returns the halfedge that is rotated clockwise around the
    /// start vertex of \c h. it is the next halfedge of the opposite halfedge of \c h.
    Halfedge cw_rotated_halfedge(Halfedge h) const
    {
        return next_halfedge(opposite_halfedge(h));
    }

    /// return the edge that contains halfedge \c h as one of its two halfedges.
    Edge edge(Halfedge h) const
    {
        return Edge(h.idx() >> 1);
    }

    /// returns whether h is a boundary halfege, i.e., if its face does not exist.
    bool is_boundary(Halfedge h) const
    {
        return !face(h).is_valid();
    }


    /// returns the \c i'th halfedge of edge \c e. \c i has to be 0 or 1.
    Halfedge halfedge(Edge e, unsigned int i) const
    {
        assert(i<=1);
        return Halfedge((e.idx() << 1) + i);
    }

    /// returns the \c i'th vertex of edge \c e. \c i has to be 0 or 1.
    Vertex vertex(Edge e, unsigned int i) const
    {
        assert(i<=1);
        return to_vertex(halfedge(e, i));
    }

    /// returns whether \c e is a boundary edge, i.e., if one of its
    /// halfedges is a boundary halfedge.
    bool is_boundary(Edge e) const
    {
        return (is_boundary(halfedge(e, 0)) || is_boundary(halfedge(e, 1)));
    }

    /// returns a halfedge of face \c f
    Halfedge halfedge(Face f) const
    {
        return fconn_[f].halfedge_;
    }

    /// sets the halfedge of face \c f to \c h
    void set_halfedge(Face f, Halfedge h)
    {
        fconn_[f].halfedge_ = h;
    }

    /// returns whether \c f is a boundary face, i.e., it one of its edges is a boundary edge.
    bool is_boundary(Face f) const
    {
        Halfedge h  = halfedge(f);
        Halfedge hh = h;
        do
        {
            if (is_boundary(opposite_halfedge(h)))
                return true;
            h = next_halfedge(h);
        }
        while (h != hh);
        return false;
    }

    //@}




public: //--------------------------------------------------- property handling

    /// \name Property handling
    //@{

    /** add a vertex property of type \c T with name \c name and default value \c t.
     fails if a property named \c name exists already, since the name has to be unique.
     in this case it returns an invalid property */
    template <class T> Vertex_property<T> add_vertex_property(const std::string& name, const T t=T())
    {
        return Vertex_property<T>(vprops_.add<T>(name, t));
    }
    /** add a halfedge property of type \c T with name \c name and default value \c t.
     fails if a property named \c name exists already, since the name has to be unique.
     in this case it returns an invalid property */
    template <class T> Halfedge_property<T> add_halfedge_property(const std::string& name, const T t=T())
    {
        return Halfedge_property<T>(hprops_.add<T>(name, t));
    }
    /** add an edge property of type \c T with name \c name and default value \c t.
     fails if a property named \c name exists already, since the name has to be unique.
     in this case it returns an invalid property */
    template <class T> Edge_property<T> add_edge_property(const std::string& name, const T t=T())
    {
        return Edge_property<T>(eprops_.add<T>(name, t));
    }
    /** add a face property of type \c T with name \c name and default value \c t.
     fails if a property named \c name exists already, since the name has to be unique.
     in this case it returns an invalid property */
    template <class T> Face_property<T> add_face_property(const std::string& name, const T t=T())
    {
        return Face_property<T>(fprops_.add<T>(name, t));
    }


    /** get the vertex property named \c name of type \c T. returns an invalid
     Vertex_property if the property does not exist or if the type does not match. */
    template <class T> Vertex_property<T> get_vertex_property(const std::string& name) const
    {
        return Vertex_property<T>(vprops_.get<T>(name));
    }
    /** get the halfedge property named \c name of type \c T. returns an invalid
     Vertex_property if the property does not exist or if the type does not match. */
    template <class T> Halfedge_property<T> get_halfedge_property(const std::string& name) const
    {
        return Halfedge_property<T>(hprops_.get<T>(name));
    }
    /** get the edge property named \c name of type \c T. returns an invalid
     Vertex_property if the property does not exist or if the type does not match. */
    template <class T> Edge_property<T> get_edge_property(const std::string& name) const
    {
        return Edge_property<T>(eprops_.get<T>(name));
    }
    /** get the face property named \c name of type \c T. returns an invalid
     Vertex_property if the property does not exist or if the type does not match. */
    template <class T> Face_property<T> get_face_property(const std::string& name) const
    {
        return Face_property<T>(fprops_.get<T>(name));
    }


    /** if a vertex property of type \c T with name \c name exists, it is returned.
     otherwise this property is added (with default value \c t) */
    template <class T> Vertex_property<T> vertex_property(const std::string& name, const T t=T())
    {
        return Vertex_property<T>(vprops_.get_or_add<T>(name, t));
    }
    /** if a halfedge property of type \c T with name \c name exists, it is returned.
     otherwise this property is added (with default value \c t) */
    template <class T> Halfedge_property<T> halfedge_property(const std::string& name, const T t=T())
    {
        return Halfedge_property<T>(hprops_.get_or_add<T>(name, t));
    }
    /** if an edge property of type \c T with name \c name exists, it is returned.
     otherwise this property is added (with default value \c t) */
    template <class T> Edge_property<T> edge_property(const std::string& name, const T t=T())
    {
        return Edge_property<T>(eprops_.get_or_add<T>(name, t));
    }
    /** if a face property of type \c T with name \c name exists, it is returned.
     otherwise this property is added (with default value \c t) */
    template <class T> Face_property<T> face_property(const std::string& name, const T t=T())
    {
        return Face_property<T>(fprops_.get_or_add<T>(name, t));
    }


    /// remove the vertex property \c p
    template <class T> void remove_vertex_property(Vertex_property<T>& p)
    {
        vprops_.remove(p);
    }
    /// remove the halfedge property \c p
    template <class T> void remove_halfedge_property(Halfedge_property<T>& p)
    {
        hprops_.remove(p);
    }
    /// remove the edge property \c p
    template <class T> void remove_edge_property(Edge_property<T>& p)
    {
        eprops_.remove(p);
    }
    /// remove the face property \c p
    template <class T> void remove_face_property(Face_property<T>& p)
    {
        fprops_.remove(p);
    }


    /** get the type_info \c T of vertex property named \c name. returns an typeid(void)
     if the property does not exist or if the type does not match. */
    const std::type_info& get_vertex_property_type(const std::string& name)
    {
        return vprops_.get_type(name);
    }
    /** get the type_info \c T of halfedge property named \c name. returns an typeid(void)
     if the property does not exist or if the type does not match. */
    const std::type_info& get_halfedge_property_type(const std::string& name)
    {
        return hprops_.get_type(name);
    }
    /** get the type_info \c T of edge property named \c name. returns an typeid(void)
     if the property does not exist or if the type does not match. */
    const std::type_info& get_edge_property_type(const std::string& name)
    {
        return eprops_.get_type(name);
    }
    /** get the type_info \c T of face property named \c name. returns an typeid(void)
     if the property does not exist or if the type does not match. */
    const std::type_info& get_face_property_type(const std::string& name)
    {
        return fprops_.get_type(name);
    }


    /// returns the names of all vertex properties
    std::vector<std::string> vertex_properties() const
    {
        return vprops_.properties();
    }
    /// returns the names of all halfedge properties
    std::vector<std::string> halfedge_properties() const
    {
        return hprops_.properties();
    }
    /// returns the names of all edge properties
    std::vector<std::string> edge_properties() const
    {
        return eprops_.properties();
    }
    /// returns the names of all face properties
    std::vector<std::string> face_properties() const
    {
        return fprops_.properties();
    }
    /// prints the names of all properties
    void property_stats() const;

    //@}




public: //--------------------------------------------- iterators & circulators

    /// \name Iterators & Circulators
    //@{

    /// returns start iterator for vertices
    Vertex_iterator vertices_begin() const
    {
        return Vertex_iterator(Vertex(0), this);
    }

    /// returns end iterator for vertices
    Vertex_iterator vertices_end() const
    {
        return Vertex_iterator(Vertex(vertices_size()), this);
    }

    /// returns start iterator for halfedges
    Halfedge_iterator halfedges_begin() const
    {
        return Halfedge_iterator(Halfedge(0), this);
    }

    /// returns end iterator for halfedges
    Halfedge_iterator halfedges_end() const
    {
        return Halfedge_iterator(Halfedge(halfedges_size()), this);
    }

    /// returns start iterator for edges
    Edge_iterator edges_begin() const
    {
        return Edge_iterator(Edge(0), this);
    }

    /// returns end iterator for edges
    Edge_iterator edges_end() const
    {
        return Edge_iterator(Edge(edges_size()), this);
    }

    /// returns start iterator for faces
    Face_iterator faces_begin() const
    {
        return Face_iterator(Face(0), this);
    }

    /// returns end iterator for faces
    Face_iterator faces_end() const
    {
        return Face_iterator(Face(faces_size()), this);
    }

    /// returns circulator for vertices around vertex \c v
    Vertex_around_vertex_circulator vertices(Vertex v) const
    {
        return Vertex_around_vertex_circulator(this, v);
    }

    /// returns circulator for outgoing halfedges around vertex \c v
    Halfedge_around_vertex_circulator halfedges(Vertex v) const
    {
        return Halfedge_around_vertex_circulator(this, v);
    }

    /// returns circulator for faces around vertex \c v
    Face_around_vertex_circulator faces(Vertex v) const
    {
        return Face_around_vertex_circulator(this, v);
    }

    /// returns circulator for vertices of face \c f
    Vertex_around_face_circulator vertices(Face f) const
    {
        return Vertex_around_face_circulator(this, f);
    }

    /// returns circulator for halfedges of face \c f
    Halfedge_around_face_circulator halfedges(Face f) const
    {
        return Halfedge_around_face_circulator(this, f);
    }

    //@}





public: //--------------------------------------------- higher-level operations

    /// \name Higher-level Topological Operations
    //@{

    /// returns whether the mesh a triangle mesh. this function simply tests
    /// each face, and therefore is not very efficient.
    /// \sa trianglate(), triangulate(Face)
    bool is_triangle_mesh() const;

    /// triangulate the entire mesh, by calling triangulate(Face) for each face.
    /// \sa trianglate(Face)
    void triangulate();

    /// triangulate the face \c f
    /// \sa trianglate()
    void triangulate(Face f);


    /// returns whether collapsing the halfedge \c h is topologically legal.
    /// \attention This function is only valid for triangle meshes.
    bool is_collapse_ok(Halfedge h);

    /** Collapse the halfedge \c h by moving its start vertex into its target
     vertex. For non-boundary halfedges this function removes one vertex, three
     edges, and two faces. For boundary halfedges it removes one vertex, two
     edges and one face.
     \attention This function is only valid for triangle meshes.
     \attention Halfedge collapses might lead to invalid faces. Call
     is_collapse_ok(Halfedge) to be sure the collapse is legal.
     \attention The removed items are only marked as deleted. You have
     to call garbage_collection() to finally remove them.
     */
    void collapse(Halfedge h);


    /** Split the face \c f by first adding point \c p to the mesh and then
     inserting edges between \c p and the vertices of \c f. For a triangle
     this is a standard one-to-three split.
     \sa split(Face, Vertex)
     */
    Vertex split(Face f, const Point& p) { Vertex v=add_vertex(p); split(f,v); return v; }

    /** Split the face \c f by inserting edges between \c p and the vertices
     of \c f. For a triangle this is a standard one-to-three split.
     \sa split(Face, const Point&)
     */
    void split(Face f, Vertex v);


    /** Split the edge \c e by first adding point \c p to the mesh and then
     connecting it to the two vertices of the adjacent triangles that are
     opposite to edge \c e.
     \attention This function is only valid for triangle meshes.
     \sa split(Edge, Vertex)
     */
    Vertex split(Edge e, const Point& p) { Vertex v=add_vertex(p); split(e,v); return v; }

    /** Split the edge \c e by connecting vertex \c v it to the two vertices
     of the adjacent triangles that are opposite to edge \c e.
     \attention This function is only valid for triangle meshes.
     \sa split(Edge, Point)
     */
    void split(Edge e, Vertex v);


    /** Subdivide the edge \c e = (v0,v1) by splitting it into the two edge
     (v0,p) and (p,v1). Note that this function does not introduce any other
     edge or faces. It simply splits the edge.
     \sa insert_vertex(Edge, Vertex)
     */
    Vertex insert_vertex(Edge e, const Point& p) { Vertex v=add_vertex(p); insert_vertex(e,v); return v; }

    /** Subdivide the edge \c e = (v0,v1) by splitting it into the two edge
     (v0,v) and (v,v1). Note that this function does not introduce any other
     edge or faces. It simply splits the edge.
     \sa insert_vertex(Edge, Point)
     */
    void insert_vertex(Edge e, Vertex v);


    /// insert edge between the to-vertices of h0 and h1.
    /// \attention h0 and h1 have to belong to the same face
    void insert_edge(Halfedge h0, Halfedge h1);


    /** Check whether flipping edge \c e is topologically
     \attention This function is only valid for triangle meshes.
     \sa flip(Edge)
     */
    bool is_flip_ok(Edge e) const;

    /** Flip edge \c e: Remove edge \c e and add an edge between the two vertices
     opposite to edge \c e of the two incident triangles.
     \attention This function is only valid for triangle meshes.
     \sa is_flip_ok(Edge)
     */
    void flip(Edge e);


    /** returns the valence (number of incident edges or neighboring vertices)
     of vertex \c v. */
    unsigned int valence(Vertex v) const;

    /// returns the valence of face \c f (its number of vertices)
    unsigned int valence(Face f) const;

    //@}




public: //------------------------------------------ geometry-related functions

    /// \name Geometry-related Functions
    //@{

    /// position of a vertex (read only)
    const Point& position(Vertex v) const { return vpoint_[v]; }

    /// position of a vertex
    Point& position(Vertex v) { return vpoint_[v]; }

    /// compute face normals by calling compute_face_normal(Face) for each face.
    void update_face_normals();

    /// compute normal vector of face \c f.
    Normal compute_face_normal(Face f) const;

    /// compute vertex normals by calling compute_vertex_normal(Vertex) for each vertex.
    void update_vertex_normals();

    /// compute normal vector of vertex \c v.
    Normal compute_vertex_normal(Vertex v) const;

    /// compute the length of edge \c e.
    Scalar edge_length(Edge e) const;

    //@}




private: //---------------------------------------------- allocate new elements

    /// allocate a new vertex, resize vertex properties accordingly.
    Vertex new_vertex()
    {
        vprops_.push_back();
        return Vertex(vertices_size()-1);
    }

    /// allocate a new edge, resize edge and halfedge properties accordingly.
    Halfedge new_edge(Vertex start, Vertex end)
    {
        assert(start != end);

        eprops_.push_back();
        hprops_.push_back();
        hprops_.push_back();

        Halfedge h0(halfedges_size()-2);
        Halfedge h1(halfedges_size()-1);

        set_vertex(h0, end);
        set_vertex(h1, start);

        return h0;
    }

    /// allocate a new face, resize face properties accordingly.
    Face new_face()
    {
        fprops_.push_back();
        return Face(faces_size()-1);
    }




private: //--------------------------------------------------- helper functions

    /// find the halfedge from start to end
    Halfedge find_halfedge(Vertex start, Vertex end) const;

    /** make sure that the outgoing halfedge of vertex v is a boundary halfedge
     if v is a boundary vertex. */
    void adjust_outgoing_halfedge(Vertex v);

    /// Helper for halfedge collapse
    void remove_edge(Halfedge h);

    /// Helper for halfedge collapse
    void remove_loop(Halfedge h);

    /// are there deleted vertices, edges or faces?
    bool garbage() const { return garbage_; }



private: //------------------------------------------------------- private data

    friend bool read_poly(Surface_mesh& mesh, const std::string& filename);

    Property_container vprops_;
    Property_container hprops_;
    Property_container eprops_;
    Property_container fprops_;

    Vertex_property<Vertex_connectivity>      vconn_;
    Halfedge_property<Halfedge_connectivity>  hconn_;
    Face_property<Face_connectivity>          fconn_;

    Vertex_property<bool>  vdeleted_;
    Edge_property<bool>    edeleted_;
    Face_property<bool>    fdeleted_;

    Vertex_property<Point>   vpoint_;
    Vertex_property<Normal>  vnormal_;
    Face_property<Normal>    fnormal_;

    unsigned int deleted_vertices_;
    unsigned int deleted_edges_;
    unsigned int deleted_faces_;
    bool garbage_;
};


//------------------------------------------------------------ output operators


inline std::ostream& operator<<(std::ostream& os, Surface_mesh::Vertex v)
{
    return (os << 'v' << v.idx());
}

inline std::ostream& operator<<(std::ostream& os, Surface_mesh::Halfedge h)
{
    return (os << 'h' << h.idx());
}

inline std::ostream& operator<<(std::ostream& os, Surface_mesh::Edge e)
{
    return (os << 'e' << e.idx());
}

inline std::ostream& operator<<(std::ostream& os, Surface_mesh::Face f)
{
    return (os << 'f' << f.idx());
}


//=============================================================================
/// @}
//=============================================================================
#endif // SURFACE_MESH_H
//=============================================================================
