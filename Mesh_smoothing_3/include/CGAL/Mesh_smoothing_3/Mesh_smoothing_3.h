#pragma once
#include <cassert>
#include <unordered_map>
#include <functional>
#include <Eigen/Eigen>
#include "mesh_representations.h"
#include "CGAL/Mesh_smoothing_3/internal/Tetrahedral_mesh_smoother.h"


namespace Mesh_smoothing_3 {

namespace Parameters {
    /*!
        Recommended parameters for boundary weights: 
        DEFAULT: 1, for standard mesh smoothing
        STRONG: 10, when surface deviation should be avoided
        SOFT: 1e-3, when only slight influence of the boundary is required
    */
    enum BOUNDARY_WEIGHTING_MODE {DEFAULT, STRONG, SOFT};

    enum OPTIMIZATION_MODE { CONFORMAL, MIXED_WITH_SIZE, LAPLACIAN, MIN_MAXING_CONFORMAL };

    enum PREDICATES_MODE { NO_CHECK, STATUS_UPDATE, CHECK_IN_LBFGS, STRONG_ENFORCEMENT }; 

}

/*!
* \ingroup pkgMeshSmoothing3Classes
* 
* \brief smooth a tetrahedral mesh with optional constraints on the boundary and along a curve network
*
* The class `Mesh_smoother` provides an interface to optimize a tetrahedral mesh with an energy, 
* while allowing the user to specify soft-constraints on the boundary and along a curve network. 
* This class is designed for linking with various mesh representations.
* For usage with CGAL common mesh representations, refer to `C3t3_smoother`.
* 
* \tparam TetrahedralMesh encodes the volumetric mesh and must be a model of `MeshDataStructure`.
* 
* \tparam BoundaryMesh represents a set of selected polygons with points from Tetrahedral mesh 
* that will be used to define surface constraints on `TetrahedralMesh`. It must be a model of `SurfaceDataStructure`.
* 
* \tparam EdgeNetwork represents a set of selected segments with points from Tetrahedral mesh 
* that will be used to add curve constraints on `TetrahedralMesh`. It must be a model of `PolylinesDataStructure`.
* 
* 
* \sa `CGAL::Mesh_smoothing_3::C3t3_smoother`
* 
*/
template<
    typename TetrahedralMesh = default_structures::Empty_mesh,
    typename BoundaryMesh = default_structures::Empty_boundary<typename TetrahedralMesh::Vertex_descriptor>,
    typename EdgeNetwork = default_structures::Empty_edge_network<typename TetrahedralMesh::Vertex_descriptor>
>
class Mesh_smoother {
public:
    /*!

    */
    using Vertex_descriptor = typename TetrahedralMesh::Vertex_descriptor;
    
    /*!

    */
    using Cell_descriptor = typename TetrahedralMesh::Cell_descriptor;
    
    /*!

    */
    using Point_3 = typename TetrahedralMesh::Point_3;
    
    /*!

    */
    using Normal_3 = typename BoundaryMesh::Normal_3;
    
    /*!

    */
    using Surface_patch_index = typename BoundaryMesh::Surface_patch_index;

        /*!

    */
    using Curve_index = typename EdgeNetwork::Curve_index;

    /*!

    */
    template<typename T> using Vertex_descriptor_map = std::unordered_map<Vertex_descriptor, T>; // todo: manage to template that?
    template<typename T> using Cell_descriptor_map = std::unordered_map<Cell_descriptor, T>; // internal usage

    using Tetrahedral_mesh_smoother = Mesh_smoothing_3_internal::Tetrahedral_mesh_smoother<Surface_patch_index, Curve_index>;

    /*!
        Constructor of the class, it does not perform any operation
    */
    Mesh_smoother(TetrahedralMesh &mesh, BoundaryMesh const &boundary = BoundaryMesh(), EdgeNetwork const &edge_network = EdgeNetwork());

    /*!
        Locks (or unlock) all vertices contained in the input BoundaryMesh. 
        It is independent from the other locking functions, and unlocking will not interfere with manually set locks. 
        Set by default as true if no `set_boundary_query` is called. 
    */
    void set_locked_boundary(bool locked = true);

    /*!

    */
    void set_verbose(bool verbose = true);
    
    /*!
        Max number of high level iterations done by the smoother. Each iteration can be seen as equivalent to a Newton step. 
        Untangling may require a high number of iterations to converge (100), simple mesh improvement will usually require less (from 1 to 10). 
        Default is 1000. 
    */
    void set_max_number_of_iteration(unsigned number_of_iterations);

    // todo: locks are relatively inefficient because they use a map. Should it be improved?


    /*!
        Hard-constraining a certain vertex to its current coordinates
    */
    void set_vertex_lock(Vertex_descriptor vertex, bool locked = true);

    /*!
        Hard-constraining only 1 dimension for the given vertex
    */
    void set_vertex_dim_lock(Vertex_descriptor vertex, unsigned dimension, bool locked = true);

    /*!
        Setting all of the locked vertices at once. 
        Use this function to maximize performance when having tricky locks configuration. 
        Note that this will override all previously add locks, except those from `set_locked_boundary`. 
    */
    void set_vertices_dim_locks(Vertex_descriptor_map<std::array<bool, 3>> const &vertex_dimension_locks);

    /*!
        Helper function to set locked vertices from a container. 
    */
    template <typename Container>
    void set_locked_vertices(Container const &vertices) {
        for (Vertex_descriptor vertex : vertices) {
            set_vertex_lock(vertex, true);
        }
    }

    /*!
        Remove all previously set locks.
    */
    void clear_locks();

    // QUERIES FOR BOUNDARY PROJECTION    
    
    /*!
        Plane type to which the boundary polygons should align.
        the last parameter is a weight that is used to indicate the importance of the constraint.
        By default, use the value 1. 
    */
    using Plane = std::tuple<Point_3, Normal_3, double>;

    /*!
        Projection function. 
        
        \param coord center of polygon 
        \param surface_id patch id of the polygon as defined in `BoundaryMesh`
        \param radius indicator of the average edge length of the polygon
        \return Plane to which the polygon should align (usually closest tangential plane of the target geometry)

        \warning Must be thread safe -- it will be called in an OpenMP context. 
    */
    using Boundary_point_query = std::function<Plane (Point_3 const &coord, Surface_patch_index surface_id, double radius)>;

    /*!
        Projection function. 
        
        \param polygon vector containing the current points of the polygon 
        \param surface_id patch id of the polygon as defined in `BoundaryMesh`
        \return Plane to which the polygon should align (usually closest tangential plane of the target geometry)

        \warning Must be thread safe -- it will be called in an OpenMP context. 
    */
    using Boundary_polygon_query = std::function<Plane (std::vector<Point_3> const &polygon, Surface_patch_index surface_id)>;

    /*!
        Batch projection function. Equivalent of `Boundary_point_query`. 
        Will be called in more systemic manner, so will produce slower results. Use if single calls are slower, or to avoid thread-safe requirements. 
        
        \param coords vector with center of each polygon 
        \param surface_ids vector with patch ids of each polygon as defined in `BoundaryMesh`
        \param radii vector with indicator of the average edge length for each polygon
        \param[out] results vector with target plane for each polygon
    */
    using Boundary_point_batch_query = std::function<void (std::vector<Point_3> const &coords, std::vector<Surface_patch_index> &surface_ids, std::vector<double> &radii, std::vector<Plane> &results)>;

    /*!
        Batch projection function. Equivalent of `Boundary_polygon_query`. 
        Will be called in more systemic manner, so will produce slower results. Use if single calls are slower, or to avoid thread-safe requirements. 
        
        \param polygons vector with a vector of points for each polygon 
        \param surface_ids vector with patch ids of each polygon as defined in `BoundaryMesh`
        \param[out] results vector with target plane for each polygon
    */
    using Boundary_polygon_batch_query = std::function<void (std::vector<std::vector<Point_3>> const &polygons, std::vector<Surface_patch_index> &surface_id, std::vector<Plane> &results)>;


    /*!
        set constraint query for surface as a `Boundary_point_query`
    */
    void set_boundary_query(Boundary_point_query boundary_query);

    /*!
        set constraint query for surface as a `Boundary_polygon_query`
    */
    void set_boundary_query(Boundary_polygon_query boundary_query);
    
    /*!
        set constraint query for surface as a `Boundary_point_batch_query`
    */
    void set_boundary_query(Boundary_point_batch_query boundary_query);
    
    /*!
        set constraint query for surface as a `Boundary_polygon_batch_query`
    */
    void set_boundary_query(Boundary_polygon_batch_query boundary_query);


    // QUERIES FOR CURVE NETWORK PROJECTION
    /*!
        Tangent vector and point to which the curve segment should align.
        Last parameter is a weight that is used to indicate the importance of the constraint.
        By default, use the value 1.
    */
    using Curve_tangent = std::tuple<Point_3, Normal_3, double>;

    /*!
        Projection function. 
        
        \note Thread safe is not required here -- subject to change if specific cases require it.
        
        \param coord center of edge 
        \param curve_id id of the curve as defined in `EdgeNetwork`
        \param radius edge length
        \return Curve_tangent to which the edge should align (usually closest tangential direction of the target geometry)
    */
    using Curve_point_query = std::function<Curve_tangent (Point_3 const &coord, Curve_index curve_id, double radius)>;
    
    /*!
        Projection function. 

        \note Thread safe is not required here -- subject to change if specific cases require it.
        
        \param edge array containing the two points of the edge
        \param curve_id id of the curve as defined in `EdgeNetwork`
        \return Curve_tangent to which the edge should align (usually closest tangential direction of the target geometry)
    */
    using Curve_segment_query = std::function<Curve_tangent (std::array<Point_3, 2> const &edge, Curve_index curve_id)>;
    
    /*!
        Batch projection function. Equivalent of `Curve_point_query`. 
        Will be called in more systemic manner, so will produce slower results. Use if single calls are slower. 
        
        \param coords vector with center of each edge 
        \param curve_ids vector with curve id of each edge as defined in `EdgeNetwork`
        \param radii vector with edge length of each edge
        \param[out] results vector with target curve tangent for each edge
    */
    using Curve_point_batch_query = std::function<void (std::vector<Point_3> const &coord, std::vector<Curve_index> &curve_ids, std::vector<double> &radius, std::vector<Curve_tangent> &results)>;
    
    /*!
        Batch projection function. Equivalent of `Curve_segment_query`. 
        Will be called in more systemic manner, so will produce slower results. Use if single calls are slower. 
        
        \param edges vector with couple of points for each edge 
        \param curve_ids vector with curve ids of each edge as defined in `EdgeNetwork`
        \param[out] results vector with target curve tangent for each edge
    */
    using Curve_segment_batch_query = std::function<void (std::vector<std::array<Point_3, 2>> const &edges, std::vector<Curve_index> &curve_ids, std::vector<Curve_tangent> &results)>;


    /*!
        set constraint query for curves as a `Curve_point_query`
    */
    void set_curves_query(Curve_point_query curve_query);

    /*!
        set constraint query for curves as a `Curve_segment_query`
    */
    void set_curves_query(Curve_segment_query curve_query);
    
    /*!
        set constraint query for curves as a `Curve_point_batch_query`
    */
    void set_curves_query(Curve_point_batch_query curve_query);
    
    /*!
        set constraint query for curves as a `Curve_segment_batch_query`
    */
    void set_curves_query(Curve_segment_batch_query curve_query);

    // ADDING QUADRATIC TARGETS FOR VERTICES

    /*!
        set soft constraint for a given vertex
    */
    void set_vertex_target_position(Vertex_descriptor v, Point_3 const &target_positions, double weight = 1.0);

    /*!
        set the list of soft constraints for vertices. Note that it overrides previously set vertex constraints. 
    */
    void set_vertex_target_positions(std::vector<std::tuple<Vertex_descriptor, Point_3, double>> const &target_positions);

    /*!
        set soft constraints from a container of [vertex, coord] pairs.  
    */
    template <typename Container>
    void set_vertex_target_positions(Container const &target_positions) {
        for (auto const & [v, pt] : target_positions) {
           set_vertex_target_position(v, pt);
        }
    }

    /*!
        Clear previously set vertex soft constraints
    */
    void clear_vertex_target_positions();

    /*!
        Set the weight on the boundary term of the energy.
        Default is 1. 10 will strongly enforce the boundary matching, 1e-3 will only have a slight effect, favoring the inside meshing.
        Warning: large value can lead to convergence issues if the boundary is not initialized on its contraint. 
    */
    void set_boundary_weight(double weight); // large values can lead to convergence issues


    
    /*!
        Set boundary with pre-set recommended parameters.
    */
    void set_boundary_weight(Parameters::BOUNDARY_WEIGHTING_MODE mode);

    /*!
        Start the optimization procedure
    */
    bool run();

public: // for advanced usage. Do not touch if you do not know what you are doing.


    void set_optimization_mode(Parameters::OPTIMIZATION_MODE mode);

    void set_minimum_valid_edge_size(double edge_size);  // should be a minimum bound on the valid edge size of the mesh, used as a reference for untangling

    unsigned get_total_number_of_lbfgs_iterations() const;

    // return a negative value if you don't want enforce a sizing on a given vertex. 
    using Point_target_sizing_query = std::function<double (Vertex_descriptor vertex, Point_3 const &coord)>;
    using Point_target_sizing_batch_query = std::function<void (std::vector<Vertex_descriptor> const &vertices, std::vector<Point_3> const &coords, std::vector<double> &sizings)>;

    void set_vertex_target_sizing_query(Point_target_sizing_query query); // todo: implementation of that is not finished
    void set_vertex_target_sizing_query(Point_target_sizing_batch_query query);

    

    /*
        ONLY USE IF YOUR MESH NEEDED PREDICATES TO BE GENERATED VALID IN THE FIRST PLACE
        i.e. the code will not break already "good" meshes, 
        and enforcing strong orientation may significantly throttle minimization
        Parameter is set to STATUS_UPDATE if invalid elements are in the input
        Default is NO_CHECK
    */
    void set_predicates_mode(Parameters::PREDICATES_MODE mode);

    // first unsigned is number of det <= 0, second is with exact predicates
    // optional std::vector gives the det and exact predicate computation of each cell.
    // Warning: will have to create a full instance of smoother, so it is not a cheap function to call
    std::pair<unsigned, unsigned> get_nb_of_inverted_cells(Cell_descriptor_map<std::pair<double, bool>> *cell_determinants = nullptr);

    // only when CHECK_IN_LBFGS is enabled. 
    unsigned get_number_of_invalid_steps_measured_by_predicates() const;

    // user callback to enforce "infinite" energy for given coordinates changes
    // It can be used to have particular guarantees to preserve
    // but I do not recommend to rely on it, it is likely to block the optimization all together
    using Update_validator = std::function<bool (Vertex_descriptor_map<Point_3> const &coords)>;

    void set_update_validator(Update_validator func);

public: // for advanced monitoring

    using Iteration_status = Tetrahedral_mesh_smoother::Iteration_status;
    using Vertex_status = Tetrahedral_mesh_smoother::Vertex_status;
    using Cell_status = Tetrahedral_mesh_smoother::Tetrahedron_status;

    using Callback_function = std::function<bool (Iteration_status const &status,
                                                  Vertex_descriptor_map<Vertex_status> const &vertex_data,
                                                  Cell_descriptor_map<Cell_status> const &cell_data
                                                 )>;
    using Callback_setting = Tetrahedral_mesh_smoother::DEBUG_CALLBACK_SETTING;
    void set_callback_function(Callback_function callback_function, Callback_setting setting = Callback_setting::OUTER_ITER);

private:
    // inputs
    TetrahedralMesh &_mesh;
    Vertex_descriptor_map<std::array<bool, 3>> _user_locks;
    BoundaryMesh const &_boundary;
    EdgeNetwork const &_edge_network;
    std::vector<std::tuple<Vertex_descriptor, Point_3, double>> _vertex_target_positions;

    // options
    bool _verbose = false;
    bool _lock_boundary = true;
    Parameters::OPTIMIZATION_MODE _optimization_mode = Parameters::CONFORMAL;
    Parameters::PREDICATES_MODE _predicates_mode = Parameters::NO_CHECK;

    unsigned _max_number_of_iteration = 1000;
    double _min_valid_edge_size = 1e-6;
    double _boundary_weight = 1.;

    bool is_vert_locked(Vertex_descriptor v) const;

    template<typename T>
    Eigen::Vector3d convert_to_eigen(T const &vector) const { return Eigen::Vector3d(vector[0], vector[1], vector[2]); }

    inline Point_3 convert_to_user(Eigen::Vector3d point) const { point = _scale*point + _shift; return Point_3(point[0], point[1], point[2]); }

    inline Eigen::Vector3d convert_to_inner(Point_3 const &point) const { return (convert_to_eigen(point) - _shift) / _scale; }


    // internal working data
    void check_refs();
    void clear_internal_data();
    void initialize_boundary();
    void initialize_curve_network();
    void create_compress_sorted_data();
    void initialise_point_targets();
    Eigen::VectorXd _compressed_coords;
    std::vector<bool> _compressed_locks;
    Vertex_descriptor_map<unsigned> _vertex_original_to_compressed;
    Cell_descriptor_map<unsigned> _cell_original_to_compressed;
    std::vector<std::array<unsigned, 4>> _tetrahedra;
    std::vector<std::array<Eigen::Vector3d, 4>> _tetrahedron_refs;
    std::vector<std::vector<unsigned>> _vert2tet_corner;


    enum QUERY_TYPE {NONE, POINT_QUERY, ELEMENT_QUERY, BATCH_POINT_QUERY, BATCH_ELEMENT_QUERY};

    std::vector<std::vector<unsigned>> _bnd_faces;
    std::vector<std::pair<unsigned, std::vector<std::array<unsigned, 2>>>> _vert_and_face_corners;
    std::vector<Surface_patch_index> _face_surface_id;
    QUERY_TYPE _boundary_query_type = NONE;
    Boundary_point_query _boundary_point_query = nullptr;
    Boundary_polygon_query _boundary_polygon_query = nullptr;
    Boundary_point_batch_query _boundary_point_batch_query = nullptr;
    Boundary_polygon_batch_query _boundary_polygon_batch_query = nullptr;
    std::vector<Point_3> _boundary_batch_info_points;
    std::vector<double> _boundary_batch_info_radii;
    std::vector<std::vector<Point_3>> _boundary_batch_info_polygons;
    std::vector<Plane> _boundary_batch_planes;

    void initialise_boundary_query(Tetrahedral_mesh_smoother &);


    std::vector<std::array<unsigned, 2>> _curve_edges;
    std::vector<Curve_index> _curve_ids;
    QUERY_TYPE _curve_query_type = NONE;
    Curve_point_query _curve_point_query = nullptr;
    Curve_segment_query _curve_segment_query = nullptr;
    Curve_point_batch_query _curve_point_batch_query = nullptr;
    Curve_segment_batch_query _curve_segment_batch_query = nullptr;
    std::vector<Point_3> _curve_batch_info_points;
    std::vector<double> _curve_batch_info_radii;
    std::vector<std::array<Point_3, 2>> _curve_batch_info_edges;
    std::vector<Curve_tangent> _curve_batch_tangents;

    void initialise_curve_queries(Tetrahedral_mesh_smoother &);


    std::vector<std::tuple<unsigned, Eigen::Vector3d, double>> _point_targets;


    QUERY_TYPE _target_size_query_type = NONE;
    Point_target_sizing_query _target_size_query = nullptr; 
    Point_target_sizing_batch_query _target_size_batch_query = nullptr;

    Update_validator _update_validator = nullptr;
    Vertex_descriptor_map<Point_3> _current_coords_to_check;


    double _scale = 1.;
    Eigen::Vector3d _shift = Eigen::Vector3d::Zero();
    double const _rescaled_bbox_scale = 0.5;
    void rescale_geometry();
    void update_mesh_coordinates();


    Callback_function _callback_function = nullptr;
    Callback_setting _callback_setting = Callback_setting::NOTHING;
    bool _callback_initialized = false;
    Iteration_status _callback_status;
    Vertex_descriptor_map<Vertex_status> _callback_vertex_map_data;
    Cell_descriptor_map<Cell_status> _callback_cell_map_data;
    void initialize_callback();
    bool run_callback(Iteration_status const &, std::vector<Vertex_status> const &, std::vector<Cell_status> const&);

    unsigned _nb_lbfgs_iterations = 0;
    unsigned _nb_predicates_invalid_steps = 0;

    void initialise_smoother(Tetrahedral_mesh_smoother &);
};


using cgal_types::C3t3_wrapper; 



/*!
* \ingroup pkgMeshSmoothing3Classes
* 
* \brief smooth a c3t3 mesh structure
*
* The class `C3t3_smoother` provides a wrapper around `Mesh_smoother` to use it with a CGAL `C3t3` mesh structure.
* Facets of the `C3t3` mesh structure are defined as boundary, and edges of the `C3t3` mesh structure are used as curve network.
* It is still required to define the wished boundary queries. 
* 
* \tparam Type of the input mesh complex, model of `MeshComplex_3InTriangulation_3`
* 
* \sa `CGAL::Mesh_smoothing_3::Mesh_smoother`
* 
*/
template<
    typename C3t3
>
class C3t3_smoother : public Mesh_smoother <C3t3_wrapper<C3t3>, C3t3_wrapper<C3t3>, C3t3_wrapper<C3t3>> {
private:
    C3t3_wrapper<C3t3> mesh_wrapper;

public:
    C3t3_smoother(C3t3 &c3t3)
    : mesh_wrapper{c3t3}
    , Mesh_smoother<C3t3_wrapper<C3t3>, C3t3_wrapper<C3t3>, C3t3_wrapper<C3t3>>(mesh_wrapper, mesh_wrapper, mesh_wrapper)
    {}

};


} // namespace Mesh_smoothing_3



#include "internal/Mesh_smoothing_3_impl.hpp"
