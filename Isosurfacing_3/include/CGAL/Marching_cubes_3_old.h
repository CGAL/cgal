#ifndef CGAL_MARCHING_CUBES_3_OLD_H
#define CGAL_MARCHING_CUBES_3_OLD_H

#include <CGAL/Isosurfacing_3/internal/Marching_cubes_3_internal.h>
#include <CGAL/Isosurfacing_3/internal/Tables_old.h>

#include <mutex>
#include <unordered_map>

namespace CGAL {
namespace Isosurfacing {

template <class Domain_, class PointRange, class PolygonRange>
void marching_cubes_cell_old(const std::size_t x, const std::size_t y, const std::size_t z, const Domain_& domain,
                             const typename Domain_::FT iso_value, PointRange& points, PolygonRange& polygons,
                             std::mutex& mutex) {

    typedef std::array<std::size_t, 3> Idx_3;
    typedef typename Domain_::FT FT;
    typedef typename Domain_::Point_3 Point_3;

    const Idx_3 idx0 = {x + 0, y + 1, z + 0};
    const Idx_3 idx1 = {x + 1, y + 1, z + 0};
    const Idx_3 idx2 = {x + 1, y + 0, z + 0};
    const Idx_3 idx3 = {x + 0, y + 0, z + 0};
    const Idx_3 idx4 = {x + 0, y + 1, z + 1};
    const Idx_3 idx5 = {x + 1, y + 1, z + 1};
    const Idx_3 idx6 = {x + 1, y + 0, z + 1};
    const Idx_3 idx7 = {x + 0, y + 0, z + 1};

    const Point_3 pos0 = domain.position(idx0[0], idx0[1], idx0[2]);
    const Point_3 pos1 = domain.position(idx1[0], idx1[1], idx1[2]);
    const Point_3 pos2 = domain.position(idx2[0], idx2[1], idx2[2]);
    const Point_3 pos3 = domain.position(idx3[0], idx3[1], idx3[2]);
    const Point_3 pos4 = domain.position(idx4[0], idx4[1], idx4[2]);
    const Point_3 pos5 = domain.position(idx5[0], idx5[1], idx5[2]);
    const Point_3 pos6 = domain.position(idx6[0], idx6[1], idx6[2]);
    const Point_3 pos7 = domain.position(idx7[0], idx7[1], idx7[2]);

    const FT dist0 = domain.value(idx0[0], idx0[1], idx0[2]);
    const FT dist1 = domain.value(idx1[0], idx1[1], idx1[2]);
    const FT dist2 = domain.value(idx2[0], idx2[1], idx2[2]);
    const FT dist3 = domain.value(idx3[0], idx3[1], idx3[2]);
    const FT dist4 = domain.value(idx4[0], idx4[1], idx4[2]);
    const FT dist5 = domain.value(idx5[0], idx5[1], idx5[2]);
    const FT dist6 = domain.value(idx6[0], idx6[1], idx6[2]);
    const FT dist7 = domain.value(idx7[0], idx7[1], idx7[2]);

    unsigned int cubeindex = 0;
    // set bit if corresponding corner is below iso
    cubeindex |= (dist0 < iso_value) << 0;
    cubeindex |= (dist1 < iso_value) << 1;
    cubeindex |= (dist2 < iso_value) << 2;
    cubeindex |= (dist3 < iso_value) << 3;
    cubeindex |= (dist4 < iso_value) << 4;
    cubeindex |= (dist5 < iso_value) << 5;
    cubeindex |= (dist6 < iso_value) << 6;
    cubeindex |= (dist7 < iso_value) << 7;

    Point_3 vertlist[12];
    // bitmap of edges that intersect the iso-plane(s)
    const int edges = edgeTable[cubeindex];

    if (edges == 0) {
        return;
    }

    // interpolation of vertices incident to the edge, according to diagram
    if (edges & (1 << 0)) {
        vertlist[0] = vertex_interpolation(pos0, pos1, dist0, dist1, iso_value);
    }
    if (edges & (1 << 1)) {
        vertlist[1] = vertex_interpolation(pos1, pos2, dist1, dist2, iso_value);
    }
    if (edges & (1 << 2)) {
        vertlist[2] = vertex_interpolation(pos2, pos3, dist2, dist3, iso_value);
    }
    if (edges & (1 << 3)) {
        vertlist[3] = vertex_interpolation(pos3, pos0, dist3, dist0, iso_value);
    }
    if (edges & (1 << 4)) {
        vertlist[4] = vertex_interpolation(pos4, pos5, dist4, dist5, iso_value);
    }
    if (edges & (1 << 5)) {
        vertlist[5] = vertex_interpolation(pos5, pos6, dist5, dist6, iso_value);
    }
    if (edges & (1 << 6)) {
        vertlist[6] = vertex_interpolation(pos6, pos7, dist6, dist7, iso_value);
    }
    if (edges & (1 << 7)) {
        vertlist[7] = vertex_interpolation(pos7, pos4, dist7, dist4, iso_value);
    }
    if (edges & (1 << 8)) {
        vertlist[8] = vertex_interpolation(pos0, pos4, dist0, dist4, iso_value);
    }
    if (edges & (1 << 9)) {
        vertlist[9] = vertex_interpolation(pos1, pos5, dist1, dist5, iso_value);
    }
    if (edges & (1 << 10)) {
        vertlist[10] = vertex_interpolation(pos2, pos6, dist2, dist6, iso_value);
    }
    if (edges & (1 << 11)) {
        vertlist[11] = vertex_interpolation(pos3, pos7, dist3, dist7, iso_value);
    }

    // std::lock_guard<std::mutex> lock(mutex);

    for (int i = 0; triTable[cubeindex][i] != -1; i += 3) {

        const Point_3& p0 = vertlist[triTable[cubeindex][i + 0]];
        const Point_3& p1 = vertlist[triTable[cubeindex][i + 1]];
        const Point_3& p2 = vertlist[triTable[cubeindex][i + 2]];

        const std::size_t p0_idx = points.size();  // TODO: not allowed

        points.push_back(p0);
        points.push_back(p1);
        points.push_back(p2);

        polygons.push_back({});
        auto& triangle = polygons.back();

        triangle.push_back(p0_idx + 0);
        triangle.push_back(p0_idx + 1);
        triangle.push_back(p0_idx + 2);
    }
}

template <class Domain_, class PointRange, class PolygonRange>
void marching_cubes_cell_RG(const std::size_t x, const std::size_t y, const std::size_t z, const Domain_& domain,
                            const typename Domain_::FT iso_value, PointRange& points, PolygonRange& polygons,
                            std::mutex& mutex, std::unordered_map<std::size_t, std::size_t>& vertex_map) {

    typedef typename Domain_::FT FT;
    typedef typename Domain_::Point_3 Point;

    /// The structure _Vertex_ represents a vertex by giving its unique global index and the unique index of the edge.
    struct Vertex {
        std::size_t g_idx;  //<! Index indicating the position in vertex array, used final shared vertex list.
        std::size_t g_edg;  //<! Unique global index used a key to find unique vertex list in the map.
    };

    // we need to compute up to 3 vertices at the interior of a cell, therefore
    // the cell shift factor is set to 3+3 = 6, i.e. 3 edges assigned to a cell for global numberig
    // and 3 vertices in the interior of the cell
    const int cell_shift_factor = 3;

    // there can be at most 12 intersections
    std::array<Vertex, 12> vertices;

    // slice hex
    // collect function values and build index
    FT values[8];
    Point corners[8];

    int vi = 0;
    std::bitset<8> index = 0;
    for (int kl = 0; kl <= 1; kl++) {
        for (int jl = 0; jl <= 1; jl++) {
            for (int il = 0; il <= 1; il++) {
                // collect scalar values and computex index
                corners[vi] = domain.position(x + il, y + jl, z + kl);
                values[vi] = domain.value(x + il, y + jl, z + kl);

                if (values[vi] >= iso_value) {
                    // index.set(VertexMapping[vi]);
                    index.set(vi);
                }
                // next cell vertex
                vi++;
            }
        }
    }

    // collect edges from table and
    // interpolate triangle vertex positon
    const int i_case = int(index.to_ullong());
    // compute for this case the vertices
    ushort flag = 1;
    for (int edge = 0; edge < 12; edge++) {
        if (flag & Cube_table::intersected_edges[i_case]) {
            // the edge global index is given by the vertex global index + the edge offset
            const std::size_t ix = x + Cube_table::global_edge_id[edge][0];
            const std::size_t iy = y + Cube_table::global_edge_id[edge][1];
            const std::size_t iz = z + Cube_table::global_edge_id[edge][2];
            const std::size_t global_index = (iz * domain.size_y() * domain.size_x() + iy * domain.size_x() + ix);
            vertices[edge].g_edg = cell_shift_factor * global_index + Cube_table::global_edge_id[edge][3];
            // generate vertex here, do not care at this point if vertex already exist
            // interpolation weight
            const int v0 = Cube_table::edge_to_vertex[edge][0];
            const int v1 = Cube_table::edge_to_vertex[edge][1];
            const FT l = (iso_value - values[v0]) / (values[v1] - values[v0]);
            // interpolate vertex
            const FT px = (1 - l) * corners[v0][0] + l * corners[v1][0];
            const FT py = (1 - l) * corners[v0][1] + l * corners[v1][1];
            const FT pz = (1 - l) * corners[v0][2] + l * corners[v1][2];
            const Point position(px, py, pz);

            // std::lock_guard<std::mutex> lock(mutex);
            // set vertex index
            // const auto s_index = vertex_map.find(vertices[edge].g_edg);
            // if (s_index == vertex_map.end()) {
            // index not found! Add index to hash map
            const std::size_t g_idx = points.size();
            // vertex_map[vertices[edge].g_edg] = g_idx;
            vertices[edge].g_idx = g_idx;
            points.push_back(position);
            //} else {
            //    vertices[edge].g_idx = s_index->second;  // this is vertex global index g_idx
            //}
        }
        flag <<= 1;
    }

    // std::lock_guard<std::mutex> lock(mutex);
    // construct triangles
    for (int t = 0; t < 16; t += 3) {
        const int t_index = i_case * 16 + t;
        // if (e_tris_list[t_index] == 0x7f)
        if (Cube_table::triangle_cases[t_index] == -1) break;

        const int eg0 = Cube_table::triangle_cases[t_index];
        const int eg1 = Cube_table::triangle_cases[t_index + 1];
        const int eg2 = Cube_table::triangle_cases[t_index + 2];

        // insert new triangle in list
        polygons.push_back({});
        auto& triangle = polygons.back();

        triangle.push_back(vertices[eg0].g_idx);
        triangle.push_back(vertices[eg1].g_idx);
        triangle.push_back(vertices[eg2].g_idx);
    }
}

template <class Domain_, class PointRange, class PolygonRange>
void make_triangle_mesh_using_marching_cubes_old(const Domain_& domain, const typename Domain_::FT iso_value,
                                                 PointRange& points, PolygonRange& polygons) {

    std::mutex mutex;

    const std::size_t size_k = domain.size_z();
    const std::size_t size_j = domain.size_y();
    const std::size_t size_i = domain.size_x();

    const std::size_t blocking_size = 100;

    // TODO: look at polygon mesh processing for tbb (also linking)

    // compute a unique global index for vertices
    // use as key the unique edge number
    std::unordered_map<std::size_t, std::size_t> v_map;


    for (std::size_t bj = 0; bj < size_j - 1; bj += blocking_size) {
        //#pragma omp parallel for
        for (std::size_t k = 0; k < size_k - 1; k++) {

            const std::size_t j_start = bj;
            const std::size_t j_end = std::min(size_j - 1, bj + blocking_size);

            for (std::size_t j = j_start; j < j_end; j++) {
                for (std::size_t i = 0; i < size_i - 1; i++) {
                    // internal::marching_cubes_cell_old(i, j, k, domain, iso_value, points, polygons, mutex);
                    internal::marching_cubes_cell_RG(i, j, k, domain, iso_value, points, polygons, mutex, v_map);
                }
            }
        }
    }
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_MARCHING_CUBES_3_OLD_H
