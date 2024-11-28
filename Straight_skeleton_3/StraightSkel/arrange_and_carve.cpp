#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <fstream>

namespace PMP = ::CGAL::Polygon_mesh_processing;
namespace pred = ::CGAL::Polygon_mesh_processing::Corefinement;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point = K::Point_3;
using Vector = K::Vector_3;
using Triangle = K::Triangle_3;
using Plane = K::Plane_3;

using Mesh = CGAL::Surface_mesh<Point>;

using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

using PID = std::size_t;
using TID = std::size_t;
using VID = std::size_t;

enum class Volume_orientation {
  UNKNOWN = 0,
  UNREACHABLE, // 1
  INWARD, // 2
  OUTWARD, // 3
  INCONSISTENT // 4
};

auto sq_distance = [](const TID tid,
                      const std::vector<Point> pts,
                      const std::vector<std::vector<PID> > tris,
                      const Mesh& mesh) -> FT
{
  // brute force but whatever
  FT sq_dist = std::numeric_limits<double>::max();
  for(face_descriptor f : faces(mesh))
  {
    Triangle tr { pts[tris[tid][0]],
                  pts[tris[tid][1]],
                  pts[tris[tid][2]] };

    halfedge_descriptor h = halfedge(f, mesh);
    Triangle f_tr = { get(CGAL::vertex_point, mesh, source(h, mesh)),
                      get(CGAL::vertex_point, mesh, target(h, mesh)),
                      get(CGAL::vertex_point, mesh, target(next(h, mesh), mesh)) };

    FT sq_ld = CGAL::squared_distance(tr, f_tr);
    if(sq_ld < sq_dist)
      sq_dist = sq_ld;
  }

  return sq_dist;
};

Mesh coplanar_merge(const Mesh& mesh)
{
  CGAL::IO::write_polygon_mesh("results/coplanar_merge_before.off", mesh,
                                CGAL::parameters::stream_precision(17));

  Mesh out;
  PMP::remesh_planar_patches(mesh, out);

  CGAL::IO::write_polygon_mesh("results/coplanar_merge_after.off", out,
                                CGAL::parameters::stream_precision(17));

  return out;
}

template <typename ValueType,
          typename VisitorBase = CGAL::Polygon_mesh_processing::Autorefinement::Default_visitor>
struct Range_updating_autoref_visitor : public VisitorBase {
    Range_updating_autoref_visitor(const std::vector<ValueType>& old_range,
                                   std::vector<ValueType>& new_range,
                                   const VisitorBase& base = VisitorBase{})
        : VisitorBase(base), old_range_(old_range), new_range_(new_range) {
        new_range.reserve(old_range.size());
    }

    void verbatim_triangle_copy(std::size_t tgt_id, std::size_t src_id) {
        // std::cout << "verbatim_triangle_copy " << tgt_id << " from " << src_id << std::endl;
        VisitorBase::verbatim_triangle_copy(tgt_id, src_id);
        new_range_.resize(tgt_id + 1);
        new_range_[tgt_id] = old_range_[src_id];
    }

    void new_subtriangle (std::size_t tgt_id, std::size_t src_id) {
        // std::cout << "new_subtriangle " << tgt_id << " from " << src_id << std::endl;
        VisitorBase::new_subtriangle(tgt_id, src_id);
        new_range_.resize(tgt_id + 1);
        new_range_[tgt_id] = old_range_[src_id];
    }

private:
    const std::vector<ValueType>& old_range_;
    std::vector<ValueType>& new_range_;
};

template <typename ValueType,
          typename BaseVisitor = CGAL::Polygon_mesh_processing::internal::Default_repair_PS_visitor>
struct Range_updating_repair_PS_visitor : public BaseVisitor {
    Range_updating_repair_PS_visitor(std::vector<ValueType>& range,
                                     const BaseVisitor& base_visitor = BaseVisitor{})
        : BaseVisitor(base_visitor), range_(range) { }

    void swap(std::size_t pos_1, std::size_t pos_2) {
        BaseVisitor::swap(pos_1, pos_2);
        std::swap(range_[pos_1], range_[pos_2]);
    }
    void resize(std::size_t new_size) {
        BaseVisitor::resize(new_size);
        range_.resize(new_size);
    }
private:
    std::vector<ValueType>& range_;
};

void generate_arrangement(const Mesh& sm,
                          const FT delta,
                          std::vector<Point>& points,
                          std::vector<std::vector<PID> >& triangles,
                          std::vector<std::array<VID, 2> >& face_volume_IDs)
{
  int piece_id = 0;

  auto add_face_to_soup = [&](const face_descriptor f,
                              const Mesh& sm2) -> void
  {
    std::vector<PID> face;

    // For each vertex in the face, add the point to the points if it’s not already there,
    // and get the index of the point in the points to store in 'face'.
    for(vertex_descriptor v : CGAL::vertices_around_face(halfedge(f, sm2), sm2))
    {
      const Point& p = get(CGAL::vertex_point, sm2, v);

      // Check if this point is already in points to avoid duplicates
      auto it = std::find(points.begin(), points.end(), p);
      if(it == points.end()) {
          points.push_back(p);
          face.push_back(points.size() - 1); // index of new point
      } else {
          face.push_back(std::distance(points.begin(), it)); // index of existing point
      }
    }

    // the purpose is to avoid duplicates AND give priority to the boundary faces
    // as to have the proper 0-1 markers
    // @todo obviously everything has terrible complexity
    std::vector<PID> canonical_face = PMP::internal::construct_canonical_polygon<K>(points, face);
    auto polygon_equal = [&](const std::vector<PID>& other_face)
    {
      return (canonical_face == PMP::internal::construct_canonical_polygon<K>(points, other_face));
    };

    if(std::find_if(std::cbegin(triangles), std::cend(triangles), polygon_equal) != std::cend(triangles)) {
      std::cerr << "Warning: not inserting face as it is a duplicate" << std::endl;
      return;
    }

    triangles.push_back(face);
  };

  // add the (triangle) face to the output
  for(face_descriptor fd : faces(sm))
  {
    add_face_to_soup(fd, sm);
    face_volume_IDs.emplace_back();
    face_volume_IDs.back()[0] = 1; // down, inside
    face_volume_IDs.back()[1] = 0; // up, outside

    CGAL::IO::write_OFF("results/contributions_" + std::to_string(piece_id++) + ".off",
                        points, std::vector<std::vector<PID> >{triangles.back()},
                        CGAL::parameters::stream_precision(17));
  }

  for(face_descriptor f : faces(sm))
  {
    halfedge_descriptor h = halfedge(f, sm);
    vertex_descriptor v = source(h, sm);
    const Point& p = get(CGAL::vertex_point, sm, v);

    Vector normal = CGAL::Polygon_mesh_processing::compute_face_normal(f, sm);
    Plane plane(p, normal);

    // Shift the plane along its normal by the offset delta
    Plane shifted_plane(p - delta * normal, normal);

    // Clip the mesh with the shifted plane
    Mesh clipped_sm = sm;
    PMP::clip(clipped_sm, shifted_plane, CGAL::parameters::clip_volume(true));

    CGAL::IO::write_OFF("results/clipped_" + std::to_string(piece_id) + ".off", clipped_sm,
                        CGAL::parameters::stream_precision(17));

    // Extract faces that lie exactly on the shifted plane from the clipped mesh
    // @todo maybe can do something smarter with the visitor
    std::size_t old_triangles_count = triangles.size();
    for(face_descriptor cf : faces(clipped_sm)) {
      bool all_vertices_on_plane = true;
      for(vertex_descriptor cv : vertices_around_face(halfedge(cf, clipped_sm), clipped_sm)) {
        if(!shifted_plane.has_on(get(CGAL::vertex_point, clipped_sm, cv))) {
          all_vertices_on_plane = false;
          break;
        }
      }

      if(!all_vertices_on_plane) {
        continue;
      }

      add_face_to_soup(cf, clipped_sm);
    }

    std::size_t new_triangles_count = triangles.size();
    std::cout << new_triangles_count - old_triangles_count << " faces in contribution" << std::endl;

    for(TID tid=old_triangles_count; tid<new_triangles_count; ++tid) {
      // duplicate check means we don't override existing face_volume_IDs
      face_volume_IDs.emplace_back();
      face_volume_IDs.back()[0] = 1; // down, inside
      face_volume_IDs.back()[1] = 1; // up, inside
    }

    std::vector<std::vector<PID> > contributions(std::cbegin(triangles) + old_triangles_count,
                                                 std::cbegin(triangles) + new_triangles_count);
    CGAL::IO::write_OFF("results/contributions_" + std::to_string(piece_id++) + ".off",
                        points, contributions,
                        CGAL::parameters::stream_precision(17));
  }

  // for(TID tid=0; tid<triangles.size(); ++tid)
  // {
  //   std::cout << "Triangle #" << tid << " markers: " << face_volume_IDs[tid][0] << " " << face_volume_IDs[tid][1] << std::endl;
  // }

  std::cout << "== autoref ==" << std::endl;

  // Now, autoref and use visitors to maintain the outside / inside tag
  std::vector<std::array<VID, 2> > updated_face_volume_IDs;
  Range_updating_autoref_visitor<std::array<VID, 2> > autoref_visitor(face_volume_IDs, updated_face_volume_IDs);

  PMP::autorefine_triangle_soup(points, triangles,
                                CGAL::parameters::visitor(autoref_visitor)
                                                  .concurrency_tag(CGAL::Parallel_if_available_tag()));

  face_volume_IDs = std::move(updated_face_volume_IDs);

  CGAL::IO::write_OFF("results/autorefed.off", points, triangles,
                      CGAL::parameters::stream_precision(17));

  std::cout << "== purge duplicates ==" << std::endl;

  // purge duplicate faces
  Range_updating_repair_PS_visitor<std::array<VID, 2> > repair_ps_visitor(face_volume_IDs);
  PMP::merge_duplicate_polygons_in_polygon_soup(points, triangles,
                                                CGAL::parameters::visitor(repair_ps_visitor)
                                                                 .erase_all_duplicates(false) /*keep one*/
                                                                 .require_same_orientation(false));

  CGAL::IO::write_OFF("results/autorefed_cleaned.off", points, triangles,
                      CGAL::parameters::stream_precision(17));
}

void carve_arrangement(const Mesh& sm,
                       const FT delta,
                       const std::vector<Point>& points,
                       const std::vector<std::vector<PID> >& triangles,
                       std::vector<std::array<VID, 2> >& face_volume_IDs)
{
  auto fill_edge_map = [](const std::vector<Point>& points,
                          const std::vector<std::vector<PID> >& triangles,
                          std::vector<std::unordered_map<PID, std::vector<TID> > >& edge_map) {
      CGAL_precondition(edge_map.size() == points.size());

      // collect duplicated edges
      for (TID ti=0; ti<triangles.size(); ++ti) {
          for (std::size_t j=0; j<3; ++j) {
              std::pair<PID, PID> e_pids = CGAL::make_sorted_pair(triangles[ti][j],
                                                                  triangles[ti][(j+1)%3]);
              edge_map[e_pids.first][e_pids.second].push_back(ti);
          }
      }

      namespace pred = CGAL::Polygon_mesh_processing::Corefinement;

      for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
          for (auto& pid1_and_edges : edge_map[pid0]) {
              std::vector<TID>& inc_triangles = pid1_and_edges.second;

              if (inc_triangles.size() == 2) { // edge is only incident to a single SS3 face
                  continue;
              }

              const PID pid1 = pid1_and_edges.first;
              CGAL_assertion(pid0 != pid1);

              // std::cout << "processing a non-manifold edge: " << points[pid0] << " -- " << points[pid1] << std::endl;

              auto get_third_point_id = [&triangles, pid0, pid1](TID tid) -> PID
              {
                  std::size_t third;

                  // need to be careful that the orientation of the edge might not match the orientation of the triangle
                  if (triangles[tid][0] == pid0 || triangles[tid][0] == pid1) {
                      if (triangles[tid][1] == pid0 || triangles[tid][1] == pid1) {
                          third = triangles[tid][2];
                      } else {
                          third = triangles[tid][1];
                      }
                  } else {
                      third = triangles[tid][0];
                  }

                  CGAL_postcondition(third != pid0 && third != pid1);
                  return third;
              };

              const Point& ref_pt = points.at(get_third_point_id(inc_triangles[0]));
              auto less = [&ref_pt, &points, pid0, pid1, get_third_point_id](TID tid1, TID tid2)
              {
                  return pred::sorted_around_edge<K>(points.at(pid0), points.at(pid1),
                                                     ref_pt,
                                                     points.at(get_third_point_id(tid1)),
                                                     points.at(get_third_point_id(tid2)));
              };

              std::sort(inc_triangles.begin()+1, inc_triangles.end(), less);

              // std::cout << "Around edge [" << pid0 << " " << pid1 << "], faces are sorted: ";
              // for(TID tid : inc_triangles)
              //   std::cout << " " << tid;
              // std::cout << std::endl;
          }
      }
  }; // lambda 'fill_edge_map'

  auto build_volume_CC = [](const TID seed_tid,
                            const VID CC_ID,
                            const bool start_from_inverted_face,
                            const bool ignore_facet_with_incompatible_orientations,
                            const bool ignore_dangling_outside_facets,
                            const std::vector<Point>& points,
                            const std::vector<std::vector<PID> >& triangles,
                            const auto& edge_map,
                            auto& volume_CCs,
                            auto& volume_orientations,
                            auto& face_volume_IDs) {

#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
      std::cout << "Building volume #" << CC_ID << " from seed face " << seed_tid << std::endl;
#endif

      volume_CCs.emplace_back();
      volume_orientations.emplace_back();

      std::stack<std::pair<TID, bool> > to_visit;
      to_visit.emplace(seed_tid, start_from_inverted_face);

      while (!to_visit.empty())
      {
          TID current_tid;
          bool invert_face;
          std::tie(current_tid, invert_face) = to_visit.top();
          to_visit.pop();

#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
          std::cout << "At face " << current_tid << " [" << triangles[current_tid][0]
                                                  << ", " << triangles[current_tid][1]
                                                  << ", " << triangles[current_tid][2] << "], ";
          std::cout << "invert: " << invert_face << ", ";
          std::cout << "VIDS: " << face_volume_IDs[current_tid][0] << " " << face_volume_IDs[current_tid][1] << std::endl;
#endif

          std::size_t pos = invert_face ? 0 : 1;
          if(face_volume_IDs[current_tid][pos] == CC_ID) {
              // already visited this facet during the flooding of this volume's boundary
              continue;
          }

          volume_CCs.back().push_back(current_tid);

          // mark face as visited
          face_volume_IDs[current_tid][pos] = CC_ID;

          // flood through the edges
          for (int j=0; j<3; ++j) {
              std::pair<PID, PID> e_pids = CGAL::make_sorted_pair(triangles[current_tid][j],
                                                                  triangles[current_tid][(j+1)%3]);
              const std::vector<TID>& inc_triangles = edge_map.at(e_pids.first).at(e_pids.second);
              CGAL_assertion(!inc_triangles.empty());

              // The faces are ordered CCW while looking from pid0.
              // So the walking while looking from [j] depends on whether [j] is pid0 or not
              int iter_direction = (e_pids.first == triangles[current_tid][j]) ? 1 : -1;

              // and it also depends on whether we are walking above or below the face
              iter_direction *= invert_face ? 1 : -1;

#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
              std::cout << "  ~~ Crossing edge [" << e_pids.first << ", " << e_pids.second << "]" << std::endl;
              std::cout << "    pos: " << points[e_pids.first] << " " << points[e_pids.second] << std::endl;
              std::cout << "    iter_direction = " << iter_direction << std::endl;
              std::cout << "    invert_face = " << invert_face << std::endl;
#endif

              TID next_tid = current_tid;
              for (;;) {
                  if(inc_triangles.size() == 1) {
#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
                      std::cerr << "Warning: dangling triangle..." << std::endl;
                      std::cout << "    over the edge, the triangle is ITSELF " << current_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                      to_visit.emplace(current_tid, !invert_face);
#endif
                      break;
                  } else if(inc_triangles.size() == 2) {
                      // we should only be there once, meaning if we do not ignore orientations,
                      // then the faces MUST be compatible
                      CGAL_assertion(next_tid == current_tid);

                      next_tid = (inc_triangles[0] == current_tid) ? inc_triangles[1] : inc_triangles[0];
#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
                      std::cout << "    over the edge, the triangle is TRIVIALLY " << next_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                      std::cout << "VIDS " << face_volume_IDs[next_tid][0] << " " << face_volume_IDs[next_tid][1] << std::endl;
#endif
                      CGAL_assertion(next_tid != current_tid);
                  } else {
                      // tricky part, now
                      auto tid_it = std::find(std::begin(inc_triangles), std::end(inc_triangles), next_tid /*updates on every iteration*/);
                      CGAL_assertion(tid_it != inc_triangles.end());

                      if(iter_direction == 1) { // CCW
#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
                          std::cout << "    CCW walk" << std::endl;
#endif
                          auto next_it = std::next(tid_it);
                          next_tid = (next_it == inc_triangles.end()) ? inc_triangles[0] : *next_it;
                      } else { // CW
#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
                          std::cout << "    CW walk" << std::endl;
#endif
                          next_tid = (tid_it == inc_triangles.begin()) ? inc_triangles.back() : *(std::prev(tid_it));
                      }

#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
                      std::cout << "    over the edge, the triangle is " << next_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                      std::cout << "VIDS " << face_volume_IDs[next_tid][0] << " " << face_volume_IDs[next_tid][1] << std::endl;
#endif
                      CGAL_assertion(next_tid != current_tid);
                  }

                  // If the next face is incident to the outside (CC_ID == 0) on both sides, ignore it
                  if(ignore_dangling_outside_facets) {
                      // CC id #0 is the outside marker for volume building of facet prisms
                      bool is_dangling = (face_volume_IDs[next_tid][0] == 0 && face_volume_IDs[next_tid][1] == 0);
                      if(is_dangling) {
                          CGAL_assertion(inc_triangles.size() != 2);
#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
                          std::cout << "  ignoring next TID because it is dangling" << std::endl;
#endif
                          continue;
                      }
                  }

                  // If the edge has the same direction in both faces (aka, the orientation changes),
                  // then we have to flip the direction of turning around the edge)
                  const auto j_it = std::find(std::begin(triangles[next_tid]),
                                            std::end(triangles[next_tid]),
                                            triangles[current_tid][j]);
                  CGAL_assertion(j_it != std::end(triangles[next_tid]));
                  const std::size_t pos = std::distance(std::begin(triangles[next_tid]), j_it);
                  CGAL_assertion(triangles[next_tid][pos] == triangles[current_tid][j]);

                  const bool flip_side = (triangles[next_tid][(pos+1)%3] == triangles[current_tid][(j+1)%3]);
#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
                  std::cout << "    flipping? " << flip_side
                            << " (N: " << triangles[next_tid][(pos+1)%3] << " C: " << triangles[current_tid][(j+1)%3] << ")" << std::endl;
#endif
                  if(ignore_facet_with_incompatible_orientations && flip_side) {
                      if(inc_triangles.size() != 2) { // @tmp
#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
                          std::cout << "  ignoring next TID because of incompatible orientation" << std::endl;
#endif
                          continue; // keep turning
                      } else {
                          std::cerr << "Warning: we should not have ignored incompatible face orientations, but there are only 2 incident triangles" << std::endl;
                      }
                  }

#ifdef CGAL_SS3_DEBUG_VOLUME_BUILDING
                  std::cout << "Final TID = " << next_tid << std::endl;
#endif
                  if(flip_side) {
                      CGAL_warning(!ignore_facet_with_incompatible_orientations);
                      volume_orientations[CC_ID] = Volume_orientation::INCONSISTENT;
                      to_visit.emplace(next_tid, !invert_face);
                  } else {
                      to_visit.emplace(next_tid, invert_face);
                  }

                  break;
              }
          }
      }
  }; // lambda 'build_volume_CC'

  std::vector<std::unordered_map<PID, std::vector<TID> > > edge_map(points.size());
  fill_edge_map(points, triangles, edge_map);

  // remove all volumes incident to the input
  std::vector<std::vector<TID> > unused_volume_CCs;
  std::vector<Volume_orientation> unused_volume_orientations;

  std::vector<TID> triangles_to_treat;
  for (TID tid=0; tid<triangles.size(); ++tid)
  {
    CGAL_assertion(face_volume_IDs[tid][0] != 0); // down shouldn't be outside
    if(face_volume_IDs[tid][1] == 0) // up == outside
      triangles_to_treat.push_back(tid);
  }

  for (TID tid : triangles_to_treat)
  {
    if(face_volume_IDs[tid][0]) // already visited
      continue;

    std::cout << "remove volume incident to face " << tid << std::endl;
    build_volume_CC(tid, 0 /*outside*/,
                    true /*start from inveted face (down)*/,
                    false /*do not ignore faces with incompatible orientations*/,
                    false /*do not ignore dangling facets*/,
                    points,
                    triangles,
                    edge_map,
                    unused_volume_CCs,
                    unused_volume_orientations,
                    face_volume_IDs);

    CGAL_postcondition(face_volume_IDs[tid][0] == 0 && face_volume_IDs[tid][1] == 0);
  }

  // debug
  std::vector<std::vector<PID> > not_outside_faces;
  for(TID tid=0; tid<triangles.size(); ++tid)
  {
    if(face_volume_IDs[tid][0] == face_volume_IDs[tid][1] && face_volume_IDs[tid][0] == 0)
      continue;

    not_outside_faces.push_back(triangles[tid]);
  }

  std::cout << not_outside_faces.size() << " outside faces after border purge" << std::endl;
  CGAL::IO::write_OFF("results/carving_-1.off", points, not_outside_faces,
                      CGAL::parameters::stream_precision(17));
  // -- debug (end)

  unsigned int carving_iter = 0;
  for(;;) // till there is nothing left to do
  {
    bool did_something = false;

    // remove boundary volumes whose boundary face is at distance < offset than the input
    for (TID tid=0; tid<triangles.size(); ++tid)
    {
      if(face_volume_IDs[tid][0] == face_volume_IDs[tid][1]) // fully outside of fully inside
        continue;

      for(std::size_t i=0; i<2; ++i)
      {
        if(face_volume_IDs[tid][i] != 0)
          continue;

        FT sq_d = sq_distance(tid, points, triangles, sm);
        std::cout << "Distance of Triangle #" << tid << " = " << sq_d << std::endl;

        if(sq_d > CGAL::square(delta) - 1e-10) // @fixme hardcoded bound
          continue;

        std::cout << "Carving through " << tid << " at dist " << CGAL::approximate_sqrt(sq_d) << std::endl;

        build_volume_CC(tid, 0 /*outside*/,
                        bool(i), // if outside is down ('0'), start upwards ('false'), and inversely
                        false /*do not ignore faces with incompatible orientations*/,
                        false /*do not ignore dangling facets*/,
                        points,
                        triangles,
                        edge_map,
                        unused_volume_CCs,
                        unused_volume_orientations,
                        face_volume_IDs);

        did_something = true;
      }
    }

    // remove boundary volumes whose boundary face is not facing outwards
    // @todo:
    // - can this happen?
    // - what does it mean in terms of distance?
    for (TID tid=0; tid<triangles.size(); ++tid)
    {
      if(face_volume_IDs[tid][0] == face_volume_IDs[tid][1]) // fully outside of fully inside
        continue;

      for(std::size_t i=0; i<2; ++i)
      {
        // configuration of down outside, up not outside (0 == 0 && 1 != 0)
        if(face_volume_IDs[tid][0] != 0 || face_volume_IDs[tid][1] == 0)
          continue;

        std::cout << "Carving through inverted " << tid << std::endl;

        build_volume_CC(tid, 0 /*outside*/,
                        false, // start upwards
                        false /*do not ignore faces with incompatible orientations*/,
                        false /*do not ignore dangling facets*/,
                        points,
                        triangles,
                        edge_map,
                        unused_volume_CCs,
                        unused_volume_orientations,
                        face_volume_IDs);

        did_something = true;
      }
    }

    // debug
    not_outside_faces.clear();
    for(TID tid=0; tid<triangles.size(); ++tid)
    {
      if(face_volume_IDs[tid][0] == face_volume_IDs[tid][1] && face_volume_IDs[tid][0] == 0)
        continue;

      not_outside_faces.push_back(triangles[tid]);
    }

    std::cout << not_outside_faces.size() << " NOT outside faces at carving iteration #" << carving_iter << std::endl;
    CGAL::IO::write_OFF("results/carving_" + std::to_string(carving_iter++) + ".off",
                        points, not_outside_faces,
                        CGAL::parameters::stream_precision(17));
    // debug (end)

    if(!did_something)
      break;
  }

  for(TID tid=0; tid<triangles.size(); ++tid)
  {
    std::cout << "Triangle #" << tid << " markers: " << face_volume_IDs[tid][0] << " " << face_volume_IDs[tid][1] << std::endl;
  }
}

void extract_outer_hull(const std::vector<Point>& points,
                        const std::vector<std::vector<PID> >& triangles,
                        const std::vector<std::array<VID, 2> >& face_volume_IDs,
                        const Mesh& sm)
{
  // Take all the faces that are incident to a boundary
  std::vector<std::vector<PID> > boundary_triangles;
  for(TID tid=0; tid<triangles.size(); ++tid)
  {
    if(face_volume_IDs[tid][0] == face_volume_IDs[tid][1])
      continue;

    std::cout << tid << " is boundary" << std::endl;
    boundary_triangles.push_back(triangles[tid]);
  }

  for(std::size_t tid=0; tid<boundary_triangles.size(); ++tid)
  {
    std::cout << "Final Face #" << tid << " at distance " << CGAL::approximate_sqrt(sq_distance(tid, points, boundary_triangles, sm)) << std::endl;
  }

  CGAL::IO::write_OFF("results/hull.off", points, boundary_triangles,
                      CGAL::parameters::stream_precision(17));
}

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: "
              << argv[0] << "\n"
              << "\tinput_filename\n"
              << "\t[output_filename] (PLY)\n"
              << std::endl;
    std::cerr << "Input format: any format readable with CGAL::IO::read_polygon_mesh()" << std::endl;
    std::cerr << "Output format: PLY" << std::endl;
    return EXIT_FAILURE;
  }

  const char* input_filename = argv[1];
  const FT delta = (argc > 2) ? std::stod(argv[2]) : 1;

  std::cout << "in: " << input_filename << std::endl;

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(input_filename, sm)) {
    std::cerr << "Error: failed to read input" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh: " << num_vertices(sm) << " NV " << num_faces(sm) << " NF" << std::endl;

  if(!CGAL::is_closed(sm)) {
    std::cerr << "Error: input mesh must be closed" << std::endl;
    return EXIT_FAILURE;
  }

  // sm = coplanar_merge(sm);
  PMP::triangulate_faces(sm); // autoref wants triangulated inputs

  CGAL::IO::write_OFF("results/input.off", sm, CGAL::parameters::stream_precision(17));

  std::cout << "Input mesh: " << num_vertices(sm) << " NV " << num_faces(sm) << " NF (post remesh)" << std::endl;

  std::vector<Point> points;
  std::vector<std::vector<PID> > faces;
  std::vector<std::vector<TID> > volume_CCs; // range of range (volume) of triangle IDs
  std::vector<std::array<VID, 2> > face_volume_IDs;


  generate_arrangement(sm, delta, points, faces, face_volume_IDs);

  carve_arrangement(sm, delta, points, faces, face_volume_IDs);

  extract_outer_hull(points, faces, face_volume_IDs, sm);

  CGAL::IO::write_PLY("results/final.ply", sm,
                      CGAL::parameters::use_binary_mode(false)
                                       .stream_precision(17));

  std::cout << "Done." << std::endl;

  return EXIT_SUCCESS;
}
