// Copyright (c) 2023-2026 GeometryFactory and Claudio Mancinelli.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Claudio Mancinelli and Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_REFINE_MESH_ALONG_PATHS_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_REFINE_MESH_ALONG_PATHS_H

#include <CGAL/license/Vector_graphics_on_surfaces.h>

#include <CGAL/Polygon_mesh_processing/locate.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

namespace CGAL {
namespace Vector_graphics_on_surfaces {

/*!
 * \ingroup VGSFunctions
 * refines `tmesh` so that each path in `paths` corresponds to a set edges of `tmesh` after the call.
 * Note that each path must be such that for two consecutive face locations, there exists a face in `tmesh` containing the two corresponding points.
 * \tparam TriangleMesh a model of `MutableFaceGraph`
 * \tparam K a model of `Kernel` with `K::FT` being a floating point number type (float or double)
 * \tparam VNM a model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `Kernel::Vector_3` as value type.
 * \tparam FNM a model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type and `Kernel::Vector_3` as value type.
 * \tparam OutputIterator an output iterator accepting `boost::graph_traits<TriangleMesh>::%halfedge_descriptor` to be put
 * \param tmesh the triangle mesh to be refined
 * \param paths a path described as a range of edge locations, with the property that
 *              for two consecutive edge locations, there exists a face in `tmesh` containing the two corresponding points.
 * \param vnm property map associating a normal to each vertex of `tmesh` that is updated by this function
 * \param fnm property map associating a normal to each face of `tmesh` that is updated by this function
 * \param out output iterator where created halfedges are put
 * \todo add named parameters
 * \todo vnm and fnm are optional
 * \todo intersection between path or self-intersections are not handle. Should it be? If not what do we do?
 * \todo out should contain edges, and also existing edges already part of a path
 * \todo shall we also have the edges in the order of the input rather than all at once
 */
template <class K, class TriangleMesh, class VNM, class FNM, class OutputIterator>
// std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>
OutputIterator
refine_mesh_along_paths(const std::vector<std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>>>& paths,
                        TriangleMesh& tmesh,
                        VNM vnm,
                        FNM fnm,
                        OutputIterator out)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  // TODO: nothing is done here to identify identical points

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using Face_loc = PMP::Face_location<TriangleMesh, typename K::FT>;
  using dscrptr_vrnt = PMP::descriptor_variant<TriangleMesh>;
  using Point_3 = typename K::Point_3;
  using EK = CGAL::Exact_predicates_exact_constructions_kernel;

  // 3 types of points: on edges, on faces, and existing vertices
  typename boost::property_map<TriangleMesh, CGAL::dynamic_face_property_t<std::size_t>>::type
    fid_map = get(CGAL::dynamic_face_property_t<std::size_t>(), tmesh, std::size_t(-1));
  typename boost::property_map<TriangleMesh, CGAL::dynamic_edge_property_t<std::size_t>>::type
    eid_map = get(CGAL::dynamic_edge_property_t<std::size_t>(), tmesh, std::size_t(-1));
  typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<std::size_t>>::type
    vid_map = get(CGAL::dynamic_vertex_property_t<std::size_t>(), tmesh, std::size_t(-1));

  std::vector< std::vector<std::pair<std::size_t, std::size_t>> > edges_per_face;
  std::vector< std::vector<std::size_t> > points_per_face;
  std::vector< std::vector<std::size_t> > points_per_edge;
  std::vector< edge_descriptor > edges_to_refine;
  std::vector< face_descriptor > faces_to_refine;
  std::vector< std::pair<std::size_t, vertex_descriptor> > input_vertices;

  std::vector<Face_loc> polyline_locations;

  auto register_point = [&](const dscrptr_vrnt& var, std::size_t id)
  {
    switch(var.index())
    {
      case 1:
      {
        halfedge_descriptor hd=std::get<halfedge_descriptor>(var);
        edge_descriptor ed = edge(hd, tmesh);
        if (get(eid_map, ed)==std::size_t(-1))
        {
          put(eid_map, ed, edges_to_refine.size());
          points_per_edge.emplace_back();
          edges_to_refine.push_back(ed);
        }
        points_per_edge[get(eid_map,ed)].push_back(id);
        return;
      }
      case 2:
      {
        face_descriptor fd=std::get<face_descriptor>(var);
        if (get(fid_map, fd)==std::size_t(-1))
        {
          put(fid_map, fd, faces_to_refine.size());
          points_per_face.emplace_back();
          edges_per_face.emplace_back();
          faces_to_refine.push_back(fd);
        }
        points_per_face[get(fid_map,fd)].push_back(id);
        return;
      }
      default:
        input_vertices.emplace_back(id, std::get<vertex_descriptor>(var));
    }
  };

  auto register_segment = [&](const dscrptr_vrnt& v1, const dscrptr_vrnt& v2,
                              const Face_loc& loc1, const Face_loc& loc2,
                              std::size_t i1, std::size_t i2)
  {
    if (v1.index()==0 && v2.index()==0) return;
    if (v1.index()==1 && v2.index()==1 && v1==v2) return;
    if (v1.index()==0 && v2.index()==1)
    {
      vertex_descriptor vd = std::get<vertex_descriptor>(v1);
      halfedge_descriptor hd = std::get<halfedge_descriptor>(v2);
      if (vd==source(hd, tmesh) || vd==target(hd, tmesh))
        return;
    }
    if (v2.index()==0 && v1.index()==1)
    {
      vertex_descriptor vd = std::get<vertex_descriptor>(v2);
      halfedge_descriptor hd = std::get<halfedge_descriptor>(v1);
      if (vd==source(hd, tmesh) || vd==target(hd, tmesh))
        return;
    }
    Face_loc copy1=loc1, copy2=loc2;
    CGAL_assertion_code(bool OK=)
    PMP::locate_in_common_face(copy1, copy2, tmesh);
    CGAL_assertion(OK);
    if (get(fid_map,copy1.first)==std::size_t(-1))
    {
      put(fid_map, copy1.first, faces_to_refine.size());
      points_per_face.emplace_back();
      edges_per_face.emplace_back();
      faces_to_refine.push_back(copy1.first);
    }
    edges_per_face[get(fid_map, copy1.first)].emplace_back(i1,i2);
  };

  for (const std::vector<Face_loc>& path : paths)
  {
    bool closed = PMP::are_locations_identical(path.front(), path.back(), tmesh);
    std::size_t nbp = path.size();

    if (closed)  nbp-=1;

    std::size_t prev_id=polyline_locations.size(), first_id=prev_id;
    dscrptr_vrnt prev_var = PMP::get_descriptor_from_location(path.front(), tmesh),
                 first_var = prev_var;

    polyline_locations.push_back( path.front() );
    register_point(prev_var, prev_id);

    for (std::size_t i=1; i<nbp; ++i)
    {
      std::size_t id=polyline_locations.size();
      polyline_locations.push_back(path[i]);
      dscrptr_vrnt var = PMP::get_descriptor_from_location(path[i], tmesh);
      register_point(var, id);
      register_segment(prev_var, var,
                       path[i-1], path[i],
                       prev_id, id);
      prev_id=id;
      prev_var=var;
    }

    if (closed)
      register_segment(prev_var, first_var,
                       path[path.size()-2], path.front(),
                       prev_id, first_id);
  }

  // precompute normals per face
  std::vector<typename K::Vector_3> face_normals(faces_to_refine.size());
  std::vector<std::array<vertex_descriptor, 3>> face_input_vertices(faces_to_refine.size());
  for (std::size_t fid=0; fid<faces_to_refine.size(); ++fid)
  {
    face_descriptor fd = faces_to_refine[fid];
    halfedge_descriptor hd = halfedge(fd, tmesh);
    face_normals[fid]=PMP::compute_face_normal(fd, tmesh); // TODO: vpm + compute using EK
    face_input_vertices[fid] = make_array(target(hd, tmesh),
                                          target(next(hd, tmesh), tmesh),
                                          target(prev(hd, tmesh), tmesh));
  }

  //vertices
  std::vector<EK::Point_3> exact_points(polyline_locations.size());
  std::vector<vertex_descriptor> polyline_vertices(
    polyline_locations.size(), boost::graph_traits<TriangleMesh>::null_vertex());

  // collect split point on edges
  std::vector<halfedge_descriptor> hedges_to_refine(edges_to_refine.size());
  std::vector< std::vector<std::pair<Point_3, std::size_t> > > points_on_hedges(edges_to_refine.size());
  for (std::size_t eid=0; eid<edges_to_refine.size(); ++eid)
  {
    edge_descriptor e = edges_to_refine[eid];
    halfedge_descriptor href = halfedge(e, tmesh);
    hedges_to_refine[eid]=href;

    std::vector< std::pair<typename K::FT, typename K::FT> > coordinates;
    coordinates.reserve(points_per_edge[eid].size());
    for (std::size_t pid : points_per_edge[eid])
    {
      halfedge_descriptor hd = halfedge(polyline_locations[pid].first, tmesh);
      int src_id=0;
      if (edge(hd, tmesh)!=e)
      {
        hd=next(hd, tmesh); ++src_id;
        if (edge(hd, tmesh)!=e)
        {
          hd=next(hd, tmesh);
          ++src_id;
        }
      }
      CGAL_assertion(edge(hd, tmesh)==e);
      coordinates.emplace_back(polyline_locations[pid].second[src_id], polyline_locations[pid].second[(src_id+1)%3]);
      if (hd!=href)
      {
        CGAL_assertion( opposite(hd, tmesh)==href );
        std::swap( coordinates.back().first, coordinates.back().second );
      }
    }

    // now sort coordinates
    std::vector<std::size_t> ids(coordinates.size());
    std::iota(ids.begin(), ids.end(), 0);

    std::sort(ids.begin(), ids.end(), [&coordinates](std::size_t i, std::size_t j)
                                      { return coordinates[i].first < coordinates[j].first; });

    points_on_hedges[eid].reserve(ids.size());
    for (std::size_t id : ids)
    {
      points_on_hedges[eid].emplace_back(
        PMP::construct_point(polyline_locations[points_per_edge[eid][id]], tmesh),
        points_per_edge[eid][id]);

      exact_points[points_per_edge[eid][id]] =
        PMP::construct_point(polyline_locations[points_per_edge[eid][id]], tmesh,
                        parameters::vertex_point_map(make_cartesian_converter_property_map<EK::Point_3>(tmesh.points())));

    }
  }

  CGAL::Cartesian_converter<K, EK> to_exact;

  // add new vertices per face
  for (std::size_t fid=0; fid<faces_to_refine.size(); ++fid)
  {
    typename K::Vector_3 face_normal=get(fnm, faces_to_refine[fid]);
    for(std::size_t vid : points_per_face[fid])
    {
      vertex_descriptor vd = add_vertex(tmesh);
      put(vnm, vd, face_normal);
      put(vid_map, vd, vid);
      tmesh.point(vd)=PMP::construct_point(polyline_locations[vid], tmesh);

      CGAL_assertion( polyline_vertices[vid]==boost::graph_traits<TriangleMesh>::null_vertex() );
      polyline_vertices[vid] = vd;
      exact_points[vid] = PMP::construct_point(polyline_locations[vid], tmesh,
                            parameters::vertex_point_map(make_cartesian_converter_property_map<EK::Point_3>(tmesh.points())));
    }
  }

  // set existing vertex
  for (const std::pair<std::size_t, vertex_descriptor>& id_and_vd : input_vertices)
  {
    polyline_vertices[id_and_vd.first]=id_and_vd.second;
    put(vid_map, id_and_vd.second, id_and_vd.first);
    exact_points[id_and_vd.first]=to_exact(tmesh.point(id_and_vd.second));
  }

  //Now split edges
  for (std::size_t eid=0; eid<edges_to_refine.size(); ++eid)
  {
    halfedge_descriptor href = hedges_to_refine[eid];
    CGAL_assertion( !is_border_edge(href, tmesh) ); // TODO shall we throw if we reach boundary?
    typename K::Vector_3 edge_normal = 0.5 * (get(fnm, face(href,tmesh))+get(fnm, face(opposite(href, tmesh), tmesh)));

    for (const std::pair<Point_3, std::size_t>& pt_and_id : points_on_hedges[eid])
    {
      CGAL_assertion( polyline_vertices[pt_and_id.second]==boost::graph_traits<TriangleMesh>::null_vertex() );
      href = ::CGAL::Euler::split_edge(href, tmesh);
      polyline_vertices[pt_and_id.second]=target(href, tmesh);
      // TODO: use VPM
      tmesh.point(target(href,tmesh)) = std::move(pt_and_id.first);
      put(vid_map, target(href,tmesh), pt_and_id.second);
      put(vnm, target(href,tmesh), edge_normal);
    }
  }

  // triangulate faces
  using CDT_traits = Projection_traits_3<EK>;
  using Vb = Triangulation_vertex_base_with_info_2<vertex_descriptor,CDT_traits>;
  using Fb = Constrained_triangulation_face_base_2<CDT_traits>;
  using TDS_2 = Triangulation_data_structure_2<Vb,Fb>;
  using CDT = Constrained_Delaunay_triangulation_2<CDT_traits,TDS_2>;
  using CDT_Vertex_handle = typename CDT::Vertex_handle;

  std::size_t nb_verts = polyline_vertices.size();
  polyline_vertices.resize(polyline_vertices.size()+3); // for input triangle vertices
  exact_points.resize(exact_points.size()+3);

  std::vector<CDT_Vertex_handle> id_to_cdt_vhandles(nb_verts);

  for (std::size_t fid=0; fid<faces_to_refine.size(); ++fid)
  {
    CDT_traits traits(to_exact(face_normals[fid]));
    CDT cdt(traits);

    // turn around the face to init the convex hull
    std::array<vertex_descriptor,3> tri_verts =
      make_array(face_input_vertices[fid][0], face_input_vertices[fid][1], face_input_vertices[fid][2]);

    std::array<CDT_Vertex_handle, 3> vhandles;
    vhandles[0]=cdt.insert_outside_affine_hull(to_exact(tmesh.point(tri_verts[0])));
    vhandles[1]=cdt.insert_outside_affine_hull(to_exact(tmesh.point(tri_verts[1])));
    vhandles[2] = cdt.tds().insert_dim_up(cdt.infinite_vertex(), false);
    vhandles[2]->set_point(to_exact(tmesh.point(tri_verts[2])));

    // assign ids to input triangle vertices
    std::size_t offset=nb_verts;
    for (int i=0;i<3;++i)
    {
      vhandles[i]->info()=tri_verts[i];
      if (get(vid_map, tri_verts[i])>=nb_verts) // and not simply -1 as previous triangulation could have set another dummy value
      {
        polyline_vertices[offset]=tri_verts[i];
        put(vid_map, tri_verts[i], offset);
        exact_points[offset]=vhandles[i]->point();
        ++offset;
      }
      else
      {
        id_to_cdt_vhandles[get(vid_map, tri_verts[i])]=vhandles[i];
      }
    }


    face_descriptor fd = faces_to_refine[fid];
    std::array<halfedge_descriptor,3> face_sides;
    for (int i=0; i<3; ++i)
    {
      for (halfedge_descriptor h : halfedges_around_target(face_input_vertices[fid][i], tmesh))
      {
        if (face(h, tmesh) == fd)
        {
          face_sides[i]=h;
          break;
        }
      }
    }

    // insert points on edges
    CDT_Vertex_handle vinf = cdt.infinite_vertex();
    for (int k=0; k<3; ++k)
    {
      if (next(face_sides[k], tmesh) != face_sides[(k+1)%3])
      {
        CDT_Vertex_handle src = vhandles[k], tgt = vhandles[(k+1)%3];
        halfedge_descriptor hcurr=next(face_sides[k], tmesh);
        CGAL_assertion(src->info()==source(hcurr, tmesh));
        CGAL_assertion(tgt->info()==target(face_sides[(k+1)%3], tmesh));

        CDT_Vertex_handle prev = src;
        while(hcurr!=face_sides[(k+1)%3])
        {
          typename CDT::Face_handle fh;
          CGAL_assertion_code(bool ok =)
          cdt.is_face(prev, tgt, vinf, fh);
          CGAL_assertion(ok);
          CGAL_assertion( get(vid_map, target(hcurr,tmesh)) != std::size_t(-1) );
          prev=cdt.insert_in_edge(exact_points[get(vid_map, target(hcurr,tmesh))], fh, fh->index(vinf));
          id_to_cdt_vhandles[get(vid_map, target(hcurr,tmesh))]=prev;
          prev->info() = target(hcurr,tmesh);
          cdt.restore_Delaunay(prev); // TODO maybe not each time but one global?
          CGAL_assertion(cdt.is_valid());
          hcurr=next(hcurr, tmesh);
        }
      }
    }

    // insert points in the face (TODO: we probably don't need spatial sorting)
    for (std::size_t vid : points_per_face[fid])
    {
      id_to_cdt_vhandles[vid] = cdt.insert(exact_points[vid]);
      id_to_cdt_vhandles[vid]->info()=polyline_vertices[vid];
    }



    // insert constrained edges
    for (const std::pair<std::size_t, std::size_t>& pi : edges_per_face[fid])
    {
      cdt.insert_constraint(id_to_cdt_vhandles[pi.first],id_to_cdt_vhandles[pi.second]);
    }

    // register cdt edge -> halfedge
    halfedge_descriptor hd = halfedge(fd, tmesh);
    std::map<std::pair<std::size_t, std::size_t>, halfedge_descriptor> edge_to_hedge;
    // triangle boundary first
    for (halfedge_descriptor h : halfedges_around_face(hd, tmesh))
    {
      edge_to_hedge[std::make_pair(get(vid_map,source(h,tmesh)), get(vid_map,target(h,tmesh)))]=h;
    }
    //grab edges that are not on the convex hull (these have already been created)
    for (typename CDT::Finite_edges_iterator it=cdt.finite_edges_begin();
                                             it!=cdt.finite_edges_end(); ++it)
    {
      typename CDT::Vertex_handle cdt_v0=it->first->vertex( cdt.ccw(it->second) );
      typename CDT::Vertex_handle cdt_v1=it->first->vertex( cdt.cw(it->second) );

      // consider edges not on the convex hull (not on the boundary of the face)
      // and create the corresponding halfedges
      if ( !cdt.is_infinite(it->first->vertex(it->second)) &&
           !cdt.is_infinite(cdt.mirror_vertex(it->first,it->second)) )
      {
        edge_descriptor e=add_edge(tmesh);
        halfedge_descriptor h=halfedge(e,tmesh), h_opp=opposite(h,tmesh);
        std::size_t i0=get(vid_map,cdt_v0->info()), i1=get(vid_map,cdt_v1->info());
        vertex_descriptor v0=polyline_vertices[i0], v1=polyline_vertices[i1];

        set_target(h,v0,tmesh);
        set_target(h_opp,v1,tmesh);
        set_halfedge(v0,h,tmesh);
        set_halfedge(v1,h_opp,tmesh);

        edge_to_hedge[std::make_pair(i0,i1)]=h_opp;
        edge_to_hedge[std::make_pair(i1,i0)]=h;

        if (it->first->is_constrained(it->second))
          *out++=h;
      }
    }

    //grab triangles.
    face_descriptor current_face = fd;
    typename K::Vector_3 face_normal = get(fnm, fd);
    for (typename CDT::Finite_faces_iterator it=cdt.finite_faces_begin(),
                                             it_end=cdt.finite_faces_end();;)
    {
      typename CDT::Vertex_handle cdt_v0=it->vertex(0);
      typename CDT::Vertex_handle cdt_v1=it->vertex(1);
      typename CDT::Vertex_handle cdt_v2=it->vertex(2);

      std::size_t i0=get(vid_map,cdt_v0->info()), i1=get(vid_map,cdt_v1->info()), i2=get(vid_map,cdt_v2->info());

      CGAL_assertion(edge_to_hedge.count(std::make_pair(i0,i1))!= 0);
      CGAL_assertion(edge_to_hedge.count(std::make_pair(i1,i2))!= 0);
      CGAL_assertion(edge_to_hedge.count(std::make_pair(i2,i0))!= 0);

      halfedge_descriptor h01=edge_to_hedge[std::make_pair(i0,i1)];
      halfedge_descriptor h12=edge_to_hedge[std::make_pair(i1,i2)];
      halfedge_descriptor h20=edge_to_hedge[std::make_pair(i2,i0)];

      set_next(h01,h12,tmesh);
      set_next(h12,h20,tmesh);
      set_next(h20,h01,tmesh);

      //update face halfedge
      set_halfedge(current_face,h01,tmesh);

      //update face of halfedges
      set_face(h01,current_face,tmesh);
      set_face(h12,current_face,tmesh);
      set_face(h20,current_face,tmesh);

      if ( ++it!=it_end )
      {
        current_face=add_face(tmesh);
        put(fnm, current_face, face_normal);
      }
      else
        break;
    }
  }

  return out;
}

} } // end of CGAL::Vector_graphics_on_surfaces namespace

#endif // CGAL_POLYGON_MESH_PROCESSING_BSURF_REFINE_MESH_ALONG_PATHS_H
