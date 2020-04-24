// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj

#ifndef CGAL_TETRAHEDRAL_REMESHING_FMLS_H
#define CGAL_TETRAHEDRAL_REMESHING_FMLS_H

#include <CGAL/license/Tetrahedral_remeshing.h>

// -------------------------------------------
// FMLS
// A Fast Moving Least Square operator for 3D
// points sets.
// -------------------------------------------

#include <CGAL/number_utils.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>

#include <cmath>
#include <vector>
#include <unordered_set>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <boost/functional/hash.hpp>


namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
// --------------------------------------------------------------
//  MLS Projection
// --------------------------------------------------------------

inline float wendland(float x, float h)
{
  x = CGAL::abs(x);
  if (x < h)
    return CGAL::square(CGAL::square(1 - x / h)) * (4 * x / h + 1);
  else
    return 0.0;
}

inline void setPNSample(std::vector<float>& p,
                        const std::size_t& i,
                        const float x, const float y, const float z,
                        const float nx, const float ny, const float nz)
{
  p[6 * i] = x;
  p[6 * i + 1] = y;
  p[6 * i + 2] = z;
  p[6 * i + 3] = nx;
  p[6 * i + 4] = ny;
  p[6 * i + 5] = nz;
}

template<typename Gt>
inline CGAL::Vector_3<Gt> projectOn(const CGAL::Vector_3<Gt>& x,
                                    const CGAL::Vector_3<Gt>& N,
                                    const CGAL::Vector_3<Gt>& P)
{
  typename Gt::Compute_scalar_product_3 scalar_product
    = Gt().compute_scalar_product_3_object();

  typename Gt::FT w = scalar_product((x - P), N);
  return x - (N * w);
}

template<typename Gt>
inline typename Gt::FT length(const CGAL::Vector_3<Gt>& v)
{
  return CGAL::approximate_sqrt(v.squared_length());
}

template<typename Gt>
inline typename Gt::FT distance(const CGAL::Vector_3<Gt>& v1,
                                const CGAL::Vector_3<Gt>& v2)
{
  CGAL::Vector_3<Gt> diff(v1.x() - v2.x(),
                          v1.y() - v2.y(),
                          v1.z() - v2.z());
  return length(diff);
}

template<typename Gt>
inline void weightedPointCombination(const CGAL::Vector_3<Gt>& x,
                                     const CGAL::Vector_3<Gt>& pi,
                                     const CGAL::Vector_3<Gt>& ni,
                                     float sigma_s, bool bilateral, float sigma_r,
                                     bool hermite,
                                     CGAL::Vector_3<Gt>& c, CGAL::Vector_3<Gt>& nc, float& sumW)
{
  float w = wendland(distance(x, pi), sigma_s);
  if (bilateral)
    w *= wendland(length(x - projectOn(x, ni, pi)), sigma_r);
  if (hermite)
    c += w * projectOn(x, ni, pi);
  else
    c += w * pi;
  nc += w * ni;
  sumW += w;
}

template<typename Gt>
class FMLS
{
  typedef typename Gt::Vector_3 Vector_3;
  typedef typename Gt::FT       FT;

public:
  // FMLS provide MLS projection and filtering from a point set.
  // The underlying data structure is a simple list of float in the PN format
  // A PN object is a list of 6xfloat32 chunk :
  // x0,y0,z0,nx0,ny0,nz0,x1,y1,z1,nx1,ny1,nz1,...
  // with {xi,yi,zi} the position and {nxi,nyi,nzi} the normal vector of the
  // i-th point sample. A PN can be read from and write to file directly
  // (identity serialization) and is handled as a simple float32 pointer.
  //
  // Use the 'fast*' methods. Brute force methods are inserted only for comparison.
  //
  // Memory policy: a FMLS class manages itself al its belonging objects.
  // Therefore, PN and filtered PN lists are the property of this class, and should
  // be copied if modified outside the class.
  FMLS()
  {
    PNSize = 0;
    PNScale = 1.0f;
    MLSRadius = 0.1f;
    bilateralRange = 0.05f;
    bilateral = false;
    hermite = false;
    numIter = 1;
  }

  // --------------------------------------------------------------
  //  Main Interface
  // --------------------------------------------------------------

  std::vector<float> createPN(unsigned int size)
  {
    return std::vector<float>(size * SURFEL_SIZE);
  }

  void setPN(const std::vector<float>& newPN,
             const unsigned int newPNSize,
             const float pointSpacing)
  {
    freeCPUMemory();
    PN = newPN;
    PNSize = newPNSize;
    computePNScale();
    MLSRadius = 3 * pointSpacing / PNScale;
    grid.clear();
    grid.init(PN, PNSize, MLSRadius * PNScale);
  }

  // Compute, according to the current point sampling stored in FMLS, the MLS projection
  // of p and store the resulting position in q and normal in n.
  void fastProjectionCPU(const Vector_3& p, Vector_3& q, Vector_3& n) const
  {
    float sigma_s = PNScale * MLSRadius;
    float sigma_r = bilateralRange;

    Vector_3 g = (p - Vector_3(grid.getMinMax()[0], grid.getMinMax()[1], grid.getMinMax()[2])) / sigma_s;
    std::array<FT, 3> gxyz = { g.x(), g.y(), g.z() };

    for (unsigned int j = 0; j < 3; j++) {
      gxyz[j] = floor(gxyz[j]);
      if (gxyz[j] < 0.f)
        gxyz[j] = 0.f;
      if (gxyz[j] >= grid.getRes()[j])
        gxyz[j] = grid.getRes()[j] - 1;
    }
    unsigned int minIt[3], maxIt[3];
    for (unsigned int j = 0; j < 3; j++) {
      if (((unsigned int)gxyz[j]) == 0)
        minIt[j] = 0;
      else
        minIt[j] = ((unsigned int)gxyz[j]) - 1;
      if (((unsigned int)gxyz[j]) == (grid.getRes()[j] - 1))
        maxIt[j] = (grid.getRes()[j] - 1);
      else
        maxIt[j] = ((unsigned int)gxyz[j]) + 1;
    }
    Vector_3 c = CGAL::NULL_VECTOR;
    float sumW = 0.f;
    unsigned int it[3];
    for (it[0] = minIt[0]; it[0] <= maxIt[0]; it[0]++)
      for (it[1] = minIt[1]; it[1] <= maxIt[1]; it[1]++)
        for (it[2] = minIt[2]; it[2] <= maxIt[2]; it[2]++) {
          unsigned int gridIndex = grid.getLUTElement(it[0], it[1], it[2]);
          if (gridIndex == 2 * PNSize)
            continue;
          unsigned int neigh = grid.getCellIndicesSize(it[0], it[1], it[2]);
          for (unsigned int j = 0; j < neigh; j++) {
            unsigned int k = grid.getIndicesElement(it[0], it[1], it[2], j);
            Vector_3 pk(PN[6 * k], PN[6 * k + 1], PN[6 * k + 2]);
            Vector_3 nk(PN[6 * k + 3], PN[6 * k + 4], PN[6 * k + 5]);
            weightedPointCombination(p, pk, nk, sigma_s, bilateral, sigma_r, hermite, c, n, sumW);
          }
        }
    if (sumW == 0.f) {
      n = Vector_3(1.f, 0.f, 0.f);
      q = p;
    }
    else {
      c /= sumW;
      normalize(n, Gt());
      q = projectOn(p, n, c);
    }
  }

  // Compute the MLS projection of the list of point stored in pv and store the resulting
  // positions and normal in qv. qv must be preallocated to stroe 6*pvSize float32.
  // The strid indicates the offsets in qv (the defautl value of 3 means that the qv
  // is compact: pv={x0,y0,z0,x1,y1,z1...}. If pv contains also normals for instance,
  // the stride should be set to 6.
  void fastProjectionCPU(const std::vector<float>& pv,
                         const std::size_t pvSize,
                         std::vector<float>& qv,
                         std::size_t stride = 3) const
  {
    for (std::size_t i = 0; i < pvSize; i++)
    {
      Vector_3 p(pv[stride * i], pv[stride * i + 1], pv[stride * i + 2]);
      Vector_3 q, n;
      for (unsigned int j = 0; j < numIter; j++)
      {
        q = CGAL::NULL_VECTOR;
        n = CGAL::NULL_VECTOR;
        fastProjectionCPU(p, q, n);
        p = q;
      }
      setPNSample(qv, i, q[0], q[1], q[2], n[0], n[1], n[2]);
    }
  }

  // Brute force version. O(PNSize) complexity. For comparison only.
  void projectionCPU(const Vector_3& x, Vector_3& q, Vector_3& n)
  {
    float sigma_s = MLSRadius * PNScale;
    float sigma_r = bilateralRange;
    Vector_3 p(x);
    for (unsigned int k = 0; k < numIter; k++) {
      Vector_3 c = CGAL::NULL_VECTOR;
      n = CGAL::NULL_VECTOR;
      float sumW = 0.f;
      for (unsigned int j = 0; j < PNSize; j++) {
        Vector_3 pj(PN[6 * j], PN[6 * j + 1], PN[6 * j + 2]);
        Vector_3 nj(PN[6 * j + 3], PN[6 * j + 4], PN[6 * j + 5]);
        weightedPointCombination(p, pj, nj, sigma_s, bilateral, sigma_r, hermite, c, n, sumW);
      }
      c /= sumW;
      n.normalize();
      q = projectOn(p, n, c);
      p = q;
    }

  }
  // Brute force version. O(pvSize*PNSize) complexity. For comparison only.
  void projectionCPU(const std::vector<float>& pv,
                     unsigned int pvSize,
                     std::vector<float>& qv,
                     unsigned int stride = 3)
  {
    for (int i = 0; i < int(pvSize); i++) {
      Vector_3 p(pv[stride * i], pv[stride * i + 1], pv[stride * i + 2]);
      Vector_3 q, n;
      for (unsigned int j = 0; j < numIter; j++) {
        q = CGAL::NULL_VECTOR;
        n = CGAL::NULL_VECTOR;
        projectionCPU(p, q, n);
        p = q;
      }
      setPNSample(qv, i, q[0], q[1], q[2], n[0], n[1], n[2]);
    }
  }

  // --------------------------------------------------------------
  //  Accessors
  // --------------------------------------------------------------

  // Number of elements of the PN. One elemnt is a 6-float32 chunk.
  inline unsigned int getPNSize() const { return PNSize; }
  inline std::vector<float>& getPN() { return PN; }
  inline const std::vector<float>& getPN() const { return PN; }

  // Min/Max corners of PN's bounding volume
  inline const float* getMinMax() const { return grid.getMinMax(); }
  // Radius of the bounding sphere of the PN
  inline float getPNScale() const { return PNScale; }
  // Normalized MLS support size
  inline float getMLSRadius() const { return MLSRadius; }
  inline void setMLSRadius(float s) { MLSRadius = s; grid.clear(); grid.init(PN, PNSize, MLSRadius * PNScale); }
  // Bilateral weighting for feature preservation (inspired by [Jones 2003]).
  inline bool isBilateral() const { return bilateral; }
  inline void toggleBilateral(bool b) { bilateral = b; }
  // Bilateral support size for the range weight
  inline float getBilateralRange() const { return bilateralRange; }
  inline void setBilateralRange(float r) { bilateralRange = r; }
  // Hermite interpolation [Alexa 2009]
  inline bool isHermite() const { return hermite; }
  inline void toggleHermite(bool b) { hermite = b; }
  // Fix number of iterations of the MLS projection
  inline unsigned int getNumOfIter() const { return numIter; }
  inline void setNumOfIter(unsigned int i) { numIter = i; }

  // --------------------------------------------------------------
  //  Misc.
  // --------------------------------------------------------------

  // Size of a point sample in bytes (6xfloat32: 3 for position and normal
  static const unsigned int SURFEL_SIZE = 24;

private:

  void computePNScale()
  {
    Vector_3 c = CGAL::NULL_VECTOR;
    for (std::size_t i = 0; i < PNSize; i++)
      c += Vector_3(PN[6 * i], PN[6 * i + 1], PN[6 * i + 2]);
    c /= PNSize;
    PNScale = 0.f;
    for (std::size_t i = 0; i < PNSize; i++) {
      float r = distance(c, Vector_3(PN[6 * i], PN[6 * i + 1], PN[6 * i + 2]));
      if (r > PNScale)
        PNScale = r;
    }
  }

  // --------------------------------------------------------------
  //  3D Grid Structure
  // --------------------------------------------------------------
  // --------------------------------------------------------------
  //  Grid data structure for fast r-ball neighborhood query
  // --------------------------------------------------------------

  class Grid
  {
  public:
    Grid()
      : cellSize(1.f)
    {}
    ~Grid()
    {
      clear();
    }

    void init(const std::vector<float>& PN, unsigned int PNSize, float sigma_s)
    {
      cellSize = sigma_s;
      for (unsigned int i = 0; i < 3; i++) {
        minMax[i] = PN[i];
        minMax[3 + i] = PN[i];
      }
      for (unsigned int i = 0; i < PNSize; i++)
        for (unsigned int j = 0; j < 3; j++) {
          if (PN[6 * i + j] < minMax[j])
            minMax[j] = PN[6 * i + j];
          if (PN[6 * i + j] > minMax[3 + j])
            minMax[3 + j] = PN[6 * i + j];
        }
      for (unsigned int i = 0; i < 3; i++) {
        minMax[i] -= 0.001f;
        minMax[3 + i] += 0.001f;
      }
      for (unsigned int i = 0; i < 3; i++)
        res[i] = (unsigned int)ceil((minMax[3 + i] - minMax[i]) / cellSize);
      unsigned int LUTSize = res[0] * res[1] * res[2];
      LUT.resize(LUTSize);
      LUT.assign(LUTSize, 0);

      unsigned int nonEmptyCells = 0;
      Vector_3 gMin(minMax[0], minMax[1], minMax[2]);
      Vector_3 gMax(minMax[3], minMax[4], minMax[5]);
      for (unsigned int i = 0; i < PNSize; i++) {
        unsigned int index = getLUTIndex(Vector_3(PN[6 * i], PN[6 * i + 1], PN[6 * i + 2]));
        if (LUT[index] == 0)
          nonEmptyCells++;
        LUT[index]++;
      }
      unsigned int indicesSize = PNSize + nonEmptyCells;
      indices.reserve(indicesSize);
      indices.assign(indicesSize, 0);

      unsigned int cpt = 0;
      for (unsigned int i = 0; i < res[0]; i++)
        for (unsigned int j = 0; j < res[1]; j++)
          for (unsigned int k = 0; k < res[2]; k++) {
            unsigned int index = getLUTIndex(i, j, k);
            if (LUT[index] != 0) {
              indices[cpt] = LUT[index];
              LUT[index] = cpt;
              cpt += indices[cpt] + 1;
              indices[cpt - 1] = 0; // local iterator for subsequent filling
            }
            else
              LUT[index] = 2 * PNSize;
          }
      for (unsigned int i = 0; i < PNSize; i++) {
        Vector_3 p = Vector_3(PN[6 * i], PN[6 * i + 1], PN[6 * i + 2]);
        unsigned int indicesIndex = getLUTElement(p);
        unsigned int totalCount = indices[indicesIndex];
        unsigned int countIndex = indicesIndex + totalCount;
        unsigned int currentCount = indices[countIndex];
        if (currentCount < indices[indicesIndex])
          indices[countIndex]++;
        unsigned int pIndex = indicesIndex + 1 + currentCount;
        indices[pIndex] = i;
      }
    }

    void clear()
    {
      cellSize = 1.f;
    }

    // Accessors

    inline const std::array<float, 6> getMinMax() const { return minMax; }
    inline const std::array<unsigned int, 3> getRes() const { return res; }
    inline float getCellSize() const { return cellSize; }
    inline std::vector<unsigned int>& getLUT() { return LUT; }
    inline const std::vector<unsigned int>& getLUT() const { return LUT; }
    inline unsigned int getLUTIndex(unsigned int i,
                                    unsigned int j,
                                    unsigned int k) const
    {
      return k * res[0] * res[1] + j * res[0] + i;
    }
    inline unsigned int getLUTElement(unsigned int i,
                                      unsigned int j,
                                      unsigned int k) const
    {
      return LUT[getLUTIndex(i, j, k)];
    }
    unsigned int getLUTIndex(const Vector_3& x) const
    {
      Vector_3 vp = (x - Vector_3(minMax[0], minMax[1], minMax[2])) / cellSize;
      std::array<FT, 3> p = { vp.x(), vp.y(), vp.z() };
      for (unsigned int j = 0; j < 3; j++) {
        p[j] = floor(p[j]);
        if (p[j] < 0)
          p[j] = 0.f;
        if (p[j] >= res[j])
          p[j] = res[j] - 1;
      }
      unsigned index = ((unsigned int)floor(p[2])) * res[0] * res[1]
                       + ((unsigned int)floor(p[1])) * res[0]
                       + ((unsigned int)floor(p[0]));
      return index;
    }
    inline unsigned int getLUTElement(const Vector_3& x) const {
      return LUT[getLUTIndex(x)];
    }
    inline std::vector<unsigned int>& getIndices() { return indices; }
    inline const std::vector<unsigned int>& getIndices() const { return indices; }
    inline unsigned int getCellIndicesSize(unsigned int i,
                                           unsigned int j,
                                           unsigned int k) const {
      return indices[getLUTElement(i, j, k)];
    }
    inline unsigned int getIndicesElement(unsigned int i,
                                          unsigned int j,
                                          unsigned int k,
                                          unsigned int e) const {
      return indices[getLUTElement(i, j, k) + 1 + e];
    }

  private:
    std::array<float, 6> minMax;
    float cellSize;
    std::array<unsigned int, 3> res;
    std::vector<unsigned int> LUT; // 3D Index Look-Up Table
    std::vector<unsigned int> indices; // 3D Grid data
  };


  // --------------------------------------------------------------
  //  Memory Managment
  // --------------------------------------------------------------

  void freeCPUMemory()
  {
    PN.clear();
    PNSize = 0;
  }

  // --------------------------------------------------------------
  //  CPU Data
  // --------------------------------------------------------------

  std::vector<float> PN;
  unsigned int PNSize;
  float PNScale; // size of the bounding sphere radius
  float MLSRadius;
  float bilateralRange;
  bool bilateral;
  bool hermite;
  unsigned int numIter;
  Grid grid;
};

template<typename Subdomain__FMLS,
         typename Subdomain__FMLS_indices,
         typename VerticesNormalsMap,
         typename VerticesSurfaceIndices,
         typename C3t3>
void createMLSSurfaces(Subdomain__FMLS& subdomain_FMLS,
                       Subdomain__FMLS_indices& subdomain_FMLS_indices,
                       const VerticesNormalsMap& vertices_normals,
                       const VerticesSurfaceIndices& vertices_surface_indices,
                       const C3t3& c3t3,
                       const int upsample = 2) // can be 0, 1 or 2
{
  typedef typename C3t3::Surface_patch_index Surface_index;
  typedef typename C3t3::Triangulation       Tr;
  typedef typename Tr::Edge                  Edge;
  typedef typename Tr::Vertex_handle         Vertex_handle;
  typedef typename Tr::Geom_traits           Gt;
  typedef typename Gt::Point_3               Point_3;
  typedef typename Gt::Vector_3              Vector_3;

  typedef typename VerticesSurfaceIndices::mapped_type    VertexSurfaces;
  typedef typename VerticesSurfaceIndices::const_iterator VerticesSurfaceIterator;

  const Tr& tr = c3t3.triangulation();

  //createAreaWeightedUpSampledMLSSurfaces(0);
  //return ;
  subdomain_FMLS.clear();
  subdomain_FMLS_indices.clear();

  typedef boost::unordered_map<Surface_index, unsigned int> SurfaceIndexMap;

  SurfaceIndexMap current_subdomain_FMLS_indices;
  SurfaceIndexMap subdomain_sample_numbers;

  //Count the number of vertices for each boundary surface (i.e. one per label)
  for (const Vertex_handle vit : tr.finite_vertex_handles())
  {
    VerticesSurfaceIterator sit = vertices_surface_indices.find(vit);
    if (sit == vertices_surface_indices.end())
      continue;

    const VertexSurfaces& v_surface_indices = vertices_surface_indices.at(vit);
    CGAL_assertion(vit->in_dimension() <= 2);

    for(const Surface_index& si : v_surface_indices)
    {
      subdomain_sample_numbers[si]++;
    }
  }

  if (upsample > 0)
  {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Up sampling MLS " << upsample << std::endl;
#endif
    for (typename C3t3::Facets_in_complex_iterator
      fit = c3t3.facets_in_complex_begin();
      fit != c3t3.facets_in_complex_end(); ++fit)
    {
      const Surface_index surf_i = c3t3.surface_patch_index(*fit);
      if (upsample == 1)
        subdomain_sample_numbers[surf_i] ++;
      else if (upsample == 2)
        subdomain_sample_numbers[surf_i] += 4;
    }
  }

  std::vector< std::vector<float> > pns;

  int count = 0;
  //Memory allocation for the point plus normals of the point samples
  for (typename SurfaceIndexMap::iterator it = subdomain_sample_numbers.begin();
       it != subdomain_sample_numbers.end(); ++it)
  {
    current_subdomain_FMLS_indices[it->first] = count;
    pns.push_back(std::vector<float>(it->second * 6, 0.f));
    count++;
  }

  std::vector<int> current_v_count(count, 0);
  std::vector<float> point_spacing(count, 0.f);
  std::vector<int> point_spacing_count(count, 0);

  //Allocation of the PN
  for (Vertex_handle vit : tr.finite_vertex_handles())
  {
    VerticesSurfaceIterator sit = vertices_surface_indices.find(vit);
    if (sit == vertices_surface_indices.end())
      continue;

    const VertexSurfaces& v_surface_indices = vertices_surface_indices.at(vit);
    CGAL_assertion(vit->in_dimension() <= 2);

    for (const Surface_index& surf_i : v_surface_indices)
    {
      const int fmls_id = current_subdomain_FMLS_indices[surf_i];

      const Point_3& p = point(vit->point());

      pns[fmls_id][6 * current_v_count[fmls_id]] = p.x();
      pns[fmls_id][6 * current_v_count[fmls_id] + 1] = p.y();
      pns[fmls_id][6 * current_v_count[fmls_id] + 2] = p.z();

      const Vector_3& normal = vertices_normals.at(vit).at(surf_i);

      pns[fmls_id][6 * current_v_count[fmls_id] + 3] = normal.x();
      pns[fmls_id][6 * current_v_count[fmls_id] + 4] = normal.y();
      pns[fmls_id][6 * current_v_count[fmls_id] + 5] = normal.z();

      current_v_count[fmls_id]++;
    }
  }

  typedef std::pair<Vertex_handle, Vertex_handle> Edge_vv;
  if (upsample == 0)
  {
    std::unordered_set<Edge_vv, boost::hash<Edge_vv> > edgeMap;

    for (typename C3t3::Facets_in_complex_iterator
         fit = c3t3.facets_in_complex_begin();
         fit != c3t3.facets_in_complex_end(); ++fit)
    {
      for (int i = 0; i < 2; i++)
      {
        for (int j = i + 1; j < 3; j++)
        {
          Edge edge(fit->first, indices(fit->second,i), indices(fit->second,j));

          Vertex_handle vh0 = edge.first->vertex(edge.second);
          Vertex_handle vh1 = edge.first->vertex(edge.third);
          Edge_vv e = make_vertex_pair(vh0, vh1);
          if ( vertices_surface_indices.find(vh0) != vertices_surface_indices.end()
               && vertices_surface_indices.find(vh1) != vertices_surface_indices.end()
               && edgeMap.find(e) == edgeMap.end())
          {
            edgeMap.insert(e);

            const Surface_index surf_i = c3t3.surface_patch_index(*fit);
            const int fmls_id = current_subdomain_FMLS_indices[surf_i];

            point_spacing[fmls_id] += CGAL::approximate_sqrt(
                                        CGAL::squared_distance(point(vh0->point()), point(vh1->point())));
            point_spacing_count[fmls_id] ++;
          }
        }
      }
    }
  }

  if (upsample > 0)
  {
    for (typename C3t3::Facets_in_complex_iterator
         fit = c3t3.facets_in_complex_begin();
         fit != c3t3.facets_in_complex_end(); ++fit)
    {
      const Surface_index surf_i = c3t3.surface_patch_index(*fit);

      const int fmls_id = current_subdomain_FMLS_indices[surf_i];

      Vertex_handle vhs[3] = { fit->first->vertex(indices(fit->second, 0)),
                               fit->first->vertex(indices(fit->second, 1)),
                               fit->first->vertex(indices(fit->second, 2)) };
      Vector_3 points[3] = { Vector_3(CGAL::ORIGIN, point(vhs[0]->point())),
                             Vector_3(CGAL::ORIGIN, point(vhs[1]->point())),
                             Vector_3(CGAL::ORIGIN, point(vhs[2]->point())) };
      Vector_3 normals[3] = { vertices_normals.at(vhs[0]).at(surf_i),
                              vertices_normals.at(vhs[1]).at(surf_i),
                              vertices_normals.at(vhs[2]).at(surf_i) };

      std::vector<Vector_3> points_to_add;
      std::vector<Vector_3> n_points_to_add;

      //Add the barycenter of the facet
      Vector_3 barycenter = (points[0] + points[1] + points[2]) / 3.;
      Vector_3 n_barycenter = (normals[0] + normals[1] + normals[2]);

      n_barycenter = n_barycenter / CGAL::sqrt((n_barycenter * n_barycenter));

      points_to_add.push_back(barycenter);
      n_points_to_add.push_back(n_barycenter);

      if (upsample == 1)
      {
        for (int i = 0; i < 3; i++)
        {
          Vector_3 space_1 = barycenter - points[i];

          point_spacing[fmls_id] += CGAL::sqrt(space_1 * space_1);
          point_spacing_count[fmls_id] ++;
        }
      }
      else if (upsample == 2)
      {
        for (int i = 0; i < 3; i++)
        {
          int i1 = (i + 1) % 3;
          int i2 = (i + 2) % 3;

          Vector_3 p = (barycenter + points[i1] + points[i2]) / 3.;
          Vector_3 n = (n_barycenter + normals[i1] + normals[i2]);

          n = n / CGAL::sqrt(n * n);

          points_to_add.push_back(p);
          n_points_to_add.push_back(n);

          Vector_3 space_1 = p - barycenter;
          Vector_3 space_2 = p - points[i1];
          Vector_3 space_3 = p - points[i2];

          point_spacing[fmls_id] += CGAL::sqrt(space_1 * space_1);
          point_spacing[fmls_id] += CGAL::sqrt(space_2 * space_2);
          point_spacing[fmls_id] += CGAL::sqrt(space_3 * space_3);

          point_spacing_count[fmls_id] += 3;
        }
      }
      for (unsigned int i = 0; i < points_to_add.size(); i++)
      {
        Vector_3& point = points_to_add[i];

        pns[fmls_id][6 * current_v_count[fmls_id]] = point.x();
        pns[fmls_id][6 * current_v_count[fmls_id] + 1] = point.y();
        pns[fmls_id][6 * current_v_count[fmls_id] + 2] = point.z();

        Vector_3& normal = n_points_to_add[i];

        pns[fmls_id][6 * current_v_count[fmls_id] + 3] = normal.x();
        pns[fmls_id][6 * current_v_count[fmls_id] + 4] = normal.y();
        pns[fmls_id][6 * current_v_count[fmls_id] + 5] = normal.z();

        current_v_count[fmls_id]++;
      }
    }
  }


  int nb_of_mls_to_create = 0;
  double average_point_spacing = 0;

  //Cretaing the actual MLS surfaces
  for (typename SurfaceIndexMap::iterator it = current_subdomain_FMLS_indices.begin();
       it != current_subdomain_FMLS_indices.end(); ++it)
  {
    if (current_v_count[it->second] > 3)
    {
      nb_of_mls_to_create++;

      float current_point_spacing = point_spacing[it->second] / point_spacing_count[it->second];
      point_spacing[it->second] = current_point_spacing;

      average_point_spacing += current_point_spacing;
    }
  }

  average_point_spacing = average_point_spacing / nb_of_mls_to_create;

  subdomain_FMLS.resize(nb_of_mls_to_create, FMLS<Gt>());

  count = 0;
  //Creating the actual MLS surfaces
  for (typename SurfaceIndexMap::iterator it = current_subdomain_FMLS_indices.begin();
       it != current_subdomain_FMLS_indices.end(); ++it)
  {
    if (current_v_count[it->second] > 3)
    {
      float current_point_spacing = point_spacing[it->second];

      //subdomain_FMLS[count].toggleHermite(true);
      subdomain_FMLS[count].setPN(pns[it->second],
                                  current_v_count[it->second],
                                  current_point_spacing);
      //  subdomain_FMLS[count].toggleHermite(true);
      subdomain_FMLS_indices[it->first] = count;

      count++;
    }
    else {
      std::cout << "Problem of number for MLS : " << current_v_count[it->second] << std::endl;
    }
  }
}

} // internal
} // Tetrahedral_remeshing
} // CGAL

#endif //CGAL_TETRAHEDRAL_REMESHING_FMLS_H
