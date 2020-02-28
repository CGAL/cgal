#ifndef CGAL_TETRAHEDRAL_REMESHING_FMLS_H
#define CGAL_TETRAHEDRAL_REMESHING_FMLS_H

// -------------------------------------------
// FMLS
// A Fast Moving Least Square operator for 3D
// points sets.
// 
// Copyright (C) 2006-2011 Tamy Boubekeur
// All rights reserved.
// -------------------------------------------

// -------------------------------------------
// Disclaimer: this code is dirty in the
// meaning that there is no attention paid to
// proper class attribute access, memory
// management or optimisation of any kind. It
// is designed for quick-and-dirty testing
// purpose.
// -------------------------------------------

#include <CGAL/number_utils.h>
#include <CGAL/Point_3.h>

#include <cmath>
#include <string>
#include <vector>

#include <boost/unordered_map.hpp>

#include <CGAL/Tetrahedral_remeshing/internal/Tetrahedral_remeshing_helpers.h>

#include "Vec3D.h"


namespace CGAL
{
  namespace Tetrahedral_remeshing
  {
    namespace internal
    {
      // --------------------------------------------------------------
      //  CPU Memory Code
      // --------------------------------------------------------------

      template<typename T>
      void freeCPUResource(T** res) {
        if (*res != NULL) {
          free(*res);
          *res = NULL;
        }
      }

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

      inline void setPNSample(float* p, unsigned int i,
                              float x, float y, float z,
                              float nx, float ny, float nz)
      {
        p[6 * i] = x;
        p[6 * i + 1] = y;
        p[6 * i + 2] = z;
        p[6 * i + 3] = nx;
        p[6 * i + 4] = ny;
        p[6 * i + 5] = nz;
      }

      inline void weightedPointCombination(const Vec3Df& x, const Vec3Df& pi, const Vec3Df& ni,
                                           float sigma_s, bool bilateral, float sigma_r,
                                           bool hermite,
                                           Vec3Df& c, Vec3Df& nc, float& sumW)
      {
        float w = wendland(Vec3Df::distance(x, pi), sigma_s);
        if (bilateral)
          w *= wendland((x - x.projectOn(ni, pi)).getLength(), sigma_r);
        if (hermite)
          c += w * x.projectOn(ni, pi);
        else
          c += w * pi;
        nc += w * ni;
        sumW += w;
      }

      class FMLS
      {
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
          PN = NULL;
          PNSize = 0;
          PNScale = 1.0f;
          MLSRadius = 0.1f;
          bilateralRange = 0.05f;
          bilateral = false;
          hermite = false;
          numIter = 1;
        }
        ~FMLS()
        {
          freeCPUMemory();
        }

        // --------------------------------------------------------------
        //  Main Interface 
        // --------------------------------------------------------------

        float* createPN(unsigned int size)
        {
          return  (float*)malloc(size * SURFEL_SIZE);
        }

        float* clonePN()
        {
          float* clone = createPN(PNSize);
          memcpy(clone, PN, PNSize * SURFEL_SIZE);
          return clone;
        }

        void loadPN(const char* filename)
        {
          freeCPUMemory();
          FILE* file = fopen(filename, "r");
          if (!file)
            throw Exception("Cannot read file" + std::string(filename));
          fseek(file, 0, SEEK_END);
          unsigned int numOfByte = ftell(file);
          fseek(file, 0, SEEK_SET);
          PNSize = numOfByte / SURFEL_SIZE;
          PN = createPN(PNSize);
          fread(PN, SURFEL_SIZE, PNSize, file);
          fclose(file);

          computePNScale();
          grid.clear();
          grid.init(PN, PNSize, MLSRadius * PNScale);
        }

        void setPN(float* newPN, unsigned int newPNSize)
        {
          freeCPUMemory();
          PN = newPN;
          PNSize = newPNSize;
          computePNScale();
          grid.clear();
          grid.init(PN, PNSize, MLSRadius * PNScale);
        }

        void setPN(float* newPN, unsigned int newPNSize, float pointSpacing)
        {
          freeCPUMemory();
          PN = newPN;
          PNSize = newPNSize;
          computePNScale();
          MLSRadius = 3 * pointSpacing / PNScale;
          grid.clear();
          grid.init(PN, PNSize, MLSRadius * PNScale);
        }

        static void savePN(float* pn, unsigned int size, const char* filename)
        {
          FILE* file = fopen(filename, "w");
          if (!file)
            throw Exception("Cannot write to file" + std::string(filename));
          fwrite(pn, SURFEL_SIZE, size, file);
          fclose(file);
        }

        template<typename K>
        void fastProjectionCPU(const CGAL::Point_3<K>& vp,
                               CGAL::Point_3<K>& vq,
                               CGAL::Vector_3<K>& vn)
        {
          Vec3Df p(vp.x(), vp.y(), vp.z());
          Vec3Df q(vq.x(), vq.y(), vq.z());
          Vec3Df n(vn.x(), vn.y(), vn.z());
          fastProjectionCPU(p, q, n);

          vq = CGAL::Point_3<K>(q[0], q[1], q[2]);
          vn = CGAL::Vector_3<K>(n[0], n[1], n[2]);
        }

        // Compute, according to the current point sampling stored in FMLS, the MLS projection
        // of p and store the resulting position in q and normal in n.
        void fastProjectionCPU(const Vec3Df& p, Vec3Df& q, Vec3Df& n)
        {
          float sigma_s = PNScale * MLSRadius;
          float sigma_r = bilateralRange;
          Vec3Df g = (p - Vec3Df(grid.getMinMax()[0], grid.getMinMax()[1], grid.getMinMax()[2])) / sigma_s;
          for (unsigned int j = 0; j < 3; j++) {
            g[j] = floor(g[j]);
            if (g[j] < 0.f)
              g[j] = 0.f;
            if (g[j] >= grid.getRes()[j])
              g[j] = grid.getRes()[j] - 1;
          }
          unsigned int minIt[3], maxIt[3];
          for (unsigned int j = 0; j < 3; j++) {
            if (((unsigned int)g[j]) == 0)
              minIt[j] = 0;
            else
              minIt[j] = ((unsigned int)g[j]) - 1;
            if (((unsigned int)g[j]) == (grid.getRes()[j] - 1))
              maxIt[j] = (grid.getRes()[j] - 1);
            else
              maxIt[j] = ((unsigned int)g[j]) + 1;
          }
          Vec3Df c;
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
                  Vec3Df pk(PN[6 * k], PN[6 * k + 1], PN[6 * k + 2]);
                  Vec3Df nk(PN[6 * k + 3], PN[6 * k + 4], PN[6 * k + 5]);
                  weightedPointCombination(p, pk, nk, sigma_s, bilateral, sigma_r, hermite, c, n, sumW);
                }
              }
          if (sumW == 0.f) {
            n = Vec3Df(1.f, 0.f, 0.f);
            q = p;
          }
          else {
            c /= sumW;
            n.normalize();
            q = p.projectOn(n, c);
          }
        }

        // Compute the MLS projection of the list of point stored in pv and store the resulting
        // positions and normal in qv. qv must be preallocated to stroe 6*pvSize float32.
        // The strid indicates the offsets in qv (the defautl value of 3 means that the qv
        // is compact: pv={x0,y0,z0,x1,y1,z1...}. If pv contains also normals for instance,
        // the stride should be set to 6.
        void fastProjectionCPU(const float* pv, unsigned int pvSize,
          float* qv, unsigned int stride = 3)
        {
#pragma omp parallel for
          for (int i = 0; i < int(pvSize); i++) {
            Vec3Df p(pv[stride * i], pv[stride * i + 1], pv[stride * i + 2]);
            Vec3Df q, n;
            for (unsigned int j = 0; j < numIter; j++) {
              q = Vec3Df();
              n = Vec3Df();
              fastProjectionCPU(p, q, n);
              p = q;
            }
            setPNSample(qv, i, q[0], q[1], q[2], n[0], n[1], n[2]);
          }

        }

        // Brute force version. O(PNSize) complexity. For comparison only.
        void projectionCPU(const Vec3Df& x, Vec3Df& q, Vec3Df& n)
        {
          float sigma_s = MLSRadius * PNScale;
          float sigma_r = bilateralRange;
          Vec3Df p(x);
          for (unsigned int k = 0; k < numIter; k++) {
            Vec3Df c;
            n = Vec3Df();;
            float sumW = 0.f;
            for (unsigned int j = 0; j < PNSize; j++) {
              Vec3Df pj(PN[6 * j], PN[6 * j + 1], PN[6 * j + 2]);
              Vec3Df nj(PN[6 * j + 3], PN[6 * j + 4], PN[6 * j + 5]);
              weightedPointCombination(p, pj, nj, sigma_s, bilateral, sigma_r, hermite, c, n, sumW);
            }
            c /= sumW;
            n.normalize();
            q = p.projectOn(n, c);
            p = q;
          }

        }
        // Brute force version. O(pvSize*PNSize) complexity. For comparison only.
        void projectionCPU(const float* pv, unsigned int pvSize,
          float* qv, unsigned int stride = 3)
        {
#pragma omp parallel for
          for (int i = 0; i < int(pvSize); i++) {
            Vec3Df p(pv[stride * i], pv[stride * i + 1], pv[stride * i + 2]);
            Vec3Df q, n;
            for (unsigned int j = 0; j < numIter; j++) {
              q = Vec3Df();
              n = Vec3Df();
              projectionCPU(p, q, n);
              p = q;
            }
            setPNSample(qv, i, q[0], q[1], q[2], n[0], n[1], n[2]);
          }
        }


       // --------------------------------------------------------------------
       // Filtering by applying MLS projection on the input point set itself.
       // --------------------------------------------------------------------
        // The 'filter*_*' methods apply the MLS projection to the PN samples themselves,
        // providing a low pass (or feature preserving, dependeing on the options)
        // version which can be gathered using 'getFilteredPN ()' afterwards.
        void fastFilterCPU(float* fPN)
        {
          fastProjectionCPU(PN, PNSize, fPN, 6);
        }

        void filterCPU(float* fPN) // Brute force method. O(PNSize^2) complexity. For comparison only.
        {
          projectionCPU(PN, PNSize, fPN, 6);
        }

        // --------------------------------------------------------------
        //  Accessors
        // --------------------------------------------------------------

        // Number of elements of the PN. One elemnt is a 6-float32 chunk.
        inline unsigned int getPNSize() const { return PNSize; }
        inline float* getPN() { return PN; }
        inline const float* getPN() const { return PN; }

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

        class Exception {
        private:
          std::string msg;
        public:
          inline Exception(const std::string& msg) : msg(msg) {}
          virtual ~Exception() {}
          inline const std::string getMessage() const { return std::string("[FMLS][Error]: ") + msg; }
        };

      private:

        void computePNScale()
        {
          Vec3Df c;
          for (unsigned int i = 0; i < PNSize; i++)
            c += Vec3Df(PN[6 * i], PN[6 * i + 1], PN[6 * i + 2]);
          c /= PNSize;
          PNScale = 0.f;
          for (unsigned int i = 0; i < PNSize; i++) {
            float r = Vec3Df::distance(c, Vec3Df(PN[6 * i], PN[6 * i + 1], PN[6 * i + 2]));
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
          {
            cellSize = 1.f;
            LUTSize = 0;
            LUT = NULL;
            indicesSize = 0;
            indices = NULL;
          }
          ~Grid()
          {
            clear();
          }

          void init(float* PN, unsigned int PNSize, float sigma_s)
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
            LUTSize = res[0] * res[1] * res[2];
            unsigned int gridLUTNumOfByte = LUTSize * sizeof(unsigned int);
            LUT = (unsigned int*)malloc(gridLUTNumOfByte);
            memset(LUT, 0, gridLUTNumOfByte);
            unsigned int nonEmptyCells = 0;
            Vec3Df gMin(minMax[0], minMax[1], minMax[2]);
            Vec3Df gMax(minMax[3], minMax[4], minMax[5]);
            for (unsigned int i = 0; i < PNSize; i++) {
              unsigned int index = getLUTIndex(Vec3Df(PN[6 * i], PN[6 * i + 1], PN[6 * i + 2]));
              if (LUT[index] == 0)
                nonEmptyCells++;
              LUT[index]++;
            }
            indicesSize = PNSize + nonEmptyCells;
            indices = (unsigned int*)malloc(indicesSize * sizeof(unsigned int));
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
              Vec3Df p = Vec3Df(PN[6 * i], PN[6 * i + 1], PN[6 * i + 2]);
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
            if (LUT != NULL)
              free(LUT);
            if (indices != NULL)
              free(indices);
            cellSize = 1.f;
            LUTSize = 0;
            LUT = NULL;
            indicesSize = 0;
            indices = NULL;
          }

          // Accessors

          inline const float* getMinMax() const { return minMax; }
          inline const unsigned int* getRes() const { return res; }
          inline float getCellSize() const { return cellSize; }
          inline unsigned int* getLUT() { return LUT; }
          inline const unsigned int* getLUT() const { return LUT; }
          inline unsigned int getLUTSize() const { return LUTSize; }
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
          unsigned int getLUTIndex(const Vec3Df& x) const
          {
            Vec3Df p = (x - Vec3Df(minMax[0], minMax[1], minMax[2])) / cellSize;
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
          inline unsigned int getLUTElement(const Vec3Df& x) const {
            return LUT[getLUTIndex(x)];
          }
          inline unsigned int* getIndices() { return indices; }
          inline const unsigned int* getIndices() const { return indices; }
          inline unsigned int getIndicesSize() const { return indicesSize; }
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
          float minMax[6];
          float cellSize;
          unsigned int  res[3];
          unsigned int LUTSize;
          unsigned int* LUT; // 3D Index Look-Up Table
          unsigned int indicesSize;
          unsigned int* indices; // 3D Grid data
        };


        // --------------------------------------------------------------
        //  Memory Managment
        // --------------------------------------------------------------

        void freeCPUMemory()
        {
          if (PN != NULL)
            freeCPUResource(&PN);
          PNSize = 0;
        }

        // --------------------------------------------------------------
        //  CPU Data
        // --------------------------------------------------------------

        float* PN;
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
               typename C3t3>
      void createMLSSurfaces(Subdomain__FMLS& subdomain_FMLS,
                             Subdomain__FMLS_indices& subdomain_FMLS_indices,
                             const VerticesNormalsMap& vertices_normals,
                             const C3t3& c3t3,
                             int upsample = 0)
      {
//        upsample = 0;
        typedef typename C3t3::Surface_patch_index Surface_index;
        typedef typename C3t3::Subdomain_index     Subdomain_index;
        typedef typename C3t3::Triangulation       Tr;
        typedef typename Tr::Edge                  Edge;
        typedef typename Tr::Vertex_handle         Vertex_handle;
        typedef typename Tr::Geom_traits           Gt;
        typedef typename Gt::Point_3               Point_3;
        typedef typename Gt::Vector_3              Vector_3;

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
          if (vit->in_dimension() == 2)
          {
            const Surface_index si = surface_patch_index(vit, c3t3);
            subdomain_sample_numbers[si]++;
          }
        }

        //if (upsample > 0) {
        //  std::cout << "Up sampling MLS " << upsample << std::endl;
        //  for (C3t3_with_info::Facet_iterator fit = c3t3_with_info.facets_begin(); fit != c3t3_with_info.facets_end(); ++fit) {
        //    Surface_index surf_i = triangulated_domain.make_surface_index(fit->first->subdomain_index(),
        //      fit->first->neighbor(fit->second)->subdomain_index());
        //    if (upsample == 1)
        //      subdomain_sample_numbers[surf_i] ++;
        //    else if (upsample == 2)
        //      subdomain_sample_numbers[surf_i] += 4;
        //  }
        //}

        std::vector< float* > pns;

        int count = 0;
        //Memory allocation for the point plus normals of the point samples
        for (typename SurfaceIndexMap::iterator it = subdomain_sample_numbers.begin();
             it != subdomain_sample_numbers.end(); ++it)
        {
          current_subdomain_FMLS_indices[it->first] = count;
          pns.push_back(new float[it->second * 6]);
          count++;
        }

        std::vector<int> current_v_count(count, 0);
        std::vector<double> point_spacing(count, 0);
        std::vector<int> point_spacing_count(count, 0);

        //Allocation of the PN
        for (Vertex_handle vit : tr.finite_vertex_handles())
        {
          if (vit->in_dimension() == 2)
          {
            const Surface_index surf_i = surface_patch_index(vit, c3t3);

            int fmls_id = current_subdomain_FMLS_indices[surf_i];

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
        typedef boost::unordered_map<Edge_vv, unsigned int> EdgeMapIndex;
        if (upsample == 0)
        {
          EdgeMapIndex edgeMap;

          for (typename C3t3::Facet_iterator fit = c3t3.facets_begin();
               fit != c3t3.facets_end(); ++fit)
          {
            for (int i = 0; i < 2; i++)
            {
              for (int j = i + 1; j < 3; j++)
              {
                Edge edge(fit->first, indices(fit->second,i), indices(fit->second,j));

                Vertex_handle vh0 = edge.first->vertex(edge.second);
                Vertex_handle vh1 = edge.first->vertex(edge.third);
                Edge_vv e = make_vertex_pair(vh0, vh1);
                if ( vh0->in_dimension() == 2
                  && vh1->in_dimension() == 2
                  && edgeMap.find(e) == edgeMap.end())
                {
                  edgeMap[e] = 0;
                  Surface_index surf_i = c3t3.surface_patch_index(*fit);
                  int fmls_id = current_subdomain_FMLS_indices[surf_i];

                  point_spacing[fmls_id] += CGAL::approximate_sqrt(
                    CGAL::squared_distance(point(vh0->point()), point(vh1->point())));
                  point_spacing_count[fmls_id] ++;
                }
              }
            }
          }
        }

//        if (upsample > 0) {
//
//          for (C3t3_with_info::Facet_iterator fit = c3t3_with_info.facets_begin(); fit != c3t3_with_info.facets_end(); ++fit) {
//
//            Surface_index surf_i = triangulated_domain.make_surface_index(fit->first->subdomain_index(),
//              fit->first->neighbor(fit->second)->subdomain_index());
//
//            int fmls_id = current_subdomain_FMLS_indices[surf_i];
//
//            Vertex_handle vhs[3] = { fit->first->vertex(indices[fit->second][0]),
//                                      fit->first->vertex(indices[fit->second][1]),
//                                      fit->first->vertex(indices[fit->second][2]) };
//            K::Vector_3 points[3] = { PointToVector(vhs[0]->point()), PointToVector(vhs[1]->point()), PointToVector(vhs[2]->point()) };
//            K::Vector_3 normals[3] = { vertices_normals[vhs[0]->info()][surf_i], vertices_normals[vhs[1]->info()][surf_i], vertices_normals[vhs[2]->info()][surf_i] };
//
//            std::vector<K::Vector_3> points_to_add;
//            std::vector<K::Vector_3> n_points_to_add;
//
//            //Add the barycenter of the facet
//            K::Vector_3 barycenter = (points[0] + points[1] + points[2]) / 3.;
//            K::Vector_3 n_barycenter = (normals[0] + normals[1] + normals[2]);
//
//            barycenter = (points[0] + points[1] + points[2]) / 3.;
//            n_barycenter = (normals[0] + normals[1] + normals[2]);
//
//            n_barycenter = n_barycenter / CGAL::sqrt((n_barycenter * n_barycenter));
//
//            points_to_add.push_back(barycenter);
//            n_points_to_add.push_back(n_barycenter);
//
//            if (upsample == 1) {
//              for (int i = 0; i < 3; i++) {
//                K::Vector_3 space_1 = barycenter - points[i];
//
//                point_spacing[fmls_id] += CGAL::to_double(CGAL::sqrt(space_1 * space_1));
//                point_spacing_count[fmls_id] ++;
//              }
//            }
//            else if (upsample == 2) {
//              for (int i = 0; i < 3; i++) {
//
//                int i1 = (i + 1) % 3;
//                int i2 = (i + 2) % 3;
//
//                K::Vector_3 p = (barycenter + points[i1] + points[i2]) / 3.;
//                K::Vector_3 n = (n_barycenter + normals[i1] + normals[i2]);
//
//                n = n / CGAL::sqrt(n * n);
//
//                points_to_add.push_back(p);
//                n_points_to_add.push_back(n);
//
//                K::Vector_3 space_1 = p - barycenter;
//                K::Vector_3 space_2 = p - points[i1];
//                K::Vector_3 space_3 = p - points[i2];
//
//                point_spacing[fmls_id] += CGAL::to_double(CGAL::sqrt(space_1 * space_1));
//                point_spacing[fmls_id] += CGAL::to_double(CGAL::sqrt(space_2 * space_2));
//                point_spacing[fmls_id] += CGAL::to_double(CGAL::sqrt(space_3 * space_3));
//
//                point_spacing_count[fmls_id] += 3;
//              }
//            }
//            for (unsigned int i = 0; i < points_to_add.size(); i++) {
//              K::Vector_3& point = points_to_add[i];
//
//              pns[fmls_id][6 * current_v_count[fmls_id]] = point.x();
//              pns[fmls_id][6 * current_v_count[fmls_id] + 1] = point.y();
//              pns[fmls_id][6 * current_v_count[fmls_id] + 2] = point.z();
//
//              K::Vector_3& normal = n_points_to_add[i];
//
//              pns[fmls_id][6 * current_v_count[fmls_id] + 3] = normal.x();
//              pns[fmls_id][6 * current_v_count[fmls_id] + 4] = normal.y();
//              pns[fmls_id][6 * current_v_count[fmls_id] + 5] = normal.z();
//
//              current_v_count[fmls_id]++;
//            }
//          }
//        }


        int nb_of_mls_to_create = 0;
        double average_point_spacing = 0;

        //Cretaing the actual MLS surfaces
        for (typename SurfaceIndexMap::iterator it = current_subdomain_FMLS_indices.begin();
             it != current_subdomain_FMLS_indices.end(); ++it)
        {
          if (current_v_count[it->second] > 3)
          {
            nb_of_mls_to_create++;

            double current_point_spacing = point_spacing[it->second] / point_spacing_count[it->second];
            point_spacing[it->second] = current_point_spacing;

            average_point_spacing += current_point_spacing;
          }
        }

        average_point_spacing = average_point_spacing / nb_of_mls_to_create;

        subdomain_FMLS.resize(nb_of_mls_to_create, FMLS());

        count = 0;
        //Creating the actual MLS surfaces
        for (typename SurfaceIndexMap::iterator it = current_subdomain_FMLS_indices.begin();
             it != current_subdomain_FMLS_indices.end(); ++it)
        {
          if (current_v_count[it->second] > 3)
          {
            double current_point_spacing = point_spacing[it->second];

            //subdomain_FMLS[count].toggleHermite(true);
            subdomain_FMLS[count].setPN(pns[it->second], current_v_count[it->second], current_point_spacing);
            //  subdomain_FMLS[count].toggleHermite(true);
            subdomain_FMLS_indices[it->first] = count;

            count++;
          }
          else {
            std::cout << "Problem of number for MLS : " << current_v_count[it->second] << std::endl;
          }
        }
      }
    }
  }
}

#endif //CGAL_TETRAHEDRAL_REMESHING_FMLS_H
