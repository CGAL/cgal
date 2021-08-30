//== INCLUDES =================================================================
#include "performance_3.h"
#include <iostream>

#include "Topology/generic/parameters.h"
#include "Topology/map/embeddedMap3.h"
#include "Topology/generic/traversor3.h"
#include "Topology/generic/traversorCell.h"
#include "Topology/generic/cellmarker.h"

#include "Geometry/vector_gen.h"

#include "Algo/Modelisation/tetrahedralization.h"

#include "Algo/Import/import.h"
#include "Algo/Export/exportVol.h"

#include "Container/fakeAttribute.h"


//== CLASS DEFINITION =========================================================

using namespace CGoGN ;

#define TEST_QUICK 0

#define TEST_MR 1

struct PFP: public PFP_STANDARD
{
  // definition of the myMap
  typedef EmbeddedMap3 MAP;
};


class CGoGN_performance_3 : public Performance_test_3
{

private:
  PFP::MAP myMap;
  VertexAttribute<PFP::VEC3> position;
  VertexAttribute<PFP::VEC3> position_smooting;
  VertexAttribute<int>       boundary_vertex;

public:

  CGoGN_performance_3() : Performance_test_3()
  {
    position_smooting = myMap.addAttribute<PFP::VEC3, VERTEX>("positionSmoothing");
  }


private:
  void display_info()
  {
    std::cout << "#Darts=" << myMap.getNbDarts();
    std::cout << ", #0-cells=" << myMap.getNbOrbits<VERTEX>();
    std::cout << ", #1-cells=" << myMap.getNbOrbits<EDGE>();
    std::cout << ", #2-cells=" << myMap.getNbOrbits<FACE>();
    std::cout << ", #3-cells=" << myMap.getNbOrbits<VOLUME>();
    std::cout << "\t" << std::endl;
  }

  Dart getShortestEdge()
  {
    double weight = (std::numeric_limits<double>::max)();
    Dart dart = NIL;
    bool boundary=false;

    TraversorE<PFP::MAP> te(myMap);
    for(Dart dit = te.begin() ; dit != te.end() ; dit = te.next())
    {
      Dart dit1 = myMap.phi1(dit);

      boundary=false;
      if (boundary_vertex[dit]==0)
      {
        if (myMap.isBoundaryVertex(dit))
        {
          boundary=true;
          boundary_vertex[dit]=1;
        }
        else
        {
          boundary_vertex[dit]=2;
        }
      }
      else
      {
        boundary = (boundary_vertex[dit]==1);
      }
      if (boundary_vertex[dit1]==0)
      {
        if (myMap.isBoundaryVertex(dit1))
        {
          boundary=true;
          boundary_vertex[dit1]=1;
        }
        else
        {
          boundary_vertex[dit1]=2;
        }
      }
      else
      {
        boundary = (boundary_vertex[dit1]==1);
      }

      if (boundary) continue;

      PFP::VEC3 p1 = position[dit];
      PFP::VEC3 p0 = position[dit1];
      PFP::VEC3 v = (p1 - p0);

      double w = sqrt(v.norm2());
      if(w < weight)
      {
        weight = w;
        dart = dit;
      }
    }

    return dart;
  }

private:
  virtual bool read_mesh(const char* _filename)
  {
    std::vector<std::string> attrNames ;
    Algo::Volume::Import::importMesh<PFP>(myMap, _filename, attrNames);
    position = myMap.getAttribute<PFP::VEC3, VERTEX>(attrNames[0]);

    myMap.enableQuickTraversal<VERTEX>();
    myMap.enableQuickTraversal<VOLUME>();

    // if enabled, don't forget activate disable functions in split_tet
    myMap.enableQuickIncidentTraversal<PFP::MAP, VERTEX, VOLUME>();
    myMap.enableQuickIncidentTraversal<PFP::MAP, VOLUME, VERTEX>();

    myMap.enableQuickAdjacentTraversal<PFP::MAP, VERTEX, EDGE>();
    myMap.enableQuickAdjacentTraversal<PFP::MAP, VERTEX, VOLUME>();

    return true;
  }

  virtual bool write_mesh(const char* _filename)
  {
    //    Algo::Volume::Export::exportTetmesh<PFP>(myMap, position, _filename);
    return true;
  }

  virtual int circulator_test()
  {
    int counter = 0;


    //for each vertex enumerate its incident volumes
    TraversorV<PFP::MAP> tv(myMap);
    for(Dart dit = tv.begin() ; dit != tv.end() ; dit = tv.next())
    {
      Traversor3VW<PFP::MAP> tvw(myMap,dit);
      for(Dart ditvw = tvw.begin() ; ditvw != tvw.end() ; ditvw = tvw.next())
      {
        ++counter;
      }
    }

    //for each volumes enumerate its vertices
    TraversorW<PFP::MAP> tw(myMap);
    for(Dart dit = tw.begin() ; dit != tw.end() ; dit = tw.next())
    {
      Traversor3WV<PFP::MAP> twv(myMap,dit);
      for(Dart ditwv = twv.begin() ; ditwv != twv.end() ; ditwv = twv.next())
      {
        --counter;
      }
    }
    return counter;
  }


  virtual int circulator2_test()
  {
    int counter = 0;
    //for each vertex enumerate its incident volumes
    TraversorV<PFP::MAP> tv(myMap);
    for(Dart dit = tv.begin() ; dit != tv.end() ; dit = tv.next())
    {
      Traversor3VVaW<PFP::MAP> twv(myMap,dit);
      for(Dart ditwv = twv.begin() ; ditwv != twv.end() ; ditwv = twv.next())
      {
        ++counter;
      }
    }
    return counter;
  }

  virtual void barycenter_test(bool draw)
  {
    PFP::VEC3 sum(0.0, 0.0, 0.0);

    TraversorW<PFP::MAP> tw(myMap);
    for(Dart dit = tw.begin() ; dit != tw.end() ; dit = tw.next())
    {
      PFP::VEC3 p = Algo::Surface::Geometry::volumeCentroid<PFP>(myMap,dit,position);
      sum += p;
    }

    if ( draw ) std::cout<<"CGoGn::barycenter: "<<sum<<std::endl;
  }

  virtual void smoothing_test()
  {
    //laplacian smoothing
    TraversorV<PFP::MAP> tv(myMap);
    for(Dart dit = tv.begin() ; dit != tv.end() ; dit = tv.next())
    {
      PFP::VEC3 p(0.0);
      unsigned int c = 0;

      Traversor3VVaE<typename PFP::MAP> trav3VVaE(myMap, dit);
      for(Dart dit3VVaF = trav3VVaE.begin() ; dit3VVaF != trav3VVaE.end() ; dit3VVaF = trav3VVaE.next())
      {
        p += position[dit3VVaF];
        ++c;
      }
      p /= double(c);

      position_smooting[dit] = p;
    }

    myMap.swapAttributes(position, position_smooting);
  }

  virtual void split_tet_test()
  {
    myMap.disableQuickTraversal<VOLUME>();
    myMap.disableQuickTraversal<VERTEX>();

    myMap.disableQuickIncidentTraversal<VERTEX, VOLUME>();
    myMap.disableQuickIncidentTraversal<VOLUME, VERTEX>();

    myMap.disableQuickAdjacentTraversal<VERTEX, EDGE>();
    myMap.disableQuickAdjacentTraversal<VOLUME, VERTEX>();

    TraversorW<typename PFP::MAP> tW(myMap);
    for(Dart dit = tW.begin() ; dit != tW.end() ; dit = tW.next())
    {
      typename PFP::VEC3 volCenter(0.0);
      volCenter += position[dit];
      volCenter += position[myMap.phi1(dit)];
      volCenter += position[myMap.phi_1(dit)];
      volCenter += position[myMap.phi_1(myMap.phi2(dit))];
      volCenter /= 4;

      Dart dres = Algo::Volume::Modelisation::Tetrahedralization::flip1To4<PFP>(myMap, dit);
      position[dres] = volCenter;
    }
  }

  virtual void collapse_test(unsigned int n)
  {
    boundary_vertex = myMap.addAttribute<int, VERTEX>("boundaryVertex");
    TraversorE<PFP::MAP> te(myMap);
    for(Dart dit = te.begin() ; dit != te.end() ; dit = te.next())
      boundary_vertex[dit]=0;

    for(unsigned int i = 0; i < n; ++i)
    {
      Dart dit = getShortestEdge();

      if(dit == NIL)
      {
        std::cerr << "No valid edge anymore, aborting at step "<<i << std::endl;
        return;
      }

      boundary_vertex[dit]=0;
      boundary_vertex[myMap.phi_1(dit)]=0;

      myMap.collapseEdge(dit);
    }
  }

};
//=============================================================================
