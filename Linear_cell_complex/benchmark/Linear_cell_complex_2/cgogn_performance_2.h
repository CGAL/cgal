//== INCLUDES =================================================================


#include "performance_2.h"
#include <iostream>

#include "Topology/generic/genericmap.h"
#include "Topology/generic/parameters.h"
#include "Topology/map/embeddedMap2.h"
#include "Topology/generic/traversor2.h"
#include "Topology/generic/traversorCell.h"
#include "Topology/generic/cellmarker.h"

#include "Geometry/vector_gen.h"

#include "Algo/Import/import.h"
#include "Algo/Export/export.h"

#include "Algo/Modelisation/subdivision.h"


//== CLASS DEFINITION =========================================================

using namespace CGoGN ;

struct PFP: public PFP_STANDARD
{
  // definition of the map
  typedef EmbeddedMap2 MAP;
};


class CGoGN_performance_2 : public Performance_test_2
{

private:
  PFP::MAP myMap;
  VertexAttribute<PFP::VEC3> position;
  VertexAttribute<PFP::VEC3> normalF;
  VertexAttribute<PFP::VEC3> normalV;

  unsigned int nbv;

  void display_info()
  {
    std::cout << "#Darts=" << myMap.getNbDarts();
    std::cout << ", #0-cells=" << myMap.getNbOrbits<VERTEX>();
    std::cout << ", #1-cells=" << myMap.getNbOrbits<EDGE>();
    std::cout << ", #2-cells=" << myMap.getNbOrbits<FACE>();
    std::cout << "\t" << std::endl;
  }

public:

  CGoGN_performance_2() : Performance_test_2()
  {
    normalF = myMap.addAttribute<PFP::VEC3, FACE>("normalsF");
    normalV = myMap.addAttribute<PFP::VEC3, VERTEX>("normalsV");
  }

private:
  virtual bool read_mesh(const char* _filename)
  {
    std::vector<std::string> attrNames ;
    if(!Algo::Surface::Import::importMesh<PFP>(myMap, _filename, attrNames))
      return false;

    position = myMap.getAttribute<PFP::VEC3, VERTEX>(attrNames[0]);

    myMap.enableQuickTraversal<VERTEX>();
    myMap.enableQuickTraversal<FACE>();

    /*myMap.enableQuickLocalIncidentTraversal<PFP::MAP, VERTEX, FACE>(myMap);
      myMap.enableQuickLocalIncidentTraversal<PFP::MAP, FACE, VERTEX>(myMap);
      myMap.enableQuickLocalAdjacentTraversal<PFP::MAP, VERTEX, EDGE>(myMap);
    */
    nbv = position.nbElements();

    return true;
  }

  virtual bool write_mesh(const char* _filename)
  {
    return Algo::Surface::Export::exportOFF<PFP>(myMap, position, _filename);
  }

  virtual int circulator_test()
  {
    int counter = 0;

    TraversorV<PFP::MAP> tv(myMap);
    for(Dart dit = tv.begin(); dit != tv.end(); dit = tv.next())
    {
      Traversor2VF<PFP::MAP> trav(myMap, dit);
      for(Dart ditF = trav.begin(); ditF != trav.end(); ditF = trav.next())
        ++counter;
    }

    TraversorF<PFP::MAP> tf(myMap);
    for(Dart dit = tf.begin(); dit != tf.end(); dit = tf.next())
    {
      Traversor2FV<PFP::MAP> trav(myMap, dit);
      for(Dart ditV = trav.begin(); ditV != trav.end(); ditV = trav.next())
        --counter;
    }

    return counter;
  }

  virtual void barycenter_test(bool draw)
  {
    PFP::VEC3 p(0.0,0.0,0.0);

    unsigned int size = position.end();
    unsigned int count = 0;
    for(unsigned int i = position.begin(); i != size; position.next(i))
    {
      p += position[i];
      ++count;
    }
    p /= PFP::REAL(count);

    for(unsigned int i = position.begin(); i != size; position.next(i))
      position[i] -= p;

    if ( draw ) std::cout<<"Barycenter: "<<p<<std::endl;
  }

  virtual void normal_test()
  {
    TraversorF<PFP::MAP> tf(myMap);
    for(Dart dit = tf.begin(); dit != tf.end(); dit = tf.next())
    {
      const typename PFP::VEC3& p1 = position[dit];
      const typename PFP::VEC3& p2 = position[myMap.phi1(dit)];
      const typename PFP::VEC3& p3 = position[myMap.phi_1(dit)];
      normalF[dit] = ((p2 - p1) ^ (p3 - p1)).normalize();
    }

    TraversorV<PFP::MAP> tv(myMap);
    for(Dart dit = tv.begin(); dit != tv.end(); dit = tv.next())
    {
      typename PFP::VEC3 n(0.0);

      Traversor2VF<PFP::MAP> trav(myMap, dit);
      for(Dart ditF = trav.begin(); ditF != trav.end(); ditF = trav.next())
        n += normalF[ditF];

      normalV[dit] = n.normalize();
    }
  }

  virtual void smoothing_test()
  {
    TraversorV<PFP::MAP> tv(myMap);
    for(Dart dit = tv.begin(); dit != tv.end(); dit = tv.next())
    {
      typename PFP::VEC3 p(0.0,0.0,0.0);
      unsigned int c = 0;

      Traversor2VVaE<PFP::MAP> trav2VVaE(myMap, dit);
      for(Dart ditVV = trav2VVaE.begin() ; ditVV != trav2VVaE.end() ; ditVV = trav2VVaE.next())
      {
        p += position[ditVV];
        ++c;
      }
      p /= PFP::REAL(c);
      position[dit] = p;
    }
  }

  virtual void subdivision_test()
  {
    VertexAutoAttribute<PFP::VEC3> new_position(myMap);

    //compute new positions of old vertices
    TraversorV<PFP::MAP> tv(myMap);
    for(Dart d = tv.begin(); d != tv.end(); d = tv.next())
    {
      unsigned int n = 0;
      typename PFP::VEC3 p(0.0,0.0,0.0);
      Traversor2VVaE<PFP::MAP> trav2VVaE(myMap, d);
      for(Dart ditVV = trav2VVaE.begin() ; ditVV != trav2VVaE.end() ; ditVV = trav2VVaE.next())
      {
        p += position[ditVV];
        ++n;
      }
      typename PFP::REAL alpha = (4.0 - 2.0*cos(2.0*M_PI/PFP::REAL(n))) / 9.0;

      new_position[d] = (1.0f-alpha)*position[d] + alpha/PFP::REAL(n)*p;
    }

    Dart last = myMap.end();

    //split faces

    AttributeHandler<Dart, FACE> qtF(&myMap, myMap.getQuickTraversal<FACE>());
    unsigned int nbf = qtF.nbElements();

    for(unsigned int i = qtF.begin(); i != nbf; ++i)
    {
      Dart d = qtF[i];

      unsigned int c = 0;
      typename PFP::VEC3 p(0.0,0.0,0.0);

      //triangule face
      Traversor2FV<PFP::MAP> trav(myMap, d);
      for(Dart ditV = trav.begin(); ditV != trav.end(); ditV = trav.next())
      {
        p += position[ditV];
        ++c;
      }

      p /= PFP::REAL(c);
      Dart cd = Algo::Surface::Modelisation::trianguleFace<PFP>(myMap, d);

      new_position[cd] = p;
    }

    myMap.swapAttributes(position, new_position);

    // flip old edges
    for(Dart d = myMap.begin(); d != last; myMap.next(d))
    {
      if(d < myMap.phi2(d))
        myMap.flipEdge(d);
    }
  }

  virtual void collapse_test()
  {
    myMap.updateQuickTraversal<FACE>();

    //split faces
    DartMarkerNoUnmark me(myMap);

    typename PFP::VEC3 p(0,0,0);

    AttributeHandler<Dart, FACE> qtF(&myMap, myMap.getQuickTraversal<FACE>());
    unsigned int nbf = qtF.end();

    for(unsigned int i = qtF.begin(); i != nbf; ++i)
    {
      Dart d = qtF[i];
      Dart cd = Algo::Surface::Modelisation::trianguleFace<PFP>(myMap, d);
      position[cd] = p;
      me.markOrbit<VERTEX>(cd);
    }

    myMap.disableQuickTraversal<VERTEX>();

    TraversorV<typename PFP::MAP> travV(myMap);
    for(Dart d = travV.begin() ; d != travV.end() ; d = travV.next())
    {
      if(me.isMarked(d))
        myMap.deleteVertex(d);
    }
  }

  virtual void remesh_test()
  {
    unsigned int nbOps = 1000;

    TraversorF<typename PFP::MAP> travF(myMap);
    for(Dart d = travF.begin(); d != travF.end(); d = travF.next())
    {
      Dart cd = Algo::Surface::Modelisation::trianguleFace<PFP>(myMap, d);
      --nbOps;
      if(nbOps == 0)
        break;
    }

    nbOps = 1000;

    TraversorE<typename PFP::MAP> travE(myMap);
    for(Dart d = travE.begin(); d != travE.end(); d = travE.next())
    {
      if(myMap.edgeCanCollapse(d))
      {
        myMap.collapseEdge(d);
        --nbOps;
      }
      if(nbOps == 0)
        break;
    }
  }

};
//=============================================================================
