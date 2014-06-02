#ifndef MY_MINSK_HEADER
#define MY_MINSK_HEADER

#include <CGAL/basic.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Origin.h>
#include <CGAL/Arrangement_with_history_2.h> 
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/IO/Arr_with_history_iostream.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include "Arr_SegmentData_traits.h"
#include <fstream>
//#include "Graphics.h"

#include <list>
#include <set>
#include <utility>
#include <algorithm>
#include <iterator>
#include <valarray>
#include <boost/unordered_set.hpp>
#include <queue>

#include <ostream>
#include "ICollisionDetector.h"
#include "NaiveCollisionDetector.h"
#include "SweepCollisionDetection.h"
#include "AABB_Collision_detector.h"
#include <boost/unordered_map.hpp>
#include <boost/timer.hpp>

// OMP
//#include <omp.h>


#define HE_WRITE 0

//#ifdef _DEBUG
#define BLA
#ifdef _DEBUG
	#define SHOW_STAGES 1
#else
	#define SHOW_STAGES 0
#endif

#define WRITE_ARR 0

namespace CGAL{

struct Less_than_handle {
  template <typename Type>
  bool operator()(Type s1, Type s2) const { return (&(*s1) < &(*s2)); }
};	


/*
struct Convseg_Less_than{
	bool operator()(ConvS s1, Type s2) const { return (&(*s1) < &(*s2)); }

};
*/

template <class HalfedgeBase_>
	class Arr_map_halfedge : public HalfedgeBase_{
	public:
		//std::string name, type;
		bool visited;
		bool isDegenerate;
		int loopNumber;
	};


template <class Traits_,
	  class HalfedgeBase_ = Arr_halfedge_base<typename Traits_::X_monotone_curve_2> > class Arr_my_extended_dcel : 
	public Arr_dcel_base<Arr_vertex_base<typename Traits_::Point_2>,
		Arr_map_halfedge<HalfedgeBase_>,
		Arr_face_base>
	{ 
	public:

		 template<typename T>
		 class rebind
		  {
			//typedef typename VertexBase_::template rebind
			//					<typename T::Point_2>              Rebind_vertex;
			//typedef typename Rebind_vertex::other                  Vertex_base;
			typedef typename HalfedgeBase_::template rebind
								<typename T::X_monotone_curve_2>   Rebind_halfedge;
			typedef typename Rebind_halfedge::other                Halfedge_base;
   
		  public:

			typedef Arr_my_extended_dcel<T,Halfedge_base>          other;
		  };

	};

template <class Kernel_, class Container_>
class Minkowski_sum_by_convolution_lien_2{
public:

  typedef Kernel_                                        Kernel;
  typedef CGAL::Polygon_2<Kernel, Container_>            Polygon_2;

public:

  // Kernel types:
  typedef typename Kernel::Point_2                       Point_2;
  typedef typename Kernel::Vector_2                      Vector_2;
  typedef typename Kernel::Direction_2                   Direction_2;
  
  // Kernel functors:
  typedef typename Kernel::Equal_2                       Equal_2;
  typedef typename Kernel::Construct_translated_point_2  Translate_point_2;
  typedef typename Kernel::Construct_vector_2            Construct_vector_2;
  typedef typename Kernel::Construct_direction_2         Construct_direction_2;
  typedef typename Kernel::Construct_opposite_line_2     Opposite_line_2;
  typedef typename Kernel::Orientation_2                 Compute_orientation_2;
  typedef typename Kernel::Compare_xy_2                  Compare_xy_2;
  typedef typename Kernel::Counterclockwise_in_between_2 Ccw_in_between_2;
  typedef typename Kernel::Angle_2						 Compute_Angle_2;
  typedef typename Kernel::Compare_x_2					 Compare_x_2;
  typedef typename Kernel::Is_vertical_2					 Is_vertical_2;
  typedef typename Kernel::Compute_x_2					 Compute_x_2;
  typedef typename Kernel::Compute_y_2					 Compute_y_2;

  // Polygon-related types:
  typedef typename Polygon_2::Vertex_circulator          Vertex_circulator;
  typedef std::pair<Vertex_circulator, unsigned int>     Vertex_ref;
  typedef std::pair<Vertex_ref, Vertex_ref>              Anchor;
  typedef std::list<Anchor>                              Anchors_queue;

  // Traits-related types:
  //typedef Arr_segment_traits_2<Kernel>                    Traits_2;//Segment_traits_2;
  typedef Arr_segment_traits_2<Kernel>                    Traits_2_A;//Segment_traits_2;
  typedef Arr_SegmentData_traits<Traits_2_A>				Traits_2;
  //typedef Arr_curve_data_traits_2<Traits_2_B,CGAL::Comparison_result> Traits_2;
  // used
  // typedef Arr_curve_data_traits_2<Traits_2_A,CGAL::Comparison_result> Traits_2;

   //typedef Arr_SegmentData_traits<Traits_2_B>				Traits_2;
//  typedef Arr_labeled_traits_2<Segment_traits_2>          Traits_2; 
  
  typedef typename Traits_2_A::Segment_2		  Base_Segment_2;
  typedef typename Traits_2::X_monotone_curve_2   Segment_2;
  //typedef typename Traits_2::X_monotone_curve_2           Labeled_segment_2;
  typedef std::list<Segment_2>                    Segments_list;
  //typedef std::list<Segments_list::iterator>				  Segments_itr_list;


  
  //typedef CGAL::Arr_default_dcel<Traits_2> Dcel;
	typedef CGAL::Arr_my_extended_dcel<Traits_2> Dcel;
  
//  typedef Arr_face_extended_dcel<Traits_2, int>          Dcel;
  typedef CGAL::Arrangement_with_history_2<Traits_2, Dcel>            Arrangement_history_2;
  typedef typename Arrangement_history_2::Halfedge 					Halfedge;
  typedef typename Arrangement_history_2::Vertex 					Vertex;
  typedef typename Arrangement_history_2::Vertex_iterator 					Vertex_iterator;
  typedef typename Arrangement_history_2::Halfedge_iterator 					Halfedge_iterator;
  typedef typename Arrangement_history_2::Edge_iterator 					Edge_iterator;
  typedef typename Arrangement_history_2::Halfedge_handle 					Halfedge_handle; 
   typedef typename Arrangement_history_2::Vertex_handle					Vertex_handle;
  typedef typename Arrangement_history_2::Face_iterator 					Face_iterator;
  typedef typename Arrangement_history_2::Face_handle 					Face_handle;
  //typedef typename Arrangement_history_2::Face_const 					Face_handle;
  typedef typename Arrangement_history_2::Hole_iterator 					Hole_iterator;
  typedef typename Arrangement_history_2::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
  typedef typename Arrangement_history_2::Ccb_halfedge_circulator			Ccb_halfedge_circulator;

  typedef typename Arrangement_history_2::Originating_curve_iterator  Originating_curve_iterator;
  typedef std::pair<int,int> StatePair;

  typedef std::set<Halfedge_handle,Less_than_handle> Edges_set; 
  typedef std::set<Face_handle,Less_than_handle> Faces_set; 
//  typedef Union_of_segment_cycles_2<Traits_2, Polygon_2>  Union_2;

  // Data members:
  Equal_2                 f_equal;
  Translate_point_2       f_add;
  Construct_vector_2      f_vector;
  Construct_direction_2   f_direction;
  Opposite_line_2         f_opp_line;
  Compute_orientation_2   f_orientation;
  Compare_xy_2            f_compare_xy;
  Ccw_in_between_2        f_ccw_in_between;
  Compute_Angle_2		  f_angle;
  Is_vertical_2			  f_is_vertical;
  Compute_x_2			  f_compute_x;
  Compute_y_2			  f_compute_y;

 // Compare_x_2			  f_compare_x;
  //Traits_2::Compare_endpoints_xy_2  f_compare_endpoints_xy;
  typename Traits_2::Compare_endpoints_xy_2 f_compare_endpoints_xy;
  typename Traits_2::Compare_y_at_x_2 f_compare_y_at_x;
  typename Traits_2::Compare_x_2 f_compare_x;
	
  friend class ConvSegMapper;
  struct ConvSegment
  {
		Halfedge_handle _he;
		ConvSegment(Halfedge_handle& he):_he(he){}
		ConvSegment(){}
		bool getVisited() const
		{
			return _he->visited;
		}
		
		bool getDegenerate() const
		{
			return _he->isDegenerate;
		}

		int getLoopNum()
		{
			return _he->loopNumber;
		}

		Vertex_handle getSrc()
		{
			return (_he->source());
		}

		Vertex_handle getDst()
		{
			//return Vertex();//_he->target();
			return (_he->target());
		}

		bool operator<(const ConvSegment& rhs) const
		{
			//Less_than_handle o();
			//return o.operator()(_he,rhs._he);
			return Less_than_handle()(_he,rhs._he);
		}

		bool operator==(const ConvSegment& rhs) const
		{
			return !Less_than_handle()(_he,rhs._he) && !Less_than_handle()(rhs._he,_he);
		}

  };

#define FIRST_LOOP 0
  struct ConvSegMapper{
	  Arrangement_history_2* _arr;
	  Minkowski_sum_by_convolution_lien_2* _mink;
	  ConvSegMapper(Arrangement_history_2 *arr,Minkowski_sum_by_convolution_lien_2* mink):_arr(arr),_mink(mink){
		  
	  }

	  ConvSegment getSegment(const Halfedge_handle& he)
	  {
		  return ConvSegment(_mink->getDirAgreeingHalfedge(*_arr,he));
	  }

	  Direction_2 getConvSegDir(const ConvSegment& seg) const
	  {
		  return _mink->getHalfedgeDir(seg._he);
	  }

	  void markVisited(ConvSegment& convSeg,int id)
	  {
		  _mink->setEdgeVisited(*convSeg._he,true,id);
	  }

	  bool isBBiggerThenAWithReagrdToC(const ConvSegment& a,const ConvSegment& b,const ConvSegment& c) const
	  {
		  Direction_2 dir_a = getConvSegDir(a);
		  Direction_2 dir_b = getConvSegDir(b);
		  Direction_2 dir_c = getConvSegDir(c);
		  return _mink->isDirImproving(dir_a,dir_c,dir_b);
	  }

	  void getNeighbouringSegments(Vertex_handle v,std::list<ConvSegment>& outSegments,std::list<ConvSegment>& inSegmets)
	  {
		  std::list<Halfedge_handle> inList,outList;
		  _mink->getEdgesFromVertex(*_arr,v,inList,outList);
		  typename std::list<Halfedge_iterator>::const_iterator itr;

		  for (itr = inList.begin();itr!= inList.end();++itr)
		  {
			  inSegmets.push_back(getSegment(*itr));
		  }

		  for (itr = outList.begin();itr!= outList.end();++itr)
		  {
			  outSegments.push_back(getSegment(*itr));
		  }
	  }

	  static bool getSegVisited(const ConvSegment& seg)
	  {
		  return (seg.getVisited()==true || seg.getDegenerate());
	  }

	  static bool getSegNotVisited(const ConvSegment& seg)
	  {
		  return (seg.getVisited()==false && !seg.getDegenerate());
	  }

	  double getSignedAngle(const ConvSegment& enter,const ConvSegment& exit)
	  {
			return _mink->getSignedAngle(enter._he,exit._he);
	  }

	  void filterNonVisitedSegments(const std::list<ConvSegment>& inputList,std::list<ConvSegment>& outList, std::list<ConvSegment>& visitedSegmentsList)
	  {
		  int out_size = count_if(inputList.begin(),inputList.end(),&ConvSegMapper::getSegNotVisited);
		  outList.resize(out_size);
		  remove_copy_if(inputList.begin(),inputList.end(),outList.begin(),&ConvSegMapper::getSegVisited);
		  visitedSegmentsList.resize(inputList.size() - out_size);
		  remove_copy_if(inputList.begin(),inputList.end(),visitedSegmentsList.begin(),&ConvSegMapper::getSegNotVisited);
	  }

	  struct SegCompare{
		  ConvSegment _incomingSeg;
		  ConvSegMapper* _mapperInstance;
		  SegCompare(ConvSegment incomingSeg,ConvSegMapper* mapperInstance):_incomingSeg(incomingSeg),_mapperInstance(mapperInstance){}
		  bool operator()(ConvSegment a,ConvSegment b)
		  {
			  return _mapperInstance->isBBiggerThenAWithReagrdToC(a,b,_incomingSeg);
		  }
	  };

	  ConvSegment getMaximalEdge(std::list<ConvSegment> outgoingSegs,ConvSegment incomingSeg)
	  {
		  SegCompare seg_comp(incomingSeg,this);
		  return *max_element(outgoingSegs.begin(),outgoingSegs.end(),seg_comp);
	  }

	  bool checkLoopClosed(const Vertex_handle& v,ConvSegment& startingSeg,int id)
	  {
		  Halfedge_handle h;
		  bool notFound = _mink->checkOutgoingNotVisited(*_arr,*v,h,id);
		  if (h != Halfedge_handle())
			startingSeg = getSegment(h);
		  return !notFound;
	  }

	  template<typename T> void fillEdgesSet(T& edges_set)
	  {
		Edge_iterator  itr;
		//Edges_set edges_set;
		for (itr = _arr->edges_begin();itr!=_arr->edges_end();++itr){
			//setEdgeVisited(*itr,false,-1);
			_mink->setEdgeVisited(*itr,false,-1);

			if (!itr->isDegenerate)
			//*oitr++ = getSegment(itr);
				edges_set.insert(getSegment(itr));
		}
	  }


	  ConvSegment getOuterSegment()
	  {
			Face_iterator startFace = _arr->unbounded_face();
			Halfedge_iterator perimiterFace = *(startFace -> holes_begin());
			return getSegment(perimiterFace);
	  }

	  void removeSegFromArr(const ConvSegment& seg)
	  {
		  _arr->remove_edge(seg._he);
	  }

	  void removeRangeFromArr(std::list<ConvSegment>& segsToRemove)
	  {
		  //for_each(segsToRemove.begin(),segsToRemove.end(),&ConvSegMapper::removeSegFromArr);
		  for (typename std::list<ConvSegment>::iterator itr = segsToRemove.begin();itr!=segsToRemove.end();++itr)
		  {
				removeSegFromArr(*itr);
		  }
	  }

	  
  };

  class TraversalManager;
  friend class ConvMovement;
  struct ConvMovement{
		ConvMovement(ConvSegMapper* mapperInstance,int id,ConvSegment& startedge,TraversalManager* managerInstance):_id(id),_mapperInstance(mapperInstance),_managerInstance(managerInstance)
		{
			_currEdge = startedge;
			_traversedEdges.push_back(_currEdge);
			_mapperInstance->markVisited(startedge,id);
			_traversing = true;
			_loop = NULL;
		}

		void Traverse()
		{
			
			while (_traversing)
			{

				Vertex_handle dst = _currEdge.getDst();
				std::list<ConvSegment> inList;
				std::list<ConvSegment> outList;
				_mapperInstance->getNeighbouringSegments(dst,outList,inList);
				std::list<ConvSegment> filteredSegments,visitedSegments;
				_mapperInstance->filterNonVisitedSegments(outList,filteredSegments,visitedSegments);
			/*	if (outList.size()>1)
				{
					drawTraversal();

				}*/
				if (checkCloseLoop())
					closeLoopEvent(filteredSegments);

				if (!_traversing)
					break;

				// Check if movement may proceed.
				if (filteredSegments.size() >0 )
				{
					ConvSegment c = _mapperInstance->getMaximalEdge(filteredSegments,_currEdge);
					// Check that there is no visited outgoing edge, with diffrenet loop id (ie closed loop),
					// which is better than c.
					typename std::list<ConvSegment>::iterator itr = visitedSegments.begin();
					for (;itr != visitedSegments.end();++itr)
					{
						if (_mapperInstance->isBBiggerThenAWithReagrdToC(c,*itr,_currEdge) && itr->getLoopNum() != _id && itr->getLoopNum() != FIRST_LOOP)
						{
							closeEvent();
							_traversing = false;
							return;
						}
					}

					double angle = _mapperInstance->getSignedAngle(_currEdge,c);
					_angles.push_back(angle);
					_anglesSum += angle;
					_currEdge = c;
					_mapperInstance->markVisited(c,_id);
					_traversedEdges.push_back(c);
				}
				else
				{
					closeEvent();
					_traversing = false;
				}
			}
		}
// #define SHOW_LOOPS_CONST
		void closeLoopEvent(std::list<ConvSegment>& outGoingOptions)
		{
			#ifdef SHOW_LOOPS_CONST
				drawTraversal();
			#endif
			_managerInstance->handleCloseLoopEvent(*_loop,outGoingOptions);
			_angles.pop_back();
		}

		void closeEvent()
		{
			#ifdef SHOW_LOOPS_CONST
				drawTraversal();
			#endif
			_managerInstance->handleCloseEvent();
		}

		bool checkCloseLoop()
		{
			ConvSegment loopStart;
			bool found = _mapperInstance->checkLoopClosed(_currEdge.getDst(),loopStart,_id);
			if (found)
			{
				if (_loop!=NULL)
					delete _loop;
				_loop =  new std::pair<ConvSegment,ConvSegment>(loopStart,_currEdge);
				double angle = _mapperInstance->getSignedAngle(_currEdge,loopStart);
				_angles.push_back(angle);
			}
			return found;
		}

		void stopTraverse()
		{
			_traversing = false;
		}

    /*
		void drawTraversal()
		{
		
			global_graphics->clear();

			for (typename std::list<ConvSegment>::iterator itr = _traversedEdges.begin();itr != _traversedEdges.end(); ++itr)
			{
				global_graphics->draw_edge<Kernel>((itr->_he)->curve(),QColor(0,255,0),true);
			}
			global_graphics->display();
		
		}
        */

		std::list<ConvSegment>& getTraversedEdges()
		{
			return _traversedEdges;
		}

		std::list<double>& getAngles()
		{
			return _angles;
		}
  private:
		int _id;
		ConvSegment _currEdge;
		std::list<ConvSegment> _traversedEdges;
		std::list<double> _angles;
		ConvSegMapper* _mapperInstance;
		TraversalManager* _managerInstance;
		bool _traversing;
		double _anglesSum;
		std::pair<ConvSegment,ConvSegment>* _loop;
  };

  friend class TraversalManager;
  class EdgesStore;
  class LoopsTracker;
  struct TraversalManager
  {
	  TraversalManager(ConvSegMapper* mapperInstance):_mapperInstance(mapperInstance)
	  {
		  _loopId = 0;
		  _activeMovment = NULL;
	  }
	  void traverseLoops()
	  {
		  EdgesStore edges_db(_mapperInstance);
		  _currEdgesStore = &edges_db;
		  int s = edges_db.getSize();
		  ConvSegment startEdge = _mapperInstance->getOuterSegment();
		  traceLoop(edges_db,startEdge,true);

		  while (!edges_db.isEmpty())
		  {
			  int s = edges_db.getSize();
			  ++_loopId;
			  startEdge = edges_db.getEdge();
			  traceLoop(edges_db,startEdge,false);
		  }

	  }

	  void traceLoop(EdgesStore& edgesDB,ConvSegment& startSeg,bool isOuter)
	  {
		 /* if (_activeMovment != NULL)
			  delete _activeMovment;*/
		  _activeMovment = new ConvMovement(_mapperInstance,_loopId,startSeg,this);

		  /*if (_activeLoopTracker != NULL)
			  delete _activeLoopTracker;*/
		  _activeLoopTracker = new LoopsTracker(&(_activeMovment->getTraversedEdges()),&(_activeMovment->getAngles()),isOuter);

		  _activeMovment->Traverse();

		  // Delete objects
		  delete _activeMovment;
		  delete _activeLoopTracker;
	  }
	
	  void handleCloseLoopEvent(std::pair<ConvSegment,ConvSegment>& loop,std::list<ConvSegment>& outGoingOptions)
	  {
		  _activeLoopTracker->addLoop(loop);

		  // check for improvment: (TODO: maybe refactor)
		  // This checks if the first segment in loop (actually second edge goes into first) is a worse choice then *itr , w.r.t second.
		  // This means the loop might be improved.
		  if (_activeLoopTracker->hasLoops())
		  {
			  bool stopItrating = true;
			  for (typename std::list<ConvSegment>::iterator itr = outGoingOptions.begin();itr != outGoingOptions.end();++itr)
			  {
				  if (_mapperInstance->isBBiggerThenAWithReagrdToC(loop.first,*itr,loop.second))
				  {
					  stopItrating = false;
					  break;
				  }
			  }

			  if (stopItrating)
			  {
				  handleCloseEvent();
				  _activeMovment->stopTraverse();
			  }
		  }
	  }

	  void handleCloseEvent()
	  {
		  // remove all traversed edges from segments list
		  
		  std::list<ConvSegment>& traversedEdges = _activeMovment->getTraversedEdges();
		  _currEdgesStore->removeRange(traversedEdges.begin(),traversedEdges.end());

		  // remove anything but loop from arrangment.
		  std::list<ConvSegment>& toRemoveFromArr = _activeLoopTracker->getNonLoopSegmentsList();

		  _mapperInstance->removeRangeFromArr(toRemoveFromArr);
	  }


  private:
	  ConvSegMapper* _mapperInstance;
	  ConvMovement* _activeMovment;
	  LoopsTracker* _activeLoopTracker;
	  EdgesStore* _currEdgesStore;
	  int _loopId;
  };

  struct LoopsTracker
  {
		LoopsTracker(std::list<ConvSegment>* segmentsList,std::list<double>* anglesList,bool outerLoop):_segmentsList(segmentsList),_anglesList(anglesList),_outerLoop(outerLoop)
		{
			_hasLoop = false;
			_loopsCreated = false;
		}

		void addLoop(std::pair<ConvSegment,ConvSegment>& loop)
		{
			typename std::list<ConvSegment>::iterator loop_begin = find(_segmentsList->begin(),_segmentsList->end(),(loop.first));
			typename std::list<ConvSegment>::iterator loop_end = find(_segmentsList->begin(),_segmentsList->end(),(loop.second));
/*			if (!_hasLoop)
			{
				_loopBegin = loop_begin;
				_loopEnd = loop_end;
				_hasLoop = true;
			}
			else
			{*/
				double angleSum = sumLoopAngles(loop_begin,loop_end);
				if ((angleSum >0 && _outerLoop) || (angleSum <0 && !_outerLoop))
				{// update loop
					_hasLoop = true;
					_loopBegin = loop_begin;
					_loopEnd = loop_end;
				}
			//}
		}

		double sumLoopAngles(typename std::list<ConvSegment>::iterator& loopBegin,typename std::list<ConvSegment>::iterator& loopEnd)
		{
			int begin = std::distance(_segmentsList->begin(),loopBegin);
			int end = std::distance(_segmentsList->begin(),loopEnd);
			std::list<double>::iterator itr_begin = _anglesList->begin();
			advance(itr_begin,begin);
			std::list<double>::iterator itr_end = _anglesList->begin();
			advance(itr_end,end + 1);
			double sum = 0;
			for (std::list<double>::iterator itr = itr_begin;itr != itr_end;++itr)
			{
				sum += *itr;
			}
			return sum;
		}
		
		std::list<ConvSegment>& getLoopSegmentsList()
		{
			//_loopsSegmentsList;
			//copy(_loopBegin,_loopEnd+1,_loopSegmentsList.begin());
			if (!_loopsCreated)
				createLoopsSegmentsLists();
			return _loopSegmentsList;
		}

		std::list<ConvSegment>& getNonLoopSegmentsList()
		{
			/*if (_loopsSegmentsList.size()==0)
				getLoopSegmentsList();
			//_loopsSegmentsList;
			sort(_loopsSegmentsList.begin(),_loopsSegmentsList.end());
			std::list<ConvSegment> segmentsListCopy(*_segmentsList);
			sort(segmentsListCopy.begin(),segmentsListCopy.end());
			set_difference(segmentsListCopy.begin(),segmentsListCopy.end(),_loopsSegmentsList.begin(),_loopsSegmentsList.end(),_nonLoopSegmentsList.begin());
			return _nonLoopSegmentsList;
			*/
			if (!_loopsCreated)
				createLoopsSegmentsLists();
			return _nonLoopSegmentsList;
		}

		void createLoopsSegmentsLists()
		{
			copy(_segmentsList->begin(),_segmentsList->end(),back_inserter(_nonLoopSegmentsList));
			if (_hasLoop)
			{
				/*std::list<ConvSegment>*/ 
				typename std::list<ConvSegment>::iterator loop_begin = find(_nonLoopSegmentsList.begin(),_nonLoopSegmentsList.end(),*_loopBegin);
				typename std::list<ConvSegment>::iterator loop_end = find(_nonLoopSegmentsList.begin(),_nonLoopSegmentsList.end(),*_loopEnd);
				//std::list<ConvSegment> _nonLoopSegmentsList(*_segmentsList);
				_loopSegmentsList.splice(_loopSegmentsList.begin(),_nonLoopSegmentsList,loop_begin,++loop_end);
			}
			
		}

		bool hasLoops()
		{
			return _hasLoop;
		}

  private:
		typename std::list<ConvSegment>::iterator _loopBegin,_loopEnd; 
		std::list<ConvSegment>* _segmentsList;
		std::list<double>* _anglesList;
		std::list<ConvSegment> _loopSegmentsList;
		std::list<ConvSegment> _nonLoopSegmentsList;
		bool _outerLoop;
		bool _hasLoop;
		bool _loopsCreated;
  };

  struct EdgesStore
  {
	  EdgesStore(ConvSegMapper* mapper)
	  {
		_edgesSet = new std::set<ConvSegment>();
		mapper->fillEdgesSet(*_edgesSet);
	//	sort(_edgesSet.begin(),_edgesSet.end());
	  }

	  template <typename InputIterator> void removeRange(InputIterator start,InputIterator end)
	  {
		  //std::vector<InputIterator::value_type> vec(start,end);
		  //sort(vec.begin(),vec.end());
		  //sort(start,end);
		  //std::vector<InputIterator::value_type>* seg_set = new std::vector<InputIterator::value_type>();
		  //set_difference(_edgesSet->begin(),_edgesSet->end(),vec.begin(),vec.end(),back_inserter(*seg_set));
		  //delete _edgesSet;
		  //_edgesSet = new std::set<ConvSegment>(seg_set->begin(),seg_set->end());
		  
		  for (InputIterator itr = start; itr != end;++itr)
		  {
			  _edgesSet->erase(*itr);
		  }
		  
	  }

	  template <typename InputIterator> void removeRange_bak(InputIterator start,InputIterator end)
	  {
		  std::vector<typename InputIterator::value_type> vec(start,end);
		  sort(vec.begin(),vec.end());
		  //sort(start,end);
		  std::vector<typename InputIterator::value_type>* seg_set = new std::vector<typename InputIterator::value_type>();
		  set_difference(_edgesSet->begin(),_edgesSet->end(),vec.begin(),vec.end(),back_inserter(*seg_set));
		  delete _edgesSet;
		  _edgesSet = new std::set<ConvSegment>(seg_set->begin(),seg_set->end());
	  }

	  ConvSegment getEdge()
	  {
		  return *(_edgesSet->begin());
	  }

	  bool isEmpty()
	  {
		  return (_edgesSet->size()==0);
	  }

	  int getSize()
	  {
		  return _edgesSet->size();
	  }

  private:
	  std::set<ConvSegment>* _edgesSet;
  };

  friend class DegenerateCassesManager;
  struct DegenerateCassesManager
  {
	  DegenerateCassesManager(Arrangement_history_2* arr,Minkowski_sum_by_convolution_lien_2* mink,Polygon_2* poly1,Polygon_2* poly2,bool isActive):_arr(arr),_mink(mink),_poly1(poly1),_poly2(poly2),_active(isActive)
	  {
		  
	  }

	  /*
		Go over edges and check if they are degenerate
	  */
	  void markDegenerateEdges()
	  {
		  Edge_iterator itr = _arr->edges_begin();
		  for(;itr != _arr->edges_end();++itr)
		  {
			  _mink->setEdgeVisited(*itr,false,-1);
			  if (_active)
			  {
				  if (_mink->checkDegenerateEdgeOppositeSegments(*_arr,itr))
				  {
					  if (!_mink->checkSegmentCollisionDetection(*_arr,itr->curve(),*_poly1,*_poly2))
					  {
						  _mink->setEdgeDegenerate(*itr,true);
					  }
					  else
						  _mink->setEdgeDegenerate(*itr,false);
				  }
				  else
					  _mink->setEdgeDegenerate(*itr,false);
			  }
			  else
			  {
				  _mink->setEdgeDegenerate(*itr,false);
			  }
			  
		  }
	  }


	  /*
		Goes over each vertex and checks if it is degenerate.
	  */
	  void findDegenerateBorderVertices()
	  {
		  if (!_active)
			  return;
		  Vertex_iterator itr = _arr->vertices_begin();
		/*  int vertices_num = _arr->number_of_vertices();
		  std::valarray<Vertex_iterator> vertices_arr(vertices_num);
		  int i=0;
		  for (;itr != _arr->vertices_end();++itr)
		  {
			  vertices_arr[i] = itr;
			  ++i;
		  }
		  
		  Point_2 p_end;
			#pragma omp parallel for default(shared),private(p_end)*/
		  for (;itr != _arr->vertices_end();++itr)
		  //for (i=0;i<vertices_num;++i)
		  {
			  //itr = vertices_arr[i];
			//#pragma omp single nowait
			
			if (_mink->checkDegenarateVertexIsIntersectionOfThreeSegments(*_arr,itr))
		//		if (_mink->checkDegenarateVertexIsIntersectionOfThreeSegments(*_arr,vertices_arr[i]))
			{
				Point_2 p_end = itr->point();
			//	p_end = vertices_arr[i]->point();
	/*			double x2 = CGAL::to_double(p_end.x());
				double y2 = CGAL::to_double(p_end.y());
				std::cout << x2 << "," << y2 << std::endl;*/
				if (!_mink->checkCollisionDetection(*_arr,itr->point(),*_poly1,*_poly2))
				//if (!_mink->checkCollisionDetection(*_arr,vertices_arr[i]->point(),_collision_detector,*_poly1,*_poly2))
				{
					//std::cout << "in" << std::endl;
					//#pragma omp critical
					_degenerate_points_list.push_back(itr->point());
					//_degenerate_points_list.push_back(vertices_arr[i]->point());
				}
			}
			//}
		  }
	  }


	  void addDegenerateVerticesToArr()
	  {
		  typename std::list<Point_2>::iterator itr = _degenerate_points_list.begin();
		  for (;itr != _degenerate_points_list.end();++itr)
		  {
				CGAL::insert_point(*_arr,*itr);
				
		  }
	  }

  private:
	  Arrangement_history_2* _arr;
	  Minkowski_sum_by_convolution_lien_2* _mink;
	  std::list<Point_2> _degenerate_points_list;
	  Polygon_2* _poly1;
	  Polygon_2* _poly2;
	  //SweepCollisionDetector<Kernel,Container_> _collision_detector;
	  //NaiveCollisionDetector<Kernel,Container_> _collision_detector;
	  bool _active;
  };

public:
	
   /*! Default constructor. */
  Minkowski_sum_by_convolution_lien_2   ()
  {
    // Obtain kernel functors.
    Kernel                ker;

    f_equal = ker.equal_2_object();
    f_add = ker.construct_translated_point_2_object(); 
    f_vector = ker.construct_vector_2_object();
    f_direction = ker.construct_direction_2_object();
    f_opp_line = ker.construct_opposite_line_2_object();
    f_orientation = ker.orientation_2_object();
    f_compare_xy = ker.compare_xy_2_object();
    f_ccw_in_between = ker.counterclockwise_in_between_2_object();
	f_angle = ker.angle_2_object();
	f_compare_endpoints_xy = Traits_2().compare_endpoints_xy_2_object();
	f_is_vertical = ker.is_vertical_2_object();
	f_compare_x =  Traits_2().compare_x_2_object();
	f_compare_y_at_x = Traits_2().compare_y_at_x_2_object();
	f_compute_x  = ker.compute_x_2_object();
	f_compute_y  = ker.compute_y_2_object();
  }  

   /*!
   * Compute the Minkowski sum of two simple polygons.
   * Note that as the input polygons may not be convex, the Minkowski sum may
   * not be a simple polygon. The result is therefore represented as
   * the outer boundary of the Minkowski sum (which is always a simple polygon)
   * and a container of simple polygons, representing the holes inside this
   * polygon.
   * \param pgn1 The first polygon.
   * \param pgn2 The second polygon.
   * \param sum_bound Output: A polygon respresenting the outer boundary
   *                          of the Minkowski sum.
   * \param sum_holes Output: An output iterator for the holes in the sum,
   *                          represented as simple polygons.
   * \pre Both input polygons are simple.
   * \return A past-the-end iterator for the holes in the sum.
   */
  template <class OutputIterator>
  OutputIterator old_version (const Polygon_2& pgn1,
                             const Polygon_2& pgn2,
                             Polygon_2& sum_bound,
                             OutputIterator sum_holes) 
  {
	CGAL_precondition (pgn1.is_simple());
    CGAL_precondition (pgn2.is_simple());
	CGAL_precondition (pgn1.orientation() == CGAL::COUNTERCLOCKWISE);
    CGAL_precondition (pgn2.orientation() == CGAL::COUNTERCLOCKWISE);
	Polygon_2 revP1 = revPoly(pgn1);
	Polygon_2 p2 = pgn2;
	_aabb_collision_detector = new AABBCollisionDetector<Kernel_,Container_>(p2,revP1);
	
	Segments_list reduced_conv;
	buildReducedConvolution(pgn1,pgn2,reduced_conv);

	Arrangement_history_2 arr;
	buildArrangementFromConv(reduced_conv,arr);

	const Minkowski_sum_by_convolution_lien_2* ptr = this;
	//	ConvSegMapper mapper( &arr,const_cast <Minkowski_sum_by_convolution_lien_2*>(ptr));*/
	DegenerateCassesManager degHandler(&arr,const_cast <Minkowski_sum_by_convolution_lien_2*>(ptr),const_cast <Polygon_2*>(&pgn1),const_cast <Polygon_2*>(&pgn2),true);
	degHandler.findDegenerateBorderVertices();
	degHandler.markDegenerateEdges();
	constructOrientableLoops(arr);

/*
	if (SHOW_STAGES)
	{
		global_graphics->clear();
		draw_arr(arr);
		global_graphics->display();
	}
    */
	Polygon_2 reverse_pgn1 = transform(typename Kernel::Aff_transformation_2(CGAL::Rotation(), 0, -1), pgn1);
	nestedLoopsFilter(arr,reverse_pgn1,pgn2);
	degHandler.addDegenerateVerticesToArr();
    /*
	if (SHOW_STAGES)
	{
		global_graphics->clear();	
		draw_arr(arr);
		global_graphics->display();
	}
    */
	//collisionDetectionFilter(arr);

	//global_graphics->clear();
	
	//draw_arr(arr);
	//global_graphics->display();
	//for (Segments_list::const_iterator itr = reduced_conv.begin();itr != reduced_conv.end();++itr)
	//{
	//	Segment_2 temp = *itr;
	//	Point_2 p_source = temp.source();
	//	Point_2 p_end = temp.target();
	//	double x1 = CGAL::to_double(p_source.x());
	//	double y1 = CGAL::to_double(p_source.y());
	//	double x2 = CGAL::to_double(p_end.x());
	//	double y2 = CGAL::to_double(p_end.y());
	////	CGAL::insert(arr,temp);
	//}
	delete _aabb_collision_detector;
	return (sum_holes);
  }


   template <class OutputIterator>
  OutputIterator operator() (const Polygon_2& pgn1,
                             const Polygon_2& pgn2,
                             Polygon_2& sum_bound,
                             OutputIterator sum_holes) 
  {
	CGAL_precondition (pgn1.is_simple());
    CGAL_precondition (pgn2.is_simple());
	CGAL_precondition (pgn1.orientation() == CGAL::COUNTERCLOCKWISE);
    CGAL_precondition (pgn2.orientation() == CGAL::COUNTERCLOCKWISE);
	Polygon_2 revP1 = revPoly(pgn1);
	Polygon_2 p2 = pgn2;
	boost::timer t_abb;
	_aabb_collision_detector = new AABBCollisionDetector<Kernel_,Container_>(p2,revP1);
	double aabb_build_time = t_abb.elapsed();
	//std::cout << "aabb tree build : " <<aabb_build_time << std::endl;
	Segments_list reduced_conv;
	//buildReducedConvolution(pgn1,pgn2,reduced_conv);
	boost::timer t1;
	buildReducedConvolutionFiberGrid(pgn1,pgn2,reduced_conv);
	//rc_time = t1.elapsed();
	
	//_my_global_counter = reduced_conv.size();
	//std::cout << "Number of convolution segments : " << _my_global_counter << std::endl;
	//std::cout << "reduced conv build : " <<rc_time << std::endl;
	Arrangement_history_2 arr;
	boost::timer t2;
	buildArrangementFromConv(reduced_conv,arr);
	//arr_build_time =  t2.elapsed();
	//std::cout << "buildArrangementFromConv : " << arr_build_time << std::endl; 

	//LogMyMemoryUsage();
	//std::cout << "sizeOfReducedConvArrangement : " << arr.number_of_edges() << std::endl;

	const Minkowski_sum_by_convolution_lien_2* ptr = this;
	//	ConvSegMapper mapper( &arr,const_cast <Minkowski_sum_by_convolution_lien_2*>(ptr));*/
	boost::timer t4;
	DegenerateCassesManager degHandler(&arr,const_cast <Minkowski_sum_by_convolution_lien_2*>(ptr),const_cast <Polygon_2*>(&pgn1),const_cast <Polygon_2*>(&pgn2),true);
	degHandler.findDegenerateBorderVertices();
	degHandler.markDegenerateEdges();
	//degenerate_stage = t4.elapsed();

	boost::timer t3;
    /*
	if (SHOW_STAGES)
	{
		global_graphics->clear();
		draw_arr(arr);
		global_graphics->display();
	}
    */
	Polygon_2 reverse_pgn1 = transform(typename Kernel::Aff_transformation_2(CGAL::Rotation(), 0, -1), pgn1);

	// trace outer loop
	markOutsideLoop(arr,sum_bound);

	// trace holes
//	int num_faces = arr.number_of_faces();
//	std::vector<Face_iterator> faces_itrs(num_faces);
//	int index = 0;
//	for (Face_iterator itr = arr.faces_begin();itr!=arr.faces_end();++itr )
//	{
//		faces_itrs[index] = itr;
//		++index;
//	}
//	omp_set_num_threads(4);
//#pragma omp threadprivate(_aabb_collision_detector)
//#pragma omp parallel for 
//		for (int i = 0;i<num_faces;++i ){
//			handleFace(arr,faces_itrs[i],reverse_pgn1,pgn2,sum_holes);
//		}
	
	
//	Original code here !
	for (Face_iterator itr = arr.faces_begin();itr!=arr.faces_end();++itr ){
			handleFace(arr,itr,reverse_pgn1,pgn2,sum_holes);
	}

	
	std::list<Halfedge_handle> removeList;

	// remove all non marked edges
	for (Edge_iterator itr = arr.edges_begin();itr != arr.edges_end();++itr)
	{
		if ((!itr->visited) && (!itr->isDegenerate))
			removeList.push_back(itr);		
	}

	for (typename std::list<Halfedge_handle>::iterator itr = removeList.begin();itr!=removeList.end();++itr)
		arr.remove_edge(*itr);

	degHandler.addDegenerateVerticesToArr();
	//final_stage_time = t3.elapsed();
	//std::cout << "degenerate stage : " <<degenerate_stage << std::endl; 
	//std::cout << "final stage : " <<final_stage_time << std::endl; 
    /*
	if (SHOW_STAGES)
	{
		global_graphics->clear();
		draw_arr(arr);
		global_graphics->display();
	}
    */

	/*constructOrientableLoops(arr);

	if (SHOW_STAGES)
	{
		global_graphics->clear();
		draw_arr(arr);
		global_graphics->display();
	}
	Polygon_2 reverse_pgn1 = transform(Kernel::Aff_transformation_2(CGAL::Rotation(), 0, -1), pgn1);
	nestedLoopsFilter(arr,reverse_pgn1,pgn2);
	degHandler.addDegenerateVerticesToArr();
	if (SHOW_STAGES)
	{
		global_graphics->clear();	
		draw_arr(arr);
		global_graphics->display();
	}
	*/
	//collisionDetectionFilter(arr);

	//global_graphics->clear();
	
	//draw_arr(arr);
	//global_graphics->display();
	//for (Segments_list::const_iterator itr = reduced_conv.begin();itr != reduced_conv.end();++itr)
	//{
	//	Segment_2 temp = *itr;
	//	Point_2 p_source = temp.source();
	//	Point_2 p_end = temp.target();
	//	double x1 = CGAL::to_double(p_source.x());
	//	double y1 = CGAL::to_double(p_source.y());
	//	double x2 = CGAL::to_double(p_end.x());
	//	double y2 = CGAL::to_double(p_end.y());
	////	CGAL::insert(arr,temp);
	//}
	delete _aabb_collision_detector;
	return (sum_holes);
  }

  void markOutsideLoop(Arrangement_history_2& arr){
	  Face_iterator ub_face = arr.unbounded_face();
	  Hole_iterator holes_itr = ub_face->holes_begin ();
	  Ccb_halfedge_circulator circ_start = *holes_itr;
	  Ccb_halfedge_circulator circ = circ_start;
	  do{
		  setEdgeVisited(*circ,true,0);
		  ++circ;
	  }while(circ!=circ_start);
	  
  }

  void markOutsideLoop(Arrangement_history_2& arr, Polygon_2& out_bound){
	  Face_iterator ub_face = arr.unbounded_face();
	  Hole_iterator holes_itr = ub_face->holes_begin ();
	  Ccb_halfedge_circulator circ_start = *holes_itr;
	  Ccb_halfedge_circulator circ = circ_start;
	  do{
		  setEdgeVisited(*circ,true,0);
		  out_bound.push_back (circ->source()->point());
		  --circ;
	  }while(circ!=circ_start);

  }
  /*
  void handleFace(Arrangement_history_2& arr,Face_handle itr,const Polygon_2& reverse_pgn1,const Polygon_2& pgn2) 
  {
	  if (itr->holes_begin() != itr->holes_end())
		  return;

	  Ccb_halfedge_circulator start = itr->outer_ccb();
	  Ccb_halfedge_circulator circ = start;
	  

	  // orientation check
	  do{
		  if (!checkTripNotSameDirWithSegment(arr,circ))
			  return;
		  ++circ;
	  }while (circ!=start);
	  //checkTripNotSameDirWithSegment

	  // collision detection check.
	  bool coll_detect = checkCollisionDetection(arr,start,reverse_pgn1,pgn2);
	  if (coll_detect)
		  return;

	  // mark as hole
	  circ = start;
	 // Polygon_2 pgn_hole;
	  do{
		  setEdgeVisited(*circ,true,0);
		 // pgn_hole.push_back (circ->source()->point());
		  --circ;
	  }while (circ!=start);

  }*/

  template <class OutputIterator>
  void handleFace(Arrangement_history_2& arr,Face_handle itr,const Polygon_2& reverse_pgn1,const Polygon_2& pgn2,OutputIterator holes) 
  {
	  if (itr->holes_begin() != itr->holes_end())
		  return;

	  Ccb_halfedge_circulator start = itr->outer_ccb();
	  Ccb_halfedge_circulator circ = start;


	  // orientation check
	  do{
		  if (!checkTripNotSameDirWithSegment(arr,circ))
			  return;
		  ++circ;
	  }while (circ!=start);
	  //checkTripNotSameDirWithSegment

	  // collision detection check.
	  bool coll_detect = checkCollisionDetection(arr,start,reverse_pgn1,pgn2);
	  if (coll_detect)
		  return;

	  // mark as hole
	  circ = start;
	  Polygon_2 pgn_hole;
	  do{
		  setEdgeVisited(*circ,true,0);
		  pgn_hole.push_back (circ->source()->point());
		  --circ;
	  }while (circ!=start);
//#pragma omp critical
//	  {
		  *holes = pgn_hole;
		  ++holes;
//	  }
	  
  }

  /*!
   * Compute the boundery Minkowski sum of two simple polygons.
   * The result is represented as
   * the outer boundary of the Minkowski sum (which is always a simple polygon).
   * \param pgn1 The first polygon.
   * \param pgn2 The second polygon.
   * \param sum_bound Output: A polygon respresenting the outer boundary
   *                          of the Minkowski sum.
   * 
   * \pre Both input polygons are simple.
   * \return A past-the-end iterator for the holes in the sum.
   */
  template <class OutputIterator>
  OutputIterator operator() (const Polygon_2& pgn1,
                             const Polygon_2& pgn2,
                             Polygon_2& sum_bound) const
  {
	CGAL_precondition (pgn1.is_simple());
    CGAL_precondition (pgn2.is_simple());
	
	Segments_list reduced_conv;
	buildReducedConvolution(pgn1,pgn2,reduced_conv);

/*	Arrangement_history_2 arr;
	buildArrangementFromConv(reduced_conv,arr);
*/
	const Minkowski_sum_by_convolution_lien_2* ptr = this;
	//	ConvSegMapper mapper( &arr,const_cast <Minkowski_sum_by_convolution_lien_2*>(ptr));*/
/*	DegenerateCassesManager degHandler(&arr,const_cast <Minkowski_sum_by_convolution_lien_2*>(ptr),const_cast <Polygon_2*>(&pgn1),const_cast <Polygon_2*>(&pgn2));
	degHandler.findDegenerateBorderVertices();
	degHandler.markDegenerateEdges();
	constructOrientableLoops(arr);

	if (SHOW_STAGES)
	{
		global_graphics->clear();
		draw_arr(arr);
		global_graphics->display();
	}
	Polygon_2 reverse_pgn1 = transform(Kernel::Aff_transformation_2(CGAL::Rotation(), 0, -1), pgn1);
	nestedLoopsFilter(arr,reverse_pgn1,pgn2);
	degHandler.addDegenerateVerticesToArr();
	if (SHOW_STAGES)
	{
		global_graphics->clear();	
		draw_arr(arr);
		global_graphics->display();
	}
	*/
  }


	void fillPolyDirs(const Polygon_2& pgn1,std::vector<Direction_2>& outVec) const 
	{	
		unsigned int n1 = pgn1.size();
		for (int i=0;i<(n1-1);++i)
		{
			outVec[i] = f_direction(f_vector(pgn1[i],pgn1[i+1]));

		}
		outVec[n1-1] = f_direction(f_vector(pgn1[n1-1],pgn1[0]));
	}
  	void buildReducedConvolution(const Polygon_2& pgn1, const Polygon_2& pgn2, Segments_list& reduced_conv) const
	{
		unsigned int n1 = pgn1.size();
		unsigned int n2 = pgn2.size();
		Vertex_circulator vert_p1,vert_p2,prev_p1,prev_p2,next_p1,next_p2;
		vert_p1 = prev_p1 = next_p1 = pgn1.vertices_circulator();
		vert_p2 = prev_p2 = next_p2 = pgn2.vertices_circulator();

		--prev_p1;
		--prev_p2;
		++next_p1;
		++next_p2;
		bool is_end_coincide;
		bool is_start_coincide;

		boost::unordered_map<std::pair<int,int>,Point_2> points_map;
		for (unsigned int i1 = 0;i1<n1;++i1)
		{
			for (unsigned int i2 = 0;i2<n2;++i2)
			{
				points_map[std::pair<int,int>(i1,i2)] = f_add(*vert_p1,Vector_2(Point_2(ORIGIN),*vert_p2));
				++vert_p2;
			}
			++vert_p1;
		}

		std::vector<Direction_2> p1_dirs(n1);
		std::vector<Direction_2> p2_dirs(n2);
		
		fillPolyDirs(pgn1,p1_dirs);
		fillPolyDirs(pgn2,p2_dirs);

		vert_p1 = pgn1.vertices_circulator();
		vert_p2 = pgn2.vertices_circulator();

		for (unsigned int i1 = 0;i1<n1;++i1)
		{
			for (unsigned int i2 = 0;i2<n2;++i2)
			{
				
				//Point_2 start_point =  f_add(*vert_p1,Vector_2(Point_2(ORIGIN),*vert_p2));
				Point_2 start_point =  points_map[std::pair<int,int>(i1,i2)];
				int prev_i1 = i1-1;
				if (prev_i1 == -1)
					prev_i1 = n1-1;
				int prev_i2 = i2-1;
				if (prev_i2 == -1)
					prev_i2 = n2-1;
				//if (!checkReflex(*prev_p1,*vert_p1,*next_p1) && checkSwept(*prev_p1,*vert_p1,*next_p1,*vert_p2,*next_p2,is_start_coincide,is_end_coincide)){
				if (!checkReflex(*prev_p1,*vert_p1,*next_p1) && checkSwept(p1_dirs[prev_i1],p1_dirs[i1],p2_dirs[i2],is_start_coincide,is_end_coincide)){
					//Point_2 end_point = f_add(start_point,Vector_2(*vert_p2,*next_p2));
					int cyc_ind = i2;
					++cyc_ind; 
					if (cyc_ind==n2)
						cyc_ind = 0;
					Point_2 end_point = points_map[std::pair<int,int>(i1,cyc_ind)];

			/*		double x1 = CGAL::to_double(start_point.x());
					double y1 = CGAL::to_double(start_point.y());
					double x2 = CGAL::to_double(end_point.x());
					double y2 = CGAL::to_double(end_point.y());
					double x3 = CGAL::to_double((*vert_p2).x());
					double y3 = CGAL::to_double((*vert_p2).y());
					double x4 = CGAL::to_double((*vert_p1).x());
					double y4 = CGAL::to_double((*vert_p1).y());*/
					CGAL::Comparison_result cres = f_compare_xy(start_point,end_point);
					Segment_2 conv_seg = Segment_2(Traits_2_A::Segment_2(start_point,end_point),cres);
					//Segment_2 conv_seg;
					//if (cres == CGAL::SMALLER)
					//{
					//	//Traits_2_B::X_monotone_curve_2 tempseg = Traits_2_B::X_monotone_curve_2(Traits_2_A::Segment_2(start_point,end_point),Segment_Data_Label(i1+i2*n1,i1+cyc_ind*n1));
					//	conv_seg = Segment_2(Traits_2_A::Segment_2(start_point,end_point),Segment_Data_Label(i1+i2*n1,i1+cyc_ind*n1,cres));
					//}
					//else
					//	conv_seg = Segment_2(Traits_2_A::Segment_2(start_point,end_point),Segment_Data_Label(i1+cyc_ind*n1,i1+i2*n1,cres));
					//CGAL::Comparison_result aaa = conv_seg.data();
					if (!is_end_coincide)
						reduced_conv.push_back(conv_seg);
				}
				//if (!checkReflex(*prev_p2,*vert_p2,*next_p2) && checkSwept(*prev_p2,*vert_p2,*next_p2,*vert_p1,*next_p1,is_start_coincide,is_end_coincide)){
				if (!checkReflex(*prev_p2,*vert_p2,*next_p2) && checkSwept(p2_dirs[prev_i2],p2_dirs[i2],p1_dirs[i1],is_start_coincide,is_end_coincide)){
					//Point_2 end_point = f_add(start_point,Vector_2(*vert_p1,*next_p1));
					int cyc_ind = i1;
					++cyc_ind; 
					if (cyc_ind==n1)
						cyc_ind = 0;
					Point_2 end_point = points_map[std::pair<int,int>(cyc_ind,i2)];


					//Segment_2 conv_seg = Segment_2(Traits_2_A::Segment_2(start_point,end_point),f_compare_xy(start_point,end_point));
					CGAL::Comparison_result cres = f_compare_xy(start_point,end_point);
					Segment_2 conv_seg = Segment_2(Traits_2_A::Segment_2(start_point,end_point),cres);
					/*Segment_2 conv_seg;
					if (cres == CGAL::SMALLER)
						conv_seg = Segment_2(Traits_2_A::Segment_2(start_point,end_point),Segment_Data_Label(i1+i2*n1,cyc_ind+i2*n1,cres));
					else
						conv_seg = Segment_2(Traits_2_A::Segment_2(start_point,end_point),Segment_Data_Label(cyc_ind+i2*n1,i1+i2*n1,cres));*/


					//Segment_2 conv_seg = Segment_2(start_point,end_point,f_compare_xy(start_point,end_point));
					if (!is_start_coincide)
						reduced_conv.push_back(conv_seg);
				}
				
				prev_p2 = vert_p2;
				vert_p2 = next_p2;
				++next_p2;
				
			}
			prev_p1 = vert_p1;
			vert_p1 = next_p1;
			++next_p1;
		}

	}

	// Increse a cyclic integer counter with limit lim.
	int cyclicInc(int i,int lim) const
	{
		i = i+1;
		if (i>=lim)
			i= 0;
		return i;
	}

	// Decrease a cyclic integer counter with limit lim.
	int cyclicDec(int i,int lim) const
	{
		i = i-1;
		if (i<0)
			i= lim-1;
		return i;
	}

	// Gets point corresponding to a state (i,j) if exists, creates this point if asked for first time.
	Point_2 addGetPoint(int i1,int i2,boost::unordered_map<std::pair<int,int>,Point_2>& points_map,const Polygon_2& pgn1, const Polygon_2& pgn2) const
	{
		Point_2 result;
		if (points_map.count(StatePair(i1,i2))==0)
		{
				result = f_add(pgn1[i1],Vector_2(Point_2(ORIGIN),pgn2[i2]));
				points_map[StatePair(i1,i2)] = result;
		}
		else
			result = points_map[StatePair(i1,i2)];
		return result;
	}

	// Builds the reduced convolution using the fiber grid approach. for each starting vertex, try to add out-going next states(two states).
	// If a visited vertex is reached then do not explore. This is a BFS like iteration beggining from each vertex in the first column of the 
	// fiber grid.
	void buildReducedConvolutionFiberGrid(const Polygon_2& pgn1, const Polygon_2& pgn2, Segments_list& reduced_conv) const
	{
		unsigned int n1 = pgn1.size();
		unsigned int n2 = pgn2.size();
		int i1 = 0;
		int i2 = 0;

		bool is_end_coincide;
		bool is_start_coincide;

		// Init the direcions of both polygons.
		std::vector<Direction_2> p1_dirs(n1);
		std::vector<Direction_2> p2_dirs(n2);
		
		fillPolyDirs(pgn1,p1_dirs);
		fillPolyDirs(pgn2,p2_dirs);

		
		boost::unordered_set<StatePair > visited_vertices_set;
		std::queue<StatePair > state_queue;
		boost::unordered_map<std::pair<int,int>,Point_2> points_map;

		// init the queue with vertices from the first column
		for (int i=n1-1;i>=0;--i){
			state_queue.push(StatePair(i,0));
		}

		
		while (state_queue.size() > 0){
			StatePair curr_state = state_queue.front();
			state_queue.pop();

			i1 = curr_state.first;
			i2 = curr_state.second;

			if (visited_vertices_set.count(curr_state) > 0)
				continue;

			visited_vertices_set.insert(curr_state);

			// add two outgoing edges:
			int next_p1 = cyclicInc(i1,n1);
			int next_p2 = cyclicInc(i2,n2);
			int prev_p1 = cyclicDec(i1,n1);
			int prev_p2 = cyclicDec(i2,n2);


			StatePair next_state_p1 = StatePair(next_p1,i2);
			StatePair next_state_p2 = StatePair(i1,next_p2);
			
			

			// add geometric entites of the transition from state (i,j) to (i+1,j) and (i,j+1), if they are in the reduced convolution.
			
		
			// Add an edge from Q
			if (checkSwept(p1_dirs[prev_p1],p1_dirs[i1],p2_dirs[i2],is_start_coincide,is_end_coincide) && !is_end_coincide)
			{
				state_queue.push(next_state_p2);
				if (!checkReflex(pgn1[prev_p1],pgn1[i1],pgn1[next_p1]))
				{
					//if (!is_end_coincide) // symbolic rotation
					//{
						Point_2 start_point = addGetPoint(i1,i2,points_map,pgn1,pgn2);
						Point_2 end_point = addGetPoint(i1,next_p2,points_map,pgn1,pgn2);

						CGAL::Comparison_result cres = f_compare_xy(start_point,end_point);
						//Segment_2 conv_seg = Segment_2(Traits_2_A::Segment_2(start_point,end_point),Segment_Data_Label(state(i1,i2),state(i1,next_p2),cres));
						Segment_2 conv_seg = Segment_2(typename Traits_2_A::Segment_2(start_point,end_point),Segment_Data_Label(state(i1,i2),state(i1,next_p2),cres,1));
						//Segment_2 conv_seg = Segment_2(Traits_2_A::Segment_2(start_point,end_point),cres);
						
						reduced_conv.push_back(conv_seg);
					//}
				}
			}
			
			// Add an edge from P
			if (checkSwept(p2_dirs[prev_p2],p2_dirs[i2],p1_dirs[i1],is_start_coincide,is_end_coincide) && !is_start_coincide)
			{
				state_queue.push(next_state_p1);
				if (!checkReflex(pgn2[prev_p2],pgn2[i2],pgn2[next_p2])	){
						//if () // symbolic rotation
						//{
							Point_2 start_point = addGetPoint(i1,i2,points_map,pgn1,pgn2);
							Point_2 end_point = addGetPoint(next_p1,i2,points_map,pgn1,pgn2);

							CGAL::Comparison_result cres = f_compare_xy(start_point,end_point);
							Segment_2 conv_seg = Segment_2(typename Traits_2_A::Segment_2(start_point,end_point),Segment_Data_Label(state(i1,i2),state(next_p1,i2),cres,0));
							//Segment_2 conv_seg = Segment_2(Traits_2_A::Segment_2(start_point,end_point),cres);
							reduced_conv.push_back(conv_seg);
						//}
				}
			}
		}

	}

private:
	SweepCollisionDetector<Kernel,Container_> collision_detector;
	AABBCollisionDetector<Kernel,Container_>* _aabb_collision_detector;
	/*
		Perform final stage of filtering : collision detection.
	*/
/*	void collisionDetectionFilter(Arrangement_history_2& arr) const
	{
		ICollisionDetector& collision_detector = new NaiveCollisionDetector();
		Face_iterator startFace = arr.unbounded_face();
		Hole_iterator hi;
		Ccb_halfedge_circulator perimiterFace;
		for (hi = startFace->holes_begin(); hi != startFace->holes_end(); ++hi) {
			 perimiterFace = *hi;
		}

		Halfedge_handle inside_face_edge = hi->twin();                                                                                                                                  
		Face_handle container_face = inside_face_edge->face();

	}
*/	
	/*
		Performs the stage 3 of filtering: removing nested loops which are not in the correct orientation.
	*/
	void nestedLoopsFilter(Arrangement_history_2& arr,const Polygon_2& pgn1,const Polygon_2& pgn2) const
	{
		Face_iterator startFace = arr.unbounded_face();
		Hole_iterator hi;
		Ccb_halfedge_circulator perimiterFace;
		for (hi = startFace->holes_begin(); hi != startFace->holes_end(); ++hi) {
			 perimiterFace = *hi;
		}

		//std::cout << "nested loop start\n";
		
		// now for each hole in main face we will determine it's orientation.
		nestedLoopsFilterRec(arr,perimiterFace,false,pgn1,pgn2);
	}

		/*
		For each loop we represent it by the halfedge which is the twin of the face loop. 
		ie this edge is part of the clockwise oriented halfedges loop outside the face.
	*/
	void nestedLoopsFilterRec(Arrangement_history_2& arr,Halfedge_handle& handle,bool inwards,const Polygon_2& pgn1,const Polygon_2& pgn2) const
	{
		std::list<Halfedge_handle> holesEdges;
		std::list<Halfedge_handle> after_removal_hole_edges;
		std::list<Face_handle> faces_list;

		Face_iterator startFace = arr.unbounded_face();
		Halfedge_handle inside_face_edge = handle->twin();                                                                                                                                  
		Face_handle container_face = inside_face_edge->face();

		Face_iterator face_itr=arr.faces_begin();
		for (;face_itr!=arr.faces_end();++face_itr)
		{
			if (!(face_itr->is_unbounded()) && face_itr!= container_face)
			{
				Halfedge_handle h_e = face_itr->outer_ccb()->twin();
				holesEdges.push_back(h_e);
			}
		}

		while(holesEdges.size() > 0) // remove loops from faces
		{
			/*int bla = holesEdges.size();
			std::cout << "in holes: " << bla << "\n";*/
			Halfedge_handle he = holesEdges.front();
		//	printHe(he);
			holesEdges.pop_front();

			if (!he->isDegenerate)
			{
				if (!removeAllNonConformingLoops(arr,he,inwards,holesEdges,pgn1,pgn2))
				{
					
					
					/*for (hi = startFace->holes_begin(); hi != startFace->holes_end(); ++hi) {
					 perimiterFace = *hi;
					}*/
					/*std::cout << "not removed";
					int bla = holesEdges.size();
					std::cout << "in holes: " << bla << "\n";*/
					after_removal_hole_edges.push_back(he);
				}
			}
			/*else
			{
				std::cout << "degen not removed";
				int bla = holesEdges.size();
				std::cout << "in holes: " << bla << "\n";
			}*/
		
		}
		Hole_iterator hi;
		//Ccb_halfedge_circulator perimiterFace;
		hi = startFace->holes_begin();
		container_face = (*hi)->twin()->face();
	/*	global_graphics->clear();
	
		draw_arr(arr);
		global_graphics->display();*/
//		Hole_iterator hi;
		Ccb_halfedge_circulator outside_face_itr;
		// push initial list of holes to list.
		for (hi =container_face->holes_begin();hi!= container_face->holes_end(); ++hi ) // remove degenrate cases which are false
		{
			outside_face_itr = *hi;
			Halfedge_handle h_e= outside_face_itr;//->twin();
			holesEdges.push_back(h_e);
		}

		while(holesEdges.size() > 0)
		{
		/*	int bla = holesEdges.size();
			std::cout << "in holes: " << bla << "\n";*/
			Halfedge_handle he = holesEdges.front();
		//	printHe(he);
			holesEdges.pop_front();

			if (!he->isDegenerate)
			{
				if (!removeAllNonConformingLoops(arr,he,inwards,holesEdges,pgn1,pgn2))
				{
					Hole_iterator hi;
					//Ccb_halfedge_circulator perimiterFace;
			/*		hi = startFace->holes_begin();
					container_face = (*hi)->twin()->face();*/
				/*	std::cout << "not removed";
					int bla = holesEdges.size();
					std::cout << "in holes: " << bla << "\n";*/
					after_removal_hole_edges.push_back(he);
				}
			}
			/*else
			{
				std::cout << "degen not removed";
				int bla = holesEdges.size();
				std::cout << "in holes: " << bla << "\n";
			}*/
		
		}

	//	std::list<Halfedge_handle> semi_holes;
		//findSemiHoles(arr,handle,semi_holes);
	//	holesEdges.insert(holesEdges.begin(),semi_holes.begin(),semi_holes.end());

		// while we still have holes , leave only loops orientable with current direction.
		
	/*	std::cout << "number of hole edges " << holesEdges.size() << "\n";
		std::cout << arr.number_of_faces() << " faces:" << std::endl;
		

		std::cout << arr.number_of_faces() << " faces:" << std::endl;*/
/*		global_graphics->clear();
	
		draw_arr(arr);
		global_graphics->display();
	*/
		
/*		for (std::list<Halfedge_handle>::iterator itr = after_removal_hole_edges.begin();itr != after_removal_hole_edges.end() ; ++itr)
		{
			nestedLoopsFilterRec(arr,*itr, !inwards,pgn1,pgn2);
		}
		*/
		/*container_face = inside_face_edge->face();
		// now all the holes are in correct orientation, recursivaly call this func:
		for (hi =container_face->holes_begin();hi!= container_face->holes_end(); ++hi )
		{
			outside_face_itr = *hi;
			nestedLoopsFilterRec(arr,outside_face_itr, !inwards,pgn1,pgn2);
		//	Halfedge_handle h_e= outside_face_itr.twin();
		//	holesEdges.push_back(h_e);
		}*/
	}



//	/*
//		For each loop we represent it by the halfedge which is the twin of the face loop. 
//		ie this edge is part of the clockwise oriented halfedges loop outside the face.
//	*/
//	void nestedLoopsFilterRec(Arrangement_history_2& arr,Halfedge_handle& handle,bool inwards,const Polygon_2& pgn1,const Polygon_2& pgn2) const
//	{
//		std::list<Halfedge_handle> holesEdges;
//		std::list<Halfedge_handle> after_removal_hole_edges;
//	//	std::list<Face_handle> faces_list;
//
//
//		Halfedge_handle inside_face_edge = handle->twin();                                                                                                                                  
//		Face_handle container_face = inside_face_edge->face();
//
//	
//
//		Hole_iterator hi;
//		Ccb_halfedge_circulator outside_face_itr;
//		// push initial list of holes to list.
//		for (hi =container_face->holes_begin();hi!= container_face->holes_end(); ++hi )
//		{
//			outside_face_itr = *hi;
//			Halfedge_handle h_e= outside_face_itr;//->twin();
//			holesEdges.push_back(h_e);
//		}
//
//	//	std::list<Halfedge_handle> semi_holes;
//		//findSemiHoles(arr,handle,semi_holes);
//	//	holesEdges.insert(holesEdges.begin(),semi_holes.begin(),semi_holes.end());
//
//		// while we still have holes , leave only loops orientable with current direction.
//		
//		std::cout << "number of hole edges " << holesEdges.size() << "\n";
//		std::cout << arr.number_of_faces() << " faces:" << std::endl;
//		while(holesEdges.size() > 0)
//		{
//		/*	int bla = holesEdges.size();
//			std::cout << "in holes: " << bla << "\n";*/
//			Halfedge_handle he = holesEdges.front();
//		//	printHe(he);
//			holesEdges.pop_front();
//
//			if (!he->isDegenerate)
//			{
//				if (!removeAllNonConformingLoops(arr,he,inwards,holesEdges,pgn1,pgn2))
//				{
//					std::cout << "not removed";
//					int bla = holesEdges.size();
//					std::cout << "in holes: " << bla << "\n";
//					after_removal_hole_edges.push_back(he);
//				}
//			}
//			else
//			{
//				std::cout << "degen not removed";
//				int bla = holesEdges.size();
//				std::cout << "in holes: " << bla << "\n";
//			}
//		
//		}
//
//		std::cout << arr.number_of_faces() << " faces:" << std::endl;
///*		global_graphics->clear();
//	
//		draw_arr(arr);
//		global_graphics->display();
//	*/
//		
//		for (std::list<Halfedge_handle>::iterator itr = after_removal_hole_edges.begin();itr != after_removal_hole_edges.end() ; ++itr)
//		{
//			nestedLoopsFilterRec(arr,*itr, !inwards,pgn1,pgn2);
//		}
//
//		/*container_face = inside_face_edge->face();
//		// now all the holes are in correct orientation, recursivaly call this func:
//		for (hi =container_face->holes_begin();hi!= container_face->holes_end(); ++hi )
//		{
//			outside_face_itr = *hi;
//			nestedLoopsFilterRec(arr,outside_face_itr, !inwards,pgn1,pgn2);
//		//	Halfedge_handle h_e= outside_face_itr.twin();
//		//	holesEdges.push_back(h_e);
//		}*/
//	}

	// find all faces who has handle's face sorrounding them, but are not holes.
	void findSemiHoles(Arrangement_history_2& arr,Halfedge_handle& handle,std::list<Halfedge_handle>& semi_holes) const
	{
		Faces_set faces_in_face;
		Face_handle outside_face = handle->face();
		Ccb_halfedge_circulator circ = handle->twin()->ccb();
		Ccb_halfedge_circulator curr = circ;
		do{
			// check if the twin edge is the regular outside face or a different one. if so we are in half island.
			Face_handle curr_sec_face = curr->twin()->face();
			if (curr_sec_face != outside_face)
			{
				typename Faces_set::iterator itr = faces_in_face.find(curr_sec_face);
			
				if (itr == faces_in_face.end())
				{
					faces_in_face.insert(curr_sec_face);
				}
			}
		}while(++curr!=circ);
		for (typename Faces_set::iterator itr = faces_in_face.begin();itr!= faces_in_face.end();++itr)
		{
			Face_handle h = *itr;
			semi_holes.push_back((h->outer_ccb()->twin()));
		}
		/*Halfedge_handle inside_face_edge = handle->twin();                                                                                                                                  
		Face_handle outside_face = handle->face();
		
		Face_handle inside_face = inside_face_edge->face();
		Faces_set visited_faces;*/
		
	}

	bool removeAllNonConformingLoops(Arrangement_history_2& arr,Halfedge_handle& handle,bool inwards,std::list<Halfedge_handle>& holesEdges,const Polygon_2& pgn1,const Polygon_2& pgn2) const
	{
		// check if loop is in the right direction ((a || b) && !(a && b))
		bool a =!checkTripSameDirWithSegment(arr, handle);
		bool b = inwards;
		// nor a,b
		bool conforming_loop =  !((a || b) && !(a && b));

		if (!conforming_loop)
		{
			//std::cout << "start remove\n";

	/*		Halfedge_handle inside_face_edge = handle->twin();
			Face_handle container_face = inside_face_edge->face();
			Hole_iterator hi;
			Ccb_halfedge_circulator outside_face_itr;
			// push initial list of holes to list.
			for (hi =container_face->holes_begin();hi!= container_face->holes_end(); ++hi )
			{
				outside_face_itr = *hi;
				Halfedge_handle h_e= outside_face_itr->twin();
				holesEdges.push_back(h_e);
			}
			*/
			removeFaceLoop(arr,handle->twin());
			//std::cout << "end remove\n";
			return true;
		}

		// Check collision detection criterion
		bool coll_detect = checkCollisionDetection(arr,handle,pgn1,pgn2);
		if (coll_detect)
		{
			removeFaceLoop(arr,handle->twin());
			//std::cout << "end remove\n";
			return true;
		}
		return false;

	}

	Polygon_2 revPoly(const Polygon_2& input) const
	{
		Polygon_2 out;
		typename Polygon_2::Vertex_iterator itr = input.vertices_begin();
		for (;itr!= input.vertices_end();++itr)
		{
			out.push_back(Point_2(-itr->x(),-itr->y()));
		}
		if (out.orientation()==CGAL::CLOCKWISE)
			out.reverse_orientation();

/*
		QColor c1(0,255,0);
		QColor c2(0,0,255);
        */
	
		
		return out;

	}

	/*void buildCollisionDetector(const Polygon_2& pgn1,const Polygon_2& pgn2)
	{
		q = revPoly(pgn2);
		
	}*/

	/*ICollisionDetector<Kernel,Container_>* getColDetect() const
	{
		return (ICollisionDetector<Kernel,Container_>*)&collision_detector;			
	}*/

	AABBCollisionDetector<Kernel,Container_>* getColDetect() const
	{
		return _aabb_collision_detector;			
	}

	/*
	This version assumes poly1 is reflected through origin. (as called from nested loops filter)
	*/
	bool checkCollisionDetection(Arrangement_history_2& arr,Halfedge_handle& handle,const Polygon_2& pgn1,const Polygon_2& pgn2) const
	{
		//ICollisionDetector<Kernel,Container_>& collision_detector = SweepCollisionDetector<Kernel,Container_>();
		//ICollisionDetector<Kernel,Container_>* collision_detector = getColDetect();
		AABBCollisionDetector<Kernel,Container_>* collision_detector = getColDetect();
		//ICollisionDetector<Kernel,Container_>& collision_detector = NaiveCollisionDetector<Kernel,Container_>();
		//handle = handle->twin();
		
		/*Point_2 p = handle->source()->point();
		Point_2 p2 = handle->target()->point();
		Point_2 mid_point = CGAL::midpoint(p,p2);*/
		Point_2 mid_point = findInsidePoint(arr,handle);
		//Polygon_2 r_pgn1 = CGAL::transform(Kernel::Aff_transformation_2(CGAL::Rotation(),0,-1), pgn1);
		Polygon_2 t_pgn1 = transform(typename Kernel::Aff_transformation_2(CGAL::Translation(), Vector_2(CGAL::ORIGIN,mid_point)), pgn1);
		collision_detector->setTranslationPoint(mid_point);		
		//QColor c1(0,255,0);
		//QColor c2(0,0,255);
	
		//global_graphics->draw_polygon(t_pgn1,c1);
		//global_graphics->draw_polygon(pgn2,c2);
		//global_graphics->display();
		//global_graphics->clear();
		return collision_detector->checkCollision(t_pgn1,pgn2);
	}

	Point_2 findInsidePoint(Arrangement_history_2& arr,Halfedge_handle& handle) const
	{
		
		Ccb_halfedge_circulator currHandle = handle->ccb();
		Ccb_halfedge_circulator nextHandle = currHandle;
		++nextHandle;
		//while ((currHandle->curve().is_vertical()) || currHandle->direction() != nextHandle->direction())
		while (currHandle->direction() != nextHandle->direction())
		{
//			printHe(currHandle);
//			printHe(nextHandle);
			++currHandle;
			++nextHandle;
			if (checkReflex(currHandle->source()->point(),currHandle->target()->point(),nextHandle->target()->point()))
				break;
		}

		Point_2 p = currHandle->source()->point();
		Point_2 p2 = currHandle->target()->point();
		Point_2 work_point = p2;

		Ccb_halfedge_circulator best_edge= handle;
		bool has_some_point = false;
		Ccb_halfedge_circulator circ = nextHandle;
		Ccb_halfedge_circulator end = handle;

		bool shoot_upwards = (currHandle->direction() == ARR_LEFT_TO_RIGHT);
		if (nextHandle->curve().is_vertical())
			work_point = CGAL::midpoint(p,p2);
		if (currHandle->curve().is_vertical())
		{
			p = nextHandle->source()->point();
		    p2 = nextHandle->target()->point();
			work_point = CGAL::midpoint(p,p2);
			++best_edge;
			++circ;
			++end;
		}
		//Point_2 mid_point = CGAL::midpoint(p,p2);
		
		++circ;
		while (circ!=end)
		{ 
			Base_Segment_2 circ_curve = circ->curve();
			
			if (f_compare_x(work_point,circ_curve.min()) != f_compare_x(work_point,circ_curve.max()))
			{ // we have an edge with same x range as endpoint of 
				bool above_first = (f_compare_y_at_x(work_point,circ_curve) == SMALLER);
				if (has_some_point)
				{
					bool under_best;
					Base_Segment_2 best_edge_curve = best_edge->curve();
					if (f_compare_x(best_edge_curve.min(),circ_curve.min()) != f_compare_x(best_edge_curve.max(),circ_curve.min())) 
					{
						under_best = f_compare_y_at_x(circ_curve.min(),best_edge_curve) == SMALLER;
						
					}
					else
					{
						under_best = f_compare_y_at_x(best_edge_curve.min(),circ_curve) != SMALLER;
					}
					if ( (shoot_upwards && above_first && under_best) ||(!shoot_upwards && !above_first && !under_best))
						best_edge = circ;
				}
				else
				{
					has_some_point = true;
					best_edge = circ;
				}
			}
			++circ;
		}

		if (best_edge->curve().is_vertical())
		{
			Base_Segment_2 best_edge_curve = best_edge->curve();
			typename Kernel::FT x0 =  f_compute_x(work_point);
			typename Kernel::FT y_point = f_compute_y(work_point); 
			
			if (shoot_upwards){
				typename Kernel::FT y_best = f_compute_y(best_edge_curve.min());
				typename Kernel::FT y = (y_best-y_point)/2 + y_point;
				return Point_2(x0,y);
			}
			else
			{
				typename Kernel::FT y_best = f_compute_y(best_edge_curve.min());
				typename Kernel::FT y = (y_point - y_best)/2 + y_best;
				return Point_2(x0,y);
			}
			return work_point;
		}
		Base_Segment_2 best_edge_curve = best_edge->curve();
		typename Kernel::FT x0 =  f_compute_x(work_point);
		typename Kernel::FT x1 =  f_compute_x(best_edge_curve.min());
		typename Kernel::FT x2 =  f_compute_x(best_edge_curve.max());
		typename Kernel::FT alpha = (x0-x2)/(x1-x2);

		typename Kernel::FT y_best = alpha*f_compute_y(best_edge_curve.min())+ (1-alpha)*f_compute_y(best_edge_curve.max());
		typename Kernel::FT y_point = f_compute_y(work_point); 
		typename Kernel::FT y = (y_best-y_point)/2 + y_point;


		return Point_2(x0,y);
		//return work_point;
	}

	/*
		This version reflects poly 1.
	*/
	bool checkCollisionDetection(Arrangement_history_2& arr,Point_2& point,const Polygon_2& pgn1,const Polygon_2& pgn2) const
	{
		//ICollisionDetector<Kernel,Container_>& collision_detector = NaiveCollisionDetector<Kernel,Container_>();
		//ICollisionDetector<Kernel,Container_>* collision_detector = getColDetect();
		AABBCollisionDetector<Kernel,Container_>* collision_detector = getColDetect();
		//handle = handle->twin();
		Point_2 p = point;
		//Polygon_2 r_pgn1 = CGAL::transform(Kernel::Aff_transformation_2(CGAL::Rotation(),0,-1), pgn1);
		Polygon_2 r_pgn1 = revPoly(pgn1);
		Polygon_2 t_pgn1 = transform(typename Kernel::Aff_transformation_2(CGAL::Translation(), Vector_2(CGAL::ORIGIN,p)), r_pgn1);
		collision_detector->setTranslationPoint(p);		
	/*	QColor c1(0,255,0);
		QColor c2(0,0,255);
	
		global_graphics->draw_polygon(t_pgn1,c1);
		global_graphics->draw_polygon(pgn2,c2);
		global_graphics->display();
		global_graphics->clear(); */
		return collision_detector->checkCollision(t_pgn1,pgn2);
	}

	bool checkSegmentCollisionDetection(Arrangement_history_2& arr,Segment_2& seg,const Polygon_2& pgn1,const Polygon_2& pgn2) const
	{
		Point_2 mid_point= CGAL::midpoint(seg.source(),seg.target());
		return checkCollisionDetection(arr,mid_point,pgn1,pgn2);
	}

	void removeFaceLoop(Arrangement_history_2& arr,Halfedge_handle& handle) const
	{
		std::list<Halfedge_handle> remove_list;
		Ccb_halfedge_circulator circ = handle->ccb();
		remove_list.push_front(handle);
		++circ;
		//Halfedge_iterator startEdge = circ.begin();
		while ( circ != handle)
		{
			remove_list.push_front(circ);
			++circ;
		}

		for (typename std::list<Halfedge_handle>::iterator itr = remove_list.begin();itr != remove_list.end();++itr)
		{
			arr.remove_edge(*itr);
		}
	}
	void buildArrangementFromConv(const Segments_list& reduced_conv,Arrangement_history_2& arr) const
	{
		CGAL_precondition (arr.is_empty());
		CGAL::insert (arr, reduced_conv.begin(),reduced_conv.end());
	}

	void constructOrientableLoops(Arrangement_history_2& arr) const 
	{
	
        /*
		if (SHOW_STAGES)
		{
			draw_arr(arr);
			global_graphics->display();
		}
        */

		if (WRITE_ARR)
		{
            std::ofstream file("arr.txt");
			file << arr;
		}
		const Minkowski_sum_by_convolution_lien_2* ptr = this;
		ConvSegMapper mapper( &arr,const_cast <Minkowski_sum_by_convolution_lien_2*>(ptr));
		TraversalManager manager(&mapper);
		manager.traverseLoops();

	/* OLD Version 10/7/11
		Edge_iterator  itr;
		Edges_set edges_set;
		for (itr = arr.edges_begin();itr!=arr.edges_end();++itr){
			setEdgeVisited(*itr,false,-1);
			edges_set.insert(itr);

		}


		// trace orientable loops:
		
		while (edges_set.size() != 0)
		{			
			traceOrientableLoops( arr,edges_set);
		}
		*/
		/*std::cout << "arr" << std::endl;
		for (itr = arr.edges_begin();itr!=arr.edges_end();++itr){
			
			printHe(itr);
		}*/
		//typename iterator_traits<InputIterator>::value_type seg;
		
	}

		
	bool checkDegenerateEdgeOppositeSegments(Arrangement_history_2& arr,Halfedge_handle he) const
	{
		Originating_curve_iterator segment_itr;// = arr.originating_curves_begin ( *he);
		//bool found = false;
		//segment = NULL;

		std::list<Direction_2> segments_dir_list;
	
		for (segment_itr = arr.originating_curves_begin(he);segment_itr != arr.originating_curves_end(he);++segment_itr)
		{
			Segment_2 segment = *segment_itr;
			/*std::cout << "checktrip: " << std::endl;
			printSegment(segment);
			printHe(he);*/
			Direction_2 seg_dir = f_direction(f_vector(segment.source(),segment.target()));
			segments_dir_list.push_back(seg_dir);
			//Direction_2 start_he_dir = f_direction((f_vector(he->source()->point(),he->target()->point())));
				
			
		}

		//Direction_2 start_dir = segments_dir_list.begin();
		segments_dir_list.sort();
		typename std::list<Direction_2>::iterator end = unique(segments_dir_list.begin(),segments_dir_list.end());
		int i =distance(segments_dir_list.begin(),end);
		return i>1;
	
	}

	bool checkDegenarateVertexIsIntersectionOfThreeSegments(Arrangement_history_2& arr,Vertex_handle vh)
	{
		Halfedge_around_vertex_circulator itr = vh->incident_halfedges();
		Halfedge_around_vertex_circulator start = itr;
		int count_degree =0;
		do{
			++count_degree;
		}
		while(++itr!=start);

		if (count_degree <=2)
			return false;
		
		Originating_curve_iterator segment_itr;
		std::list<Segment_2*> orig_segments_list;
		std::list<Direction_2> segments_dir_list;
		//std::list<> originating_curves_list;
	
		// Handle the standard case where we have two intersecting convolution segments.
		if (count_degree == 4){
			do{
				for (segment_itr = arr.originating_curves_begin(itr);segment_itr != arr.originating_curves_end(itr);++segment_itr){
					Segment_2 segment = *segment_itr;
					orig_segments_list.push_back(&(*segment_itr));
				}
			} while(++itr!=start);
			orig_segments_list.sort();
			typename std::list<Segment_2*>::iterator end = unique(orig_segments_list.begin(),orig_segments_list.end());
			int i =distance(orig_segments_list.begin(),end);
			if (i==2) // this is two curves crossing case. 
				return false;
		}


	/* DEBUG
	if ( != 4)
		{
			std::cout << "Degree not 4 : " << count_degree << std::endl;
		}
*/
		
		do{
			
			for (segment_itr = arr.originating_curves_begin(itr);segment_itr != arr.originating_curves_end(itr);++segment_itr)
			{
				Segment_2 segment = *segment_itr;
				//printSegment(*segment_itr);
				//orig_segments_list.push_back(&(*segment_itr));
				Direction_2 seg_dir = f_direction(f_vector(segment.source(),segment.target()));
				segments_dir_list.push_back(seg_dir);
				
			}
			
		}
		while(++itr!=start);
		
		

		segments_dir_list.sort();
		typename std::list<Direction_2>::iterator end = unique(segments_dir_list.begin(),segments_dir_list.end());
		int i =distance(segments_dir_list.begin(),end);
		return i>2;
		//orig_segments_list.sort();


		/*std::list<Segment_2*>::iterator end = unique(orig_segments_list.begin(),orig_segments_list.end());
		int i =distance(orig_segments_list.begin(),end);

		return (distance(orig_segments_list.begin(),end) > 2);*/
	}

	// Gets the he that agrees in direction with the convolution segment.
	Halfedge_handle getDirAgreeingHalfedge(Arrangement_history_2& arr,const Halfedge_handle& he) const
	{
		Halfedge_handle curr_halfedge = he;
		if (!checkTripSameDirWithSegment(arr,he))
			curr_halfedge = curr_halfedge->twin();
		return curr_halfedge;
	}

	// Gets list of incoming and outgoing edges(as defined by directions of segments in the convolution) from the vertex.
	void getEdgesFromVertex(Arrangement_history_2& arr,Vertex_handle v_src,std::list<Halfedge_handle>& inList,std::list<Halfedge_handle>& outList) const
	{
		outList.clear();
		inList.clear();
		Halfedge_around_vertex_circulator itr = v_src->incident_halfedges();
		Halfedge_around_vertex_circulator start = itr;
		do{
			Halfedge_handle curr_edge = getDirAgreeingHalfedge(arr,itr);
			if ((curr_edge->source()) == v_src)
			{
				outList.push_back(curr_edge);
			}else
			{
				inList.push_back(curr_edge);
			}

		}
		while(++itr!=start);
	}

	// Returns the direction of a half edge
	Direction_2 getHalfedgeDir(const Halfedge_handle& he) const
	{
		Direction_2 dir = f_direction((f_vector(he->source()->point(),he->target()->point())));
		return dir;
	}

	double getSignedAngle(const Halfedge_handle& h_enter,const Halfedge_handle& h_exit) const
	{
		Direction_2 dir_enter = getHalfedgeDir(h_enter);
		Direction_2 dir_exit = getHalfedgeDir(h_exit);
		Vector_2 vec_enter =  dir_enter.vector();
		Vector_2 vec_exit =  dir_exit.vector();
		Point_2 org(CGAL::ORIGIN);
		Vector_2 origin_vec(org,org);
		Orientation sign_or = f_orientation(vec_enter,vec_exit);
		float sign = 0.f;
		if (sign_or == CGAL::LEFT_TURN)
			sign = 1;
		else 
			if (sign_or == CGAL::RIGHT_TURN)
				sign = -1;
			else 
				sign = 0;
		
		double prod = CGAL::to_double(vec_enter*vec_exit);
		
		/*if ((CGAL::is_zero(vec_enter.squared_length()) == CGAL::Tag_true) || (CGAL::is_zero(vec_exit.squared_length()) == CGAL::Tag_true)) 
			return 0;*/
		if (f_equal(vec_enter,origin_vec) || f_equal(vec_exit,origin_vec))
			return 0;
		
		double len1 = sqrt(CGAL::to_double(vec_enter.squared_length()));
		double len2 = sqrt(CGAL::to_double(vec_exit.squared_length()));
		//if (abs(len1 -eps()) == 0 || len2 == 0)
		double p = prod/(len1*len2);
		p = min((double)(1),p);
		p = max((double)(-1),p);
		double ang = acos(p);
		return sign*ang;
	}

    /*
	void traceOrientableLoopNew(Arrangement_history_2& arr,Halfedge_handle& start_halfedge,bool is_outer_loop,int loop_counter,Edges_set& edges_set) const
	{
		std::list<Halfedge_handle> temp_segments;
		Halfedge_handle trace_start = getDirAgreeingHalfedge(arr,start_edge);
		Halfedge_handle curr_halfedge = trace_start;
		double sum_angles = 0;
		bool done = false;
		bool improving = false;
		while (!done){
			bool has_loop_closed = false;

			Halfedge_handle next_halfedge = traverseNextHalfedge(Arrangement_history_2& arr,curr_halfedge,has_next_edge,loop_counter,has_loop_closed,loop_start_handle,loop_end_handle);
			
			if (has_loop_closed)
			{
				double temp_sum = sum + getSignedAngle(curr_edge,loop_closed_handle);
				// check if loop is in correct orientation.
				if ((temp_sum > 0 && is_outer_loop) || (temp_sum < 0 && !is_outer_loop))
				{
					
					improving = true;
					sum =0;
					// loop was closed and nowhere to continue
					if (!has_next_edge)
					{
						// a loop closes with part of the edges. we have to remove edges till next_halfedge from arrangment,
						// and clear the list.
						removeHalfEdgesArrPart(arr,edges_set,temp_segments,loop_start_handle);
						//removeHalfEdgesArrPart(arr,edges_set,temp_segments,curr_halfedge);
						nextLoop(temp_segments, arr,edges_set, curr_halfedge,start_halfedge,done,loop_counter);

					}
				}
			}
			else
			{

			}
		}
	}
    */

	Halfedge_handle traverseNextHalfedge()
	{

	}



	void traceOrientableLoops(Arrangement_history_2& arr,Edges_set& edges_set) const
	{
		//#define SHOW_LOOPS_CONST
		int loop_counter = 0;
		std::list<Halfedge_iterator> temp_segments;

		double angles_sum = 0;

		// get an edge that is surly on the outside border. we will start from the external loop.
		Face_iterator startFace = arr.unbounded_face();
		Halfedge_iterator perimiterFace = *(startFace -> holes_begin());
		//Ccb_halfedge_circulator perimiterFace = *(startFace -> holes_begin());
		/*for (hi = startFace->holes_begin(); hi != startFace->holes_end(); ++hi) {
			 perimiterFace = *hi;
		}*/
		// change we trace outer loop first
		Halfedge_iterator curr_halfedge =  perimiterFace;//(*edges_set.begin());
		
		//Originating_curve_iterator first_segment_itr = arr.originating_curves_begin ( *curr_halfedge);
		//Segment_2 first_segment = *first_segment_itr;
		//Direction_2 start_seg_dir = f_direction(f_vector(first_segment.source(),first_segment.target()));
		//Direction_2 start_he_dir = f_direction((f_vector(curr_seg.source(),curr_seg.target()));

		if (!checkTripSameDirWithSegment(arr,curr_halfedge))
			curr_halfedge = curr_halfedge->twin();

		
		Halfedge_iterator start_halfedge = curr_halfedge;
		temp_segments.push_back(curr_halfedge);
		setEdgeVisited(*curr_halfedge,true,loop_counter);
		Halfedge_handle close_loop_handle;
		Halfedge_handle before_close_loop_handle;
		bool isOrientable = true;
		bool close_loop_found = false;
		bool tracing_holes = false;
		//int loop_counter;

		while(isOrientable )
		{
			// after first loop we are starting to trace the holes.
			if (loop_counter >0)
				tracing_holes = true;
			printHe(curr_halfedge);
			/*for (Originating_curve_iterator itr = arr.originating_curves_begin ( curr_halfedge);itr!=arr.originating_curves_end( curr_halfedge);++itr ){
				Segment_2 segment = *itr;

			}*/
			/////////////////////////DEBUG 
			//Segment_2 temp = curr_halfedge->curve();
			Point_2 p_source = curr_halfedge->source()->point();
			Point_2 p_end = curr_halfedge->target()->point();
			double x1 = CGAL::to_double(p_source.x());
			double y1 = CGAL::to_double(p_source.y());
			double x2 = CGAL::to_double(p_end.x());
			double y2 = CGAL::to_double(p_end.y());
			int size_edges_set = edges_set.size();

			/////////////////////////DEBUG 
			Vertex_iterator v_target = curr_halfedge->target();
			 p_source = v_target->point();
			double x = CGAL::to_double(p_source.x());
			double y = CGAL::to_double(p_source.y());
			Halfedge_around_vertex_circulator itr = v_target->incident_halfedges();
			bool next_edge_found;
			bool temp_close_loop_found =false;
			
			Halfedge_iterator next_halfedge = getLargestExitingClockwiseEdge(arr,curr_halfedge,*v_target,next_edge_found,loop_counter,temp_close_loop_found,close_loop_handle);
			if (!next_halfedge->visited && next_edge_found)
				setEdgeVisited(*next_halfedge,true,loop_counter);
			// if any ending of loop was found remember it.
			if (temp_close_loop_found)
			{
				before_close_loop_handle = curr_halfedge;
				close_loop_found = true;
			}
		//	printHe(next_halfedge);
			//Halfedge_handle next_halfedge = getLargestExitingClockwiseEdge();
			
			// Check if we are stuck or we have met a visited edge which closes the loop:
			if (!next_edge_found ){
				Halfedge_handle temp_twin = next_halfedge->twin();
				if ((next_halfedge->loopNumber == loop_counter) && (next_halfedge != curr_halfedge)  && (temp_twin != curr_halfedge)){ 
					#ifdef SHOW_LOOPS_CONST
					//global_graphics->clear();

					for (std::list<Halfedge_iterator>::iterator itr = temp_segments.begin();itr != temp_segments.end(); ++itr)
					{
						//global_graphics->draw_edge<Kernel>((*itr)->curve(),QColor(0,255,0),true);
					}
					//global_graphics->display();
					#endif
					// a loop closes with part of the edges. we have to remove edges till next_halfedge from arrangment,
					// and clear the list.
					removeHalfEdgesArrPart(arr,edges_set,temp_segments,next_halfedge);
					//removeHalfEdgesArrPart(arr,edges_set,temp_segments,curr_halfedge);
					nextLoop(temp_segments, arr,edges_set, curr_halfedge,start_halfedge,isOrientable,loop_counter);
				}

				else 
				{
					if (close_loop_found)
					{
						#ifdef SHOW_LOOPS_CONST
						//global_graphics->clear();

						for (std::list<Halfedge_iterator>::iterator itr = temp_segments.begin();itr != temp_segments.end(); ++itr)
						{
							//global_graphics->draw_edge<Kernel>((*itr)->curve(),QColor(0,255,0),true);
						}
						//global_graphics->display();
						#endif
						removeHalfEdgesArrPartEnd(arr,edges_set,temp_segments,close_loop_handle,before_close_loop_handle);
						nextLoop(temp_segments, arr,edges_set, curr_halfedge,start_halfedge,isOrientable,loop_counter);
						close_loop_found = false;
					}
					else
					{// we are stuck.
						#ifdef SHOW_LOOPS_CONST
						//global_graphics->clear();

						for (std::list<Halfedge_iterator>::iterator itr = temp_segments.begin();itr != temp_segments.end(); ++itr)
						{
							//global_graphics->draw_edge<Kernel>((*itr)->curve(),QColor(0,255,0),true);
						}
						
						//global_graphics->display();
						#endif
						removeHalfEdgesArr(arr,edges_set,temp_segments);
						nextLoop(temp_segments, arr,edges_set, curr_halfedge,start_halfedge,isOrientable,loop_counter);
					}
				}
			}
			else{
			
				if (close_loop_found)
					{  // We have closed the loop, attempted to expand it and got stuck so remove everything but the loop.
						// THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????
						// THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????
						// THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????
						// THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????
						// THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????
						// THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????

						// THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????
						// THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????

						// THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????
						// THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????


						#ifdef SHOW_LOOPS_CONST
						//global_graphics->clear();

						for (std::list<Halfedge_iterator>::iterator itr = temp_segments.begin();itr != temp_segments.end(); ++itr)
						{
							//global_graphics->draw_edge<Kernel>((*itr)->curve(),QColor(0,255,0),true);
						}
						//global_graphics->display();
						#endif
						removeHalfEdgesArrPartEnd(arr,edges_set,temp_segments,close_loop_handle,before_close_loop_handle);
						nextLoop(temp_segments, arr,edges_set, curr_halfedge,start_halfedge,isOrientable,loop_counter);
						close_loop_found = false;
					}
				else
				{
					curr_halfedge = next_halfedge;
					temp_segments.push_back(curr_halfedge);
					int size_temp_segments = temp_segments.size();
					// check if loop has closed.
				//	if (next_halfedge->source() == start_halfedge->source() && next_halfedge->target() == start_halfedge->target())
					if (next_halfedge== start_halfedge)
					{ 
						#ifdef SHOW_LOOPS_CONST
						//global_graphics->clear();

						for (std::list<Halfedge_iterator>::iterator itr = temp_segments.begin();itr != temp_segments.end(); ++itr)
						{
							//global_graphics->draw_edge<Kernel>((*itr)->curve(),QColor(0,255,0),true);
						}

						//global_graphics->display();
						#endif
						removeHalfedgesSet(arr,edges_set,temp_segments);
						nextLoop(temp_segments, arr,edges_set, curr_halfedge,start_halfedge,isOrientable,loop_counter);					
					}
				}
			}

			

		}
		//printHe(curr_halfedge);

	}

	void nextLoop(std::list<Halfedge_iterator>& temp_segments,Arrangement_history_2& arr,Edges_set& edges_set,Halfedge_iterator& curr_halfedge,Halfedge_iterator& start_halfedge,bool& isOrientable,int& loop_counter) const
	{
		temp_segments.clear();
		if (edges_set.size()>0){
			++loop_counter;
			curr_halfedge =  (*edges_set.begin());
			while (curr_halfedge->visited && edges_set.size()>0 )
			{
				typename Edges_set::iterator testItem = edges_set.find(curr_halfedge);
				if (testItem == edges_set.end()){
					testItem = edges_set.find(((curr_halfedge)->twin()));
					if (testItem != edges_set.end()){
						edges_set.erase(testItem);	
					}
				}
				else
					edges_set.erase(testItem);	
				curr_halfedge =  (*edges_set.begin());
			}
			if (!checkTripSameDirWithSegment(arr,curr_halfedge))
				curr_halfedge = curr_halfedge->twin();
			setEdgeVisited(*curr_halfedge,true,loop_counter);
			temp_segments.push_back(curr_halfedge);
			start_halfedge = curr_halfedge;
		}
		else{
			isOrientable = false;
		}
	}

	void setEdgeVisited(Halfedge& he,bool value,int id) const
	{
		he.visited = value;
		he.twin()->visited = value;
		he.loopNumber = id;
		he.twin()->loopNumber = id;
	}

	void setEdgeDegenerate(Halfedge& he,bool value) const
	{
		he.isDegenerate = value;
		he.twin()->isDegenerate = value;
	}

	bool getEdgeDegenerate(Halfedge_handle& he) const
	{
		return he.isDegenerate;
	}

	bool getEdgeVisited(Halfedge_handle& he) const
	{
		return he.visited;
	}



	// Removes the edges matching to halfedges from the set
	void removeHalfedgesSet(Arrangement_history_2& arr,Edges_set&  edges_set,std::list<Halfedge_iterator>& temp_segments) const
	{
		removeHalfEdgesInner(arr,edges_set,temp_segments,false,false,false,*(temp_segments.begin()),*(temp_segments.begin()));
	}
	
	// Removes the edges matching to halfedges from the set and arrangment
	void removeHalfEdgesArr(Arrangement_history_2& arr,Edges_set& edges_set,std::list<Halfedge_iterator>& temp_segments) const
	{
		removeHalfEdgesInner(arr,edges_set,temp_segments,true,false,false,*(temp_segments.begin()),*(temp_segments.begin()));
	}

	// Remove edges before loop and keep loop 
	void removeHalfEdgesArrPart(Arrangement_history_2& arr,Edges_set& edges_set,std::list<Halfedge_iterator>& temp_segments,Halfedge_handle& partition) const
	{
		removeHalfEdgesInner(arr,edges_set,temp_segments,false,true,false,partition,*(temp_segments.begin()));
	}

	// Remove edges before loop and keep loop and remove edges after loop.
	void removeHalfEdgesArrPartEnd(Arrangement_history_2& arr,Edges_set& edges_set,std::list<Halfedge_iterator>& temp_segments,Halfedge_handle& partition,Halfedge_handle& end_partition) const
	{
		removeHalfEdgesInner(arr,edges_set,temp_segments,false,true,true,partition,end_partition);
	}

	// remove half edges from set or set and arrangment.
	void removeHalfEdgesInner(Arrangement_history_2& arr,Edges_set& edges_set,std::list<Halfedge_iterator>& temp_segments,bool remove_arr,bool range,bool end_range,Halfedge_handle& partition_itr,Halfedge_handle& partition_end_itr ) const
	{
		bool part_reached = false;
		 
		typename Edges_set::iterator partItem = edges_set.find(partition_itr);
		if (partItem == edges_set.end()){
			partItem = edges_set.find(((*partition_itr).twin()));
		}

		typename Edges_set::iterator part_item_end = edges_set.find(partition_end_itr);
		if (part_item_end == edges_set.end()){
			part_item_end = edges_set.find(((*partition_end_itr).twin()));
		}

		Halfedge_handle h_e_start; 
		if (partItem != edges_set.end()){
			h_e_start = *partItem;
		}
		Halfedge_handle h_e_finish; 
		if (part_item_end != edges_set.end())
		{
			h_e_finish = *part_item_end;		 
		}

		printHe(h_e_start);

		bool mark_change_part = false;

		for (typename std::list<Halfedge_iterator>::iterator itr = temp_segments.begin();itr != temp_segments.end();++itr)
		{
			//int a =  edges_set.size();
			//int b = temp_segments.size();
			//Halfedge_handle h1 = *itr;
			//Halfedge_handle h2 = *partItem;
			printHe(*itr);
			
			if (range && !part_reached && ((*itr == h_e_start) || (*itr == h_e_start->twin())) )
				part_reached = true;

			// If we have to close the loop check if range end has arrived
			if (range && end_range &&part_reached && (((*itr) == h_e_finish) || (*itr == h_e_finish->twin()) ) )
				mark_change_part = true;

			
			
			typename Edges_set::iterator testItem = edges_set.find(*itr);
			if (testItem == edges_set.end()){
				testItem = edges_set.find(((*itr)->twin()));
				if (testItem != edges_set.end()){
					edges_set.erase(testItem);
					if (remove_arr || (range && !part_reached))
						arr.remove_edge(*itr);
				}
			}
			else{
				edges_set.erase(testItem);
					if (remove_arr || (range && !part_reached))
						arr.remove_edge(*itr);
			}

			if (mark_change_part)
				part_reached = false;
		}
	}

	
	Halfedge_handle getLargestExitingClockwiseEdge(Arrangement_history_2& arr,Halfedge_iterator& curr_halfedge,Vertex& v_target,bool& next_edge_found,int loop_number,bool& close_loop_found,Halfedge_handle& handle_close_loop) const
	{
		Point_2 p_source = v_target.point();
		double x1 = CGAL::to_double(p_source.x());
		double y1 = CGAL::to_double(p_source.y());
		Halfedge_around_vertex_circulator itr = v_target.incident_halfedges();
		Halfedge_around_vertex_circulator start = itr;
		Halfedge_around_vertex_circulator tmp = itr;
		++itr;
		//unsigned int number_of_incident_edges = itr.size();
		if (itr == tmp){
			next_edge_found = false;
			return itr;
		}
//		--itr;
		bool found = false;
		Direction_2 entering_dir = f_direction((f_vector(curr_halfedge->source()->point(),curr_halfedge->target()->point())));

		// find first edge which is valid in direction with respect to original edges.
		//
		int count =0;
		bool res;
		
		while ((!(res = checkTripSameDirWithSegment(arr,((itr->twin())))) || (curr_halfedge == itr) || (itr->visited)) && (itr!=start) )
		{
			++itr;
			++ count;
		}
			Halfedge_handle maybe_visited;
		if (count == 1)
		{
		
			close_loop_found = false;
			if	(!checkOutgoingNotVisited(arr,v_target,maybe_visited,loop_number))
			{
				close_loop_found = true;
				handle_close_loop = maybe_visited;
			}
		}

		if (!res || itr->visited ){
			next_edge_found = false;
			return itr->twin();
		}

		if (count != 1)
		{
			//Halfedge_handle maybe_visited;
			close_loop_found = false;
			if	(!checkOutgoingNotVisited(arr,v_target,maybe_visited,loop_number))
			{
				close_loop_found = true;
				handle_close_loop = maybe_visited;
			}
		}
		
		//Halfedge& min_he = itr->twin();

		// Now we have the first edge, check if we can improve 
		Halfedge_around_vertex_circulator next_edge_itr = itr;
		++next_edge_itr;
		//Halfedge_iterator best_he =  
		count =0;
		Direction_2 min_edge_dir = f_direction((f_vector(itr->twin()->source()->point(),itr->twin()->target()->point())));
		while (next_edge_itr != itr) 
		{
			if (checkTripSameDirWithSegment(arr,((next_edge_itr->twin())))){
				Direction_2 new_edge_dir = f_direction((f_vector(next_edge_itr->twin()->source()->point(),next_edge_itr->twin()->target()->point())));
				// If the new edge improves the old one, ie it is larger and satisfies the improvment rule, which is under consideration.
				// currently, being the anticlocwise most if we are in right halfplane of the entering direction, and then choose anti clockwise most from left plane.
				if ( isDirImproving(min_edge_dir,entering_dir,new_edge_dir) && 
					  (!next_edge_itr->visited)
					  ){
					//(f_ccw_in_between(new_edge_dir,-entering_dir,min_edge_dir) &&  f_ccw_in_between(new_edge_dir,min_edge_dir,entering_dir))){
						itr = next_edge_itr;
						min_edge_dir = f_direction((f_vector(itr->twin()->source()->point(),itr->twin()->target()->point())));
				}
			}
			++count;
			++next_edge_itr;
		}
		if (close_loop_found)
		{
			Direction_2 new_edge_dir = f_direction((f_vector(maybe_visited->twin()->source()->point(),maybe_visited->twin()->target()->point())));
			if (isDirImproving(min_edge_dir,entering_dir,new_edge_dir))
			{
				next_edge_found = true;
				return maybe_visited->twin();
			}
		}
		
		next_edge_found = true;
		return itr->twin();
		
		
	}

	//bool isDirImproving(Direction_2& min_edge_dir,Direction_2& entering_dir,Direction_2& new_edge_dir) const
	//{
	///*	return (( ((entering_dir == new_edge_dir ) || (f_ccw_in_between(entering_dir,new_edge_dir,min_edge_dir) || (min_edge_dir == entering_dir) )) && (!f_ccw_in_between(min_edge_dir,-entering_dir,new_edge_dir)))  ||
	//				  (f_ccw_in_between(min_edge_dir,new_edge_dir,-entering_dir) &&  !f_ccw_in_between(entering_dir,min_edge_dir,new_edge_dir) && (min_edge_dir!=entering_dir)));*/
	//	// 
	//	if ((entering_dir == min_edge_dir ) || f_ccw_in_between(min_edge_dir,-entering_dir,entering_dir)) // if minimal dir is to the right
	//	{
	//		return f_ccw_in_between(new_edge_dir,-entering_dir,min_edge_dir);
	//	}
	//	else // min dir is to the left.
	//	{
	//		if (f_ccw_in_between(new_edge_dir,-entering_dir,entering_dir) || (entering_dir == new_edge_dir ) ) // if new dir is to the right it is better.
	//			return true;
	//		return f_ccw_in_between(new_edge_dir,min_edge_dir,-entering_dir); // else check improvment on left side.
	//	}
	//}

	bool isDirImproving(Direction_2& min_edge_dir,Direction_2& entering_dir,Direction_2& new_edge_dir) const
	{
	/*	return (( ((entering_dir == new_edge_dir ) || (f_ccw_in_between(entering_dir,new_edge_dir,min_edge_dir) || (min_edge_dir == entering_dir) )) && (!f_ccw_in_between(min_edge_dir,-entering_dir,new_edge_dir)))  ||
					  (f_ccw_in_between(min_edge_dir,new_edge_dir,-entering_dir) &&  !f_ccw_in_between(entering_dir,min_edge_dir,new_edge_dir) && (min_edge_dir!=entering_dir)));*/
		// 
		Direction_2 opp_enter =  -entering_dir;
		if ((opp_enter == min_edge_dir ) ) // if minimal dir equals -entering dir 
		{
			return  true;
		}
		else // 
		{
			return f_ccw_in_between(new_edge_dir,opp_enter,min_edge_dir);
		/*	if (f_ccw_in_between(new_edge_dir,-entering_dir,entering_dir) || (entering_dir == new_edge_dir ) ) // if new dir is to the right it is better.
				return true;
			return f_ccw_in_between(new_edge_dir,min_edge_dir,-entering_dir); // else check improvment on left side.*/
		}
	}


	/* 
		Checks that the edge leads to a vertex which an outgoing visited edge has been ie returns false if we close a loop.
		h returns the edge which begins the loop
	*/
	bool checkOutgoingNotVisited(Arrangement_history_2& arr,Vertex& v_target,Halfedge_handle& h,int loop_number) const
	{
		Point_2 p_source = v_target.point();
		
		Halfedge_around_vertex_circulator itr = v_target.incident_halfedges();
		Halfedge_around_vertex_circulator start = itr;
		
		do{
			if (checkTripSameDirWithSegment(arr,((itr->twin()))) && itr->visited && itr->loopNumber == loop_number)
			{
				h = itr;
				return false;
			}
		}while (++itr != start );
		return true;
	}

	bool checkTripSameDirWithSegment_bak(Arrangement_history_2& arr,Halfedge_handle he) const
	{
		Originating_curve_iterator segment_itr;// = arr.originating_curves_begin ( *he);
		//bool found = false;
		//segment = NULL;
		Direction_2 start_he_dir = f_direction((f_vector(he->source()->point(),he->target()->point())));
		
		for (segment_itr = arr.originating_curves_begin (he);segment_itr != arr.originating_curves_end(he);++segment_itr)
		{
			Segment_2 segment = *segment_itr;
			/*std::cout << "checktrip: " << std::endl;
			printSegment(segment);
			printHe(he);*/
			Direction_2 start_seg_dir = f_direction(f_vector(segment.source(),segment.target()));
			
			
			if (start_seg_dir == start_he_dir){
				return true;
				
			}
		}

		return false;
	}

	bool checkTripSameDirWithSegment(Arrangement_history_2& arr,Halfedge_handle he) const
	{
		Originating_curve_iterator segment_itr;// = arr.originating_curves_begin ( *he);
		//bool found = false;
		//segment = NULL;
		//Vector_2 start_he_dir = (f_vector(he->source()->point(),he->target()->point()));
		
		for (segment_itr = arr.originating_curves_begin (he);segment_itr != arr.originating_curves_end(he);++segment_itr)
		{
			Segment_2 segment = *segment_itr;
			/*std::cout << "checktrip: " << std::endl;
			printSegment(segment);
			printHe(he);*/
			//Vector_2 start_seg_dir = f_vector(segment.source(),segment.target());
			
			
			//if (CGAL::angle(start_he_dir,start_seg_dir)==CGAL::ACUTE){
			//if (f_angle(start_he_dir,start_seg_dir)==CGAL::ACUTE){
			
			/*
			if (CGAL::collinear_are_ordered_along_line(segment.source(),he->source()->point(),he->target()->point())){
				return true;}
			*/
			
			
			CGAL::Comparison_result c1 = f_compare_endpoints_xy(segment);
			CGAL::Comparison_result c2 = (CGAL::Comparison_result)he->direction();//f_compare_xy(he->source()->point(),he->target()->point());
			bool same_dir = (c1==c2);
			if (same_dir)
				return true;
				
			
		}

		return false;
	}

	bool checkTripNotSameDirWithSegment(Arrangement_history_2& arr,Halfedge_handle he) const
	{
		Originating_curve_iterator segment_itr;// = arr.originating_curves_begin ( *he);
		//bool found = false;
		//segment = NULL;
		//Vector_2 start_he_dir = (f_vector(he->source()->point(),he->target()->point()));
		
		for (segment_itr = arr.originating_curves_begin (he);segment_itr != arr.originating_curves_end(he);++segment_itr)
		{
			Segment_2 segment = *segment_itr;
			/*std::cout << "checktrip: " << std::endl;
			printSegment(segment);
			printHe(he);*/
			//Vector_2 start_seg_dir = f_vector(segment.source(),segment.target());
			
			
			//if (CGAL::angle(start_he_dir,start_seg_dir)==CGAL::ACUTE){
			//if (f_angle(start_he_dir,start_seg_dir)==CGAL::ACUTE){
			
			/*
			if (CGAL::collinear_are_ordered_along_line(segment.source(),he->source()->point(),he->target()->point())){
				return true;}
			*/
			
			
			CGAL::Comparison_result c1 = segment.label()._orientation;//f_compare_endpoints_xy(segment);
			//CGAL::Comparison_result c1 = segment.data();//f_compare_endpoints_xy(segment);
			CGAL::Comparison_result c2 = (CGAL::Comparison_result)he->direction();//f_compare_xy(he->source()->point(),he->target()->point());
			bool same_dir = (c1!=c2);
			if (same_dir)
				return true;
				
			
		}

		return false;
	}



	bool checkReflex(const Point_2& prev,const Point_2& curr,const Point_2& next) const
	{
		CGAL::Orientation res_ori = f_orientation (prev, curr, next);
		//if (f_orientation (prev, curr, next) == COLLINEAR)
		
	  return ((res_ori == RIGHT_TURN) || (res_ori == COLLINEAR));
	}

	/*bool checkSwept(const Point_2& prev,const Point_2& curr,const Point_2& next,const Point_2& start,const Point_2& end) const
	{
		Direction_2 dir_start = f_direction(f_vector(prev,curr));
		Direction_2 dir_end = f_direction(f_vector(curr,next));
		Direction_2 dir_new = f_direction(f_vector(start,end));
		return (dir_new == dir_start) || f_ccw_in_between(dir_new,dir_start,dir_end) || (dir_end == dir_new);
	}*/

	bool checkSwept(const Point_2& prev,const Point_2& curr,const Point_2& next,const Point_2& start,const Point_2& end, bool& isStartConcide,bool& isEndConcide) const
	{
		Direction_2 dir_start = f_direction(f_vector(prev,curr));
		Direction_2 dir_end = f_direction(f_vector(curr,next));
		Direction_2 dir_new = f_direction(f_vector(start,end));
		isStartConcide = dir_new == dir_start;
		isEndConcide = dir_end == dir_new;

		return isStartConcide || f_ccw_in_between(dir_new,dir_start,dir_end) ||isEndConcide;
	}

	bool checkSwept(Direction_2& dir_start,Direction_2& dir_end,Direction_2& dir_new, bool& isStartConcide,bool& isEndConcide) const
	{
	/*	Direction_2 dir_start = f_direction(f_vector(prev,curr));
		Direction_2 dir_end = f_direction(f_vector(curr,next));
		Direction_2 dir_new = f_direction(f_vector(start,end));*/
		isStartConcide = dir_new == dir_start;
		isEndConcide = dir_end == dir_new;

		return isStartConcide || f_ccw_in_between(dir_new,dir_start,dir_end) ||isEndConcide;
	}

	void printHe(Halfedge_handle he) const
	{
		//Segment_2 temp = he->curve();
		if (!HE_WRITE)
			return;
		Point_2 p_source = he->source()->point();
		Point_2 p_end = he->target()->point();
		double x1 = CGAL::to_double(p_source.x());
		double y1 = CGAL::to_double(p_source.y());
		double x2 = CGAL::to_double(p_end.x());
		double y2 = CGAL::to_double(p_end.y());

		std::cout << x1 << "," << y1 << "," << x2 << "," << y2 << std::endl;
	}

	void printSegment(Segment_2& seg) const
	{
		//Segment_2 temp = he->curve();
		if (!HE_WRITE)
			return;
		Point_2 p_source = seg.source();
		Point_2 p_end = seg.target();
		double x1 = CGAL::to_double(p_source.x());
		double y1 = CGAL::to_double(p_source.y());
		double x2 = CGAL::to_double(p_end.x());
		double y2 = CGAL::to_double(p_end.y());

		std::cout << x1 << "," << y1 << "," << x2 << "," << y2 << std::endl;
	}

	void draw_arr(Arrangement_history_2& arr) const
	{
		//QColor c(0,255,0);
		for (Edge_iterator itr = arr.edges_begin();itr!=arr.edges_end();++itr){
			 printSegment(itr->curve());
			//global_graphics->draw_edge<Kernel>(itr->curve(),c,true);
			
		}

		// draw isolated vertices
	/*	Face_handle unbounded_face = arr.unbounded_face();
	
		Hole_iterator hi;
		Ccb_halfedge_circulator perimiterFace;
		for (hi = unbounded_face->holes_begin(); hi != unbounded_face->holes_end(); ++hi) {
			 perimiterFace = *hi;
		}
		hi->*/
		Vertex_iterator vit;
		//QColor c2(0,255,255);
		for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
		
		  if (vit->is_isolated())
		  {
			  //global_graphics->draw_cross<Kernel>(vit->point(),c2,0.1);
		  }
    
		}

		
	}

	


};

template <class Kernel, class Container>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_2_ (const Polygon_2<Kernel,Container>& pgn1,
                 const Polygon_2<Kernel,Container>& pgn2)
{
  Minkowski_sum_by_convolution_lien_2<Kernel, Container>  mink_sum;
  Polygon_2<Kernel,Container>                        sum_bound;
  std::list<Polygon_2<Kernel,Container> >            sum_holes;

  if (pgn1.size() > pgn2.size())
    mink_sum (pgn1, pgn2, sum_bound, std::back_inserter(sum_holes));
  else
    mink_sum (pgn2, pgn1, sum_bound, std::back_inserter(sum_holes));

  return (Polygon_with_holes_2<Kernel,Container> (sum_bound,
                                                  sum_holes.begin(),
                                                  sum_holes.end()));
}

};

#endif
