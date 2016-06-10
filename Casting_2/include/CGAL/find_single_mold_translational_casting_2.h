// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efif@gmail.com>

#ifndef CGAL_FIND_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H
#define CGAL_FIND_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H
#include <CGAL/Polygon_2.h>
#include <CGAL/enum.h>

#include <iostream>
#include<list>

template <typename Kernel>
class CircleArrangment
{

	typedef typename Kernel::Direction_2 Direction_2;

	inline static bool isOpenDirectionContainedInSegment(Direction_2 d,bool isDIsLeftClockwise, std::pair<Direction_2,Direction_2> s)
	{
		if((isDIsLeftClockwise && d==s.second )||(!isDIsLeftClockwise && d==s.first ))
			return false;
		return !d.counterclockwise_in_between(s.first,s.second);

	}
	//only one segment can be closed so if B is closed and A isnt, we dont care
	inline static bool isAcontainedInB(bool isAStartClosed,bool isAEndClosed,std::pair<Direction_2,Direction_2> A, std::pair<Direction_2,Direction_2> B)
	{
		//A is closed, B is open and they share an vertex -> A not contained in B
		if( ( isAStartClosed &&A.first == B.first )|| (isAEndClosed && A.second == B.second))
			return false;
		if((A.first == B.second ||B.first == A.second)&& A.first != A.second)
			return false;

		return (!A.first.counterclockwise_in_between(B.first,B.second)&&
				!A.second.counterclockwise_in_between(B.first,B.second)&&
				!A.first.counterclockwise_in_between(B.first,A.second));
	}

	class CircleArrangmentEdge
	{
	public :
		bool startIsClosed;
		Direction_2 edgeStartEngle; //the end is the start of the next edge
		char count; //number of outer circles that cover the edge (0/1/2+)
		size_t edgeIndex;
		CircleArrangmentEdge(Direction_2 edgeStartEngle,size_t edgeIndex,bool startIsClosed,bool setCountToZero=false)
		{
			this->startIsClosed=startIsClosed;
			this->edgeStartEngle= edgeStartEngle;
			if(setCountToZero)
			{
				this->count = 0;
			}
			else
			{
				this->count = 1;
				this->edgeIndex = edgeIndex;
			}
		}
		void plusplus(size_t edgeIndex)
		{
			if(this->count ==0)
			{
				this->edgeIndex = edgeIndex;
				this->count = 1;
			}
			else if(this->count ==1)
			{
				this->count =2;
			}

		}
		bool isCovered(){return count ==2;}

	};
	typedef typename std::list<struct CircleArrangmentEdge >::iterator CircleEdgeIterator;

	std::list<CircleArrangmentEdge >edges;

	void insertIfLegal(const CircleEdgeIterator curIt,const CircleEdgeIterator nextIt,const struct CircleArrangmentEdge &edge)
	{
		if(( edge.startIsClosed && !nextIt->startIsClosed) || edge.edgeStartEngle != nextIt->edgeStartEngle)
			if(( curIt->startIsClosed && !edge.startIsClosed) || edge.edgeStartEngle != curIt->edgeStartEngle)
				edges.insert(nextIt,edge);
	}

public :
	CircleArrangment(std::pair<Direction_2,Direction_2> firstSegmentOuterCircle)
	{
		edges.push_back(CircleArrangmentEdge(firstSegmentOuterCircle.first,0,false));
		edges.push_back(CircleArrangmentEdge(firstSegmentOuterCircle.second,0,true,true));
	}

	void addSegmentOuterCircle(std::pair<Direction_2,Direction_2> segmentOuterCircle, size_t edge_index)
	{
		std::pair<Direction_2,Direction_2> edge;
		bool isStartClosedSegment=edges.begin()->startIsClosed;
		bool isEndClosedSegment=edges.begin()->startIsClosed;
		edge.first = edges.begin()->edgeStartEngle;
		edge.second = edges.begin()->edgeStartEngle;
		CircleEdgeIterator nextIt=edges.begin();
		CircleEdgeIterator it=edges.begin();
		bool done = false;
		while(!done)
		{
			it = nextIt;
			nextIt=it;
			nextIt++;
			if(nextIt==edges.end()) {
				done =true;
				nextIt= edges.begin();
			}

			isStartClosedSegment = !isEndClosedSegment;
			isEndClosedSegment = !nextIt->startIsClosed;
			edge.first = edge.second;
			edge.second =nextIt->edgeStartEngle;

			if(it->count == 2)
			{

				continue;
			}
			if( isAcontainedInB(isStartClosedSegment,isEndClosedSegment,edge,segmentOuterCircle))
			{
				it->plusplus(edge_index);
				continue;
			}
			bool isStartContained = isOpenDirectionContainedInSegment(segmentOuterCircle.first,true,edge);
			bool isEndContained = isOpenDirectionContainedInSegment(segmentOuterCircle.second,false,edge);
			// o~~~~~~~~~~~~o  = new segment
			// ?------------?  = "old" segment (the edge from the array
			if(isStartContained)
			{
				if(isEndContained)
				{
					bool isordered =  !segmentOuterCircle.second.counterclockwise_in_between(segmentOuterCircle.first,edge.second);
					if(isordered)
					{
						// 		o~~~~~~~~~~~~o
						// ?-----------------------?
						// __________________________
						// ?----c
						// 		o~-~-~-~-~-~-o
						// 					 c-----?
						struct CircleArrangmentEdge edge2 = *it;
						edge2.startIsClosed=false;
						edge2.edgeStartEngle=segmentOuterCircle.first;
						edge2.plusplus(edge_index);
						this->insertIfLegal(it,nextIt,edge2);
						struct CircleArrangmentEdge edge3 = *it;
						edge3.startIsClosed=true;
						edge3.edgeStartEngle=segmentOuterCircle.second;
						this->insertIfLegal(it,nextIt,edge3);


					}
					else
					{
						// ...~~~~~~~~~o	o~~~~~~~~~~... (round)
						//		   ?------------?
						// __________________________
						// 		   ?~-~o
						//			   c----c
						// 		        	o-~-?

						struct CircleArrangmentEdge edge2 = *it;
						edge2.startIsClosed=true;
						edge2.edgeStartEngle=segmentOuterCircle.second;
						this->insertIfLegal(it,nextIt,edge2);
						struct CircleArrangmentEdge edge3 = *it;
						edge3.startIsClosed=false;
						edge3.edgeStartEngle=segmentOuterCircle.first;
						edge3.plusplus(edge_index);
						this->insertIfLegal(it,nextIt,edge3);
						it->plusplus(edge_index);
					}
				}
				else
				{
					// 		o~~~~~~~~~~~~o
					// ?-----------?
					//_____________________
					// ?----c
					// 		o-~-~-~?
					struct CircleArrangmentEdge edge2 = *it;
					edge2.startIsClosed=false;
					edge2.edgeStartEngle=segmentOuterCircle.first;
					edge2.plusplus(edge_index);
					this->insertIfLegal(it,nextIt,edge2);

				}
			}
			else
			{
				if(isEndContained)
				{
					// o~~~~~~~~~~~~o
					// 		?------------?
					//_____________________
					// 		?-~-~-~-o
					// 				c----?
					struct CircleArrangmentEdge edge2 = *it;
					edge2.startIsClosed=true;
					edge2.edgeStartEngle=segmentOuterCircle.second;
					it->plusplus(edge_index);
					this->insertIfLegal(it,nextIt,edge2);
				}
				//else -  no intersection, do noting


			}
		}

	}


	void mergeAdjacent2EdgesAndRemoveEmpty()
	{
		bool inTwoEdge=false;
		for(CircleEdgeIterator it=edges.begin();it!=edges.end();)
		{


			if(it->isCovered())
			{
				if(inTwoEdge)
				{
					it= edges.erase(it);
					continue;
				}
				inTwoEdge = true;
			}
			else
			{
				inTwoEdge = false;
			}
			it++;
		}
	}
	// debug function
	void printArrangement()
	{
		for(CircleEdgeIterator it=edges.begin();it!=edges.end();++it)
		{
			if(it->startIsClosed)
				std::cout<<")[";
			else
				std::cout<<"](";
			std::cout<<it->edgeStartEngle;
			std::cout<<","<<(int)it->count;

		}
		std::cout<<"\n\n";

	}
	template< typename OutputIterator>
	void getAll1Edges(OutputIterator oi)
	{
		for(CircleEdgeIterator it=edges.begin();it!=edges.end();)
		{

			if((*it).count == 1)
			{
				std::pair<size_t, 	std::pair<Direction_2,Direction_2> > edge;
				edge.first = (*it).edgeIndex;
				edge.second.first = (*it).edgeStartEngle;
				it++;
				edge.second.second = (*((it == edges.end())?edges.begin():(it))).edgeStartEngle;
				oi=(edge);
			}
			else
			{
				it++;
			}
		}
	}
	bool allIsCoveredTwice(){return edges.size()==1;}

};
template <typename Kernel>

inline std::pair<typename Kernel::Direction_2,typename Kernel::Direction_2> getSegmentOuterCircle(typename Kernel::Segment_2 seg,	CGAL::Orientation orientation)
{
	typename Kernel::Direction_2 forward( seg);
	typename Kernel::Direction_2 backward(-forward);
	return (orientation ==CGAL::Orientation::CLOCKWISE)?
			std::make_pair(backward, forward) : std::make_pair(forward, backward);

}
template <typename Kernel, typename OutputIterator>
OutputIterator
find_single_mold_translational_casting_2(const CGAL::Polygon_2<Kernel>& pgn, OutputIterator oi)
{
	typename CGAL::Polygon_2<Kernel>::Edge_const_iterator e_it = pgn.edges_begin();
	size_t edge_index =0;
	CGAL::Orientation poly_orientation =pgn.orientation();
	std::pair<typename Kernel::Direction_2,typename Kernel::Direction_2> segmentOuterCircle = getSegmentOuterCircle<Kernel>(*e_it++,poly_orientation);
	CircleArrangment<Kernel> circleArrangment(segmentOuterCircle);
	edge_index++;
	//circleArrangment.printArrangement();

	for(;e_it!= pgn.edges_end();e_it++,edge_index++)
	{
		segmentOuterCircle = getSegmentOuterCircle<Kernel>(*e_it,poly_orientation);
		circleArrangment.addSegmentOuterCircle(segmentOuterCircle,edge_index);
		circleArrangment.mergeAdjacent2EdgesAndRemoveEmpty();
		//	circleArrangment.printArrangement();
		if(circleArrangment.allIsCoveredTwice())
			return oi;
	}
	circleArrangment.getAll1Edges(oi);
	return oi;
}
#endif
