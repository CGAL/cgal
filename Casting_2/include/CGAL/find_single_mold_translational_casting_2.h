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

/*
 * Terms:
 * 		point = Represented as  Direction_2. It is the intersection between the fitting Direction_2 and the unit circle
 *
 *		arc = Represented as A pair of point. clockwise arc between the first point and the second point.
 *		(each of its sides might be open or closed)
 *
 *		SegmentOuterCircle  = Arc that represent all the directions that points out from the polygon if it start from the
 *		fitting segment. This arc is always open half circle.
 */

/*! \CircleArrangment
    \brief This class represents an partition of the unit-circle in to cells of depth 0,1,2+
    \ where depth the number of inserted open half-circles inserted that covers this cell

	in addition this class contains some static functions that are in this class inorder of sharing its typedefs

	all the circle is always covered by some cell. there can't be an hole
 */
template <typename Kernel>
class CircleArrangment
{

	typedef typename Kernel::Direction_2 point;
	typedef std::pair<point,point>  arc;
	/*
	 * Reminder:
	 * 		point = Represented as  Direction_2. It is the intersection between the fitting Direction_2 and the unit circle
	 *
	 *		arc = Represented as A pair of point. clockwise arc between the first point and the second point.
	 *		(each of its sides might be open or closed)
	 */

	/*! \fn bool isOpenDirectionContainedInArc(point p,bool isCounterClockwise, arc a)
	    \brief checks whether an open epsilon area clockwise/counterclockwise from a point p is contained in an arc s
	    \param p -  a point .
	    \param isCounterClockwise - true: we care about the counterclockwise epsilon area of p. false: same with clockwise
	    \param a - an arc that should contain the epsilon area
	 */
	inline static bool isOpenDirectionContainedInArc(point p,bool isCounterClockwise, arc a)
	{
		if((isCounterClockwise && p==a.second )||(!isCounterClockwise && p==a.first ))
			return false;
		return !p.counterclockwise_in_between(a.first,a.second);

	}
	/*! \fn bool isAcontainedInB(bool isAStartClosed,bool isAEndClosed,arc A, arc B)
	    \brief checks whether an arc A is contained in an arc B
	    \param isAStartClosed - do A contains its start point (clockwise)
	    \param isAEndClosed - do A contains its end point (clockwise)
	    \param A - an arc
	    \param B - an *open* arc
	 */
	inline static bool isAcontainedInB(bool isAStartClosed,bool isAEndClosed,arc A, arc B)
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

	/*! \CircleArrangmentEdge
	    \brief This class represents a cells (a point or an arc) of depth 0,1,2+ in the CircleArrangment
	    \ where depth the number of inserted open half-circles inserted that cover this cell

		This edge (cell) is described by the first point of the edge (clockwise).
		The last point can be deduced by the next instance of CircleArrangmentEdge in the list in CircleArrangment

		this class also keeps the cell depth.
	 */
	class CircleArrangmentEdge
	{
	public :
		bool startIsClosed;
		point edgeStartEngle; //the end is the start of the next edge
		char count; //number of outer circles that cover the edge (0/1/2+)
		size_t edgeIndex; // the index of the polygon edge who's open half-circle covers this cell - only relevant if count ==1

		/*! \ctor CircleArrangmentEdge(point edgeStartEngle,size_t edgeIndex,bool startIsClosed,bool setCountToOne=true)
			    \brief create a new edge (Arc), this edge count must be 0 or 1
			    \param edgeStartEngle - the first point of the arc (clockwise)
			    \param edgeIndex - the index of the polygon edge who's open half-circle covers this cell - only relevant if count == 1
			    \param startIsClosed - is the point edgeStartEngle contained in this cell
			    \param setCountToOne - to set the count to one (or zero if this var is false)
		 */
		CircleArrangmentEdge(point edgeStartEngle,size_t edgeIndex,bool startIsClosed,bool setCountToOne=true)
		{
			this->startIsClosed=startIsClosed;
			this->edgeStartEngle= edgeStartEngle;
			this->count = (int) setCountToOne;
			this->edgeIndex = edgeIndex;
		}

		/*! \fn void plusplus(size_t edgeIndex)
				\brief add new polygon edge who's open half-circle covers this cell
			    \param edgeIndex - the index of this edge
			    increase the edge count by one (if it is 2+, it will stay 2+)
			    set this new edge to be the one covers the cell if the count was zero before. (only relevant if now count == 1)
		 */
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


	/*! \fn void insertIfLegal(const CircleEdgeIterator curIt,const CircleEdgeIterator nextIt,const struct CircleArrangmentEdge &edge)
			\brief add new edge to the arrangement if it won't create some empty edges
		    \param curIt - iterator to the edge before where the new edge should be inserted
		    \param nextIt - iterator to the edge after where the new edge should be inserted
		    \param edge - the new edge that should be inserted

		    Notice that nextIt is redundant since it can be deduced from curIt.
		    But it was easier for me to just send it as well.

	 */
	void insertIfLegal(const CircleEdgeIterator curIt,const CircleEdgeIterator nextIt,const struct CircleArrangmentEdge &edge)
	{
		if(( edge.startIsClosed && !nextIt->startIsClosed) || edge.edgeStartEngle != nextIt->edgeStartEngle)
			if(( curIt->startIsClosed && !edge.startIsClosed) || edge.edgeStartEngle != curIt->edgeStartEngle)
				edges.insert(nextIt,edge);
	}
	/*! \fn void mergeAdjacent2EdgesAndRemoveEmpty()
			\brief merge all the arcs that are adjacent and of depth 2+

			it doesn't merge the first and last ones since it is easier this way.
	 */
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

public :
	/*! \ctor CircleArrangment(arc firstSegmentOuterCircle)
			\brief creates an arrangement on circle with two edges the one covered by firstSegmentOuterCircle and the other one
		    \param firstSegmentOuterCircle - the outer circle of the first segment of the polygon.

		    Notice that you might consider implementing the ctor as an full circle of depth 0,
		    but it was much easier for me to ignore the case where the all circle is a single arc,
		    so I choose this implementation.

	 */
	CircleArrangment(arc firstSegmentOuterCircle)
	{
		edges.push_back(CircleArrangmentEdge(firstSegmentOuterCircle.first,0,false));
		edges.push_back(CircleArrangmentEdge(firstSegmentOuterCircle.second,0,true,false));
	}
	/*! \fn addSegmentOuterCircle(arc segmentOuterCircle, size_t edge_index)
			\brief updates the arrangement in respect to a new segment outer open circle
		    \param segmentOuterCircle - the outer circle of the current segment of the polygon.
		    \param edge_index - this segment id

			This is the main funtion of this code.
			It separates the cells in which the endpoints of the new arc is contained to two parts and increase count
			for all the cells that the new arc covers.

			In the end the function mergeAdjacent2EdgesAndRemoveEmpty is called to remove redundant cells

	 */
	void addSegmentOuterCircle(arc segmentOuterCircle, size_t edge_index)
	{
		arc edge;
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
			bool isStartContained = isOpenDirectionContainedInArc(segmentOuterCircle.first,true,edge);
			bool isEndContained = isOpenDirectionContainedInArc(segmentOuterCircle.second,false,edge);
			// o~~~~~~~~~~~~o  = new arc
			// ?------------?  = "old" arc (the edge from the array)
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
		mergeAdjacent2EdgesAndRemoveEmpty();
	}

#if 0
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
#endif


	/*! \fn void getAll1Edges(OutputIterator oi)
			\brief insert to oi all the cells in depth 1 i.e. the cells that represent legal pullout directions
			\param oi - the output iterator to put the cells in

			Puts in oi var of type pair<size_t, std::pair<Kernel::Direction_2 ,Kernel::Direction_2 > >
			foreach valid top edge.
			Should only be called after all of the polygon edges where inserted.
     */

	template< typename OutputIterator>
	void getAll1Edges(OutputIterator oi)
	{
		for(CircleEdgeIterator it=edges.begin();it!=edges.end();)
		{

			if((*it).count == 1)
			{
				std::pair<size_t, 	arc > edge;
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
	/*! \fn bool allIsCoveredTwice()
			\brief returns if all of the arrangement is in depth 2+

			Before running this run mergeAdjacent2EdgesAndRemoveEmpty() or addSegmentOuterCircle() which calls
			mergeAdjacent2EdgesAndRemoveEmpty().

			The funtions checks that the whole circle is a single cell, which can happen only if this cell is of depth 2,
			so there is no need to check the depth as well.

     */
	bool allIsCoveredTwice(){return edges.size()==1;}

};

/*! \fn std::pair<typename Kernel::Direction_2,typename Kernel::Direction_2> getSegmentOuterCircle(typename Kernel::Segment_2 seg,	CGAL::Orientation orientation)
 * 	\fn arc getSegmentOuterCircle(typename Kernel::Segment_2 seg,	CGAL::Orientation orientation)
		\brief returns the open outer half-circle of the edge.
		\param seg - the polygon segment
		\param orientation - the orientation of the segment (and the polygon).
							if CLOCKWISE then the outer half circle is to the left.
 */

template <typename Kernel>
inline std::pair<typename Kernel::Direction_2,typename Kernel::Direction_2> getSegmentOuterCircle(typename Kernel::Segment_2 seg,	CGAL::Orientation orientation)
{
	typename Kernel::Direction_2 forward( seg);
	typename Kernel::Direction_2 backward(-forward);
	return (orientation ==CGAL::Orientation::CLOCKWISE)?
			std::make_pair(backward, forward) : std::make_pair(forward, backward);

}



/*! \fn OutputIterator find_single_mold_translational_casting_2(const CGAL::Polygon_2<Kernel>& pgn, OutputIterator oi)
 	\brief returns all the possible top edges of the polygon and there pullout direction (with no rotation)
		\param pgn - the input polygon that we want to check if is castable or not.
		\param oi - the output iterator to put the top edges in
 */
template <typename Kernel, typename OutputIterator>
OutputIterator
find_single_mold_translational_casting_2(const CGAL::Polygon_2<Kernel>& pgn, OutputIterator oi)
{
	typedef typename CGAL::Polygon_2<Kernel>::Edge_const_iterator  PolygonEdgeIterator;
	typedef typename Kernel::Direction_2 point;
	typedef std::pair<point,point>  arc;
	/*
	 * Reminder:
	 * 		point = Represented as  Direction_2. It is the intersection between the fitting Direction_2 and the unit circle
	 *
	 *		arc = Represented as A pair of point. clockwise arc between the first point and the second point.
	 *		(each of its sides might be open or closed)
	 */

	PolygonEdgeIterator e_it = pgn.edges_begin();
	size_t edge_index =0;
	CGAL::Orientation poly_orientation =pgn.orientation();
	arc segmentOuterCircle = getSegmentOuterCircle<Kernel>(*e_it++,poly_orientation);
	CircleArrangment<Kernel> circleArrangment(segmentOuterCircle);
	edge_index++;
	for(;e_it!= pgn.edges_end();e_it++,edge_index++)
	{
		segmentOuterCircle = getSegmentOuterCircle<Kernel>(*e_it,poly_orientation);
		circleArrangment.addSegmentOuterCircle(segmentOuterCircle,edge_index);
		if(circleArrangment.allIsCoveredTwice())
			return oi;
	}
	circleArrangment.getAll1Edges(oi);
	return oi;
}
#endif
