// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
#ifndef CGAL_ARR_LANDMARKS_ANN_H
#define CGAL_ARR_LANDMARKS_ANN_H

/*! \file
* Definition of the Arr_landmarks_ann<Arrangement> template.
*/
#include <CGAL/basic.h>
#include <CGAL/Arr_point_location/ANN.h>"
#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>

//#define CGAL_LM_ANN_DEBUG
#ifdef CGAL_LM_ANN_DEBUG
	#define PRINT_ANN_DEBUG(expr)   std::cout << expr << std::endl
#else
	#define PRINT_ANN_DEBUG(expr)
#endif

CGAL_BEGIN_NAMESPACE

/*! \class
* A class that answers nearest neighbor queries.
* It recieves a set of points, and builds a kd-tree for them.
* Given a query point, it finds the closest point to the query.
*/
template <class Traits_>
class Arr_landmarks_ann 
{
public:
	typedef Traits_										                    Traits_2;
  typedef typename Traits_2::Approximate_number_type	  ANT;
	typedef typename Traits_2::Point_2					          Point_2;
	typedef Arr_landmarks_ann<Traits_2>	                  Self;

	class  NN_Point_2 {
	public:
		NN_Point_2() 
		{ 
			vec[0]= vec[1]  = 0; 
		}

    NN_Point_2 (Point_2 &pnt)
		{ 
			m_point = pnt;
      Traits_2 traits; 
			vec[0]= traits.approximate_2_object()(pnt, 0);
      vec[1]= traits.approximate_2_object()(pnt, 1);
		}

		NN_Point_2 (Point_2 &pnt, Object &object) 
		{ 
			m_object = object;
			m_point = pnt;
      Traits_2 traits; 
			vec[0]= traits.approximate_2_object()(pnt, 0);
      vec[1]= traits.approximate_2_object()(pnt, 1);
		}

		Object &object() { return m_object; }

		Point_2 &point() { return m_point; }

		ANT& x()
		{
		  return vec[ 0 ];
		}
		
		ANT& y()
		{
		  return vec[ 1 ];
		}

		bool operator==(const NN_Point_2& p) const 
		{	
			return (x() == p.x()) && (y() == p.y())  ;
		}

		bool  operator!=(const NN_Point_2& p) const 
		{ 
			return ! (*this == p); 
		}

		ANT		    vec[2];
		Object		m_object;
		Point_2		m_point;
	}; //end of class


	//class Construct_coord_iterator {
	//public:
	//	const ANT* operator()(const NN_Point_2& p) const 
	//	{ return static_cast<const ANT*>(p.vec); }

	//	const ANT* operator()(const NN_Point_2& p, int)  const
	//	{ return static_cast<const ANT*>(p.vec+2); }
	//};

  typedef std::list<NN_Point_2>                           Point_list;	
	typedef typename Point_list::const_iterator             Input_iterator;

protected:
	int					    nPts;					  // actual number of data points
	ANNpointArray		dataPts;				// data points
	ANNkd_tree*			kdTree;					// search structure

	bool            b_allocs;       //if the arrays were allocated
  Object *        objects_array;  //array of objects correspond to the dataPts
  Point_2 *       point_2_array;  //array of Point_2 correspond to the dataPts

private:

  /*! Copy constructor - not supported. */
  Arr_landmarks_ann (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

public:
	  /*! Default constructor. */
	  Arr_landmarks_ann () : 
      nPts(0), 
      dataPts(NULL), 
      kdTree(NULL),
      b_allocs(false),
      objects_array(NULL), 
      point_2_array(NULL)
	  {
	  }

	  /*! Distructor - must be because we allocated a tree */
	  ~Arr_landmarks_ann() 
	  {
		  clean();
	  }

	  /*! init should allocate the tree and initialize it with all points.
	  it will be better if we could create a blank tree, 
	  and afterwards add points to it
	  */
	  void init(Input_iterator begin, Input_iterator beyond) 
	  {
		  if (b_allocs || kdTree) {
			  std::cerr << "ERROR: init - tree exists" << std::endl;
			  return;
		  }

		  if (begin != beyond) {
        call_ann_to_build_tree(begin, beyond);
		  }
	  }

protected:
	  void call_ann_to_build_tree(Input_iterator begin, Input_iterator beyond) 
	  {
      Input_iterator iter;      
      int dim = 2;
      int i = 0; 
      for (iter = begin; iter != beyond; iter++)
      {
        nPts++;
      }
	    dataPts = annAllocPts(nPts, dim);			// allocate data points
      objects_array = new Object[nPts];
      point_2_array = new Point_2[nPts];
		  b_allocs = true;

      for (i=0, iter = begin; (i< nPts) && (iter != beyond); i++, iter++)
      {
        NN_Point_2 np = *iter;
        dataPts[i][0] = np.x();
        dataPts[i][1] = np.y();
        objects_array[i] = np.object();
        point_2_array[i] = np.point();
        PRINT_ANN_DEBUG ("i = "<< i
          << " dataPts[i][0] = "<< dataPts[i][0] 
          << " dataPts[i][1] = "<< dataPts[i][1] 
          << " point_2_array[i] = "<< point_2_array[i] );
      }

      kdTree = new ANNkd_tree(					// build search structure
        dataPts,					// the data points
        nPts,						// number of points
        dim);						// dimension of space

    }

public:

	  /*! clean - deletes the tree in order to create a new one later
	  */
	  void clean() 
	  {
      PRINT_ANN_DEBUG(" in clean function ...");
		  if (b_allocs){
        if (kdTree)
          delete kdTree;
        if (objects_array)
          delete [] objects_array;
        if (point_2_array)
          delete [] point_2_array;
		  }
		  b_allocs = false;
		  kdTree = NULL;
	    annClose();									// done with ANN
	  }

	  /*! find the nearest point to the query (and its location object)
	  */
	  Point_2 find_nearest_neighbor(Point_2 query, Object &obj) const
	  {
		  //create NN_Point_2 from Point_2 
		  NN_Point_2 nn_query(query);

      //create ANNPoint from nn_query
      ANNcoord ann_coord[2];
      ann_coord[0] = nn_query.x();
      ann_coord[1] = nn_query.y();
      ANNpoint queryPt = ann_coord;

      PRINT_ANN_DEBUG("queryPt = (" << queryPt[0]<< ", " 
                      << queryPt[1]    << ")");

      int	k	= 1;			// number of nearest neighbors
      double eps = 0;	// error bound

	    //ANNidxArray	nnIdx;					// near neighbor indices
	    //ANNdistArray dists;					// near neighbor distances
      ANNidx ann_idx;
      //nnIdx = &ann_idx;
      ANNdist ann_dist;
      //dists = &ann_dist;

    	kdTree->annkSearch(		// search
				queryPt,						// query point
				k,								  // number of near neighbors
				&ann_idx,						// nearest neighbors (returned)
				&ann_dist,					// distance (returned)
				eps);							  // error bound

      ann_dist = sqrt (ann_dist);
    //  std::cout << "\tNN:\tIndex\tDistance\n";
		  //for (int i = 0; i < k; i++) {			// print summary
			 // dists[i] = sqrt(dists[i]);			// unsquare distance
    //    std::cout << "\t" << i << "\t" << nnIdx[i] << "\t" << dists[i] << "\n";
		  //}

		  //get the object 
      //int nn_index = nnIdx[0];
      obj = objects_array[ann_idx];	
      PRINT_ANN_DEBUG("index = "<< ann_idx
                << ", nearest neighbor = "<< point_2_array[ann_idx]
                << ", dist = "<< ann_dist);
		  return (point_2_array[ann_idx]);
	  }
	 
};

CGAL_END_NAMESPACE


#endif
