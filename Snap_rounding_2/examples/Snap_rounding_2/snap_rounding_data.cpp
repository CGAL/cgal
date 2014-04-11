// Copyright 2009,2014 Max-Planck-Institute Saarbruecken (Germany).
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
// author(s)     : Waqar Khan <wkhan@mpi-inf.mpg.de>


/* Usage
 *
 *  This example converts arbitrary-precision arrangment into fixed-precision using Snap Rounding and by using INPUT DATA FROM A USER SPECIFIED FILE.
 *  <Argument 1> (Mandatory) path to the input file containing the arrangment information.
 *  <Argument 2> (Optional)  path to the output file where the results of snap rounding will be stored.
 *							Not providing this argument will print the result on standard output.
 *
 *  Input file format
 *  Line # 1: 		Number of line-segments present in the file.
 *  Line # 2 to N+1:	segment_start_point_x <space> segment_start_point_y <space> segment_end_point_x <space> segment_end_point_y
 *
 *  Each line should contain information about just one segment.
*/

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>
#include <fstream>

typedef CGAL::Quotient<CGAL::MP_Float>           Number_type;
typedef CGAL::Cartesian<Number_type>             Kernel;
typedef CGAL::Snap_rounding_traits_2<Kernel>     Traits;
typedef Kernel::Segment_2                        Segment_2;
typedef Kernel::Point_2                          Point_2;
typedef std::list<Segment_2>                     Segment_list_2;
typedef std::list<Point_2>                       Polyline_2;
typedef std::list<Polyline_2>                    Polyline_list_2;

int main(int argc, char* argv[])
{
	//if(argc > 3 || argc < 2)
	if(argc > 3)
	{
		std::cout<< "Incorrect input. <Arg 1> path to the INPUT file. <Arg 2> (optional) path to the OUTPUT file. No arguments to choose the default data file" << std::endl;
		return -1;
	}

	Segment_list_2 seg_list;
	Polyline_list_2 output_list;

	std::ifstream my_read_file;
	std::ofstream my_write_file;

	if(argc > 1)
		my_read_file.open(argv[1]);
	else
		my_read_file.open("data/snap_rounding_data");

	if(!my_read_file.is_open())
	{
		std::cout<< "Error opening the input file"<< std::endl;
		return -1;
	}

	if(argc==3)
	{
		my_write_file.open(argv[2]);

		if(!my_read_file.is_open())
		{
			std::cout<< "Error opening the output file"<< std::endl;
			return -1;
		}
	}

	unsigned int number_of_lines = 0;
	my_read_file >> number_of_lines;

	for(unsigned int i=0; i<number_of_lines; i++)
	{
		double point_start_x, point_start_y, point_end_x, point_end_y;
		my_read_file >> point_start_x;
		my_read_file >> point_start_y;
		my_read_file >> point_end_x;
		my_read_file >> point_end_y;

		seg_list.push_back(Segment_2(Point_2(point_start_x, point_start_y), Point_2(point_end_x, point_end_y)));
	}

	// Generate an iterated snap-rounding representation, where the centers of
	// the hot pixels bear their original coordinates, using 1 kd trees:
	CGAL::snap_rounding_2<Traits,Segment_list_2::const_iterator,Polyline_list_2>
	  									(seg_list.begin(), seg_list.end(), output_list, 1.0, true, false, 1);

 	int counter = 0;
	Polyline_list_2::const_iterator iter1;

	if(argc == 3) //output to the file
	{
		for (iter1 = output_list.begin(); iter1 != output_list.end(); ++iter1)
		{
		    my_write_file << "Polyline number " << ++counter << ":\n";
		    Polyline_2::const_iterator iter2;

		    for (iter2 = iter1->begin(); iter2 != iter1->end(); ++iter2)
		      my_write_file << "    (" << iter2->x() << ":" << iter2->y() << ")\n";
		}

		my_write_file.close();
	}
	else //output to std output
	{
		for (iter1 = output_list.begin(); iter1 != output_list.end(); ++iter1)
		{
		    std::cout << "Polyline number " << ++counter << ":\n";
		    Polyline_2::const_iterator iter2;

		    for (iter2 = iter1->begin(); iter2 != iter1->end(); ++iter2)
		      std::cout << "    (" << iter2->x() << ":" << iter2->y() << ")\n";
		}
	}

	my_read_file.close();

	return(0);
}
