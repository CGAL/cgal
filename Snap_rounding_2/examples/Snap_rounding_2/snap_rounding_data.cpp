#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>
#include <fstream>
#include <CGAL/Timer.h>

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
	if(argc > 3 || argc < 2)
	{
		std::cout<< "Incorrect input. please provide the file path of the data file only. (optionally) enter the output file name" << std::endl;
		return -1;
	}

	Segment_list_2 seg_list;
	Polyline_list_2 output_list;

	std::ifstream my_read_file;
	std::ofstream my_write_file;

	my_read_file.open(argv[1]);

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

	

	CGAL::Timer segment_creation_time, snap_rounding_time;

	unsigned int number_of_lines = 0;
	my_read_file >> number_of_lines;
	
	segment_creation_time.start();
	
	for(int i=0; i<number_of_lines; i++)
	{
		double point_start_x, point_start_y, point_end_x, point_end_y;
		my_read_file >> point_start_x;
		my_read_file >> point_start_y;
		my_read_file >> point_end_x;
		my_read_file >> point_end_y;
		
		seg_list.push_back(Segment_2(Point_2(point_start_x, point_start_y), Point_2(point_end_x, point_end_y)));
	}
	
	segment_creation_time.stop();


	snap_rounding_time.start();
	// Generate an iterated snap-rounding representation, where the centers of
	// the hot pixels bear their original coordinates, using 1 kd trees:
	CGAL::snap_rounding_2<Traits,Segment_list_2::const_iterator,Polyline_list_2>
	  									(seg_list.begin(), seg_list.end(), output_list, 1.0, true, false, 1);

	snap_rounding_time.stop();  

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
		
		my_write_file << "\n\nSegment creation of " << number_of_lines << " took: " << segment_creation_time.time() << " sec" <<std::endl;
		my_write_file << "\n\nSnap rounding took: " << snap_rounding_time.time() << " sec" <<std::endl;
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
	
	std::cerr << "\n\nSegment creation of " << number_of_lines << " took: " << segment_creation_time.time() << " sec" <<std::endl;
	std::cerr << "\n\nSnap rounding took: " << snap_rounding_time.time() << " sec" <<std::endl;
	}
	
	my_read_file.close();

	return(0);
}