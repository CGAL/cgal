// Author     : Nader Salman

#ifndef _Gyroviz_segmentation2_
#define _Gyroviz_segmentation2_


#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
using namespace std;



//#define cimg_convert_path "C:\Users\nsalman\Documents\ImageMagick\ImageMagick-6.4.0\VisualMagick\bin\convert"
#include <CImg.h>
using namespace cimg_library;


inline
CImg <unsigned char> to_grayscale(const CImg <unsigned char>& image)
{
	CImg <unsigned char> grayscale_img(image.width, image.height,1,3,0);
	//3 channels (for instance RGB)
	cimg_forXY(image,i,j) 
	{ 
		// The effective luminance of a pixel is calculated with the following formula:
		// Y=0.3RED+0.59GREEN+0.11BLUE
		grayscale_img(i,j,0,0) =  (char)(image(i,j,0,0)*0.3 + image(i,j,0,2)*0.59 + image(i,j,0,1)*0.11);
		grayscale_img(i,j,0,1) =  grayscale_img(i,j,0,0);
		grayscale_img(i,j,0,2) =  grayscale_img(i,j,0,0); 
	}
	return grayscale_img;
} 


/* Spatial filters operate over a locall neighborhood
* about a pixel in an image.
* Two kinds of Spatial filters are exposed in here:
*
*  1° Low-pass (or Smoothing filters):
*     Such filters consist of a weighted summation of
*     all pixels in a neighborhood, divided by some scaling amount.
*
*  2° High-pass (or Sharpening edge crispening)
*     Such filters will bring out the details in an image.
*/

//Spatial filters
//---------------

//3x3 Gaussian kernel 
inline
CImg <unsigned char> gauss3(const CImg <unsigned char>& image) {

  CImg <unsigned char> filtered (image.dimx(),image.dimy(),1,3,0);

  cimg_for_insideXY(image,row,column,1){

    filtered(row,column,0,0) = (4*image(row,column,0,0) + image(row-1,column-1,0,0) + 2*(row-1,column,0,0) + image(row-1,column+1,0,0) + 
      2*image(row,column-1,0,0) + 2*image(row,column+1,0,0) + image(row+1,column-1,0,0) + 2*image(row+1,column,0,0) + 
      image( row+1,column+1,0,0))/16 ;  

    filtered(row,column,0,1) = filtered(row,column,0,0);
    filtered(row,column,0,2) = filtered(row,column,0,0);
  }
  return filtered;
}


//3x3 Average kernel 
inline
CImg <unsigned char> average3(const CImg <unsigned char>& image) {

  CImg <unsigned char> filtered (image.dimx(),image.dimy(),1,3,0);

  cimg_for_insideXY(image,row,column,1){

    filtered(row,column,0,0) = (image(row,column,0,0) + image(row-1,column-1,0,0) + (row-1,column,0,0) + image(row-1,column+1,0,0) + 
      image(row,column-1,0,0) + image(row,column+1,0,0) + image(row+1,column-1,0,0) + image(row+1,column,0,0) + 
      image( row+1,column+1,0,0))/9 ;

    filtered(row,column,0,1) = filtered(row,column,0,0);
    filtered(row,column,0,2) = filtered(row,column,0,0);
  }
  return filtered;
}


//High-pass filters
//-----------------

//3x3 High-pass filter 
inline
CImg <unsigned char> highP1(const CImg <unsigned char>& image) {

	CImg <unsigned char> filtered (image.dimx(),image.dimy(),1,3,0); // displayed pixels are normalized

	cimg_for_insideXY(image,row,column,1){

		filtered(row,column,0,0) = 10*image(row,column,0,0) - image(row-1,column,0,0) - image(row,column-1,0,0) - image(row,column+1,0,0) - image(row+1,column,0,0) ;
		filtered(row,column,0,1) = filtered(row,column,0,0);
		filtered(row,column,0,2) = filtered(row,column,0,0);
	}
	return filtered;//.normalize(0,256);
}



//Gradient edge operators
//-----------------------
//Note: For all the gradient operators a Standard Threshold will be used


//Frei-Chen Mask
inline
CImg <unsigned char> grad_freiChen(const CImg <unsigned char>& image,int threshold) {

	float grad_x, grad_y;  
	CImg <unsigned char> grad (image.dimx(),image.dimy(),1,3,0);

	float sqrt_2 = sqrt(2.0);

	cimg_for_insideXY(image,row,column,1){
		grad_x = grad_y = 0; 

		/*The row and column gradients are normalized to provide unit-gain
		positive and negativee weighted averages about a separated edge position.
		Normalized by : 1/(2+sqrt(2))
		*/

		grad_x = (1/(2+sqrt_2))*(sqrt_2*image(row-1,column,0,0) + image(row-1,column-1,0,0) + image(row-1,column+1,0,0) -
			sqrt_2*image(row+1,column,0,0) - image(row+1,column+1,0,0) - image(row+1,column-1,0,0)); 
		grad_y = (1/(2+sqrt_2))*(sqrt_2*image(row,column-1,0,0) + image(row-1,column-1,0,0) + image(row+1,column-1,0,0) -
			sqrt_2*image(row,column+1,0,0) - image(row-1,column+1,0,0) - image(row+1,column+1,0,0));

		//Square-root gradient
		grad(row,column,0,0)= (char)sqrt(grad_x*grad_x + grad_y*grad_y);


		//Threshold selection
		if(grad(row,column,0,0) < threshold)
		{ 
			// Non edge is black (0,0,0)       // Temporarily gray for Opengl visualisation problems
			grad(row,column,0,0) = 128;
			grad(row,column,0,1) = 128;
			grad(row,column,0,2) = 128;
		}
		else 
		{
			// Edge is white
			grad(row,column,0,0)= 255;
			grad(row,column,0,1)= grad(row,column,0,0);
			grad(row,column,0,2)= grad(row,column,0,0);
		}
		//grad(row,column,0,0) *= 128.0/threshold;
		//grad(row,column,0,1)= grad(row,column,0,0);
		//grad(row,column,0,2)= grad(row,column,0,0);

	}

	return grad;
}



//Pixel difference Mask
//inline
//CImg <unsigned char> grad_pdiff(const CImg <unsigned char>& image,int threshold) {
//
//  unsigned char grad_x, grad_y;  
//  CImg <unsigned char> grad (image.dimx(),image.dimy(),1,1,0);
//
//  cimg_for_insideXY(image,row,column,1){
//    grad_x = grad_y = 0;  
//
//    grad_x = image(row,column) - image(row+1,column); 
//    grad_y = image(row,column) - image(row,column+1);
//
//    //Square-root gradient
//    grad(row,column)= sqrt(grad_x*grad_x + grad_y*grad_y);
//
//    //Threshold selection
//    if(grad(row,column) < threshold)
//      grad(row,column) = 0;
//  }
//
//  return grad;
//}


//Separated pixel difference Mask
//inline
//CImg <unsigned char> grad_spdiff(const CImg <unsigned char>& image,int threshold) {
//
//  unsigned char grad_x, grad_y;  
//  CImg <unsigned char> grad (image.dimx(),image.dimy(),1,1,0);
//
//  cimg_for_insideXY(image,row,column,1){
//    grad_x = grad_y = 0;  
//
//    grad_x = image(row-1,column) - image(row+1,column); 
//    grad_y = image(row,column-1) - image(row,column+1);
//
//    //Square-root gradient
//    grad(row,column)= sqrt(grad_x*grad_x + grad_y*grad_y);
//
//    //Threshold selection
//    if(grad(row,column) < threshold)
//      grad(row,column) = 0;
//  }
//
//  return grad;
//}


//Roberts Mask
//inline
//CImg <unsigned char> grad_roberts(const CImg <unsigned char>& image,int threshold) {
//
//  unsigned char grad_x, grad_y;  
//  CImg <unsigned char> grad (image.dimx(),image.dimy(),1,1,0);
//
//  cimg_for_insideXY(image,row,column,1){
//    grad_x = grad_y = 0;  
//
//    grad_x = image(row,column) - image(row+1,column+1); 
//    grad_y = image(row,column) - image(row-1,column+1);
//
//    //Square-root gradient
//    grad(row,column)= sqrt(grad_x*grad_x + grad_y*grad_y);
//
//    //Threshold selection
//    if(grad(row,column) < threshold)
//      grad(row,column) = 0;
//  }
//
//  return grad;
//}



//Prewitt Mask
//inline
//CImg<unsigned char> grad_prewitt(CImg<unsigned char> image,int threshold) {
//
//  unsigned char grad_x, grad_y;  
//  CImg <unsigned char> grad (image.dimx(),image.dimy(),1,1,0);
//
//  cimg_for_insideXY(image,row,column,1){
//    grad_x = grad_y = 0;     
//
//    /*The row and column gradients are normalized to provide unit-gain
//    positive and negative weighted averages about a separated edge position.
//    Normalized by : 1/3
//    */
//
//    grad_x = (1/3.0)*(image(row-1,column+1) + image(row-1,column) + image(row-1,column-1) -
//      image(row+1,column+1) - image(row+1,column) - image(row+1,column-1));
//    grad_y = (1/3.0)*(image(row-1,column-1) + image(row,column-1) + image(row+1,column-1) -
//      image(row-1,column+1) - image(row,column+1) - image(row+1,column+1));
//
//    //Square-root gradient
//    grad(row,column)= sqrt(grad_x*grad_x + grad_y*grad_y);
//
//    //Threshold selection
//    if(grad(row,column) < threshold)
//      grad(row,column) = 0;
//  }
//
//  return grad;
//}



//Sobel Mask
//inline
//CImg<unsigned char> grad_sobel(CImg<unsigned char> image,int threshold) {
//
//  unsigned char grad_x, grad_y;  
//  CImg <unsigned char> grad (image.dimx(),image.dimy(),1,1,0);
//
//  cimg_for_insideXY(image,row,column,1){
//    grad_x = grad_y = 0; 
//
//    /*The row and column gradients are normalized to provide unit-gain
//    positive and negative weighted averages about a separated edge position.
//    Normalized by : 1/4
//    */
//
//    grad_x = (1/4.0)*(2*image(row-1,column) + image(row-1,column-1) + image(row-1,column+1) -
//      2*image(row+1,column) - image(row+1,column+1) - image(row+1,column-1)); 
//    grad_y = (1/4.0)*(2*image(row,column-1) + image(row-1,column-1) + image(row+1,column-1) -
//      2*image(row,column+1) - image(row-1,column+1) - image(row+1,column+1));
//
//    //Square-root gradient
//    grad(row,column)= sqrt(grad_x*grad_x + grad_y*grad_y);
//
//    //Threshold selection
//    if(grad(row,column) < threshold)
//      grad(row,column) = 0;
//  }
//
//  return grad;
//}





//Laplacian edge operator
//-----------------------


//Torre and Poggio Mask
//inline
//CImg <unsigned char> lap_torPoggio(CImg<unsigned char> image,int threshold)  {
//
//  unsigned char lap_x, lap_y;
//  CImg <unsigned char> lap(image.dimx(),image.dimy(),1,1,0);
//
//  cimg_for_insideXY(image,row,column,1){
//
//    lap_x = lap_y = 0;
//
//    lap_x = -image(row-1,column-1) + 2*image(row,column) - image(row+1,column);
//    lap_y = -image(row,column+1) + 2*image(row,column) - image(row,column-1);
//
//    //Highlight the edges
//    lap(row,column)= lap_x+lap_y;
//
//    //Threshold selection
//    if(lap(row,column) < threshold)
//      lap(row,column) = 0;
//  }
//  return lap;
//}

//Prewitt Mask
//inline
//CImg <unsigned char> lap_prewitt(CImg<unsigned char> image,int threshold)  {
//
//  unsigned char lap_x, lap_y;
//  CImg <unsigned char> lap(image.dimx(),image.dimy(),1,1,0);
//
//  cimg_for_insideXY(image,row,column,1){
//
//    lap_x = lap_y = 0; 
//
//    lap_x = 2*image(row,column+1) + 2*image(row,column) + 2*image(row,column-1) - image(row-1,column+1) - image(row-1,column) - image(row-1,column-1) - image(row+1,column+1) - image(row+1,column) - image(row+1,column-1);
//
//    lap_y = 2*image(row-1,column) + 2*image(row,column) + 2*image(row+1,column) - image(row-1,column+1) - image(row,column+1) - image(row+1,column+1) - image(row-1,column-1) - image(row,column-1) - image(row+1,column-1);
//
//    //Highlight the edges
//    lap(row,column)= lap_x+lap_y;
//
//    //Threshold selection
//    if(lap(row,column) < threshold)
//      lap(row,column) = 0;
//  }
//  return lap;
//}





#endif // _Gyroviz_segmentation2_


