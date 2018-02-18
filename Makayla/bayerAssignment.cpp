//g++ -o Makayla.exe bayerAssignment.cpp -O2 -lgdi32

#include <iostream>
#include "CImg.h"

using namespace cimg_library;
using namespace std;

int main() {
  CImg<unsigned char> orig("cherryBlossoms.bmp");
  //Create empty image to fill in with colored Bayer Filter
  CImg<unsigned char> bayerFilter(orig.width(), orig.height(), 1, 3, 0);
  CImg<unsigned char> reconstructedImg(orig.width(), orig.height(), 1, 3, 0);
  
  //Copy over green intensities (every other pixel)
  for (int i = 0; i < orig.width(); i++){
	  for(int j = i%2; j < orig.height(); j+=2){
		  //if ( j == 0 )
			  //cout << "first row: " << (int)orig(i,j,0,1) << "\n";
		 //Copy over green channel intensity for this pixel 
		 bayerFilter(i,j,0,1) = orig(i,j,0,1); 
		 //red and blue are zero
		 bayerFilter(i,j,0,0) = 0.0;
		 bayerFilter(i,j,0,2) = 0.0;
	  }
  }
  
  //Copy over blue intensities (every other pixel of every other row)
  //Fill in here//
    for (int i = 0; i < orig.width(); i+=2){
	  for(int j = 1; j < orig.height(); j+=2){
		 bayerFilter(i,j,0,2) = orig(i,j,0,2); 
		// bayerFilter(i,j,0,1) = 0.0;
		 bayerFilter(i,j,0,0) = 0.0;
	  }
  }
  
  //Copy over red intensities (every other pixel of every other row)
  //Fill in here//  
    for (int i = 1; i < orig.width(); i+=2){
	  for(int j = 0; j < orig.height(); j+=2){
		 bayerFilter(i,j,0,0) = orig(i,j,0,0); 
		 bayerFilter(i,j,0,2) = 0.0;
		// bayerFilter(i,j,0,1) = 0.0;
	  }
  }
  
  //Display and save filtered image
  CImgDisplay disp(bayerFilter); 
  while (!disp.is_closed())
      disp.wait(); 
  
  bayerFilter.save("bayerFilterBlossoms.bmp");
  
  //Reconstruct image from filter
  //Fill in here//
  //GREEN INTENSITIES
 for (int i = 1; i < orig.width()-1; i++){
	  for(int j = i%2; j < orig.height()-1; j+=2){
		 if (j==0)
			 j+=2;
		reconstructedImg(i,j,0,1) = bayerFilter(i,j,0,1);
		if (i%2 == 0){
			reconstructedImg(i,j,0,0) = 0.5*bayerFilter(i+1, j, 0, 0) + 0.5*bayerFilter(i-1, j, 0, 0);
			reconstructedImg(i,j,0,2) = 0.5*bayerFilter(i, j+1, 0, 2) + 0.5*bayerFilter(i, j-1, 0, 2);
		}
		else if (i%2 == 1){
			reconstructedImg(i,j,0,0) = 0.5*bayerFilter(i, j+1, 0, 0) + 0.5*bayerFilter(i, j-1, 0, 0);	
			reconstructedImg(i,j,0,2) = 0.5*bayerFilter(i+1, j, 0, 2) + 0.5*bayerFilter(i-1, j, 0, 2);
		}
		 
		
	  }
  }
  
   //Copy over blue intensities (every other pixel of every other row)
  //BLUE INTENSITIES/
    for (int i = 0; i < orig.width(); i+=2){
	  for(int j = 1; j < orig.height(); j+=2){
		 //cout<< (int)bayerFilter(i,j,0,2) <<endl;
		 reconstructedImg(i,j,0,2) = bayerFilter(i,j,0,2);
		 reconstructedImg(i,j,0,1) = 0.25*bayerFilter(i+1,j,0,1) + 0.25*bayerFilter(i-1,j,0,1) + 0.25*bayerFilter(i,j+1,0,1) + 0.25*bayerFilter(i,j-1,0,1);
		 reconstructedImg(i,j,0,0) = 0.25*bayerFilter(i+1,j+1,0,0) + 0.25*bayerFilter(i-1,j-1,0,0) + 0.25*bayerFilter(i+1,j-1,0,0) + 0.25*bayerFilter(i-1,j+1,0,0);
	  }
  }

   //Copy over red intensities (every other pixel of every other row)
  //RED INTENSITIES//  
    for (int i = 1; i < orig.width(); i+=2){
	  for(int j = 0; j < orig.height(); j+=2){
		 //cout<< (int)bayerFilter(i,j,0,0) <<endl;
		 reconstructedImg(i,j,0,0) = bayerFilter(i,j,0,0); 
		 reconstructedImg(i,j,0,2) = 0.25*bayerFilter(i+1,j+1,0,2) + 0.25*bayerFilter(i-1,j-1,0,2) + 0.25*bayerFilter(i-1,j+1,0,2) + 0.25*bayerFilter(i+1,j-1,0,2);
		 reconstructedImg(i,j,0,1) = 0.25*bayerFilter(i+1,j,0,1) + 0.25*bayerFilter(i-1,j,0,1) + 0.25*bayerFilter(i,j+1,0,1) + 0.25*bayerFilter(i,j-1,0,1);
	  }
  }
  
  //////////////////BORDERS/////////////////
  
  //TOP GREEN//
  for (int i = 2; i < orig.width() - 1; i+=2){
	  int j = 0;
	  reconstructedImg(i,j,0,1) = bayerFilter(i,j,0,1);//green
	  reconstructedImg(i,j,0,0) = 0.5*bayerFilter(i-1,j,0,0) + 0.5*bayerFilter(i+1,j,0,0); //red
	  reconstructedImg(i,j,0,2) = bayerFilter(i,j+1,0,2); //blue
  }
  
  //TOP RED//
  //cout<< "PointA" << endl;
  for(int i = 1; i < orig.width()-1; i+=2){
		 int j = 0;
		 reconstructedImg(i,j,0,0) = bayerFilter(i,j,0,0); //red
		 reconstructedImg(i,j,0,1) = 0.333*((float)bayerFilter(i-1,j,0,1)) + 0.3333*((float)bayerFilter(i+1,j,0,1)) + 0.33333*((float)bayerFilter(i,j+1,0,1)); //green
		 //cout<< (int)reconstructedImg(i,j,0,1) << endl;
		 reconstructedImg(i,j,0,2) = 0.5*bayerFilter(i+1,j+1,0,2) + 0.5*bayerFilter(i-1,j+1,0,2); //blue
	  
  }
  
  
  //bottom//
  if((orig.height()-1)%2 == 1){
	  for(int i = 1; i < orig.width()-1; i+=2){
		  int j = orig.height()-1;
		  reconstructedImg(i,j,0,1) = bayerFilter(i,j,0,1);
		  reconstructedImg(i,j,0,0) = bayerFilter(i,j+1,0,0);
		  reconstructedImg(i,j,0,2) = 0.5*bayerFilter(i-1,j,0,2) + 0.5*bayerFilter(i+1,j,0,2);
	  }
	 for(int i = 2; i<orig.width()-1; i+=2){
		 int j = orig.height()-1;
		 reconstructedImg(i,j,0,2) = bayerFilter(i,j,0,2);
		 reconstructedImg(i,j,0,0) = 0.5*bayerFilter(i-1,j+1,0,0) + 0.5*bayerFilter(i+1,j+1,0,2);
		 reconstructedImg(i,j,0,1) = 0.333333*((float)bayerFilter(i-1,j,0,1)) + 0.33333*((float)bayerFilter(i,j+1,0,1)) + 0.33333*((float)bayerFilter(i+1,j,0,1));
	 }
  }
  
  else if((orig.height()-1)%2 == 0){
	  for(int i = 1; i<orig.width()-1; i+=2){
	  int j = orig.height()-1;
	  reconstructedImg(i,j,0,0) = bayerFilter(i,j,0,0);
	  reconstructedImg(i,j,0,1) = 0.33333*((float)bayerFilter(i-1,j,0,1)) + 0.333333*((float)bayerFilter(i,j+1,0,1)) + 0.3333333*((float)bayerFilter(i+1,j,0,1));
	  reconstructedImg(i,j,0,2) = 0.5*bayerFilter(i-1,j+1,0,2) + 0.5*bayerFilter(i+1,j+1,0,2);
  }
	  for(int i = 2; i<orig.width()-1; i+=2){
	  int j = orig.height()-1;
	  reconstructedImg(i,j,0,1) = bayerFilter(i,j,0,1);
	  reconstructedImg(i,j,0,0) = 0.5*bayerFilter(i-1,j,0,0) + 0.5*bayerFilter(i+1,j,0,0);
	  reconstructedImg(i,j,0,2) = bayerFilter(i+1,j,0,2);
	}
  }


  
  ///////////// EDGES OF THE PICTURE///////////////
  
  //LEFT EDGE BLUE//
  for(int j = 1; j<orig.height()-1; j+=2){
	  int i = 0;
	  reconstructedImg(i,j,0,2) = bayerFilter(i,j,0,2);
	  reconstructedImg(i,j,0,1) = 0.333333*((float)bayerFilter(i,j+1,0,1)) + 0.3333333*((float)bayerFilter(i,j-1,0,1)) + 0.3333333*((float)bayerFilter(i+1,j,0,1));
	  reconstructedImg(i,j,0,0) = 0.5*bayerFilter(i+1,j+1,0,0) + 0.5*bayerFilter(i+1,j-1,0,0);
  }
  
  //LEFT EDGE GREEN //
  for(int j=2; j<orig.height()-1 ; j+=2){
	  int i = 0;
	  reconstructedImg(i,j,0,1) = bayerFilter(i,j,0,1);
	  reconstructedImg(i,j,0,2) = 0.5*bayerFilter(i,j+1,0,2) + 0.5*bayerFilter(i,j-1,0,2);
	  reconstructedImg(i,j,0,0) = bayerFilter(i+1,j,0,0);
  }
  
  ////right edge////
  
  //blue//
  if (orig.width()%2 == 1){
	  for (int j =1; j<orig.height()-1;j+=2){
		  int i = orig.width()-1;
		  reconstructedImg(i,j,0,2) = bayerFilter(i,j,0,2);
		  reconstructedImg(i,j,0,1) = 0.33333*((float)bayerFilter(i,j+1,0,1)) + 0.33333*((float)bayerFilter(i,j-1,0,1)) + 0.33333*((float)bayerFilter(i-1,j,0,1));
		  reconstructedImg(i,j,0,0) = 0.5*bayerFilter(i-1,j+1,0,0) + 0.5*bayerFilter(i-1,j-1,0,0);
	  }
	  //green//
	  for (int j=2; j<orig.height()-1; j+=2){
		  int i = orig.width()-1;
		  reconstructedImg(i,j,0,1) = bayerFilter(i,j,0,1);
		  reconstructedImg(i,j,0,2) = 0.5*bayerFilter(i,j+1,0,2) + 0.5*bayerFilter(i,j-1,0,2);
		  reconstructedImg(i,j,0,0) = bayerFilter(i-1,j,0,0);
	  }
  }
  else if(orig.width()%2 == 0){
  //RIGHT EDGE RED//
	  for(int j =2; j<orig.height()-1; j+=2){
		  int i = orig.width()-1;
		  reconstructedImg(i,j,0,0) = bayerFilter(i,j,0,0);
		  reconstructedImg(i,j,0,1) = 0.33333*((float)bayerFilter(i,j+1,0,1)) + 0.33333*((float)bayerFilter(i-1,j,0,1)) + 0.3333333*((float)bayerFilter(i,j-1,0,1));
		  reconstructedImg(i,j,0,2) = 0.5*bayerFilter(i-1,j+1,0,2) + 0.5*bayerFilter(i-1,j-1,0,2);
	  }
  //RIGHT EDGE GREEN//
	  for(int j = 1; j<=orig.height()-1; j+=2){
		  int i = orig.width()-1;
		  reconstructedImg(i,j,0,1) = bayerFilter(i,j,0,1);
		  reconstructedImg(i,j,0,0) = 0.5*bayerFilter(i,j+1,0,0) + 0.5*bayerFilter(i,j-1,0,0);
		  reconstructedImg(i,j,0,2) = bayerFilter(i-1,j,0,2);
	  }
  }
  
  ////////////////CORNERS///////////////////////////
  
  //top left corner//
  reconstructedImg(0,0,0,1) = bayerFilter(0,0,0,1);
  reconstructedImg(0,0,0,0) = bayerFilter(1,0,0,0);
  reconstructedImg(0,0,0,2) = bayerFilter(0,1,0,2);
  
  
  //top right corner//
  if((orig.width()-1)%2 == 1){
	  int i = orig.width()-1;
	  int j = orig.height()-1;
	 reconstructedImg(i,0,0,0) = bayerFilter(i,0,0,0);
	 reconstructedImg(i,0,0,1) = 0.5*bayerFilter(i-1,0,0,1) + 0.5*bayerFilter(i,1,0,1);
	 reconstructedImg(i,0,0,2) = bayerFilter(i-1,j-1,0,2);
  }
  else if((orig.width()-1)%2 == 0){
	  int i = orig.width()-1;
	 reconstructedImg(i,0,0,1) = bayerFilter(i,0,0,1);//green
	 reconstructedImg(i,0,0,2) = bayerFilter(i-1,0,0,2);//red
	 reconstructedImg(i,0,0,0) = bayerFilter(i,1,0,0);//blue
  }
  
  //bottom left corner//
  if((orig.height()-1)%2 == 0){
	  int i = 0;
	  int j = orig.height()-1;
	  reconstructedImg(i,j,0,1) = bayerFilter(i,j,0,1); //green
	  reconstructedImg(i,j,0,2) = bayerFilter(i+1,j,0,2); //blue
	  reconstructedImg(i,j,0,0) = bayerFilter(i,j+1,0,0); //red
  }
  else if((orig.height()-1)%2 == 1){
	  int i = 0;
	  int j = orig.height()-1;
	  reconstructedImg(i,j,0,2) = bayerFilter(i,j,0,2); //blue
	  reconstructedImg(i,j,0,1) = 0.5*bayerFilter(i+1,j,0,1) + 0.5*bayerFilter(i,j+1,0,1); //green
	  reconstructedImg(i,j,0,0) = bayerFilter(i+1,j+1,0,0); //red
  }	  

  //bottom right corner//
  if((orig.height()-1)%2 == 0){
	  int i = orig.width()-1;
	  int j = orig.height()-1;
	  reconstructedImg(i,j,0,0) = bayerFilter(i,j,0,0);
	  reconstructedImg(i,j,0,1) = 0.5*bayerFilter(i-1,j,0,1) + 0.5*bayerFilter(i,j+1,0,1);
	  reconstructedImg(i,j,0,2) = bayerFilter(i-1,j+1,0,2);
  }
  else if((orig.height()-1)%2 == 1){
	  int i = orig.width()-1;
	  int j = orig.height()-1;
	  reconstructedImg(i,j,0,1) = bayerFilter(i,j,0,1);
	  reconstructedImg(i,j,0,0) = bayerFilter(i+1,j,0,0);
	  reconstructedImg(i,j,0,2) = bayerFilter(i,j+1,0,2);
  }
  
  
  //Display and save reconstructed image  
  //Fill in here//
	CImgDisplay disp2(reconstructedImg); 
  while (!disp2.is_closed())
      disp2.wait(); 
  
  reconstructedImg.save("reconstructedImg.bmp");
  
  return 0;
}
