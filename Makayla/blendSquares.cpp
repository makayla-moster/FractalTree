
//g++ -o blendSquares.exe blendSquares.cpp -O2 -lgdi32

#include <iostream>
#include "CImg.h"

using namespace cimg_library;
using namespace std;

int main() {
CImg<unsigned char> blendSq(100, 100, 1, 3, 0);


// Corners:
blendSq(0, 0, 0, 2) = 0.0; //top left blue
blendSq(0, 0, 0, 1) = 0.0; //top left green
blendSq(0, 0, 0, 0) = 0.0; //top left red

blendSq(0, blendSq.height()-1, 0, 1) = 150.0; //bottom left green
blendSq(0, blendSq.height()-1, 0, 2) = 150.0; //bottom left blue
blendSq(0, blendSq.height()-1, 0, 0) = 0.0; // bottom left red

blendSq(blendSq.width()-1, 0, 0, 0) = 255.0; //top right red
blendSq(blendSq.width()-1, 0, 0, 1) = 0.0; //top right green
blendSq(blendSq.width()-1, 0, 0, 2) = 150.0; //top right blue

blendSq(blendSq.width()-1, blendSq.height()-1, 0, 0) = 255.0; //bottom right red
blendSq(blendSq.width()-1, blendSq.height()-1, 0, 1) = 255.0; //bottom right green
blendSq(blendSq.width()-1, blendSq.height()-1, 0, 2) = 255.0; //bottom right blue


// GET SIDES OF IMAGE
float j = 1;
// left side
while(j != blendSq.height()-1){
	 float i = 0;
	 float blue = (((blendSq(0,0, 0, 2))*(blendSq.height() - j))/blendSq.height()) + (blendSq(0,blendSq.height()-1, 0, 2))*(j/blendSq.height());
	 float green = (((blendSq(0,0, 0, 1))*(blendSq.height() - j))/blendSq.height()) + (blendSq(0,blendSq.height()-1, 0, 1))*(j/blendSq.height());
	 float red = (((blendSq(0,0, 0, 0))*(blendSq.height() - j))/blendSq.height()) + (blendSq(0,blendSq.height()-1, 0, 0))*(j/blendSq.height());
	 //cout<<m<<endl;
	 //cout<<blendSq(0,0, 0, 1)<<endl;
	 //cout<<blendSq(0,blendSq.height()-1, 0, 1)<<endl;
	 blendSq(i, j, 0, 2) = blue;
	 blendSq(i, j, 0, 0) = red;
	 blendSq(i, j, 0, 1) = green;
	 j+=1; 
	}

float y = 1;
//right side 
while(y != blendSq.height()-1){
	 float i = blendSq.width()-1;
	 float red2 = (((blendSq(blendSq.width()-1, 0, 0, 0))*(blendSq.height() - y))/blendSq.height()) + (blendSq(blendSq.width()-1, blendSq.height()-1, 0, 0))*(y/blendSq.height());
	 float blue2 = (((blendSq(blendSq.width()-1, 0, 0, 2))*(blendSq.height() - y))/blendSq.height()) + (blendSq(blendSq.width()-1, blendSq.height()-1, 0, 2))*(y/blendSq.height());
	 float green2 = (((blendSq(blendSq.width()-1, 0, 0, 1))*(blendSq.height() - y))/blendSq.height()) + (blendSq(blendSq.width()-1, blendSq.height()-1, 0, 1))*(y/blendSq.height());
	 //cout<<red2<<endl;
	 //cout<<blue2<<endl;
	 //cout<<green2<<endl;
	 //cout<<blendSq(blendSq.width()-1, 0, 0, 1)<<endl;
	 //cout<<blendSq(blendSq.width()-1, blendSq.height()-1, 0, 1)<<endl;
	 //cout<<green2<<endl;
	 //cout<<blendSq(0,blendSq.height()-1, 0, 1)<<endl;
	 blendSq(i, y, 0, 2) = blue2;
	 blendSq(i, y, 0, 0) = red2;
	 blendSq(i, y, 0, 1) = green2;
	 y+=1; 
	}
	
float i = 1;
while(i != blendSq.width()-1){
	 float j = 0;
	 float blue = (((blendSq(0,0, 0, 2))*(blendSq.width() - i))/blendSq.width()) + (blendSq(blendSq.width()-1, 0, 0, 2))*(i/blendSq.width());
	 float green = (((blendSq(0,0, 0, 1))*(blendSq.width() - i))/blendSq.width()) + (blendSq(blendSq.width()-1, 0, 0, 1))*(i/blendSq.width());
	 float red = (((blendSq(0,0, 0, 0))*(blendSq.width() - i))/blendSq.width()) + (blendSq(blendSq.width()-1, 0, 0, 0))*(i/blendSq.width());
	 blendSq(i, j, 0, 0) = red;
	 blendSq(i, j, 0, 1) = green;
	 blendSq(i, j, 0, 2) = blue;
	 i+=1;
}

float m = 1;
while (m != blendSq.height()){
	float n = 1;
	 while(n != blendSq.width()-1){ 
		 float blue = (((blendSq(0,m, 0, 2))*(blendSq.width() - n))/blendSq.width()) + (blendSq(blendSq.width()-1, m, 0, 2))*(n/blendSq.width());
		 float green = (((blendSq(0,m, 0, 1))*(blendSq.width() - n))/blendSq.width()) + (blendSq(blendSq.width()-1, m, 0, 1))*(n/blendSq.width());
		 float red = (((blendSq(0,m, 0, 0))*(blendSq.width() - n))/blendSq.width()) + (blendSq(blendSq.width()-1, m, 0, 0))*(n/blendSq.width());
		 blendSq(n, m, 0, 0) = red;
		 blendSq(n, m, 0, 1) = green;
		 blendSq(n, m, 0, 2) = blue;
	 n+=1;
	 }
m+=1;
}
	
	
	
  //Display and save image
  CImgDisplay disp(blendSq); 
  while (!disp.is_closed())
      disp.wait(); 
  
  blendSq.save("interpolationSq.bmp");






}

