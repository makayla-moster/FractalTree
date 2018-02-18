
//g++ -o blendTriangle.exe blendTriangle.cpp -lgdi32

#include <iostream>
#include "CImg.h"

using namespace cimg_library;
using namespace std;

int main() {
CImg<unsigned char> blendTri(100, 100, 1, 3, 0);

blendTri(49, 0, 0, 2) = 150.0; //top blue
blendTri(49, 0, 0, 1) = 0.0; //top green
blendTri(49, 0, 0, 0) = 255.0; //top red

blendTri(0, 20, 0, 1) = 0.0; //bottom left green
blendTri(0, 20, 0, 2) = 255.0; //bottom left blue
blendTri(0, 20, 0, 0) = 0.0; // bottom left red

blendTri(blendTri.width()-1, 20, 0, 0) = 0.0; //bottom right red
blendTri(blendTri.width()-1, 20, 0, 1) = 255.0; //bottom right green
blendTri(blendTri.width()-1, 20, 0, 2) = 100.0; //bottom right blue

float x0 = 49;
float y0 = 0;
float x1 = blendTri.width()-1;
float y1 = 20;
float x = x0 + 1;
float y = y0 + 0.5;
float d = (y0 - y1)*x + (x1 - x0)*y + (x0*y1) - (x1*y0);
float distanceA = sqrt(pow((x1-x0), 2) + pow((y1-y0), 2));

while(x != x1){
	float h = sqrt(pow((x-x0), 2) + pow((y-y0),2));
	float red = (((blendTri(49, 0, 0, 0))*(distanceA - h)/distanceA) + ((blendTri(blendTri.width()-1, 20, 0, 0))*(h/distanceA)));
	float green = (((blendTri(49, 0, 0, 1))*(distanceA - h)/distanceA) + ((blendTri(blendTri.width()-1, 20, 0, 1))*(h/distanceA)));
	float blue = (((blendTri(49, 0, 0, 2))*(distanceA - h)/distanceA) + ((blendTri(blendTri.width()-1, 20, 0, 2))*(h/distanceA)));
	blendTri(x, y, 0, 0) = red;
	blendTri(x, y, 0, 1) = green;
	blendTri(x, y, 0, 2) = blue;
	if(d<0){
		y += 1;
		d += (y0 - y1) + (x1 - x0);
	}
	else{
		d += (y0 - y1);
	}
	x+=1;
}

float i0 = 0;
float j0 = 20;
float i1 = blendTri.width()-1;
float j1 = 20;
float i = i0+1;
float j = j0+0.5;
float g = (j0 - j1)*i + (i1 - i0)*j + (i0*j1) - (i1*j0);
float distanceB = sqrt(pow((i1-i0),2) + pow((j1-j0),2));
	
while(i != i1){
	float red = (((blendTri(0, 20, 0, 0))*(blendTri.width() - i)/blendTri.width()) + ((blendTri(blendTri.width()-1, 20, 0, 0))*(i/blendTri.width())));
	float green = (((blendTri(0, 20, 0, 1))*(blendTri.width() - i)/blendTri.width()) + ((blendTri(blendTri.width()-1, 20, 0, 1))*(i/blendTri.width())));
	float blue = (((blendTri(0, 20, 0, 2))*(blendTri.width() - i)/blendTri.width()) + ((blendTri(blendTri.width()-1, 20, 0, 2))*(i/blendTri.width())));
	blendTri(i, j, 0, 0) = red;
	blendTri(i, j, 0, 1) = green;
	blendTri(i, j, 0, 2) = blue;
	if(g<0){
		j += 1;
		g += (j0 - j1) + (i1 - i0);
	}
	else{
		g += (j0 - j1);
	}
	i+=1;
}

float m0 = 0;
float n0 = 20;
float m1 = 49;
float n1 = 0;
float m = m0+1;
float n = n0-0.5;
float z = (n0 + n1)*m + (m1 - m0)*n - (m0*n1) - (m1*n0);
float distanceC = sqrt(pow((m1-m0),2) + pow((n1-n0),2));

while(m != m1){
	float l = sqrt(pow((m-m0),2) + pow((n-n0),2));
	float red = (((blendTri(0, 20, 0, 0))*(distanceC - l)/distanceC) + ((blendTri(49, 0, 0, 0))*(l/distanceC)));
	float green = (((blendTri(0, 20, 0, 1))*(distanceC - l)/distanceC) + ((blendTri(49, 0, 0, 1))*(l/distanceC)));
	float blue = (((blendTri(0, 20, 0, 2))*(distanceC - l)/distanceC) + ((blendTri(49, 0, 0, 2))*(l/distanceC)));
	blendTri(m, n, 0, 2) = blue;
	blendTri(m, n, 0, 1) = green;
	blendTri(m, n, 0, 0) = red;
	//cout<<n<<endl;
	if (n>0){
		if(z>0){
		 n -= 1;
		 //cout<<n<<" PointB"<<endl;
		 z += (n0 - n1) + (m1 - m0);
		}
	//cout<<n<<endl;
		else{
		 z += (n0 - n1);
		 }
	}
	else if (n<0){
		n+=1;
		z+= (n0 - n1);
	}
	m+=1;
	z = (n0 + n1)*m + (m1 - m0)*n - (m0*n1) - (m1*n0);
}

//float totalArea = (20*99)/2;
//midpoint @ (49,10)

//cout << distanceA << endl;

//cout << distanceB << endl;

//cout << distanceC <<endl;
float s = (distanceA + distanceB + distanceC)/2;
//cout << s << endl;
float totalArea = sqrt(s*(s-distanceA)*(s-distanceB)*(s-distanceC));
//cout << totalArea << endl;


// Colors image red
for(int yVal = 0; yVal<blendTri.height()-1; yVal++){
	int found = -1;
	int left, right = 0;
	for(int xVal = 0; xVal<blendTri.width(); xVal++){
		if((blendTri(xVal, yVal, 0, 0)>0.0) && found==-1){
		 left = xVal;
		 found *= -1;
		}
		if((blendTri(xVal, yVal, 0, 0)>0.0) && found ==1){
			right = xVal;
		}
	}
	for(int xVal = left; xVal < right; xVal++){
			float distance_a = sqrt(pow((xVal-x0), 2)+pow((yVal-y0), 2));
			float distance_b = sqrt(pow((i0-xVal), 2)+pow((j0-yVal), 2));
			float distance_c = sqrt(pow((xVal-x1), 2)+pow((yVal-y1), 2));
			float sa = (distance_a+distanceA+distance_c)/2;
			float sb = (distance_a+distance_b+distanceC)/2;
			float sc = (distanceB+distance_b+distance_c)/2;
			float areaA = sqrt(sa*(sa-distance_a)*(sa-distance_c)*(sa-distanceA));
			float areaB = sqrt(sb*(sb-distance_a)*(sb-distance_b)*(sb-distanceC));
			float areaC = sqrt(sc*(sc-distance_b)*(sc-distance_c)*(sc-distanceB));
			float red = ((areaA*blendTri(0,20,0,0)/totalArea)+((areaB*blendTri(blendTri.width()-1,20,0,0)))/totalArea)+((areaC*blendTri(49,0,0,0))/totalArea);
			float green = ((areaA*blendTri(0,20,0,1)/totalArea)+((areaB*blendTri(blendTri.width()-1,20,0,1)))/totalArea)+((areaC*blendTri(49,0,0,1))/totalArea);
			float blue = ((areaA*blendTri(0,20,0,2)/totalArea)+((areaB*blendTri(blendTri.width()-1,20,0,2)))/totalArea)+((areaC*blendTri(49,0,0,2))/totalArea);
			//cout<<"red "<<red<<endl;
			//cout<<"green "<<green<<endl;
			//cout<<"blue "<<blue<<endl;
			blendTri(xVal, yVal, 0, 0) = red;
			blendTri(xVal, yVal, 0, 1) = green;
			blendTri(xVal, yVal, 0, 2) = blue;
	}
}

 //Display and save image
  CImgDisplay disp(blendTri); 
  while (!disp.is_closed())
      disp.wait(); 
  
  blendTri.save("interpolationTri.bmp");
  }
