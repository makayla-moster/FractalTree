//g++ -w -o Makayla.exe Fractal.cpp libglew32.dll.a libglfw3dll.a -I include -lOpenGL32 -L ./ -lglew32 -lglfw3

#include "gl_utils.h"
#include <GL/glew.h>		/* include GLEW and new version of GL on Windows */
#include <GLFW/glfw3.h>     /* GLFW helper library */
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stack>
#include <vector>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <time.h>
#include <assert.h>
#include <string>
#include <string.h>
#include <stdarg.h>   
using namespace std;

void multiply(GLfloat matrix1[], GLfloat matrix2[], GLfloat result[]){   // this is to multiply a 4x4 and a 4x1 matrix together
	result[0] = (matrix1[0]*matrix2[0]) + (matrix1[4]*matrix2[1]) + (matrix1[8]*matrix2[2]) + (matrix1[12]*matrix2[3]);
	result[1] = (matrix1[1]*matrix2[0]) + (matrix1[5]*matrix2[1]) + (matrix1[9]*matrix2[2]) + (matrix1[13]*matrix2[3]);
	result[2] = (matrix1[2]*matrix2[0]) + (matrix1[6]*matrix2[1]) + (matrix1[10]*matrix2[2]) + (matrix1[14]*matrix2[3]);
	result[3] = (matrix1[3]*matrix2[0]) + (matrix1[7]*matrix2[1]) + (matrix1[11]*matrix2[2]) + (matrix1[15]*matrix2[3]);
}

string generatePattern(){												//Generates a pattern to create a tree.
    int numIts = 4; // Number of iterations
    string pattern = "F"; //"[X]";    // Using F for the pattern 
    
    for (int i = 0; i < numIts; i++){
        string newPattern = ""; 
        for (int idx = 0; idx < pattern.length(); idx++){
            //cout << "char: " << pattern.substr(idx,1) << endl;
            if (pattern.substr(idx,1).compare("F") == 0) 
			newPattern += "F[F][-F][+F]";   //"F[-F][F][+F]"
            else if (pattern.substr(idx,1).compare("X") == 0) 
                newPattern += "F-[[X]+X]+F[+FX]-X";
            else{
                newPattern += pattern.substr(idx,1);
            }        
        }
        pattern = newPattern;
    }
    
    return pattern;    
}

int countLabel(string modelName, char label[]){
	int numLab = 0;
	
	FILE *objFile;
	objFile = fopen(modelName.c_str(),"r");

	char buf[128];
	while (fscanf(objFile, "%s", &buf) != EOF){
	    if (strcmp(buf,label) == 0)
		numLab++;
	}
	
	cout << "Model has " << numLab << " " << label << "\n";
	fclose(objFile);
	return numLab;

}

//read in vertices
void loadVertices(string modelName, GLfloat verts[]){	
	cout << "Loading vertices\n";
	int numVert = 0;

	FILE *objFile;
	objFile = fopen(modelName.c_str(),"r");

	float maxX = -100000000.0;
	float minX = 100000000.0;
	float maxY = -100000000.0;
	float minY = 100000000.0;
	float maxZ = -100000000.0;
	float minZ = 100000000.0;

	char buf[128];
	char label[] = "v";
	float a, b, c;
	while (fscanf(objFile, "%s", &buf) != EOF){ 
	    if (strcmp(buf,label) == 0){
		//cout << "vertex" << endl;
		fscanf(objFile, "%f %f %f\n", &a, &b, &c);
		//cout << "a,b,c: " << a << ", " << b << ", " << c << endl;
		if (a > maxX)
		    maxX = a;
		if (a < minX)
		    minX = a;
		if (b > maxY)
		    maxY = b;
		if (b < minY)
		    minY = b;
		if (c > maxZ)
		    maxZ = c;
		if (c < minZ)
		    minZ = c;
		verts[3*numVert + 0] = 1.0*a;
		verts[3*numVert + 1] = 1.0*b;
		verts[3*numVert + 2] = 1.0*c;
		numVert++;
	    }
	}

	float scaleX = maxX-minX;
	float scaleY = maxY-minY;
	float scaleZ = maxZ-minZ;
	float transX = 0.5*(maxX+minX);
	float transY = 0.5*(maxY+minY);
	float transZ = 0.5*(maxZ+minZ);
	cout << "scales: " << scaleX << ", " << scaleY << ", " << scaleZ << endl;

	for (int i = 0; i < numVert; i++){
	    verts[3*i+0] = (verts[3*i+0] - transX)/scaleX;
	    verts[3*i+1] = (verts[3*i+1] - transY)/scaleY;
	    verts[3*i+2] = (verts[3*i+2] - transZ)/scaleZ;
	}

	fclose(objFile);
	cout << "Done loading vertices\n";
}


void computeFaceNormals(GLfloat faceNormals[], GLfloat verts[], GLint faces[], int numFaces){
	for (int i = 0; i < numFaces; i++){
		int idx1 = faces[i*3 + 0];
		int idx2 = faces[i*3 + 1];
		int idx3 = faces[i*3 + 2];
		float ux = verts[idx1*3 + 0] - verts[idx2*3 + 0];
		float uy = verts[idx1*3 + 1] - verts[idx2*3 + 1];
		float uz = verts[idx1*3 + 2] - verts[idx2*3 + 2];
		float vx = verts[idx1*3 + 0] - verts[idx3*3 + 0];
		float vy = verts[idx1*3 + 1] - verts[idx3*3 + 1];
		float vz = verts[idx1*3 + 2] - verts[idx3*3 + 2];
		float nx = uy*vz - uz*vy;
		float ny = uz*vx - ux*vz;
		float nz = ux*vy - uy*vx;
		float mag = sqrt(nx*nx + ny*ny + nz*nz);
		
	//	cout << "avg norm" << nx << ", " << ny << ", " << nz << endl;
		//cout << "mag: " << mag << endl;
		faceNormals[3*i + 0] = nx/mag;
		faceNormals[3*i + 1] = ny/mag;
		faceNormals[3*i + 2] = nz/mag;
	}
	
}



void computeVertNormals(GLfloat normals[], GLfloat verts[], int numVerts, GLint faces[], int numFaces, GLfloat faceNormals[]){
	for (int i = 0; i < numVerts; i++){
		float avgX = 0.0;
		float avgY = 0.0;
		float avgZ = 0.0;
		int numF_vert = 0;
		
		//find all the faces that contain this vertex
		for (int j = 0; j < numFaces; j++){
			int found = 0;
			for (int k = 0; k < 3; k++)
				if (faces[j*3 + k] == i)
					found = 1;
			if (found){
				avgX += faceNormals[j*3 + 0];
				avgY += faceNormals[j*3 + 1];
				avgZ += faceNormals[j*3 + 2];				
				numF_vert++;				
			}			
		}
		
		avgX /= numF_vert;
		avgY /= numF_vert;
		avgZ /= numF_vert;
		//cout << "avg norm" << avgX << ", " << avgY << ", " << avgZ << endl;
		normals[i*3 + 0] = avgX;
		normals[i*3 + 1] = avgY;
		normals[i*3 + 2] = avgZ;
	}
}

/* void loadFaces(string modelName, GLint faces[]){	
	//cout << "Loading faces\n";
	int numFaces = 0;

	FILE *objFile;
	objFile = fopen(modelName.c_str(),"r");

	char buf[128];
	char label[] = "f";
	GLfloat a, b, c;
	while (fscanf(objFile, "%s", &buf) != EOF){ 
	    if (strcmp(buf,label) == 0){
		fscanf(objFile, "%f %f %f\n", &a, &b, &c);
		faces[3*numFaces + 0] = a-1;
		faces[3*numFaces + 1] = b-1;
		faces[3*numFaces + 2] = c-1;
		numFaces++;
	    }
	}

	fclose(objFile);
	//cout << "Done loading faces\n";
} */

void loadFaces(string modelName, GLint faces[]){    			//To read in Maya OBJ files.
    cout << "Loading new faces\n";

    FILE *objFile;
    objFile = fopen(modelName.c_str(),"r");
    int numFaces = 0;
    char buf[128];
    char label[] = "f";
    float a, b, c, a1, b1, c1, a2, b2, c2;
    while (fscanf(objFile, "%s", &buf) != EOF){ 
    //for (int z0 = 0; z0 < 100; z0++){
        //cout << "buf: " << buf << endl;
        if (strcmp(buf,label) == 0){
        //cout << "vertex" << endl;
        fscanf(objFile, "%f/%f/%f %f/%f/%f %f/%f/%f\n", &a, &b, &c, &a1, &b1, &c1, &a2, &b2, &c2);
        //cout << "a,b,c: " << a << ", " << b << ", " << c << endl;
        
                faces[3*numFaces+0] = a-1;
                faces[3*numFaces+1] = a1-1;
                faces[3*numFaces+2] = a2-1;
                numFaces++;
        }
    }

    fclose(objFile);
    cout << "Done loading faces\n";
}


int main() {	
 
	GLFWwindow *window = NULL;
	const GLubyte *renderer;
	const GLubyte *version;
	GLuint vao;
	GLuint vbo;
	
	int rotation;
	int count = 1;
	int countBracket = 1;
	string pattern = generatePattern();
	//cout << "pattern: " << pattern << endl;
	//cout << endl;
	for (int idx = 0; idx < pattern.length(); idx++){
			if (pattern.substr(idx, 1).compare("F") == 0){
				count ++;
				//cout << "F" << endl;
			}
			else if (pattern.substr(idx, 1).compare("]") == 0){
				countBracket ++;
				//cout << "]" << endl;
			}
		}
	//cout << "Bracket " << countBracket << endl << endl;
	int totalCount = count + countBracket;							//Total amount of points, including the backtracking points that are added for the lines.
	GLfloat branchPoints[totalCount*3]; 							//List of points to make the branches. Includes extra points for lines.
	GLfloat leafPoints[countBracket*3];								//List of points to place the leaves - only on the ends of branches though.
	GLfloat* result = new float[4];									//This is a new 4x1 matrix to acquire the new heading.	
	GLfloat* result2 = new float[4];
	stack<float> PositionStack;
	stack<float> HeadingStack;
	float rz = (90 * 3.14159) / 180;	//For the first rotation, so the trunk is 90 degrees from the bottom of the screen. Converts degrees to radians.
	float rx = (0 * 3.14159) / 180;
	float ry = (90 * 3.14159) / 180;
	int pointsCount = 0;
	int leafCount = 0;
	GLfloat sx = 1;
	GLfloat sy = 1;
	GLfloat sz = 1;
	GLfloat currentPosition[] = {0.0f, -0.25f, 0.0f, 1.0f};
	GLfloat currentHeading[] = {0.0f, 0.5f, 0.0f, 0.0f};
	GLfloat rotateZ[] = 
		{cos(rz),sin(rz),0,0,
		-sin(rz),cos(rz),0,0,
		0,0,1,0,
		0,0,0,1};
	GLfloat scale[] =
		{sx,0,0,0,
		 0,sy,0,0,
		 0,0,sz,0,
		 0,0,0,1};
	GLfloat rotateX[] = 
		{1,0,0,0,
		 0,cos(rx),-sin(rx),0,
		 0,sin(rx),cos(rx),0,
		 0,0,0,1};
	GLfloat rotateY[] = 
		{cos(ry),0,sin(ry),0,
		0,1,0,0,
		-sin(ry),0,cos(ry),0,
		0,0,0,1};

	branchPoints[pointsCount] = currentPosition[0];					//These lines add the first set of points to the list of points.
	//cout << "Points 1: " << branchPoints[pointsCount] << endl;
	pointsCount++;
	branchPoints[pointsCount] = currentPosition[1];
	//cout << "Points 2: " << branchPoints[pointsCount] << endl;
	pointsCount++;
	branchPoints[pointsCount] = currentPosition[2];
	//cout << "Points 3: " << branchPoints[pointsCount] << endl;
	pointsCount++;
	//branchPoints[pointsCount] = currentPosition[3];
	//cout << "Points 4: " << branchPoints[pointsCount] << endl;
	//pointsCount++;
	//cout << endl;
	
	for (int idx = 0; idx < pattern.length(); idx++){
		
		//cout << "STRING: " << pattern.substr(idx, 1) << endl << endl;
		
		if (pattern.substr(idx,1).compare("[") == 0){
			//cout << "Printing currentPosition before push " << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			PositionStack.push(currentPosition[3]);					//Pushes the currentPosition onto the PositionStack.
			PositionStack.push(currentPosition[2]);
			PositionStack.push(currentPosition[1]);
			PositionStack.push(currentPosition[0]);
			//cout << "Printing currentPosition after push " << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			//cout << endl;
			
			//cout << "Printing currentHeading before push " << currentHeading[0] << " " << currentHeading[1] << " " << currentHeading[2] << " " << currentHeading[3] << endl;
			HeadingStack.push(currentHeading[3]);					//Pushes the currentHeading onto the HeadingStack.
			HeadingStack.push(currentHeading[2]);
			HeadingStack.push(currentHeading[1]);
			HeadingStack.push(currentHeading[0]);
			//cout << "Printing currentHeading after push " << currentHeading[0] << " " << currentHeading[1] << " " << currentHeading[2] << " " << currentHeading[3] << endl;
			//cout << endl;
		}
		
		else if (pattern.substr(idx, 1).compare("]") == 0){
			//cout << "Printing currentPosition before pop " << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			currentPosition[0] = PositionStack.top();				//Sets the current position back to the top of the stack.
			leafPoints[leafCount] = currentPosition[0];				//Adds the point to the array which will be where I put my leaf OBJ.
			leafCount++;
			PositionStack.pop();									//Pops the current position from the top of the stack.
			
			currentPosition[1] = PositionStack.top();
			leafPoints[leafCount] = currentPosition[1];
			leafCount++;
			PositionStack.pop();
			
			currentPosition[2] = PositionStack.top();
			leafPoints[leafCount] = currentPosition[2];
			leafCount++;
			PositionStack.pop();
			
			currentPosition[3] = PositionStack.top();
			//leafPoints[leafCount] = currentPosition[3];
			//leafCount++;
			PositionStack.pop();
			//cout << "Printing currentPosition after pop " << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			//cout << endl;
			
			branchPoints[pointsCount] = currentPosition[0];			//Adds the currentPosition to the list of points.
			//cout << "Points 1: " << branchPoints[pointsCount] << endl;
			pointsCount++;
			branchPoints[pointsCount] = currentPosition[1];
			//cout << "Points 2: " << branchPoints[pointsCount] << endl;
			pointsCount++;
			branchPoints[pointsCount] = currentPosition[2];
			//cout << "Points 3: " << branchPoints[pointsCount] << endl;
			pointsCount++;
			//branchPoints[pointsCount] = currentPosition[3];
			//cout << "Points 4: " << branchPoints[pointsCount] << endl;
			//pointsCount++;
			//cout << endl; 
			
			//cout << "Printing currentHeading before pop " << currentHeading[0] << " " << currentHeading[1] << " " << currentHeading[2] << " " << currentHeading[3] << endl;			
			currentHeading[0] = HeadingStack.top();					//Sets the currentHeading to the top of the HeadingStack.
			HeadingStack.pop();										//Pops the currentHeading from the top of the stack.
			currentHeading[1] = HeadingStack.top();
			HeadingStack.pop();
			currentHeading[2] = HeadingStack.top();
			HeadingStack.pop();
			currentHeading[3] = HeadingStack.top();
			HeadingStack.pop();
			//cout << "Printing currentHeading after pop " << currentHeading[0] << " " << currentHeading[1] << " " << currentHeading[2] << " " << currentHeading[3] << endl;
			//cout << endl;
		}
		
		else if (pattern.substr(idx, 1).compare("F") == 0){
			
			//cout << " F before adding" << endl;
			//cout << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			
			currentPosition[0] += currentHeading[0]*.2;				//Changes the height of the tree, I like .2.
			currentPosition[1] += currentHeading[1]*.2;
			currentPosition[2] += currentHeading[2]*.2;
			currentPosition[3] += currentHeading[3];
			
			//cout << " F after adding" << endl;
			//cout << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			//cout << endl;
			
			branchPoints[pointsCount] = currentPosition[0];			//Adds the currentPosition to the list of points.
			//cout << "Points 1: " << branchPoints[pointsCount] << endl;
			pointsCount++;
			branchPoints[pointsCount] = currentPosition[1];
			//cout << "Points 2: " << branchPoints[pointsCount] << endl;
			pointsCount++;
			branchPoints[pointsCount] = currentPosition[2];
			//cout << "Points 3: " << branchPoints[pointsCount] << endl;
			pointsCount++;
			//branchPoints[pointsCount] = currentPosition[3];
			//cout << "Points 4: " << branchPoints[pointsCount] << endl;
			//pointsCount++;
			//cout << endl;
		}
		
		else if (pattern.substr(idx, 1).compare("+") == 0){
			//cout <<  "+ rz before " << rz << endl;
			rotation = rand() % 65 + 1;							//Chooses a random number to rotate by from 0 to 65.
			float rz =-  ((rotation * 3.14159) / 180);			//Converts degrees of the rotation to radians.
			rotateZ[0] = cos(rz);
			rotateZ[1] = sin(rz);
			rotateZ[4] = -sin(rz);
			rotateZ[5] = cos(rz);
			//cout <<  "+ rz after " << rz << endl;
			//cout << endl;
			multiply(rotateZ, currentHeading, result);		//Multiplies the rotateZ matrix by the currentHeading to get the new currentHeading.
			float magnitude = sqrt((result[0]*result[0]) + (result[1]*result[1]) + (result[2]*result[2]) + (result[3]*result[3]));	//Finds magnitude for normalization of heading.
			currentHeading[0] = result[0] / magnitude;		//Normalizes the currentHeading vector.
			currentHeading[1] = result[1] / magnitude;
			currentHeading[2] = result[2] / magnitude;
			currentHeading[3] = result[3] / magnitude;
			//cout << "Current Heading after Normalization: " << currentHeading[0] << " " << currentHeading[1] << " " << currentHeading[2] << " " << currentHeading[3] << endl;
			//cout << endl;
		}
		
		else if (pattern.substr(idx, 1).compare("-") == 0){
			//cout <<  "- rz before " << rz << endl;
			rotation = rand() % 65 + 1;							//Chooses a random number to rotate by from 0 to 65.
			float rz =+ ((rotation * 3.14159) / 180);			//Converts degrees of the rotation to radians.
			rotateZ[0] = cos(rz);
			rotateZ[1] = sin(rz);
			rotateZ[4] = -sin(rz);
			rotateZ[5] = cos(rz);
			//cout <<  "- rz after " << rz << endl;
			//cout << endl;
			multiply(rotateZ, currentHeading, result);			//Multiplies the rotateZ matrix by the currentHeading to get the new currentHeading.
			float magnitude = sqrt((result[0]*result[0]) + (result[1]*result[1]) + (result[2]*result[2]) + (result[3]*result[3]));	//Finds magnitude for normalization of heading.
			currentHeading[0] = result[0] / magnitude;			//Normalizes the currentHeading vector.
			currentHeading[1] = result[1] / magnitude;
			currentHeading[2] = result[2] / magnitude;
			currentHeading[3] = result[3] / magnitude;
			/*cout << "Current Heading after Normalization: " << currentHeading[0] << " " << currentHeading[1] << " " << currentHeading[2] << " " << currentHeading[3] << endl;
			cout << endl;*/
		}
	}

	string modelName = "Leaf.obj";								//Name of the OBJ to load.
	
	int numVert = countLabel(modelName, "v");
	GLfloat* verts = new GLfloat[3*numVert];
	loadVertices(modelName, verts);
	
	int numFaces = countLabel(modelName, "f");
	GLint* faces = new GLint[3*numFaces];
	loadFaces(modelName, faces);
	
	GLfloat* faceNormals = new GLfloat[3*numFaces];
	computeFaceNormals(faceNormals, verts, faces, numFaces);
	
	GLfloat* vertNormals = new GLfloat[3*numVert];
	computeVertNormals(vertNormals, verts, numVert, faces, numFaces, faceNormals);
	
	GLfloat* points = new GLfloat[9*numFaces];
	GLfloat* normals = new GLfloat[9*numFaces];
	for (int i = 0; i < numFaces; i++){
	    int idx1 = faces[3*i + 0];
	    int idx2 = faces[3*i + 1];
	    int idx3 = faces[3*i + 2];
	    points[i*9 + 0] = verts[3*idx1+0];
	    points[i*9 + 1] = verts[3*idx1+1];
	    points[i*9 + 2] = verts[3*idx1+2];
	    points[i*9 + 3] = verts[3*idx2+0];
	    points[i*9 + 4] = verts[3*idx2+1];
	    points[i*9 + 5] = verts[3*idx2+2];
	    points[i*9 + 6] = verts[3*idx3+0];
	    points[i*9 + 7] = verts[3*idx3+1];
	    points[i*9 + 8] = verts[3*idx3+2];
	    normals[i*9 + 0] = vertNormals[3*idx1+0];
	    normals[i*9 + 1] = vertNormals[3*idx1+1];
	    normals[i*9 + 2] = vertNormals[3*idx1+2];
	    normals[i*9 + 3] = vertNormals[3*idx2+0];
	    normals[i*9 + 4] = vertNormals[3*idx2+1];
	    normals[i*9 + 5] = vertNormals[3*idx2+2];
	    normals[i*9 + 6] = vertNormals[3*idx3+0];
	    normals[i*9 + 7] = vertNormals[3*idx3+1];
	    normals[i*9 + 8] = vertNormals[3*idx3+2];
	}
	int numPoints = 3*numFaces;
	
	/* these are the strings of code for the shaders
	the vertex shader positions each vertex point */
	const char *vertex_shader = "#version 410\n"
		"attribute vec3 vp;"
		"void main () {"
		"  gl_Position = vec4 (vp, 1.0);"
		"}";
	/* the fragment shader colours each fragment (pixel-sized area of the
	triangle) */
	const char *fragment_shader = "#version 410\n"
		"out vec4 frag_colour;"
		"void main () {"
		"  frag_colour = vec4 (0.80, 0.5215, 0.247, 1.0);"
		"}";
	/* GL shader objects for vertex and fragment shader [components] */
	GLuint vert_shader, frag_shader;
	/* GL shader program object [combined, to link] */
	GLuint shader_programme;
	
	//-----------------------------------------------------------------------------
	/* these are the strings of code for the shaders
	the vertex shader positions each vertex point */
	const char *vertex_shader2 = "#version 410\n"
		"attribute vec3 vp;"
		"uniform mat4 rotateX;"
		"void main () {"
		"  gl_Position = rotateX * vec4 (vp, 1.0);"
		"}";
	/* the fragment shader colours each fragment (pixel-sized area of the
	triangle) */
	const char *fragment_shader2 = "#version 410\n"
		"out vec4 frag_colour;"
		"void main () {"
		"  frag_colour = vec4 (0.0, 1.0, 0.5, 1.0);"
		"}";
	/* GL shader objects for vertex and fragment shader [components] */
	GLuint vert_shader2, frag_shader2;
	/* GL shader program object [combined, to link] */
	GLuint shader_programme2;
	//----------------------------------------------------------------------------

	/* start GL context and O/S window using the GLFW helper library */
	if ( !glfwInit() ) {
		fprintf( stderr, "ERROR: could not start GLFW3\n" );
		return 1;
	}

	window = glfwCreateWindow(640, 480, "Fractal Window", NULL, NULL);
	if ( !window ) {
		fprintf( stderr, "ERROR: could not open window with GLFW3\n" );
		glfwTerminate();
		return 1;
	}
	
	glfwMakeContextCurrent( window );
	/* start GLEW extension handler */
	glewExperimental = GL_TRUE;
	glewInit();

	/* get version info */
	renderer = glGetString( GL_RENDERER ); /* get renderer string */
	version = glGetString( GL_VERSION );	 /* version as a string */
	printf( "Renderer: %s\n", renderer );
	printf( "OpenGL version supported %s\n", version );

	/* tell GL to only draw onto a pixel if the shape is closer to the viewer
	than anything already drawn at that pixel */
	glEnable( GL_DEPTH_TEST ); /* enable depth-testing */
	/* with LESS depth-testing interprets a smaller depth value as meaning "closer" */
	glDepthFunc( GL_LESS );
	/* a vertex buffer object (VBO) is created here. this stores an array of
	data on the graphics adapter's memory. in our case - the vertex points */
	glGenBuffers( 1, &vbo );
	glBindBuffer( GL_ARRAY_BUFFER, vbo );
	glBufferData( GL_ARRAY_BUFFER, (totalCount*3) * sizeof( GLfloat ), branchPoints, GL_STATIC_DRAW ); // count*3 to not lose points 

	/* the vertex array object (VAO) is a little descriptor that defines which
	data from vertex buffer objects should be used as input variables to vertex
	shaders. in our case - use our only VBO, and say 'every three floats is a
	variable' */
	glGenVertexArrays( 1, &vao );
	glBindVertexArray(vao);
	// "attribute #0 should be enabled when this vao is bound"
	glEnableVertexAttribArray(0);
	// this VBO is already bound, but it's a good habit to explicitly specify which
	// VBO's data the following
	// vertex attribute pointer refers to
	glBindBuffer( GL_ARRAY_BUFFER, vbo );
	// "attribute #0 is created from every 3 variables in the above buffer, of type
	// float (i.e. make me vec3s)"
	glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, NULL );
	
	//----------------------------------------------------------------------------------------------------
	GLuint points_vbo;
	glGenBuffers (1, &points_vbo);
	glBindBuffer (GL_ARRAY_BUFFER, points_vbo);
	glBufferData (GL_ARRAY_BUFFER, 3 * numPoints * sizeof (GLfloat), points, GL_STATIC_DRAW);
	
	GLuint normals_vbo;
	glGenBuffers (1, &normals_vbo);
	glBindBuffer (GL_ARRAY_BUFFER, normals_vbo);
	glBufferData (GL_ARRAY_BUFFER, 3 * numPoints * sizeof (GLfloat), normals, GL_STATIC_DRAW);
	
	GLuint vao2;
	glGenVertexArrays (1, &vao2);
	glBindVertexArray (vao2);
	glBindBuffer (GL_ARRAY_BUFFER, points_vbo);
	glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glBindBuffer (GL_ARRAY_BUFFER, normals_vbo);
	glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray (0);
	glEnableVertexAttribArray (1);
	//----------------------------------------------------------------------------------------------------

	/* here we copy the shader strings into GL shaders, and compile them. we
	then create an executable shader 'program' and attach both of the compiled
	shaders. we link this, which matches the outputs of the vertex shader to
	the inputs of the fragment shader, etc. and it is then ready to use */
	vert_shader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vert_shader, 1, &vertex_shader, NULL);
	glCompileShader(vert_shader);
	frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(frag_shader, 1, &fragment_shader, NULL);
	glCompileShader(frag_shader);
	shader_programme = glCreateProgram();
	glAttachShader(shader_programme, frag_shader);
	glAttachShader(shader_programme, vert_shader);
	glLinkProgram(shader_programme);
	glPointSize(5.0);
	
	//------------------------------------------------------------------------------------
	vert_shader2 = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vert_shader2, 1, &vertex_shader2, NULL);
	glCompileShader(vert_shader2);
	frag_shader2 = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(frag_shader2, 1, &fragment_shader2, NULL);
	glCompileShader(frag_shader2);
	shader_programme2 = glCreateProgram();
	glAttachShader(shader_programme2, frag_shader2);
	glAttachShader(shader_programme2, vert_shader2);
	glLinkProgram(shader_programme2);
	//------------------------------------------------------------------------------------

	/* this loop clears the drawing surface, then draws the geometry described
			by the VAO onto the drawing surface. we 'poll events' to see if the window
			was closed, etc. finally, we 'swap the buffers' which displays our drawing
			surface onto the view area. we use a double-buffering system which means
			that we have a 'currently displayed' surface, and 'currently being drawn'
			surface. hence the 'swap' idea. in a single-buffering system we would see
			stuff being drawn one-after-the-other */
	while ( !glfwWindowShouldClose( window ) ) {
		/* wipe the drawing surface clear */
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glUseProgram(shader_programme);
		glBindVertexArray(vao);
		/* draw points 0-3 from the currently bound VAO with current in-use shader */
		glDrawArrays(GL_LINES, 0, totalCount);
	//------------------------------------------------------------------------------------	
		glUseProgram(shader_programme2);
		glBindVertexArray(vao2);
		glDrawArrays(GL_TRIANGLES, 0, numPoints);
	//------------------------------------------------------------------------------------
		
		/* update other events like input handling */
		glfwPollEvents();
		/* put the stuff we've been drawing onto the display */
		glfwSwapBuffers( window );
	}

	/* close GL context and any other GLFW resources */
	glfwTerminate();
	return 0;
}