//g++ -w -o Makayla.exe objTest.cpp libglew32.dll.a libglfw3dll.a -I include -lOpenGL32 -L ./ -lglew32 -lglfw3

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


int main(){
	GLFWwindow *window = NULL;
	const GLubyte *renderer;
	const GLubyte *version;
	GLuint vao;
	GLuint vbo;
	float rz = (90 * 3.14159) / 180;	//For the first rotation, so the trunk is 90 degrees from the bottom of the screen. Converts degrees to radians.
	float rx = (50 * 3.14159) / 180;
	GLfloat rotateZ[] = 
		{cos(rz),sin(rz),0,0,
		-sin(rz),cos(rz),0,0,
		0,0,1,0,
		0,0,0,1};
	//GLfloat scale[] =
	//	{sx,0,0,0,
	//	 0,sy,0,0,
	//	 0,0,sz,0,
	//	 0,0,0,1};
	GLfloat rotateX[] = 
		{1,0,0,0,
		 0,cos(rx),-sin(rx),0,
		 0,sin(rx),cos(rx),0,
		 0,0,0,1};
	
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
	const char *vertex_shader2 = "#version 410\n"
		"in vec3 vp;"
		"uniform mat4 rotateX;"
		"void main () {"
		"  gl_Position = rotateX * vec4 (vp, 1.0);"
		"}";
	/* the fragment shader colours each fragment (pixel-sized area of the
	triangle) */
	const char *fragment_shader2 = "#version 410\n"
		"out vec4 frag_colour;"
		"void main () {"
		"  frag_colour = vec4 (0.0, 1.0, 0.0, 1.0);"
		"}";
	/* GL shader objects for vertex and fragment shader [components] */
	GLuint vert_shader2, frag_shader2;
	/* GL shader program object [combined, to link] */
	GLuint shader_programme2;
	
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
	
	renderer = glGetString( GL_RENDERER ); /* get renderer string */
	version = glGetString( GL_VERSION );	 /* version as a string */
	printf( "Renderer: %s\n", renderer );
	printf( "OpenGL version supported %s\n", version );

	/* tell GL to only draw onto a pixel if the shape is closer to the viewer
	than anything already drawn at that pixel */
	glEnable( GL_DEPTH_TEST ); /* enable depth-testing */
	/* with LESS depth-testing interprets a smaller depth value as meaning "closer" */
	glDepthFunc( GL_LESS );
	
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
	
	vert_shader2 = glCreateShader( GL_VERTEX_SHADER );
	glShaderSource( vert_shader2, 1, &vertex_shader2, NULL );
	glCompileShader( vert_shader2 );
	frag_shader2 = glCreateShader( GL_FRAGMENT_SHADER );
	glShaderSource( frag_shader2, 1, &fragment_shader2, NULL );
	glCompileShader( frag_shader2 );
	shader_programme2 = glCreateProgram();
	glAttachShader( shader_programme2, frag_shader2 );
	glAttachShader( shader_programme2, vert_shader2 );
	glLinkProgram( shader_programme2 );
	
	while ( !glfwWindowShouldClose( window ) ) {
		/* wipe the drawing surface clear */
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
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