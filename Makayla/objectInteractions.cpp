//g++ -o Makayla.exe objectInteractions.cpp libglew32.dll.a libglfw3dll.a -I include -lOpenGL32 -L ./ -lglew32 -lglfw3


/******************************************************************************\
| OpenGL 4 Example Code.                                                       |
| Accompanies written series "Anton's OpenGL 4 Tutorials"                      |
| Email: anton at antongerdelan dot net                                        |
| First version 27 Jan 2014                                                    |
| Copyright Dr Anton Gerdelan, Trinity College Dublin, Ireland.                |
| See individual libraries for separate legal notices                          |
|******************************************************************************|
| "Hello Triangle". Just the basics.                                           |
| If you're on Apple un-comment the version number code at the beginning. It   |
| will give you the latest, even if you say 3.2!                               |
| This uses the libraries GLEW and GLFW3 to start GL. Download and compile     |
| these first. Linking them might be a pain, but you'll need to master this.   |
\******************************************************************************/
#include <GL/glew.h>		/* include GLEW and new version of GL on Windows */
#include <GLFW/glfw3.h> /* GLFW helper library */
#include <stdio.h>
#include <math.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
using namespace std;


void multiply(GLfloat matrix1[], GLfloat matrix2[], GLfloat result[]){
	//GLfloat* result = new float[16];
	result[0] = (matrix1[0]*matrix2[0])+(matrix1[4]*matrix2[1])+(matrix1[8]*matrix2[2])+(matrix1[12]*matrix2[3]);
	result[4] = (matrix1[0]*matrix2[4])+(matrix1[4]*matrix2[5])+(matrix1[8]*matrix2[6])+(matrix1[12]*matrix2[7]);
	result[8] = (matrix1[0]*matrix2[8])+(matrix1[4]*matrix2[9])+(matrix1[8]*matrix2[10])+(matrix1[12]*matrix2[11]);
	result[12] = (matrix1[0]*matrix2[12])+(matrix1[4]*matrix2[13])+(matrix1[8]*matrix2[14])+(matrix1[12]*matrix2[15]);
	result[1] = (matrix1[1]*matrix2[0])+(matrix1[5]*matrix2[1])+(matrix1[9]*matrix2[2])+(matrix1[13]*matrix2[3]);
	result[5] = (matrix1[1]*matrix2[4])+(matrix1[5]*matrix2[5])+(matrix1[9]*matrix2[6])+(matrix1[13]*matrix2[7]);
	result[9] = (matrix1[1]*matrix2[8])+(matrix1[5]*matrix2[9])+(matrix1[9]*matrix2[10])+(matrix1[13]*matrix2[11]);
	result[13] = (matrix1[1]*matrix2[12])+(matrix1[5]*matrix2[13])+(matrix1[9]*matrix2[14])+(matrix1[13]*matrix2[15]);
	result[2] = (matrix1[2]*matrix2[0])+(matrix1[6]*matrix2[1])+(matrix1[10]*matrix2[2])+(matrix1[14]*matrix2[3]);
	result[6] = (matrix1[2]*matrix2[4])+(matrix1[6]*matrix2[5])+(matrix1[10]*matrix2[6])+(matrix1[14]*matrix2[7]);
	result[10] = (matrix1[2]*matrix2[8])+(matrix1[6]*matrix2[9])+(matrix1[10]*matrix2[10])+(matrix1[14]*matrix2[11]);
	result[14] = (matrix1[2]*matrix2[12])+(matrix1[6]*matrix2[13])+(matrix1[10]*matrix2[14])+(matrix1[14]*matrix2[15]);
	result[3] = (matrix1[3]*matrix2[0])+(matrix1[7]*matrix2[1])+(matrix1[11]*matrix2[2])+(matrix1[15]*matrix2[3]);
	result[7] = (matrix1[3]*matrix2[4])+(matrix1[7]*matrix2[5])+(matrix1[11]*matrix2[6])+(matrix1[15]*matrix2[7]);
	result[11] = (matrix1[3]*matrix2[8])+(matrix1[7]*matrix2[9])+(matrix1[11]*matrix2[10])+(matrix1[15]*matrix2[11]);
	result[15] = (matrix1[3]*matrix2[12])+(matrix1[7]*matrix2[13])+(matrix1[11]*matrix2[14])+(matrix1[15]*matrix2[15]);
	//return result;
}

float x,y,z =0;
float sx = 1;
float sy = 1;
float sz = 1;
float d = 1.55;
float rx,ry,rz = 0;
float fov = 67 * 3.14159 /180.0;
float aspect = 1.0;
float near = 0.01;
float far = 100.0;
float range = tan(fov*0.5)*near;
float Sx = (2*near)/((range*aspect)+(range*aspect));
float Sy = near/range;
float Sz = -(far+near)/(far-near);
float Pz = -(2*far*near)/(far-near);

GLfloat proMat[] = {Sx, 0.0f, 0.0f, 0.0f,
					0.0f, Sy, 0.0f, 0.0f,
					0.0f, 0.0f, Sz, -1.0,
					0.0f, 0.0f, Pz, 1.0f};
					
GLfloat orthagonalPro[] =  {1,0,0,0,
							0,1,0,0,
							0,0,0,0,
							0,0,0,1};
					
GLfloat lookAt[] = {1.0f, 0.0f, -0.0f, 0.0f,
					0.0f, 1.0f, -0.0f, 0.0f,
					0.0f, 0.0f, 1.0f, 0.0f,
					-0.0f, -0.0f, -1.0f, 1.0f};

GLfloat view[] =  { 1,0,0,0,
					0,1,0,0,
					0,0,1,0,
					0,0,0,1};
GLfloat* viewResult = new float[16];
GLfloat* viewMat = new float[16];

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){
	
	//cout<<"\n key \n";
	if(key == GLFW_KEY_Q && action == GLFW_PRESS){
		//cout<<"Q pressed\n";
		x -= 0.1;
	}
	else if(key == GLFW_KEY_W && action == GLFW_PRESS){
		//cout<<"W pressed\n";
		x += 0.1;
	}
	else if(key == GLFW_KEY_A && action == GLFW_PRESS){
		//cout<<"A pressed\n";
		y -= 0.1;
	}
	else if(key == GLFW_KEY_S && action == GLFW_PRESS){
		//cout<<"S pressed\n";
		y += 0.1;
	}
	else if(key == GLFW_KEY_E && action == GLFW_PRESS){
		//cout<<"E pressed\n";
		z -= 0.1;
	}
	else if(key == GLFW_KEY_D && action == GLFW_PRESS){
		//cout<<"D pressed\n";
		z += 0.1;
	}
	else if(key == GLFW_KEY_R && action == GLFW_PRESS){
		//cout<<"R pressed\n";
		sx -= 0.1;
	}
	else if(key == GLFW_KEY_T && action == GLFW_PRESS){
		//cout<<"T pressed\n";
		//cout<<"SX  "<<sx<<"\n";
		sx += 0.1;
	}
	else if(key == GLFW_KEY_F && action == GLFW_PRESS){
		//cout<<"F pressed\n";
		sy -= 0.1;
	}
	else if(key == GLFW_KEY_G && action == GLFW_PRESS){
		//cout<<"F pressed\n";
		//cout<<"SY  "<<sy<<"\n";
		sy += 0.1;
	}
	else if(key == GLFW_KEY_Y && action == GLFW_PRESS){
		//cout<<"Y pressed\n";
		sz -= 0.1;
	}
	else if(key == GLFW_KEY_H && action == GLFW_PRESS){
		//cout<<"H pressed\n";
		//cout<<"SZ  "<<sz<<"\n";
		sz += 0.1;
	}
	else if(key == GLFW_KEY_U && action == GLFW_PRESS){
		//cout<<"Y pressed\n";
		rx -= 1.0;
	}
	else if(key == GLFW_KEY_I && action == GLFW_PRESS){
		//cout<<"H pressed\n";
		rx += 1.0;
	}
	else if(key == GLFW_KEY_J && action == GLFW_PRESS){
		//cout<<"J pressed\n";
		ry -= 1.0;
	}
	else if(key == GLFW_KEY_K && action == GLFW_PRESS){
		//cout<<"K pressed\n";
		ry += 1.0;
	}
	else if(key == GLFW_KEY_O && action == GLFW_PRESS){
		//cout<<"J pressed\n";
		rz -= 1.0;
	}
	else if(key == GLFW_KEY_L && action == GLFW_PRESS){
		//cout<<"K pressed\n";
		rz += 1.0;
	}
	else if(key == GLFW_KEY_B && action == GLFW_PRESS){
		//cout<<"J pressed\n";
		d -= 0.1;
	}
	else if(key == GLFW_KEY_N && action == GLFW_PRESS){
		//cout<<"K pressed\n";
		//cout<<d<<"\n";
		d += 0.1;
	}
	else if(key == GLFW_KEY_Z && action == GLFW_PRESS){
		//cout<<"J pressed\n";
		fov -= 1.0;
	}
	else if(key == GLFW_KEY_X && action == GLFW_PRESS){
		//cout<<"K pressed\n";
		fov += 1.0;
	}
	else if(key == GLFW_KEY_C && action == GLFW_PRESS){
		//cout<<"J pressed\n";
		aspect -= 0.1;
	}
	else if(key == GLFW_KEY_V && action == GLFW_PRESS){
		//cout<<"K pressed\n";
		aspect += 0.1;
	}
	else if(key == GLFW_KEY_P && action == GLFW_PRESS){
		//cout<<"P pressed\n";
		//Perspective view
		viewMat=proMat;
		
	}
	else if(key == GLFW_KEY_M && action == GLFW_PRESS){
		//cout<<"M pressed\n";
		viewMat=orthagonalPro;
		
	}
}



int main() {
	GLFWwindow *window = NULL;
	const GLubyte *renderer;
	const GLubyte *version;
	GLuint vao;
	GLuint vbo;
	float numPoints = 12;
	
	
	/* geometry to use. these are 3 xyz points (9 floats total) to make a triangle */
	GLfloat points[] = { 0.0f,  0.5f, -0.5f, //front triangle point
						 0.5f, -0.5f, 0.0f, 
						-0.5f, -0.5f, 0.0f,
						
						0.0f,  0.5f, -0.5f, //right triangle
						0.5f, -0.5f,  0.0f, //point
						0.5f, -0.5f, -1.0f,
						
						-0.5f, -0.5f, -1.0f, //back triangle
						 0.5f, -0.5f, -1.0f, 
						 0.0f,  0.5f, -0.5f,
						
						  0.0f,  0.5f, -0.5f, //left triangle point
						 -0.5f, -0.5f,  0.0f, 
						 -0.5f, -0.5f, -1.0f
						 
						 /*-0.5f, -0.5f,  0.0f,//bottom triangle 1
						 -0.5f, -0.5f, -1.0f,
						  0.5f, -0.5f,  0.0f,
						  
						  0.5f, -0.5f,  0.0f, //bottom triangle 2
						  0.5f, -0.5f, -1.0f,
						 -0.5f, -0.5f, -1.0f*/};
								 
						
	GLfloat translate[]={1,0,0,0,
						 0,1,0,0,
						 0,0,1,0,
						 x,y,z,1};
	
	GLfloat scale[] =  {sx,0,0,0,
						0,sy,0,0,
						0,0,sz,0,
						0,0,0,1};
						
	GLfloat skew[] ={1,0,0,0,
					(cos(d)/sin(d)),1,0,0,
					 0,0,1,0,
					 0,0,0,1};
					 
	GLfloat rotateX[]= {1,0,0,0,
						0,cos(rx),sin(rx),0,
						0,-sin(rx),cos(rx),0,
						0,0,0,1};		
	
	GLfloat rotateY[]= {cos(ry),0,-sin(ry),0,
						0,1,0,0,
						sin(ry),0,cos(ry),0,
						0,0,0,1};
					 
	GLfloat rotateZ[]= {cos(rz),sin(rz),0,0,
						-sin(rz),cos(rz),0,0,
						0,0,1,0,
						0,0,0,1};
	
	GLfloat trans[] = { 1,0,0,0,
						0,1,0,0,
						0,0,1,0,
						0,0,0,1};
						
	GLfloat* result = new float[16];
	
	/* these are the strings of code for the shaders
	the vertex shader positions each vertex point */
	/*const char *vertex_shader = "#version 410\n"
								"in vec3 vp;"
								"void main () {"
								"	gl_Position = vec4 (vp, 1.0);"
								"}";*/
	const char *vertex_shader = "#version 410\n"
								"attribute vec3 vp;"
								"uniform mat4 model, view, proj;"
								"void main () {"
								"gl_Position = proj * view * model * vec4(vp, 1.0);"								
								//"gl_Position = view* vec4(vp, 1.0);"
								"}";
	
	/* the fragment shader colours each fragment (pixel-sized area of the
	triangle) */
	const char *fragment_shader = "#version 410\n"
								"out vec4 frag_colour;"
								"void main () {"
								"	frag_colour = vec4 (0.0, 1.0, 0.5, 1.0);"
								"}";
	/* GL shader objects for vertex and fragment shader [components] */
	GLuint vert_shader, frag_shader;
	/* GL shader programme object [combined, to link] */
	GLuint shader_programme;

	/* start GL context and O/S window using the GLFW helper library */
	if ( !glfwInit() ) {
		fprintf( stderr, "ERROR: could not start GLFW3\n" );
		return 1;
	}

/* We must specify 3.2 core if on Apple OS X -- other O/S can specify
anything here. I defined 'APPLE' in the makefile for OS X */
#ifdef APPLE
	glfwWindowHint( GLFW_CONTEXT_VERSION_MAJOR, 3 );
	glfwWindowHint( GLFW_CONTEXT_VERSION_MINOR, 2 );
	glfwWindowHint( GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE );
	glfwWindowHint( GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE );
#endif

	window = glfwCreateWindow( 640, 480, "Hello Triangle", NULL, NULL );
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
	glBufferData( GL_ARRAY_BUFFER, 3 * numPoints * sizeof( GLfloat ), points, GL_STATIC_DRAW );
	

	/* the vertex array object (VAO) is a little descriptor that defines which
	data from vertex buffer objects should be used as input variables to vertex
	shaders. in our case - use our only VBO, and say 'every three floats is a
	variable' */
	glGenVertexArrays( 1, &vao );
	glBindVertexArray( vao );
	// "attribute #0 should be enabled when this vao is bound"
	glEnableVertexAttribArray( 0 );
	// this VBO is already bound, but it's a good habit to explicitly specify which
	// VBO's data the following
	// vertex attribute pointer refers to
	glBindBuffer( GL_ARRAY_BUFFER, vbo );
	// "attribute #0 is created from every 3 variables in the above buffer, of type
	// float (i.e. make me vec3s)"
	glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, NULL );

	/* here we copy the shader strings into GL shaders, and compile them. we
	then create an executable shader 'program' and attach both of the compiled
			shaders. we link this, which matches the outputs of the vertex shader to
	the inputs of the fragment shader, etc. and it is then ready to use */
	vert_shader = glCreateShader( GL_VERTEX_SHADER );
	glShaderSource( vert_shader, 1, &vertex_shader, NULL );
	glCompileShader( vert_shader );
	frag_shader = glCreateShader( GL_FRAGMENT_SHADER );
	glShaderSource( frag_shader, 1, &fragment_shader, NULL );
	glCompileShader( frag_shader );
	shader_programme = glCreateProgram();
	glAttachShader( shader_programme, frag_shader );
	glAttachShader( shader_programme, vert_shader );
	glLinkProgram( shader_programme );

	/* this loop clears the drawing surface, then draws the geometry described
			by the VAO onto the drawing surface. we 'poll events' to see if the window
	was closed, etc. finally, we 'swap the buffers' which displays our drawing
			surface onto the view area. we use a double-buffering system which means
			that we have a 'currently displayed' surface, and 'currently being drawn'
			surface. hence the 'swap' idea. in a single-buffering system we would see
			stuff being drawn one-after-the-other */
	//multiply(view, proMat, viewResult);
	while ( !glfwWindowShouldClose( window ) ) {
		//Multiplying translation matrices
		multiply(trans, rotateX, result);
		multiply(result, rotateY, result);
		multiply(result, rotateZ, result);
		multiply(result, skew, result);
		multiply(result, scale, result);
		multiply(result, translate, result);
		
		//multiply(trans, scale, result);
		//multiply(result, translate, result);
		
		//Change projection matrices
		//multiply(view, proMat, viewResult);
	
		//View matrix info
		int trans_mat_location = glGetUniformLocation (shader_programme, "model");
		glUseProgram( shader_programme );
		glUniformMatrix4fv (trans_mat_location, 1, GL_FALSE, result);
		int view_mat_location = glGetUniformLocation (shader_programme, "view");
		glUseProgram( shader_programme );
		glUniformMatrix4fv (view_mat_location, 1, GL_FALSE, lookAt);
		int proj_mat_location = glGetUniformLocation (shader_programme, "proj");
		glUseProgram( shader_programme );
		glUniformMatrix4fv (proj_mat_location, 1, GL_FALSE, viewResult);
		//glUniformMatrix4fv (proj_mat_location, 1, GL_FALSE, proMat);
		
		/* wipe the drawing surface clear */
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		glUseProgram( shader_programme );
		glBindVertexArray( vao );
		/* draw points 0-3 from the currently bound VAO with current in-use shader */
		glDrawArrays( GL_TRIANGLES, 0, numPoints );
		/* update other events like input handling */
		glfwPollEvents();
		//
		glfwSetKeyCallback(window, key_callback);
		
		translate[12]=x;
		translate[13]=y;
		translate[14]=z;
		scale[0]=sx;
		scale[5]=sy;
		scale[10]=sz;
		rotateX[5]=cos(rx);
		rotateX[6]=sin(rx);
		rotateX[9]=-sin(rx);
		rotateX[10]=cos(rx);
		rotateY[0]=cos(ry);
		rotateY[2]=-sin(ry);
		rotateY[8]=sin(ry);
		rotateY[10]=cos(ry);
		rotateZ[0]=cos(rz);
		rotateZ[1]=sin(rz);
		rotateZ[4]=-sin(rz);
		rotateZ[5]=cos(rz);
		skew[4]=(cos(d)/sin(d));
		float range = tan(fov*0.5)*near;
		proMat[0]=(2*near)/((range*aspect)+(range*aspect));
		proMat[5]=near/range;
		proMat[10]=-(far+near)/(far-near);
		proMat[14]=-(2*far*near)/(far-near);
		
		multiply(view, viewMat, viewResult);
		//multiply(view, orthagonalPro, viewResult);
		
		/* put the stuff we've been drawing onto the display */
		glfwSwapBuffers( window );
	}

	/* close GL context and any other GLFW resources */
	glfwTerminate();
	return 0;
}
