
//g++ -o Makayla.exe Fractal.cpp libglew32.dll.a libglfw3dll.a -I include -lOpenGL32 -L ./ -lglew32 -lglfw3

#include <GL/glew.h>		/* include GLEW and new version of GL on Windows */
#include <GLFW/glfw3.h>     /* GLFW helper library */
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stack>
#include <vector>
using namespace std;

/*void  multiply(GLfloat matrix1[], GLfloat matrix2[], GLfloat result[]){ // this is to multiply a 4x4 with another 4x4 
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
} */

void multiply(GLfloat matrix1[], GLfloat matrix2[], GLfloat result[]){   // this is to multiply a 4x4 and a 4x1 matrix together
	result[0] = (matrix1[0]*matrix2[0]) + (matrix1[4]*matrix2[1]) + (matrix1[8]*matrix2[2]) + (matrix1[12]*matrix2[3]);
	result[1] = (matrix1[1]*matrix2[0]) + (matrix1[5]*matrix2[1]) + (matrix1[9]*matrix2[2]) + (matrix1[13]*matrix2[3]);
	result[2] = (matrix1[2]*matrix2[0]) + (matrix1[6]*matrix2[1]) + (matrix1[10]*matrix2[2]) + (matrix1[14]*matrix2[3]);
	result[3] = (matrix1[3]*matrix2[0]) + (matrix1[7]*matrix2[1]) + (matrix1[11]*matrix2[2]) + (matrix1[15]*matrix2[3]);
}

string generatePattern(){
    int numIts = 1; // Number of iterations
    string pattern = "F"; //"[X]";    // Using F for the pattern 
    
    for (int i = 0; i < numIts; i++){
        string newPattern = ""; 
        for (int idx = 0; idx < pattern.length(); idx++){
            //cout << "char: " << pattern.substr(idx,1) << endl;
            if (pattern.substr(idx,1).compare("F") == 0) 
                newPattern += "F[-F][F][+F]"; 
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

int main() {
	
	//GLfloat matrix1[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
	//GLfloat matrix2[] = {2, 3, 4, 5};
	//GLfloat* result = new float[4];
	//multiply(matrix1, matrix2, result);
	//cout << matrix1[14] << " " << matrix2[3] << " " << result[2] << endl;
	/* string topofStack;
	stack<string> mystack;
	mystack.push("H");
	mystack.push("W");
	topofStack = mystack.top();
	mystack.pop();
	cout << topofStack << endl; */
	 
	GLFWwindow *window = NULL;
	const GLubyte *renderer;
	const GLubyte *version;
	GLuint vao;
	GLuint vbo;
	
	int count = 1;
	string pattern = generatePattern();
	cout << "pattern: " << pattern << endl;
	cout << endl;
	for (int idx = 0; idx < pattern.length(); idx++){
			if (pattern.substr(idx, 1).compare("F") == 0){
				count ++;
				//cout << "F" << endl;
			}
		}
		
	GLfloat points[count*4];
	GLfloat* result = new float[4];
	stack<float> PositionStack;
	stack<float> HeadingStack;
	float rz = (90 * 3.14159) / 180;
	int pointsCount = 0;
	GLfloat currentPosition[] = {-0.25f, -0.25f, 0.0f, 1.0f};
	GLfloat currentHeading[] = {0.0f, 0.5f, 0.0f, 0.0f};
	GLfloat rotateZ[] = 
		{cos(rz),sin(rz),0,0,
		-sin(rz),cos(rz),0,0,
		0,0,1,0,
		0,0,0,1};

	points[pointsCount] = currentPosition[0];
	cout << "Points 1: " << points[pointsCount] << endl;
	pointsCount++;
	points[pointsCount] = currentPosition[1];
	cout << "Points 2: " << points[pointsCount] << endl;
	pointsCount++;
	points[pointsCount] = currentPosition[2];
	cout << "Points 3: " << points[pointsCount] << endl;
	pointsCount++;
	points[pointsCount] = currentPosition[3];
	cout << "Points 4: " << points[pointsCount] << endl;
	pointsCount++;
	cout << endl;
	
	for (int idx = 0; idx < pattern.length(); idx++){
		
		
		if (pattern.substr(idx,1).compare("[") == 0){
			cout << "Printing currentPosition before push " << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			PositionStack.push(currentPosition[3]);
			PositionStack.push(currentPosition[2]);
			PositionStack.push(currentPosition[1]);
			PositionStack.push(currentPosition[0]);
			cout << "Printing currentPosition after push " << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			cout << endl;
			
			HeadingStack.push(currentHeading[3]);
			HeadingStack.push(currentHeading[2]);
			HeadingStack.push(currentHeading[1]);
			HeadingStack.push(currentHeading[0]);
		}
		
		else if (pattern.substr(idx, 1).compare("]") == 0){
			cout << "Printing currentPosition before pop " << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			currentPosition[0] = PositionStack.top();
			PositionStack.pop();
			currentPosition[1] = PositionStack.top();
			PositionStack.pop();
			currentPosition[2] = PositionStack.top();
			PositionStack.pop();
			currentPosition[3] = PositionStack.top();
			PositionStack.pop();
			cout << "Printing currentPosition after pop " << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			
			currentHeading[0] = HeadingStack.top();
			HeadingStack.pop();
			currentHeading[1] = HeadingStack.top();
			HeadingStack.pop();
			currentHeading[2] = HeadingStack.top();
			HeadingStack.pop();
			currentHeading[3] = HeadingStack.top();
			HeadingStack.pop();
			cout << endl;
		}
		
		else if (pattern.substr(idx, 1).compare("F") == 0){
			
			cout << " F before adding" << endl;
			cout << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			
			currentPosition[0] += currentHeading[0]*.02;
			currentPosition[1] += currentHeading[1]*.02;
			currentPosition[2] += currentHeading[2]*.02;
			currentPosition[3] += currentHeading[3];
			
			cout << "F after adding" << endl;
			cout << currentPosition[0] << " " << currentPosition[1] << " " << currentPosition[2] << " " << currentPosition[3] << endl;
			cout << endl;
			
			points[pointsCount] = currentPosition[0];
			cout << "Points 1: " << points[pointsCount] << endl;
			pointsCount++;
			points[pointsCount] = currentPosition[1];
			cout << "Points 2: " << points[pointsCount] << endl;
			pointsCount++;
			points[pointsCount] = currentPosition[2];
			cout << "Points 3: " << points[pointsCount] << endl;
			pointsCount++;
			points[pointsCount] = currentPosition[3];
			cout << "Points 4: " << points[pointsCount] << endl;
			pointsCount++;
			cout << endl;
		}
		
		else if (pattern.substr(idx, 1).compare("+") == 0){
			cout <<  "+ rz before " << rz << endl;
			float rz =-  ((10 * 3.14159) / 180);
			rotateZ[0] = cos(rz);
			rotateZ[1] = sin(rz);
			rotateZ[4] = -sin(rz);
			rotateZ[5] = cos(rz);
			cout <<  "+ rz after " << rz << endl;
			cout << endl;
			multiply(rotateZ, currentHeading, result);
			float magnitude = sqrt((result[0]*result[0]) + (result[1]*result[1]) + (result[2]*result[2]) + (result[3]*result[3]));
			currentHeading[0] = result[0] / magnitude;
			currentHeading[1] = result[1] / magnitude;
			currentHeading[2] = result[2] / magnitude;
			currentHeading[3] = result[3] / magnitude;
			cout << "Current Heading after Normalization: " << currentHeading[0] << " " << currentHeading[1] << " " << currentHeading[2] << " " << currentHeading[3] << endl;
			cout << endl;
		}
		
		else if (pattern.substr(idx, 1).compare("-") == 0){
			cout <<  "- rz before " << rz << endl;
			float rz =+ ((10 * 3.14159) / 180);
			rotateZ[0] = cos(rz);
			rotateZ[1] = sin(rz);
			rotateZ[4] = -sin(rz);
			rotateZ[5] = cos(rz);
			cout <<  "- rz after " << rz << endl;
			cout << endl;
			multiply(rotateZ, currentHeading, result);	
			float magnitude = sqrt((result[0]*result[0]) + (result[1]*result[1]) + (result[2]*result[2]) + (result[3]*result[3]));
			currentHeading[0] = result[0] / magnitude;
			currentHeading[1] = result[1] / magnitude;
			currentHeading[2] = result[2] / magnitude;
			currentHeading[3] = result[3] / magnitude;
			cout << "Current Heading after Normalization: " << currentHeading[0] << " " << currentHeading[1] << " " << currentHeading[2] << " " << currentHeading[3] << endl;
			cout << endl;
		}
	}

	/* geometry to use. these are 3 xyz points (9 floats total) to make a triangle */
	//GLfloat points[] = { 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, -1.0f, -1.0f, 0.0f };
	/* these are the strings of code for the shaders
	the vertex shader positions each vertex point */
	const char *vertex_shader = "#version 410\n"
		"in vec3 vp;"
		"void main () {"
		"  gl_Position = vec4 (vp, 1.0);"
		"}";
	/* the fragment shader colours each fragment (pixel-sized area of the
	triangle) */
	const char *fragment_shader = "#version 410\n"
		"out vec4 frag_colour;"
		"void main () {"
		"  frag_colour = vec4 (0.0, 1.0, 0.5, 1.0);"
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

	window = glfwCreateWindow( 640, 480, "Fractal Window", NULL, NULL );
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
	glBufferData( GL_ARRAY_BUFFER, 9 * sizeof( GLfloat ), points, GL_STATIC_DRAW );

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
	glPointSize(5.0);

	/* this loop clears the drawing surface, then draws the geometry described
			by the VAO onto the drawing surface. we 'poll events' to see if the window
			was closed, etc. finally, we 'swap the buffers' which displays our drawing
			surface onto the view area. we use a double-buffering system which means
			that we have a 'currently displayed' surface, and 'currently being drawn'
			surface. hence the 'swap' idea. in a single-buffering system we would see
			stuff being drawn one-after-the-other */
	while ( !glfwWindowShouldClose( window ) ) {
		/* wipe the drawing surface clear */
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		glUseProgram( shader_programme );
		glBindVertexArray( vao );
		/* draw points 0-3 from the currently bound VAO with current in-use shader */
		glDrawArrays( GL_POINTS , 0, count);
		/* update other events like input handling */
		glfwPollEvents();
		/* put the stuff we've been drawing onto the display */
		glfwSwapBuffers( window );
	}

	/* close GL context and any other GLFW resources */
	glfwTerminate();
	return 0;
}