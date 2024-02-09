#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <map>
#include <conio.h>
#include <unordered_set>
#include <queue>


#include <windows.h>
#include <glut.h>

using namespace std;

#define pi (2*acos(0.0))
#define eps 0.0000001
#define inf 10000.0

int cr;
int nump;
int T = 100;
int t;

FILE* fin;

class point
{
public:
    double x,y;
    point(double xx=0.0, double yy=0.0)
    {
        x=xx;
        y=yy;
    }
};

vector< vector<int> > G;
vector<point> locs;
vector<point> s;
vector<point> d;

double sqr_dist(point a, point b)
{
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

void generate_graph()
{
    int i, j;
    G.clear();

    for (i = 0; i < locs.size(); i++)
    {
        vector<int> tmp;
        tmp.clear();
        for (j = 0; j < locs.size(); j++)
        {
            if (j != i && sqr_dist(locs[i], locs[j]) <= cr * cr) tmp.push_back(j);
        }
        G.push_back(tmp);
    }
}


void input_file()
{
    fin = fopen("output.txt","r");

	fscanf(fin,"%d %d\n", &nump, &cr);

    s = vector<point> (nump, point());
    d = vector<point> (nump, point());

	for (int i = 0; i < nump; i++)
    {
        fscanf(fin,"%lf %lf\n", &s[i].x, &s[i].y);
        fscanf(fin,"%lf %lf\n", &d[i].x, &d[i].y);

        cout << s[i].x << " " <<s[i].y << endl;

    }

	fclose(fin);
}


void drawSquare(double x, double y, double s)
{
	glBegin(GL_QUADS);{
		glVertex3f( x + s, y + s, 0);
		glVertex3f( x + s, y - s, 0);
		glVertex3f( x - s, y - s, 0);
		glVertex3f( x - s, y + s, 0);
	}glEnd();
}


void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			//drawgrid = 1 - drawgrid;
			break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			//cameraHeight -= 3.0;
			break;
		case GLUT_KEY_UP:		// up arrow key
			//cameraHeight += 3.0;
			break;

		case GLUT_KEY_RIGHT:
			//cameraAngle += 0.03;
			break;
		case GLUT_KEY_LEFT:
			//cameraAngle -= 0.03;
			break;

		case GLUT_KEY_PAGE_UP:
			break;
		case GLUT_KEY_PAGE_DOWN:
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP

			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
	gluLookAt(0,0,0,	0,0,-1,	0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects


    int i,j;

    glColor3f(1, 0, 0);
    for (i = 0; i < locs.size(); i++)
    {
        drawSquare(locs[i].x, locs[i].y, 3);
    }

    glColor3f(0, 1, 0);
    for (i = 0; i < locs.size(); i++)
    {
        for (j = 0; j < G[i].size(); j++)
        {
            glBegin(GL_LINES);
            {
                glVertex3f(locs[i].x, locs[i].y, 0);
                glVertex3f(locs[G[i][j]].x, locs[G[i][j]].y, 0);
            }
            glEnd();
        }
    }




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(int val){
	//codes for any changes in Models, Camera
	locs.clear();
	locs = vector<point> (nump, point());
	//cout << locs.size() << endl;

	for (int i = 0; i < nump; i++)
    {
        locs[i].x = (d[i].x * t + s[i].x * (T - t)) / ((double) T);
        locs[i].y = (d[i].y * t + s[i].y * (T - t)) / ((double) T);
    }
    generate_graph();

    if (t < T) t++;

    glutPostRedisplay();
	glutTimerFunc(50, animate, 0);
}

void init(){
	//codes for initialization

	input_file();

    t = 0;



	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluOrtho2D(0, 500, 0, 500);
	//gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutTimerFunc(50, animate, 0);		//what you want to do in the idle time (when no drawing is occuring)

	//glutIdleFunc(animate);

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}



//void genobs()
//{
//    obsdata.clear();
//
//    obs o;
//    int i,ri,rii;
//    double rd;
//    srand(time(NULL));
//    int cnt=numobs;
//    while(cnt!=0)
//    {
//        ri=rand()%(rrange-osa-osd);
//        rii=rand()%1000;
//        o.xl=ri+(double)rii/1000.0;
//        ri=rand()%(2*osd)+(osa-osd);
//        rii=rand()%1000;
//        o.xh=o.xl+ri+(double)rii/1000.0;
//
//        ri=rand()%(rrange-osa-osd);
//        rii=rand()%1000;
//        o.yl=ri+(double)rii/1000.0;
//        ri=rand()%(2*osd)+(osa-osd);
//        rii=rand()%1000;
//        o.yh=o.yl+ri+(double)rii/1000.0;
//
//
//        if(doesintn(obsdata,o)==false)
//        {
//            obsdata.push_back(o);
//            cnt--;
//        }
//    }
//}






