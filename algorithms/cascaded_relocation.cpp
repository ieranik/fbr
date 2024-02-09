#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <time.h>
#include <ctime>
#include <stdlib.h>
#include <map>
#include <conio.h>
#include <unordered_set>
#include <queue>
#include <sys/time.h>


#include <windows.h>
#include <glut.h>

using namespace std;

#define pi (2*acos(0.0))
#define eps 0.0000001
#define inf 10000.0
#define numtests 1


FILE* fin;
FILE* fout;

double cr;
int numpoints;


double max2(double a,double b)
{
    return (a>b)?a:b;
}

double min2(double a,double b)
{
    return (a<b)?a:b;
}

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
vector<point> allpoints;
vector<point> saved;
vector<point> start;
vector<bool> isboun;
int fn;
int bci = 0;
point bcp = point(0, 0);


void printlocs()
{
    for (int i = 0; i < locs.size(); i++) cout << locs[i].x << " " << locs[i].y << endl;
    cout << endl;
}


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
            if (j != i && sqr_dist(locs[i], locs[j]) <= cr * cr + eps) tmp.push_back(j);
        }
        G.push_back(tmp);
    }
}


int first_edge(int cur)
{
    double amin = 91.0;
    int idx = -1;
    for (int i = 0; i < G[cur].size(); i++)
    {
        double angle = atan2(locs[G[cur][i]].y - locs[cur].y, locs[G[cur][i]].x - locs[cur].x);
        angle *= (180 / pi);
        if (angle < amin)
        {
            amin = angle;
            idx = G[cur][i];
        }
    }
    //cout << idx << endl;
    return idx;
}


int next_edge(int cur, int prev)
{
    double refangle = atan2(locs[cur].y - locs[prev].y, locs[cur].x - locs[prev].x);
    refangle *= (180 / pi);
    double amin = 200.0;
    int idx = -1;
    for (int i = 0; i < G[cur].size(); i++)
    {
        if (G[cur][i] == prev) continue;
        double angle = atan2(locs[G[cur][i]].y - locs[cur].y, locs[G[cur][i]].x - locs[cur].x);
        angle *= (180 / pi);
        angle = angle - refangle;
        if (angle > 180) angle -= 360.0;
        else if (angle < -180) angle += 360.0;

        if (angle < amin)
        {
            amin = angle;
            idx = G[cur][i];
        }
    }
    return idx;
}

void boundary_points()
{
    isboun.clear();
    isboun = vector<bool>(G.size(), false);
    double xmin = 1000;
    int mi = -1;
    for (int i = 0; i < locs.size(); i++)
    {
        if (locs[i].x < xmin)
        {
            xmin = locs[i].x;
            mi = i;
        }
    }
    //cout << locs[0].x << endl;
    //cout << mi << endl;
    int prev = mi;
    int cur = first_edge(prev);
    //cout << cur << endl;
    isboun[cur] = true;
    isboun[prev] = true;
    while (true)
    {
        int next = next_edge(cur, prev);

        isboun[next] = true;
        if (next == mi) break;
        prev = cur;
        cur = next;
    }
}



point repos2(int p, int c)
{
    double len = sqrt(sqr_dist(locs[p], locs[c]));
    if (len < cr) return locs[c];
    double dx = locs[c].x - locs[p].x;
    double dy = locs[c].y - locs[p].y;
    dx *= (cr / len);
    dy *= (cr / len);
    return point(locs[p].x + dx, locs[p].y +  dy);
}

bool all_con(point p)
{
    for (int i = 0; i < G[fn].size(); i++)
    {
        double dis = sqr_dist(p, locs[G[fn][i]]);
        if (dis > cr * cr) return false;
    }
    return true;
}

point subdivide(int d, int s, double f)
{
    double dx = locs[d].x - locs[s].x;
    double dy = locs[d].y - locs[s].y;
    dx *= f;
    dy *= f;
    return point(locs[s].x + dx, locs[s].y +  dy);
}

point go_close(int d, int s)
{
    point tp;
    int cnt = 15;
    double f = 0.25;
    double rt = 0.5;
    tp = subdivide(d, s, rt);
    while (cnt--)
    {
        if (all_con(subdivide(d, s, rt)) == true) rt -= f;
        else rt += f;

        f /= 2.0;
    }
    tp = subdivide(d, s, rt);

    return tp;
}

void find_bc()
{
    double mindis = 1000000;

    for (int i = 0; i < G[fn].size(); i++)
    {
        double dis = sqr_dist(locs[fn], locs[G[fn][i]]);
        if (dis < mindis)
        {
            mindis = dis;
            bci = G[fn][i];
        }
    }

    bcp = go_close(fn, bci);



}

void repos(int p, int c)
{
    double len = sqrt(sqr_dist(locs[p], locs[c]));
    double dx = locs[c].x - locs[p].x;
    double dy = locs[c].y - locs[p].y;
    dx *= (cr / len);
    dy *= (cr / len);
    locs[c].x = locs[p].x + dx;
    locs[c].y = locs[p].y + dy;
}


void cascade(int r, point rloc)
{
    locs[r] = rloc;

    bool *vis = new bool[G.size()];
    for (int i = 0; i < G.size(); i++) vis[i] = false;

    queue<int> Q;

    Q.push(r);
    vis[r] = true;

    while(Q.empty() == false)
    {
        int t = Q.front();
        Q.pop();
        for (int i = 0; i < G[t].size(); i++)
        {
            int v = G[t][i];
            if (vis[v] == true) continue;
            vis[v] = true;
            if (sqr_dist(locs[t], locs[v]) > cr * cr) repos(t, v);

            Q.push(v);
        }
    }
}



void input_data()
{
    fin = fopen("errortest.txt","r");
	fscanf(fin, " %d", &numpoints);
	fscanf(fin, " %lf", &cr);

	allpoints.clear();
	for (int i = 0; i < numpoints * numtests; i++)
    {
        point tmp;
        fscanf(fin, " %lf %lf", &tmp.x, &tmp.y);
        allpoints.push_back(tmp);
    }

    //cout << " ok" << locs[1].x << endl;
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
//			find_bc();
//            cascade(bci, bcp);
//            generate_graph();
//
//            cout << fn << " " << bci << endl;
//            cout << sqr_dist(locs[11], locs[13]) << endl;
			break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
            //cascade(0, point(locs[0].x, locs[0].y - 5));
            //generate_graph();
			break;
		case GLUT_KEY_UP:		// up arrow key
			//cameraHeight += 3.0;
			//cascade(0, point(locs[0].x, locs[0].y + 5));
            //generate_graph();
			break;

		case GLUT_KEY_RIGHT:
		    //cascade(0, point(locs[0].x + 5, locs[0].y));
            //generate_graph();
			//cameraAngle += 0.03;
			break;
		case GLUT_KEY_LEFT:
		    //cascade(0, point(locs[0].x - 5, locs[0].y));
            //generate_graph();
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

                find_bc();
                cascade(bci, bcp);
                generate_graph();

                cout << fn << " " << bci << endl;


			}
			break;

		case GLUT_RIGHT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP

                locs = saved;
                generate_graph();
                fn++;

			}
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
	/ set-up camera f
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
	gluLookAt(0,0,0,    0,0,-1, 0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects


    int i,j;


    for (i = 0; i < locs.size(); i++)
    {
        glColor3f(0, 0, 1);
        drawSquare(locs[i].x, locs[i].y, 3);
    }

    glColor3f(1, 0, 0);
    drawSquare(locs[fn].x, locs[fn].y, 5);

    glColor3f(1, 1, 0);
    drawSquare(locs[bci].x, locs[bci].y, 5);

    glColor3f(0, 1, 1);
    drawSquare(bcp.x + 2, bcp.y + 2, 5);



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

    glColor3f(1, 0, 1);

    for (i = 0; i < start.size() - 1; i++)
    {
        drawSquare(start[i].x, start[i].y, 5);
    }






	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}



void animate(int val){
	//codes for any changes in Models, Camera



    glutPostRedisplay();
	glutTimerFunc(50, animate, 0);
}



void init(){
	//codes for initialization


    double sum = 0;

	input_data();

	fout = fopen("out256s.txt","w");

	srand(time(NULL));

	struct timeval tps;
    struct timeval tpe;
    long int st;
    long int et;

    gettimeofday(&tps, NULL);
    st = tps.tv_sec * 1000000 + tps.tv_usec ;


    for (int i = 0; i < numtests; i++)
    {
        locs = vector<point>(allpoints.begin() + i * numpoints, allpoints.begin() + (i + 1) * numpoints);
        start = locs;

        fn = numpoints - 1;
        generate_graph();
        cout << "here1" << endl;
        boundary_points();
        cout << "here2" << endl;
        find_bc();
        cout << bci << endl;
        cout << bcp.x << " " << bcp.y << endl;
        cascade(bci, bcp);

        double maxdis = -1;
        for (int i = 0; i < locs.size(); i++)
        {
            if (sqr_dist(locs[i], start[i]) > maxdis)
                maxdis = sqr_dist(locs[i], start[i]);
        }

        fprintf(fout, "%f\n", sqrt(maxdis));
        sum += sqrt(maxdis);
        //cout << i << " " << sqrt(maxdis) << endl;


    }

    gettimeofday(&tpe, NULL);
    et = tpe.tv_sec * 1000000 + tpe.tv_usec ;

    long int oldt = (et - st);
    cout << "Time: " << oldt / 100 << endl;
    cout << "Distance: " << sum / (100 * cr) << endl;

    fclose(fout);














    //generate_graph();



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


