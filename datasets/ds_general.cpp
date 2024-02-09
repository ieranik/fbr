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
#define cr 200.0
#define dsize 16


FILE* fin;
FILE* fout;


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

class bcnode
{
public:
    bool iscv;
    vector<int> nodelist;
    int depth;
    bcnode()
    {
        iscv = false;
        nodelist.clear();
        depth = 0;
    }
};

class amnode
{
public:
    int src;
    int dst;
    double weight;
    amnode()
    {
        src = dst = -1;
        weight = inf;
    }
};

void printbcn(bcnode n)
{
    cout << n.iscv << " ";
    for (int i = 0; i < n.nodelist.size(); i++) cout << n.nodelist[i] << " ";
    cout << endl;
}




vector< vector<int> > G;
vector<point> locs;
vector<bool> isap;
vector<int> leaves;
int numcv;
bool is2con;
bool iscon;


vector<bcnode> bcl;
vector<int> blocknum;
vector<int> vtb;
vector< pair<int, int> > bridgelist;

vector< vector<int> > T;
vector< vector< pair<int, double> > > DT;
vector< pair<int, int> > edges;
vector< pair<int, int> > add_edges;

void printlocs()
{
    for (int i = 0; i < locs.size(); i++) cout << locs[i].x << " " << locs[i].y << endl;
    cout << endl;
}


double sqr_dist(point a, point b)
{
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

vector<int> removedup(vector<int> v)
{
    vector<int> tmp = vector<int> (G.size(), 0);
    for (int i = 0; i < v.size(); i++) tmp[v[i]] = 1;
    vector<int> ret;
    for (int i = 0; i < tmp.size(); i++) if (tmp[i] == 1) ret.push_back(i);
    return ret;
}

vector<int> removedup1(vector<int> v)
{
    if (v.size() == 1 && isap[v[0]] == true) return v;
    vector<int> tmp = vector<int> (G.size(), 0);
    for (int i = 0; i < v.size(); i++) if (isap[v[i]] == false) tmp[v[i]] = 1;
    vector<int> ret;
    for (int i = 0; i < tmp.size(); i++) if (tmp[i] == 1) ret.push_back(i);
    return ret;
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





void get_ap(int u, bool visited[], int disc[],  int low[], int parent[], bool ap[])
{
    // A static variable is used for simplicity, we can avoid use of static
    // variable by passing a pointer.
    static int time = 0;

    // Count of children in DFS Tree
    int children = 0;

    // Mark the current node as visited
    visited[u] = true;

    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;

    // Go through all vertices aadjacent to this
    for (int i = 0; i < G[u].size(); i++)
    {
        int v = G[u][i];  // v is current adjacent of u

        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (!visited[v])
        {
            children++;
            parent[v] = u;
            get_ap(v, visited, disc, low, parent, ap);

            // Check if the subtree rooted with v has a connection to
            // one of the ancestors of u
            low[u]  = min(low[u], low[v]);

            if (low[v] > disc[u]) bridgelist.push_back(pair<int, int> (u, v));

            // u is an articulation point in following cases

            // (1) u is root of DFS tree and has two or more chilren.
            if (parent[u] == -1 && children > 1)
               ap[u] = true;

            // (2) If u is not root and low value of one of its child is more
            // than discovery value of u.
            if (parent[u] != -1 && low[v] >= disc[u])
               ap[u] = true;
        }

        // Update low value of u for parent function calls.
        else if (v != parent[u])
            low[u]  = min(low[u], disc[v]);
    }
}

void art_points()
{
    bool *visited = new bool[G.size()];
    int *disc = new int[G.size()];
    int *low = new int[G.size()];
    int *parent = new int[G.size()];
    bool *ap = new bool[G.size()]; // To store articulation points

    // Initialize parent and visited, and ap(articulation point) arrays
    for (int i = 0; i < G.size(); i++)
    {
        parent[i] = -1;
        visited[i] = false;
        ap[i] = false;
    }

    bridgelist.clear();
    // Call the recursive helper function to find articulation points
    // in DFS tree rooted with vertex 'i'
    for (int i = 0; i < G.size(); i++)
        if (visited[i] == false)
            get_ap(i, visited, disc, low, parent, ap);

    // Now ap[] contains articulation points, print them


    isap = vector<bool>(G.size(), false);
    is2con = true;
    for (int i = 0; i < G.size(); i++)
    {
        isap[i] = ap[i];
        if (ap[i] == true) is2con = false;
    }

}



int dfsvis(int u, bool visited[])
{
    visited[u] = true;
    int cnt = 1;

    for (int i = 0; i < G[u].size(); i++)
    {
        int v = G[u][i];
        if (!visited[v]) cnt += dfsvis(v, visited);
    }
    return cnt;
}



void dfs()
{
    bool *visited = new bool[G.size()];

    for (int i = 0; i < G.size(); i++) visited[i] = false;

    int cnt = dfsvis(0, visited);

    if (cnt < (int)G.size()) iscon = false;
    else iscon = true;
}


vector<int> get_bn(int u, bool visited[])
{
    visited[u] = true;
    vector<int> ret;
    ret.clear();
    ret.push_back(u);
    if (isap[u] == true) return ret;



    for (int i = 0; i < G[u].size(); i++)
    {
        int v = G[u][i];  // v is current adjacent of u

        if (!visited[v] || isap[v] == true)
        {
            vector<int> tmp = get_bn(v, visited);
            for (int j = 0; j < tmp.size(); j++) ret.push_back(tmp[j]);
        }
    }
    return ret;
}

void bnums()
{
    bool *visited = new bool[G.size()];

    for (int i = 0; i < G.size(); i++)
    {
        visited[i] = false;
    }

    for (int i = 0; i < G.size(); i++)
    {
        if (visited[i] == false && isap[i] == false)
        {
            vector<int> ret = get_bn(i, visited);
            if (ret.size() != 0)
            {
                bcnode tmp;
                tmp.iscv = false;
                tmp.nodelist = ret;
                bcl.push_back(tmp);
            }
        }
    }
}


void create_bcl()
{
    bcl.clear();

	for (int i = 0; i < G.size(); i++)
    {
        if (isap[i] == true)
        {
            bcnode tmp;
            tmp.iscv = true;
            vector<int> vt;
            vt.clear();
            vt.push_back(i);
            tmp.nodelist = vt;
            bcl.push_back(tmp);
        }
    }

    bnums();

    for (int i = 0; i < bridgelist.size(); i++)
    {
        if (isap[bridgelist[i].first] == true && isap[bridgelist[i].second] == true)
        {
            bcnode tmp;
            tmp.iscv = false;
            vector<int> vt;
            vt.clear();
            vt.push_back(bridgelist[i].first);
            vt.push_back(bridgelist[i].second);
            tmp.nodelist = vt;
            bcl.push_back(tmp);
        }
    }

    for (int i = 0; i < bcl.size(); i++) bcl[i].nodelist = removedup(bcl[i].nodelist);

    cout << "BC nodes" << endl;
    for (int i = 0; i < bcl.size(); i++) printbcn(bcl[i]);

}


void create_tree()
{
    vtb.clear();
    vtb = vector<int> (G.size(), 0);
    for (int i = 0; i < bcl.size(); i++)
    {
        if (bcl[i].iscv == true)
        {
            vtb[bcl[i].nodelist[0]] = i;
        }
        else
        {
            for (int j = 0; j < bcl[i].nodelist.size(); j++)
            {
                int v = bcl[i].nodelist[j];
                if (isap[v] == false) vtb[v] = i;
            }
        }
    }



    T.clear();
    for (int i = 0; i < bcl.size(); i++)
    {
        vector<int> tmp;
        tmp.clear();
        T.push_back(tmp);
    }
    for (int i = 0; i < bcl.size(); i++)
    {
        if (bcl[i].iscv == false)
        {
            for (int j = 0; j < bcl[i].nodelist.size(); j++)
            {
                int v = bcl[i].nodelist[j];
                if (isap[v] == true)
                {
                    int u = vtb[v];
                    T[i].push_back(u);
                    T[u].push_back(i);
                }
            }
        }
    }
    cout << endl;

    cout << "Tree edges" << endl;
    for (int i = 0; i < T.size(); i++)
    {
        cout << i << ": ";
        for (int j = 0; j < T[i].size(); j++)
        {
            cout << T[i][j] << " ";
        }
        cout << endl;
    }
}


void create_dtree()
{
    DT.clear();
    for (int i = 0; i < bcl.size(); i++)
    {
        vector<pair<int, double> > tmp;
        tmp.clear();
        DT.push_back(tmp);
    }
    queue<int> Q;


    int maxx = -1;
    numcv = -1;
    for (int i = 0; i < bcl.size(); i++)
    {
        //cout << bcl[i].nodelist.size() << endl;
        if ((int)bcl[i].nodelist.size() > maxx)
        {
            maxx = (int)bcl[i].nodelist.size();
            numcv = i;
        }
    }

    //cout << numcv << endl;

//    numcv = 0;
//    while (bcl[numcv].iscv == true) numcv++;

    bool *visited = new bool[T.size()];
    for (int i = 0; i < T.size(); i++) visited[i] = false;

    bcl[numcv].depth = 0;
    Q.push(numcv);
    visited[numcv] = true;

    while(Q.empty() == false)
    {
        int t = Q.front();
        Q.pop();
        for (int i = 0; i < T[t].size(); i++)
        {
            int v = T[t][i];
            if (visited[v] == true) continue;
            visited[v] = true;
            Q.push(v);
            bcl[v].depth = bcl[t].depth + 1;
            DT[v].push_back(pair<int, double> (t, 0));
            //am[min(v,t)][max(v,t)].weight = 0;

        }
    }
    cout << endl;

    cout << "Parent" << endl;

    for (int i = 0; i < DT.size(); i++)
    {
        cout << i << ": ";
        for (int j = 0; j < DT[i].size(); j++)
        {
            cout << DT[i][j].first << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Depth" << endl;

    for (int i = 0; i < bcl.size(); i++) cout << i << " " << bcl[i].depth << endl;

    cout << endl;

}


void find_leaves()
{
    leaves.clear();
    vector<int> indeg = vector<int>(DT.size(), 0);
    for (int i = 0; i < DT.size(); i++)
    {
        if ((int)DT[i].size() > 0) indeg[DT[i][0].first]++;
    }
    //cout << "leaves" << endl;
    for (int i = 0; i < indeg.size(); i++)
    {
        if (indeg[i] == 0)
        {
            //cout << i << endl;
            leaves.push_back(i);
        }
    }
}


point scale(point p)
{
    double len = sqrt(p.x * p.x + p.y * p.y);
    p.x = p.x * ((len - cr) / cr);
    p.y = p.y * ((len - cr) / cr);
    return p;
}

void process_leaves()
{
    for (int i = 0; i < leaves.size(); i++)
    {
        //cout << "loop" << endl;
        int lb = leaves[i];
        int cv = DT[lb][0].first;

        int pb = DT[cv][0].first;
        int cvn = bcl[cv].nodelist[0];

        //cout << lb << " " << cv << " " << pb << " " << cvn << endl;

        double mindis = 1000000.0;
        int lmin = -1;
        int pmin = -1;
        for (int j = 0; j < bcl[lb].nodelist.size(); j++)
        {
            for (int k = 0; k < bcl[pb].nodelist.size(); k++)
            {
                if (bcl[lb].nodelist[j] != cvn && bcl[pb].nodelist[k] != cvn)
                {
                    int ln = bcl[lb].nodelist[j];
                    int pn = bcl[pb].nodelist[k];
                    if (sqr_dist(locs[ln], locs[pn]) < mindis)
                    {
                        mindis = sqr_dist(locs[ln], locs[pn]);
                        lmin = ln;
                        pmin = pn;
                    }
                }
            }
        }
        double dx = locs[pmin].x - locs[lmin].x;
        double dy = locs[pmin].y - locs[lmin].y;
        point tp(dx, dy);
        tp = scale(tp);
        dx = tp.x;
        dy = tp.y;
        for (int j = 0; j < bcl[lb].nodelist.size(); j++)
        {
            if (bcl[lb].nodelist[j] != cvn)
            {
                 locs[bcl[lb].nodelist[j]].x += dx;
                 locs[bcl[lb].nodelist[j]].y += dy;
            }
        }
    }
    return;
}

void add_to_file()
{
	for (int i = 0; i < locs.size(); i++)
    {
        fprintf(fout,"%f %f\n", locs[i].x, locs[i].y);
    }
}


void input_data()
{
    fin = fopen("points0.txt","r");
	int numpoints;
	fscanf(fin, " %d", &numpoints);
	locs.clear();
	for (int i = 0; i < numpoints; i++)
    {
        point tmp;
        fscanf(fin, " %lf %lf", &tmp.x, &tmp.y);
        locs.push_back(tmp);
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
			if(state == GLUT_DOWN)
            {
                point p;
                int ri, rii;

                while (true)
                {

                    ri = rand() % 500;
                    rii = rand() % 100;
                    double x = ri + (double)rii / 100.0;
                    ri = rand() % 500;
                    rii = rand() % 100;
                    double y = ri + (double)rii / 100.0;

                    cout << x  << " " << y << endl;

                    bool flag = false;
                    for (int i = 0; i < locs.size(); i++)
                    {
                        if (abs(locs[i].x - x) < eps && abs(locs[i].y - y) < eps) flag = true;
                    }

                    if (flag == false)
                    {
                        locs.push_back(point(x, y));
                        break;
                    }
                }
                generate_graph();
                dfs();
                art_points();

                cout << iscon << is2con << endl;

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
	gluLookAt(0,0,0,	0,0,-1,	0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects


    int i,j;


    for (i = 0; i < locs.size(); i++)
    {
//        if (isap[i] == false) glColor3f(1, 0, 0);
//        else glColor3f(0, 0, 1);
        glColor3f(0, 0, 1);
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

    glutPostRedisplay();
	glutTimerFunc(50, animate, 0);
}

void dsgen()
{
    fout = fopen("data16.txt","w");

    fprintf(fout,"%d\n", dsize);
	fprintf(fout,"%f\n", cr);

	int cnt = 0;
    point p;
    int ri, rii;

    while (cnt < 100)
    {
        locs.clear();
        while ((int)locs.size() < dsize)
        {
            ri = rand() % 500;
            rii = rand() % 100;
            double x = ri + (double)rii / 100.0;
            ri = rand() % 500;
            rii = rand() % 100;
            double y = ri + (double)rii / 100.0;

            //cout << x  << " " << y << endl;

            bool flag = false;
            for (int i = 0; i < locs.size(); i++)
            {
                if (abs(locs[i].x - x) < eps && abs(locs[i].y - y) < eps) flag = true;
            }

            if (flag == false)
            {
                locs.push_back(point(x, y));
            }

        }

        generate_graph();
        dfs();
        art_points();
        //cout << locs.size() << endl;

        //cout << iscon << is2con << endl;
        if (iscon == true && is2con == false)
        {
            add_to_file();
            cnt++;
            cout << "yes" << endl;
        }

    }

    fclose(fout);

}

void init(){
	//codes for initialization

	srand(time(NULL));

	dsgen();






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

