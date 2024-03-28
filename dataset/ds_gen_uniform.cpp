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

using namespace std;

#define pi (2*acos(0.0))
#define eps 0.0000001
#define inf 10000.0
#define cr 150.0
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


vector< vector<int> > G;
vector<point> locs;
vector<bool> isap;
vector<int> leaves;
int numcv;
bool is2con;
bool iscon;


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


void add_to_file()
{
	for (int i = 0; i < locs.size(); i++)
    {
        fprintf(fout,"%f %f\n", locs[i].x, locs[i].y);
    }
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
      
        if (iscon == true && is2con == false)
        {
            add_to_file();
            cnt++;
        }

    }

    fclose(fout);

}

int main(){
	srand(time(NULL));
	dsgen();
    return 0;
}


