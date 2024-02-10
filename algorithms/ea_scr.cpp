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
#include <sys/time.h>
#include <time.h>



#include <windows.h>
#include <glut.h>

using namespace std;

#define pi (2*acos(0.0))
#define eps 0.0000001
#define inf 10000.0
#define numtestcases 1

double cr;
int numpoints;

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
vector<point> allpoints;
vector<bool> isap;
int numcv;

vector<bcnode> bcl;
vector<int> blocknum;
vector<int> vtb;
vector< pair<int, int> > bridgelist;

vector< vector<int> > T;
vector< vector< pair<int, double> > > DT;
amnode am[100][100];
amnode amd[100][100];
vector< pair<int, int> > edges;
vector< pair<int, int> > add_edges;


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

    //bridgelist.clear();

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
    for (int i = 0; i < G.size(); i++) isap[i] = ap[i];
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

//    cout << "here" << endl;
//    for (int i = 0; i < bcl.size(); i++) printbcn(bcl[i]);

}


void create_tree()
{
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

//    for (int i = 0; i < T.size(); i++)
//    {
//        cout << i << ": ";
//        for (int j = 0; j < T[i].size(); j++)
//        {
//            cout << T[i][j] << " ";
//        }
//        cout << endl;
//    }
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

    numcv = 0;
    while (bcl[numcv].iscv == true) numcv++;

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

//    cout << "hereee" << endl;
//
//    for (int i = 0; i < DT.size(); i++)
//    {
//        cout << i << ": ";
//        for (int j = 0; j < DT[i].size(); j++)
//        {
//            cout << DT[i][j].first << " ";
//        }
//        cout << endl;
//    }
//
//    for (int i = 0; i < bcl.size(); i++) cout << i << " " << bcl[i].depth << endl;

}



void super_impose()
{
    int n = bcl.size();
    for (int i = 0; i < n; i++) bcl[i].nodelist = removedup1(bcl[i].nodelist);

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            for (int k = 0; k < bcl[i].nodelist.size(); k++)
            {
                for (int l = 0; l < bcl[j].nodelist.size(); l++)
                {
                    int u = bcl[i].nodelist[k];
                    int v = bcl[j].nodelist[l];
                    double dis = sqrt(sqr_dist(locs[u], locs[v]));
                    if (dis < am[i][j].weight)
                    {
                        am[i][j].weight = dis;
                        am[i][j].src = u;
                        am[i][j].dst = v;
                    }
                }
            }
        }
    }


//    for (int i = 0; i < T.size(); i++)
//    {
//        for (int j = 0; j < T.size(); j++)
//        {
//            cout << am[i][j].weight << "\t";
//        }
//        cout << endl;
//    }
//    cout << endl;


}

//v is desc of u
bool is_desc(int u, int v)
{
    if (bcl[u].depth > bcl[v].depth) return false;

    int vd = bcl[v].depth;
    int ud = bcl[u].depth;

    int dif = vd - ud;

    while (dif--) v = DT[v][0].first;

    return (u == v);
}

int onedesc(int u, int v)
{
    int vd = bcl[v].depth;
    int ud = bcl[u].depth;

    int dif = vd - ud - 1;

    while (dif--) v = DT[v][0].first;

    return v;
}


int lca(int u, int v)
{
    int vd = bcl[v].depth;
    int ud = bcl[u].depth;

    int dif = vd - ud;

    if (dif > 0)
    {
        while (dif--) v = DT[v][0].first;
    }
    else
    {
        dif = -dif;
        while (dif--) u = DT[u][0].first;
    }
    while (u != v)
    {
        u = DT[u][0].first;
        v = DT[v][0].first;
    }
    return u;
}




void create_aug_tree()
{
    int n = bcl.size();

    cout << "n: " << n << endl;

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (am[i][j].weight == inf) continue;
            int id = bcl[i].depth;
            int jd = bcl[j].depth;


            if (id > jd && is_desc(j, i) == true)
            {
                if (amd[j][i].weight > am[i][j].weight)
                {
                    amd[j][i].src = j;
                    amd[j][i].dst = i;
                    amd[j][i].weight = am[i][j].weight;
                }
                continue;
            }

            if (id < jd && is_desc(i, j) == true)
            {
                if (amd[i][j].weight > am[i][j].weight)
                {
                    amd[i][j].src = j;
                    amd[i][j].dst = i;
                    amd[i][j].weight = am[i][j].weight;
                }
                continue;
            }

            int t = lca(i, j);

            if (amd[t][i].weight > am[i][j].weight)
            {
                amd[t][i].src = j;
                amd[t][i].dst = i;
                amd[t][i].weight = am[i][j].weight;
            }

            if (amd[t][j].weight > am[i][j].weight)
            {
                amd[t][j].src = j;
                amd[t][j].dst = i;
                amd[t][j].weight = am[i][j].weight;
            }

            if (amd[j][i].weight > am[i][j].weight)
            {
                amd[j][i].src = j;
                amd[j][i].dst = i;
                amd[j][i].weight = am[i][j].weight;
            }

            if (amd[i][j].weight > am[i][j].weight)
            {
                amd[i][j].src = j;
                amd[i][j].dst = i;
                amd[i][j].weight = am[i][j].weight;
            }
        }
    }
//    cout << "amd" << endl;
//    for (int i = 0; i < T.size(); i++)
//    {
//        for (int j = 0; j < T.size(); j++)
//        {
//            cout << amd[i][j].weight << "\t";
//        }
//        cout << endl;
//    }
//
//    cout << endl;
//
//    for (int i = 0; i < T.size(); i++)
//    {
//        for (int j = 0; j < T.size(); j++)
//        {
//            cout << amd[i][j].src << "\t";
//        }
//        cout << endl;
//    }
//
//    cout << endl;
//
//
//    for (int i = 0; i < T.size(); i++)
//    {
//        for (int j = 0; j < T.size(); j++)
//        {
//            cout << amd[i][j].dst << "\t";
//        }
//        cout << endl;
//    }
//
//    cout << "amd over" << endl;






    for (int i = 0; i < T.size(); i++)
    {
        if (bcl[i].iscv == true)
        {
            for (int j = 0; j < T.size(); j++)
            {
                if (j != i && is_desc(i, j) == true ) //&& amd[i][j].weight < inf
                {
                    double pw = amd[i][j].weight;
                    amd[i][j].weight = inf;

                    int od = onedesc(i, j);
                    if (od != j && amd[od][j].weight > pw)
                    {
                        amd[od][j].weight = pw;
                        amd[od][j].src = i;
                        amd[od][j].dst = j;
                    }
                }
            }
        }
    }

//    cout << endl;
//
//    for (int i = 0; i < T.size(); i++)
//    {
//        for (int j = 0; j < T.size(); j++)
//        {
//            cout << amd[i][j].weight << "\t";
//        }
//        cout << endl;
//    }
//
//    cout << endl;
//
//    for (int i = 0; i < T.size(); i++)
//    {
//        for (int j = 0; j < T.size(); j++)
//        {
//            cout << amd[i][j].src << "\t";
//        }
//        cout << endl;
//    }
//
//    cout << endl;
//
//    for (int i = 0; i < T.size(); i++)
//    {
//        for (int j = 0; j < T.size(); j++)
//        {
//            cout << amd[i][j].dst << "\t";
//        }
//        cout << endl;
//    }
//
//    cout << "amd2 over" << endl;

    for (int i = 0; i < DT.size(); i++)
    {
        if (DT[i].size() > 0)
        {
            amd[i][DT[i][0].first].weight = 0;
            amd[i][DT[i][0].first].src = i;
            amd[i][DT[i][0].first].dst = DT[i][0].first;
        }

    }

//    for (int i = 0; i < DT.size(); i++)
//    {
//        amd[i][11].weight = inf;
//        amd[11][i].weight = inf;
//
//    }

//    cout << endl;
//
//    for (int i = 0; i < T.size(); i++)
//    {
//        for (int j = 0; j < T.size(); j++)
//        {
//            cout << amd[i][j].weight << "\t";
//        }
//        cout << endl;
//    }
//
//    cout << endl;
//
//    for (int i = 0; i < T.size(); i++)
//    {
//        for (int j = 0; j < T.size(); j++)
//        {
//            cout << amd[i][j].src << "\t";
//        }
//        cout << endl;
//    }
//
//    cout << endl;
//
//    for (int i = 0; i < T.size(); i++)
//    {
//        for (int j = 0; j < T.size(); j++)
//        {
//            cout << amd[i][j].dst << "\t";
//        }
//        cout << endl;
//    }

}




int findmin(vector<int> v, double cost[])
{
    double mn = inf;
    int ret;
    for (int i = 0; i < v.size(); i++)
    {
        if (cost[v[i]] < mn)
        {
            ret = v[i];
            mn = cost[v[i]];
        }
    }
    return ret;
}

vector<int> delete_elem(vector<int> v, int e)
{
    vector<int> ret;
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] != e)
        {
            ret.push_back(v[i]);
        }
    }
    return ret;
}

bool is_elem(vector<int> F, int e)
{
    for (int i = 0; i < F.size(); i++) if (F[i] == e) return true;
    return false;
}


void aug_edges()
{
    bool *visited = new bool[T.size()];
    double *cost = new double[T.size()];
    int *parent = new int[T.size()];
    for (int i = 0; i < T.size(); i++)
    {
        cost[i] = inf;
        visited[i] = false;
        parent[i] = -1;
    }
    vector<int> F;


    cost[numcv] = 0;
    F.push_back(numcv);

    int loop = T.size();

    while (loop--)
    {

//        for (int i = 0; i < T.size(); i++)
//            cout << i << "\t" << visited[i] << "\t" << cost[i] << "\t" << parent[i] << "\n";
//        cout << endl;
        int v = findmin(F, cost);
        visited[v] = true;
        F = delete_elem(F,v);
        for (int w = 0; w < T.size(); w++)
        {
            if (amd[v][w].weight < inf && visited[w] == false && is_elem(F, w) == false)
            {
                F.push_back(w);
                cost[w] = amd[v][w].weight;
                parent[w] = v;
            }
            else if (is_elem(F, w) == true && amd[v][w].weight < cost[w])
            {
                cost[w] = amd[v][w].weight;
                parent[w] = v;
            }
        }
    }

//    for (int i = 0; i < T.size(); i++)
//        cout << i << "\t" << visited[i] << "\t" << cost[i] << "\t" << parent[i] << "\n";
//    cout << endl;
    edges.clear();

    for (int i = 0; i < T.size(); i++)
    {
        if (i != numcv) edges.push_back(pair<int, int> (amd[parent[i]][i].src, amd[parent[i]][i].dst));
    }

    //for (int i = 0; i < edges.size(); i++) cout << edges[i].first << " " << edges[i].second << endl;
}


bool is_adj(int u, int v)
{
    if (u == -1 || v == -1) return true;
    for (int i = 0; i < T[u].size(); i++) if (T[u][i] == v) return true;
    return false;
}

pair<int, int> minp(int u, int v)
{
    double mindis = inf;
    pair<int, int> ret = pair<int, int>(-1, -1);
    for (int i = 0; i < bcl[u].nodelist.size(); i++)
    {
        for (int j = 0; j < bcl[v].nodelist.size(); j++)
        {
            int uu = bcl[u].nodelist[i];
            int vv = bcl[v].nodelist[j];
            double dis = sqrt(sqr_dist(locs[uu], locs[vv]));
            if (dis < mindis)
            {
                mindis = dis;
                ret.first = uu;
                ret.second = vv;
            }
        }
        //cout << mindis << endl;
    }
    return ret;
}


void final_edges()
{
    add_edges.clear();
    for (int i = 0; i < edges.size(); i++)
    {
        //cout << edges[i].first << " " << edges[i].second << endl;

        if (is_adj(edges[i].first, edges[i].second) == false)
        {
            pair<int, int> piit = minp(edges[i].first, edges[i].second);
            //cout << piit.first << " " << piit.second << endl;
            if (piit.first != -1) add_edges.push_back(piit);
            //cout << piit.first << " " << piit.second << endl;

        }
    }
}

void generate_file()
{
	fprintf(fout,"%d\n", locs.size());
	for (int i = 0; i < locs.size(); i++)
    {
        fprintf(fout,"%f %f\n", locs[i].x, locs[i].y);
    }

    int edgecount = 0;
    for (int i = 0; i < locs.size(); i++)
    {
        edgecount += G[i].size();
    }
    edgecount /= 2;
    edgecount += add_edges.size();

    fprintf(fout,"%d\n", edgecount);


    for (int i = 0; i < locs.size(); i++)
    {
        for (int j = 0; j < G[i].size(); j++)
        {
            if (i < G[i][j]) fprintf(fout,"%d %d\n", i, G[i][j]);
        }
    }


    for (int i = 0; i < add_edges.size(); i++)
    {
        fprintf(fout,"%d %d\n", add_edges[i].first, add_edges[i].second);
    }




	//fclose(fout);



    //fin=fopen("input.txt","r");
	//fscanf(fin,"%d",&i);

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

point reposhalf(int p, int c)
{
    double len = sqrt(sqr_dist(locs[p], locs[c]));
    point tp;
    double dif = len - cr;
    dif = (dif / 2.0) + cr;
    double dx = locs[c].x - locs[p].x;
    double dy = locs[c].y - locs[p].y;
    dx *= (dif / len);
    dy *= (dif / len);
    tp.x = locs[p].x + dx;
    tp.y = locs[p].y + dy;
    return tp;
}

point reposhalf2(int p, int c)
{
    double len = sqrt(sqr_dist(locs[p], locs[c]));
    point tp;
    double dif = cr;
    double dx = locs[c].x - locs[p].x;
    double dy = locs[c].y - locs[p].y;
    dx *= (dif / len);
    dy *= (dif / len);
    tp.x = locs[p].x + dx;
    tp.y = locs[p].y + dy;
    return tp;
}

double adjust_pos()
{
    vector<point> start = locs;
    for (int i = 0; i < add_edges.size(); i++)
    {
        cascade(add_edges[i].first, reposhalf(add_edges[i].second, add_edges[i].first));
        cascade(add_edges[i].second, reposhalf2(add_edges[i].first, add_edges[i].second));
    }
    double maxdis = -1;
    for (int i = 0; i < locs.size(); i++)
    {
        if (sqr_dist(locs[i], start[i]) > maxdis)
            maxdis = sqr_dist(locs[i], start[i]);
    }

    fprintf(fout, "%f\n", sqrt(maxdis));
    return sqrt(maxdis);
}









void input_data()
{
    fin = fopen("errortest2.txt","r");
	fscanf(fin, " %d", &numpoints);
	fscanf(fin, " %lf", &cr);

	allpoints.clear();
	for (int i = 0; i < numpoints * numtestcases; i++)
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


    for (i = 0; i < locs.size(); i++)
    {
        if (isap[i] == false) glColor3f(1, 0, 0);
        else glColor3f(0, 0, 1);
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

    glColor3f(0, 1, 1);
    for (i = 0; i < add_edges.size(); i++)
    {
        glBegin(GL_LINES);
        {
            glVertex3f(locs[add_edges[i].first].x, locs[add_edges[i].first].y, 0);
            glVertex3f(locs[add_edges[i].second].x, locs[add_edges[i].second].y, 0);
        }
        glEnd();
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

	input_data();

	fout = fopen("out256g.txt","w");
	fprintf(fout, "%f\n", cr);

	srand(time(NULL));

	struct timeval tps;
    struct timeval tpe;
    long int st;
    long int et;

    gettimeofday(&tps, NULL);
    st = tps.tv_sec * 1000000 + tps.tv_usec ;
    double sum = 0;

    for (int i = 0; i < numtestcases; i++)
    {
        locs = vector<point>(allpoints.begin() + i * numpoints, allpoints.begin() + (i + 1) * numpoints);


        generate_graph();
        art_points();
        create_bcl();
        create_tree();
        create_dtree();
        super_impose();
        create_aug_tree();
        aug_edges();
        final_edges();
        //generate_file();
        sum += adjust_pos();

    }

    gettimeofday(&tpe, NULL);
    et = tpe.tv_sec * 1000000 + tpe.tv_usec ;

    long int oldt = (et - st);
    cout << "Distance: " << sum / (100 * cr) << endl;
    cout << "Time: " << oldt / 100 << endl;

    fclose(fout);















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

