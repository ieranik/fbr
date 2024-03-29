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

using namespace std;

#define eps 1e-6
#define inf 1e12
#define num_test_cases 100
#define K 3

double com_range;
int num_points;

FILE* fin;

class point {
public:
    double x,y;
    point(double xx=0.0, double yy=0.0)
    {
        x=xx;
        y=yy;
    }
};

struct edge {
    int u, v;
    double w;
};

bool compare_edge(edge e1, edge e2) {
    return (e1.w < e2.w);
}

double sqr_dist(point a, point b) {
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}


vector< vector<int> > G;
vector<point> locs;
vector<point> allpoints;
vector<bool> isap;
bool iscon, is2con;

vector< pair<int, int> > bridgelist;
vector< pair<int, int> > add_edges;


// generate graph from positions (locs)
// edge exists if distance is leq communication range (com_range)
void generate_graph() {
    int i, j;
    G.clear();

    for (i = 0; i < locs.size(); i++) {
        vector<int> tmp;
        tmp.clear();
        for (j = 0; j < locs.size(); j++) {
            if (j != i && sqr_dist(locs[i], locs[j]) <= com_range * com_range + eps) tmp.push_back(j);
        }
        G.push_back(tmp);
    }
}


// helper function to determine articulation points
// articulation points are the nodes removal of which disconnect the graph
void get_ap(int u, vector<bool>& visited, vector<int>& disc,  vector<int>& low, vector<int>& parent, vector<bool>& ap) {
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
    for (int i = 0; i < G[u].size(); i++) {
        int v = G[u][i];  // v is current adjacent of u

        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (!visited[v]) {
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


// determine articulation points
// articulation points are the nodes removal of which disconnect the graph
void art_points() {
    int sz = (int) G.size();

    vector<bool> visited(sz);
    vector<bool> ap(sz);
    vector<int> disc(sz);
    vector<int> low(sz);
    vector<int> parent(sz);
  
    // Initialize parent and visited, and ap(articulation point) arrays
    for (int i = 0; i < G.size(); i++) {
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

    //set is2con true if there is no articulation points, and false otherwise
    isap = vector<bool>(sz, false);
    is2con = true;
    for (int i = 0; i < G.size(); i++) {
        isap[i] = ap[i];
        if (ap[i] == true) is2con = false;
    }
}


// dfs helper function
// returns the number of nodes visited
int dfsvis(int u, bool visited[]) {
    visited[u] = true;
    int cnt = 1;

    for (int i = 0; i < G[u].size(); i++) {
        int v = G[u][i];
        if (!visited[v]) cnt += dfsvis(v, visited);
    }
    return cnt;
}

// vanilla dfs
// determines if graph is connected or not
// sets the iscon flag accordingly
void dfs()
{
    bool *visited = new bool[G.size()];

    for (int i = 0; i < G.size(); i++) visited[i] = false;

    int cnt = dfsvis(0, visited);

    if (cnt < (int)G.size()) iscon = false;
    else iscon = true;
}

// modified dfs that does not visit nodes in v
// determines if graph G minus the nodes in v is connected or not
// sets the iscon flag accordingly
void dfs_mod(vector<int> v)
{
    bool *visited = new bool[G.size()];

    for (int i = 0; i < G.size(); i++) visited[i] = false;
    for (int i = 0; i < v.size(); i++) visited[v[i]] = true;

    int startnode = 0;
    while (true) {
        if (visited[startnode] == true) startnode++;
        else break;
    }

    int cnt = dfsvis(startnode, visited);

    if (cnt + K - 1 < (int)G.size()) iscon = false;
    else iscon = true;
}

// test is graph is kconneced
bool is_kconnected(){
    if (K == 2) {
        dfs();
        art_points();
        return iscon and is2con;
    }
    else if (K == 3) {
        for (int i = 0; i < num_points; i++) {
            for (int j = i + 1; j < num_points; j++) {
                vector<int> v = {i, j};
                dfs_mod(v);
                if (iscon == false) return false;
            }
        }
        return true;
    }
    else if (K == 4) {
        for (int i = 0; i < num_points; i++) {
            for (int j = i + 1; j < num_points; j++) {
                for (int k = j + 1; k < num_points; k++) {
                    vector<int> v = {i, j, k};
                    dfs_mod(v);
                    if (iscon == false) return false;
                }
            }
        }
        return true;
    }
    return true;
}

// moves node c towards node p such that their distance equals com_range
void repos(int p, int c) {
    double len = sqrt(sqr_dist(locs[p], locs[c]));
    double dx = locs[c].x - locs[p].x;
    double dy = locs[c].y - locs[p].y;
    dx *= (com_range / len);
    dy *= (com_range / len);
    locs[c].x = locs[p].x + dx;
    locs[c].y = locs[p].y + dy;
}


// moves node r to position rloc
// if the movement disconnects the graph, initiate cascade of movements 
// similar to BFS
void cascade(int r, point rloc) {
    locs[r] = rloc;

    bool *vis = new bool[G.size()];
    for (int i = 0; i < G.size(); i++) vis[i] = false;

    queue<int> Q;

    Q.push(r);
    vis[r] = true;

    while(Q.empty() == false) {
        int t = Q.front();
        Q.pop();
        for (int i = 0; i < G[t].size(); i++) {
            int v = G[t][i];
            if (vis[v] == true) continue;
            vis[v] = true;
            if (sqr_dist(locs[t], locs[v]) > com_range * com_range) repos(t, v);

            Q.push(v);
        }
    }
}

// return where to move node c 
// move c towards node p half the distance to reconnect the edge
point reposhalf(int p, int c) {
    double len = sqrt(sqr_dist(locs[p], locs[c]));
    point tp;
    double dif = len - com_range;
    dif = (dif / 2.0) + com_range;
    double dx = locs[c].x - locs[p].x;
    double dy = locs[c].y - locs[p].y;
    dx *= (dif / len);
    dy *= (dif / len);
    tp.x = locs[p].x + dx;
    tp.y = locs[p].y + dy;
    return tp;
}


// return where to move node c 
// move c towards node p such that edge between them reconnects
point reposhalf2(int p, int c)
{
    double len = sqrt(sqr_dist(locs[p], locs[c]));
    point tp;
    double dif = com_range;
    double dx = locs[c].x - locs[p].x;
    double dy = locs[c].y - locs[p].y;
    dx *= (dif / len);
    dy *= (dif / len);
    tp.x = locs[p].x + dx;
    tp.y = locs[p].y + dy;
    return tp;
}

// move nodes in pairs to reconnect edges between them
// returns sum and maximum of robot movements
pair<double, double> adjust_pos() {
    vector<point> start = locs;
    for (int i = 0; i < add_edges.size(); i++) {
        cascade(add_edges[i].first, reposhalf(add_edges[i].second, add_edges[i].first));
        cascade(add_edges[i].second, reposhalf2(add_edges[i].first, add_edges[i].second));
    }

    double maxdis = -1;
    double sum = 0;

    for (int i = 0; i < locs.size(); i++) {
        double dist = sqrt(sqr_dist(locs[i], start[i]));
        sum += dist;
        maxdis = max(dist, maxdis);
    }

    return {sum, maxdis};
}

// calculate number of edges in graph G
int count_edges() {
    int cnt = 0;
    for (int i = 0; i < num_points; i++) cnt += G[i].size();
    return cnt / 2;
}

// determine which edges to reconnect to restore k-connectivity
// put the edges in add_edges vector
void augment_edges() {
    add_edges.clear();

    //cout << count_edges() << endl;

    // determine all edges currently not in G and put those in com_edges
    vector<edge> com_edges;
    for (int i = 0; i < num_points; i++) {
        for (int j = i + 1; j < num_points; j++) {
            if (sqr_dist(locs[i], locs[j]) > com_range * com_range + eps) {
                edge tmp;
                tmp.u = i;
                tmp.v = j;
                tmp.w = sqr_dist(locs[i], locs[j]);
                com_edges.push_back(tmp);
            }
        }
    }

    // sort the edges in increasing order of length
    sort(com_edges.begin(), com_edges.end(), compare_edge);
    
    // keep adding edges until graph becomes k-connected
    int cnt = 0;
    while (true) {
        edge tmp = com_edges[cnt];
        int u = tmp.u;
        int v = tmp.v;
        G[u].push_back(v);
        G[v].push_back(u);
        if (is_kconnected() == true) break;
        cnt++; 
    }

    //cout << count_edges() << endl;

    add_edges.push_back({com_edges[cnt].u, com_edges[cnt].v});

    // remove unnecessary edges without compromising k-connectivity
    for (int i = cnt - 1; i >= 0; i--) {
        edge tmp = com_edges[i];
        int u = tmp.u;
        int v = tmp.v;
        
        vector<int> nu;
        vector<int> nv;
        for (int j = 0; j < G[u].size(); j++) if (G[u][j] != v) nu.push_back(G[u][j]);
        for (int j = 0; j < G[v].size(); j++) if (G[v][j] != u) nv.push_back(G[v][j]);
        
        G[u] = nu;
        G[v] = nv;

        if (is_kconnected() == true) continue;

        G[u].push_back(v);
        G[v].push_back(u);
        add_edges.push_back({u, v});
    }

    //cout << count_edges() << endl << endl;
 
    return;
}

// take input from file
void input_data() {
    fin = fopen("data32.txt","r");
	fscanf(fin, " %d", &num_points);
	fscanf(fin, " %lf", &com_range);

	allpoints.clear();
	for (int i = 0; i < num_points * num_test_cases; i++) {
        point tmp;
        fscanf(fin, " %lf %lf", &tmp.x, &tmp.y);
        allpoints.push_back(tmp);
    }

    fclose(fin);
}


int main(){

	input_data();

	srand(time(NULL));

	struct timeval tps;
    struct timeval tpe;
    long int st;
    long int et;

    gettimeofday(&tps, NULL);
    st = tps.tv_sec * 1000000 + tps.tv_usec ;
    
    double sum_sum = 0;
    double max_sum = 0;

    for (int i = 0; i < num_test_cases; i++)
    {
        // process each test case
        locs = vector<point>(allpoints.begin() + i * num_points, allpoints.begin() + (i + 1) * num_points);

        generate_graph();
        augment_edges();
        pair<double, double> ret = adjust_pos();
        sum_sum += ret.first;
        max_sum += ret.second;
        // generate_graph();
        // cout << is_kconnected() << endl;
    }

    gettimeofday(&tpe, NULL);
    et = tpe.tv_sec * 1000000 + tpe.tv_usec ;

    long int oldt = (et - st);
    cout << "Sum of Movements: " << sum_sum / (num_test_cases * com_range) << endl;
    cout << "Maximum Movement: " << max_sum / (num_test_cases * com_range) << endl;
    cout << "Time: " << oldt / num_test_cases << endl;

    return 0;
}
