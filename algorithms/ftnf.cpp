#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

#define INF 1e18
#define eps 1e-12
 
struct point {
    double x, y;
};
 
struct circle {
    point c;
    double r;
};
 
double dist(const point& a, const point& b)
{
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}
 
bool is_inside(const circle& c, const point& p)
{
    return dist(c.c, p) <= c.r;
}
 
 
// Helper method to get a circle defined by 3 points, point A is origin
point get_circle_center(double bx, double by, double cx, double cy)
{
    double B = bx * bx + by * by;
    double C = cx * cx + cy * cy;
    double D = bx * cy - by * cx;
    return { (cy * B - by * C) / (2 * D), (bx * C - cx * B) / (2 * D) };
}
 
// Function to return a unique circle that intersects three points
circle circle_from(const point& A, const point& B, const point& C)
{
    point I = get_circle_center(B.x - A.x, B.y - A.y, C.x - A.x, C.y - A.y);
 
    I.x += A.x;
    I.y += A.y;
    return {I, dist(I, A)};
}
 
// Function to return the smallest circle that intersects 2 points
circle circle_from(const point& A, const point& B)
{
    point C = { (A.x + B.x) / 2.0, (A.y + B.y) / 2.0 };
     return { C, dist(A, B) / 2.0 };
}
 
// Function to check whether a circle encloses the given points
bool is_valid_circle(const circle& c, const vector<point>& P)
{
    for (const point& p : P)
        if (!is_inside(c, p))
            return false;
    return true;
}
 
// Function to return the minimum enclosing circle for N <= 3
circle min_circle_trivial(vector<point>& P)
{
    if (P.empty()) {
        return { { 0, 0 }, 0 };
    }
    else if (P.size() == 1) {
        return { P[0], 0 };
    }
    else if (P.size() == 2) {
        return circle_from(P[0], P[1]);
    }
 
    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
 
            circle c = circle_from(P[i], P[j]);
            if (is_valid_circle(c, P))
                return c;
        }
    }
    return circle_from(P[0], P[1], P[2]);
}
 
// Welzl's algorithm
circle welzl_helper(vector<point>& P, vector<point> R, int n)
{
    // Base case when all points processed or |R| = 3
    if (n == 0 || R.size() == 3) {
        return min_circle_trivial(R);
    }
 
    // Pick a random point randomly
    int idx = rand() % n;
    point p = P[idx];
 
    // Put the picked point at the end of P
    swap(P[idx], P[n - 1]);
 
    // Get the MEC circle d from the set of points P - {p}
    circle d = welzl_helper(P, R, n - 1);
 
    // If d contains p, return d
    if (is_inside(d, p)) {
        return d;
    }
 
    // Otherwise, must be on the boundary of the MEC
    R.push_back(p);
 
    // Return the MEC for P - {p} and R U {p}
    return welzl_helper(P, R, n - 1);
}
 
circle welzl(const vector<point>& P)
{
    vector<point> P_copy = P;
    random_shuffle(P_copy.begin(), P_copy.end());
    return welzl_helper(P_copy, {}, P_copy.size());
}

FILE* fin;
FILE* fout;

int num_points;
double com_range;
vector<point> locs;
vector<point> new_locs;
vector<bool> isap;
vector< vector<int> > G;
vector< pair<int, int> > bridgelist;
vector<bool> is_processed;

bool iscon;
bool is2con;


void input_data()
{
    fin = fopen("data.txt","r");
	fscanf(fin, " %d", &num_points);
	fscanf(fin, " %lf", &com_range);

    //cout << num_points << " " << com_range << endl;

    locs.clear();
	for (int i = 0; i < num_points; i++)
    {
        point tmp;
        fscanf(fin, " %lf %lf", &tmp.x, &tmp.y);
        //cout << tmp.x << " " << tmp.y << endl;
        locs.push_back(tmp);
    }

    fclose(fin);
}

int closest_to_center_idx(point c, vector<point> p) {
    return 0;
    
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
            if (j != i && sqr_dist(locs[i], locs[j]) <= com_range * com_range + eps) tmp.push_back(j);
        }
        G.push_back(tmp);
    }
}



void get_ap(int u, vector<bool>& visited, vector<int>& disc,  vector<int>& low, vector<int>& parent, vector<bool>& ap)
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
    int sz = (int) G.size();

    vector<bool> visited(sz);
    vector<bool> ap(sz);
    vector<int> disc(sz);
    vector<int> low(sz);
    vector<int> parent(sz);
  

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


    isap = vector<bool>(sz, false);
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


bool is_biconnected(){
    dfs();
    art_points();
    return iscon and is2con;
}

int closest_to_center(point c) {
    int ret_idx = -1;
    double min_dis = INF;
    for (int i = 0; i < locs.size(); i++) {
        double dis = dist(c, locs[i]);
        if (dis < min_dis) {
            min_dis = dis;
            ret_idx = i;
        } 
    }
    return ret_idx;
}

pair<int, int> two_nearest_neibors(int ci) {
    int idx1 = -1;
    int idx2 = -1;
    double min_dis = INF;
    double min_dis2 = INF;
    for (int i = 0; i < locs.size(); i++) {
        if (i == ci) continue;
        double dis = dist(locs[i], locs[ci]);
        if (dis < min_dis) {
            min_dis2 = min_dis;
            min_dis = dis;
            idx2 = idx1;
            idx1 = i;
        }
        else if (dis < min_dis2) {
            min_dis2 = dis;
            idx2 = i;
        }
    }
    return pair<int, int>(idx1, idx2);
}

double farthest_pair_dist(int i1, int i2, int i3) {
    double d = -INF;
    d = max(d, dist(locs[i1], locs[i2]));
    d = max(d, dist(locs[i1], locs[i3]));
    d = max(d, dist(locs[i3], locs[i2]));
    return d;
}

point find_pos(point a, point b, point c) {
    double disb = dist(b, a);
    double disc = dist(c, a);
    if (disb < disc) b = c;
    disb = dist(b, a);
    point ret;
    ret.x = b.x + (a.x - b.x) * (com_range / disb);
    ret.y = b.y + (a.y - b.y) * (com_range / disb);
    return ret;
}

void print_point(point p) {
    cout << p.x << " " << p.y << endl;
}

void print_vec(vector<int> v){ 
    for (int i = 0; i < v.size(); i++) {
        cout << v[i] << " ";
    }
    cout << endl;
}


int main() {
    input_data();
    generate_graph();

    if (is_biconnected()) return 0;

    new_locs = vector<point>(num_points);
    is_processed = vector<bool>(num_points, false);

    circle mec = welzl(locs);
    point c = mec.c;

    //cout << c.x << " " << c.y << endl;

    int ci = closest_to_center(c);
    is_processed[ci] = true;
    new_locs[ci] = locs[ci];

    pair<int, int> tmp = two_nearest_neibors(ci);
    //cout << ci << " " << tmp.first << " " << tmp.second << endl;

    double fd = farthest_pair_dist(ci, tmp.first, tmp.second);

    double factor = min(1.0, com_range / fd);

    //cout << factor << endl;
    is_processed[tmp.first] = true;
    is_processed[tmp.second] = true;

    // cout << ci << " " << tmp.first << " " << tmp.second << endl;

    new_locs[tmp.first].x = locs[ci].x + (locs[tmp.first].x - locs[ci].x) * factor; 
    new_locs[tmp.first].y = locs[ci].y + (locs[tmp.first].y - locs[ci].y) * factor; 

    new_locs[tmp.second].x = locs[ci].x + (locs[tmp.second].x - locs[ci].x) * factor; 
    new_locs[tmp.second].y = locs[ci].y + (locs[tmp.second].y - locs[ci].y) * factor; 

    // cout << new_locs[tmp.first].x << " " << new_locs[tmp.first].y << endl;
    // cout << new_locs[tmp.second].x << " " << new_locs[tmp.second].y << endl;

    int loop = num_points - 3;


    while (loop--) {
        vector<int> in_idx;
        vector<int> out_idx;
        for (int i = 0; i < num_points; i++) {
            if (is_processed[i] == true) in_idx.push_back(i);
            else out_idx.push_back(i);
        }

        // print_vec(in_idx);
        // print_vec(out_idx);

        double min_dis = INF;
        int in_id = -1;
        int out_id = -1;
        for (int i = 0; i < in_idx.size(); i++) {
            for (int j = 0; j < out_idx.size(); j++) {
                double dis = dist(new_locs[in_idx[i]], locs[out_idx[j]]);
                if (dis < min_dis) {
                    min_dis = dis;
                    in_id = in_idx[i]; 
                    out_id = out_idx[j];
                }
            }
        }
        min_dis = INF;
        int s_id = -1;
        for (int i = 0; i < in_idx.size(); i++) {
            if (in_idx[i] == in_id) continue;
            double dis = dist(new_locs[in_id], new_locs[in_idx[i]]);
            if (dis < min_dis) {
                min_dis = dis;
                s_id = in_idx[i];
            }
        }
        point new_pos = find_pos(locs[out_id], new_locs[in_id], new_locs[s_id]);

        // cout << out_id << endl;
        // print_point(locs[out_id]);
        // print_point(new_pos);
        // cout << endl;

        is_processed[out_id] = true;
        new_locs[out_id] = new_pos;

    }

    for (int i = 0; i < num_points; i++) print_point(new_locs[i]);
    
    
    locs = new_locs;
    generate_graph();
    if (is_biconnected()) cout << "it works" << endl;
    
 
    return 0;
}