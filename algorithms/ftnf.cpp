#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <time.h>

using namespace std;

#define INF 1e18
#define eps 1e-6
#define num_test_cases 100
#define K 4

 
struct point {
    double x, y;
};
 
struct circle {
    point c;
    double r;
};

struct point_dis {
    int idx;
    double dis;
};

bool compare_point_dis(point_dis e1, point_dis e2) {
    return (e1.dis < e2.dis);
}

 
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
    // fin = fopen("data16.txt","r");
	// fscanf(fin, " %d", &num_points);
	// fscanf(fin, " %lf", &com_range);

    //cout << num_points << " " << com_range << endl;

    locs.clear();
	for (int i = 0; i < num_points; i++)
    {
        point tmp;
        fscanf(fin, " %lf %lf", &tmp.x, &tmp.y);
        //cout << tmp.x << " " << tmp.y << endl;
        locs.push_back(tmp);
    }

    //fclose(fin);
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


vector<int> k_nearest_neibors(int ci) {
    vector<point_dis> v;
    for (int i = 0; i < locs.size(); i++) {
        if (i == ci) continue;
        double dis = dist(locs[i], locs[ci]);
        point_dis tmp;
        tmp.dis = dis;
        tmp.idx = i;
        v.push_back(tmp);
    }
    sort(v.begin(), v.end(), compare_point_dis);
    vector<int> ret;
    for (int i = 0; i < K; i++) ret.push_back(v[i].idx);
    return ret;
}


double farthest_pair_dist(int ci, vector<int> v) {
    v.push_back(ci);
    double d = -INF;
    for (int i = 0; i < v.size(); i++) {
        for (int j = i + 1; j < v.size(); j++) {
            d = max(d, dist(locs[v[i]], locs[v[j]]));
        }
    }
    return d;
}


point find_pos(int oi, vector<int> inc) {
    point p = locs[oi];
    double max_dis = -INF;
    int id = -1;
    for (int i = 0; i < inc.size(); i++) {
        double dis = dist(p, new_locs[inc[i]]);
        if (dis > max_dis) {
            id = inc[i];
            max_dis = dis;
        }
    }

    point ret;
    point b = new_locs[id];
    ret.x = b.x + (p.x - b.x) * (com_range / max_dis);
    ret.y = b.y + (p.y - b.y) * (com_range / max_dis);
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
    fin = fopen("data16.txt","r");
	fscanf(fin, " %d", &num_points);
	fscanf(fin, " %lf", &com_range);

    srand(time(NULL));

	struct timeval tps;
    struct timeval tpe;
    long int st;
    long int et;

    gettimeofday(&tps, NULL);
    st = tps.tv_sec * 1000000 + tps.tv_usec ;


    int iter = num_test_cases;

    double sum_sum = 0;
    double max_sum = 0;

    while (iter--) {
        input_data();
        generate_graph();

        if (is_kconnected()) continue;

        new_locs = vector<point>(num_points);
        is_processed = vector<bool>(num_points, false);

        circle mec = welzl(locs);
        point c = mec.c;

        //cout << c.x << " " << c.y << endl;

        int ci = closest_to_center(c);
        is_processed[ci] = true;
        new_locs[ci] = locs[ci];

        vector<int> tmp = k_nearest_neibors(ci);
        //cout << ci << " " << tmp.first << " " << tmp.second << endl;

        double fd = farthest_pair_dist(ci, tmp);
        double factor = min(1.0, com_range / fd);

        //cout << factor << endl;
        for (int i = 0; i < tmp.size(); i++) {
            is_processed[tmp[i]] = true;
            new_locs[tmp[i]].x = locs[ci].x + (locs[tmp[i]].x - locs[ci].x) * factor; 
            new_locs[tmp[i]].y = locs[ci].y + (locs[tmp[i]].y - locs[ci].y) * factor; 
        }

        // cout << new_locs[tmp.first].x << " " << new_locs[tmp.first].y << endl;
        // cout << new_locs[tmp.second].x << " " << new_locs[tmp.second].y << endl;

        int loop = num_points - K - 1;


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

            vector<point_dis> v;

            for (int i = 0; i < in_idx.size(); i++) {
                if (in_idx[i] == in_id) continue;
                double dis = dist(new_locs[in_id], new_locs[in_idx[i]]);
                point_dis tmp;
                tmp.dis = dis;
                tmp.idx = in_idx[i];
                v.push_back(tmp);
            }
            
            sort(v.begin(), v.end(), compare_point_dis);

            vector<int> inc;
            inc.push_back(in_id);
            for (int i = 0; i < K -1; i++) inc.push_back(v[i].idx);


            point new_pos = find_pos(out_id, inc);

            // cout << out_id << endl;
            // print_point(locs[out_id]);
            // print_point(new_pos);
            // cout << endl;

            is_processed[out_id] = true;
            new_locs[out_id] = new_pos;

        }

        //for (int i = 0; i < num_points; i++) print_point(new_locs[i]);
        double max_dis = -INF;
        double sum = 0;
        for (int i = 0; i < num_points; i++) {
            double distance = dist(locs[i], new_locs[i]);
            max_dis = max(max_dis, distance);
            sum += distance;
        }
        //cout << max_dis << endl;
        max_sum += max_dis;
        sum_sum += sum;
        
        
        // locs = new_locs;
        // generate_graph();
        // if (is_kconnected()) cout << "ok" << endl;
        // else cout << "not ok" << endl;


    }

    gettimeofday(&tpe, NULL);
    et = tpe.tv_sec * 1000000 + tpe.tv_usec ;

    long int oldt = (et - st);
    cout << "Sum of Movements: " << sum_sum / (num_test_cases * com_range) << endl;
    cout << "Maximum Movement: " << max_sum / (num_test_cases * com_range) << endl;
    cout << "Time: " << oldt / num_test_cases << endl;
    
    
    // cout << wcnt << endl;
    // cout << sum / (100 * com_range) << endl;

    fclose(fin);
    
    
 
    return 0;
}