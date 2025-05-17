#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef long double ld;
typedef unsigned int ui;
typedef unsigned long long ull;

#define Rep(i,n) for(int i = 0; i < (n); ++i)
#define Repd(i,n) for(int i = (n)-1; i >= 0; --i)
#define For(i,a,b) for(int i = (a); i <= (b); ++i)
#define Ford(i,a,b) for(int i = (a); i >= (b); --i)
#define Fit(i,v) for(__typeof((v).begin()) i = (v).begin(); i != (v).end(); ++i)
#define Fitd(i,v) for(__typeof((v).rbegin()) i = (v).rbegin(); i != (v).rend(); ++i)
#define mp make_pair
#define pb push_back
#define fi first
#define se second
#define sz(a) ((int)(a).size())
#define all(a) (a).begin(), (a).end()
#define ms(a,x) memset(a, x, sizeof(a))

template<class F, class T> T convert(F a, int p = -1) { stringstream ss; if (p >= 0) ss << fixed << setprecision(p); ss << a; T r; ss >> r; return r; }
template<class T> T gcd(T a, T b){ T r; while (b != 0) { r = a % b; a = b; b = r; } return a;}
template<class T> T lcm(T a, T b) { return a / gcd(a, b) * b; }
template<class T> T sqr(T x) { return x * x; }
template<class T> T cube(T x) { return x * x * x; }
template<class T> int getbit(T s, int i) { return (s >> i) & 1; }
template<class T> T onbit(T s, int i) { return s | (T(1) << i); }
template<class T> T offbit(T s, int i) { return s & (~(T(1) << i)); }
template<class T> int cntbit(T s) { return s == 0 ? 0 : cntbit(s >> 1) + (s & 1); }

typedef pair<int, int> ii;
typedef vector<int> vi;
typedef vector<ii> vii;

const ld PI = acos(-1.0);
const ld eps = 1e-9;
const int dr[] = {-1, 0, +1, 0};
const int dc[] = {0, +1, 0, -1};
const int inf = (int)1e9 + 5;
const ll linf = (ll)1e16 + 5;
const ll mod = (ll)1e9 + 7;
const int UNVISITED = -1;
const int VISITED   =  1;
const int MaxN      = 10000;
const int INF       = 1000000000;

//int dr[] = {1,1,0,-1,-1,-1, 0, 1}; // trick to explore an implicit 2D grid
//int dc[] = {0,1,1, 1, 0,-1,-1,-1}; // S,SE,E,NE,N,NW,W,SW neighbors

int V, E, VL, VR;
vi dfs_num;
vector<vii> AdjList;

void enter()
{
    int u, v, c;
    cin >> V >> E;
    Rep(i, V) AdjList[i].clear();
    Rep(i, E) {
        cin >> u >> v >> c;
        AdjList[u-1].pb(mp(v-1, c));
        AdjList[v-1].pb(mp(u-1, c));
    }
}

void dfs(int u)
{
    dfs_num[u] = VISITED;
    for (int j = 0; j < (int)AdjList[u].size(); j++)
    {
        ii v = AdjList[u][j];
        if (dfs_num[v.first] == UNVISITED)
            dfs(v.first);
    }
}

void bfs(int s) {
    vi d(V, INF);
    d[s] = 0; // distance from source s to s is 0
    queue<int> q; q.push(s); // start from source

    while (!q.empty()) {
        int u = q.front(); q.pop(); // queue: layer by layer!
        for (int j = 0; j < (int)AdjList[u].size(); j++) {
            ii v = AdjList[u][j]; // for each neighbor of u
            if (d[v.first] == INF) { // if v.first is unvisited + reachable
                d[v.first] = d[u] + 1; // make d[v.first] != INF to flag it
                q.push(v.first); // enqueue v.first for the next iteration
            }
        }
    }
}

int connectedComponent() {
    int numCC = 0;
    dfs_num.assign(V, UNVISITED); // sets all vertices’ state to UNVISITED
    for (int i = 0; i < V; i++) // for each vertex i in [0..V-1]
        if (dfs_num[i] == UNVISITED) // if vertex i is not visited yet
        {
            ++numCC;
            dfs(i);
        }

    return numCC;
}

bool isBipartiteGraph(int s) {
    // inside int main()
    queue<int> q; q.push(s);
    vi color(V, INF); color[s] = 0;
    bool isBipartite = true; // addition of one boolean flag, initially true
    while (!q.empty() & isBipartite)  // similar to the original BFS routine
    {
        int u = q.front(); q.pop();
        for (int j = 0; j < (int)AdjList[u].size(); j++)
        {
            ii v = AdjList[u][j];
            if (color[v.first] == INF)  // but, instead of recording distance,
            {
                color[v.first] = 1 - color[u]; // we just record two colors {0, 1}
                q.push(v.first);
            }
            else if (color[v.first] == color[u])  // u & v.first has same color
            {
                isBipartite = false; break;
            }
        }
    } // we have a coloring conflict

    return isBipartite;
}

int dfsNumberCounter, dfsRoot, rootChildren;
vi dfs_low, dfs_parent, articulation_vertex;

void articulationPointAndBridge(int u)
{
    dfs_low[u] = dfs_num[u] = dfsNumberCounter++; // dfs_low[u] <= dfs_num[u]
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ii v = AdjList[u][j];
        if (dfs_num[v.first] == UNVISITED) { // a tree edge
            dfs_parent[v.first] = u;
            if (u == dfsRoot) rootChildren++; // special case if u is a root
            articulationPointAndBridge(v.first);
            if (dfs_low[v.first] >= dfs_num[u]) // for articulation point
                articulation_vertex[u] = true; // store this information first
            if (dfs_low[v.first] > dfs_num[u]) // for bridge
                printf(" Edge (%d, %d) is a bridge\n", u, v.first);
            dfs_low[u] = min(dfs_low[u], dfs_low[v.first]); // update dfs_low[u]
        }
        else if (v.first != dfs_parent[u]) // a back edge and not direct cycle
            dfs_low[u] = min(dfs_low[u], dfs_num[v.first]); // update dfs_low[u]
    }
}

int numSCC;
vi S, visited;

void findArticulationPointAndBridge()
{
    dfsNumberCounter = 0; dfs_num.assign(V, UNVISITED); dfs_low.assign(V, 0);
    dfs_parent.assign(V, 0); articulation_vertex.assign(V, 0);
    printf("Bridges:\n");
    for (int i = 0; i < V; i++)
        if (dfs_num[i] == UNVISITED) {
            dfsRoot = i; rootChildren = 0; articulationPointAndBridge(i);
            articulation_vertex[dfsRoot] = (rootChildren > 1);
        } // special case

    printf("Articulation Points:\n");
    for (int i = 0; i < V; i++)
        if (articulation_vertex[i])
            printf(" Vertex %d\n", i);
}

void tarjanSCC(int u)
{
    dfs_low[u] = dfs_num[u] = dfsNumberCounter++; // dfs_low[u] <= dfs_num[u]
    S.push_back(u); // stores u in a vector based on order of visitation
    visited[u] = 1;
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ii v = AdjList[u][j];
        if (dfs_num[v.first] == UNVISITED)
            tarjanSCC(v.first);
        if (visited[v.first]) // condition for update
            dfs_low[u] = min(dfs_low[u], dfs_low[v.first]);
    }

    if (dfs_low[u] == dfs_num[u]) { // if this is a root (start) of an SCC
        printf("SCC %d:", ++numSCC); // this part is done after recursion
        while (1) {
            int v = S.back(); S.pop_back(); visited[v] = 0;
            printf(" %d", v);
            if (u == v) break; }
        printf("\n");
    }
}

void findStrongConnectedComponents() //Directed Graph
{
    dfs_num.assign(V, UNVISITED); dfs_low.assign(V, 0); visited.assign(V, 0);
    dfsNumberCounter = numSCC = 0;
    for (int i = 0; i < V; i++)
        if (dfs_num[i] == UNVISITED)
            tarjanSCC(i);
}

class UnionFind {
private:
    vi p, rank;
public:
    UnionFind(int N) {
        rank.assign(N, 0);
        p.assign(N, 0); for (int i = 1; i < N; i++) p[i] = i;
    }

    int findSet(int i) { return (p[i] == i) ? i : (p[i] = findSet(p[i])); }
    bool isSameSet(int i, int j) { return findSet(i) == findSet(j); }
    void unionSet(int i, int j) {
        if (!isSameSet(i, j)) { // if from different set
            int x = findSet(i), y = findSet(j);
            if (rank[x] > rank[y]) p[y] = x; // rank keeps the tree short
            else { p[x] = y;
                if (rank[x] == rank[y]) rank[y]++; }
        }
    }
};

void krukalAlgorithm()
{
    vector< pair<int, ii> > EdgeList; // (weight, two vertices) of the edge
    int u, v, w;
    for (int i = 0; i < E; i++) {
        scanf("%d %d %d", &u, &v, &w); // read the triple: (u, v, w)
        EdgeList.push_back(make_pair(w, ii(u, v)));
    } // (w, u, v)

    sort(EdgeList.begin(), EdgeList.end()); // sort by edge weight O(E log E)
    // note: pair object has built-in comparison function
    int mst_cost = 0;
    UnionFind UF(V); // all V are disjoint sets initially
    for (int i = 0; i < E; i++) { // for each edge, O(E)
        pair<int, ii> front = EdgeList[i];
        if (!UF.isSameSet(front.second.first, front.second.second)) { // check
            mst_cost += front.first; // add the weight of e to MST
            UF.unionSet(front.second.first, front.second.second); // link them
        }
    } // note: the runtime cost of UFDS is very light

    // note: the number of disjoint sets must eventually be 1 for a valid MST
    printf("MST cost = %d (Kruskal’s)\n", mst_cost);
}

vi taken; // global boolean flag to avoid cycle
priority_queue<ii> pq; // priority queue to help choose shorter edges
// note: default setting for C++ STL priority_queue is a max heap

void process(int vtx) { // so, we use -ve sign to reverse the sort order
    taken[vtx] = 1;
    for (int j = 0; j < (int)AdjList[vtx].size(); j++) {
        ii v = AdjList[vtx][j];
        if (!taken[v.first]) pq.push(ii(-v.second, -v.first));
    }
} // sort by (inc) weight then by (inc) id

void primAlgorithm()
{
    taken.assign(V, 0); // no vertex is taken at the beginning
    process(0); // take vertex 0 and process all edges incident to vertex 0
    int mst_cost = 0, u, w;
    while (!pq.empty())
    { // repeat until V vertices (E=V-1 edges) are taken
        ii front = pq.top(); pq.pop();
        u = -front.second, w = -front.first; // negate the id and weight again
        if (!taken[u]) // we have not connected this vertex yet
            mst_cost += w, process(u); // take u, process all edges incident to u
    } // each edge is in pq only once!

    printf("MST cost = %d (Prim’s)\n", mst_cost);
}

// SSSP: Single-Source Shortest Paths
void dijkstraAlgorithm(int s) // SSSP on Weighted Graph
{
    vi dist(V, INF); dist[s] = 0; // INF = 1B to avoid overflow
    priority_queue< ii, vector<ii>, greater<ii> > pq; pq.push(ii(0, s));
    while (!pq.empty()) { // main loop
        ii front = pq.top(); pq.pop(); // greedy: get shortest unvisited vertex
        int d = front.first, u = front.second;
        if (d > dist[u]) continue; // this is a very important check
        for (int j = 0; j < (int)AdjList[u].size(); j++) {
            ii v = AdjList[u][j]; // all outgoing edges from u
            if (dist[u] + v.second < dist[v.first]) {
                dist[v.first] = dist[u] + v.second; // relax operation
                pq.push(ii(dist[v.first], v.first));
            }
        }
    } // this variant can cause duplicate items in the priority queue
}

vi dist(V, INF);

void bellmanAlgorithm(int s) // SSSP on Graph with Negative Weight Cycle
{
    dist[s] = 0;
    for (int i = 0; i < V - 1; i++) // relax all E edges V-1 times
        for (int u = 0; u < V; u++) // these two loops = O(E), overall O(VE)
            for (int j = 0; j < (int)AdjList[u].size(); j++) {
                ii v = AdjList[u][j]; // record SP spanning here if needed
                dist[v.first] = min(dist[v.first], dist[u] + v.second); // relax
            }
}

bool hasNegativeCycle() {
    // after running the O(VE) Bellman Ford’s algorithm shown above
    bool hasNegativeCycle = false;
    for (int u = 0; u < V; u++) // one more pass to check
        for (int j = 0; j < (int)AdjList[u].size(); j++) {
            ii v = AdjList[u][j];
            if (dist[v.first] > dist[u] + v.second) // if this is still possible
                hasNegativeCycle = true; // then negative cycle exists!
        }
    printf("Negative Cycle Exist? %s\n", hasNegativeCycle ? "Yes" : "No");
    return hasNegativeCycle;
}

list<int> cyc; // we need list for fast insertion in the middle
void EulerTour(list<int>::iterator i, int u)
{
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ii v = AdjList[u][j];
        if (v.second) { // if this edge can still be used/not removed
            v.second = 0; // make the weight of this edge to be 0 (‘removed’)
            for (int k = 0; k < (int)AdjList[v.first].size(); k++) {
                ii uu = AdjList[v.first][k]; // remove bi-directional edge
                if (uu.first == u && uu.second) {
                    uu.second = 0;
                    break;
                } }
            EulerTour(cyc.insert(i, u), v.first);
        }
    }
}

void findEulerTour(int A) {
    cyc.clear();
    EulerTour(cyc.begin(), A); // cyc contains an Euler tour starting at A
    for (list<int>::iterator it = cyc.begin(); it != cyc.end(); it++)
        printf("%d\n", *it); // the Euler tour
}

vi match, vis; // global variables
vector<vi> _AdjList;

int Aug(int l) { // return 1 if an augmenting path is found
    if (vis[l]) return 0; // return 0 otherwise
    vis[l] = 1;
    for (int j = 0; j < (int)_AdjList[l].size(); j++) {
        int r = _AdjList[l][j]; // edge weight not needed -> vector<vi> AdjList
        if (match[r] == -1 || Aug(match[r])) {
            match[r] = l; return 1; // found 1 matching
        } }
    return 0; // no matching
}

int findMaximumBipartite() {
    // build unweighted bipartite graph with directed edge left->right set
    int MCBM = 0;
    match.assign(V, -1); // V is the number of vertices in bipartite graph
    for (int l = 0; l < VL; l++) { // n = size of the left set
        vis.assign(VL, 0); // reset before each recursion
        MCBM += Aug(l);
    }
    printf("Found %d matchings\n", MCBM);
    return MCBM;
}

struct MincostFlow {
    int n, s, t, E, adj[maxe], next[maxe], last[maxv], which[maxv];
    ll totalCost, totalFlow, cap[maxe], flow[maxe], cost[maxe], pot[maxe], dist[maxv];

    void init(int _n, int _s, int _t) {
        n = _n; s = _s; t = _t;
        ms(last, -1); E = 0;
    }

    void add(int u, int v, ll ca, ll co) {
        adj[E] = v; cap[E] = ca; flow[E] = 0; cost[E] = +co; next[E] = last[u]; last[u] = E++;
        adj[E] = u; cap[E] =  0; flow[E] = 0; cost[E] = -co; next[E] = last[v]; last[v] = E++;
    }

    void bellman() {
        bool stop = false;
        ms(pot, 0);

        while (!stop) {
            stop = true;
            Rep(u, n) for (int e = last[u]; e != -1; e = next[e]) if (flow[e] < cap[e]) {
                int v = adj[e];
                if (pot[v] > pot[u] + cost[e]) {
                    pot[v] = pot[u] + cost[e];
                    stop = false;
                }
            }
        }
    }

    bool dijkstra() {
        typedef pair<ll, int> node;
        priority_queue<node, vector<node>, greater<node> > que;

        Rep(u, n) dist[u] = linf;
        dist[s] = 0;
        que.push(mp(0, s));

        while (!que.empty()) {
            ll dnow = que.top().fi;
            int u = que.top().se;
            que.pop();

            if (dnow > dist[u]) continue;

            for (int e = last[u]; e != -1; e = next[e]) if (flow[e] < cap[e]) {
                int v = adj[e];
                ll dnext = dnow + cost[e] + pot[u] - pot[v];

                if (dist[v] > dnext) {
                    dist[v] = dnext;
                    which[v] = e;
                    que.push(mp(dnext, v));
                }
            }
        }

        return dist[t] < linf;
    }

    bool maxflow(ll desireFlow = linf) {
        totalCost = totalFlow = 0;
        bellman();

        while (totalFlow < desireFlow) {
            if (!dijkstra()) return false;

            ll delta = desireFlow - totalFlow;
            for (int v = t, e = which[v]; v != s; v = adj[e ^ 1], e = which[v]) delta = min(delta, cap[e] - flow[e]);
            for (int v = t, e = which[v]; v != s; v = adj[e ^ 1], e = which[v]) {
                flow[e] += delta;
                flow[e ^ 1] -= delta;
            }

            totalFlow += delta;
            totalCost += delta * (dist[t] - pot[s] + pot[t]);

            Rep(u, n) pot[u] += dist[u];
        }

        return true;
    }

} mcf;

void process() {

}

void printResult() {

}


int main() {
    freopen ("input.txt" , "r", stdin);
    freopen ("output.txt" , "w", stdout);

    ios::sync_with_stdio(false); cin.tie(0);
    int T;
    cin >> T;
    while (T--) {

    }

    fclose(stdin); fclose(stdout);
    return 0;
}
