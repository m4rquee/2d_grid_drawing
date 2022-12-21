#include "gurobi_c++.h"
#include <array>
#include <cmath>
#include <sstream>

using namespace std;
const double EPS = 1E-1;
const long unsigned SEED = 42;// the random number generator seed

inline bool side(double p[2], const double u[2], double v[2]) {
    // Calculates the vector from u to v:
    v[0] -= u[0];
    v[1] -= u[1];
    // Calculates the vector from u to p:
    p[0] -= u[0];
    p[1] -= u[1];
    auto cross_prod = v[0] * p[1] - v[1] * p[0];
    return cross_prod >= 0;
}

// Intersection elimination callback. Whenever a feasible solution is found,
// find the edges intersections, and add an intersection elimination constraint
// if some edges intersect.
class IntersectElim : public GRBCallback {
public:
    int m;
    array<int, 2> *E;
    GRBVar *X, *Y, *edge_len;

    IntersectElim(int _m, array<int, 2> *_E, GRBVar *_X, GRBVar *_Y, GRBVar *_edge_len)
        : m(_m), E(_E), X(_X), Y(_Y), edge_len(_edge_len) {}

protected:
    void callback() override {
        try {
            // Check if all the model variables are integer:
            if (where != GRB_CB_MIPSOL) return;// as this code do not take advantage of the other options

            // Found an integer feasible solution - does it have an edge intersection?

            return;
            bool itersect;
            for (int e = 0; e < m; e++)
                for (int f = e + 1; f < m; f++) {
                    double a[2] = {getSolution(X[E[e][0]]), getSolution(Y[E[e][0]])},
                           b[2] = {getSolution(X[E[e][1]]), getSolution(Y[E[e][1]])},
                           c[2] = {getSolution(X[E[f][0]]), getSolution(Y[E[f][0]])},
                           d[2] = {getSolution(X[E[f][1]]), getSolution(Y[E[f][1]])};

                    itersect = !(side(a, c, d) == side(b, c, d) || side(a, b, c) == side(a, b, d));

                    if (itersect)
                        cout << E[e][0] + 1 << "-" << E[e][1] + 1 << " intersects with " << E[f][0] + 1 << "-"
                             << E[f][1] + 1 << endl;
                }
        } catch (GRBException &e) {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) { cout << "Error during callback" << endl; }
    }
};

string itos(int i) {
    stringstream s;
    s << i;
    return s.str();
}

array<int, 2> *read_graph(int *n, int *m) {
    cin >> *n >> *m;
    auto ret = new array<int, 2>[*m];
    for (int i = 0; i < *m; i++) cin >> ret[i][0] >> ret[i][1];
    return ret;
}

void method() {
    // Read the graph edges data: -----------------------------------------
    int n, m;
    array<int, 2> *E = read_graph(&n, &m);

    // Sort the edges in lexicographic order:
    auto comp = [](const array<int, 2> &u, const array<int, 2> &v) {
        return u[0] != v[0] ? u[0] < v[0] : u[1] < v[1];
    };
    sort(E, E + m, comp);

    auto X = new GRBVar[n], Y = new GRBVar[n];// coordinate for each vertex
    GRBVar width, height;                     // dimensions of the drawing
    auto edge_len = new GRBVar[m];

    int M = (int) pow(n, 2);// a huge number
    double time_limit = 180;// 3 min
    try {
        auto *env = new GRBEnv();
        env->set(GRB_IntParam_Seed, SEED);
        env->set(GRB_DoubleParam_TimeLimit, time_limit);

        GRBModel model = GRBModel(*env);
        model.set(GRB_StringAttr_ModelName, "2D Grid Drawing");
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.set(GRB_IntParam_NonConvex, 2);

        // Must set LazyConstraints parameter when using lazy constraints:
        // model.set(GRB_IntParam_LazyConstraints, 1);

        // Set callback function:
        // IntersectElim cb = IntersectElim(m, E, X, Y, edge_len);
        // model.setCallback(&cb);

        // Focus primarily on feasibility of the relaxation:
        // model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
        // model.set(GRB_IntParam_Cuts, GRB_CUTS_AGGRESSIVE);
        // model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);

        // Variable creation: _____________________________________________________________________
        width = model.addVar(ceil(sqrt(n)) - 1, M, 0.0, GRB_INTEGER, "w");
        height = model.addVar(1.0, M, 0.0, GRB_INTEGER, "h");
        model.setObjective(width * height);
        for (int i = 0; i < n; i++) {
            X[i] = model.addVar(0.0, M, 0.0, GRB_INTEGER, "x_" + itos(i));
            Y[i] = model.addVar(0.0, M, 0.0, GRB_INTEGER, "y_" + itos(i));
        }
        for (int e = 0; e < m; e++)
            edge_len[e] = model.addVar(1.0, sqrt(2 * M * M), 0.0, GRB_CONTINUOUS, "el_" + itos(e));
        model.update();// run update to use model inserted variables

        // Constraint creation: ___________________________________________________________________
        model.addConstr(height <= width, "lied_down_box");

        // All coordinates must be withing the bounding box:
        for (int i = 0; i < n; i++) {
            model.addConstr(X[i] <= width, "x_limit_" + itos(i));
            model.addConstr(Y[i] <= height, "y_limit_" + itos(i));
        }

        // All coordinates must be different from each other:
        for (int u = 0; u < n; u++)
            for (int v = u + 1; v < n; v++)
                model.addQConstr((X[u] - X[v]) * (X[u] - X[v]) + (Y[u] - Y[v]) * (Y[u] - Y[v]) >= 1);

        // Relate the coordinates with the edges' length:
        for (int e = 0; e < m; e++) {
            int u = E[e][0], v = E[e][1];
            GRBQuadExpr edge_len2 = (X[u] - X[v]) * (X[u] - X[v]) + (Y[u] - Y[v]) * (Y[u] - Y[v]);
            model.addQConstr(edge_len[e] * edge_len[e] == edge_len2);
        }

        // Prohibit intersections between edges with a common endpoint:
        for (int e = 0, e_final; e < m; e = e_final) {
            auto U_x = X[E[e][0]], U_y = Y[E[e][0]];
            for (e_final = e; e_final < m && E[e_final][0] == E[e][0]; e_final++) {}

            for (int vi = e; vi < e_final; vi++)
                for (int wi = vi + 1; wi < e_final; wi++) {
                    auto V_x = X[E[vi][1]], V_y = Y[E[vi][1]], W_x = X[E[wi][1]], W_y = Y[E[wi][1]];
                    // Ensure that they are not collinear (i.e., cos = v.w/(|v||w|) < 1):
                    auto dot_prod = (V_x - U_x) * (W_x - U_x) + (V_y - U_y) * (W_y - U_y);
                    auto name = "int_" + itos(vi) + "_" + itos(wi);
                    model.addQConstr(dot_prod <= edge_len[vi] * edge_len[wi] - EPS, name);
                }
        }

        model.update();// run update before optimize
        model.optimize();
        if (model.get(GRB_IntAttr_SolCount) == 0) throw GRBException("Could not obtain a solution!", -1);

        /*cout << "Edges' length:" << endl;
        for (int e = 0; e < m; e++)
            cout << e + 1 << "th edge length - " << edge_len[e].get(GRB_DoubleAttr_X) << endl;*/
        cout << "Points coordinates:" << endl;
        for (int i = 0; i < n; i++)
            cout << i + 1 << "th vtx - (" << (int) X[i].get(GRB_DoubleAttr_X) << "; "
                 << (int) Y[i].get(GRB_DoubleAttr_X) << ")" << endl;
        cout << "- Width " << width.get(GRB_DoubleAttr_X) << endl;
        cout << "- Height " << height.get(GRB_DoubleAttr_X) << endl;
        cout << "- Area Cost " << model.getObjective().getValue() << endl;
    } catch (GRBException &e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) { cout << "Error during callback" << endl; }
    delete[] X;
    delete[] Y;
}

int main() {
    method();
    return 0;
}
