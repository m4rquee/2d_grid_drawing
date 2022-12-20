#include "gurobi_c++.h"
#include <array>
#include <cmath>
#include <sstream>

using namespace std;

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
class InterElim : public GRBCallback {
public:
    int m;
    GRBVar *X, *Y;
    array<int, 2> *E;

    InterElim(int _m, array<int, 2> *_E, GRBVar *_X, GRBVar *_Y) : m(_m), X(_X), Y(_Y), E(_E) {}

protected:
    bool itersect(int e, int f) const {
        double a[2] = {X[E[e][0]].get(GRB_DoubleAttr_Obj), Y[E[e][0]].get(GRB_DoubleAttr_Obj)},
               b[2] = {X[E[e][1]].get(GRB_DoubleAttr_Obj), Y[E[e][1]].get(GRB_DoubleAttr_Obj)},
               c[2] = {X[E[f][0]].get(GRB_DoubleAttr_Obj), Y[E[f][0]].get(GRB_DoubleAttr_Obj)},
               d[2] = {X[E[f][1]].get(GRB_DoubleAttr_Obj), Y[E[f][1]].get(GRB_DoubleAttr_Obj)};

        if (side(a, c, d) == side(b, c, d) || side(a, b, c) == side(a, b, d)) return false;
        return true;
    }

    void callback() override {
        try {
            if (where == GRB_CB_MIPSOL) {
                // Found an integer feasible solution - does it have no edge intersection?
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

    auto XorYdiff = new GRBVar *[n];// if the [u][v] vertices coordinates differ in X or Y
    auto Xu_gt_Xv = new GRBVar *[n];// given that u differs from v in the X coordinate, then Xu < Xv?
    auto Yu_gt_Yv = new GRBVar *[n];// given that u differs from v in the Y coordinate, then Yu < Yv?

    int M = (int) pow(n, 3);// a huge number
    double time_limit = 180;// 3 min
    try {
        auto *env = new GRBEnv();
        env->set(GRB_DoubleParam_TimeLimit, time_limit);

        GRBModel model = GRBModel(*env);
        model.set(GRB_StringAttr_ModelName, "2D Grid Drawing");
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.set(GRB_IntParam_NonConvex, 2);

        // Must set LazyConstraints parameter when using lazy constraints:
        model.set(GRB_IntParam_LazyConstraints, 0);

        // Set callback function:
        InterElim cb = InterElim(m, E, X, Y);
        model.setCallback(&cb);

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
        for (int i = 0; i < n; i++) {
            XorYdiff[i] = new GRBVar[n], Xu_gt_Xv[i] = new GRBVar[n], Yu_gt_Yv[i] = new GRBVar[n];
            for (int j = i + 1; j < n; j++) {
                XorYdiff[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_or_y_diff_" + itos(i) + itos(j));
                Xu_gt_Xv[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "xu_less_xv_" + itos(i) + itos(j));
                Yu_gt_Yv[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "yu_less_yv_" + itos(i) + itos(j));
            }
        }
        model.update();// run update to use model inserted variables

        // Constraint creation: ___________________________________________________________________
        model.addConstr(height <= width, "lied_down_box");

        // All coordinates must be withing the bounding box:
        for (int i = 0; i < n; i++) {
            model.addConstr(X[i] <= width, "x_limit_" + itos(i));
            model.addConstr(Y[i] <= height, "y_limit_" + itos(i));
        }

        // All coordiates must be different from each other:
        for (int u = 0; u < n; u++)
            for (int v = u + 1; v < n; v++) {
                model.addConstr(X[u] - X[v] <= XorYdiff[u][v] * M + Xu_gt_Xv[u][v] * M - 1);
                model.addConstr(X[u] - X[v] >= -XorYdiff[u][v] * M + (Xu_gt_Xv[u][v] - 1) * M + 1);

                model.addConstr(Y[u] - Y[v] <= (1 - XorYdiff[u][v]) * M + Yu_gt_Yv[u][v] * M - 1);
                model.addConstr(Y[u] - Y[v] >= -(1 - XorYdiff[u][v]) * M + (Yu_gt_Yv[u][v] - 1) * M + 1);
            }

        model.update();// run update before optimize
        model.optimize();
        if (model.get(GRB_IntAttr_SolCount) == 0) throw GRBException("Could not obtain a solution!", -1);

        cout << "Points coordinates:" << endl;
        for (int i = 0; i < n; i++)
            cout << i + 1 << "th vtx - (" << (int) X[i].get(GRB_DoubleAttr_X) << "; "
                 << (int) Y[i].get(GRB_DoubleAttr_X) << ")" << endl;
        cout << "- Width " << width.get(GRB_DoubleAttr_X) << endl;
        cout << "- Height " << height.get(GRB_DoubleAttr_X) << endl;
        cout << "- Area Cost " << model.getObjective().getValue() << endl;

        bool itersect;
        for (int e = 0; e < m; e++)
            for (int f = e + 1; f < m; f++) {
                double a[2] = {X[E[e][0]].get(GRB_DoubleAttr_X), Y[E[e][0]].get(GRB_DoubleAttr_X)},
                       b[2] = {X[E[e][1]].get(GRB_DoubleAttr_X), Y[E[e][1]].get(GRB_DoubleAttr_X)},
                       c[2] = {X[E[f][0]].get(GRB_DoubleAttr_X), Y[E[f][0]].get(GRB_DoubleAttr_X)},
                       d[2] = {X[E[f][1]].get(GRB_DoubleAttr_X), Y[E[f][1]].get(GRB_DoubleAttr_X)};

                itersect = !(side(a, c, d) == side(b, c, d) || side(a, b, c) == side(a, b, d));

                if (itersect)
                    cout << E[e][0] + 1 << "-" << E[e][1] + 1 << " intersects with " << E[f][0] + 1 << "-"
                         << E[f][1] + 1 << endl;
            }
    } catch (GRBException &e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) { cout << "Error during callback" << endl; }
    for (int i = 0; i < n; i++) {
        delete[] XorYdiff[i];
        delete[] Xu_gt_Xv[i];
        delete[] Yu_gt_Yv[i];
    }
    delete[] X;
    delete[] Y;
    delete[] XorYdiff;
    delete[] Xu_gt_Xv;
    delete[] Yu_gt_Yv;
}

int main() {
    method();
    return 0;
}
