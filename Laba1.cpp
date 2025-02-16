#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <cmath>
#include <limits>
#include <algorithm>
#define PRINT
using namespace std;

struct Constraint {
    vector<double> coefficients;
    double rhs;
    char type;
};

struct CanonicalConstraint {
    vector<double> coefficients;
    double rhs;
};

struct SimplexResult {
    vector<double> solution;
    double optimal;
    bool unbounded;
    bool infeasible;
};

// Преобразование задачи в каноническую форму
auto CanonForm() {
    vector<double> original_obj = { -2, -3, 1, -4, -1 };
    vector<bool> is_free = { false, false, false, false, true };
    vector<Constraint> original_constraints = {
        {{1, 1, 1, 0, 0}, 10, '='},
        {{2, 0, 0, -1, 1}, 5, '='},
        {{0, 0, 1, 2, -1}, 3, '='},
        {{1, -2, 1, 0, 0}, 8, '<'},
        {{0, 0, 0, 1, 1}, 2, '>'}
    };

    map<int, pair<int, int>> split_vars;
    vector<double> canonical_obj;
    vector<string> var_names;

    for (size_t i = 0; i < original_obj.size(); ++i) {
        if (!is_free[i]) {
            canonical_obj.push_back(original_obj[i]);
            var_names.push_back("x" + to_string(i + 1));
        }
        else {
            int plus_idx = canonical_obj.size();
            canonical_obj.push_back(original_obj[i]);
            var_names.push_back("x" + to_string(i + 1) + "+");

            int minus_idx = canonical_obj.size();
            canonical_obj.push_back(-original_obj[i]);
            var_names.push_back("x" + to_string(i + 1) + "-");

            split_vars[i] = { plus_idx, minus_idx };
        }
    }

    vector<CanonicalConstraint> canonical_constraints;
    int slack_counter = 1;
    for (const auto& con : original_constraints) {
        CanonicalConstraint new_con;
        new_con.rhs = con.rhs;
        new_con.coefficients.resize(canonical_obj.size(), 0.0);

        for (size_t i = 0; i < con.coefficients.size(); ++i) {
            if (!is_free[i]) {
                new_con.coefficients[i] += con.coefficients[i];
            }
            else {
                auto& indices = split_vars[i];
                new_con.coefficients[indices.first] += con.coefficients[i];
                new_con.coefficients[indices.second] -= con.coefficients[i];
            }
        }

        if (con.type != '=') {
            string var_name = "s" + to_string(slack_counter++);
            var_names.push_back(var_name);
            canonical_obj.push_back(0.0);
            new_con.coefficients.push_back(con.type == '<' ? 1.0 : -1.0);
        }
        canonical_constraints.push_back(new_con);
    }
#ifdef PRINT
    // Вывод канонической формы
    cout << "Canin Form:\n";
    cout << "Aim Func: Z = ";
    for (size_t i = 0; i < canonical_obj.size(); ++i) {
        if (i > 0) cout << " + ";
        cout << canonical_obj[i] << "*" << var_names[i];
    }
    cout << "\n\nOgr:\n";

    for (const auto& con : canonical_constraints) {
        for (size_t i = 0; i < con.coefficients.size(); ++i) {
            if (i > 0) cout << " + ";
            cout << con.coefficients[i] << "*" << var_names[i];
        }
        cout << " = " << con.rhs << endl;
    }

    cout << "\nAll variables >= 0\n";
#endif
    return make_tuple(canonical_constraints, canonical_obj, var_names);
}

void Pivot(vector<vector<double>>& T, int m, int n, int enter, int leave) {
    double pivot = T[leave][enter];
    for (int j = 0; j <= n; j++) {
        T[leave][j] /= pivot;
    }
    for (int i = 0; i <= m; i++) {
        if (i != leave) {
            double factor = T[i][enter];
            for (int j = 0; j <= n; j++) {
                T[i][j] -= factor * T[leave][j];
            }
        }
    }
}

SimplexResult SimplexMethod(vector<vector<double>> T, int m, int n) {
    const double eps = 1e-9;
    while (true) {
        int enter = -1;
        for (int j = 0; j < n; j++) {
            if (T[m][j] < -eps) {
                enter = j;
                break;
            }
        }
        if (enter == -1) break;

        int leave = -1;
        double minRatio = numeric_limits<double>::infinity();
        for (int i = 0; i < m; i++) {
            if (T[i][enter] > eps) {
                double ratio = T[i][n] / T[i][enter];
                if (ratio < minRatio || (fabs(ratio - minRatio) < eps && i < leave)) {
                    minRatio = ratio;
                    leave = i;
                }
            }
        }
        if (leave == -1) {
            return { {}, 0.0, true, false };
        }
        Pivot(T, m, n, enter, leave);
    }

    vector<double> solution(n, 0.0);
    for (int j = 0; j < n; j++) {
        bool is_basic = false;
        int basic_row = -1;
        for (int i = 0; i < m; i++) {
            if (fabs(T[i][j]) > eps) {
                if (!is_basic) {
                    is_basic = true;
                    basic_row = i;
                }
                else {
                    is_basic = false;
                    break;
                }
            }
        }
        if (is_basic && basic_row != -1) {
            solution[j] = T[basic_row][n];
        }
    }
    double optimal = T[m][n];
    return { solution, optimal, false, false };
}

int main() {
    auto [canonical_constraints, canonical_obj, var_names] = CanonForm();
    vector<vector<double>> A;
    vector<double> b;
    for (const auto& con : canonical_constraints) {
        A.push_back(con.coefficients);
        b.push_back(con.rhs);
    }

    int m = A.size();
    int n = canonical_obj.size();
    vector<vector<double>> T(m + 1, vector<double>(n + 1, 0.0));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            T[i][j] = (j < A[i].size()) ? A[i][j] : 0.0;
        }
        T[i][n] = b[i];
    }
    for (int j = 0; j < n; j++) {
        T[m][j] = -canonical_obj[j];
    }

    SimplexResult result = SimplexMethod(T, m, n);
    if (result.unbounded) {
        cout << "The problem is unbounded." << endl;
    }
    else if (result.infeasible) {
        cout << "The problem is infeasible." << endl;
    }
    else {
        cout << "Optimal value: " << result.optimal << endl;
        cout << "Solution: ";
        for (double val : result.solution) {
            cout << val << " ";
        }
        cout << endl;
    }

    return 0;
}
