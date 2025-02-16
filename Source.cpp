#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;

struct Constraint {
    vector<double> coefficients;
    double rhs;
    char type; // '<' for <=, '>' for >=, '=' for =
};

struct DualConstraint {
    vector<double> coefficients;
    double rhs;
    char type; // '>' for >=
};

struct SimplexResult {
    vector<double> solution;
    double optimal;
    bool unbounded;
    bool infeasible;
};

// Выполняет шаг симплекс-метода, преобразуя таблицу
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

// Реализует симплекс-метод для поиска оптимального решения
SimplexResult SimplexMethod(vector<vector<double>> T, int m, int n) {
    const double eps = 1e-9;
    int iteration = 0;
    while (true) {
        cout << "Iteration " << ++iteration << ":\n";
        for (const auto& row : T) {
            for (double val : row) {
                cout << val << "\t";
            }
            cout << endl;
        }
        cout << endl;

        int enter = -1;
        for (int j = 0; j < n; j++) {
            if (T[m][j] < -eps) {  // Для максимизации ищем отрицательный коэффициент
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
                if (ratio < minRatio) {
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
    double optimal = -T[m][n]; // Для максимизации меняем знак
    return { solution, optimal, false, false };
}

// Преобразует задачу в двойственную форму и выводит её
void DualForm(const vector<double>& obj, const vector<Constraint>& constraints) {
    int m = constraints.size();
    int n = obj.size();

    vector<double> dual_obj(m);
    vector<DualConstraint> dual_constraints(n);

    // Формирование целевой функции двойственной задачи
    for (int i = 0; i < m; ++i) {
        dual_obj[i] = constraints[i].rhs;
    }

    // Формирование ограничений двойственной задачи
    for (int j = 0; j < n; ++j) {
        dual_constraints[j].coefficients.resize(m);
        for (int i = 0; i < m; ++i) {
            dual_constraints[j].coefficients[i] = constraints[i].coefficients[j];
        }
        dual_constraints[j].rhs = obj[j];
        dual_constraints[j].type = '>';  // Ограничение двойственной задачи становится >=
    }

    // Вывод двойственной задачи
    cout << "Dual Form:\n";
    cout << "Minimize W = ";
    for (int i = 0; i < m; ++i) {
        if (i > 0) cout << " + ";
        cout << dual_obj[i] << "*y" << (i + 1);
    }
    cout << "\n\nSubject to:\n";
    for (const auto& con : dual_constraints) {
        for (size_t i = 0; i < con.coefficients.size(); ++i) {
            if (i > 0) cout << " + ";
            cout << con.coefficients[i] << "*y" << (i + 1);
        }
        cout << " >= " << con.rhs << endl;
    }
    cout << "\nAll variables y >= 0\n";

    // Создание симплекс-таблицы для двойственной задачи
    vector<vector<double>> T(m + 1, vector<double>(n + 1, 0.0));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            T[i][j] = dual_constraints[j].coefficients[i];
        }
        T[i][n] = dual_obj[i];
    }
    for (int j = 0; j < n; j++) {
        T[m][j] = -obj[j];  // Для максимизации знак меняется
    }

    // Запуск симплекс-метода
    SimplexResult result = SimplexMethod(T, m, n);
    if (result.unbounded) {
        cout << "The dual problem is unbounded." << endl;
    }
    else if (result.infeasible) {
        cout << "The dual problem is infeasible." << endl;
    }
    else {
        cout << "Optimal value of dual problem: " <<-1* result.optimal << endl;
        cout << "Dual solution: ";
        for (double val : result.solution) {
            cout << val << " ";
        }
        cout << endl;
    }
}

int main() {
    vector<double> original_obj = { 2, 3, -1, 4, 1 };
    vector<Constraint> original_constraints = {
        {{1, 1, 1, 0, 0}, 10, '>'},
        {{2, 0, 0, -1, 1}, 5, '>'},
        {{0, 0, 1, 2, -1}, 3, '>'},
        {{1, -2, 1, 0, 0}, 8, '>'},
        {{0, 0, 0, 1, 1}, 2, '>'}
    };

    DualForm(original_obj, original_constraints);

    return 0;
}
