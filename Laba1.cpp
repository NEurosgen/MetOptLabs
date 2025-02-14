#include <iostream>
#include <vector>
#include <string>
#include <map>

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

auto CreateMatrix(std::vector<CanonicalConstraint> Mat,int max) {
    std::vector<std::vector<double>>  ans(Mat.size(), std::vector<double>(max));
    std::vector<double> vectorans(Mat.size());
    for (int i = 0; i < Mat.size(); ++i) {
        vectorans[i] = Mat[i].rhs;
        for (int j = 0; j < max; ++j) {
            if (j >= Mat[i].coefficients.size()) {
                ans[i][j] = 0;
                
            }
            else {
                ans[i][j] = Mat[i].coefficients[j];
            }

        }

    }
    return std::make_pair( ans, vectorans);

}

auto CanonForm(){
   // Исходные данные пока в таком виде 
    vector<double> original_obj = { 2, 3, -1, 4, 1 };
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

    for (int i = 0; i < original_obj.size(); ++i) {
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

        // Обработка коэффициентов
        for (int i = 0; i < con.coefficients.size(); ++i) {
            if (!is_free[i]) {
                new_con.coefficients[i] += con.coefficients[i];
            }
            else {
                auto& indices = split_vars[i];
                new_con.coefficients[indices.first] += con.coefficients[i];
                new_con.coefficients[indices.second] -= con.coefficients[i];
            }
        }

        // Добавление slack/surplus переменных
        if (con.type != '=') {
            string var_name = "s" + to_string(slack_counter++);
            var_names.push_back(var_name);
            canonical_obj.push_back(0.0);

            if (con.type == '<') {
                new_con.coefficients.push_back(1.0);
            }
            else {
                new_con.coefficients.push_back(-1.0);
            }
        }

        canonical_constraints.push_back(new_con);
    }

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



    auto [MatrixA, VectorB] = CreateMatrix(canonical_constraints, var_names.size());
    return std::make_pair(MatrixA, VectorB);

   
}

int main() {
    auto [MatrixA, VectorB] = CanonForm();
}