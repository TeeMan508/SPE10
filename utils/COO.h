#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>

class COO {
private:
    std::vector<int> ia = {};
    std::vector<int> ja = {};
    std::vector<double> a = {};


public:
    COO() {
    }

    int len_mat() {
        if (ia[ia.size() - 1] > ja[ja.size() - 1]) { return ia[ia.size() - 1] + 1; }
        else return ja[ja.size() - 1] + 1;
    }


    int insert_val(int row, int col, double value) {
        if ((ia.empty()) ||
            ((row >= ia[ia.size() - 1]) && (col > ja[ja.size() - 1]))) {
            ia.push_back(row);
            ja.push_back(col);
            a.push_back(value);
        } else {
/*
 * Find the elements of the same row as the element being inserted into the matrix.
 *
 * If this is the first element in a specific row, the two iterators returned by equal_range()
 * point to the first element of the next larger row.
 *
 * If there are already other elements in the same row, the two iterators returned by equal_range()
 * point to the first element of the row and the first element of the next larger row.
 *
 * Using the iterators also calculate indices to the elements returned by equal_range().
 * These are used to index the corresponding elements in the other two vectors representing
 * the sparse matrix (ja and values).
 */
            const auto p = std::equal_range(ia.begin(), ia.end(), row);
            const auto index_of_first = p.first - ia.begin();
            const auto index_of_last = p.second - ia.begin();

/*
 * Create iterators to point to the corresponding elements in ja.
 */
            const auto first = next(ja.begin(), index_of_first);
            const auto last = next(ja.begin(), index_of_last);

/*
 * Find the correct position where the new element must be inserted and perform the corresponding
 * insertions into the three vectors representing the sparse matrix.
 */
            auto col_pos_it = upper_bound(first, last, col);
            auto pos = col_pos_it - ja.begin();
            ja.insert(col_pos_it, col);

            auto row_pos_it = next(ia.begin(), pos);
            ia.insert(row_pos_it, row);

            auto val_pos_it = next(a.begin(), pos);
            a.insert(val_pos_it, value);
        }
    }

    double operator()(int index) {
        return a[index];
    }


    double operator()(int row, int col) {
        if (row > len_mat() - 1 || row > len_mat() - 1) {
            std::cerr << "out of range";
            throw;
        }
        else {
            int i = 0;
            while (ia[i] != row && i <= ia.size() - 1) {
                i++;
            }
            if (i == ia.size()) { return 0; }
            while (ja[i] != col && ia[i] == row) {
                i++;
            }
            if (ia[i] != row) { return 0; }
            return a[i];
        }
    }

    void print_coo() {
        std::cout << "ia = [";
        for (int i = 0; i < ia.size() - 1; ++i) {
            std::cout << ia[i] << ", ";
        }
        std::cout << ia[ia.size() - 1] << "]" << std::endl;

        std::cout << "ja = [";
        for (int i = 0; i < ja.size() - 1; ++i) {
            std::cout << ja[i] << ", ";
        }
        std::cout << ja[ja.size() - 1] << "]" << std::endl;

        std::cout << "a = [";
        for (int i = 0; i < a.size() - 1; ++i) {
            std::cout << a[i] << ", ";
        }
        std::cout << a[a.size() - 1] << "]" << std::endl;
    }

    void print_mat() {
//        std::cout << "len = " << len_mat() << std::endl;
        for (int i = 0; i < len_mat(); ++i) {
            for (int j = 0; j < len_mat(); ++j) {

                auto k = operator()(i,j);
                std::cout<<std::setw(5)<<this->operator()(i,j);
            }
            std::cout<<std::endl;
        }
    }
};