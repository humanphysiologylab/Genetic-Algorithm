#include "optimization_problem.h"
#include <iomanip>

BlockOfTable::BlockOfTable(int * model_indices, Block xblock)
    : Base(xblock), model_indices(model_indices)
{}
int BlockOfTable::get_model_pos(int col) const
{
    return model_indices[col];
}
Table::Table(int m, int n, const std::vector<std::string> & header)
    : Base(Base::Constant(m, n, 0)), header(header)
{    
}

template<typename T>
std::vector<T> operator+(const std::vector<T> & v1, const std::vector<T> & v2)
{
    std::vector<T> res(v1.size() + v2.size());
    res.insert(res.end(), v1.begin(), v1.end());
    res.insert(res.end(), v2.begin(), v2.end());
    return res;
}

TableSplit::TableSplit(int m, int n, std::vector<int> states_model_indices,
            std::vector<int> alg_model_indices,
            const std::vector<std::string> & states_names,
            const std::vector<std::string> & alg_names)
        : Base(m, n, states_names + alg_names)
{
    states_cols = states_model_indices.size();
    alg_cols = alg_model_indices.size();
    assert(n == states_cols + alg_cols);
    model_indices = new int [states_model_indices.size() + alg_model_indices.size()];
    for (unsigned i = 0; i < states_model_indices.size(); i++) {
        model_indices[i] = states_model_indices[i];
    }
    for (unsigned i = 0; i < alg_model_indices.size(); i++) {
        model_indices[i + states_model_indices.size()] = alg_model_indices[i];
    }
}
TableSplit::~TableSplit()
{
    delete [] model_indices;
}
BlockOfTable TableSplit::get_algebraic()
{
    return BlockOfTable(model_indices + states_cols, this->rightCols(alg_cols));
}
BlockOfTable TableSplit::get_states()
{
    return BlockOfTable(model_indices, this->leftCols(states_cols));
}
void Table::export_csv(const std::string & filename)
{
    std::ofstream file;
    file.open(filename + ".csv");
    if (!file.is_open())
        throw(filename + " cannot be opened");
    for (const auto & x: header)
        file << x << " ";
    file << std::endl;
    file << std::scientific << std::setprecision(12);
    for (int i = 0; i < Base::rows(); i++) {
        for(int j = 0; j < Base::cols(); j++) {
            file << Base::operator()(i, j) << " ";
        }
        file << std::endl;
    }
}



int halfheight_index(const std::vector<double> & v)
{
    /* return depolarization halfheight position in v 
     * It is assumed that action potential is completely inside the single period
     * represented by v
     */
    
    assert(v.size() > 0);
    double max = v[0], min = v[0];
    //find max and min
    for (int i = 0; i < v.size(); i++) {
        if (v[i] < min) {
            min = v[i];
        } else if (v[i] > max) {
            max = v[i];
        }
    }
    const double hh = (max + min) / 2;
    
    //find i : v[i] <= hh <= v[i + 1]
    for (int i = 0; i < v.size(); i++) {
        if (v[i] < hh)
            continue;
        //so v[i] >= hh
        if (i == 0)
            break;
        //find closest among i-1 and i
        const double a1 = std::fabs(v[i-1] - hh);
        const double a2 = std::fabs(v[i] - hh);
        if (a1 < a2)
            return i - 1;
        else
            return i;
    }
    return -1;//smth wrong, or just v[-1] <= hh <= v[0] xx
}
