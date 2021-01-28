#include "optimization_problem.h"

BlockOfTable::BlockOfTable(int * model_indices, Block xblock)
    : Base(xblock), model_indices(model_indices)
{}
int BlockOfTable::get_model_pos(int col) const
{
    return model_indices[col];
}
Table::Table(int m, int n, std::vector<int> states_model_indices,
            std::vector<int> alg_model_indices,
            const std::vector<std::string> & states_names,
            const std::vector<std::string> & alg_names)
        : Base(m, n)
{
    states_cols = states_model_indices.size();
    alg_cols = alg_model_indices.size();
    assert(n == states_cols + alg_cols);
    model_indices = new int [states_model_indices.size() + alg_model_indices.size()];
    for (unsigned i = 0; i < states_model_indices.size(); i++) {
        model_indices[i] = states_model_indices[i];
        header.push_back(states_names[i]);
    }
    for (unsigned i = 0; i < alg_model_indices.size(); i++) {
        model_indices[i + states_model_indices.size()] = alg_model_indices[i];
        header.push_back(alg_names[i]);
    }
}
Table::~Table()
{
    delete [] model_indices;
}
BlockOfTable Table::get_algebraic()
{
    return BlockOfTable(model_indices + states_cols, this->rightCols(alg_cols));
}
BlockOfTable Table::get_states()
{
    return BlockOfTable(model_indices, this->leftCols(states_cols));
}
void Table::export_csv(const std::string & filename)
{
    std::ofstream file;
    file.open(filename);
    if (!file.is_open())
        throw(filename + " cannot be opened");
    for (const auto & x: header)
        file << x << " ";
    file << std::endl;
   
    for (int i = 0; i < Base::rows(); i++) {
        for(int j = 0; j < Base::cols(); j++) {
            file << Base::operator()(i, j) << " ";
        }
        file << std::endl;
    }
}
