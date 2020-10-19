#include "tournament_selection.h"
#include <forward_list>
#include <cstdlib>

void tournament_basic(int *mpool, int mpool_size, std::vector<std::pair<double, int>> & sd_n_index, int number_of_ignored_losers)
{
    int list_size = sd_n_index.size() - number_of_ignored_losers;
    std::forward_list<int> list_of_indices(list_size);
    
    auto a = list_of_indices.begin();
    for (int i = 0; i < list_size; i++) {
        *a = sd_n_index[i].second;
        a++;
    }
    
    a = list_of_indices.begin();
    for (int i = 0; i < mpool_size; i++) {
        assert(a != list_of_indices.end());
        mpool[i] = *a;
        
        auto eraser = a;
        const int eraser_index = rand() % (list_size - 1 - i); //it works but needs explanation
        for (int j = 0; j < eraser_index; j++)
            eraser++;
    
        list_of_indices.erase_after(eraser);
        list_size--;
        a++;
    }
    
}


void tournament_selection(int *mpool, int mutant_number, std::vector<std::pair<double, int>> & sd_n_index, int number_of_ignored_losers)
{
    assert(mutant_number % 2 == 0);
    tournament_basic(mpool                    , mutant_number / 2, sd_n_index, number_of_ignored_losers);
    tournament_basic(mpool + mutant_number / 2, mutant_number / 2, sd_n_index, number_of_ignored_losers);
}
