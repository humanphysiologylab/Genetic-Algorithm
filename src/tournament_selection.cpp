#include "tournament_selection.h"
#include <forward_list>
#include <cstdlib>

void tournament_basic_half(int *mpool, int mpool_size, const std::vector<std::pair<double, int>> & sd_n_index, int number_of_ignored_losers)
{
    /* sd_n_index is expected to be sorted by sd in increasing order.
     * tournament_basic_half fills in mpool with indices of organisms in mating pool according to the tournament selection without replacement.
     * At the end of sd_n_index the worst organisms are stored, and bottom number_of_ignored_losers are totally excluded 
     * and will not be in the mpool.
     * 
     */
     
    /* We use forward_list as an array we will run through
     * picking winners and kicking out random losers.
     */
    int list_size = sd_n_index.size() - number_of_ignored_losers;
    std::forward_list<int> list_of_indices(list_size);
    
    /* First, copy indices from sd_n_index
     * to the list wrt the order
     */
    auto a = list_of_indices.begin();
    for (int i = 0; i < list_size; i++) {
        *a = sd_n_index[i].second;
        a++;
    }
    
    /* 'a' points at the current winner, since our list is sorted,
     * every organism after 'a' is inferior to it.
     */
    a = list_of_indices.begin();
    for (int i = 0; i < mpool_size; i++) {
        assert(a != list_of_indices.end());
        /* save current winner
         */
        mpool[i] = *a;
        
        /* now find a loser to the right 
         */
        auto eraser = a;
        const int eraser_index = rand() % (list_size - 1 - i); //it works correctly
        for (int j = 0; j < eraser_index; j++)
            eraser++;
        /* Now eraser points at the element in the list which
         * is to the left for the picked loser.
         * We need to erase the loser from the list
         * and move to the winner of the next pair.
         */
        list_of_indices.erase_after(eraser);
        list_size--;
        a++;
    }
    
}


void tournament_selection(int *mpool, int mpool_size, const std::vector<std::pair<double, int>> & sd_n_index, int number_of_ignored_losers)
{
    /* sd_n_index is expected to be sorted by sd in increasing order.
     * tournament_selection fills in mpool with indices of organisms in mating pool.
     * At the end of sd_n_index the worst organisms are stored, and bottom number_of_ignored_losers are totally excluded 
     * and will not be in the mpool.
     * 
     * each call of tournament_basic_half fills half of mating pool
     */
    assert(mpool_size % 2 == 0);
    tournament_basic_half(mpool                 , mpool_size / 2, sd_n_index, number_of_ignored_losers);
    tournament_basic_half(mpool + mpool_size / 2, mpool_size / 2, sd_n_index, number_of_ignored_losers);
}
