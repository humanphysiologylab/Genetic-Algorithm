#ifndef GA_TOURNAMENT_SELECTION_H
#define GA_TOURNAMENT_SELECTION_H

#include <vector>
#include <forward_list>
#include <random>
#include <cassert>
#include <algorithm>

template<typename InitializedRandomGenerator>
class TournamentSelection
{
    InitializedRandomGenerator rg;

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
            const int eraser_index = std::uniform_int_distribution<int>(0, list_size - 2 - i)(rg); //it works correctly
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

public:
    TournamentSelection(InitializedRandomGenerator rg)
    : rg(rg)
    {}

    void operator()(int *mpool, int mpool_size, std::vector<std::pair<double, int>> & sd_n_index, int number_of_ignored_losers)
    {
        //sort by error increasing
        std::sort(sd_n_index.begin(), sd_n_index.end(),
                      [](const std::pair<double, int> &left_element, const std::pair<double, int> &right_element) {
                          return left_element.first < right_element.first;
                      });
        /* sd_n_index is expected to be sorted by sd in increasing order.
         * tournament selection fills in mpool with indices of organisms in mating pool.
         * At the end of sd_n_index the worst organisms are stored, and bottom number_of_ignored_losers are totally excluded 
         * and will not be in the mpool.
         * 
         * each call of tournament_basic_half fills half of mating pool
         */
        assert(mpool_size % 2 == 0);
        tournament_basic_half(mpool                 , mpool_size / 2, sd_n_index, number_of_ignored_losers);
        tournament_basic_half(mpool + mpool_size / 2, mpool_size / 2, sd_n_index, number_of_ignored_losers);
        
        //we also need to shuffle mpool because it is sorted!
        for (int i = 0; i < mpool_size; i++) {
            int j = std::uniform_int_distribution<int>(i, mpool_size - 1)(rg);
            std::swap(mpool[i], mpool[j]);
        }
    }
};





template<typename InitializedRandomGenerator>
class TournamentSelectionFast
{
    InitializedRandomGenerator rg;

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
            const int eraser_index = std::uniform_int_distribution<int>(0, list_size - 2 - i)(rg); //it works correctly
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

public:
    TournamentSelectionFast(InitializedRandomGenerator rg)
    : rg(rg)
    {}

    void operator()(int *mpool, int mpool_size, std::vector<std::pair<double, int>> & sd_n_index, int number_of_ignored_losers)
    {
        TODO!
        /* FAST version!
         * tournament selection fills in mpool with indices of organisms in mating pool.
         * each call of tournament_basic_half fills half of mating pool.
         * no guarantees of ignoring losers
         */
        assert(mpool_size % 2 == 0);
        tournament_basic_half(mpool                 , mpool_size / 2, sd_n_index, number_of_ignored_losers);
        tournament_basic_half(mpool + mpool_size / 2, mpool_size / 2, sd_n_index, number_of_ignored_losers);
        
        //we also need to shuffle mpool because it is sorted!
        for (int i = 0; i < mpool_size; i++) {
            int j = std::uniform_int_distribution<int>(i, mpool_size - 1)(rg);
            std::swap(mpool[i], mpool[j]);
        }
    }
};



#endif //GA_TOURNAMENT_SELECTION_H
