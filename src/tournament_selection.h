#ifndef GA_TOURNAMENT_SELECTION_H
#define GA_TOURNAMENT_SELECTION_H

#include <vector>
#include <forward_list>
#include <random>
#include <cassert>
#include <iostream>

template<typename InitializedRandomGenerator>
class TournamentSelection
{
    InitializedRandomGenerator rg;

    void tournament_basic_half(int *mpool, size_t mpool_size, const std::vector<std::pair<double, int>> & sd_n_index, int number_of_ignored_losers)
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
        size_t list_size = sd_n_index.size() - number_of_ignored_losers;
        std::forward_list<int> list_of_indices(list_size);
        
        /* First, copy indices from sd_n_index
         * to the list wrt the order
         */
        auto a = list_of_indices.begin();
        for (size_t i = 0; i < list_size; i++) {
            *a = sd_n_index[i].second;
            a++;
        }
        
        /* 'a' points at the current winner, since our list is sorted,
         * every organism after 'a' is inferior to it.
         */
        a = list_of_indices.begin();
        for (size_t i = 0; i < mpool_size; i++) {
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

    void tournament_basic_half_fast(int *mpool, size_t mpool_size, std::vector<std::pair<double, int>> & sd_n_index)
    {
        /* tournament_basic_half_fast fills in mpool with indices of organisms in mating pool according to the tournament selection without replacement.
         */

        //first, shuffle sd_n_index
        for (unsigned i = 0; i < sd_n_index.size(); i++) {
            int j = std::uniform_int_distribution<int>(i,  sd_n_index.size() - 1)(rg);
            std::swap(sd_n_index[i], sd_n_index[j]);
        }
        
        //now, for each following pair, find a winner and put it into mpool
        assert(mpool_size * 2 <= sd_n_index.size());
        assert(sd_n_index.size() % 2 == 0);
        for (size_t i = 0; i < mpool_size; i++) {
            const size_t one = 2 * i;
            const size_t two = 2 * i + 1;
            if (sd_n_index[one].first < sd_n_index[two].first)
                mpool[i] = sd_n_index[one].second;
            else
                mpool[i] = sd_n_index[two].second;
        }
    }

public:
    TournamentSelectionFast(InitializedRandomGenerator rg)
    : rg(rg)
    {}

    void operator()(int *mpool, size_t mpool_size, std::vector<std::pair<double, int>> & sd_n_index, int number_of_ignored_losers, double threshold)
    {
        /* FAST version!
         * sd_n_index is expected to be partially sorted by sd
         * in increasing order with a tail of worst number_of_ignored_losers organisms.
         * tournament selection fills in mpool with indices of organisms in mating pool.
         * each call of tournament_basic_half fills half of mating pool.
         * sd_n_index will be permutated!
         * any specie with error larger than threshold is filtered 
         */
        
        assert(mpool_size % 2 == 0);
        //trim losers
        assert(sd_n_index.size() - number_of_ignored_losers >= mpool_size);
        sd_n_index.resize(sd_n_index.size() - number_of_ignored_losers);

        //replace the rest of really bad losers by random fine species
        int first_bad_index;
        for (first_bad_index = sd_n_index.size(); first_bad_index > 0; first_bad_index--) {
            if (sd_n_index[first_bad_index - 1].first < threshold)
                break;
        }
        std::cout << "really bad filtered by tournament (sd > " << threshold << ")" << ": " << sd_n_index.size() - first_bad_index << std::endl;
        if (first_bad_index == 0)
            throw("Tournament: all species are bad!");
        for (int i = first_bad_index; i < sd_n_index.size(); i++) {
            int j = std::uniform_int_distribution<int>(0,  first_bad_index - 1)(rg);
            sd_n_index[i] = sd_n_index[j];
        }
        //probably sd_n_index is not sorted now
        //but anyway we gonna shuffle it
        tournament_basic_half_fast(mpool                 , mpool_size / 2, sd_n_index);
        tournament_basic_half_fast(mpool + mpool_size / 2, mpool_size / 2, sd_n_index);
    }
};

#endif //GA_TOURNAMENT_SELECTION_H
