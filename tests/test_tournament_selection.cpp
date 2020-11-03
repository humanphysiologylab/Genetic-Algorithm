#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "tournament_selection.h"
#include "pcg_random.hpp"



TEST(tournament_tests, nullInput)
{
    TournamentSelection tournament_selection{pcg64()};
    std::vector<std::pair<double, int>> sd_n_index;
    tournament_selection(0, 0, sd_n_index, 0);
}

TEST(tournament_tests, pool_is_complete)
{
    TournamentSelection tournament_selection{pcg64()};
    const int pool_size = 10;
    int mpool[pool_size];
    std::vector<std::pair<double, int>> sd_n_index(pool_size);
    //initialize
    for (int i = 0; i < pool_size; i++) {
        mpool[i] = -1;
        sd_n_index[i].first = i;
        sd_n_index[i].second = i;
    }
    tournament_selection(mpool, pool_size, sd_n_index, 0);
    for (int i = 0; i < pool_size; i++) {
        EXPECT_TRUE(mpool[i] >= 0);
        EXPECT_TRUE(mpool[i] < pool_size);
    }
}

TEST(tournament_tests, pool_is_random)
{
    TournamentSelection tournament_selection{pcg64()};
    const int pool_size = 50;
    int mpool[pool_size];
    std::vector<std::pair<double, int>> sd_n_index(pool_size);
    //initialize
    for (int i = 0; i < pool_size; i++) {
        sd_n_index[i].first = i;
        sd_n_index[i].second = i;
    }
    tournament_selection(mpool, pool_size, sd_n_index, 0);
    
    int mpool2[pool_size];
    tournament_selection(mpool2, pool_size, sd_n_index, 0);
    
    bool equal = true;
    for (int i = 0; i < pool_size; i++) {
        if (mpool[i] != mpool2[i]) {
            equal = false;
            break;
        }
    }
    EXPECT_FALSE(equal);
}
