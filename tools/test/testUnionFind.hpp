#include "gtest/gtest.h"
#include "../src/UnionFind.hpp"
#include <iostream>

//todo labeling
class TestUnionFind: public testing::TestWithParam<long> {
public:
};

TEST_P(TestUnionFind,construckt){
    UnionFind uf1(GetParam());
    ASSERT_EQ(uf1.size(), GetParam());
}

TEST_P(TestUnionFind,add){
    UnionFind uf1(20);
    uf1.add(GetParam());
    if(GetParam()<20)
    ASSERT_EQ(uf1.size(), 20);
    else
    ASSERT_EQ(uf1.size(), GetParam()+1);
}

TEST_P(TestUnionFind,find_myself){
    UnionFind uf1(20);
    ASSERT_EQ(uf1.Find(GetParam()), GetParam());
}

TEST_P(TestUnionFind,union_and_find){
    UnionFind uf1(20);
    uf1.Union(5, GetParam());
    ASSERT_EQ(uf1.Find(5), uf1.Find(GetParam()));
    uf1.Union(18, GetParam());
    ASSERT_EQ(uf1.Find(18), uf1.Find(GetParam()));
    ASSERT_EQ(uf1.Find(18), uf1.Find(5));
}

TEST_P(TestUnionFind,ClusterSize){
    UnionFind uf1(20);
    uf1.Union(0, GetParam());
    uf1.Union(18, GetParam());
    std::vector<long> c=uf1.ClusterSize();
    std::sort(c.begin(), c.end());
    ASSERT_EQ(*(c.end()-1),3 );
}

INSTANTIATE_TEST_CASE_P(UnionFind,TestUnionFind,testing::Values(1,12,21,4864));
