#include "gtest/gtest.h"
#include "../src/Histogram.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdlib>

TEST(Histogram, SetRange_simple){
  Histogram h;
  std::vector<double> r={0.5,1.5,2.5,3.5};
  h.SetRange(r);
  std::vector<double> r2;
  std::vector<double> p2;
  h.pdf(r2, p2);
  EXPECT_EQ(r, r2);
}

TEST(Histogram, SetRange_unorderd){
  Histogram h;
  std::vector<double> r={3.5,0.5,2.5,1.5};
  h.SetRange(r);
  std::vector<double> r2;
  std::vector<double> p2;
  h.pdf(r2, p2);
  std::sort(r.begin(), r.end());
  EXPECT_EQ(r, r2);
}

TEST(Histogram, SetUniformRange_simple){
  Histogram h;
  std::vector<double> r={0.5,1.5,2.5,3.5};
  h.SetUniformRange(0.5, 3.5, 3);
  std::vector<double> r2;
  std::vector<double> p2;
  h.pdf(r2, p2);
  EXPECT_EQ(r, r2);
}

TEST(Histogram, SetUniformRange_wrong_input){
  Histogram h;
  EXPECT_ANY_THROW(h.SetUniformRange(3.5, 0.5, 3));
  EXPECT_ANY_THROW(h.SetUniformRange(0.5, 0.5, 3));
}

TEST(Histogram, SetLogRange_simple){
  Histogram h;
  h.SetLogRange(1.8, 1800, 3);
  std::vector<double> r={1.8,18,180,1800};
  std::vector<double> r2;
  std::vector<double> p2;
  h.pdf(r2, p2);
  for (size_t i = 0; i < r.size(); i++) {
    EXPECT_DOUBLE_EQ(r[i], r2[i]);
  }

}

TEST(Histogram, SetLogRange_wrong_input){
  Histogram h;
  EXPECT_ANY_THROW(h.SetLogRange(3.5, 0.5, 3))<<"wrong order";
  EXPECT_ANY_THROW(h.SetLogRange(0.5, 0.5, 3))<<"equal borders";
  EXPECT_ANY_THROW(h.SetLogRange(-0.5, 0.5, 3))<<"negativ border";
  EXPECT_ANY_THROW(h.SetLogRange(0, 0.5, 3))<<"0 as border";
}

TEST(Histogram, SetAdptiveLogRange_simple){
  Histogram h;
  h.SetAdptiveLogRange(1, 1000, 9, 3);
  std::vector<double> r={1,10,100,1000};
  std::vector<double> r2;
  std::vector<double> p2;
  h.pdf(r2, p2);
  for (size_t i = 0; i < r.size(); i++) {
    EXPECT_DOUBLE_EQ(r[i], r2[i]);
  }
}

TEST(Histogram, SetAdptiveLogRange_steps){
  Histogram h;
  h.SetAdptiveLogRange(1, 8900, 5, 3);

  std::vector<double> r2;
  std::vector<double> p2;
  h.pdf(r2, p2);

  for (size_t i = 0; i < r2.size()-1; i++) {
    EXPECT_EQ(std::fmod((r2[i+1]-r2[i]),5), 0);
  }
}

TEST(Histogram, SetAdptiveLogRange_wrong_input){
  Histogram h;
  EXPECT_ANY_THROW(h.SetAdptiveLogRange(3.5, 0.5, 1, 3))<<"wrong order";
  EXPECT_ANY_THROW(h.SetAdptiveLogRange(0.5, 0.5, 1, 3))<<"equal borders";
  EXPECT_ANY_THROW(h.SetAdptiveLogRange(-0.5, 0.5, 1, 3))<<"negativ border";
  EXPECT_ANY_THROW(h.SetAdptiveLogRange(0, 8.5, 1, 3))<<"0 as border";
  EXPECT_ANY_THROW(h.SetAdptiveLogRange(1, 10, 2, 5))<<"not enough spase";
}

TEST(Histogram,Accumulate_simple){
  Histogram h;
  h.SetUniformRange(0.5, 10.5, 10);
  for (size_t i = 1; i < 11; i++) {
    h(i);
  }
  std::vector<double> r = h.GetBin();
  for (size_t i = 0; i < r.size(); i++) {
    EXPECT_EQ(r[i],1);
  }
}


TEST(Histogram,Accumulate_rigth_bin){
  Histogram h;
  h.SetUniformRange(0, 10, 10);
  h(1);
  std::vector<double> r = h.GetBin();
  EXPECT_EQ(r[1],1);
  EXPECT_EQ(r[0],0);
}

TEST(Histogram,Accumulate_not_in_range){
  Histogram h;
  h.SetUniformRange(0.5, 10.5, 10);
  EXPECT_ANY_THROW(h(0))<<"lower bound";
  EXPECT_ANY_THROW(h(11))<<"upper bound";
}

TEST(Histogram, PDF_nomelised_uniform){
  Histogram h;
  h.SetUniformRange(0.5, 10.5, 10);
  for (size_t i = 0; i < 10000; i++) {
    h((rand()%10)+1);
  }
  std::vector<double> r2;
  std::vector<double> p2;
  h.pdf(r2, p2);
  double sum=0;
  for (size_t i = 0; i < p2.size(); i++) {
    sum+=(r2[i+1]-r2[i])*p2[i];
  }
  EXPECT_EQ(sum, 1);
  // std::cout<<sum<<std::endl;
}
TEST(Histogram, PDF_nomelised_log){
  Histogram h;
  h.SetLogRange(1, 1000, 3);
  for (size_t i = 0; i < 100000; i++) {
    h((rand()%999)+1);
  }
  std::vector<double> r2;
  std::vector<double> p2;
  h.pdf(r2, p2);
  double sum=0;
  for (size_t i = 0; i < p2.size(); i++) {
    sum+=(r2[i+1]-r2[i])*p2[i];
  }
  EXPECT_EQ(sum, 1);
}
