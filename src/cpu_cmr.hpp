#pragma once

#include "cpu_sys.hpp"
#include <random>
#include <fstream>
#include <string>

typedef unsigned long long spin_t;

class cpu_cmr : public cpu_sys {
private:
spin_t * s;
spin_t * J;
int * nn;
int L;
int L2;
long N;
double P;
double P_blue[32];
double P_red[32];
double B1;
double B2;
bool if_measured;
vector<float> measure_r;
long *label;
int num_label;
ofstream file;
ofstream data;
bool file_is_open=false;
unsigned long long seed=12345UL;
std::vector<long> clusters;
std::mt19937_64 gen;
std::uniform_real_distribution<double> uni_dist;
public:
cpu_cmr(int L_, double B1_, double B2_,double P_);
void set_seed(long se) override;
void init_rand() override;
void init_order() override;
void init_J() override;
void load_J(string fname) override;
void save_J(string fname) override;
void sweep() override;
void set_T(double ) override;
void set_T2(double ,double);
vector<float> measure() override;
void save_sys(string prefix) override;
long get_N() override;
double get_T() override;
void swap(cpu_sys* sys, std::unique_ptr<spin_t[]>) override;
// void swap(cpu_cmr* sys, std::unique_ptr<spin_t[]>);
ostream& save(ostream& stream) override;
istream& load(istream& stream) override;
cpu_sys * clone() const override;
void tempering();
~cpu_cmr();
};
