#include "cpu_cmr.hpp"
#include "UnionFind.hpp"
#include "bin_io.hpp"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <bitset>
#include <sstream>


/**
 * @brief bit spin to double
 * @details locks if bit is set and returns depending on it -1 or 1
 *
 * @param s spin
 * @param j selected bit
 *
 * @return -1.0 if bit j in s is set else 1.0
 */
inline double bit_to_double(spin_t s, int j){
	spin_t a=(s<<j&((spin_t)1ULL<<63))|0x3FF0000000000000ULL; //0x3FF0000000000000 in double = 1. If bit j in s is set sets sing bit in double
	void * b= &a; //casting to void pinter
	return *((double*) b); //
}

/**
 * @brief swape bits
 * @details swape bit i and j
 *
 * @param i, j selected bit positon to swap
 * @param b input to swap
 *
 * @return swaped bits
 */
inline spin_t SwapBit(unsigned int i, unsigned int j, spin_t b){
unsigned int x = ((b >> i) ^ (b >> j)) & ((1U << 1) - 1); // XOR temporary
return b ^ ((x << i) | (x << j));
}


/**
 * @brief converts xyz codinets to id

 *
 * @param x, y, z cordinets
 * @param L lienear size
 * @param L2 squer of lienear size
 *
 * @return id
 */
inline int xyz2id(int x, int y, int z, int L, int L2) {
  return ((x + L) % L) + L * ((y + L) % L) + L2 * ((z + L) % L);
}


cpu_cmr::cpu_cmr(int L_, double B1_, double B2_, double P_)
    : L(L_), P(P_), B1(B1_), B2(B2_), uni_dist(0.0, 1.0) {
  N = L * L * L;
  L2 = L * L;
	double dB=(B2-B1)/31.;
	for (size_t i = 0; i < 32; i++) {
		P_blue[i] = 1. - exp((-4.) * (B1+dB*i));
		P_red[i] = 1. - exp((-2.) * (B1+dB*i));
	}
  gen.seed(seed);
  s = new spin_t[N];
  J = new spin_t[N * 6];
  nn = new int[N * 6];
  label = new long[N];
  for (long x = 0; x < L; x++) {
    for (long y = 0; y < L; y++) {
      for (long z = 0; z < L; z++) {
        nn[6 * xyz2id(x, y, z, L, L2) + 0] = xyz2id(x - 1, y, z, L, L2);
        nn[6 * xyz2id(x, y, z, L, L2) + 1] = xyz2id(x + 1, y, z, L, L2);
        nn[6 * xyz2id(x, y, z, L, L2) + 2] = xyz2id(x, y - 1, z, L, L2);
        nn[6 * xyz2id(x, y, z, L, L2) + 3] = xyz2id(x, y + 1, z, L, L2);
        nn[6 * xyz2id(x, y, z, L, L2) + 4] = xyz2id(x, y, z - 1, L, L2);
        nn[6 * xyz2id(x, y, z, L, L2) + 5] = xyz2id(x, y, z + 1, L, L2);
      }
    }
  }
	if_measured=false;
}

void cpu_cmr::set_seed(long se) {
  seed =se;
  gen.seed(se);
  uni_dist.reset();
  if (file_is_open) {
    file.close();
    file_is_open=false;
  }
}
/**
 * @brief initialise spins randomly
 */
void cpu_cmr::init_rand() {
  for (long i = 0; i < N; i++) {
    s[i]=gen();
  }
}
/**
 * @brief initialise spins to -1
 */
void cpu_cmr::init_order() {
  for (long i = 0; i < N; i++) {
    s[i] = ~(spin_t)0ULL;
  }
}
/**
 * @brief initialise J
 */
void cpu_cmr::init_J() {
  for (long i = 0; i < N; i++) {
    for (long j = 0; j < 3; j++) {
      if (uni_dist(gen) < P) {
        J[6 * i + 2 * j] = (spin_t)0ULL;
        J[6 * nn[6 * i + 2 * j] + 2 * j + 1] = (spin_t)0ULL;

      } else {
        J[6 * i + 2 * j] = ~(spin_t)0ULL;
        J[6 * nn[6 * i + 2 * j] + 2 * j + 1] = ~(spin_t)0ULL;
      }
    }
  }
}

void cpu_cmr::set_T(double B1_) {
  B1=B1_;
	double dB=(B2-B1)/31.;
	for (size_t i = 0; i < 32; i++) {
		P_blue[i] = 1 - exp((-4.) * (B1+dB*i));
		P_red[i] 	= 1 - exp((-2.) * (B1+dB*i));
	}
 }


void cpu_cmr::set_T2(double B1_,double B2_) {
  B1=B1_;
	B1=B2_;
 	double dB=(B2-B1)/31.;
 	for (size_t i = 0; i < 32; i++) {
 		P_blue[i] = 1 - exp((-4.) * (B1+dB*i));
 		P_red[i] 	= 1 - exp((-2.) * (B1+dB*i));
 	}
}

long cpu_cmr::get_N() { return N; }

double cpu_cmr::get_T() { return B1; }

void cpu_cmr::save_J(string fname) {
  std::ofstream file(fname);
  if (!file.is_open()) {
    // error
    exit(-1);
  }
  for (long i = 0; i < N; i++) {
    for (long j = 0; j < 6; j++) {
      file << J[6 * i + j] << " ";
    }
    file << "\n";
  }
  file.flush();
  file.close();
}

void cpu_cmr::load_J(string fname) {
  std::ifstream file(fname);
  if (!file.is_open()) {
    // error
    exit(-1);
  }
  for (long i = 0; i < N; i++) {
    for (long j = 0; j < 6; j++) {
      file >> J[6 * i + j];
    }
  }
  file.close();
}
vector<float> cpu_cmr::measure() {
	if(!if_measured){
		spin_t mask1=0x5555555555555555ULL;
		vector<float> result;
		double sum[64+32]{};
		for (long i = 0; i < N; i++) {
			for (long j = 0; j < 3; j++) {
				spin_t e =J[6 * i + 2 * j] ^s[i]^s[nn[6 * i + 2 * j]];
				for (size_t l = 0; l < 64; l++) {
					sum[l]-=bit_to_double(e,63-l);
				}
			}
			spin_t q = (s[i]^(s[i]>>1))&mask1;
			// std::cerr << (std::bitset<64>) q << '\n';
			for (size_t j = 0; j < 32; j++) {
				sum[64+j]-=bit_to_double(q,2*(31-j)+1);
				// if(bit_to_double(q,2*j)!=1)
				// std::cerr << bit_to_double(q,2*j) << '\n';
			}
		}
		for (size_t i = 0; i < 64+32; i++) {
			result.push_back(sum[i]/(float)N);
		}
		if_measured=true;
		measure_r=result;
	}
  return measure_r;
}

void cpu_cmr::sweep() {
	if_measured=false;
  UnionFind uf_blue[32]{N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N};
  UnionFind uf_gray[32]{N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N};
  spin_t mask1=0x5555555555555555ULL;
  spin_t mask2=0xAAAAAAAAAAAAAAAAULL;
  spin_t * rand=new spin_t[N];
  for (long i = 0; i < N; i++) {
    rand[i]=gen();
		// if (i==1) {
		// 	std::cerr << __builtin_popcountll(rand[i]&mask1) <<'\n';
		// }
    for (long j = 0; j < 6; j += 2) {
      spin_t a =J[6 * i + j]^s[nn[6 * i + j]]^s[i];
      spin_t blue = a & (a>>1);
      spin_t red 	= a ^ (a>>1);
			spin_t pb	=	0ULL;
			spin_t pr =	0ULL;
			for (size_t i = 0; i < 32; i++) {
				pb |= uni_dist(gen) < P_blue[i]?0x1ULL<<(i*2):0ULL;
				pr |= uni_dist(gen) < P_red[i]?0x2ULL<<(i*2):0ULL;
			}
      blue	= blue&pb;
      red		= red&pr;
      for (size_t k = 0; k < 32; k++) {
        if ((blue>>(2*k))&0x1ULL) {
					// if (k==1) {
					//   std::cerr << i << "\t" <<  nn[6 * i + j] <<"\t"<<L<< '\n';
					// }
          uf_blue[k].Union(nn[6 * i + j], i);
          uf_gray[k].Union(nn[6 * i + j], i);
        }
        if ((red >>(2*k))&0x1ULL) {
          uf_gray[k].Union(nn[6 * i + j], i);
        }
      }
    }
  }
  for (size_t i = 0; i <(unsigned)N; i++) {
		spin_t t_rand=0x0ULL;
		for (size_t j = 0; j < 32; j++) {
			if (!uf_gray[j].Singleton(i)) {
				long   c  = uf_gray[j].Find(i);
				t_rand   |= rand[c]&(0x3ULL<<(2*j));
				// if (i==1) {
				//   std::cerr << (std::bitset<64>) t_rand << '\n';
				// }
			}
		}

		spin_t r1 = (t_rand&mask1)|((t_rand&mask1)<<1);
		spin_t r2 = (t_rand&mask2)|((t_rand&mask2)>>1);
		spin_t pq = (s[i]^(s[i]>>1))&mask1;

		spin_t f  = pq^mask1^r2^s[i];
		// std::cerr << c << '\n';
		// if (i==1) {
		// 	std::cerr << "test" << '\n';
		//   std::cerr << (std::bitset<64>) (s[i]) << '\n';
		// 	std::cerr << (std::bitset<64>) (f&r1) << '\n';
		// }
		s[i]     ^= f&r1;
		// if (i==1) {
		//   std::cerr << (std::bitset<64>) (s[i]) << '\n';
		// }
  }
  delete[] rand;
}

istream &cpu_cmr::load(istream &stream) { return stream; }

ostream &cpu_cmr::save(ostream &stream) { return stream; }

void cpu_cmr::swap(cpu_sys *sys, std::unique_ptr<spin_t[]> mask) {
  cpu_cmr *sys_wolff = dynamic_cast<cpu_cmr *>(sys);
  if (sys_wolff != NULL) {
    std::cerr << "error wrong type "<<mask[0];
    // swap(sys_wolff, std::move(mask));
  } else {
		std::cerr << "error is not implemented "<< '\n';
    // LOG(LOG_ERROR) << "conversion error";
    exit(-2);
  }
}

// void cpu_cmr::swap(cpu_cmr *sys, std::unique_ptr<spin_t[]> p) {
//   if (p[0] != 0) {
//
//     // int *b = s;
//     // s = sys->s;
//     // sys->s = b;
//   }
// }
cpu_sys *cpu_cmr::clone() const { return new cpu_cmr(L, B1, P, B2); }

cpu_cmr::~cpu_cmr() {
  // std::cerr << save_cusrter_size.size() << '\n';
  if(data.is_open()){
    binary_write(data,(int)0,1);
  }
  data.close();

  delete[] s;
  delete[] J;
  delete[] nn;
  delete[] label;
}

void cpu_cmr::save_sys(string prefix) {
  std::cerr << prefix << '\n';
	for (int i = 0; i < 64; ++i) {
    stringstream str_i;
    str_i << i;
    // image setup
    ofstream file(prefix + str_i.str() + ".pbm");
    file << "P1" << endl;
    file << L << " " << L * L << endl;
    // print image
    for (int j = 0; j < L * L; ++j) {
      for (int k = 0; k < L; ++k) {
        file << (((s[j * L + k] & ((spin_t)1 << i)) == 0)
                     ? "0 "
                     : "1 ");  // muss ich noch mal überprüfen
      }
      file << endl;
    }
  }

  // exit(-3);
}

void cpu_cmr::tempering() {
	vector<float> result=measure();
	static std::uniform_int_distribution<int> tempPar(0, 30);
	static std::uniform_int_distribution<int> tempSce(0, 1);
	int i = tempPar(gen);
	int a = 2*i+tempSce(gen);
	int b = 2*(i+1)+tempSce(gen);
	double dB=((B2-B1)/31.);
	double dE=result[b]-result[a];
	std::cerr << exp(dE*dB) << '\n';
	if (uni_dist(gen) < exp(dE*dB)) {
		for (long i = 0; i < N; i++) {
			s[i]=SwapBit(a,b,s[i]);
		}
	}
}
