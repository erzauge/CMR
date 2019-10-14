#include "UnionFind.hpp"


UnionFind::UnionFind(long N_n=50){
        N=N_n;
        node.resize(N, -1);
}

void UnionFind::add(long x){
        if (N<=x) {
                N=x+1;
                node.resize(N, -1);
        }
}

long UnionFind::Find(long x){
        if(x<N) {
                long y=x;
                while (node[y]>=0) { // root search
                        y=node[y];
                }
                while (node[x]>=0) { // path compression
                        long z = node[x];
                        node[x] = y;
                        x = z;
                }
                return y;
        }
        else{
                add(x);
                return x;
        }
}

void UnionFind::Union(long a, long b){
        long r_a=Find(a);
        long r_b=Find(b);

        if (r_a!=r_b) { // join if not already joind
                if (node[r_a]>node[r_b]) { //append the smale to the biger one
                        node[r_a]+=node[r_b];
                        node[r_b]=r_a;
                }
                else{
                        node[r_b]+=node[r_a];
                        node[r_a]=r_b;
                }
        }
}

std::vector<long> UnionFind::ClusterSize(){
        std::vector<long> r;
        for (long i : node) {
                if (i<0) {
                        r.push_back(-i);
                }
        }
        return r;
}

long UnionFind::Label(long *label){
        std::vector<long> helper;
        helper.resize(node.size(),0);
        long count=0;
        for (unsigned long i = 0; i < node.size(); i++) {
                long x=Find(i);
                if (helper[x]!=0) {
                        label[i]=helper[x];
                }
                else{
                        count++;
                        helper[x]=count;
                        label[i]=count;
                }
        }
        return count;
}

long UnionFind::size(){
        return node.size();
}

bool UnionFind::Singleton(long x){
  return node[x]==-1;
}
