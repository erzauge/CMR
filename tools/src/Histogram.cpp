#include "Histogram.hpp"
#include <algorithm>
#include <iterator>
#include <limits>
#include <cmath>
#include <iostream>

Histogram::Histogram():bins(1){
  range.push_back(-std::numeric_limits<float>::infinity());
  range.push_back(std::numeric_limits<float>::infinity());
  bin.push_back(0);
}

Histogram::Histogram(std::vector<double> _ragne){
  SetRange(_ragne);
}

void Histogram::operator()(double x){
  auto upper=std::upper_bound(range.begin(),range.end(),x);
  unsigned int bin_pos=std::distance(range.begin(), upper);
	if(bin_pos==0||bin_pos==bins+1) throw("out of range");
	bin[bin_pos-1]++;
  normalise++;
}

void Histogram::SetRange(std::vector<double> _range){
  if (_range.size()<2) throw("size of range to smale");
  range = _range;
  std::sort(range.begin(), range.end());
  bins  = range.size()-1;
  bin.clear();
  bin.resize(bins, 0);
}

void Histogram::SetUniformRange(double min, double max, int _bins){
  if (min>=max) throw("swape min and max");
  range.clear();
  bins=_bins;
  double step=(max-min)/(bins);
  for (size_t i = 0; i < bins+1; i++) {
    range.push_back(min+i*step);
  }
  bin.clear();
  bin.resize(bins, 0);
}

void Histogram::SetLogRange(double min, double max, int _bins){
  if (min>=max) throw("swape min and max");
  if (min <= 0) throw("smale or equel then 0");
  range.clear();
  bins=_bins;
  double ratio=std::pow((max/min),1./bins);
  for (size_t i = 0; i < bins+1; i++) {
    range.push_back(min*std::pow(ratio,i));
  }
  bin.clear();
  bin.resize(bins, 0);
}

void Histogram::SetAdptiveLogRange(double min, double max, double step_size, int _bins){
  if (min >= max) throw("swape min and max");
  if (min <= 0) throw("log not defind for valuse smaler 0");
  if ((max-min)/step_size<_bins) throw("not enough spase");
  range.clear();
  bins=_bins;
  int n=bins;
  range.push_back(min);
  double ratio=std::pow((max/min),1./bins);
  for (size_t i = 1; i < bins+1; i++) {
    double d = range.back()*ratio-range.back();
    d = ceil(d/step_size)*step_size;
    range.push_back(range.back()+d);
    n--;
    ratio=std::pow((max/range.back()),1./n);
  }
  bin.clear();
  bin.resize(bins, 0);
}

std::ostream& operator<<(std::ostream& os, const Histogram& H)
{
    for (size_t i = 0; i < H.bins; i++) {
      double d=H.range[i+1]-H.range[i];
      double p=H.bin[i]/(H.normalise*d);
      os << H.range[i] << "\t" << H.range[i+1] << "\t" << p << std::sqrt(p*d*(1-p*d)/(H.normalise))/d << "\n";
    }
    return os;
}

void Histogram::pdf(std::vector<double> &_range,std::vector<double> &_bin){
  _range.clear();
  _bin.clear();
  _bin.resize(bins);
  _range = range;
  double sum=0;
  for (size_t i = 0; i < bin.size(); i++) {
    sum+=bin[i];
  }
  // if(sum == 0) sum=1;
  for (size_t i = 0; i < bin.size(); i++) {
    _bin[i]=(bin[i])/((range[i+1]-range[i])*sum);
  }
}

std::vector<double> Histogram::GetBin(){
  return bin;
}
