#pragma once
#include <vector>
#include <iostream>

class Histogram {
private:
  unsigned int bins;
  std::vector<double> range;
  std::vector<double> bin;
  double normalise;
public:
  Histogram ();
  Histogram (std::vector<double> _range);
  void SetRange (std::vector<double> _range);
  void SetUniformRange (double min, double max,int _bins);
  void SetLogRange (double min, double max,int _bins);
  void SetAdptiveLogRange (double min, double max, double step_size,int _bins);
  void operator() (double x);
  void pdf(std::vector<double> &_range,std::vector<double> &_bin);
  std::vector<double> GetBin();
  friend std::ostream& operator<<(std::ostream& os, const Histogram& dt);
  // ~Histogram ();
};
