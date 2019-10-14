#include "cpu_cmr.hpp"
#include <iostream>

#include <tclap/CmdLine.h>

namespace TCLAP {
template <> struct ArgTraits<block> {
  typedef ValueLike ValueCategory;
};
} // namespace TCLAP

std::string dbl2str(double d)
{
    size_t len = std::snprintf(0, 0, "%.1f", d);
    std::string s(len+1, 0);
    // technically non-portable, see below
    std::snprintf(&s[0], len+1, "%.1f", d);
    // remove nul terminator
    s.pop_back();
    return s;
}

int main(int argc, char const *argv[]) {
  TCLAP::CmdLine cmd("", ' ', "", true);
  TCLAP::ValueArg<double> PArg("p", "p", "fraction antiferomagentic copling ",
                               false, 0, "double");
  cmd.add(PArg);
  TCLAP::ValueArg<double> TArg("T", "T", "invers Temperaturs", false, 0.2, "double");
  cmd.add(TArg);

  TCLAP::ValueArg<double> dTArg("", "dT", "invers Temperaturs range", false, 0.1, "double");
  cmd.add(dTArg);

  TCLAP::ValueArg<int> LArg("L", "L", "system size", false, 20, "int");
  cmd.add(LArg);

  TCLAP::ValueArg<unsigned long> seedArg("s", "seed", "seed", false,
                                  12345, "string");
  cmd.add(seedArg);
  TCLAP::ValueArg<long> numArg("r", "", "numper of steps", false,
                                  1000, "string");
  cmd.add(numArg);

  TCLAP::SwitchArg orderdSwitch("","orderd","intitialies orderd", cmd, false);

  cmd.parse(argc, argv);

  cpu_cmr sys_cmr(LArg.getValue(), TArg.getValue(), TArg.getValue()+.1, PArg.getValue());
  sys_cmr.set_seed(seedArg.getValue());
  sys_cmr.init_J();
  if (orderdSwitch.isSet()) {
    sys_cmr.init_order();
  }
  else{
    sys_cmr.init_rand();
  }
  vector<float> result;
  for (size_t i = 0; i < (unsigned)numArg.getValue(); i++) {
    sys_cmr.sweep();
    result=sys_cmr.measure();
    for (size_t j = 0; j < result.size(); j++) {
      std::cout << result[j]<<"\t";
    }
    std::cout << std::endl;
    if(i%5==0){
      sys_cmr.tempering();
    }
  }
  // sys_cmr.save_sys("test_crm");
  exit(0);


  return 0;
}
