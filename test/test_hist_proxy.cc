#include "XSecAna/HistProxyish.h"
#include <iostream>
#include "test_utils.h"

using namespace xsec;

class HistWrap {
public:
    HistWrap() = default;
    HistWrap(std::vector<HistWrap*> & registry)
    : fHist(10, 0, 10){
        registry.push_back(this);
    }
    void SetHist(Hist<double, 10> h) { fHist = h; }
    Hist<double, 10> GetHist() { return fHist; }
private:
    Hist<double, 10> fHist;
};

// TODO this might be a problem. Would like to do partial instantiation
template<>
Hist<double, 10> HistProxyish<HistWrap, double, 10>::GetHist() {
    return fWrapped.GetHist();
}

int main(int argc, char ** argv)
{
    bool verbose = false;
    if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
    bool pass = true;
    bool test;

    std::vector<HistWrap*> registry;
    HistProxyish<HistWrap, double, 10> wrapper(registry);

    // wrapped histogram starts off as empty
    TEST_ARRAY("starts empty",
               (wrapper.GetHist().Contents()),
               (Eigen::Array<double, 1, 10>::Zero(10)),
               0);

    // externally fill the wrapped hist
    registry[0]->SetHist(test::utils::get_simple_data<double, 10>());

    TEST_ARRAY("filled",
               (wrapper.GetHist().Contents()),
               (test::utils::get_simple_data<double, 10>().Contents()),
               0)

    bool caught_exception = false;
    try {
        wrapper.Normalize("area");
    }
    catch (xsec::exceptions::HistProxyishError & e) {
        caught_exception = true;
    }
    pass &= caught_exception;

    // TODO test that HistProxyish objects will work with analysis components
    return !pass;
}