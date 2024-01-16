#pragma once
#include <vector>
#include <string>

namespace caf
{
    class SRSystParamHeader
    {
    public:
        SRSystParamHeader();
        ~SRSystParamHeader();

        int nshifts;
        std::string name;
        int id;
    };

    class SRWeightGlobal
    {

    public:
    SRWeightGlobal();
    ~SRWeightGlobal();

    std::vector<SRSystParamHeader> params;
    };

    class SRGlobal
    {
    public:
    SRGlobal();
    ~SRGlobal();

    SRWeightGlobal wgts;
    };


  

} // end namespace
