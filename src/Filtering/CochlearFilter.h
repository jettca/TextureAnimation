#pragma once

#include "Filtering/Filter.h"

namespace TextureSynthesis
{
    class CochlearFilter : public Filter
    {
    public:
        CochlearFilter(double centerFrequency);
        using Filter::filter;
        void filter(Aquila::SpectrumType& spectrum, double sampleRate) const;

    private:
        double _centerFrequency;
    };
}
