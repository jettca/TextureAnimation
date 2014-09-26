#pragma once

#include "Filter.h"

namespace TextureSynthesis
{
    class LowpassFilter : public Filter
    {
    public:
        LowpassFilter(double lowpassFrequency);
        using Filter::filter;
        void filter(Aquila::SpectrumType& spectrum, double sampleRate);

    private:
        double _lowpassFrequency;
    };
}
