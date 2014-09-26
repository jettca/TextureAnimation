#pragma once

#include <vector>
#include "aquila/aquila.h"

namespace TextureSynthesis
{
    class Signal
    {
    public:
        Signal(int length, double sampleRate)
            : _samples(length), _sampleRate(sampleRate)
        {}

        std::vector<Aquila::SampleType> _samples;
        double _sampleRate;
    };
}
