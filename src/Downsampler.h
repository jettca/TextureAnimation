#pragma once

#include "Signal.h"

namespace TextureSynthesis
{
    class Downsampler
    {
    public:
        Downsampler(Signal& signal, double newSampleRate);
        void revert(Signal& signal);
        double actualSampleRate();
    private:
        Signal _fineStructure;
        double _newSampleRate;
    };
}
