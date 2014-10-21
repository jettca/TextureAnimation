#pragma once

#include "Filtering/Signal.h"

namespace TextureSynthesis
{
    class Downsampler
    {
    public:
        Downsampler(Signal& signal, double targetRate);
        void revert(Signal& signal);
        double actualSampleRate();

        static std::pair<int, double> newSizeAndRate(int oldSize, double oldRate,
                double targetRate);

    private:
        Signal _fineStructure;
        double _newSampleRate;
    };
}
