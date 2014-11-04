#pragma once

#include "Filtering/Signal.h"

namespace TextureSynthesis
{
    class Downsampler
    {
    public:

        /* Downsample signal while maintaining power-of-two length
         */
        Downsampler(Signal& signal, double targetRate);

        /* Revert signal to pre-downsampled state
         */
        void revert(Signal& signal) const;

        /* Returns the actual downsampled rate to maintain the signal's length
         */
        double actualSampleRate() const;

        /* Returns the actual size and rate a signal would be if it was attempted
         * to be downsampled
         */
        static std::pair<int, double> newSizeAndRate(int oldSize, double oldRate,
                double targetRate);

    private:
        /* Fine structure used for reconstructing signals
         */
        Signal _fineStructure;

        /* Downsampled rate
         */
        double _newSampleRate;
    };
}
