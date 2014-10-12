#pragma once

#include "Signal.h"
#include "FilterBank.h"
#include "Downsampler.h"

namespace TextureSynthesis
{
    // TODO: Make this class not static!
    class TextureFilterer
    {
    public:
        static void auditoryFilter(const Signal& signal,
                std::vector<Signal>& cochlearSignals,
                std::vector<std::vector<Signal>>& modulationSignals);

        static void recombine(std::vector<Signal>& cochlearEnvelopes,
                Signal& combinedSignal);
    private:
        static void generateCochlearBank(FilterBank& cochlearBank);
        static void generateModulationBank(FilterBank& modulationBank);
        static double erbsInverse(int numFilters);

        static std::vector<Signal> _someModulationSignals;
        static std::vector<Downsampler> _downsamplers;
        static std::vector<Signal> _subbandPhases;
        static double _cochlearExponent;
    };
}
