#pragma once

#include "Filtering/Signal.h"
#include "Filtering/FilterBank.h"
#include "Filtering/Downsampler.h"

namespace TextureSynthesis
{
    class TextureFilterer
    {
    public:
        TextureFilterer();

        void auditoryFilter(const Signal& signal,
                std::vector<Signal>& cochlearEnvelopes,
                std::vector<std::vector<Signal>>& modulationSignals,
                bool makeRecombinable); 

        void modulationFilter(const std::vector<Signal>& cochlearEnvelopes,
                std::vector<std::vector<Signal>>& modulationSignals);

        void recombine(std::vector<Signal>& cochlearEnvelopes,
                Signal& combinedSignal);

        static double targetDownsampleRate;
        static int numCochlearEnvelopes;
        static int numModulationSignals;
    private:
        void generateCochlearBank(FilterBank& cochlearBank);
        void generateModulationBank(FilterBank& modulationBank);
        double erbsInverse(int numFilters);

        std::vector<Downsampler> _downsamplers;
        std::vector<Signal> _subbandPhases;
        double _cochlearExponent;
    };
}
