#pragma once

#include "Filtering/Signal.h"
#include "Filtering/FilterBank.h"
#include "Filtering/Downsampler.h"

namespace TextureSynthesis
{
    class TextureFilterer
    {
    public:

        /* Constructs a new filterer
         */
        TextureFilterer();

        /* Takes an input signal and filters it into cochlear envelopes and
         * modulation signals.  If makeRecombinable is true, it will store
         * all lost data so that the cochlear envelopes can be recombined
         * into the original signal.
         */
        void auditoryFilter(const Signal& signal,
                std::vector<Signal>& cochlearEnvelopes,
                std::vector<std::vector<Signal>>& modulationSignals,
                bool makeRecombinable);

        /* Recombines the cochlear envelopes into the combined signal using data
         * from the last call to auditoryFilter for which makeRecombinable was
         * set to true.
         */
        void recombine(std::vector<Signal> cochlearEnvelopes,
                Signal& combinedSignal) const;

        /* Filters a set of cochlear envelopes into their modulation signals
         */
        void modulationFilter(const std::vector<Signal>& cochlearEnvelopes,
                std::vector<std::vector<Signal>>& modulationSignals) const;

        static void modulationFilter(const Signal& cochlearEnvelope,
                std::vector<Signal>& modulationSignals);

        // Various static getters
        static int numCochlearEnvelopes();
        static int numModulationSignals();
        static double targetDownsampleRate();
        static const FilterBank& modulationBank();

    private:

        // Data for recombining envelopes into signals
        std::vector<Downsampler> _downsamplers;
        std::vector<Signal> _subbandPhases;

        /* Generates the filter bank used to convert a signal into its cochlear envelopes.
         */
        static FilterBank generateCochlearBank();

        /* Generates the filter bank used to convert a cochlear envelope into its
         * modulation signals.
         */
        static FilterBank generateModulationBank();

        /* Computes the spacing between bandpass filter frequencies.
         */
        static double erbsInverse(int numFilters);

        // Static filter banks
        static FilterBank _cochlearBank;
        static FilterBank _modulationBank;

        // Fixed filtering parameters
        static int _numCochlearEnvelopes;
        static int _numModulationSignals;
        static double _targetDownsampleRate;
        static double _cochlearExponent;
        static double _oneOverCochlearExponent;
    };
}
