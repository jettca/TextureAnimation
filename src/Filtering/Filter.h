#pragma once

#include "aquila/aquila.h"
#include "Filtering/Signal.h"

namespace TextureSynthesis
{
    class Filter
    {
    public:
        virtual void filter(Aquila::SpectrumType& spectrum,
                double sampleRate) = 0;
        
        void filter(const Signal& signal, Aquila::SpectrumType& spectrum) const;
        void filter(Signal& signal) const;
    };
}
