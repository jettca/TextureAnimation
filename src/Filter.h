#pragma once

#include "aquila/aquila.h"

namespace TextureSynthesis
{
    class Filter
    {
    public:
        virtual void filter(Aquila::SpectrumType& spectrum,
                double sampleRate) = 0;
        
        Aquila::SpectrumType filter(Aquila::SignalSource signal);

        virtual ~Filter() { }
    };
}
