#pragma once

#include <vector>
#include "aquila/aquila.h"

namespace TextureSynthesis
{
    class Signal
    {
    public:
        Signal(int length, double sampleRate)
            : _signal(length), _sampleRate(sampleRate)
        { }

        Aquila::SignalType _signal;
        double _sampleRate;
    };
}
