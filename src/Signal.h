#pragma once

#include <vector>
#include "aquila/aquila.h"

namespace TextureSynthesis
{
    class Signal
    {
    public:
        Signal(int length, double sampleRate);
        Signal(Aquila::SignalSource source);
        void makeAnalytic();
        void makeEnvelope();
        void pow(double a);
        std::vector<double> realPart();
        std::vector<double> imaginaryPart();

        Aquila::SignalType _signal;
        double _sampleRate;
    };
}
