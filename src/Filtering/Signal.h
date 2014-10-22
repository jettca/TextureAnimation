#pragma once

#include <vector>
#include "aquila/aquila.h"

namespace TextureSynthesis
{
    class Signal
    {
    public:
        Signal(int length, double sampleRate);
        Signal(const Aquila::SignalSource& source);
        void makeAnalytic();
        void makeEnvelope();
        void makeEnvelope(Signal& phase);
        void pow(double a);
        std::vector<double> realPart() const;
        std::vector<double> imaginaryPart() const;

        void set(const Signal& signal);

        Aquila::SignalType _signal;
        double _sampleRate;
    };
}
