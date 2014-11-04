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

        void set(const Signal& signal);

        std::complex<double>& operator[](const int index);
        const std::complex<double>& operator[](const int index) const;
        void resize(size_t size);
        size_t size() const;
        Aquila::SignalType& samples();
        const Aquila::SignalType& samples() const;


        void makeAnalytic();
        void makeEnvelope();
        void makeEnvelope(Signal& phase);

        void pow(double a);

        std::vector<double> realPart() const;
        std::vector<double> imaginaryPart() const;

        double sampleRate;

    private:
        Aquila::SignalType _signal;
    };
}
