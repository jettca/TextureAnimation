#include "aquila/aquila.h"
#include "FilterBank.h"
#include "LowpassFilter.h"
#include "CochlearFilter.h"
#include "AudioDevice.h"

#include <complex>
#include <iostream>

using namespace TextureSynthesis;

int main()
{
    const std::size_t SIZE = pow(2, 16);
    const Aquila::FrequencyType sampleRate = 2000, f1 = 180, f2 = 20;

    Aquila::SineGenerator sine1(sampleRate);
    sine1.setAmplitude(16).setFrequency(f1).generate(SIZE);
    Aquila::SineGenerator sine2(sampleRate);
    sine2.setAmplitude(1).setFrequency(f2).generate(SIZE);
    Signal sum = Signal(sine1 * sine2);

    return 0;
}
