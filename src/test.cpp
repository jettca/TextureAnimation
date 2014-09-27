#include "aquila/aquila.h"
#include "FilterBank.h"
#include "LowpassFilter.h"
#include "CochlearFilter.h"
#include "AudioDevice.h"

using namespace TextureSynthesis;

static bool _donePlaying = false;

void donePlayingCallback()
{
    _donePlaying = true;
}

int main()
{
    const std::size_t SIZE = 64;//pow(2, 16);
    const Aquila::FrequencyType sampleRate = 2000, f1 = 96, f2 = 813;

    Aquila::SineGenerator sine1(sampleRate);
    sine1.setAmplitude(32).setFrequency(f1).generate(SIZE);
    Aquila::SineGenerator sine2(sampleRate);
    sine2.setAmplitude(8).setFrequency(f2).setPhase(0.75).generate(SIZE);
    Aquila::SignalSource sum = sine1 + sine2;

    CochlearFilter *filter1P = new CochlearFilter(96);
    LowpassFilter *filter2P = new LowpassFilter(500);

    FilterBank filterBank;
    filterBank.addFilter(std::shared_ptr<Filter>((Filter*)filter1P));
    filterBank.addFilter(std::shared_ptr<Filter>((Filter*)filter2P));

    std::vector<Signal> signals = filterBank.apply(sum);

    return 0;
} 
