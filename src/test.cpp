#include "aquila/aquila.h"
#include "Synthesis/AudioDevice.h"
#include "Synthesis/TextureSynthesizer.h"

#include <iostream>
#include <cmath>
#include <thread>
#include <chrono>
#include <fenv.h>

using namespace TextureSynthesis;

static bool stillPlaying = true;

void donePlayingCallback()
{
    stillPlaying = false;
}

int main(int argc, char **argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    bool makeWaveFile = false;
    std::string outfile;
    if(argc < 2)
    {
        std::cerr << "Usage: TextureAnimation <input file> [<output file>]\n";
        return 1;
    }
    if(argc > 2)
    {
        makeWaveFile = true;
        outfile = argv[2];
    }

    Aquila::WaveFile input(argv[1]);
    int sourceLen = pow(2, (int)(log(input.getSamplesCount()) / log(2)));
    sourceLen = std::min(sourceLen, (int)pow(2, 12));
    double sampleRate = input.getSampleFrequency();
    
    int numChannels = input.getChannelsNum();
    double channelAvg;

    Signal sourceSignal(sourceLen, sampleRate);
    for(int i = 0; i < sourceLen; i++)
    {
        channelAvg = 0;
        for(int j = 0; j < numChannels; j++)
            channelAvg += input.sample(i);
        sourceSignal._signal[i] = channelAvg / numChannels;
    }

    TextureSynthesizer synthesizer(sourceSignal);
    Signal outSignal(sourceLen, sampleRate);
    synthesizer.synthesize(outSignal);

    if(makeWaveFile)
    {
        Aquila::SignalSource output(outSignal.realPart(), sampleRate);
        Aquila::WaveFile::save(output, outfile);
    }
    else
    {
        AudioDevice audioDevice;
        audioDevice.play(outSignal, donePlayingCallback);

        while(stillPlaying)
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    return 0;
}
