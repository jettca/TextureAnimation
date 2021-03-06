#include "Synthesis/AudioDevice.h"

using namespace TextureSynthesis;

AudioDevice::AudioDevice()
    : _curUserDataP(NULL)
{ }

AudioDevice::~AudioDevice()
{ }

void AudioDevice::play(const Signal& signal, std::function<void()> callback)
{
    _curUserDataP = std::shared_ptr<UserData>(new UserData(&signal, callback));

    SDL_AudioSpec desiredSpec, obtainedSpec;
    desiredSpec.freq = signal.sampleRate;
    desiredSpec.format = AUDIO_S16;
    desiredSpec.channels = 1;
    desiredSpec.samples = 1024;
    desiredSpec.callback = writeAudio;
    desiredSpec.userdata = _curUserDataP.get();

    if(SDL_OpenAudio(&desiredSpec, &obtainedSpec) != 0)
        std::cerr << "Error starting audio device: " << SDL_GetError() << "\n\n";
    if(desiredSpec.format != obtainedSpec.format)
        std::cerr << "Couldn't output audio with the given format.\n";
    if(desiredSpec.freq != obtainedSpec.freq)
        std::cerr << "Couldn't output audio with at the signal's sample rate.\n";

    SDL_PauseAudio(0);
}

void AudioDevice::writeAudio(void *userDataVP, Uint8 *stream, int len)
{
    UserData *userDataP = (UserData*)userDataVP;
    const Signal *signalP = userDataP->_signalP;
    std::function<void()> callback = userDataP->_callback;
    int *signalPosP = &(userDataP->_signalPos);

    int signalLen = signalP->size();
    int streamPos;

    short sample;
    for(streamPos = 0; streamPos < len && *signalPosP < signalLen; streamPos+=2)
    {
        sample = (short)std::real((*signalP)[*signalPosP]);
        stream[streamPos] = sample & 0xFF;
        stream[streamPos + 1] = sample >> 8;
        (*signalPosP)++;
    }

    if(*signalPosP >= signalLen)
    {
        for(; streamPos < len; streamPos++)
            stream[streamPos] = 0;
        callback();
    }
}

UserData::UserData(const Signal *signalP, std::function<void()> callback)
    : _signalP(signalP), _callback(callback), _signalPos(0), _isPlaying(true)
{ }
