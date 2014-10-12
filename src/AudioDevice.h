#include <SDL2/SDL.h>
#include <SDL2/SDL_mixer.h>

#include "Signal.h"

namespace TextureSynthesis
{
    class UserData
    {
    public:
        UserData(const Signal *signal, std::function<void()> callback);
        const Signal *_signalP;
        std::function<void()> _callback;
        int _signalPos;
        bool _isPlaying;
    };

    class AudioDevice
    {
    public:
        AudioDevice();
        void play(const Signal& signal, std::function<void()> callback);
        ~AudioDevice();

    private:
        std::shared_ptr<UserData> _curUserDataP;

        static void writeAudio(void *userDataVP, Uint8 *stream, int len);
    };
}
