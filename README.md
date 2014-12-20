TextureAnimation
================

The results of a research project on texture animation, under the supervision
of Professor Changxi Zheng at Columbia University.

Table of Contents
I.      Summary
II.     Results
III.    Running the code
IV.     Navigating the code



I. Summary
This project was a C++ implementation of the audio texture synthesis work done by
Josh McDermott and Eero Simoncelli in their paper which can be found here:
http://mcdermottlab.mit.edu/papers/McDermott_Simoncelli_2011_sound_texture_synthesis.pdf

There you can find a full summary of the work being done.  At a high level, the aim
is to take an audio file as input, which stores a sound texture.  This includes
sounds like rain, fire, crickets, murmuring, streams, and so forth.  The output should
be a sound of the same nature, but which is distinct from the original.  A listener
could hear the output and recognize it as coming from the same source as the input
sound.

My C++ implementation is based in part on the paper and in part on the Matlab
implementation which can be found on their website along with the paper.  The purpose
of a C++ implementation is greater efficiency and a result that can be included in
a wider variety of practical projects, such as the sound synthesis work being done
by this class's professor.



II. Results
You can find results in the samples/ directory.  The only results presented is a short
clip of a stream, which is missing many of the low frequencies from the input, but is
otherwise recognizable.  I tried using other inputs and longer lengths, but in both
cases my implementation failed to converge and/or produced garbage.



III. Running the code
This part is a bit tricky.  You'll have to install cmake, GSL, SDL 2.0, SDL Mixer, and
my modified version of the open-source Aquilla digital signal processing library.
Note also that the tutorial presented has only been tested on Ubuntu 14.04.
To install GSL, visit http://gnu.mirrors.pair.com/gnu/gsl/ and grab gsl-latest.tar.gz
to compile from source.  To install SDL and SDL Mixer, do the same from
https://www.libsdl.org/download-2.0.php and http://www.libsdl.org/projects/SDL_mixer/
respectively.  Finally, on to Aquilla.  Navigate to aquilla-src/ and run "cmake ." then
"make" and finally "sudo make install".  After this is done, there is one last library
to move into place.  Copy aquilla-src/lib/libOoura_fft.a into a directory in your
LD_LIBRARY_PATH environment variable (I chose /usr/local/lib).  You may have to
add the directory to the variable yourself.  Once all this is done, you can navigate
to the top level directory, run "cmake ." and "make" to compile.  To run the code,
use "bin/TextureSynthesis samples/stream_source.wav samples/stream_out.wav" and it
will begin!



IV. Navigating the code
The root directory for the source code is, unsurprisingly, src/.  If you want to
explore the code yourself a bit, the main function is in test.cpp. The meat of the
code really comes from the Synthesis/TextureSynthesizer class, where everything comes
together.  It's quite a large project, admittedly, so it will take some time to work
through.  Email me at jettca1@gmail.com if you have any questions!
