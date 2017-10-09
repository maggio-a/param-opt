#ifndef TIMER_H
#define TIMER_H

#include <chrono>

struct Timer {
    using hrc = std::chrono::high_resolution_clock;

    hrc::time_point start;
    hrc::time_point last;

    Timer() : start(hrc::now()) { last = start; }

    float TimeElapsed() { last = hrc::now(); return std::chrono::duration<float>(last - start).count(); }

    float TimeSinceLastCheck() { auto t = last; last = hrc::now(); return std::chrono::duration<float>(last - t).count(); }

    void Reset() { start = last = hrc::now(); }

};

#endif // TIMER_H

