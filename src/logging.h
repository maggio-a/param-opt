#ifndef LOGGING_H
#define LOGGING_H

#include <vcg/math/histogram.h>

#include <string>
#include <sstream>
#include <memory>
#include <thread>
#include <mutex>
#include <unordered_map>

#include <cstdio>


#include "gl_utils.h"

struct ParameterizationStrategy;
struct MeshGraph;
struct RasterizedParameterizationStats;
struct Args;

void LogStrategy(const ParameterizationStrategy& strategy, double tol);
void PercentilePlot(vcg::Distribution<double>& d);
void LogDistortionStats(std::shared_ptr<MeshGraph> graph);
void LogParameterizationStats(std::shared_ptr<MeshGraph> graph, const std::vector<RasterizedParameterizationStats>& statsVec);

void LogAggregateStats(const std::string& str, std::shared_ptr<MeshGraph> graph, TextureObjectHandle textureObject);

void LogExecutionParameters(const Args& args, const ParameterizationStrategy& strategy);

// logging macros

#define LOG_INIT(level) (logging::Logger::Init(level))
#define LOG_SET_THREAD_NAME(name) (logging::Logger::RegisterName(name))
#define LOG(level) (level > logging::Logger::GetLogLevel()) ? ((void) 0) : logging::V_() & logging::Buffer(level)

#define LOG_ERR     LOG(logging::Level::Error)
#define LOG_WARN    LOG(logging::Level::Warning)
#define LOG_INFO    LOG(logging::Level::Info)
#define LOG_VERBOSE LOG(logging::Level::Verbose)
#define LOG_DEBUG   LOG(logging::Level::Debug)

namespace logging {

enum Level {
    Error   = -2,
    Warning = -1,
    Info    =  0,
    Verbose =  1,
    Debug   =  2
};

class Buffer {

    std::ostringstream os;

public:

    Buffer(int level);
    ~Buffer();

    template <typename T>
    Buffer& operator<<(const T& t)
    {
        os << t;
        return *this;
    }

    /*
    Buffer& operator<<(std::ostream& (*f)(std::ostream&))
    {
        f(os);
        return *this;
    }
    */

};

struct V_ {
    void operator &(const Buffer&) { }
};

class Logger {

    static int logLevel;
    static std::vector<std::ostream *> streamVec;
    static std::unordered_map<std::thread::id, std::string> threadNames;

    static std::mutex singletonMtx;

public:

    static void Init(int level);

    static int GetLogLevel();
    static std::string GetName();
    static void RegisterStream(std::ostream *os);
    static void RegisterName(const std::string& threadName);
    static void Log(const std::string& s);

};

}

#endif // LOGGING_H

