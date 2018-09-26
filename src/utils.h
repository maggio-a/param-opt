#ifndef UTILS_H
#define UTILS_H

#include <cstdlib>
#include <iostream>

#define ensure_condition(expr) \
    ((expr)\
     ? (void) (0) \
     : ensure_condition_fail(#expr, __FILE__, __LINE__))

inline void ensure_condition_fail(const char *expr, const char *filename, unsigned int line)
{
    std::cerr << filename << " (line " << line << "): Failed check `" << expr << "'" << std::endl;
    std::abort();
}

struct Args {
    std::string filename;
    int regionCount;
    bool gui;
    bool filter;
    bool fixContextCapture;

    Args() : filename{}, regionCount{20}, gui{false}, filter{true} {}
};

/*
 * List of valid options:
 *
 *   --gui                 Invoke the GUI
 *
 *   --nofilter            Do not use the push-pull filter when generating the texture
 *
 *   --regioncount=NUM     Use the integer value NUM as the target atlas size in the clustering procedure
 *                         Note that this value is summed to the number of connected components...
 *                         default = 20
 *
 *   --fixcontextcapture   Use a simple topological filter to try to clean troublesome regions
 *
 */
inline Args parse_args(int argc, char *argv[])
{
    Args args;
    for (int i = 1; i < argc; ++i) {
        std::string s(argv[i]);
        if (s.substr(0, 2) == std::string("--")) {
            std::string arg = s.substr(2);
            if (arg == std::string("gui"))
                args.gui = true;
            else if (arg == std::string("nofilter"))
                args.filter = false;
            else if (arg == "regioncount")
                args.regionCount = std::stoi(s.substr(s.find_first_of("=") + 1));
            else if (arg == "fixcontextcapture")
                args.fixContextCapture = true;
            else
                ensure_condition(0 && "Invalid flag option");
        } else {
            ensure_condition(args.filename.empty() && "Invalid argument");
            args.filename = s;
        }
    }
    return args;
}


#endif // UTILS_H

