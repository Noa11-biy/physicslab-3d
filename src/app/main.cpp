#include <cstdio>
#include <cstring>

#include "Application.hpp"

int main(int argc, char** argv) {
    pl::AppOptions options;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--smoke-test") == 0) {
            options.smokeTest = true;
        } else {
            std::fprintf(stderr, "Usage : %s [--smoke-test]\n", argv[0]);
            return std::strcmp(argv[i], "--help") == 0 ? 0 : 2;
        }
    }
    return pl::runApplication(options);
}
