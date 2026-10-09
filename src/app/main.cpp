#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "Application.hpp"
#include "physicslab/core/Level.hpp"

int main(int argc, char** argv) {
    pl::AppOptions options;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--smoke-test") == 0) {
            options.smokeTest = true;
        } else if (std::strcmp(argv[i], "--gpu-test") == 0) {
            options.gpuTest = true;
        } else if (std::strcmp(argv[i], "--gpu-max-n") == 0 && i + 1 < argc) {
            options.gpuMaxN = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--level") == 0 && i + 1 < argc) {
            options.level = std::atoi(argv[++i]);
            if (options.level < 1 || options.level > pl::kLevelCount) {
                std::fprintf(stderr, "Le niveau doit être compris entre 1 et %d\n", pl::kLevelCount);
                return 2;
            }
        } else if (std::strcmp(argv[i], "--sim") == 0 && i + 1 < argc) {
            options.simulation = std::atoi(argv[++i]);
        } else {
            std::fprintf(stderr, "Usage : %s [--level 1..6] [--sim 1..10] [--smoke-test] [--gpu-test [--gpu-max-n N]]\n", argv[0]);
            return std::strcmp(argv[i], "--help") == 0 ? 0 : 2;
        }
    }
    return pl::runApplication(options);
}
