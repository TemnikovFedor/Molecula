#include "particles.h"

int main(int argc, char const *argv[]) {
    std::cout << "Before start" << std::endl;
    Particles p;
    p.moving();
    std::cout << "Before end" << std::endl;
    return 0;
}