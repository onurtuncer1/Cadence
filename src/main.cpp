#include "Cadence/Vector.h"
#include <iostream>

int main() {
    using namespace Cadence;

    constexpr Vector3<double> a{1.0, 2.0, 3.0};
    constexpr Vector3<double> b{4.0, 5.0, 6.0};
    constexpr auto c = a + b;

    std::cout << "c = (" << c.x() << ", " << c.y() << ", " << c.z() << ")\n";
    std::cout << "dot = " << a.dot(b) << "\n";
    std::cout << "norm(a) = " << a.norm() << "\n";

    return 0;
}
