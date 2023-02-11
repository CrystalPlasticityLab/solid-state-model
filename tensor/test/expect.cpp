#include "expect.h"
int expect(bool expr, std::string msg) {
    if (expr){
        // 0 for background Color(Black)
        // A for text color(Green)
        system("Color 0A");
        std::cout << "    PASSED: " << msg  << std::endl;
        return 1;
    } else {
        system("Color 0A");
        Color::Modifier red(Color::FG_RED);
        std::cout << "NOT PASSED: " << msg << std::endl;
        return 0;
    }
}