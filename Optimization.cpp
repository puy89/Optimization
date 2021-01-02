
#include <iostream>
#include "Func.h"

int main()
{
    Var var(1);
    double cs1[] = { 1, 1 };
    Poly p1(var, cs1, 1);
    double cs2[] = { 6, -5, 1 };
    Poly p2(var, cs2, 2);
    //f = (x+1)(x2-5x+6)
    Func& f = p1 * p2;
    double a = 1, x0 = -10;
    Array root(x0);
    std::cout << "root1 = " << f.newton_raphson(root, 10).toString() << "\n";
    root = 1;
    std::cout << "root2 = " << f.newton_raphson(root, 10).toString() << "\n";
    root = 10;
    std::cout << "root3 = " << f.newton_raphson(root, 10).toString() << "\n";
}
