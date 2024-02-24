#include <iostream>
#include <vector>

int main()
{

    std::vector<double> v(5);

    for(int i = 0; i <= 5; i++)
        v.at(i) = i;
    


    return 0;
}