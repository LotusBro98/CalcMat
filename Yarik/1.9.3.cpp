#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdint.h>

double calcC(double eps)
{
	int64_t n = 1;
	double CCur = 100500;
	double CPrev;
	double sumPrev = 0;
	double diff;
	do
	{
		CPrev = CCur;

		CCur = 0;
		for (int64_t k = n; k > n / 2; k--)
			CCur += 1.0 / k;
		CCur += sumPrev;
		sumPrev = CCur;
		
		CCur -= std::log(n);
		
		diff = std::abs(CCur - CPrev);
		std::cout << CCur << " " << diff << " " << n << "\n";
		
		n *= 2;
	}
	while (diff >= eps);

	return CCur;
}

int main()
{
	std::cout << std::setprecision(10) << calcC(1e-10) << "\n";
	
	return 0;
}
