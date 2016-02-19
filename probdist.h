#include "config.h"
#include <algorithm>
#include <math.h>

#ifndef _PROBDIST_H_
#define _PROBDIST_H_

class poissondd {
	private:
		whole factorial(whole n) {return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}
		
	public:
		real lambda, invstart = 0, invinc = 1;
	
		poissondd (real expval): lambda(expval) {}
	 
		void set(real expval) {lambda = expval;}
		real prob(real val) {return pow(lambda, val) * exp(-lambda) / factorial((whole) val);}
		real randval()
		{
			real sum = 0, suml, samp = (rand() + 1.0) / (RAND_MAX + 1.0), val = invstart;
			
			while (sum < samp)
			{
				suml = sum;
				sum += invinc * pow(lambda, val) * exp(-lambda) / factorial((whole) val);
				val += invinc;
			}
			
			if (fabs(suml - samp) < fabs(sum - samp)) val -= invinc;
			
			return val;
		}
};

class binomialdd {
	private:
		whole factorial(whole n) {return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}
		whole choose(whole n, whole k) {return factorial(n) / (factorial(k) * factorial(n - k));}
		
	public:
		real n, p, invstart = 0, invinc = 1;
	
		binomialdd (real num, real prob): n(num), p(prob) {}
	 
		void set(real num, real prob) {n = num; p = prob;}
		real prob(real val) {return choose(n, val) * pow(p, val) * pow(1 - p, n - val);}
		real randval()
		{
			real sum = 0, suml, samp = (rand() + 1.0) / (RAND_MAX + 1.0), val = invstart;
			
			while (sum < samp)
			{
				suml = sum;
				sum += invinc * choose(n, val) * pow(p, val) * pow(1 - p, n - val);
				val += invinc;
			}
			
			if (fabs(suml - samp) < fabs(sum - samp)) val -= invinc;

			return val;
		}
};

class normalcd {
	public:
		real mu, sigma, invstart = -10, invinc = 0.01;
	
		normalcd (real mean, real stdev): mu(mean), sigma(stdev) {invstart = -5 * stdev;}
	 
		void set(real mean, real stdev) {mu = mean; sigma = stdev; invstart = -5 * stdev;}
		
		real prob(real val1, real val2)
		{
			real sum = 0, suml, val = val1;
			
			while (val < val2)
			{
				suml = sum;
				sum += invinc * 1/(sigma * sqrt(2 * PI)) * exp((val - mu) * (mu - val) / (2 * sigma * sigma));
				val += invinc;
			}
			
			return sum;
		}
		
		real randval()
		{
			real sum = 0, suml, samp = (rand() + 1.0) / (RAND_MAX + 1.0), val = invstart;
			
			while (sum < samp)
			{
				suml = sum;
				sum += invinc * 1/(sigma * sqrt(2 * PI)) * exp((val - mu) * (mu - val) / (2 * sigma * sigma));
				val += invinc;
				//printf("%f %f\n", sum, samp);
			}
			
			if (fabs(suml - samp) < fabs(sum - samp)) val -= invinc;
			
			return val;
		}
};

class gammacd {
	private:
		real gamma(real arg)
		{
			real sum = 0, suml = -1, val = 0;
			whole pastmax = 0;
			
			while (sum - suml > 0.000001 || !pastmax)
			{
				suml = sum;
				sum += 0.01 * pow(val, arg - 1) * exp(-val);
				val += 0.01;
				if (pow(val, arg - 1) * exp(-val) > pow(val + 0.01, arg - 1) * exp(-val - 0.01)) pastmax = 1;
			}

			return sum;
		}
		
	public:
		real k, theta, gammak, invstart = 0, invinc = 0.01;
	
		gammacd (real shape, real scale): k(shape), theta(scale) {gammak = gamma(k);}
	 
		void set(real shape, real scale) {k = shape; theta = scale; gammak = gamma(k);}
		
		real prob(real val1, real val2)
		{
			real sum = 0, suml, val = val1;
			
			while (val < val2)
			{
				suml = sum;
				sum += invinc * 1 / (gammak * pow(theta, k)) * pow(val, k - 1) * exp(-val / theta);
				val += invinc;
			}
			
			return sum;
		}
		
		real randval()
		{
			real sum = 0, suml, samp = (rand() + 1.0) / (RAND_MAX + 1.0), val = invstart;
			
			while (sum < samp)
			{
				suml = sum;
				sum += invinc * 1 / (gammak * pow(theta, k)) * pow(val, k - 1) * exp(-val / theta);
				val += invinc;
			}
			
			if (fabs(suml - samp) < fabs(sum - samp)) val -= invinc;

			return val;
		}
};

class lognormcd {
	public:
		real mu, sigma, invstart = 0.01, invinc = 0.01;
	
		lognormcd (real mean, real stdev): mu(mean), sigma(stdev) {}
	 
		void set(real mean, real stdev) {mu = mean; sigma = stdev;}
		
		real prob(real val1, real val2)
		{
			real sum = 0, suml, val = val1;
			
			while (val < val2)
			{
				suml = sum;
				sum += invinc * 1 / (sigma * val * sqrt(2 * PI)) * exp((log(val) - mu) * (mu - log(val)) / (2 * sigma * sigma));
				val += invinc;
			}
			
			return sum;
		}
		
		real randval()
		{
			real sum = 0, suml, samp = (rand() + 1.0) / (RAND_MAX + 1.0), val = invstart;
			
			while (sum < samp)
			{
				suml = sum;
				sum += invinc * 1 / (sigma * val * sqrt(2 * PI)) * exp((log(val) - mu) * (mu - log(val)) / (2 * sigma * sigma));
				val += invinc;
			}
			
			if (fabs(suml - samp) < fabs(sum - samp)) val -= invinc;
			
			return val;
		}
};

class weibullcd {
	public:
		real lambda, k, invstart = 0, invinc = 0.01;
	
		weibullcd (real scale, real shape): k(shape), lambda(scale) {}
	 
		void set(real scale, real shape) {k = shape; lambda = scale;}
		
		real prob(real val1, real val2)
		{
			real sum = 0, suml, val = val1;
			
			while (val < val2)
			{
				suml = sum;
				sum += invinc * k / lambda * pow(val / lambda, k - 1) * exp(-pow(val / lambda, k));
				val += invinc;
			}
			
			return sum;
		}
		
		real randval()
		{
			real sum = 0, suml, samp = (rand() + 1.0) / (RAND_MAX + 1.0), val = invstart;
			
			while (sum < samp)
			{
				suml = sum;
				sum += invinc * k / lambda * pow(val / lambda, k - 1) * exp(-pow(val / lambda, k));
				val += invinc;
			}
			
			if (fabs(suml - samp) < fabs(sum - samp)) val -= invinc;
			
			return val;
		}
};

#endif // #ifndef _PROBDIST_H_
