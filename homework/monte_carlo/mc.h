#include <vector>
#include <random>
#include <cmath>

using vec = std::vector<double>;

namespace pp
{

    std::vector<int> prime_numbers(int n){
        std::vector<int> primes; int candidate = 2;
        while ((int)primes.size() < n){
            bool candidate_is_prime = true;
            for (int i=0; i<(int)primes.size(); i++){
                if (primes[i]*primes[i]>candidate){break;}
                if (candidate%primes[i]==0) {candidate_is_prime = false; break;}
            }
            if (candidate_is_prime == true) {primes.push_back(candidate);}
            candidate += 1;
        }
        return primes;
    };

    std::pair<double,double> plainmc(auto f, vec a, vec b, int N, auto rng){
        int dim = a.size();
        double V = 1.0; 
        for(int i=0; i<dim; i++){V*=b[i]-a[i];};
        double sum1 = 0; double sum2 = 0;

        vec x(dim);
        for(int i=0; i<N; i++){
            for(int k=0; k<dim; k++){
                x[k] = a[k] + rng.next() * (b[k] - a[k]);
            }
            double fx = f(x); sum1 += fx; sum2 += fx*fx;
        }
        double mean = sum1 / N;
        double sigma = std::sqrt(sum2 / N - mean * mean);
        return std::pair<double,double>(mean * V, sigma * V / std::sqrt(N));
    };

    double quasimc(auto f, vec a, vec b, int N){
        int dim = a.size();
        double V = 1.0; 
        for(int i=0; i<dim; i++){V*=b[i]-a[i];};
        double sum1 = 0;

        vec x(dim);
        std::vector<int> base = prime_numbers(dim);

        for(int i=0; i<N; i++){
            for(int k=0; k<dim; k++){
                
                int n = i; double q = 0.0; int b0 = base[k]; double bk = 1.0/b0;
                while(n>0){ q+=n % b0 * bk; n /= b0; bk/=b0; };

                x[k] = a[k] + q * (b[k] - a[k]);
            }
            double fx = f(x); sum1 += fx;
        }
        double mean = sum1 / N;
        return mean * V;
    };

} // namespace pp
