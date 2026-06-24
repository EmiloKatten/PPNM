#include"run_timing.h"
#include"evd.h"

void run_timing(){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);

    std::ofstream timingfile("timing.data");
    
    for(int N = 30; N <= 300; N += 30){
        pp::matrix A(N,N);
        for(int i = 0; i < N; i++){
            A(i,i) = dist(gen);
            for(int j = i+1; j < N; j++){
                double value = dist(gen);
                A(i,j) = value;
                A(j,i) = value;
            }
        }
        double lambda_start = dist(gen);
        pp::vector v_start(N);
        for (int i=0; i<N; i++){
            v_start[i] = dist(gen);
        }

        // help from chatbot with timing
        auto start = std::chrono::high_resolution_clock::now();

        auto [lambda, v, steps] = pp::eigen_newton(A, lambda_start, v_start);

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = stop - start;

        timingfile << N << " " << elapsed.count() << "\n";
        std::cout << "\tN = " << N << ",   time = " << elapsed.count() << " sec";
        std::cout<<",   # of calls ="<<steps<<"\n";
    }
    timingfile.close();
}

void run_compare_timing(){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);

    std::ofstream timingfile("comp_timing.data");   

    for(int N = 10; N <= 100; N += 10){
        pp::matrix A(N,N);
        for(int i = 0; i < N; i++){
            A(i,i) = dist(gen);
            for(int j = i+1; j < N; j++){
                double value = dist(gen);
                A(i,j) = value;
                A(j,i) = value;
            }
        }
        double lambda_start = dist(gen);
        pp::vector v_start(N);
        for (int i=0; i<N; i++){
            v_start[i] = dist(gen);
        }

        // help from chatbot with timing
        auto start = std::chrono::high_resolution_clock::now();
        auto [lambda, v, steps] = pp::eigen_newton(A, lambda_start, v_start);
        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = stop - start;

        auto start_EVD = std::chrono::high_resolution_clock::now();
        pp::EVD evdT(A);
        auto stop_EVD = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_EVD = stop_EVD - start_EVD;

        timingfile << N << " " << elapsed.count() << " " << elapsed_EVD.count() << "\n";
        std::cout << "\tN = " << N << ".   Newton: " << elapsed.count() << "s" \
                  << ".   EVD: " << elapsed_EVD.count() << "s\n";
    }
    timingfile.close();
}