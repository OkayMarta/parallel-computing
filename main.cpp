#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <limits>
#include <chrono>
#include <thread>

#ifdef _WIN32
#include <windows.h>
#endif

struct Point {
    double x, y;
};

struct SearchResult {
    Point min_p1, min_p2;
    Point max_p1, max_p2;
    double min_dist;
    double max_dist;
};

// ---------------------------------------------------------
// Допоміжні функції
// ---------------------------------------------------------

// Квадрат відстані
double distSq(const Point& a, const Point& b) {
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

// Генерація випадкових точок
std::vector<Point> generatePoints(int n) {
    std::vector<Point> points(n);
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);
    for (int i = 0; i < n; ++i) {
        points[i].x = dist(gen);
        points[i].y = dist(gen);
    }
    return points;
}

// ---------------------------------------------------------
// РОБОТА З ФАЙЛАМИ
// ---------------------------------------------------------

// Збереження вхідних даних у файл
void saveToFile(const std::vector<Point>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    file << points.size() << "\n";
    file << std::fixed << std::setprecision(6);
    for (const auto& p : points) {
        file << p.x << " " << p.y << "\n";
    }
    file.close();
}

// Зчитування вхідних даних із файлу
std::vector<Point> loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<Point> points;
    if (!file.is_open()) return points;
    int n;
    file >> n;
    points.resize(n);
    for (int i = 0; i < n; ++i) {
        file >> points[i].x >> points[i].y;
    }
    file.close();
    return points;
}

// ---------------------------------------------------------
// 1. ПОСЛІДОВНИЙ АЛГОРИТМ
// ---------------------------------------------------------
SearchResult findPairsSequential(const std::vector<Point>& points) {
    int n = points.size(); // Дізнаємося, скільки всього точок (75 000)
    
	double min_dist_sq = std::numeric_limits<double>::max(); // Початково задаю величезне число
    double max_dist_sq = -1.0; // Початково -1, бо квадрат відстані завжди ≥0
    int min_i = 0, min_j = 0, max_i = 0, max_j = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
			// Рахуємо квадрат відстані між точкою i та точкою j
            double d_sq = distSq(points[i], points[j]);

			// Перевіряємо, чи це новий рекорд мінімуму
            if (d_sq < min_dist_sq) {
                min_dist_sq = d_sq;
                min_i = i; min_j = j;
            }
			// Перевіряємо, чи це новий рекорд максимуму
            if (d_sq > max_dist_sq) {
                max_dist_sq = d_sq;
                max_i = i; max_j = j;
            }
        }
    }

    SearchResult res;
    res.min_p1 = points[min_i]; res.min_p2 = points[min_j];
    res.max_p1 = points[max_i]; res.max_p2 = points[max_j];
    res.min_dist = std::sqrt(min_dist_sq);
    res.max_dist = std::sqrt(max_dist_sq);
    return res;
}

// ---------------------------------------------------------
// 2. ПАРАЛЕЛЬНИЙ АЛГОРИТМ: Функція для одного потоку
// ---------------------------------------------------------
void workerThread(const std::vector<Point>& points, int start_idx, int step, SearchResult& local_res) {
    int n = points.size();

	// Ініц. локальні рекорди для цього потоку: кожен потік шукає "свій" min і max незалежно від інших
    double local_min_sq = std::numeric_limits<double>::max();
    double local_max_sq = -1.0;
    int min_i = 0, min_j = 0, max_i = 0, max_j = 0;

    // Циклічний розподіл (Interleaved) для балансування навантаження
    for (int i = start_idx; i < n; i += step) {
        for (int j = i + 1; j < n; ++j) {
            double d_sq = distSq(points[i], points[j]);
            if (d_sq < local_min_sq) {
                local_min_sq = d_sq;
                min_i = i; min_j = j;
            }
            if (d_sq > local_max_sq) {
                local_max_sq = d_sq;
                max_i = i; max_j = j;
            }
        }
    }

    // Зберігаємо ЛОКАЛЬНІ квадрати відстаней (корінь беру в кінці)
    local_res.min_dist = local_min_sq;
    local_res.max_dist = local_max_sq;
    local_res.min_p1 = points[min_i]; local_res.min_p2 = points[min_j];
    local_res.max_p1 = points[max_i]; local_res.max_p2 = points[max_j];
}

// ---------------------------------------------------------
// ПАРАЛЕЛЬНИЙ АЛГОРИТМ: Головна функція
// ---------------------------------------------------------
SearchResult findPairsParallel(const std::vector<Point>& points, int num_threads) {
    std::vector<std::thread> threads; // Список для самих потоків
    std::vector<SearchResult> local_results(num_threads); // Масив результатів які повернуть потоки

    // Запускаємо потоки
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back(workerThread, std::ref(points), t, num_threads, std::ref(local_results[t]));
    }

    // Чекаємо завершення всіх потоків
    for (auto& th : threads) {
        th.join();
    }

    // Об'єднуємо результати з усіх потоків
    SearchResult global_res = local_results[0];
    for (int t = 1; t < num_threads; ++t) {
        if (local_results[t].min_dist < global_res.min_dist) {
            global_res.min_dist = local_results[t].min_dist;
            global_res.min_p1 = local_results[t].min_p1;
            global_res.min_p2 = local_results[t].min_p2;
        }
        if (local_results[t].max_dist > global_res.max_dist) {
            global_res.max_dist = local_results[t].max_dist;
            global_res.max_p1 = local_results[t].max_p1;
            global_res.max_p2 = local_results[t].max_p2;
        }
    }

    // Беремо корінь з фінальних глобальних значень
    global_res.min_dist = std::sqrt(global_res.min_dist);
    global_res.max_dist = std::sqrt(global_res.max_dist);
    return global_res;
}

int main() {
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif

    int N = 75000;

    std::cout << "--- Згенерую " << N << " точок ---\n";
    std::vector<Point> points = generatePoints(N);
    
    // Зберігання точок у файл для пункту 3
    saveToFile(points, "points_data.txt");
    std::cout << "(Дані збережено у файл points_data.txt)\n";

    // ---------------------------------------------------------
    // 1. ПОСЛІДОВНИЙ ЗАМІР
    // ---------------------------------------------------------
    std::cout << "\n1. Запускаю ПОСЛІДОВНИЙ алгоритм... (очікуйте ~5 сек)\n";
    auto start_seq = std::chrono::high_resolution_clock::now();
    
    SearchResult seq_res = findPairsSequential(points);
    
    auto end_seq = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_seq = end_seq - start_seq;
    
    std::cout << "Відстань найближчих: " << seq_res.min_dist << "\n";
    std::cout << "Відстань найвіддаленіших: " << seq_res.max_dist << "\n";
    std::cout << ">> ЧАС ПОСЛІДОВНОГО: " << time_seq.count() << " секунд\n";

    // ---------------------------------------------------------
    // 2. ПАРАЛЕЛЬНИЙ ЗАМІР (БЕНЧМАРК)
    // ---------------------------------------------------------
	std::cout << "\n==========================================================\n";
    std::cout << "   ПАРАЛЕЛЬНИЙ алгоритм (Залежність від к-сті потоків)\n";
    std::cout << "==========================================================\n";
    
    // Шапка таблиці (Вирівняна пробілами вручну, щоб обійти баг C++ з кирилицею)
    std::cout << "    Потоки      Час (сек)    Прискорення      Перевірка\n";
    std::cout << "----------------------------------------------------------\n";

    // Масив з кількістю потоків для тестування
    std::vector<int> test_threads = {1, 2, 3, 4, 6, 8, 12, 16};

    // Запускаємо паралельну функцію для кожної кількості потоків
    for (int t : test_threads) {
        auto start_par = std::chrono::high_resolution_clock::now();
        
        SearchResult par_res = findPairsParallel(points, t);
        
        auto end_par = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_par = end_par - start_par;

        // Рахуємо прискорення (час послідовного ділимо на час паралельного)
        double speedup = time_seq.count() / time_par.count();

        // Перевіряємо правильність результату
        std::string check = "ПОМИЛКА";
        if (std::abs(seq_res.min_dist - par_res.min_dist) < 1e-6) {
            check = "ОК";
        }

        // Виводимо рядок таблиці (тут setw працює ідеально, бо цифри)
        std::cout << std::setw(10) << t 
                  << std::setw(15) << std::fixed << std::setprecision(4) << time_par.count() 
                  << std::setw(14) << std::fixed << std::setprecision(2) << speedup << "x"
                  << std::setw(15) << check << "\n";
    }
    std::cout << "==========================================================\n";

    return 0;
}