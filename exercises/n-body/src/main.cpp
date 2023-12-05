#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

auto main(int argc, char *argv[]) -> int {
    if (!argc) {
        std::cerr << "Needed argument\n";
        return -1;
    }

    std::ifstream file(argv[1]);

    if (!file.is_open()) {
        std::cerr << "Unable to open file\n";
        return -1;
    }

    std::string line;
    while (getline(file, line)) {
        std::cout << line << '\n';
        std::vector<std::string> words;

        // we have 10 values
        // Arrays no, Masses[i], x[i], y[i], z[i], Vx[i], Vy[i], Vz[i], softening[i], potential[i]
        int i = 0;
        std::stringstream ss(line);
        std::string word;
        while (ss >> word) {
            words.push_back(word);
            i++;
        }
        // Printing the words
        for (const auto &w : words) {
            std::cout << w << std::endl;
        }
    }
    file.close();
}
