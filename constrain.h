#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class ConstraintsData {
public:
    double net_unit_r = 0;
    double net_unit_c = 0;
    double max_net_rc = 0;
    double max_fanout = 0;
    double buffer_delay = 0;

    void display() const;
};

void readConstraints(const std::string& filename, ConstraintsData& Ddata);

void ConstraintsData::display() const{
        std::cout << "Net Unit R: " << net_unit_r << std::endl;
        std::cout << "Net Unit C: " << net_unit_c << std::endl;
        std::cout << "Max Net RC: " << max_net_rc << std::endl;
        std::cout << "Max Fanout: " << max_fanout << std::endl;
        std::cout << "Buffer Delay: " << buffer_delay << std::endl;
    }

void readConstraints(const std::string& filename, ConstraintsData& data) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Remove spaces around the line
        line.erase(0, line.find_first_not_of(" "));
        line.erase(line.find_last_not_of(" ") + 1);
        
        // Check if the line contains '='
        if (line.find('=') != std::string::npos) {
            std::istringstream iss(line);
            std::string key;
            double value;
            std::getline(iss, key, '=');
            key.erase(0, key.find_first_not_of(" ")); // Remove leading spaces
            key.erase(key.find_last_not_of(" ") + 1); // Remove trailing spaces

            if (iss >> value) {
                if (key == "net_unit_r") data.net_unit_r = value;
                else if (key == "net_unit_c") data.net_unit_c = value;
                else if (key == "max_net_rc") data.max_net_rc = value - 1;
                else if (key == "max_fanout") data.max_fanout = value;
                else if (key == "buffer_delay") data.buffer_delay = value;
            } else {
                std::cerr << "Failed to read value for key: " << key << std::endl;
            }
        } else {
            std::cerr << "Missing '=' in line: " << line << std::endl; // Debug output
        }
    }

    file.close();
}
