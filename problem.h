#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

class ProblemData {
public:
    struct Coordinate {
        double x; // X coordinate
        double y; // Y coordinate
    };

    int units_distance = 0;
    Coordinate die_area[4]; // DIEAREA vertices
    Coordinate ff_size;      // FF dimensions (length and width)
    Coordinate buf_size;     // BUF dimensions (length and width)
    Coordinate clk_position; // CLK coordinates
    std::vector<Coordinate> ff_positions; // FF coordinates
    int ff_count = 0;        // Number of FFs

    void display() const;
};

   void readProblemDef(const std::string& filename, ProblemData& data);

void ProblemData::display() const {
        std::cout << "UNITS DISTANCE MICRONS: " << units_distance << std::endl;
        std::cout << "DIEAREA Coordinates:" << std::endl;
        for (const auto& coord : die_area) {
            std::cout << "(" << coord.x << ", " << coord.y << ")" << std::endl;
        }
        std::cout << "FF Size: (" << ff_size.x << ", " << ff_size.y << ")" << std::endl;
        
        std::cout << "BUF Size: (" << buf_size.x << ", " << buf_size.y << ")" << std::endl;
        std::cout << "CLK Position: (" << clk_position.x << ", " << clk_position.y << ")" << std::endl;
        std::cout << "FF Count: " << ff_count << std::endl;
        std::cout << "FF Positions:" << std::endl;
        for (const auto& pos : ff_positions) {
            std::cout << "(" << pos.x << ", " << pos.y << ")" << std::endl;
        }
    }

void readProblemDef(const std::string& filename, ProblemData& data) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;

        // Parse UNITS
        if (line.find("UNITS DISTANCE MICRONS") != std::string::npos) {
            iss >> token >> token >> token >> data.units_distance;
        }

        // Parse DIEAREA coordinates
        else if (line.find("DIEAREA") != std::string::npos) {
         // Skip "DIEAREA" token
         iss >> token;
         for (int i = 0; i < 4; ++i) {
           char skip; // To hold and skip the '(' and ')'
           iss >> skip; // Skip '('
           iss >> data.die_area[i].x >> data.die_area[i].y; // Read the coordinates
           iss >> skip; // Skip ')'
           //将坐标转换为微米单位
           //data.die_area[i].x /= data.units_distance;
           //data.die_area[i].y /= data.units_distance;
         }
        }

        // Parse BUF dimensions
        else if (line.find("BUF (") != std::string::npos) {
            iss >> token >> token >> data.buf_size.x >> data.buf_size.y; // Read BUF dimensions
            //将BUF尺寸转换为微米单位
            //data.buf_size.x /= data.units_distance;
            //data.buf_size.y /= data.units_distance;
        }

        // Parse CLK coordinates
        else if (line.find("CLK (") != std::string::npos) {
            iss >> token >> token >> data.clk_position.x >> data.clk_position.y; // Read CLK coordinates
            //将CLK坐标转换为微米单位
            //data.clk_position.x /= data.units_distance;
            //data.clk_position.y /= data.units_distance;
        }

        // Parse FF components
        else if (line.find("- FF") != std::string::npos) {
            ProblemData::Coordinate ff_pos; // Use ProblemData::Coordinate
            iss >> token >> token >> token >> token; // Skip "- FFx"
            //读取FF左下角坐标
            if (iss >> ff_pos.x >> ff_pos.y) {
              
              //将FF坐标转换为微米单位
             // ff_pos.x /= data.units_distance;
              //ff_pos.y /= data.units_distance;
               
              //计算中心坐标：左下角坐标加上尺寸的一半
              ff_pos.x += data.ff_size.x / 2;
              ff_pos.y += data.ff_size.y / 2;

              data.ff_positions.push_back(ff_pos); 
            } else {
                std::cerr << "Failed to read FF position." << std::endl;
            }
        }

  
       // Parse FF dimensions
       else if (line.find("FF (") != std::string::npos) {
            iss >> token; // Skip "FF"
            char skip; // To hold and skip the '('
            iss >> skip; // Skip '('
            iss >> data.ff_size.x >> data.ff_size.y;
            //将FF尺寸转换为微米单位
            //data.ff_size.x /= data.units_distance;
           //data.ff_size.y /= data.units_distance;
        }

       //删除最后一行
        else if(line.find("END")!=std::string::npos) {
          continue;
        }
       


        // Parse COMPONENTS count
        else if (line.find("COMPONENTS") != std::string::npos) {
            iss >> token >> data.ff_count;
        }
    }

    file.close();
}
