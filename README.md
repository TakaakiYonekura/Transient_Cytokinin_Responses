🚧 This project **IS** under review. Please **DO NOT** cite until the corresponding paper is published.

# Transient Cytokinin Responses

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview
This project provides a simulation code for **parameter fitting of a mathematical model describing cytokinin-driven gene expression and demonstrates transient peak of cytokinin responses**.

- Developed using **Visual C++** with Eigen and NLopt libraries.
- Implements numerical optimization and mathematical modeling techniques.
- Designed for **simulating gene expression dynamics**.
- Includes a Jupyter Notebook (`Cytokinin_Plot.ipynb`) for visualizing simulation results.

## Requirements
- **OS**: Windows 10/11
- **Compiler**: Microsoft Visual Studio **2019 or later** (MSVC)
- **CMake**: Version 3.10 or later
- **Dependencies**:
  - **Eigen** (3.3.9) - [Eigen Official Repository](https://gitlab.com/libeigen/eigen)
  - **NLopt** (2.7.1) - [NLopt Official Repository](https://github.com/stevengj/nlopt)
  - Standard C++ libraries (`<iostream>`, `<fstream>`, `<vector>`, etc.)
- **Jupyter Notebook and Python packages** (installed via `requirements.txt`)


## Installation
Follow these steps to set up the project.

### 1. Install Required Libraries
It is recommended to clone NLopt and Eigen into the `third_party` directory within the project.
```sh
mkdir -p third_party
cd third_party
```

#### **Eigen Installation**
Clone the Eigen repository:
```sh
git clone https://gitlab.com/libeigen/eigen.git
```

#### **NLopt Installation**
- Download the precompiled binaries from [NLopt GitHub](https://github.com/stevengj/nlopt).
- Or build from source:
  ```sh
  git clone https://github.com/stevengj/nlopt.git
  cd nlopt
  mkdir build
  cd build
  cmake .. -G "Visual Studio 17 2022"
  cmake --build . --config Release
  ```

### **2. Build the Project**
Open the project in Visual Studio.
Apply CMakeLists.txt settings.
Build the project using Ctrl + Shift + B.

### **Install Required Libraries (Jupyter Notebook)** 
To run the Jupyter Notebook, install the necessary Python packages:
```sh
pip install -r requirements.txt
```

## **Usage**
Run the program with:
  ```sh
  ./FindSolution_Cytokinin.exe
  ```
### Visualizing Results with Jupyter Notebook
To visualize the simulation results, open the Jupyter Notebook:
  ```sh
  jupyter notebook Cytokinin_Plot.ipynb
  ```
The notebook provides:
- Plots of the simulated gene expression response
- Analysis of parameter fitting results

## **Repository Structure**
  ```markdown
  /project_root
   ├── src/          # C++ source code
   │    ├── main.cpp
   │    ├── difAve.cpp
   │    ├── DifferentialEquation.cpp
   │    ├── readCSV.cpp
   │    └── stdafx.cpp
   ├── include/      # Header files
   │    ├── difAve.h
   │    ├── DifferentialEquation.h
   │    ├── readCSV.h
   │    └── stdafx.h
   ├── Visualization/          # visualization code with jupyter notebook
   │    ├── Cytokinin_Plot.ipynb
   │    └── requirements.txt
   ├── CMakeLists.txt
   ├── README.md
   └── LICENSE
  ```

## **References**
- [Eigen Documentation](https://eigen.tuxfamily.org/)
- [NLopt Documentation](https://nlopt.readthedocs.io/en/latest/)
- [CMake Official Documentation](https://cmake.org/documentation/)
- [Microsoft Visual Studio C++ Documentation](https://learn.microsoft.com/en-us/cpp/)
- [GitHub - NLopt Repository](https://github.com/stevengj/nlopt)

## **License**
This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.

## **Citation**
Shunji Shimadzu, Takaaki Yonekura, Tomoyuki Furuya, Mikiko Kojima, Kimitsune Ishizaki, Masashi Asahina, Kyoko Ohashi-Ito, Hitoshi Sakakibara, Hidehiro Fukaki, Hiroo Fukuda, Yuki Kondo. A cytokinin response maximum induces and activates bifacial stem cells for radial growth.

🚧 Under review yet.
