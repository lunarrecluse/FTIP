install instructions:

1) build a directory named "build" and cd to it
mkdir build && cd build

2) configure
cmake -DCMAKE_BUILD_TYPE=Release ..

3) build
make

4) copy ../config.txt to current directory
cp ../config.txt .

5) edit the config.txt according to the info inside that file.

6) run the executable
e.g., when at the build directory: ./protein_registration
