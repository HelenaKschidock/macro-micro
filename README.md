# macro-micro

## Setup of the macro-micro heat simulation
1. Follow the setup instructions of the [dumux-adapter](https://github.com/precice/dumux-adapter "dumux-adapter"). You folder structure should look something like this: `dumux/dumux-adapter/`.
2. Navigate into `dumux-adapter/examples/` and git clone this repo into `dumux-adapter/examples/macro-micro`.
3. Add the line `add_subdirectory(macro-micro)` to `dumux-adapter/examples/CMakeLists.txt`.

## Building the executables
1. Rebuild the dumux-adapter by navigating into your `dumux` folder and executing `dunecontrol --only=dumux-adapter all`. Repeat this step any time you modify the CMake files.
2. Build the nutils micro executable by navigating into `dumux/dumux-adapter/build-cmake/examples/macro-micro/python-micro-heat` and executing `make copy_micro_heat_python`. Repeat this step any time you modify the micro simulation files.
3. Build the DuMuX macro executable by navigating into `dumux/dumux-adapter/build-cmake/examples/macro-micro/macro-heat` and executing `make test_macro_heat`. 

## Running the simulating
### Running only the macro simulation
* Set the `RunWithCoupling ` parameter in `macro-heat/params.input`to `false`. Rebuild afterwards.
* Navigate to `dumux/dumux-adapter/build-cmake/examples/macro-micro/macro-heat` and execute `./test_macro_heat`
### Running the coupled simulation
* Set the `RunWithCoupling ` parameter in `macro-heat/params.input`to `true`. Rebuild afterwards.
* Open two separate shells. In one of them, run the macro simulation as above.
* In the other, navigate to `dumux/dumux-adapter/build-cmake/examples/macro-micro/python-micro-heat` and execute `python3 run-micro-manager.py`.
* To run the micro code in parallel, use `mpirun -n <NThreads> python3 run-micro-problems.py` instead.