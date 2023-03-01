#include <iostream>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h> // numpy arrays
#include <pybind11/stl.h> // std::vector conversion
#include <numbers> // std::numbers

namespace py = pybind11;

class MicroSimulation
{
public:
    MicroSimulation(int sim_id);
    void initialize();
    // solve takes python dict for macro_write data, dt, and returns python dict for macro_read data
    py::dict solve(py::dict macro_write_data, double dt);
    void save_checkpoint();
    void reload_checkpoint();
    int get_dims();

private:
    const double pi_ = 3.14159265358979323846;
    int _sim_id;
    int _dims;
    //double _micro_scalar_data;
    //std::vector<double> _micro_vector_data;
    double _k_00;
    double _k_01;
    double _k_10;
    double _k_11;
    double _porosity;
    double _checkpoint;
};

// Constructor
MicroSimulation::MicroSimulation(int sim_id) : _sim_id(sim_id), _dims(3), _k_00(0), _k_01(0), _k_10(0),_k_11(0),_porosity(0), _checkpoint(0) {}

// Initialize
void MicroSimulation::initialize()
{
    std::cout << "Initialize micro problem (" << _sim_id << ")\n";
    //_micro_scalar_data = 0;
    _k_00 = 0;
    _k_01 = 0;
    _k_10 = 0;
    _k_11 = 0;
    //_micro_vector_data.clear();
    _checkpoint = 0;
}

// Solve
py::dict MicroSimulation::solve(py::dict macro_write_data, double dt)
{
    std::cout << "Solve timestep of micro problem (" << _sim_id << ")\n";

    // assert(dt != 0);
    if (dt == 0)
    {
        std::cout << "dt is zero\n";
        exit(1);
    }

    //! Here, insert your code, changing the data and casting it to the correct type
    // create double variable from macro_write_data["micro_scalar_data"]; which is a python float
    //double macro_scalar_data = macro_write_data["macro-scalar-data"].cast<double>();
    double conc = macro_write_data["concentration"].cast<double>();

    //TODO set _k_00 etc
    _k_00 = conc; //TODO remove
    // create python dict for micro_write_data
    py::dict micro_write_data;
    // add micro_scalar_data and micro_vector_data to micro_write_data
    micro_write_data["k_00"] = _k_00;
    micro_write_data["k_10"] = _k_10;
    micro_write_data["k_01"] = _k_01;
    micro_write_data["k_11"] = _k_11;
    micro_write_data["porosity"] = _porosity;
    micro_write_data["grain_size"] = std::sqrt((1-_porosity)/pi_);
    
    // return micro_write_data
    return micro_write_data;
}
// Save Checkpoint
void MicroSimulation::save_checkpoint()
{
    std::cout << "Saving state of micro problem (" << _sim_id << ")\n";
    _checkpoint = _k_00; //TODO
}

// Reload Checkpoint
void MicroSimulation::reload_checkpoint()
{
    std::cout << "Reverting to old state of micro problem (" << _sim_id << ")\n";
    _k_00 = _checkpoint;
}

int MicroSimulation::get_dims()
{
    return _dims;
}

PYBIND11_MODULE(micro_sim, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    py::class_<MicroSimulation>(m, "MicroSimulation")
        .def(py::init<int>())
        .def("initialize", &MicroSimulation::initialize)
        .def("solve", &MicroSimulation::solve)
        .def("save_checkpoint", &MicroSimulation::save_checkpoint)
        .def("reload_checkpoint", &MicroSimulation::reload_checkpoint)
        .def("get_dims", &MicroSimulation::get_dims);
}

// compile with 
// c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) cpp_dummy.cpp -o cpp_dummy$(python3-config --extension-suffix)
// then from the same directory run python3 -c "import cpp_dummy; cpp_dummy.MicroSimulation(1)