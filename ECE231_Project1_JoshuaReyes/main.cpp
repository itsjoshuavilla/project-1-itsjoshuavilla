#include <iostream>
#include <iomanip>
using namespace std;

struct _capacitor {
    double *time;
    double *voltage;
    double *current;
    double C;
};
typedef struct _capacitor Capacitor;

bool allocateCapacitor(Capacitor &cap, size_t nSteps) {
    cap.time    = new (nothrow) double[nSteps];
    cap.voltage = new (nothrow) double[nSteps];
    cap.current = new (nothrow) double[nSteps];
    return cap.time && cap.voltage && cap.current;
}

void freeCapacitor(Capacitor &cap) {
    delete[] cap.time;
    delete[] cap.voltage;
    delete[] cap.current;
    cap.time = cap.voltage = cap.current = nullptr;
}

void fillTime(Capacitor &cap, size_t nSteps, double dt) {
    for (size_t i = 0; i < nSteps; ++i) cap.time[i] = i * dt;
}

// Constant current: V(t+1) = V(t) + I(t)*dt/C, with I(t)=I_const
void simulateConstantCurrent(Capacitor &cap, size_t nSteps, double dt, double I_const) {
    cap.current[0] = I_const;
    cap.voltage[0] = 0.0;
    for (size_t i = 1; i < nSteps; ++i) {
        cap.current[i] = I_const;
        cap.voltage[i] = cap.voltage[i - 1] + cap.current[i - 1] * dt / cap.C;
    }
}

// Constant voltage (series R): I(t+1) = I(t) - I(t)*dt/(R*C); Vc from I = C dV/dt
void simulateConstantVoltage(Capacitor &cap, size_t nSteps, double dt, double R, double V0) {
    cap.current[0] = V0 / R;
    cap.voltage[0] = 0.0;
    for (size_t i = 1; i < nSteps; ++i) {
        double I_prev = cap.current[i - 1];
        cap.current[i] = I_prev - (I_prev / (R * cap.C)) * dt;
        cap.voltage[i] = cap.voltage[i - 1] + I_prev * dt / cap.C;
    }
}

void printEvery(const string &title, const Capacitor &cap, size_t nSteps, size_t stride) {
    cout << "\n=== " << title << " ===\n";
    cout << "index, time(s), voltage_Vc(V), current_I(A)\n";
    cout << fixed;
    for (size_t i = 0; i < nSteps; i += stride) {
        cout << i << ", "
             << setprecision(10) << cap.time[i] << ", "
             << setprecision(6)  << cap.voltage[i] << ", "
             << setprecision(6)  << cap.current[i] << "\n";
    }
    if ((nSteps - 1) % stride != 0) {
        size_t i = nSteps - 1;
        cout << i << ", "
             << setprecision(10) << cap.time[i] << ", "
             << setprecision(6)  << cap.voltage[i] << ", "
             << setprecision(6)  << cap.current[i] << "\n";
    }
}

int main() {
    const double dt      = 1e-10;
    const double t_final = 5e-6;
    const size_t nSteps  = static_cast<size_t>(t_final / dt); // 50000
    const double R       = 1e3;
    const double C_val   = 100e-12;
    const double I_const = 1e-2;
    const double V0      = 10.0;
    const size_t stride  = 200;

    Capacitor cc{}, rc{};
    cc.C = C_val; rc.C = C_val;

    if (!allocateCapacitor(cc, nSteps) || !allocateCapacitor(rc, nSteps)) {
        cerr << "Allocation failed.\n";
        freeCapacitor(cc); freeCapacitor(rc);
        return 1;
    }

    fillTime(cc, nSteps, dt);
    fillTime(rc, nSteps, dt);

    simulateConstantCurrent(cc, nSteps, dt, I_const);
    simulateConstantVoltage(rc, nSteps, dt, R, V0);

    printEvery("Constant Current Source (I = constant)", cc, nSteps, stride);
    printEvery("Constant Voltage Source with Series Resistor", rc, nSteps, stride);

    freeCapacitor(cc);
    freeCapacitor(rc);
    return 0;
}
