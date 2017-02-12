//
// Created by George Oblapenko (kunstmord@kunstmord.com)
//

#ifndef easy_numeric_hpp
#define easy_numeric_hpp

#include <cmath>
#include <functional>

    inline double integrate_semi_inf(std::function<double(double) > f, double a = 0, int subdivisions = 5, double* error_estimate = nullptr)  {
    // integrate the function f on the semi-infinite interval [a, +infinity)
    // it performs a change of variables: x = a + (1-t)/t and then integrates over the interval [0,1]
    // the interval is split into equal-sized subintervals (controlled by the subdivisions parameter) and the integral is computed over each
    // subinterval using the 15-point Kronrod rule
    double result = 0., tmp_error;
    double step_size = 1. / subdivisions;

    auto f_transform = [f, a](double t) {
        return f(a + (1 - t) / t) / (t * t);
    };

    for (int i = 0; i < subdivisions; i++) {
        if (error_estimate == nullptr) {
            result += integrate_interval(f_transform, i*step_size, (i + 1) * step_size);
        } else {
            result += integrate_interval(f_transform, i*step_size, (i + 1) * step_size, &tmp_error);
            *error_estimate += tmp_error;
        }
    }

    return result;
}

    // integrate a function f on the semi-infinite interval [a,+infinity)
    // it performs a change of variables and then integrates over the interval [0,1]
    // the interval is split into equal-sized subintervals and the integral is computed over each
    // subinterval, the amount of subintervals is determined by the subdivisions parameter
    // if error_estimate is not a NULL pointer, write the error_estimate to the passed variable
    //
    // default values: int subdivisions=5, double a = 0.0, double* error_estimate = nullptr
    //
    // usage: to integrate a function f(double x1, double x2, int x3)
    // over x1 from 0 to infinity (x2, x3 are fixed values defined somewhere in the code):
    //
    // auto integrand = [x2, x3](double x1) {return f(x1, x2, x3)};
    // double result = integrate_semi_inf(integrand);

    inline double integrate_interval(std::function<double (double) > f, double a, double b, double* error_estimate = nullptr)  {
    // integrate a function f on the finite interval [a, b], where b > a, using a 15-point Kronrod rule
    // if error_estimate is not a NULL pointer, write the error_estimate to the passed variable
    //
    // default values: double* error_estimate = nullptr
    //
    // usage: to integrate a function f(double x1, double x2, int x3)
    // over x1 from -1 to 1 (x2, x3 are fixed values defined somewhere in the code):
    //
    // auto integrand = [x2, x3](double x1) {return f(x1, x2, x3)};
    // double result = integrate_semi_interval(integrand, -1, 1);
    double result;
    double bma = (b - a) / 2., bpa = (a + b) / 2.;

    result = (f(-0.991455371120813 * bma + bpa) + f(0.991455371120813 * bma + bpa)) * 0.022935322010529;
    result += (f(-0.949107912342759 * bma + bpa) + f(0.949107912342759 * bma + bpa)) * 0.063092092629979;
    result += (f(-0.864864423359769 * bma + bpa) + f(0.864864423359769 * bma + bpa)) * 0.104790010322250;
    result += (f(-0.741531185599394 * bma + bpa) + f(0.741531185599394 * bma + bpa)) * 0.140653259715525;
    result += (f(-0.586087235467691 * bma + bpa) + f(0.586087235467691 * bma + bpa)) * 0.169004726639267;
    result += (f(-0.405845151377397 * bma + bpa) + f(0.405845151377397 * bma + bpa)) * 0.190350578064785;
    result += (f(-0.207784955007898 * bma + bpa) + f(0.207784955007898 * bma + bpa)) * 0.204432940075298;
    result += f(0. + bpa) * 0.209482141084728;

    result *= bma;

    if (error_estimate != nullptr) {
        *error_estimate = (f(-0.949107912342759 * bma + bpa) + f(0.949107912342759 * bma + bpa)) * 0.129484966168870;
        *error_estimate += (f(-0.741531185599394 * bma + bpa) + f(0.741531185599394 * bma + bpa)) * 0.279705391489277;
        *error_estimate += (f(-0.405845151377397 * bma + bpa) + f(0.405845151377397 * bma + bpa)) * 0.381830050505119;
        *error_estimate += f(0. + bpa) * 0.417959183673469;
        *error_estimate *= bma;
        *error_estimate = pow(200. * fabs(*error_estimate - result), 1.5);
    }
    return result;
}

#endif /* easy_numeric_hpp */
