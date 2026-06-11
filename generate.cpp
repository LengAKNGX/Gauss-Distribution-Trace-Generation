#include "necessary.h"

#include <exception>
#include <iostream>
#include <string>

namespace {

void print_usage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n"
        << "Options:\n"
        << "  --curve line|ellipse\n"
        << "  --steps N\n"
        << "  --dt value\n"
        << "  --L value\n"
        << "  --alpha value\n"
        << "  --beta value\n"
        << "  --v-exp value\n"
        << "  --sigma-theta value\n"
        << "  --sigma-v value\n"
        << "  --sigma-measure value\n"
        << "  --sample-N N\n"
        << "  --seed N\n"
        << "  --projection-delta value\n"
        << "  --newton-iters N\n"
        << "  --newton-tol value\n"
        << "  --arc exact|first|second|numerical\n"
        << "  --output path\n"
        << "  --a value --b value\n"
        << "  --x0 value --y0 value --theta0 value --v0 value --t0 value\n"
        << "  --allow-negative-speed\n"
        << "  --help\n";
}

bool has_value(int index, int argc, const std::string& option) {
    if (index + 1 >= argc) {
        throw std::invalid_argument("missing value for " + option);
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::string curve_type = "ellipse";
        double ellipse_a = 5.0;
        double ellipse_b = 2.0;

        gauss_trace::SimulationConfig config;
        config.output_csv_path = "outputs/ellipse_trace.csv";
        gauss_trace::State initial_state{{5.4, 0.0}, gauss_trace::pi / 2.0, 1.0, 0.0};

        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") {
                print_usage(argv[0]);
                return 0;
            } else if (arg == "--curve" && has_value(i, argc, arg)) {
                curve_type = argv[++i];
            } else if (arg == "--steps" && has_value(i, argc, arg)) {
                config.n_steps = std::stoi(argv[++i]);
            } else if (arg == "--dt" && has_value(i, argc, arg)) {
                config.dt = std::stod(argv[++i]);
            } else if (arg == "--L" && has_value(i, argc, arg)) {
                config.lookahead_L = std::stod(argv[++i]);
            } else if (arg == "--alpha" && has_value(i, argc, arg)) {
                config.alpha = std::stod(argv[++i]);
            } else if (arg == "--beta" && has_value(i, argc, arg)) {
                config.beta = std::stod(argv[++i]);
            } else if (arg == "--v-exp" && has_value(i, argc, arg)) {
                config.v_exp = std::stod(argv[++i]);
            } else if (arg == "--sigma-theta" && has_value(i, argc, arg)) {
                config.sigma_theta = std::stod(argv[++i]);
            } else if (arg == "--sigma-v" && has_value(i, argc, arg)) {
                config.sigma_v = std::stod(argv[++i]);
            } else if (arg == "--sigma-measure" && has_value(i, argc, arg)) {
                config.sigma_measure = std::stod(argv[++i]);
            } else if (arg == "--sample-N" && has_value(i, argc, arg)) {
                config.sample_interval_N = std::stoi(argv[++i]);
            } else if (arg == "--seed" && has_value(i, argc, arg)) {
                config.seed = std::stoull(argv[++i]);
            } else if (arg == "--projection-delta" && has_value(i, argc, arg)) {
                config.local_projection_delta = std::stod(argv[++i]);
            } else if (arg == "--newton-iters" && has_value(i, argc, arg)) {
                config.max_newton_iters = std::stoi(argv[++i]);
            } else if (arg == "--newton-tol" && has_value(i, argc, arg)) {
                config.newton_tol = std::stod(argv[++i]);
            } else if (arg == "--arc" && has_value(i, argc, arg)) {
                config.arc_length_method = gauss_trace::parse_arc_length_method(argv[++i]);
            } else if (arg == "--output" && has_value(i, argc, arg)) {
                config.output_csv_path = argv[++i];
            } else if (arg == "--a" && has_value(i, argc, arg)) {
                ellipse_a = std::stod(argv[++i]);
            } else if (arg == "--b" && has_value(i, argc, arg)) {
                ellipse_b = std::stod(argv[++i]);
            } else if (arg == "--x0" && has_value(i, argc, arg)) {
                initial_state.P.x = std::stod(argv[++i]);
            } else if (arg == "--y0" && has_value(i, argc, arg)) {
                initial_state.P.y = std::stod(argv[++i]);
            } else if (arg == "--theta0" && has_value(i, argc, arg)) {
                initial_state.theta = std::stod(argv[++i]);
            } else if (arg == "--v0" && has_value(i, argc, arg)) {
                initial_state.v = std::stod(argv[++i]);
            } else if (arg == "--t0" && has_value(i, argc, arg)) {
                initial_state.t = std::stod(argv[++i]);
            } else if (arg == "--allow-negative-speed") {
                config.speed_nonnegative = false;
            } else {
                throw std::invalid_argument("unknown option: " + arg);
            }
        }

        if (curve_type == "line" || curve_type == "vertical_line") {
            config.output_csv_path = config.output_csv_path == "outputs/ellipse_trace.csv"
                                         ? "outputs/line_trace.csv"
                                         : config.output_csv_path;
        }

        const auto curve = gauss_trace::make_curve(curve_type, ellipse_a, ellipse_b);
        gauss_trace::Simulation simulation(*curve, config, initial_state);
        const auto records = simulation.run();
        gauss_trace::write_csv_file(config.output_csv_path, *curve, config, initial_state, records);

        std::cout << "Wrote " << records.size() << " rows to " << config.output_csv_path << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
