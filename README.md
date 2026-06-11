# Gauss-Distribution-Trace-Generation
Discussing the possible to generate trace along a simple parametric curve.
---
The theoretical part is in `Trace.md` , as the realization details are in `AGENTS.md`.

## Build

```powershell
cmake -S . -B build-mingw -G "MinGW Makefiles"
cmake --build build-mingw
ctest --test-dir build-mingw --output-on-failure
```

## Run examples

Generate an ellipse trace:

```powershell
.\build-mingw\generate.exe --curve ellipse --steps 1000 --arc second --output outputs\ellipse_trace.csv
```

Generate a deterministic vertical-line trace:

```powershell
.\build-mingw\generate.exe --curve line --steps 300 --sigma-theta 0 --sigma-v 0 --sigma-measure 0 --output outputs\line_trace.csv
```

Use `--help` to see all configurable simulation parameters.
