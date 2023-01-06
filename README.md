# Cising

Cising is the C/C++ implementation of the Bornholdt Ising Model. All
background information can be found at the
[main page](https://github.com/kenokrieger/multising).

## Compiling

This is project was created using Jetbrain's CLion IDE. The easiest way to
compile it, is by loading it into CLion. Since it only consists of one file,
you may also try to compile it manually.

## Usage

The program expects a file "multising.conf" in the directory it is called from.
This file contains all the values for the simulation like grid size and parameters.
The path to the file can also be passed as the first argument in the terminal.

Example configuration file:

```
lattice_height = 2048
lattice_width = 2048
total_updates = 100000
seed = 2458242462
alpha = 15.0
j = 1.0
beta = 0.6
init_up = 0.5
```

For **alpha = 0** this model resembles the standard Ising model.

## License

This project is licensed under MIT License (see LICENSE.txt).
