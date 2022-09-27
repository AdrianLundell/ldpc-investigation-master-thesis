# Using AFF3CT as a library

This part aims at building onto AFF3CT in order to make it compatible with simulations of NAND Flash Memories and assumes that a Linux system is used.

Get the AFF3CT library:

	$ git submodule update --init --recursive

Compile the library on Linux/MacOS/MinGW:

	$ cd lib/aff3ct
	$ mkdir build
	$ cd build
	$ cmake .. -G"Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-funroll-loops -march=native" -DAFF3CT_COMPILE_EXE="OFF" -DAFF3CT_COMPILE_STATIC_LIB="ON" -DAFF3CT_COMPILE_SHARED_LIB="ON"
	$ make -j4

Now the AFF3CT library has been built in the `lib/aff3ct/build` folder.

The next step is to build the additional files that enable NAND Flash simulations. In order to do so, copy the cmake configuration files from the AFF3CT build:

	$ cd ../../../
	$ mkdir cmake mkdir cmake/Modules
	$ cp lib/aff3ct/build/lib/cmake/aff3ct-*/* cmake/Modules
	
Next, build the code:

	$ mkdir build
	$ cd build
	$ cmake .. -G"Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-funroll-loops -march=native" -DAFF3CT_COMPILE_EXE="ON" -DAFF3CT_COMPILE_TESTS="OFF"
	$ make


Run the simulations through:

	$ ./bin/aff3ct_nand_flash

If you have created some tests, you can choose to compile those instead:

	$ cmake .. -G"Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-funroll-loops -march=native" -DAFF3CT_COMPILE_EXE="OFF" -DAFF3CT_COMPILE_TESTS="ON"
	$ make
	
Run the tests through:

	$ ./tests/aff3ct_nand_flash


# Make adjustments to simulations

The simulations can be modified in the [src/main.cpp](src/main.cpp). Please refer to [aff3ct documentation](https://aff3ct.readthedocs.io/en/latest/) for more information on how to make adjustments.
