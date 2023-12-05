# Exercise 1
This folder contains code and data for the first assignment. Further info can be found [here](https://github.com/Foundations-of-HPC/Foundations_of_HPC_2022/blob/main/Assignment/exercise1/Assignment_exercise1.pdf).\
The goal of the exercise is to implement a parallel version of a variant of the [Conway’s “Game of Life”](https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life).\
Different evolution methods are implemented:
- Static Evolution
- Ordered Evolution
- Wave Evolution

## Structure of the repository
- `src` folder: here you can find the code for the different evolution methods we implemented. Additionally functions for the initialization of the playground are provided:
  - [Static evolution](./src/static_update.c)
  - [Ordered evolution](./src/ordered_update.c)
  - [Wave evolution](./src/wave_update.c)
  - [Initializing, reading and writing in parallel](./src/io_init.c)

- [`data`](./data/) folder: here are stored all the .csv files with the times obtained by running the code. The folder is subdivided according to the evolution methods used.
- [`sbatch_files`](./sbatch_files/) folder: here are stored the .sh files used to run and time the code on the cluster.
- [`files`](./files/) folder: here are stored the .pgm files with the snapshots of the playgrounds. The folder is subdivided according to the evolution methods used, plus a folder for the initial playgrounds: `init`.
- [`old_versions`](./old_versions/) folder: here are stored various attempts in tackling the problem that didn't make it to the final implementation.
## How to compile and run the code
To compile the code simply run:
```
make
```
OpenMP and OpenMPI are required to compile the code.\
If you are using the ORFEO cluster, you can load the modules with:
```
module load openMPI/4.1.5/gnu/12.2.1
```

To run the code use `mpirun` with the following arguments:
- -np *int value*: number of processes
- -i: initialize a playground
- -r: run a playground
- -k *value*: playground size
- -e [0|1|2]: type of evoultion; 0 means "ordered", 1 means "static", 2 means "wave"
- -f *string*: if used with -i, it is the name of the file where the initial playground is saved. If used with -r, it is the name of the file where the initial playground is read from.
- -n *int value*: number of steps to be calculated
- -s *int value*:every how many steps a dump of the system is saved on a file (0 meaning only at the end). This will be saved in the folder `files/` inside with name of the form *snapshot_nnnnn.pgm*, where *nnnnn* is the step at which the image is captured padded with zeros to obtain 5 digits.
