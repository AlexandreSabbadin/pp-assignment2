# Parallel Programming - Assignment 2

https://github.com/AlexandreSabbadin/pp-assignment2

How to run the program:

```bash
module load openmpi
mpicc -O3 n_body.c -o n_body -lm
srun -n 2 n_body
```

To test the program I ran it with `-n` on *dione* cluster and with `-N` on *titan*. For each message size I did 5 tries then took the smallest result. All the results are compiled into the following table. Time unit is seconds.

| N body | -n  2     | -n 4      | -n 8      | -n 16    | -N 2      | -N 4      | -N 8     | -N 16    |
|--------|-----------|-----------|-----------|----------|-----------|-----------|----------|----------|
| 64     | 0.171029  | 0.122132  | 0.075690  | 0.084980 | 0.217245  | 0.186230  | 0.262533 | 0.512607 |
| 128    | 0.708361  | 0.376637  | 0.219380  | 0.159933 | 0.756690  | 0.454817  | 0.377408 | 0.524881 |
| 256    | 2.675324  | 1.442556  | 0.769425  | 0.445679 | 2.428269  | 1.542974  | 0.925659 | 0.864886 |
| 512    | 10.719141 | 5.796739  | 2.986953  | 1.573160 | 8.982571  | 4.956951  | 2.872476 | 1.975917 |
| 1024   | 42.511683 | 22.102144 | 11.756022 | 6.087977 | 35.497059 | 18.200705 | 9.625564 | 5.448056 |

Console logs for `srun -n 4 n_body`:

![nbody_out](https://user-images.githubusercontent.com/85174595/196035710-ee03a1b3-1606-46aa-b883-7c535ac57f92.PNG)

## How it works

**Initialization:** I initialize all local body arrays in process 0, then write positions in a file, then send them to all process. For this I made a `Body` *struct* and `MPI_Datatype` with all the data (position, velocity, forces). All values are randomly generated with a seed so they're the same for all runs.

**Computation:** I used the ring method to do all the computations in parallel. Each process send its local array to the next in a loop of all process then compute with the received previous one. I used non-blocking sending method `MPI_Isend`.

**Results:** All the final local body arrays are sent back to process 0 then the final positions are written in a file.
