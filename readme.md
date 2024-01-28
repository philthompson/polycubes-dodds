# Polycubes

This program computes the values for https://oeis.org/A000162 (except for trivially small values of N):

The number of unique 3D polycubes that can be constructed with *n* cubes.

This is a rust port of Stanley Dodds's algorithm, whose original C# code is available at: https://oeis.org/A000162/a000162.cs.txt

The "v1" code was used to compute the `n=21` count used to extend that OEIS sequence.

# Thread Tasks

By changing `FORK_DEPTH`, we can control how many thread tasks the work is split among.  Thread tasks will have varying running times, but in general, as we increase the number of tasks the average running time per task decreases.

Even with a million thread tasks, memory usage and initialization time aren't an issue on my M1 Mac mini with 16GB of RAM.

| FORK_DEPTH | thread tasks |
|   :---:    |         ---: |
|     2      |           15 |
|     3      |           86 |
|     4      |          534 |
|     5      |        3,481 |
|     6      |       23,502 |
|     7      |      162,913 |
|     8      |    1,152,870 |
|     9      |    8,294,738 |

# Performance

For v2, on an M1 Mac mini, with 8 threads:
```
n=11: 000d:00h:00m:00.217s (FORK_DEPTH=2)
n=12: 000d:00h:00m:00.423s (FORK_DEPTH=2)
n=13: 000d:00h:00m:01.320s (FORK_DEPTH=2)
n=14: 000d:00h:00m:02.717s (FORK_DEPTH=3)
n=15: 000d:00h:00m:10.714s (FORK_DEPTH=3)
n=16: 000d:00h:00m:52.436s (FORK_DEPTH=4)
n=17: 000d:00h:05m:33.623s (FORK_DEPTH=4)
n=18: 000d:00h:40m:17.341s (FORK_DEPTH=4)
```

For v1, on an M1 Mac mini, with 8 threads and `FILTER_DEPTH=5`:
```
n=11: 000d:00h:00m:01.509s
n=12: 000d:00h:00m:02.168s
n=13: 000d:00h:00m:04.851s
n=14: 000d:00h:00m:07.157s
n=15: 000d:00h:00m:16.229s
n=16: 000d:00h:01m:18.004s
n=17: 000d:00h:08m:28.415s
n=18: 000d:01h:04m:18.560s 
```

In AWS EC2, on a 16-thread machine, the "v1" code computed `n=21` over a period of 10 days.  For an analysis of running n=21 in AWS EC2, see my blog post here:
https://philthompson.me/2024/Counting-Polycubes-of-Size-21.html

# Running

```
cargo run --release --bin polycubes-v2

# optional, with a checkpoint that matches N:
cargo run --release --bin polycubes-v2 -- 17 5 10306999600203458141 6036447860508

# on an Apple Silicon Mac, the nightly aarch64 toolchain must be installed:
cargo +nightly-aarch64-apple-darwin run --release --bin polycubes-v2
```

# License

This code is released under the CC-BY-SA-4.0 license, as it is derived from Stanley Dodds's C# code published at the OEIS.
