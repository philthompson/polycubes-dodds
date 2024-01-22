# Polycubes

This program computes the values for https://oeis.org/A000162 (except for trivially small values of N):

The number of unique 3D polycubes that can be constructed with *n* cubes.

This is a rust port of Stanley Dodds's algorithm, whose original C# code is available at: https://oeis.org/A000162/a000162.cs.txt

This code was used to compute the `n=21` count used to extend that OEIS sequence.

# Performance

On an M1 Mac mini, with 8 threads and `FILTER_DEPTH=5`:
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

In AWS EC2, on a 16-thread machine, `n=21` was computed over a period of 10 days.  For an analysis of running n=21 in AWS EC2, see my blog post here:
https://philthompson.me/2024/Counting-Polycubes-of-Size-21.html

# Running

```
cargo run --release

# on an Apple Silicon Mac, the nightly aarch64 toolchain must be installed:
cargo +nightly-aarch64-apple-darwin run --release
```

# License

This code is released under the CC-BY-SA-4.0 license, as it is derived from Stanley Dodds's C# code published at the OEIS.
