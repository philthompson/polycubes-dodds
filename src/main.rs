/*

this program computes the values for https://oeis.org/A000162
  (except for trivially small values of N)

this is a rust port of Stanley Dodds's algorithm, whose original
  C# code is available at:
  https://oeis.org/A000162/a000162.cs.txt

rust port by Phil Thompson in January 2024, along with changes:
- use multiple threads for computation of "nontrivial symmetries" portion
- allow use of previously-computed counts (to resume a long-running computation)
- allow specifying more threads than worker tasks
- allow reversing "trivial symmetry" worker tasks to reduce idle CPU time with a
    high number of threads
- print timestamps (elapsed time) for most messages

-----

as it stands now, the maximum number of threads needed (or AWS EC2 vCPUs) is 32,
  even for n=22, because the slowest tasks run in ~5x the time as the fastest

by starting tasks in reverse order (starting at filter=0, using RUN_JOBS_REVERSED=true)
   running with half the threads (32) will complete in the same overall runtime as
   starting with the highest filter number with 64 threads

it is possible to use even fewer threads (16) to have a somewhat longer running time
   but at lower overall AWS EC2 cost -- see this blog post for details:
   https://philthompson.me/2024/Counting-Polycubes-of-Size-21.html

for n=20, by FILTER_DEPTH:
5 -> 58 + 1 worker tasks (including 0th task)
6 -> 54 + 1 worker tasks (including 0th task)
7 -> 50 + 1 worker tasks (including 0th task)

for n=21, by FILTER_DEPTH:
5 -> 62 + 1 worker tasks (including 0th task)
6 -> 58 + 1 worker tasks (including 0th task)
7 -> 54 + 1 worker tasks (including 0th task)

for n=22, by FILTER_DEPTH:
5 -> 66 + 1 worker tasks (including 0th task)
6 -> 62 + 1 worker tasks (including 0th task)
7 -> 58 + 1 worker tasks (including 0th task)
*/

use std::io::Write;
use std::time::Duration;
use std::time::Instant;
use std::thread;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
// this import needed depending on rust toolchain (e.g. on Amazon Linux)
//use std::convert::TryInto;

const N: usize = 15; // number of polycube cells. Need n >= 4 if single threading, or n >= filterDepth >= 5 if multithreading (I think)
const FILTER_DEPTH: usize = 5; // keep this at 5 to divide the work among the most worker tasks
const THREADS: usize = 4;
const USE_PRECOMPUTED_SYMM: bool = true; // use precomputed nontrivial symmetries, if available

// the first portion of the program computes "nontrivial symmetries" which is much
//   faster than the second portion of the program.  if running a high n value, or
//   piping output to a file, then set this to false
const SHOW_NONTRIVIAL_SYMM_ETA: bool = true;
// set RUN_JOBS_REVERSED=false to start the fastest "trivial symmetry" worker tasks first
// set RUN_JOBS_REVERSED=true to start the slowest tasks first
// when reversed, this will reduce overall idle CPU time (and thus on-demand compute
//   cost) if THREADS is ~ 1/3 - 1/2 the number of "trivial symmetry" worker tasks
// for example, for n=18 or n=22, 64 threads with RUN_JOBS_REVERSED=false will
//   run in the same overall time as a 32 threads with RUN_JOBS_REVERSED=true
// (see above note and blog post link)
const RUN_JOBS_REVERSED: bool = true;

const X: i32 = (N as i32 + 5) / 4 * ((N as i32 + 5) / 4 * 3 - 2); // set X<Y<Z such that aX+bY+cZ = 0 implies a = b = c = 0 or |a|+|b|+|c| > n
const Y: i32 = X + 1; // trivial choice is X = 1, Y = n, Z = n * n. A simple reduction is X = 1, Y = n, Z = n * (n / 2) + (n + 1) / 2
const Z: i32 = X + (N as i32 + 5) / 4 * 3; // minimising Z is memory efficient. Unclear if this noticably affects performance. Z ~ 3/16 n^2 is the best I can find for arbitrary n

const X2: usize = (X + X) as usize;
const Y2: usize = (Y + Y) as usize;
const Z2: usize = (Z + Z) as usize;
const SYX: usize = (Y + X) as usize;
const SZX: usize = (Z + X) as usize;
const SZY: usize = (Z + Y) as usize;
const DYX: usize = (Y - X) as usize;
const DZX: usize = (Z - X) as usize;
const DZY: usize = (Z - Y) as usize;

const DESCRIPTIONS: [&str; 4] = [
	"orthogonal order 2 rotation",
	"orthogonal order 4 rotation",
	"short diagonal order 2 rotation",
	"long diagonal order 3 rotation"
];
const AUT_CLASS_SIZES: [usize; 4] = [3, 6, 6, 8];
const MATRIX_REPS: [[i32; 9]; 4] = [
	[-1, 0, 0, 0, -1, 0, 0, 0, 1],
	[0, -1, 0, 1, 0, 0, 0, 0, 1],
	[0, 1, 0, 1, 0, 0, 0, 0, -1],
	[0, 0, 1, 1, 0, 0, 0, 1, 0]
];
const AFFINE1: [[i32; 3]; 4] = [
	[1, 0, 0],
	[1, 0, 0],
	[1, -1, 0],
	[1, 0, -1]
];
const AFFINE2: [[i32; 3]; 4] = [
	[0, 1, 0],
	[0, 1, 0],
	[0, 0, 1],
	[0, 1, -1]
];
const BIASES: [[i32; 3]; 4] = [
	[(2 * N) as i32, (2 * N) as i32, 0],
	[(2 * N) as i32, 0, 0],
	[0, 0, 2],
	[(N - 1) as i32, 0, 1 - (N as i32)]
];

// with the USE_PRECOMPUTED_SYMM constant:
// lets us optionally use previously-calculated values here
const NONTRIVIAL_SYMMETRIES_COUNTS: [usize; 23] = [0,
	/* n=1  */             0,
	/* n=2  */             0,
	/* n=3  */             0,
	/* n=4  */             0,
	/* n=5  */             0,
	/* n=6  */             0,
	/* n=7  */             0,
	/* n=8  */             0,
	/* n=9  */             0,
	/* n=10 */             0,
	/* n=11 */        45_979,
	/* n=12 */             0,
	/* n=13 */             0,
	/* n=14 */     1_138_017,
	/* n=15 */     2_446_834,
	/* n=16 */     8_460_765,
	/* n=17 */    18_283_068,
	/* n=18 */    63_525_537,
	/* n=19 */             0,
	/* n=20 */             0,
	/* n=21 */ 1_055_564_170,
	/* n=22 */ 3_699_765_374];

// for long-running counts, we can record previously-computed worker
//   task counts here
// this allows us to resume cancelled runs while keeping everything
//   contained in this single file
// === comment these out to re-compute them! ===
const TRIVIAL_SYMMETRIES_COUNTS_BY_N_FILTER_DEPTH_FILTER: [(usize, usize, usize, usize); 67] = [
	//N, FILTER_DEPTH, filter, count
	// N=17
	(17, 5, 4, 289836494778),
	(17, 5, 8, 583673049457),
	(17, 5, 9, 646961888347),
	(17, 5, 22, 223300442187),
	// N=21
	(21, 5, 0, 309420173722300),   // took 2d:17h:12m on an AWS EC2 c7g instance
	(21, 5, 1, 345395528117841),   // took 3d:00h:52m on an AWS EC2 c7g instance
	(21, 5, 2, 512845497549304),   // took 3d:16h:02m on an AWS EC2 c7g instance
	(21, 5, 3, 696292340206626),   // took 4d:07h:12m on an AWS EC2 c7g instance
	(21, 5, 4, 863102288279808),   // took 4d:19h:21m on an AWS EC2 c7g instance
	(21, 5, 5, 1153440327911952),  // took 5d:12h:03m on an AWS EC2 c7g instance
	(21, 5, 6, 1387546375414082),  // took 5d:22h:44m on an AWS EC2 c7g instance
	(21, 5, 7, 1605849914003528),  // took 6d:05h:44m on an AWS EC2 c7g instance
	(21, 5, 8, 1787870546368474),  // took 6d:08h:08m on an AWS EC2 c7g instance
	(21, 5, 9, 1997960513986194),  // took 6d:09h:40m on an AWS EC2 c7g instance
	(21, 5, 10, 2119506128616415), // took 6d:06h:47m on an AWS EC2 c7g instance
	(21, 5, 11, 2202018997735872), // took 6d:01h:04m on an AWS EC2 c7g instance
	(21, 5, 12, 2244177905798528), // took 5d:18h:05m on an AWS EC2 c7g instance
	(21, 5, 13, 2248001799465557), // took 5d:09h:53m on an AWS EC2 c7g instance
	(21, 5, 14, 2198464880362828), // took 5d:00h:26m on an AWS EC2 c7g instance
	(21, 5, 15, 2125278404299050), // took 4d:15h:04m on an AWS EC2 c7g instance
	(21, 5, 16, 2010116122002288), // took 4d:05h:33m on an AWS EC2 c7g instance
	(21, 5, 17, 1876807006549215), // took 3d:20h:29m on an AWS EC2 c7g instance
	(21, 5, 18, 1725640003589544), // took 3d:12h:03m on an AWS EC2 c7g instance
	(21, 5, 19, 1562659895679425), // took 3d:04h:18m on an AWS EC2 c7g instance
	(21, 5, 20, 1393374663858664), // took 2d:21h:18m on an AWS EC2 c7g instance
	(21, 5, 21, 1226996979655412), // took 2d:15h:09m on an AWS EC2 c7g instance
	(21, 5, 22, 1063878243997825), // took 2d:09h:46m on an AWS EC2 c7g instance
	(21, 5, 23, 909715901464397),  // took 2d:05h:09m on an AWS EC2 c7g instance
	(21, 5, 24, 767003283534118),  // took 2d:01h:14m on an AWS EC2 c7g instance
	(21, 5, 25, 637664564613305),  // took 1d:21h:57m on an AWS EC2 c7g instance
	(21, 5, 26, 522253112618485),  // took 1d:19h:16m on an AWS EC2 c7g instance
	(21, 5, 27, 421693153196302),  // took 1d:17h:04m on an AWS EC2 c7g instance
	(21, 5, 28, 335497140472044),  // took 1d:15h:19m on an AWS EC2 c7g instance
	(21, 5, 29, 263031935220780),  // took 1d:13h:57m on an AWS EC2 c7g instance
	(21, 5, 30, 203070780532962),  // took 1d:12h:51m on an AWS EC2 c7g instance
	(21, 5, 31, 154382100269004),  // took 1d:12h:01m on an AWS EC2 c7g instance
	(21, 5, 32, 115559882167060),  // took 1d:11h:24m on an AWS EC2 c7g instance
	(21, 5, 33, 85152304485303),   // took 1d:10h:56m on an AWS EC2 c7g instance
	(21, 5, 34, 61710362959494),   // took 1d:10h:35m on an AWS EC2 c7g instance
	(21, 5, 35, 43982755199356),   // took 1d:10h:21m on an AWS EC2 c7g instance
	(21, 5, 36, 30822164177826),   // took 1d:10h:09m on an AWS EC2 c7g instance
	(21, 5, 37, 21225762905510),   // took 1d:10h:02m on an AWS EC2 c7g instance
	(21, 5, 38, 14350322250780),   // took 1d:09h:56m on an AWS EC2 c7g instance
	(21, 5, 39, 9520721822656),    // took 1d:09h:54m on an AWS EC2 c7g instance (0d:21h:55m on an M1 Mac mini)
	(21, 5, 40, 6196997682686),    // took 1d:09h:52m on an AWS EC2 c7g instance (0d:21h:53m on an M1 Mac mini)
	(21, 5, 41, 3955593817336),    // took 1d:09h:49m on an AWS EC2 c7g instance (0d:21h:52m on an M1 Mac mini)
	(21, 5, 42, 2471981286083),    // took 1d:09h:50m on an AWS EC2 c7g instance (0d:21h:51m on an M1 Mac mini)
	(21, 5, 43, 1511722821800),    // took 1d:09h:48m on an AWS EC2 c7g instance (0d:22h:30m on an M1 Mac mini)
	(21, 5, 44, 904967465059),     // took 1d:09h:49m on an AWS EC2 c7g instance (0d:22h:29m on an M1 Mac mini)
	(21, 5, 45, 529652026618),     // took 1d:09h:47m on an AWS EC2 c7g instance (0d:22h:29m on an M1 Mac mini)
	(21, 5, 46, 303040281933),     // took 1d:09h:48m on an AWS EC2 c7g instance (0d:22h:29m on an M1 Mac mini)
	(21, 5, 47, 168852583688),     // took 1d:09h:47m on an AWS EC2 c7g instance (0d:21h:50m on an M1 Mac mini)
	(21, 5, 48, 91752335576),      // took 1d:09h:47m on an AWS EC2 c7g instance (0d:21h:50m on an M1 Mac mini)
	(21, 5, 49, 48862051752),      // took 1d:09h:47m on an AWS EC2 c7g instance (0d:21h:50m on an M1 Mac mini)
	(21, 5, 50, 25216528208),      // took 1d:09h:47m on an AWS EC2 c7g instance (0d:21h:50m on an M1 Mac mini)
	(21, 5, 51, 12578787374),      // took 1d:09h:48m on an AWS EC2 c7g instance (0d:21h:51m on an M1 Mac mini)
	(21, 5, 52, 6221625934),       // took 1d:09h:48m on an AWS EC2 c7g instance (0d:21h:51m on an M1 Mac mini)
	(21, 5, 53, 2898865668),       // took 1d:09h:47m on an AWS EC2 c7g instance (0d:21h:51m on an M1 Mac mini)
	(21, 5, 54, 1307081085),       // took 1d:09h:47m on an AWS EC2 c7g instance (0d:21h:51m on an M1 Mac mini)
	(21, 5, 55, 623956362),        // took 1d:09h:49m on an AWS EC2 c7g instance (0d:22h:09m on an M1 Mac mini)
	(21, 5, 56, 254705772),        // took 1d:09h:49m on an AWS EC2 c7g instance (0d:22h:09m on an M1 Mac mini)
	(21, 5, 57, 83869444),         // took 1d:09h:46m on an AWS EC2 c7g instance (0d:22h:09m on an M1 Mac mini)
	(21, 5, 58, 46709392),         // took 1d:09h:48m on an AWS EC2 c7g instance (0d:22h:09m on an M1 Mac mini)
	(21, 5, 59, 19610164),         // took 1d:09h:47m on an AWS EC2 c7g instance (0d:21h:52m on an M1 Mac mini)
	(21, 5, 60, 2192944),          // took 1d:09h:49m on an AWS EC2 c7g instance (0d:21h:52m on an M1 Mac mini)
	(21, 5, 61, 1157758),          // took 1d:09h:48m on an AWS EC2 c7g instance (0d:21h:52m on an M1 Mac mini)
	(21, 5, 62, 1291256),          // took 1d:09h:48m on an AWS EC2 c7g instance (0d:21h:52m on an M1 Mac mini)
];

// for smaller values of N, we can greatly speed up computation by
//   lowering the sleep time used when waiting for threads to finish
// for larger values of N, we can check less often and waste fewer
//   cycles in the main thread
const N_SLEEP_MILLIS: u64 = if N > 19 { N as u64 * 100 } else if N > 17 { N as u64 * 25 } else if N > 12 { N as u64 * 2 } else { N as u64 + 5 };

fn seconds_to_dur(s: f64) -> String {
	let days = (s / 86400.0).floor();
	let hours = ((s - (days * 86400.0)) / 3600.0).floor();
	let minutes = ((s - (days * 86400.0) - (hours * 3600.0)) / 60.0).floor();
	let seconds = s - (days * 86400.0) - (hours * 3600.0) - (minutes * 60.0);
	return format!("{:0>3}d:{:0>2}h:{:0>2}m:{:0>6.3}s", days, hours, minutes, seconds);
}

fn print_w_time(start_time: Instant, message: String) {
	println!("[{}] {message}", seconds_to_dur(start_time.elapsed().as_secs_f64()));
}

// thanks to https://stackoverflow.com/a/67834588/259456
fn integer_with_thousands_separator(num: usize) -> String {
	return num.to_string()
    .as_bytes()
    .rchunks(3)
    .rev()
    .map(std::str::from_utf8)
    .collect::<Result<Vec<&str>, _>>()
    .unwrap()
    .join(",");  // separator
}

fn main() {
	let overall_start_time = Instant::now();
	let symmetry_time_elapsed: f64;
	let mut count: usize = 0; // enumerate the sum over the order 24 group of the size of the fix of each group element, and divide by 24 (Burnside's lemma)

	if USE_PRECOMPUTED_SYMM && NONTRIVIAL_SYMMETRIES_COUNTS[N] > 0 {
		count += NONTRIVIAL_SYMMETRIES_COUNTS[N];
		symmetry_time_elapsed = overall_start_time.elapsed().as_secs_f64();
		println!("used precomputed nontrivial symmetries count:");
		println!("total count for nontrivial symmetries is {} for polycubes with {} cells\nTook {}", count, N, seconds_to_dur(symmetry_time_elapsed));
	} else if THREADS == 1 {
		for sym in 0..4 {
			let mut subcount = 0;
			for i in (1 - (N as i32))..(N as i32) {
				for j in (1 - (N as i32) + i.abs())..=((N as i32) - 1 - i.abs()) {
					subcount += count_symmetric_polycubes(
						&MATRIX_REPS[sym],
						[
							i * AFFINE1[sym][0] + j * AFFINE2[sym][0] + BIASES[sym][0],
							i * AFFINE1[sym][1] + j * AFFINE2[sym][1] + BIASES[sym][1],
							i * AFFINE1[sym][2] + j * AFFINE2[sym][2] + BIASES[sym][2]
						]);
				}
			}
			println!("\n{} polycubes fixed under each {}", subcount, DESCRIPTIONS[sym]);
			println!("{} in total (symmetry occurs {} times)", AUT_CLASS_SIZES[sym] * subcount, AUT_CLASS_SIZES[sym]);
			count += AUT_CLASS_SIZES[sym] * subcount;
		}
		symmetry_time_elapsed = overall_start_time.elapsed().as_secs_f64();
		println!("total count for nontrivial symmetries is {} for polycubes with {} cells\nTook {}", count, N, seconds_to_dur(symmetry_time_elapsed));
	} else {
		let n_sleep = Duration::from_millis(N_SLEEP_MILLIS);
		println!("N={N}, n_sleep_millis={N_SLEEP_MILLIS}");
		// run a batch of worker tasks for each of the 4 symmetries
		for sym in 0..4 {
			print_w_time(overall_start_time, format!("Counting {} symmetries (set {} of 4):", DESCRIPTIONS[sym], sym+1));
			if SHOW_NONTRIVIAL_SYMM_ETA {
				println!("(note: percentage complete and ETA is just for this symmetry, not the entire overall count)");
			}
			let sw: Instant = Instant::now();
			let mut subcount = 0;

			let mut task_params: Vec<(i32, i32)> = Vec::new();
			let mut tasks = Vec::new();
			let mut total_worker_jobs: f64 = 0.0;

			for i in (1 - (N as i32))..(N as i32) {
				for j in (1 - (N as i32) + i.abs())..=((N as i32) - 1 - i.abs()) {
					task_params.push((i, j));
					total_worker_jobs += 1.0;
				}
			}
			for _ in 0..THREADS {
				if task_params.is_empty() {
					break
				}
				let (i, j) = task_params.pop().unwrap();
				let handle = thread::spawn(move || count_symmetric_polycubes_task(i, j, sym));
				tasks.push(handle);
			}

			let mut handle_index_to_join: Option<usize> = None;
			let mut last_stats = Instant::now();
			let mut compl_worker_jobs: f64 = 0.0;
			while !tasks.is_empty() || !task_params.is_empty() {
				let mut j: usize = 0;
				while j < tasks.len() {
					if tasks[j].is_finished() {
						handle_index_to_join = Some(j);
						compl_worker_jobs += 1.0;
						break;
					}
					j += 1;
				}
				if SHOW_NONTRIVIAL_SYMM_ETA && last_stats.elapsed().as_secs_f32() > 1.0 && compl_worker_jobs > 0.0 {
					last_stats = Instant::now();
					let time_elapsed = sw.elapsed();
					let seconds_per_sym = time_elapsed.as_secs_f64() / compl_worker_jobs;
					let sym_remaining = total_worker_jobs - compl_worker_jobs;
					let seconds_remaining = sym_remaining * seconds_per_sym;
					let pct_complete = (compl_worker_jobs * 100.0) / total_worker_jobs;
					print!("    {:.4}% complete, ETA:[{}], counting for n={}, sym={}, outstanding sym:[{}-{}={}]        \r",
						pct_complete,
						seconds_to_dur(seconds_remaining),
						N,
						sym,
						total_worker_jobs as usize,
						compl_worker_jobs as usize,
						(total_worker_jobs - compl_worker_jobs) as usize);
					std::io::stdout().flush().unwrap();
				}
				if handle_index_to_join.is_none() {
					thread::sleep(n_sleep);
					continue;
				}
				// .take() the handle to join and remove it from the tasks list
				subcount += tasks.remove(handle_index_to_join.take().unwrap()).join().unwrap();
				// replace the completed thread with a new one
				if let Some((i, j)) = task_params.pop() {
					let handle = thread::spawn(move || count_symmetric_polycubes_task(i, j, sym));
					tasks.push(handle);
				}
			}

			println!("\n{} polycubes fixed under each {}", subcount, DESCRIPTIONS[sym]);
			println!("{} in total (symmetry occurs {} times)\n", AUT_CLASS_SIZES[sym] * subcount, AUT_CLASS_SIZES[sym]);
			count += AUT_CLASS_SIZES[sym] * subcount;
		}

		symmetry_time_elapsed = overall_start_time.elapsed().as_secs_f64();
		println!("\ntotal count for nontrivial symmetries is {} for polycubes with {} cells\nTook {}", count, N, seconds_to_dur(symmetry_time_elapsed));
	}

	{
		print_w_time(overall_start_time, format!("\nCounting polycubes via extensions..."));
		let n_sleep = Duration::from_millis(N_SLEEP_MILLIS);
		println!("N={N}, n_sleep_millis={N_SLEEP_MILLIS}");
		let sw = Instant::now();
		let mut subcount: usize = 0;
		let mut tasks = Vec::new(); // please don't judge my multithreading - this was basically just my first attempt, to see if it helped
		let mut jobs_remaining: isize = (4 * (N as i32 - FILTER_DEPTH as i32) - 2).try_into().unwrap(); // this is the maximum left stack length, which is the value being filtered to separate work
		let total_tasks: usize = (jobs_remaining as usize) + 1;
		let mut competed_tasks: usize = 0;
		// see if any previously-recorded jobs' counts are available
		let mut recorded_counts: BTreeMap<usize, usize> = BTreeMap::new();
		for recorded_count in TRIVIAL_SYMMETRIES_COUNTS_BY_N_FILTER_DEPTH_FILTER {
			// for each tuple (N, FILTER_DEPTH, filter) keep the previously-computed count
			if recorded_count.0 == N && recorded_count.1 == FILTER_DEPTH {
				recorded_counts.insert(recorded_count.2, recorded_count.3);
			}
		}
		while tasks.len() < THREADS {
			let filter = if RUN_JOBS_REVERSED { total_tasks as isize - jobs_remaining - 1 } else { jobs_remaining }; // copy of i, since lambda expression captures the variable
			jobs_remaining -= 1;
			match recorded_counts.get(&(filter as usize)) {
				Some(prev_count) => {
					println!("using previously-recorded count [{prev_count}] for (N={N}, FILTER_DEPTH={FILTER_DEPTH}, filter={filter})");
					subcount += prev_count;
					competed_tasks += 1;
				},
				None => {
					//println!("no prev computed count for filter={filter}");
					let handle = thread::spawn(move || count_extensions_subset(filter as usize, overall_start_time));
					tasks.push(handle);
				}
			}
			// if we have more threads than jobs, stop here (jobs_remaining needs to
			//   get to -1 before we stop)
			if jobs_remaining < 0 {
				break;
			}
		}
		let mut handle_index_to_replace: Option<usize> = None;
		// this needs to continue once more when jobs_remaining == 0
		while jobs_remaining >= 0 || !tasks.is_empty() {
			let mut j: usize = 0;
			while j < tasks.len() {
				if tasks[j].is_finished() {
					handle_index_to_replace = Some(j);
					break;
				}
				j += 1;
			}
			if handle_index_to_replace.is_none() {
				thread::sleep(n_sleep);
				continue;
			}
			subcount += tasks.remove(handle_index_to_replace.take().unwrap()).join().unwrap();
			competed_tasks += 1;
			print_w_time(overall_start_time, format!("tasks remaining: [{total_tasks}-{competed_tasks}={}]", total_tasks-competed_tasks));
			// this needs to continue once more when jobs_remaining == 0
			while jobs_remaining >= 0 {
				let filter = if RUN_JOBS_REVERSED { total_tasks as isize - jobs_remaining - 1 } else { jobs_remaining }; // copy of i, since lambda expression captures the variable
				jobs_remaining -= 1;

				match recorded_counts.get(&(filter as usize)) {
					Some(prev_count) => {
						print_w_time(overall_start_time, format!("using previously-recorded count [{prev_count}] for (N={N}, FILTER_DEPTH={FILTER_DEPTH}, filter={filter})"));
						subcount += prev_count;
						competed_tasks += 1;
					},
					None => {
						//print_w_time(overall_start_time, format!("no prev computed count for filter={filter}"));
						let handle = thread::spawn(move || count_extensions_subset(filter as usize, overall_start_time));
						tasks.push(handle);
						break; // stop looping once we spawn a thread to replace the completed one
					}
				}
			}
		}
		count += subcount;
		println!("\n\n{:} polycubes with {} cells (number of polycubes fixed by trivial symmetry)\nTook {}", subcount, N, seconds_to_dur(sw.elapsed().as_secs_f64()));
	}
	println!();
	count /= 24;
	println!("{} free polycubes with {} cells", integer_with_thousands_separator(count), N);
	println!("overall time: {}", seconds_to_dur(overall_start_time.elapsed().as_secs_f64()));
}

fn count_symmetric_polycubes_task(i: i32, j: i32, sym: usize) -> usize {
	return count_symmetric_polycubes(
		&MATRIX_REPS[sym],
		[
			i * AFFINE1[sym][0] + j * AFFINE2[sym][0] + BIASES[sym][0],
			i * AFFINE1[sym][1] + j * AFFINE2[sym][1] + BIASES[sym][1],
			i * AFFINE1[sym][2] + j * AFFINE2[sym][2] + BIASES[sym][2]
		]);
}

fn count_symmetric_polycubes(linear_map: &[i32; 9], affine_shift: [i32; 3]) -> usize {
	let mut adjacency_counts = vec![vec![vec![0u8; (N + 2) as usize]; (2 * N + 1) as usize]; (2 * N + 1) as usize];
	'zloop: for z in 0..2 {
		for y in 0..(2 * N + 1) {
			for x in 0..(2 * N + 1) {
				adjacency_counts[x as usize][y as usize][z as usize] = 1;
				if z == 1 && y == N && x == N {
					break 'zloop;
				}
			}
		}
	}
	let mut required_cells = BTreeSet::new();
	let mut extension_stack: Vec<(i32, i32, i32)> = Vec::new();
	extension_stack.push((N as i32, N as i32, 1));
	let mut recovery_stack: Vec<(i32, i32, i32)> = Vec::new();
	return count_extensions(N as i32, &mut adjacency_counts, &mut required_cells, &mut extension_stack, &mut recovery_stack, linear_map, affine_shift);
}

fn count_extensions(
		mut cells_to_add: i32,
		adjacency_counts: &mut Vec<Vec<Vec<u8>>>,
		required_cells: &mut BTreeSet<(i32, i32, i32)>,
		extension_stack: &mut Vec<(i32, i32, i32)>,
		recovery_stack: &mut Vec<(i32, i32, i32)>,
		linear_map: &[i32; 9],
		affine_shift: [i32; 3]) -> usize {
	cells_to_add -= 1;
	let mut count = 0;
	let original_length = extension_stack.len();
	while let Some((x, y, z)) = extension_stack.pop() {
		recovery_stack.push((x, y, z));

		let existing_requirement = required_cells.remove(&(x, y, z));
		if !existing_requirement {
			if cells_to_add < required_cells.len() as i32 {
				continue; // number of req cells will only grow...
			}
			let mut temp_x = x;
			let mut temp_y = y;
			let mut temp_z = z;
			loop // works for general transformations of finite order
			{
				let (tx, ty, tz) = (
					linear_map[0] * temp_x + linear_map[1] * temp_y + linear_map[2] * temp_z + affine_shift[0],
					linear_map[3] * temp_x + linear_map[4] * temp_y + linear_map[5] * temp_z + affine_shift[1],
					linear_map[6] * temp_x + linear_map[7] * temp_y + linear_map[8] * temp_z + affine_shift[2],
				);
				temp_x = tx;
				temp_y = ty;
				temp_z = tz;
				if x == temp_x && y == temp_y && z == temp_z {
					break;
				}
				required_cells.insert((temp_x, temp_y, temp_z));
			}
		}

		if cells_to_add >= required_cells.len() as i32 {
			// if there are too many req cells...
			if cells_to_add == 0 {
				count += 1;
			} else {
				let inner_original_length = extension_stack.len();

				let adj_count: &mut u8 = &mut adjacency_counts[(x - 1) as usize][y as usize][z as usize];
				if *adj_count == 0 {
					extension_stack.push((x - 1, y, z));
				}
				*adj_count += 1;
				let adj_count: &mut u8 = &mut adjacency_counts[x as usize][(y - 1) as usize][z as usize];
				if *adj_count == 0 {
					extension_stack.push((x, y - 1, z));
				}
				*adj_count += 1;
				let adj_count: &mut u8 = &mut adjacency_counts[x as usize][y as usize][(z - 1) as usize];
				if *adj_count == 0 {
					extension_stack.push((x, y, z - 1));
				}
				*adj_count += 1;
				let adj_count: &mut u8 = &mut adjacency_counts[(x + 1) as usize][y as usize][z as usize];
				if *adj_count == 0 {
					extension_stack.push((x + 1, y, z));
				}
				*adj_count += 1;
				let adj_count: &mut u8 = &mut adjacency_counts[x as usize][(y + 1) as usize][z as usize];
				if *adj_count == 0 {
					extension_stack.push((x, y + 1, z));
				}
				*adj_count += 1;
				let adj_count: &mut u8 = &mut adjacency_counts[x as usize][y as usize][(z + 1) as usize];
				if *adj_count == 0 {
					extension_stack.push((x, y, z + 1));
				}
				*adj_count += 1;

				count += count_extensions(
					cells_to_add,
					adjacency_counts,
					required_cells,
					extension_stack,
					recovery_stack,
					linear_map,
					affine_shift,
				);
				adjacency_counts[(x - 1) as usize][y as usize][z as usize] -= 1;
				adjacency_counts[x as usize][(y - 1) as usize][z as usize] -= 1;
				adjacency_counts[x as usize][y as usize][(z - 1) as usize] -= 1;
				adjacency_counts[(x + 1) as usize][y as usize][z as usize] -= 1;
				adjacency_counts[x as usize][(y + 1) as usize][z as usize] -= 1;
				adjacency_counts[x as usize][y as usize][(z + 1) as usize] -= 1;
				while extension_stack.len() != inner_original_length {
					extension_stack.pop(); //should replace this with custom stack to avoid this unnecessary loop
				}
			}
		}

		if existing_requirement {
			required_cells.insert((x, y, z));
			break; //this required cell will no longer be available in the extension stack, so no more valid polycubes are possible in this branch
		} else {
			let mut temp_x = x;
			let mut temp_y = y;
			let mut temp_z = z;
			loop {
				let (tx, ty, tz) = (
					linear_map[0] * temp_x + linear_map[1] * temp_y + linear_map[2] * temp_z + affine_shift[0],
					linear_map[3] * temp_x + linear_map[4] * temp_y + linear_map[5] * temp_z + affine_shift[1],
					linear_map[6] * temp_x + linear_map[7] * temp_y + linear_map[8] * temp_z + affine_shift[2],
				);
				temp_x = tx;
				temp_y = ty;
				temp_z = tz;
				if x == temp_x && y == temp_y && z == temp_z {
					break;
				}
				required_cells.remove(&(temp_x, temp_y, temp_z));
			}
		}
	}
	while extension_stack.len() != original_length {
		extension_stack.push(recovery_stack.pop().unwrap());
	}
	return count;
}

fn count_extensions_subset(filter: usize, overall_start_time: Instant) -> usize {

	unsafe {
		let start_time = Instant::now();
		let mut byte_board_arr: [u8; (N as usize+ 2) * Z as usize] = [0; (N as usize + 2) * Z as usize];
		let mut byte_board = byte_board_arr.as_mut_ptr();
		let mut ref_stack_arr: [*mut u8; (N - 2) * 4] = [std::ptr::null_mut(); (N - 2) * 4];
		let ref_stack = ref_stack_arr.as_mut_ptr();
		byte_board = byte_board.offset(Z as isize); // seeded with first index of the byte board as the only allowed extension
		*ref_stack = byte_board;
		let mut i = byte_board.offset((N as isize + 1) * Z as isize);
		i = i.sub(1);
		// the first Z + 1 bytes are disallowed extensions; first Z are less than the minimum, last 1 due to edge case of initial polycube having no neighbours
		while i != byte_board {
			*i = 255;
			i = i.sub(1);
		}

		print_w_time(overall_start_time, format!("started  count_extensions_subset with filter={filter}, N={N}, FILTER_DEPTH={FILTER_DEPTH}"));
		let count: usize = count_extensions_subset_inner(N as usize, ref_stack.add(1), ref_stack.add((N - 2) * 4 as usize), ref_stack, filter);

		let time_elapsed = start_time.elapsed().as_secs_f64();
		print_w_time(overall_start_time, format!("finished count_extensions_subset with filter={filter}, N={N}, FILTER_DEPTH={FILTER_DEPTH}, took {} with count={count}", seconds_to_dur(time_elapsed)));
		return count;
	}
}

fn count_extensions_subset_inner(depth: usize, mut stack_top_1: *mut *mut u8, mut stack_top_2: *mut *mut u8, ref_stack: *mut *mut u8, filter: usize) -> usize {
	unsafe {
		let mut count: usize = 0;
		let stack_top_original = stack_top_1;
		let filter_stop = filter as isize;
		while stack_top_1 != ref_stack {
			stack_top_1 = stack_top_1.sub(1);
			let index = *stack_top_1;
			let mut stack_top_inner = stack_top_1;

			let incr = index.sub(X as usize);
			*incr = (*incr).wrapping_add(1); // thanks to https://stackoverflow.com/a/70776258/259456
			if *incr == 0 {
				*stack_top_inner = index.sub(X as usize);
				stack_top_inner = stack_top_inner.add(1);
			}
			let incr = index.sub(Y as usize);
			*incr = (*incr).wrapping_add(1);
			if *incr == 0 {
				*stack_top_inner = index.sub(Y as usize);
				stack_top_inner = stack_top_inner.add(1);
			}
			let incr = index.sub(Z as usize);
			*incr = (*incr).wrapping_add(1);
			if *incr == 0 {
				*stack_top_inner = index.sub(Z as usize);
				stack_top_inner = stack_top_inner.add(1);
			}
			let incr = index.add(X as usize);
			*incr = (*incr).wrapping_add(1);
			if *incr == 0 {
				*stack_top_inner = index.add(X as usize);
				stack_top_inner = stack_top_inner.add(1);
			}
			let incr = index.add(Y as usize);
			*incr = (*incr).wrapping_add(1);
			if *incr == 0 {
				*stack_top_inner = index.add(Y as usize);
				stack_top_inner = stack_top_inner.add(1);
			}
			let incr = index.add(Z as usize);
			*incr = (*incr).wrapping_add(1);
			if *incr == 0 {
				*stack_top_inner = index.add(Z as usize);
				stack_top_inner = stack_top_inner.add(1);
			}

			if depth == 4 {
				count += count_extensions_subset_final(stack_top_inner, ref_stack);
			} else if depth != FILTER_DEPTH || stack_top_1.offset_from(ref_stack) == filter_stop {
				count += count_extensions_subset_inner(depth - 1, stack_top_inner, stack_top_2, ref_stack, filter); // if multithreading is not wanted, remove "if (condition)" from this else statement
			}
			*index.sub(X as usize) = (*index.sub(X as usize)).wrapping_sub(1);
			*index.sub(Y as usize) = (*index.sub(Y as usize)).wrapping_sub(1);
			*index.sub(Z as usize) = (*index.sub(Z as usize)).wrapping_sub(1);
			*index.add(X as usize) = (*index.add(X as usize)).wrapping_sub(1);
			*index.add(Y as usize) = (*index.add(Y as usize)).wrapping_sub(1);
			*index.add(Z as usize) = (*index.add(Z as usize)).wrapping_sub(1);
			stack_top_2 = stack_top_2.sub(1);
			*stack_top_2 = index; // doing this push before the recursion would add one extra unnecessary element to the stack at each level of recursion
		}
		while stack_top_1 != stack_top_original {
			*stack_top_1 = *stack_top_2;
			stack_top_1 = stack_top_1.add(1);
			stack_top_2 = stack_top_2.add(1);
		}
		return count;
	}
}

fn count_extensions_subset_final(mut stack_top: *mut *mut u8, ref_stack: *mut *mut u8) -> usize {
	unsafe {
		let length = stack_top.offset_from(ref_stack);
		let mut count: isize = (length * (length - 1) * (length - 2) / 6).try_into().unwrap();
		let mut stack_top_temp = stack_top;
		while stack_top_temp != ref_stack {
			stack_top_temp = stack_top_temp.sub(1);
			let i = *stack_top_temp;
			let mut neighbours = 0;
			let mut subcount: usize = 128;

			let incr = i.sub(X as usize); if *incr > 127 { neighbours += 1; subcount += *i.sub(X2) as usize + *i.sub(SYX) as usize + *i.sub(SZX) as usize + *i.add(DYX) as usize + *i.add(DZX) as usize; *incr -= 1; count += *incr as isize; }
			let incr = i.sub(Y as usize); if *incr > 127 { neighbours += 1; subcount += *i.sub(Y2) as usize + *i.sub(SYX) as usize + *i.sub(SZY) as usize + *i.sub(DYX) as usize + *i.add(DZY) as usize; *incr -= 1; count += *incr as isize; }
			let incr = i.sub(Z as usize); if *incr > 127 { neighbours += 1; subcount += *i.sub(Z2) as usize + *i.sub(SZX) as usize + *i.sub(SZY) as usize + *i.sub(DZX) as usize + *i.sub(DZY) as usize; *incr -= 1; count += *incr as isize; }
			let incr = i.add(X as usize); if *incr > 127 { neighbours += 1; subcount += *i.add(X2) as usize + *i.add(SYX) as usize + *i.add(SZX) as usize + *i.sub(DYX) as usize + *i.sub(DZX) as usize; *incr -= 1; count += *incr as isize; }
			let incr = i.add(Y as usize); if *incr > 127 { neighbours += 1; subcount += *i.add(Y2) as usize + *i.add(SYX) as usize + *i.add(SZY) as usize + *i.add(DYX) as usize + *i.sub(DZY) as usize; *incr -= 1; count += *incr as isize; }
			let incr = i.add(Z as usize); if *incr > 127 { neighbours += 1; subcount += *i.add(Z2) as usize + *i.add(SZX) as usize + *i.add(SZY) as usize + *i.add(DZX) as usize + *i.add(DZY) as usize; *incr -= 1; count += *incr as isize; }

			// some of these values are sometimes negative
			count += ((subcount >> 8) as isize) + ((neighbours * (neighbours + (length << 1) - 511)) >> 1);
		}
		while stack_top != ref_stack {
			stack_top = stack_top.sub(1);
			let i = *stack_top;
			let incr = i.sub(X as usize); *incr |= *incr >> 4;
			let incr = i.sub(Y as usize); *incr |= *incr >> 4;
			let incr = i.sub(Z as usize); *incr |= *incr >> 4;
			let incr = i.add(X as usize); *incr |= *incr >> 4;
			let incr = i.add(Y as usize); *incr |= *incr >> 4;
			let incr = i.add(Z as usize); *incr |= *incr >> 4;
		}
		return count as usize;
	}
}
