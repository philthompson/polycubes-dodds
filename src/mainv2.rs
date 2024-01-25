/*

this program computes the values for https://oeis.org/A000162
  (except for trivially small values of N)

this is a rust port of Stanley Dodds's algorithm, whose original
  C# code is available at:
  https://oeis.org/A000162/a000162.cs.txt

rust port by Phil Thompson in January 2024, along with changes:
- use multiple threads for computation of "nontrivial symmetries" portion
- allow the work to be divided into many more thread tasks by changing FILTER_TASK relative to N
- allow specifying more threads than worker tasks
- print timestamps (elapsed time) for most messages



looks like we can fork in different thread tasks at some FILTER_DEPTH to control
  how many tasks we create (and thus how much memory we use):

tasks for FILTER_DEPTH=5:
N = 11 -> 23502
N = 12 -> 162913
N = 13 -> 1152870

tasks for FILTER_DEPTH=7:
N = 11 -> 534
N = 12 -> 3481
N = 13 -> 23502

tasks for FILTER_DEPTH=9:
N = 11 -> 15
N = 12 -> 86
N = 13 -> 534
N = 14 -> ?
N = 15 -> 23502

tasks for FILTER_DEPTH=12:
N = 21 -> 8294738

tasks for FILTER_DEPTH=14:
N = 17 -> 86
N = 19 -> 3481
N = 22 -> 1152870

tasks for FILTER_DEPTH=16:
N = 22 -> 23502
*/

use std::io::Write;
use std::time::Duration;
use std::time::Instant;
use std::thread;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::VecDeque;
// this import needed depending on rust toolchain (e.g. on Amazon Linux)
//use std::convert::TryInto;

const N: usize = 15; // number of polycube cells
const FILTER_DEPTH: usize = 12; // keep smaller than N to use multiple threads (lower=more shorter thread tasks)
const THREADS: usize = 8;
const USE_PRECOMPUTED_SYMM: bool = false; // use precomputed nontrivial symmetries, if available

// output checkpoint count+finished tasks info no more often than this many seconds
const CHECKPOINT_SECONDS: f32 = 15.0;

// the first portion of the program computes "nontrivial symmetries" which is much
//   faster than the second portion of the program.  if running a high n value, or
//   piping output to a file, then set this to false
const SHOW_NONTRIVIAL_SYMM_ETA: bool = true;

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

const BYTE_BOARD_LEN: usize = (N + 2) * (Z as usize);
const REF_STACK_LEN: usize = (N - 2) * 4;

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
	/* n=19 */   138_370_095,
	/* n=20 */             0,
	/* n=21 */ 1_055_564_170,
	/* n=22 */ 3_699_765_374];

const FREE_PORTION_CHECKPOINTS: [(usize, usize, usize, usize); 0] = [
	// N, FILTER_DEPTH, task num, cumulative free polycubes count
	//(15, 11, 418, 151499577788),
	//(16, 12, 152, 523474916610),
	//(16, 12, 330, 1449472908088),
	//(16, 12, 416, 1147058871532),
];

// for smaller values of N, we can greatly speed up computation by
//   lowering the sleep time used when waiting for threads to finish
// for larger values of N, we can check less often and waste fewer
//   cycles in the main thread
const N_SLEEP_MILLIS: u64 = if N > 19 { N as u64 * 100 } else if N > 17 { N as u64 * 25 } else if N > 12 { N as u64 * 2 } else { N as u64 + 5 };

static mut THREAD_TASK_COUNTER: usize = 0;

pub struct ThreadTask {
	pub num: usize,
	pub depth: usize,
	pub stack_top_inner_offset: isize,
	pub stack_top_2_offset: isize,
	pub ref_stack_offset: isize,
	pub byte_board: [u8; BYTE_BOARD_LEN],
	// can't actually store the pointers here, instead, need to
	//   store offsets (isize) to re-create this
	pub ref_stack: [isize; REF_STACK_LEN]
}

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
	.join(",");
}

fn main() {
	let overall_start_time = Instant::now();
	let symmetry_time_elapsed: f64;
	let mut count: usize = 0; // enumerate the sum over the order 24 group of the size of the fix of each group element, and divide by 24 (Burnside's lemma)

	/* ====-====-====-====-====-====-====
	first portion: symmetric polycubes
	====-====-====-====-====-====-==== */

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

	/* ====-====-====-====-====-====-====-====-====-====-====-====-====-====
	main/second portion: fixed polycubes (https://oeis.org/A001931)
	====-====-====-====-====-====-====-====-====-====-====-====-====-==== */
	{
		print_w_time(overall_start_time, format!("\nCounting polycubes via extensions..."));
		let n_sleep = Duration::from_millis(N_SLEEP_MILLIS);
		println!("N={N}, n_sleep_millis={N_SLEEP_MILLIS}");
		let sw = Instant::now();
		let mut subcount: usize = 0;

		// using a VecDeque so we can pop tasks off the front
		let mut tasks: VecDeque<ThreadTask>;

		let mut completed_task_counts: BTreeMap<usize, usize> = BTreeMap::new();
		let mut last_checkpoint_task_number: usize = 0;
		let mut last_checkpoint_time = Instant::now();
		let mut last_checkpoint_total: usize = 0;
		let mut competed_tasks: usize = 0;

		// find the highest available checkpoint for N+FILTER_DEPTH, if any
		for checkpoint in FREE_PORTION_CHECKPOINTS {
			if checkpoint.0 != N || checkpoint.1 != FILTER_DEPTH {
				continue;
			}
			completed_task_counts.insert(checkpoint.2, checkpoint.3);
			if checkpoint.2 > last_checkpoint_task_number {
				last_checkpoint_task_number = checkpoint.2;
				last_checkpoint_total = checkpoint.3;
				competed_tasks = last_checkpoint_task_number + 1;
				subcount = last_checkpoint_total;
			}
		}

		// initialize the Vec of thread tasks
		unsafe {
			let start_time = Instant::now();
			let mut byte_board_arr: [u8; BYTE_BOARD_LEN] = [0; BYTE_BOARD_LEN];
			let mut byte_board = byte_board_arr.as_mut_ptr();
			let mut ref_stack_arr: [*mut u8; REF_STACK_LEN] = [std::ptr::null_mut(); REF_STACK_LEN];
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
	
			print_w_time(overall_start_time, format!("started  initializing tasks, N={N}, FILTER_DEPTH={FILTER_DEPTH}"));
	
			// a vec of stacks, and offsets for stack_top_1
			// for each entry in the vec, with a separate thread, run count_extensions_subset_final
			let count_and_tasks: (usize, VecDeque<ThreadTask>) =
				count_extensions_subset_inner(N as usize, ref_stack.add(1), ref_stack.add((N - 2) * 4 as usize), ref_stack, &byte_board_arr, &ref_stack_arr);
	
			subcount += count_and_tasks.0;
			tasks = count_and_tasks.1;

			let time_elapsed = start_time.elapsed().as_secs_f64();
			print_w_time(overall_start_time, format!("finished initializing {} tasks, N={N}, FILTER_DEPTH={FILTER_DEPTH}, took {}", tasks.len(), seconds_to_dur(time_elapsed)));
		}

		let total_tasks: usize = tasks.len();
		let mut thread_handles = Vec::new();

		// spawn the desired number of threads
		while thread_handles.len() < THREADS {
			let task = tasks.pop_front();
			if task.is_none() {
				break;
			}
			let task = task.unwrap();
			if last_checkpoint_task_number > 0 && task.num <= last_checkpoint_task_number {
				continue;
			}
			let handle = thread::spawn(move || count_extensions_subset_outer(task, overall_start_time));
			thread_handles.push(handle);
		}
		let mut handle_index_to_replace: Option<usize> = None;

		// replace completed threads with new ones until there are none left
		while competed_tasks < total_tasks {
			let mut j: usize = 0;
			while j < thread_handles.len() {
				if thread_handles[j].is_finished() {
					handle_index_to_replace = Some(j);
					break;
				}
				j += 1;
			}
			if handle_index_to_replace.is_none() {
				thread::sleep(n_sleep);
				continue;
			}
			let (task_count, task_num) = thread_handles.remove(handle_index_to_replace.take().unwrap()).join().unwrap();
			completed_task_counts.insert(task_num, task_count);
			subcount += task_count;
			competed_tasks += 1;
			print_w_time(overall_start_time, format!("tasks remaining: [{total_tasks}-{competed_tasks}={}]", total_tasks-competed_tasks));

			// if enough time has passed, see if we can print a new checkpoint
			if competed_tasks > 0 && last_checkpoint_time.elapsed().as_secs_f32() > CHECKPOINT_SECONDS {
				last_checkpoint_time = Instant::now();
				let mut new_checkpoint_num = last_checkpoint_task_number;
				let mut new_checkpoint_total = last_checkpoint_total;
				// really should use an Option<usize> for last_checkpoint_task_number
				if last_checkpoint_task_number == 0 && !completed_task_counts.contains_key(&last_checkpoint_task_number) {
					new_checkpoint_num = total_tasks + 1;
				}
				// find the next longest contiguous run of completed task numbers
				let mut any_new = false;
				while new_checkpoint_num < total_tasks && completed_task_counts.contains_key(&new_checkpoint_num) {
					any_new = true;
					new_checkpoint_total += completed_task_counts.get(&new_checkpoint_num).unwrap();
					new_checkpoint_num += 1;
				}
				new_checkpoint_num -= 1;
				if any_new && new_checkpoint_num > last_checkpoint_task_number {
					print_w_time(overall_start_time, format!("checkpoint for N={N}, FILTER_DEPTH={FILTER_DEPTH}: thru task [{new_checkpoint_num}], we have cumulative free polycubes count={new_checkpoint_total}"));
					last_checkpoint_task_number = new_checkpoint_num;
					last_checkpoint_total = new_checkpoint_total;
				}
			}

			loop {
				let task = tasks.pop_front();
				if task.is_none() {
					break;
				}
				let task = task.unwrap();
				if last_checkpoint_task_number > 0 && task.num <= last_checkpoint_task_number {
					continue;
				}
				let handle = thread::spawn(move || count_extensions_subset_outer(task, overall_start_time));
				thread_handles.push(handle);
				break;
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


// rebuild state from ThreadTask, then continue execution with count_extensions_subset_inner()
fn count_extensions_subset_outer(task: ThreadTask, overall_start_time: Instant) -> (usize, usize) {
	unsafe {
		let start_time = Instant::now();
		let mut byte_board_arr: [u8; BYTE_BOARD_LEN] = task.byte_board;
		let byte_board = byte_board_arr.as_mut_ptr();
		let mut ref_stack_arr: [*mut u8; REF_STACK_LEN] = [std::ptr::null_mut(); REF_STACK_LEN];
		// restore the ref stack using offsets
		for i in 0..REF_STACK_LEN {
			ref_stack_arr[i] = byte_board.offset(task.ref_stack[i]);
		}

		let depth = task.depth;
		let ref_stack = ref_stack_arr.as_mut_ptr().offset(task.ref_stack_offset);
		let stack_top_inner = ref_stack.offset(task.stack_top_inner_offset);
		let stack_top_2 = ref_stack.offset(task.stack_top_2_offset);

		print_w_time(overall_start_time, format!("started  count_extensions_subset_outer() with task_number={}, N={N}, FILTER_DEPTH={FILTER_DEPTH}", task.num));

		let count_and_tasks: (usize, VecDeque<ThreadTask>) =
			count_extensions_subset_inner(depth, stack_top_inner, stack_top_2, ref_stack, &byte_board_arr, &ref_stack_arr);

		let time_elapsed = start_time.elapsed().as_secs_f64();
		print_w_time(overall_start_time, format!("finished count_extensions_subset_outer() with task_number={}, N={N}, FILTER_DEPTH={FILTER_DEPTH}, took {} with count={}", task.num, seconds_to_dur(time_elapsed), count_and_tasks.0));
		
		if !count_and_tasks.1.is_empty() {
			println!("count_extensions_subset_outer() got a non-empty set of tasks! ({} tasks)", count_and_tasks.1.len())
		}
		return (count_and_tasks.0, task.num);
	}
}

fn copy_ref_stack_as_offsets(byte_board_arr: &[u8; BYTE_BOARD_LEN], ref_stack_arr: &[*mut u8; REF_STACK_LEN]) -> [isize; REF_STACK_LEN] {
	let mut ref_stack_copy: [isize; REF_STACK_LEN] = [0; REF_STACK_LEN];
	let byte_board = (*byte_board_arr).as_ptr();
	unsafe {
		for i in 0..REF_STACK_LEN {
			// !!!!! !!!!! !!!!! !!!!! !!!!! !!!!! !!!!!
			// very odd... need to print (or at least format a string) to avoid segfauly
			// !!!!! !!!!! !!!!! !!!!! !!!!! !!!!! !!!!!
			let o = ref_stack_arr[i].offset_from(byte_board);
			let _ = format!("offset for ref_stack[{i}]: {o}");
			ref_stack_copy[i] = o;
		}
	}
	return ref_stack_copy;
}

fn count_extensions_subset_inner(
		depth: usize,
		mut stack_top_1: *mut *mut u8,
		mut stack_top_2: *mut *mut u8,
		ref_stack: *mut *mut u8,
		byte_board_arr: &[u8; BYTE_BOARD_LEN],
		ref_stack_arr: &[*mut u8; REF_STACK_LEN]) -> (usize, VecDeque<ThreadTask>) {
	unsafe {
		let mut count: usize = 0;
		let mut tasks: VecDeque<ThreadTask> = VecDeque::new();
		let stack_top_original = stack_top_1;
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

			} else if depth != FILTER_DEPTH {
				let mut deeper = count_extensions_subset_inner(depth - 1, stack_top_inner, stack_top_2, ref_stack, byte_board_arr, ref_stack_arr);
				count += deeper.0;
				tasks.append(&mut deeper.1);
			
			// at depth == FILTER_DEPTH, create a copy of the state to be later resumed by a thread
			} else {
				let task = ThreadTask {
					num: THREAD_TASK_COUNTER, 
					depth: depth - 1,
					stack_top_inner_offset: stack_top_inner.offset_from(ref_stack),
					stack_top_2_offset: stack_top_2.offset_from(ref_stack),
					ref_stack_offset: ref_stack.offset_from(ref_stack_arr.as_ptr()),
					byte_board: byte_board_arr.clone(),
					ref_stack: copy_ref_stack_as_offsets(byte_board_arr, ref_stack_arr)
				};
				tasks.push_back(task);
				THREAD_TASK_COUNTER += 1;
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
		return (count, tasks);
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
