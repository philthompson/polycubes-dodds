/*
this is code mostly generated by Microsoft Copilot, using GPT4 on Dec. 29th-30th 2023,
  when asking to port pasted C# code to Rust

i fixed and replaced a few things to make it compile, and it's initially finding
  fewer than expected polycubes for each rotational symmetry, and it's panic'ing
  in count_extensions_subset_inner() on an "attempt to add with overflow"

for initial testing, i found that the C# code produces the expected final count
  for N = 6, FILTER_DEPTH = 5, THREADS = 2, so i am initially using that to
  compare both programs' output

as of Jan. 2nd 2024, i've fixed up the ported rust code and it's now performing
   identically to the original C# code.

------

going to parallelize the initial symmetry counting portion

before parallelizing running time:
7859514470 free polycubes with 15 cells
overall time: 00h:00m:34.250s

after parallelizing running time:
7859514470 free polycubes with 15 cells
overall time: 00h:00m:11.723s

the speed of this thing is insane!  on my M1 Mac mini,
with 4 performance cores and 4 efficiency cores:
457409613979 free polycubes with 17 cells
overall time: 00h:05m:19.144s

5 minutes to calculate n=17!  the other algorithm i
  was optimizing would take 30+ days for n=17!

*/

use std::io::Write;
use std::time::Duration;
use std::time::Instant;
use std::thread;
use std::collections::HashSet;

const N: usize = 16; // number of polycube cells. Need n >= 4 if single threading, or n >= filterDepth >= 5 if multithreading (I think)
const FILTER_DEPTH: usize = 6;
const THREADS: i32 = 3;

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

pub fn seconds_to_dur(s: f64) -> String {
	let days = (s / 86400.0).floor();
	let hours = ((s - (days * 86400.0)) / 3600.0).floor();
	let minutes = ((s - (days * 86400.0) - (hours * 3600.0)) / 60.0).floor();
	let seconds = s - (days * 86400.0) - (hours * 3600.0) - (minutes * 60.0);
	let fsec = format!("{}{:.3}", if seconds < 10.0 { "0" } else { "" }, seconds);
	if days > 0.0 {
		return format!("{} days {:0>2}h:{:0>2}m:{}s", days, hours, minutes, fsec);
	}
	return format!("{:0>2}h:{:0>2}m:{}s", hours, minutes, fsec);
}

fn main() {
	let overall_start_time = Instant::now();
	let symmetry_time_elapsed: f64;
	let mut count: usize = 0; // enumerate the sum over the order 24 group of the size of the fix of each group element, and divide by 24 (Burnside's lemma)
	{
		//let total_sym_count = (4 * ((2 * N) - 1)) as f64;
		//let mut sym_count_done: f64 = 0.0;
		//for sym in 0..4 {
		//	let mut subcount = 0;
		//	for i in (1 - (N as i32))..(N as i32) {
		//		for j in (1 - (N as i32) + i.abs())..=((N as i32) - 1 - i.abs()) {
		//			let handle = thread::spawn(move || count_symmetric_polycubes_task(i, j, sym));
		//			total_worker_jobs += 1.0;
		//			tasks.push(handle);
		//			subcount += count_symmetric_polycubes(
		//				&MATRIX_REPS[sym],
		//				[
		//					i * AFFINE1[sym][0] + j * AFFINE2[sym][0] + BIASES[sym][0],
		//					i * AFFINE1[sym][1] + j * AFFINE2[sym][1] + BIASES[sym][1],
		//					i * AFFINE1[sym][2] + j * AFFINE2[sym][2] + BIASES[sym][2]
		//				]);
		//		}
		//	}
		//	println!("\n{} polycubes fixed under each {}", subcount, DESCRIPTIONS[sym]);
		//	println!("{} in total (symmetry occurs {} times)", AUT_CLASS_SIZES[sym] * subcount, AUT_CLASS_SIZES[sym]);
		//	count += AUT_CLASS_SIZES[sym] * subcount;
		//}

		// run a batch of worker tasks for each of the 4 symmetries
		for sym in 0..4 {
			println!("\nCounting {} symmetries (set {} of 4):", DESCRIPTIONS[sym], sym+1);
			println!("note: percentage complete and ETA is just for this symmetry, not the entire overall count)");
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
					} else {
						thread::sleep(Duration::from_millis(50));
					}
					j += 1;
				}
				if last_stats.elapsed().as_secs_f32() > 1.0 {
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
			println!("{} in total (symmetry occurs {} times)", AUT_CLASS_SIZES[sym] * subcount, AUT_CLASS_SIZES[sym]);
			count += AUT_CLASS_SIZES[sym] * subcount;
		}

		symmetry_time_elapsed = overall_start_time.elapsed().as_secs_f64();
		println!("\ntotal count for nontrivial symmetries is {} for polycubes with {} cells\nTook {}", count, N, seconds_to_dur(symmetry_time_elapsed));
	}

	{
		println!("\nCounting polycubes via extensions...");
		let sw = Instant::now();
		let mut subcount: usize = 0;
		let mut tasks = Vec::new(); // please don't judge my multithreading - this was basically just my first attempt, to see if it helped
		let mut jobs_remaining: isize = (4 * (N as i32 - FILTER_DEPTH as i32) - 2).try_into().unwrap(); // this is the maximum left stack length, which is the value being filtered to separate work
		let total_worker_jobs = jobs_remaining as f64;
		for _ in 0..THREADS {
			let filter = jobs_remaining; // copy of i, since lambda expression captures the variable
			jobs_remaining -= 1;
			let handle = thread::spawn(move || count_extensions_subset(filter as usize));
			tasks.push(handle);
		}
		let mut handle_index_to_replace: Option<usize> = None;
		let mut last_stats = Instant::now();
		let mut compl_worker_jobs: f64 = 0.0;
		// why should this continue once more when jobs_remaining == 0?
		while jobs_remaining >= 0 {
			let mut j: usize = 0;
			while j < tasks.len() {
				if tasks[j].is_finished() {
					handle_index_to_replace = Some(j);
					compl_worker_jobs += 1.0;
					break;
				} else {
					thread::sleep(Duration::from_millis(50));
				}
				j += 1;
			}
			if last_stats.elapsed().as_secs_f32() > 1.0 {
				last_stats = Instant::now();
				let time_elapsed = sw.elapsed();
				let seconds_per_thread = time_elapsed.as_secs_f64() / compl_worker_jobs;
				let threads_remaining = jobs_remaining as f64;
				let seconds_remaining = threads_remaining * seconds_per_thread;
				let pct_complete = (compl_worker_jobs * 100.0) / total_worker_jobs;
				let total_seconds = seconds_remaining + time_elapsed.as_secs_f64() + symmetry_time_elapsed;
				print!("    {:.4}% complete, ETA:[{}], total:[{}], counting for n={}, outstanding threads:[{}-{}={}]        \r",
					pct_complete,
					seconds_to_dur(seconds_remaining),
					seconds_to_dur(total_seconds),
					N,
					total_worker_jobs as usize,
					compl_worker_jobs as usize,
					(total_worker_jobs - compl_worker_jobs) as usize);
				std::io::stdout().flush().unwrap();
			}
			if handle_index_to_replace.is_none() {
				continue;
			}
			subcount += tasks.remove(handle_index_to_replace.take().unwrap()).join().unwrap();
			let filter = jobs_remaining; // copy of i, since lambda expression captures the variable
			jobs_remaining -= 1;
			let handle = thread::spawn(move || count_extensions_subset(filter as usize));
			tasks.push(handle);
		}
		subcount += tasks.into_iter().map(|t| t.join().unwrap()).sum::<usize>(); // collect the results from all tasks
		count += subcount;
		println!("\n\n{:} polycubes with {} cells (number of polycubes fixed by trivial symmetry)\nTook {}", subcount, N, seconds_to_dur(sw.elapsed().as_secs_f64()));
	}
	println!();
	count /= 24;
	println!("{:} free polycubes with {} cells", count, N);
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
	let mut required_cells = HashSet::new();
	let mut extension_stack: Vec<(i32, i32, i32)> = Vec::new();
	extension_stack.push((N as i32, N as i32, 1));
	let mut recovery_stack: Vec<(i32, i32, i32)> = Vec::new();
	return count_extensions(N as i32, &mut adjacency_counts, &mut required_cells, &mut extension_stack, &mut recovery_stack, linear_map, affine_shift);
}

fn count_extensions(
		mut cells_to_add: i32,
		adjacency_counts: &mut Vec<Vec<Vec<u8>>>,
		required_cells: &mut HashSet<(i32, i32, i32)>,
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
			//let mut iter: usize = 0; // iter is just for debugging and can be removed later
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

fn count_extensions_subset(filter: usize) -> usize {

	unsafe {
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

		//println!("started task {}", filter);
		let count: usize = count_extensions_subset_inner(N as usize, ref_stack.add(1), ref_stack.add((N - 2) * 4 as usize), ref_stack, filter);
		//println!("finished task {} with subcount {}", filter, count);
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
			*incr = (*incr).wrapping_add(1); // thanks to https://stackoverflow.com/a/70776258/259456
			if *incr == 0 {
				*stack_top_inner = index.sub(Y as usize);
				stack_top_inner = stack_top_inner.add(1);
			}
			let incr = index.sub(Z as usize);
			*incr = (*incr).wrapping_add(1); // thanks to https://stackoverflow.com/a/70776258/259456
			if *incr == 0 {
				*stack_top_inner = index.sub(Z as usize);
				stack_top_inner = stack_top_inner.add(1);
			}
			let incr = index.add(X as usize);
			*incr = (*incr).wrapping_add(1); // thanks to https://stackoverflow.com/a/70776258/259456
			if *incr == 0 {
				*stack_top_inner = index.add(X as usize);
				stack_top_inner = stack_top_inner.add(1);
			}
			let incr = index.add(Y as usize);
			*incr = (*incr).wrapping_add(1); // thanks to https://stackoverflow.com/a/70776258/259456
			if *incr == 0 {
				*stack_top_inner = index.add(Y as usize);
				stack_top_inner = stack_top_inner.add(1);
			}
			let incr = index.add(Z as usize);
			*incr = (*incr).wrapping_add(1); // thanks to https://stackoverflow.com/a/70776258/259456
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

			// corrected:
			let incr = i.sub(X as usize); if *incr > 127 { neighbours += 1; subcount += *i.sub(X2) as usize + *i.sub(SYX) as usize + *i.sub(SZX) as usize + *i.add(DYX) as usize + *i.add(DZX) as usize; *incr -= 1; count += *incr as isize; }
			let incr = i.sub(Y as usize); if *incr > 127 { neighbours += 1; subcount += *i.sub(Y2) as usize + *i.sub(SYX) as usize + *i.sub(SZY) as usize + *i.sub(DYX) as usize + *i.add(DZY) as usize; *incr -= 1; count += *incr as isize; }
			let incr = i.sub(Z as usize); if *incr > 127 { neighbours += 1; subcount += *i.sub(Z2) as usize + *i.sub(SZX) as usize + *i.sub(SZY) as usize + *i.sub(DZX) as usize + *i.sub(DZY) as usize; *incr -= 1; count += *incr as isize; }
			let incr = i.add(X as usize); if *incr > 127 { neighbours += 1; subcount += *i.add(X2) as usize + *i.add(SYX) as usize + *i.add(SZX) as usize + *i.sub(DYX) as usize + *i.sub(DZX) as usize; *incr -= 1; count += *incr as isize; }
			let incr = i.add(Y as usize); if *incr > 127 { neighbours += 1; subcount += *i.add(Y2) as usize + *i.add(SYX) as usize + *i.add(SZY) as usize + *i.add(DYX) as usize + *i.sub(DZY) as usize; *incr -= 1; count += *incr as isize; }
			let incr = i.add(Z as usize); if *incr > 127 { neighbours += 1; subcount += *i.add(Z2) as usize + *i.add(SZX) as usize + *i.add(SZY) as usize + *i.add(DZX) as usize + *i.add(DZY) as usize; *incr -= 1; count += *incr as isize; }

			//let incr = i.sub(X as usize);
			//if *incr > 127 {
			//	neighbours += 1;
			//	println!("neighbors incremented to {neighbours}");
			//	print!("-X: subcount goes from {subcount} to ");
			//	//subcount += *incr.sub(X2) as usize + *incr.sub(SZX) as usize + *incr.add(DYX) as usize + *incr.add(DZX) as usize;
			//	subcount += *i.sub(X2) as usize + *i.sub(SYX) as usize + *i.sub(SZX) as usize + *i.add(DYX) as usize + *i.add(DZX) as usize;
			//	//*incr.sub(X as usize) -= 1;
			//	println!("{subcount}");
			//	//println!("added {}, {}, {}, {}, {}", *incr.sub(X2), *incr.sub(SYX), *incr.sub(SZX), *incr.add(DYX), *incr.add(DZX));
			//	println!("added {}, {}, {}, {}, {}", *i.sub(X2), *i.sub(SYX), *i.sub(SZX), *i.add(DYX), *i.add(DZX));
			//	*incr -= 1;
			//	count += *incr as isize;
			//	println!("count is now {count}");
			//}
			//let incr = i.sub(Y as usize);
			//if *incr > 127 {
			//	neighbours += 1;
			//	println!("neighbors incremented to {neighbours}");
			//	print!("-Y: subcount goes from {subcount} to ");
			//	//subcount += *incr.sub(Y2) as usize + *incr.sub(SZY) as usize + *incr.sub(DYX) as usize + *incr.add(DZY) as usize;
			//	subcount += *i.sub(Y2) as usize + *i.sub(SYX) as usize + *i.sub(SZY) as usize + *i.sub(DYX) as usize + *i.add(DZY) as usize;
			//	//*incr.sub(Y as usize) -= 1;
			//	println!("{subcount}");
			//	*incr -= 1;
			//	count += *incr as isize;
			//	println!("count is now {count}");
			//}
			//let incr = i.sub(Z as usize);
			//if *incr > 127 {
			//	neighbours += 1;
			//	println!("neighbors incremented to {neighbours}");
			//	print!("-Z: subcount goes from {subcount} to ");
			//	//subcount += *incr.sub(Z2) as usize + *incr.sub(SZY) as usize + *incr.sub(DZX) as usize + *incr.sub(DZX) as usize;
			//	subcount += *i.sub(Z2) as usize + *i.sub(SZX) as usize + *i.sub(SZY) as usize + *i.sub(DZX) as usize + *i.sub(DZY) as usize;
			//	println!("{subcount}");
			//	*incr -= 1;
			//	count += *incr as isize;
			//	println!("count is now {count}");
			//}
			//let incr = i.add(X as usize);
			//if *incr > 127 {
			//	neighbours += 1;
			//	println!("neighbors incremented to {neighbours}");
			//	print!("+X: subcount goes from {subcount} to ");
			//	//subcount += *incr.add(X2) as usize + *incr.add(SZX) as usize + *incr.sub(DYX) as usize + *incr.sub(DYX) as usize;
			//	subcount += *i.add(X2) as usize + *i.add(SYX) as usize + *i.add(SZX) as usize + *i.sub(DYX) as usize + *i.sub(DZX) as usize;
			//	println!("{subcount}");
			//	println!("added {}, {}, {}, {}, {}", *i.add(X2), *i.add(SYX), *i.add(SZX), *i.sub(DYX), *i.sub(DZX));
			//	*incr -= 1;
			//	count += *incr as isize;
			//	println!("count is now {count}");
			//}
			//let incr = i.add(Y as usize);
			//if *incr > 127 {
			//	neighbours += 1;
			//	println!("neighbors incremented to {neighbours}");
			//	print!("+Y: subcount goes from {subcount} to ");
			//	//subcount += *incr.add(Y2) as usize + *incr.add(SZY) as usize + *incr.add(DYX) as usize + *incr.sub(DYX) as usize;
			//	subcount += *i.add(Y2) as usize + *i.add(SYX) as usize + *i.add(SZY) as usize + *i.add(DYX) as usize + *i.sub(DZY) as usize;
			//	println!("{subcount}");
			//	//println!("added {} + {} + {} + {}", *incr.add(Y2), *incr.add(SZY), *incr.add(DYX), *incr.sub(DYX));
			//	*incr -= 1;
			//	count += *incr as isize;
			//	println!("count is now {count}");
			//}
			//let incr = i.add(Z as usize);
			//if *incr > 127 {
			//	neighbours += 1;
			//	println!("neighbors incremented to {neighbours}");
			//	print!("+Z: subcount goes from {subcount} to ");
			//	//subcount += *incr.add(Z2) as usize + *incr.add(SZY) as usize + *incr.add(DZX) as usize + *incr.add(DZY) as usize;
			//	subcount += *i.add(Z2) as usize + *i.add(SZX) as usize + *i.add(SZY) as usize + *i.add(DZX) as usize + *i.add(DZY) as usize;
			//	println!("{subcount}");
			//	*incr -= 1;
			//	count += *incr as isize;
			//	println!("count is now {count}");
			//}

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