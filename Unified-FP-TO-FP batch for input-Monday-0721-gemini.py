import time
from datetime import datetime
from collections import Counter
import sys

# --- GLOBAL VARIABLES FOR PROGRAM 1's SOLVER (Crucial for its behavior) ---
# These must be reset for each call to _find_sequence if used outside
# the _find_sequence's encapsulating scope.
_successful_permutation_found_global = False
_found_path_nodes_global = None
_total_backtracks_for_sequence_search_global = 0

# --- Core Utilities (Unchanged from previous versions unless specified) ---

def cyclic_len(u, v, p):
    """Calculates the cyclic length (shortest distance) between two vertices u and v
    in a cycle of p vertices."""
    return min(abs(u - v), p - abs(u - v))

def edge_freq(path, p):
    """Calculates the frequency of each distinct edge length (hop) in a given path.
    'p' is the total number of vertices in the cycle/graph."""
    if len(path) < 2:
        return Counter()
    
    if len(set(path)) != p:
        return Counter()

    if not all(0 <= node < p for node in path):
        return Counter()

    return Counter(cyclic_len(path[i], path[i+1], p) for i in range(len(path)-1))

def parse_fp_tuple(fp_tuple):
    """Converts a frequency partition tuple (e.g., (2, 1, 0)) into a Counter object.
    Keys are hop lengths (1-indexed), values are their frequencies."""
    return Counter({i+1: fp_tuple[i] for i in range(len(fp_tuple))})

def satisfies_divisor_condition(fp, v):
    """Checks the BHR (Bandwidth-Hop Restriction) necessary condition for a Hamiltonian path.
    For every divisor 'd' of 'v', the sum of counts of hops that are multiples of 'd'
    must not be greater than v - d.
    """
    divisors = [d for d in range(1, v + 1) if v % d == 0]
    for d in divisors:
        count = sum(cnt for hop, cnt in fp.items() if hop % d == 0)
        if count > v - d:
            return False, d, count
    return True, None, None

def reuse_insert(H, fp2, p_curr):
    """Attempts to construct a new Hamiltonian path H2 by inserting a new vertex
    (p_curr - 1, which is len(H)) into all possible positions of an existing path H,
    and checks if its edge frequencies match fp2.
    """
    t0 = time.time()
    new_v_val = p_curr - 1
    
    if new_v_val in H:
        return None, time.time() - t0

    for pos in range(len(H) + 1):
        H2 = H[:pos] + [new_v_val] + H[pos:]
        
        if len(set(H2)) != p_curr:
            continue 

        if edge_freq(H2, p_curr) == fp2:
            return H2, time.time() - t0
    return None, time.time() - t0

def greedy_insert(H, fp2, p_curr, top_k=3):
    """Attempts to construct a new Hamiltonian path H2 by inserting a new vertex
    into an existing path H. It ranks insertions by how closely their edge
    frequencies match fp2 and checks the top_k candidates."""
    t0 = time.time()
    new_v_val = p_curr - 1
    candidates = []

    if new_v_val in H:
        return None, time.time() - t0

    for pos in range(len(H) + 1):
        H2_temp = H[:pos] + [new_v_val] + H[pos:]
        
        if len(set(H2_temp)) != p_curr:
             continue 

        current_freq = edge_freq(H2_temp, p_curr)
        score = sum(abs(current_freq.get(k, 0) - fp2.get(k, 0)) for k in set(current_freq) | set(fp2))
        candidates.append((score, H2_temp))
        
    candidates.sort(key=lambda x: x[0])
    
    for score, H2 in candidates[:top_k]:
        if score == 0:
            return H2, time.time() - t0
    return None, time.time() - t0


# --- PROGRAM 1's EXACT BACKTRACKING SOLVER LOGIC ---

def _recursive_sequence_builder(p, current_path_nodes, current_prefix_hops, remaining_multiset_hops):
    global _successful_permutation_found_global
    global _found_path_nodes_global
    global _total_backtracks_for_sequence_search_global

    if _successful_permutation_found_global: # Optimization: If a path is already found, stop
        return True

    # Base case: if a full path of `p` nodes has been built
    if len(current_path_nodes) == p:
        # Check if all hops from the multiset have been used
        if not remaining_multiset_hops: # If multiset is empty, all hops are used
            _successful_permutation_found_global = True
            _found_path_nodes_global = current_path_nodes[:] # Store a copy of the successful path
            return True
        else:
            # Path length is correct, but not all required hops were consumed
            # (e.g., ran out of available nodes before using all hops, or vice versa)
            return False 

    # If the current node count exceeds P, or if there are no hops left but path is not full
    # These are implicit checks by `len(current_path_nodes) == p` and `remaining_multiset_hops`.
    # No explicit checks needed here.

    current_node = current_path_nodes[-1]

    # Iterate through each *available hop value* in the remaining multiset
    # Creating a list for iteration allows safe modification during pop/insert later
    distinct_available_hops = sorted(list(set(remaining_multiset_hops))) # Get unique values, sorted for consistency

    for hop_val in distinct_available_hops:
        if hop_val in remaining_multiset_hops: # Check if this hop value is still in the multiset
            remaining_multiset_hops.remove(hop_val) # Temporarily remove for this branch

            # Try both directions: current_node + hop and current_node - hop
            for next_node in [(current_node + hop_val) % p, (current_node - hop_val + p) % p]:
                if next_node not in current_path_nodes: # Check if node is already visited
                    current_path_nodes.append(next_node)
                    current_prefix_hops.append(hop_val) # Record the hop used for this step

                    if _recursive_sequence_builder(p, current_path_nodes, current_prefix_hops, remaining_multiset_hops):
                        return True # Solution found
                    
                    # Backtrack: undo the choice
                    current_path_nodes.pop()
                    current_prefix_hops.pop()
                else:
                    _total_backtracks_for_sequence_search_global += 1 # Increment backtrack count if node already used

            # Backtrack: add the hop back to the multiset
            remaining_multiset_hops.append(hop_val)

    return False # No solution found from this branch


def _find_sequence(p, initial_multiset_list, hint=None):
    """
    Main entry point for Program 1's solver.
    Initializes global variables and calls the recursive builder.
    """
    global _successful_permutation_found_global
    global _found_path_nodes_global
    global _total_backtracks_for_sequence_search_global

    # Reset global variables for each new search
    _successful_permutation_found_global = False
    _found_path_nodes_global = None
    _total_backtracks_for_sequence_search_global = 0

    t0 = time.time()
    
    # Start the recursive search from node 0
    # The initial multiset must be a list of individual hop values
    # e.g., for Counter({1:2, 3:1}), this should be [1,1,3]
    # This is critical for the _recursive_sequence_builder's pop/append logic.
    
    # The `initial_multiset_list` passed to this function will be derived from the FP counter.
    
    if _recursive_sequence_builder(p, [0], [], initial_multiset_list):
        return _found_path_nodes_global, time.time() - t0, _total_backtracks_for_sequence_search_global
    else:
        return None, time.time() - t0, _total_backtracks_for_sequence_search_global

# --- Batch Job (Modified to use _find_sequence) ---

def batch_extend(HP1, FP1_tuple, num_iterations=10, mode=2):
    """
    Extends an initial Hamiltonian Path (HP1) based on an evolving Frequency Partition (FP).
    Applies BHR check and various construction methods.
    """
    base_fp_list = list(FP1_tuple)
    base_path = HP1.copy()
    results = []
    timestamp = datetime.now().strftime("%H%M%S")
    log_file = f"log_unified_{timestamp}.txt"

    print(f"\nStarting batch job. Log file: {log_file}")

    for i in range(num_iterations):
        fp_prev = parse_fp_tuple(tuple(base_fp_list))
        p_prev = len(base_path)

        # --- Evolve FP for current iteration ---
        if mode == 1:
            if len(base_fp_list) == 0:
                 print(f"⚠️ Iteration {i+1} skipped: No hops in FP to increment for P={p_prev}.")
                 continue
            
            idx = i % len(base_fp_list)
            base_fp_list[idx] += 1
            
        elif mode == 2:
            max_possible_hop_value = (p_prev + 1) // 2
            if (len(base_fp_list) + 1) > max_possible_hop_value:
                print(f"⚠️ Iteration {i+1} skipped: Adding a new hop would exceed max possible hop length ({max_possible_hop_value}) for P={p_prev + 1}.")
                t_stamp = datetime.now().strftime("%H:%M:%S")
                with open(log_file, "a", encoding="utf-8") as f:
                    f.write(f"--- Iteration {i+1} ---\n")
                    f.write(f"Timestamp: {t_stamp}\n")
                    f.write(f"Method used: SKIPPED\n")
                    f.write(f"p (vertices): {p_prev + 1} (would be)\n")
                    f.write(f"Previous FP: {fp_prev}\n")
                    f.write(f"Evolved FP:  {parse_fp_tuple(tuple(base_fp_list + [1]))} (attempted)\n")
                    f.write(f"⚠️ Skipped: Cannot add new hop as it would exceed ⌊p/2⌋ for target P.\nResult: SKIPPED\n\n")
                continue
            base_fp_list.append(1)
            
        else:
            print("❌ Invalid mode selected. Exiting batch job.")
            return

        fp_curr = parse_fp_tuple(tuple(base_fp_list))
        p_curr = sum(fp_curr.values()) + 1
        t_stamp = datetime.now().strftime("%H:%M:%S")

        print(f"\n--- Iteration {i+1} (p={p_curr}) ---")
        print(f"Evolving FP from {fp_prev} to {fp_curr}")

        # ✅ Check BHR divisor condition
        is_valid, bad_div, bad_count = satisfies_divisor_condition(fp_curr, p_curr)
        if not is_valid:
            print(f"❌ Iteration {i+1}: FP = {fp_curr}")
            print(f"❌ Failed BHR necessary condition at divisor {bad_div}: count = {bad_count} > {p_curr - bad_div}")
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"--- Iteration {i+1} ---\n")
                f.write(f"Timestamp: {t_stamp}\n")
                f.write(f"Method used: BHR Fail\n")
                f.write(f"p (vertices): {p_curr}\n")
                f.write(f"Previous FP: {fp_prev}\n")
                f.write(f"Evolved FP:  {fp_curr}\n")
                f.write(f"❌ BHR failed at divisor {bad_div}: count = {bad_count} > {p_curr - bad_div}\nResult: SKIPPED\n\n")
            continue
        else:
            print(f"✅ BHR passed for iteration {i+1}")

        # Attempt path construction methods
        H2 = None
        t_elapsed = 0
        backtracks = 0
        stage = "Unknown"

        # 1. Try Reuse-Insert
        print(f"{datetime.now().strftime('[%H:%M:%S]')} Attempting reuse-insert...")
        H2, t_elapsed_reuse = reuse_insert(base_path, fp_curr, p_curr)
        t_elapsed += t_elapsed_reuse
        stage = "reuse-insert"

        if not H2:
            print(f"{datetime.now().strftime('[%H:%M:%S]')} ❌ Reuse-insert failed. Attempting greedy-insert...")
            # 2. Try Greedy-Insert if Reuse failed
            H2, t_elapsed_greedy = greedy_insert(base_path, fp_curr, p_curr)
            t_elapsed += t_elapsed_greedy
            stage = "greedy-insert"

            if not H2:
                print(f"{datetime.now().strftime('[%H:%M:%S]')} ❌ Greedy-insert failed. Falling back to Program 1's _find_sequence...")
                # 3. Call Program 1's _find_sequence
                # Convert the Counter FP to a flat list (multiset) of hop values
                initial_multiset_for_solver = []
                for hop_val, count in fp_curr.items():
                    initial_multiset_for_solver.extend([hop_val] * count)

                H2, t_elapsed_back, backtracks = _find_sequence(p_curr, initial_multiset_for_solver)
                t_elapsed += t_elapsed_back
                stage = "backtrack (Program 1 solver)"

        elapsed_formatted = round(t_elapsed, 6)

        # Log results
        with open(log_file, "a", encoding="utf-8") as f:
            f.write(f"--- Iteration {i+1} ---\n")
            f.write(f"Timestamp: {t_stamp}\n")
            f.write(f"Method used: {stage}\n")
            f.write(f"p (vertices): {p_curr}\n")
            f.write(f"Previous FP: {fp_prev}\n")
            f.write(f"Evolved FP:  {fp_curr}\n")
            f.write("✅ BHR passed\n")

            if H2:
                f.write(f"✅ HP: {H2}\n")
                f.write(f"HP freq: {edge_freq(H2, p_curr)}\n")
                f.write(f"🔙 Backtracks: {backtracks}\n")
                f.write(f"⏱ Time: {elapsed_formatted:.2f} sec\nResult: SUCCESS\n\n")
                results.append((i + 1, p_curr, H2, fp_curr, backtracks))
                base_path = H2
                print(f"✅ HP found (Method: {stage}, Time: {elapsed_formatted:.2f}s, Backtracks: {backtracks})")
            else:
                f.write("❌ Failed to construct HP\n")
                f.write(f"🔙 Backtracks: {backtracks}\n")
                f.write(f"⏱ Time: {elapsed_formatted:.2f} sec\nResult: FAILURE\n\n")
                print(f"❌ Failed to construct HP (Method: {stage}, Time: {elapsed_formatted:.2f}s, Backtracks: {backtracks})")

    # Final Summary Log
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(f"--- Final Summary ---\n")
        f.write(f"Total Successful Constructions: {len(results)}\n")
        for idx, p, hp, fp, bktrk in results:
            f.write(f"Iteration {idx}: p = {p}, Backtracks: {bktrk}\n")
        f.write("\nEnd of log.\n")

    print(f"\n✅ Batch complete. Log saved to: {log_file}")
    return results

# --- Menu Entry Point ---

def get_mode():
    """Prompts the user to choose the mode for evolving the frequency partition."""
    while True:
        print("\nChoose Mode:")
        print("1 - Evolve with same number of parts (increment existing hops)")
        print("2 - Add one new edge length (if allowed)")
        choice = input("Your choice: ").strip()
        if choice in {"1", "2"}:
            return int(choice)
        print("❌ Invalid choice. Please enter 1 or 2.")

def get_int_list(prompt):
    """Prompts the user for a comma-separated list of integers and returns it as a list."""
    while True:
        try:
            entries = input(f"{prompt} (comma-separated): ").strip()
            nums = list(map(int, entries.split(',')))
            return nums
        except ValueError:
            print("❌ Invalid input. Please enter only comma-separated integers.")

def get_int_tuple(prompt):
    """Prompts the user for a comma-separated list of integers and returns it as a tuple."""
    return tuple(get_int_list(prompt))

def get_iteration_count():
    """Prompts the user for the number of iterations, ensuring it's within a valid range."""
    while True:
        try:
            n = int(input("Number of iterations (1–50): ").strip())
            if 1 <= n <= 50:
                return n
            print("⚠️ Must be between 1 and 50.")
        except ValueError:
            print("❌ Invalid input. Enter a number.")

if __name__ == "__main__":
    print("🧠 Hamiltonian Path Generator (Batch Mode with BHR Check)")
    mode_choice = get_mode()
    HP1 = get_int_list("Enter Hamiltonian Path")

    actual_total_edges_in_hp1 = len(HP1) - 1
    if actual_total_edges_in_hp1 < 1:
        print("❌ Error: Hamiltonian Path must have at least 2 vertices (1 edge).")
        sys.exit(1)

    while True:
        FP1_tuple = get_int_tuple("Enter Frequency Partition Tuple")
        parsed_fp1 = parse_fp_tuple(FP1_tuple)
        expected_total_edges_from_fp = sum(parsed_fp1.values())

        if expected_total_edges_from_fp != actual_total_edges_in_hp1:
            print("\n❌ Error: Mismatch between provided HP and FP.")
            print(f"    Your HP has {len(HP1)} vertices, implying {actual_total_edges_in_hp1} edges.")
            print(f"    Your FP implies a path with {expected_total_edges_from_fp} edges (sum of frequencies).")
            print("    Please re-enter your Frequency Partition Tuple to match the HP.")
        else:
            print(f"✅ FP matches HP (Total edges: {actual_total_edges_in_hp1}).")
            break

    num_iterations = get_iteration_count()

    batch_extend(HP1, FP1_tuple, num_iterations=num_iterations, mode=mode_choice)
