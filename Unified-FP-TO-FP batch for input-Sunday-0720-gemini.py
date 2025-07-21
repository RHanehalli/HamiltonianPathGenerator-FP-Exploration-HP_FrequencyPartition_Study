import time
from datetime import datetime
from collections import Counter

# --- Core Utilities ---

def cyclic_len(u, v, p):
    """Calculates the cyclic length (shortest distance) between two vertices u and v
    in a cycle of p vertices."""
    return min(abs(u - v), p - abs(u - v))

def edge_freq(path):
    """Calculates the frequency of each distinct edge length (hop) in a given path.
    The 'p' for cyclic_len is derived from the length of the path itself,
    as a Hamiltonian path on p vertices has p-1 edges.
    However, for cyclic_len, p should be the total number of vertices in the cycle,
    which for a path of length n (n vertices) is also n.
    """
    if len(path) < 2:
        return Counter() # No edges if path has less than 2 vertices
    p = len(path) # The number of vertices in the current path, used as the cycle length
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

def reuse_insert(H, fp2):
    """Attempts to construct a new Hamiltonian path H2 by inserting a new vertex
    (len(H)) into all possible positions of an existing path H, and checks if
    its edge frequencies match fp2."""
    t0 = time.time()
    new_v = len(H) # The value of the new vertex to insert
    for pos in range(len(H) + 1):
        H2 = H[:pos] + [new_v] + H[pos:]
        if edge_freq(H2) == fp2:
            return H2, time.time() - t0
    return None, time.time() - t0

def greedy_insert(H, fp2, top_k=3):
    """Attempts to construct a new Hamiltonian path H2 by inserting a new vertex
    into an existing path H. It ranks insertions by how closely their edge
    frequencies match fp2 and checks the top_k candidates."""
    t0 = time.time()
    new_v = len(H)
    candidates = []
    for pos in range(len(H) + 1):
        H2_temp = H[:pos] + [new_v] + H[pos:]
        # Calculate score: sum of absolute differences between actual and target frequencies
        current_freq = edge_freq(H2_temp)
        score = sum(abs(current_freq.get(k, 0) - fp2.get(k, 0)) for k in set(current_freq) | set(fp2))
        candidates.append((score, H2_temp))
    
    candidates.sort(key=lambda x: x[0]) # Sort by score (lower is better)
    
    for score, H2 in candidates[:top_k]:
        # Only return if there's an exact match after sorting
        if score == 0: # This means edge_freq(H2) == fp2
            return H2, time.time() - t0 # Return time elapsed
    return None, time.time() - t0 # Return None and time if no exact match found in top_k

def backtrack_search(fp2, p):
    """Performs a depth-first search (backtracking) to construct a Hamiltonian path
    with 'p' vertices that matches the given frequency partition 'fp2'."""
    multiset = []
    for hop, cnt in fp2.items():
        multiset.extend([hop] * cnt) # Create a list of all required hop lengths

    # Path starts at 0, used keeps track of visited vertices
    # back counts backtracks, t0 is start time
    path, used, back, t0 = [0], {0}, 0, time.time()

    def dfs():
        nonlocal back # Allow modification of 'back' from outer scope
        
        # Base case: If path has 'p' vertices, a Hamiltonian path is found
        if len(path) == p:
            return True
        
        # Iterate through available hops in the multiset
        for i, hop in enumerate(list(multiset)): # Create a copy for safe iteration while modifying
            # Try two possible next vertices: current_vertex + hop and current_vertex - hop (modulo p)
            current_vertex = path[-1]
            
            # Remove the current hop from multiset for this branch
            rem = multiset.pop(i) 

            for nxt in sorted(((current_vertex + hop) % p, (current_vertex - hop) % p)): # Sort to potentially standardize search order
                if nxt not in used:
                    used.add(nxt)
                    path.append(nxt)
                    
                    if dfs(): # Recursive call
                        return True # Solution found
                    
                    # Backtrack: undo the choice
                    path.pop()
                    used.remove(nxt)
                    back += 1 # Increment backtrack counter
            
            # Backtrack: add the hop back to the multiset
            multiset.insert(i, rem) # Re-insert at original position

        return False # No solution found from this branch

    found = dfs()
    return (path if found else None), time.time() - t0, back

# --- Batch Job ---

def batch_extend(HP1, FP1_tuple, num_iterations=10, mode=2):
    """
    Extends an initial Hamiltonian Path (HP1) based on an evolving Frequency Partition (FP).
    Applies BHR check and various construction methods.
    """
    base_fp_list = list(FP1_tuple) # Use a list for mutability in evolving FP
    base_path = HP1.copy()
    results = []
    timestamp = datetime.now().strftime("%H%M%S")
    log_file = f"log_unified_{timestamp}.txt"

    print(f"\nStarting batch job. Log file: {log_file}")

    for i in range(num_iterations):
        fp_prev = parse_fp_tuple(tuple(base_fp_list)) # FP of the *previous* successful path
        p_prev = sum(fp_prev.values()) + 1 # Vertices in the *previous* path

        # --- Evolve FP for current iteration ---
        if mode == 1:
            idx = i % len(base_fp_list)
            base_fp_list[idx] += 1
        elif mode == 2:
            max_parts = p_prev // 2
            # Check for p/2 condition before appending a new part
            if len(base_fp_list) >= max_parts:
                print(f"⚠️ Iteration {i+1} skipped: FP has max ⌊p/2⌋ = {max_parts} for p = {p_prev}")
                # Log the skip before continuing
                t_stamp = datetime.now().strftime("%H:%M:%S")
                with open(log_file, "a", encoding="utf-8") as f:
                    f.write(f"--- Iteration {i+1} ---\n")
                    f.write(f"Timestamp: {t_stamp}\n")
                    f.write(f"p (vertices): {p_prev + 1} (would be)\n") # p_curr if it weren't skipped
                    f.write(f"Previous FP: {fp_prev}\n")
                    f.write(f"Evolved FP:  {parse_fp_tuple(tuple(base_fp_list + [1]))} (attempted)\n")
                    f.write(f"⚠️ Skipped: FP has max ⌊p/2⌋ = {max_parts} for p = {p_prev} (current p), cannot add new hop.\nResult: SKIPPED\n\n")
                continue # Skip this iteration
            base_fp_list.append(1) # Add a new hop of length 1
        else:
            print("❌ Invalid mode selected. Exiting batch job.")
            return # Exit batch_extend

        fp_curr = parse_fp_tuple(tuple(base_fp_list)) # FP for the *current* iteration
        p_curr = sum(fp_curr.values()) + 1 # Vertices for the *current* path
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
                f.write(f"Method used: BHR Fail\n") # Indicate BHR failure as method
                f.write(f"p (vertices): {p_curr}\n")
                f.write(f"Previous FP: {fp_prev}\n")
                f.write(f"Evolved FP:  {fp_curr}\n")
                f.write(f"❌ BHR failed at divisor {bad_div}: count = {bad_count} > {p_curr - bad_div}\nResult: SKIPPED\n\n")
            continue # Skip to next iteration if BHR fails
        else:
            print(f"✅ BHR passed for iteration {i+1}")

        # Attempt path construction methods
        H2 = None
        t_elapsed = 0
        backtracks = 0
        stage = "Unknown"

        # 1. Try Reuse-Insert
        H2, t_elapsed = reuse_insert(base_path, fp_curr)
        stage = "reuse-insert"

        if not H2:
            # 2. Try Greedy-Insert if Reuse failed
            H2, t_elapsed_greedy = greedy_insert(base_path, fp_curr)
            t_elapsed += t_elapsed_greedy # Add greedy time
            stage = "heuristic" if H2 else "backtrack" # Update stage based on H2 success

            if not H2:
                # 3. Try Backtrack Search if both inserts failed
                H2, t_elapsed_back, backtracks = backtrack_search(fp_curr, p_curr)
                t_elapsed += t_elapsed_back # Add backtrack time
                # Stage is already "backtrack"

        elapsed_formatted = round(t_elapsed, 6)

        # Log results
        with open(log_file, "a", encoding="utf-8") as f:
            f.write(f"--- Iteration {i+1} ---\n")
            f.write(f"Timestamp: {t_stamp}\n")
            f.write(f"Method used: {stage}\n")
            f.write(f"p (vertices): {p_curr}\n")
            f.write(f"Previous FP: {fp_prev}\n")
            f.write(f"Evolved FP:  {fp_curr}\n")
            f.write("✅ BHR passed\n") # BHR already passed check at this point

            if H2:
                f.write(f"✅ HP: {H2}\n")
                f.write(f"HP freq: {edge_freq(H2)}\n")
                f.write(f"🔙 Backtracks: {backtracks}\n")
                f.write(f"⏱ Time: {elapsed_formatted:.2f} sec\nResult: SUCCESS\n\n")
                results.append((i + 1, p_curr, H2, fp_curr, backtracks))
                base_path = H2 # Update base_path for next iteration
                print(f"✅ HP found (Method: {stage}, Time: {elapsed_formatted:.2f}s, Backtracks: {backtracks})")
            else:
                f.write("❌ Failed to construct HP\n")
                f.write(f"🔙 Backtracks: {backtracks}\n") # Backtracks could be 0 if not backtrack search or if it returned None quickly
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

    # Calculate actual edges in HP1 once for validation
    actual_total_edges_in_hp1 = len(HP1) - 1
    if actual_total_edges_in_hp1 < 1:
        print("❌ Error: Hamiltonian Path must have at least 2 vertices (1 edge).")
        exit()

    # Loop to allow correction of FP1 if it doesn't match HP1's edge count
    while True:
        FP1_tuple = get_int_tuple("Enter Frequency Partition Tuple")

        # --- VALIDATION CHECK ---
        parsed_fp1 = parse_fp_tuple(FP1_tuple)
        expected_total_edges_from_fp = sum(parsed_fp1.values())

        if expected_total_edges_from_fp != actual_total_edges_in_hp1:
            print("\n❌ Error: Mismatch between provided HP and FP.")
            print(f"   Your HP has {len(HP1)} vertices, implying {actual_total_edges_in_hp1} edges.")
            print(f"   Your FP implies a path with {expected_total_edges_from_fp} edges (sum of frequencies).")
            print("   Please re-enter your Frequency Partition Tuple to match the HP.")
            # The loop will continue, prompting for FP1_tuple again
        else:
            print(f"✅ FP matches HP (Total edges: {actual_total_edges_in_hp1}).")
            break # Exit the loop if validation passes
        # --- END VALIDATION CHECK ---

    num_iterations = get_iteration_count()

    batch_extend(HP1, FP1_tuple, num_iterations=num_iterations, mode=mode_choice)
