from hopcroftkarp import HopcroftKarp
from collections import defaultdict
import copy

VALID_PAIRS = set(["AU", "UA", "CG", "GC", "GU", "UG"])

def parse_secondary_structure(struc, get_pair_info = False):
    """
    Parse an RNA secondary structure string and return the indices of paired and unpaired nucleotides.

    Args:
        struc (str): A string representing RNA secondary structure, where '(' and ')' denote paired nucleotides
                     and '.' denotes unpaired nucleotides.

    Returns:
        tuple: A tuple containing two lists:
            - A list of tuples representing the indices of paired nucleotides in the structure string.
            - A list of indices representing the indices of unpaired nucleotides in the structure string.

    Example:
        >>> parse_secondary_structure('((..((...)).))')
        ([(0, 11), (4, 9), (5, 8)], [2, 3, 6, 7, 12, 13])

    If the input string contains unbalanced parentheses, the function returns None and prints an error message.
    """
    stack1, stack2, stack3, stack4 = [], [], [], []
    paired_indices = []  # list of tuples: [(i1, j1), (i2, j2), ...]
    unpaired_indices = []  # list of indices: [i1, i2, ...]

    try:
        for i, x in enumerate(struc):
            if x == "(":
                stack1.append(i)
            elif x == "[":
                stack2.append(i)
            elif x == "{":
                stack3.append(i)
            elif x == "<":
                stack4.append(i)
            elif x == ")":
                u, v = stack1.pop(), i
                assert struc[u] == '(' and struc[v] == ')', "Unbalanced parenthesis in structure string"
                paired_indices.append((u, v, ('(', ')')) if get_pair_info else (u, i))
            elif x == "]":
                # paired_indices.append((stack2.pop(), i, ('[', ']')) if get_pair_info else (stack2.pop(), i))
                u, v = stack2.pop(), i
                assert struc[u] == '[' and struc[v] == ']', "Unbalanced parenthesis in structure string"
                paired_indices.append((u, v, ('[', ']')) if get_pair_info else (u, v))
            elif x == "}":
                # paired_indices.append((stack3.pop(), i, ('{', '}')) if get_pair_info else (stack3.pop(), i))
                u, v = stack3.pop(), i
                assert struc[u] == '{' and struc[v] == '}', "Unbalanced parenthesis in structure string"
                paired_indices.append((u, v, ('{', '}')) if get_pair_info else (u, v))
            elif x == ">":
                # paired_indices.append((stack4.pop(), i, ('<', '>')) if get_pair_info else (stack4.pop(), i))
                u, v = stack4.pop(), i
                assert struc[u] == '<' and struc[v] == '>', "Unbalanced parenthesis in structure string"
                paired_indices.append((u, v, ('<', '>')) if get_pair_info else (u, v))
            elif x == ".":
                unpaired_indices.append(i)
    except Exception as _:
        print("[Error] Unbalanced parenthesis in structure string")
        return None

    if stack1 or stack2 or stack3 or stack4:
        print("[Error] Unbalanced parenthesis in structure string")

    return paired_indices, unpaired_indices


def evaluate(pred_struc = "", gold_struc = "", allow_slip=False, gold_paired_pos_tuple = None, gold_unpaired_pos = None, pred_paired_pos_tuple = None, pred_unpaired_pos = None, no_assert = False):
    """
    Evaluates the predicted RNA secondary structure against a gold structure.

    Args:
        pred_struc (str): A string representing the predicted RNA secondary structure, where '(' and ')' denote paired
                          nucleotides and '.' denotes unpaired nucleotides.
        gold_struc (str): A string representing the gold RNA secondary structure, where '(' and ')' denote paired
                          nucleotides and '.' denotes unpaired nucleotides.
        allow_slip (bool): A boolean indicating whether to allow one-nucleotide slips in the predicted structure when
                           calculating precision and sensitivity. If True, a predicted paired nucleotide position is
                           considered correct if it is within one nucleotide of a gold paired nucleotide position.

    Returns:
        tuple: A tuple containing four float values:
            - Precision: The fraction of predicted paired nucleotide positions that are correct.
            - Sensitivity: The fraction of gold paired nucleotide positions that are correctly predicted.
            - F1 score: The harmonic mean of precision and sensitivity.
            - Structural distance: The difference between the length of the gold structure string and twice the number
                                   of common paired nucleotide positions plus the number of common unpaired nucleotide
                                   positions. This represents how close the predicted structure is to the gold structure.

    Raises:
        AssertionError: If the length of the predicted and gold structure strings are not equal.

    Example:
        >>> evaluate('((..((...)).))', '((..((...)).))')
        (1.0, 1.0, 1.0, 0)
        >>> evaluate('((..((...)).))', '((..(.()..).))')
        (0.75, 0.75, 0.75, 4)
    """
    if not no_assert:
        assert len(pred_struc) == len(
            gold_struc
        ), "Length of predicted and gold structure strings must be equal\nPredicted Length: {}\nGold Length: {}".format(
            len(pred_struc), len(gold_struc)
        )

    if gold_paired_pos_tuple is None or gold_unpaired_pos is None:
        gold_paired_pos_tuple, gold_unpaired_pos = parse_secondary_structure(gold_struc)
    if pred_paired_pos_tuple is None or pred_unpaired_pos is None:
        pred_paired_pos_tuple, pred_unpaired_pos = parse_secondary_structure(pred_struc)


    if allow_slip:
        graph = defaultdict(set)
        for (i, j) in pred_paired_pos_tuple:
            for (x,y) in [(i, j), (i-1,j), (i+1,j), (i,j-1), (i,j+1)]:
                if (x, y) in gold_paired_pos_tuple:
                    graph[(i,j)].add((str(x), str(y)))
        
        matching = HopcroftKarp(graph).maximum_matching()
        # only select values that are tuple of string
        common_paired = set([(int(i), int(j)) for (i, j) in matching.values() if isinstance(i, str)])

        all_paired_pos = set()
        pred_new_unpaired_pos = []
        for (i, j) in common_paired:
            all_paired_pos.add(i)
            all_paired_pos.add(j)
        for (i, j) in pred_paired_pos_tuple:
            if (i, j) not in matching:
                all_paired_pos.add(i)
                all_paired_pos.add(j)

        # get new unpaired pos
        for i in range(len(pred_struc)):
            if i not in all_paired_pos and pred_struc[i] != '*':
                pred_new_unpaired_pos.append(i)
        # get common unpaired pos
        common_unpaired = set(pred_new_unpaired_pos).intersection(gold_unpaired_pos)
    else:
        common_paired = set(pred_paired_pos_tuple).intersection(gold_paired_pos_tuple)
        common_unpaired = set(pred_unpaired_pos).intersection(gold_unpaired_pos)


    precision = len(common_paired) / (len(pred_paired_pos_tuple) + 1e-6)
    sensitivity = len(common_paired) / (len(gold_paired_pos_tuple) + 1e-6)
    f1 = 2 * precision * sensitivity / (precision + sensitivity)
    structural_distance = len(gold_struc) - (2 * len(common_paired) + len(common_unpaired))

    return precision, sensitivity, f1, structural_distance




def get_alignment_to_sequence_mapping(aligned_sequence):
    """
    Returns a dictionary mapping the indices of the aligned sequence to the indices of the unaligned sequence.

    Args:
        aligned_sequence (str): A string representing the aligned sequence, where '-' denotes a gap.

    Returns:
        dict: A dictionary mapping the indices of the aligned sequence to the indices of the unaligned sequence.

    Example:
        >>> get_alignment_to_sequence_mapping('AUCG-AUCG')
        {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 4, 6: 5, 7: 6}
    """
    mapping = {}
    j = 0
    for i, x in enumerate(aligned_sequence):
        if x != "-":
            mapping[i] = j
            j += 1
        else:
            mapping[i] = j - 1
    return mapping


def get_bpp_matrix(bpp_file_path, threshold=0.3):
    bpp_matrix = {}
    with open(bpp_file_path, 'r') as f:
        lines = f.readlines()
    seq_length = -1
    for line in lines:
        split = line.split(' ')
        if len(split) == 3 and split[0] not in ['(' or '.']:
            i, j, prob = int(split[0]) - 1, int(split[1]) - 1, float(split[2])
            if prob >= threshold:
                bpp_matrix[i, j] = prob
            seq_length = max(seq_length, i + 1, j + 1)
    return bpp_matrix, seq_length


def pairs_to_struc(pairs, seq_length):
    struc = ['.' for _ in range(seq_length)]
    for i, j in pairs:
        struc[i] = '('
        struc[j] = ')'
    return ''.join(struc)

def map_consns_struc_to_aln_seq(consns_struc, aln_seq):
    a2s_map = get_alignment_to_sequence_mapping(aln_seq)
    seq = aln_seq.replace("-", "")  # ungapped sequence

    # get the structure corresponding to the unaligned sequence
    struc = ["."] * len(seq)

    for i, j, (b1, b2) in parse_secondary_structure(consns_struc, True)[0]:
        if aln_seq[i] + aln_seq[j] not in VALID_PAIRS:
            continue
        if seq[a2s_map[i]] + seq[a2s_map[j]] not in VALID_PAIRS:
            continue

        struc[a2s_map[i]] = b1
        struc[a2s_map[j]] = b2

    return "".join(struc), seq

def map_consns_bpp_to_aln_seq(consns_bpp, aln_seq, threshold=0.01):
    a2s_map = get_alignment_to_sequence_mapping(aln_seq)
    seq = aln_seq.replace("-", "")  # ungapped sequence

    # get the structure corresponding to the unaligned sequence
    bpp = defaultdict(lambda: defaultdict(float))

    for i in range(len(aln_seq)):
        for j in consns_bpp[i].keys():
            if aln_seq[i] + aln_seq[j] not in VALID_PAIRS:
                continue
            if seq[a2s_map[i]] + seq[a2s_map[j]] not in VALID_PAIRS:
                continue
            if consns_bpp[i][j] >= threshold:
                bpp[a2s_map[i]][a2s_map[j]] = consns_bpp[i][j]
    
    return bpp, seq

    # for i, j, prob in consns_bpp:
    #     if aln_seq[i] + aln_seq[j] not in VALID_PAIRS:
    #         continue
    #     if seq[a2s_map[i]] + seq[a2s_map[j]] not in VALID_PAIRS:
    #         continue

    #     bpp[a2s_map[i]][a2s_map[j]] = prob

    # return bpp, seq


def get_struc_from_file(file_path, backward_search=False):
    """
    Get the structure string from a file.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in (lines[::-1] if backward_search else lines):
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] in ['.', '(']:
                return line
            
def get_seq_from_file(file_path, backward_search=False):
    """
    Get the sequence string from a file.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate((lines[::-1] if backward_search else lines)):
            header = lines[i-1].strip() if i > 0 else ""
            line = line.strip()
            if line[0] in set(['A', 'U', 'C', 'G', '-']):
                return header, line
    
def get_bpp_from_file(file_path, threshold=0.01, zero_indexed=False):
    """
    Get the bpp matrix from a file.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    bpp_matrix = defaultdict(lambda: defaultdict(float))
    seq_length = -1
    for line in lines:
        split = line.split(' ')
        if len(split) == 3 and len(split[0]) < 100 and float(split[2]) >= threshold:
            if zero_indexed:
                i, j, prob = int(split[0]), int(split[1]), float(split[2])
            else:
                i, j, prob = int(split[0]) - 1, int(split[1]) - 1, float(split[2])
            bpp_matrix[i][j] = prob
            seq_length = max(seq_length, j + 1)
    return bpp_matrix, seq_length

def convert_bpp_to_list(bpp_matrix, seq_length, threshold=0.01):
    consns_bpp = []
    for i in range(seq_length):
        for j in bpp_matrix[i].keys():
            if bpp_matrix[i][j] >= threshold:
                consns_bpp.append((i, j, bpp_matrix[i][j]))
    return consns_bpp


def get_unpaired_probs(paired_probs, seq_length):
    unpaired_probs = {}

    for i in range(seq_length):
        unpaired_probs[i] = 1.00

    for l in paired_probs.keys():
        for r, p in paired_probs[l].items():
            unpaired_probs[l] -= p
            unpaired_probs[r] -= p

    return unpaired_prob
    