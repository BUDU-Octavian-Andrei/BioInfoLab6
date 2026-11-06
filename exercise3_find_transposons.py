"""Exercise 3, part 2: detect the start and end of transposable elements in a DNA sequence.

A DNA transposon starts and ends with inverted repeats (the right end is the reverse
complement of the left end), and inserting it duplicates the target site, so the
element is flanked by two copies of the same short sequence (direct repeats):

    ... ATGCA  ACGTTGCAAG ...body... CTTGCAACGT  ATGCA ...
        direct inverted repeat       inverted    direct
        repeat                       repeat      repeat

The detector looks for exactly that pattern, without knowing where anything was inserted.
"""
import argparse
from pathlib import Path
from typing import NamedTuple

MIN_INVERTED_REPEAT = 8
MIN_DIRECT_REPEAT = 4
MAX_DIRECT_REPEAT = 10
MIN_ELEMENT_LENGTH = 30
MAX_ELEMENT_LENGTH = 150

DATA_DIR = Path(__file__).parent / "data"
FASTA_PATH = DATA_DIR / "artificial_sequence.fasta"
ANSWER_PATH = DATA_DIR / "artificial_elements.tsv"

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


class Element(NamedTuple):
    start: int  # 0-based, first base of the left inverted repeat
    end: int  # 0-based, last base of the right inverted repeat
    inverted_repeat: str
    direct_repeat: str

    @property
    def length(self):
        return self.end - self.start + 1


def reverse_complement(dna):
    return "".join(COMPLEMENT[base] for base in reversed(dna))


def inverted_repeat_length(sequence, start, end):
    """How many bases from start onwards pair with the bases from end backwards."""
    length = 0
    while start + length < end - length and sequence[start + length] == COMPLEMENT[sequence[end - length]]:
        length += 1
    return length


def direct_repeat_length(sequence, start, end):
    """Length of the longest identical pair of sequences just before start and just after end."""
    for length in range(MAX_DIRECT_REPEAT, MIN_DIRECT_REPEAT - 1, -1):
        if start - length >= 0 and sequence[start - length:start] == sequence[end + 1:end + 1 + length]:
            return length
    return 0


def find_transposons(sequence):
    candidates = []
    for start in range(MIN_DIRECT_REPEAT, len(sequence)):
        last_end = min(start + MAX_ELEMENT_LENGTH, len(sequence) - MIN_DIRECT_REPEAT)
        for end in range(start + MIN_ELEMENT_LENGTH - 1, last_end):
            inverted = inverted_repeat_length(sequence, start, end)
            if inverted < MIN_INVERTED_REPEAT:
                continue
            direct = direct_repeat_length(sequence, start, end)
            if direct < MIN_DIRECT_REPEAT:
                continue
            candidates.append(Element(start, end, sequence[start:start + inverted], sequence[start - direct:start]))

    # Overlapping candidates are alternatives for the same element: keep the best supported.
    # The inverted repeat ranks first, because shifting the ends inwards can trade a base
    # of it for a longer direct repeat when the flanking host bases happen to match.
    candidates.sort(key=lambda e: (len(e.inverted_repeat), len(e.direct_repeat)), reverse=True)
    elements = []
    for candidate in candidates:
        if all(candidate.end < kept.start or candidate.start > kept.end for kept in elements):
            elements.append(candidate)
    return sorted(elements)


def load_fasta(path):
    lines = Path(path).read_text().splitlines()
    return lines[0].lstrip(">"), "".join(lines[1:]).upper()


def load_answer():
    """The (start, end) pairs the generator wrote, 0-based, or None if there is no answer file."""
    if not ANSWER_PATH.exists():
        return None
    rows = ANSWER_PATH.read_text().splitlines()[1:]
    return [(int(row.split("\t")[0]) - 1, int(row.split("\t")[1]) - 1) for row in rows]


def sequence_map(sequence, elements, width=60):
    """The sequence in rows, with a line under each row marking the parts of every element."""
    marks = [" "] * len(sequence)
    for element in elements:
        inverted = len(element.inverted_repeat)
        direct = len(element.direct_repeat)
        for i in range(element.start, element.end + 1):
            marks[i] = "-"
        for i in range(inverted):
            marks[element.start + i] = ">"
            marks[element.end - i] = "<"
        for i in range(direct):
            marks[element.start - 1 - i] = "="
            marks[element.end + 1 + i] = "="

    rows = []
    for row_start in range(0, len(sequence), width):
        rows.append(f"{row_start + 1:>5}  {sequence[row_start:row_start + width]}")
        rows.append(f"{'':>5}  {''.join(marks[row_start:row_start + width])}")
    return "\n".join(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fasta", nargs="?", default=FASTA_PATH, help="sequence to search")
    args = parser.parse_args()

    description, sequence = load_fasta(args.fasta)
    print(f"Sequence: {description} ({len(sequence)} nt)\n")

    elements = find_transposons(sequence)

    # Positions are printed 1-based.
    print(f"{'#':>2}  {'Start':>5}  {'End':>5}  {'Length':>6}  {'Inverted repeat':<16}  Direct repeat")
    for number, element in enumerate(elements, start=1):
        print(f"{number:>2}  {element.start + 1:>5}  {element.end + 1:>5}  {element.length:>6}  "
              f"{element.inverted_repeat:<16}  {element.direct_repeat}")
    print(f"\nFound {len(elements)} transposable elements.\n")

    print(sequence_map(sequence, elements))
    print("\n= direct repeat (duplicated target site)   > < inverted repeats   - element body")

    answer = load_answer()
    if answer is not None and Path(args.fasta) == FASTA_PATH:
        found = [(element.start, element.end) for element in elements]
        if found == answer:
            print(f"\nAll positions match the {len(answer)} elements the generator inserted.")
        else:
            print(f"\nMismatch with the generator's answer file:")
            print(f"  inserted: {[(start + 1, end + 1) for start, end in answer]}")
            print(f"  found:    {[(start + 1, end + 1) for start, end in found]}")


if __name__ == "__main__":
    main()
