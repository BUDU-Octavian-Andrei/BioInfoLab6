"""Exercise 3, part 1: build an artificial DNA sequence holding 3-4 transposable elements.

Each element is inserted the way a DNA transposon inserts itself. The host is cut with a
staggered cut at the target site, the element goes in, and filling in the gaps copies the
target site on both sides:

    host    ... TAC ATGCA CAG ...
    after   ... TAC ATGCA [ inverted repeat ... body ... inverted repeat ] ATGCA CAG ...
"""
import argparse
import random

from exercise3_find_transposons import ANSWER_PATH, COMPLEMENT, DATA_DIR, FASTA_PATH, Element, reverse_complement

MIN_LENGTH = 200
MAX_LENGTH = 400
MIN_ELEMENTS = 3
MAX_ELEMENTS = 4
INVERTED_REPEAT_LENGTH = (8, 12)
BODY_LENGTH = (20, 40)
TARGET_SITE_LENGTH = (4, 6)
MIN_GAP = 10  # host bases kept before, between and after the insertions


def random_dna(length):
    return "".join(random.choices("ACGT", k=length))


def make_element():
    inverted_repeat = random_dna(random.randint(*INVERTED_REPEAT_LENGTH))
    body = random_dna(random.randint(*BODY_LENGTH))
    return inverted_repeat + body + reverse_complement(inverted_repeat), inverted_repeat


def make_target_site(length):
    """A target site whose first and last bases don't pair.

    If they paired, the copies on either side of the element would extend its inverted
    repeats by one base, and the element's true ends could not be told from the sequence.
    """
    while True:
        site = random_dna(length)
        if site[-1] != COMPLEMENT[site[0]]:
            return site


def build_sequence():
    """Returns the sequence and the inserted elements, with their positions in it."""
    count = random.randint(MIN_ELEMENTS, MAX_ELEMENTS)
    elements = [make_element() for _ in range(count)]
    site_lengths = [random.randint(*TARGET_SITE_LENGTH) for _ in range(count)]

    # Each insertion adds the element plus a second copy of its target site.
    inserted = sum(len(element) for element, _ in elements) + sum(site_lengths)
    host_needed = (count + 1) * MIN_GAP + sum(site_lengths)
    total_length = random.randint(max(MIN_LENGTH, inserted + host_needed), MAX_LENGTH)
    host = random_dna(total_length - inserted)

    # Share the spare host bases out at random between the gaps around the insertions.
    spare = len(host) - host_needed
    cuts = sorted(random.randint(0, spare) for _ in range(count))
    extra_gaps = [b - a for a, b in zip([0] + cuts, cuts + [spare])]

    sequence = ""
    inserted_elements = []
    position = 0  # where the next target site starts in the host
    previous = 0  # host bases up to here are already in the sequence
    for (element, inverted_repeat), site_length, extra_gap in zip(elements, site_lengths, extra_gaps):
        position += MIN_GAP + extra_gap
        site = make_target_site(site_length)
        host = host[:position] + site + host[position + site_length:]
        sequence += host[previous:position + site_length]
        inserted_elements.append(Element(len(sequence), len(sequence) + len(element) - 1, inverted_repeat, site))
        sequence += element + site
        position += site_length
        previous = position
    sequence += host[previous:]

    return sequence, inserted_elements


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, help="repeat the same sequence")
    args = parser.parse_args()
    random.seed(args.seed)

    sequence, elements = build_sequence()

    DATA_DIR.mkdir(exist_ok=True)
    description = f"artificial sequence with {len(elements)} transposable elements"
    lines = [sequence[i:i + 70] for i in range(0, len(sequence), 70)]
    FASTA_PATH.write_text(f">{description}\n" + "\n".join(lines) + "\n")

    # Positions are written 1-based, so they read the same as the detector's output.
    rows = ["start\tend\tinverted_repeat\tdirect_repeat"]
    rows += [f"{e.start + 1}\t{e.end + 1}\t{e.inverted_repeat}\t{e.direct_repeat}" for e in elements]
    ANSWER_PATH.write_text("\n".join(rows) + "\n")

    print(f"Built an {description} ({len(sequence)} nt).\n")
    print(f"{'#':>2}  {'Start':>5}  {'End':>5}  {'Length':>6}  {'Inverted repeat':<16}  Target site")
    for number, element in enumerate(elements, start=1):
        print(f"{number:>2}  {element.start + 1:>5}  {element.end + 1:>5}  {element.length:>6}  "
              f"{element.inverted_repeat:<16}  {element.direct_repeat}")
    print(f"\nSaved {FASTA_PATH} and the answer key {ANSWER_PATH}")


if __name__ == "__main__":
    main()
