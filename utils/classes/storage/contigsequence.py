from Bio.Seq import Seq


class ContigSequence:
    """ Stores the current contig's genomic sequence, read from the fasta format """
    def __init__(self, sequence: str = ""):
        self.sequence: str = sequence

    def get_sequence(self, start: int, stop: int) -> str:
        """ Returns the sequence after converting positions in O-base to slice the string properly """
        return self.sequence[start - 1: stop].lower() if self.sequence else ""

    def get_gc(self, list_positions: list[tuple[str, str, int, int]]) -> float | str:
        """ Computes the percent GC from all segments in a list of position tuples
        (contig:str, strand:str, start:int, stop:int) """
        seq = ""
        if self.sequence:
            for contig, _, start, stop in list_positions:
                if start <= stop:
                    seq += self.sequence[start - 1: stop].lower()
        return self._get_percent_gc(seq)

    @staticmethod
    def _get_percent_gc(dna: str) -> float | str:
        """ Returns the % GC of a DNA sequence (only IUPAC-compatible characters), or 0.0 if the sequence is empty """
        if dna:
            dna: str = dna.lower()
            g: float
            c: float
            s: float
            g, c, s = dna.count("g"), dna.count("c"), dna.count("s")
            return round((g + c + s) * 100 / len(dna), 2)
        return "-"

    def get_reverse_complement(self, start: int, stop: int):
        """ Returns the reverse-complementary sequence for the provided start, stop interval """
        return self._reverse_complement(self.sequence[start - 1: stop].lower()) if self.sequence else ""

    @staticmethod
    def _reverse_complement(dna: str) -> str:
        """ Returns the reverse-complement of a DNA sequence (only IUPAC-compatible characters) """
        return str(Seq(dna).reverse_complement())
