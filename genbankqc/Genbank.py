import os
import pandas as pd
from logbook import Logger

from genbankqc import Species

taxdump_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"


class Genbank:
    def __init__(self, path):
        """
        GenBank
        """
        self.path = os.path.abspath(path)
        self.info_dir = os.path.join(self.path, ".info")
        self.assembly_summary_path = os.path.join(self.info_dir, "assembly_summary.txt")
        self.log = Logger("GenBank")
        if not os.path.isdir(self.path):
            os.mkdir(self.path)
        if not os.path.isdir(self.info_dir):
            os.mkdir(self.info_dir)
        try:
            self.assembly_summary = pd.read_csv(self.assembly_summary_path, sep="\t", index_col=0)
        except FileNotFoundError:
            self.assembly_summary = pd.read_csv(assembly_summary_url, sep="\t",
                                                index_col=0, skiprows=1)
            self.assembly_summary.to_csv(self.assembly_summary_path, sep="\t")
            self.log.info("Downloaded assembly_summary.txt")

    @property
    def species(self):
        for d in os.listdir(self.path):
            species_path = os.path.join(self.path, d)
            if not os.path.isdir(species_path):
                continue
            fastas = len([f for f in os.listdir(species_path) if f.endswith('fasta')])
            if fastas < 10:
                self.log.info("Not enough genomes for {}".format(d))
                continue
            yield Species.Species(species_path, assembly_summary=self.assembly_summary)

    def qc(self):
        for species in self.species:
            species.qc()
