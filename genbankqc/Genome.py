import os.path
import re
import subprocess

from logbook import Logger
from retrying import retry
import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError
from collections import defaultdict
from genbankqc import Metadata

import pandas as pd
from Bio import SeqIO


class Genome:
    def __init__(self, genome, assembly_summary=None):
        """
        :param genome: Path to genome
        :returns: Path to genome and name of the genome
        :rtype:
        """
        self.path = os.path.abspath(genome)
        self.species_dir, self.fasta = os.path.split(self.path)
        self.name = os.path.splitext(self.fasta)[0]
        self.log = Logger(self.name)
        self.qc_dir = os.path.join(self.species_dir, "qc")
        self.msh = os.path.join(self.qc_dir, self.name + ".msh")
        self.stats_path = os.path.join(self.qc_dir, self.name + '.csv')
        if os.path.isfile(self.stats_path):
            self.stats = pd.read_csv(self.stats_path, index_col=0)
        self.assembly_summary = assembly_summary
        self.metadata = defaultdict(lambda: 'missing')
        self.xml = defaultdict(lambda: 'missing')
        try:
            self.accession_id = re.match('GCA_.*\.\d', self.name).group()
            self.metadata["accession"] = self.accession_id
        except AttributeError:
            # Raise custom exception
            self.accession_id = "missing"
            self.log.error("Invalid accession ID")
            self.log.exception()
        if isinstance(self.assembly_summary, pd.DataFrame):
            try:
                biosample = assembly_summary.loc[self.accession_id].biosample
                self.metadata["biosample_id"] = biosample
            except (AttributeError, KeyError):
                self.log.info("Unable to get biosample ID")
        self.log.info("Instantiated")

    def get_contigs(self):
        """
        Return a list of of Bio.Seq.Seq objects for fasta and calculate
        the total the number of contigs.
        """
        try:
            self.contigs = [seq.seq for seq in SeqIO.parse(self.path, "fasta")]
            self.count_contigs = len(self.contigs)
            self.log.info("Contigs: {}".format(self.count_contigs))
        except UnicodeDecodeError:
            self.log.exception()

    def get_assembly_size(self):
        """Calculate the sum of all contig lengths"""
        # TODO: map or reduce might be more elegant here
        self.assembly_size = sum((len(str(seq)) for seq in self.contigs))
        self.log.info("Assembly Size: {}".format(self.assembly_size))

    def get_unknowns(self):
        """Count the number of unknown bases, i.e. not [ATCG]"""
        # TODO: Would it be useful to allow the user to define p?
        p = re.compile("[^ATCG]")
        self.unknowns = sum((len(re.findall(p, str(seq)))
                             for seq in self.contigs))
        self.log.info("Unknowns: {}".format(self.unknowns))

    def get_distance(self, dmx_mean):
        self.distance = dmx_mean.loc[self.name]
        self.log.info("Distance: {}".format(self.distance))

    def sketch(self):
        cmd = "mash sketch '{}' -o '{}'".format(self.path, self.msh)
        if os.path.isfile(self.msh):
            self.log.info("Sketch file already exists")
        else:
            subprocess.Popen(cmd, shell="True", stderr=subprocess.DEVNULL).wait()
            self.log.info("Sketch file created")

    def get_stats(self, dmx_mean):
        if not os.path.isfile(self.stats_path):
            self.get_contigs()
            self.get_assembly_size()
            self.get_unknowns()
            self.get_distance(dmx_mean)
            data = {"contigs": self.count_contigs,
                    "assembly_size": self.assembly_size,
                    "unknowns": self.unknowns,
                    "distance": self.distance}
            self.stats = pd.DataFrame(data, index=[self.name])
            self.stats.to_csv(self.stats_path)
            self.log.info("Generated stats and wrote to disk")

    one_minute = 60000

    # Retry 3 times over a period of 3 minutes max,
    # waiting five seconds in between retries
    @retry(stop_max_attempt_number=3,
           stop_max_delay=10000,
           wait_fixed=100)
    def efetch(self, db):
        """
        Use NCBI's efetch tools to get xml for genome's biosample id or SRA id
        """
        if db == "biosample":
            db_id = db + "_id"
        elif db == "sra":
            db_id = db + "_id"
        cmd = ("esearch -db {} -query {} | "
               "efetch -format docsum".format(db, self.metadata[db_id]))
        # Make efetch timeout and retry after 30 seconds
        time_limit = 30
        if self.metadata[db_id] is not 'missing':
            try:
                p = subprocess.run(cmd, shell="True",
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.DEVNULL,
                                   timeout=time_limit)
                xml = p.stdout
                self.xml[db] = xml
                self.log.info("{} XML downloaded".format(db))
            except subprocess.TimeoutExpired:
                self.log.error("Retrying efetch after timeout")
                raise subprocess.TimeoutExpired(cmd, time_limit)
            except Exception:
                self.log.error(db)
                self.log.exception()

    def parse_biosample(self):
        """
        Get what we need to get out of the xml returned by efetch("biosample")
        Including the SRA ID and fields of interest as defined in
        Metadata.biosample_fields
        """
        try:
            tree = ET.fromstring(self.xml["biosample"])
            sra = tree.find('DocumentSummary/SampleData/BioSample/Ids/Id/[@db="SRA"]')
            self.log.info("Parsed biosample XML")
            try:
                self.metadata["sra_id"] = sra.text
            except AttributeError:
                self.metadata["sra_id"] = "missing"
            for name in Metadata.Metadata.biosample_fields:
                xp = ('DocumentSummary/SampleData/BioSample/Attributes/'
                      'Attribute/[@harmonized_name="{}"]'.format(name))
                attrib = tree.find(xp)
                try:
                    self.metadata[name] = attrib.text
                except AttributeError:
                    self.metadata[name] = "missing"

        except ParseError:
            self.log.error("Parse error for biosample XML")

    def parse_sra(self):
        try:
            tree = ET.fromstring(self.xml["sra"])
            elements = tree.iterfind("DocumentSummary/Runs/Run/[@acc]")
            self.log.info("Parsed SRA XML")
            srr_accessions = []
            for el in elements:
                    items = el.items()
                    acc = [i[1] for i in items if i[0] == 'acc']
                    acc = acc[0]
                    srr_accessions.append(acc)
            self.metadata["srr_accessions"] = ','.join(srr_accessions)
        except ParseError:
            self.log.error("Parse error for SRA XML")

    def get_metadata(self):
        self.efetch("biosample")
        self.parse_biosample()
        self.efetch("sra")
        self.parse_sra()


# make sure Genome reads in the assembly summary here
def sketch_genome(path):
    genome = Genome(path)
    genome.sketch()


def mp_stats(path, dmx_mean):
    genome = Genome(path)
    genome.get_stats(dmx_mean)
    return genome.stats
