import os
import stat
import pandas as pd
import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError

from genbankqc import Genbank


class Metadata(Genbank.Genbank):
    def __init__(self, genbank):
        super().__init__(genbank)
        self.biosample_dir = os.path.join(
            self.genbank,
            "biosample_data",
        )
        self.sra_dir = os.path.join(
            self.genbank,
            "sra",
        )

    biosample_fields = [
        "geo_loc_name",
        "collection_date",
        "strain",
        "isolation_source",
        "host",
        "collected_by",
        "sample_type",
        "sample_name",
        "host_disease",
        "isolate",
        "host_health_state",
        "serovar",
        "env_biome",
        "env_feature",
        "ref_biomaterial",
        "env_material",
        "isol_growth_condt",
        "num_replicons",
        "sub_species",
        "host_age",
        "genotype",
        "host_sex",
        "serotype",
        "host_disease_outcome",
    ]

    @property
    def srs_df(self):
        from glob import glob
        srs_all = os.path.join(self.sra_dir, "srs_all.csv")
        if os.path.isfile(srs_all):
            os.remove(srs_all)
        csvs = glob(os.path.join(self.sra_dir, "*csv"))
        with open(srs_all, "a") as f:
            for csv in csvs:
                f.write(open(csv).readlines()[1])
        srs_df = pd.read_csv(srs_all, index_col=0, header=None)
        srs_df.columns = ["SRS"]
        srs_df.to_csv(srs_all)
        return srs_df

    @property
    def biosample_df(self):
        from glob import glob
        biosample_all = os.path.join(self.biosample_dir, "biosample_all.csv")
        if os.path.isfile(biosample_all):
            os.remove(biosample_all)
        csvs = glob(os.path.join(self.biosample_dir, "*csv"))
        with open(biosample_all, "a") as f:
            for csv in csvs:
                f.write(open(csv).readlines()[1])
        # dfs = (pd.read_csv(f, index_col=0) for f in csvs)
        biosample_df = pd.read_csv(biosample_all, index_col=0, header=None)
        biosample_df.columns = [
            "scientific_name",
            "SRA",
            "accession_id",
            "geo_loc_name",
            "collection_date",
            "strain",
            "isolation_source",
            "host",
            "collected_by",
            "sample_type",
            "sample_name",
            "host_disease",
            "isolate",
            "host_health_state",
            "serovar",
            "env_biome",
            "env_feature",
            "ref_biomaterial",
            "env_material",
            "isol_growth_condt",
            " num_replicons",
            " sub_species",
            " host_age",
            " genotype",
            " host_sex",
            " serotype",
            " host_disease_outcome",
        ]
        biosample_df.to_csv(
            os.path.join(self.biosample_dir, "biosample_all.csv"))
        return biosample_df

    @property
    def biosample_ids(self):
        ids = self.assembly_summary.biosample[
            self.assembly_summary.biosample.notnull()]
        return ids

    @property
    def sra_ids(self):
        ids = self.biosample_df.SRA[
            self.biosample_df.SRA.notnull()]
        return ids

    @property
    def biosample_xml(self):
        for f in os.listdir(self.biosample_dir):
            if f.endswith(".xml"):
                yield os.path.join(self.biosample_dir, f)

    @property
    def sra_xml(self):
        for f in os.listdir(self.sra_dir):
            if not f.endswith(".xml"):
                continue
            yield os.path.join(self.sra_dir, f)

    def commands_biosample(self):
        efetch_biosample = os.path.join(
            self.biosample_dir,
            "efetch_biosample.sh",
        )
        if not os.path.isdir(self.biosample_dir):
            os.mkdir(self.biosample_dir)
        if os.path.isfile(efetch_biosample):
            os.remove(efetch_biosample)
        ext = ".xml"
        with open(efetch_biosample, "a") as f:
            for i in self.biosample_ids:
                out = os.path.join(self.biosample_dir, i + ext)
                if os.path.isfile(out):
                    continue
                cmd = ("esearch -db biosample -query {} | "
                       "efetch -format docsum "
                       "> {}\n".format(i, out))
                f.write(cmd)
        os.chmod(efetch_biosample, stat.S_IRWXU)

    def commands_sra(self):
        efetch_sra = os.path.join(
            self.sra_dir,
            "efetch_sra.sh",
        )
        if not os.path.isdir(self.sra_dir):
            os.mkdir(self.sra_dir)
        if os.path.isfile(efetch_sra):
            os.remove(efetch_sra)
        ext = ".xml"
        with open(efetch_sra, "a") as f:
            for i in self.sra_ids:
                out = os.path.join(self.sra_dir, i + ext)
                if os.path.isfile(out):
                    continue
                cmd = ("esearch -db sra -query {} | "
                       "efetch -format docsum "
                       "> {}\n".format(i, out))
                f.write(cmd)
        os.chmod(efetch_sra, stat.S_IRWXU)

    def parse_biosample_xml(self):
        for f in self.biosample_xml:
            out = os.path.splitext(os.path.basename(f))[0] + '.csv'
            out = os.path.join(self.biosample_dir, out)
            if os.path.isfile(out):
                continue
            try:
                tree = ET.ElementTree(file=f)
            except ParseError:
                continue
            accession = tree.find("DocumentSummary/Accession").text
            sra = 'DocumentSummary/SampleData/BioSample/Ids/Id/[@db="SRA"]'
            try:
                sra = tree.find(sra).text
            except AttributeError:
                sra = None
            try:
                gca = self.assembly_summary.index[
                    self.assembly_summary.biosample == accession].format()[0]
            except IndexError:
                continue
            scientific_name = self.assembly_summary.loc[gca].scientific_name
            df = pd.DataFrame()
            df.loc[accession, "scientific_name"] = scientific_name
            df.loc[accession, "SRA"] = sra
            df.loc[accession, "accession_id"] = gca
            for name in self.biosample_fields:
                xp = (
                    'DocumentSummary/SampleData/BioSample/Attributes/Attribute'
                    '[@harmonized_name="{}"]'.format(name)
                )
                try:
                    attrib = tree.find(xp).text
                    df.loc[accession, name] = attrib
                except AttributeError:
                    df.loc[accession, name] = None
            df.to_csv(out)

    def parse_sra_xml(self):
        for f in self.sra_xml:
            name = os.path.splitext(os.path.basename(f))[0]
            out = os.path.join(self.sra_dir, name + '.csv')
            df = pd.DataFrame()
            if os.path.isfile(out):
                continue
            try:
                tree = ET.ElementTree(file=f)
            except ParseError:
                continue
            elements = tree.iterfind("DocumentSummary/Runs/Run/[@acc]")
            srs_accessions = []
            for el in elements:
                    items = el.items()
                    acc = [i[1] for i in items if i[0] == 'acc']
                    acc = acc[0]
                    srs_accessions.append(acc)
            srs_accessions = ','.join(srs_accessions)
            df.loc[name, "SRS"] = srs_accessions
            df.to_csv(out)

    def metadata(self):
        biosample_df = self.biosample_df
        for s in self.species:
            biosamples_out = os.path.join(s.qc_dir, "biosamples.csv")
            srs_out = os.path.join(s.qc_dir, "SRS.csv")
            if os.path.isfile(biosamples_out):
                continue
            elif os.path.isfile(srs_out):
                continue
            # df = pd.DataFrame(columns=[''])
            gcas = s.accession_ids
            biosamples = self.assembly_summary.loc[gcas].biosample
            try:
                species_biosamples = biosample_df.loc[biosamples.values]
            except KeyError:
                continue
            try:
                srs = self.srs_df[self.srs_df.SRS.notnull()]
                srs = srs.loc[species_biosamples.SRA.tolist()]
            except KeyError:
                continue
            species_biosamples.to_csv(biosamples_out)
            srs.to_csv(srs_out)
            print(s.species)
