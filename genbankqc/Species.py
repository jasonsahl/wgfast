import os
import re
import pickle
import functools

import logbook

from subprocess import DEVNULL, Popen
from pathos.multiprocessing import ProcessingPool as Pool

import pandas as pd

from ete3 import Tree
from genbankqc import Genome
# Figure out how to supress error output from this

os.environ['QT_QPA_PLATFORM'] = 'offscreen'


class Species:
    def __init__(self, path, max_unknowns=200, contigs=3.0, assembly_size=3.0,
                 mash=3.0, assembly_summary=None, processes=1):
        """Represents a collection of genomes in `path`

        :param path: Path to the directory of related genomes you wish to analyze.
        :param max_unknowns: Number of allowable unknown bases, i.e. not [ATCG]
        :param contigs: Acceptable deviations from median number of contigs
        :param assembly_size: Acceptable deviations from median assembly size
        :param mash: Acceptable deviations from median MASH distances
        :param assembly_summary: a pandas DataFrame with assembly summary information
        """
        self.max_unknowns = max_unknowns
        self.contigs = contigs
        self.assembly_size = assembly_size
        self.mash = mash
        self.assembly_summary = assembly_summary
        self.deviation_values = [max_unknowns, contigs, assembly_size, mash]
        self.ncpus = processes
        self.path = os.path.abspath(path)
        self.name = os.path.basename(os.path.normpath(path))
        self.log = logbook.Logger(self.name)
        self.qc_dir = os.path.join(self.path, "qc")
        self.label = '-'.join(map(str, self.deviation_values))
        self.qc_results_dir = os.path.join(self.qc_dir, self.label)
        self.passed_dir = os.path.join(self.qc_results_dir, "passed")
        self.stats_path = os.path.join(self.qc_dir, 'stats.csv')
        self.nw_path = os.path.join(self.qc_dir, 'tree.nw')
        self.dmx_path = os.path.join(self.qc_dir, 'dmx.csv')
        self.failed_path = os.path.join(self.qc_results_dir, "failed.csv")
        self.tree_img = os.path.join(self.qc_results_dir, "tree.svg")
        self.summary_path = os.path.join(self.qc_results_dir, "summary.txt")
        self.allowed_path = os.path.join(self.qc_results_dir, "allowed.p")
        self.paste_file = os.path.join(self.qc_dir, 'all.msh')
        # Figure out if defining these as None is necessary
        self.tree = None
        self.stats = None
        self.dmx = None
        if os.path.isfile(self.stats_path):
            self.stats = pd.read_csv(self.stats_path, index_col=0)
        if os.path.isfile(self.nw_path):
            self.tree = Tree(self.nw_path, 1)
        if os.path.isfile(self.failed_path):
            self.failed_report = pd.read_csv(self.failed_path, index_col=0)
        if os.path.isfile(self.dmx_path):
            try:
                self.dmx = pd.read_csv(self.dmx_path, index_col=0, sep="\t")
                self.log.info("Distance matrix read succesfully")
            except pd.errors.EmptyDataError:
                self.log.exception()
        self.metadata_path = os.path.join(self.qc_dir, "{}_metadata.csv".format(self.name))
        try:
            self.metadata_df = pd.read_csv(self.metadata_path, index_col="accession")
        except FileNotFoundError:
            self.metadata_df = pd.DataFrame(columns=["accession"])
        self.criteria = ["unknowns", "contigs", "assembly_size", "distance"]
        self.tolerance = {"unknowns": max_unknowns, "contigs": contigs,
                          "assembly_size": assembly_size, "distance": mash}
        self.passed = self.stats
        self.failed = {}
        self.med_abs_devs = {}
        self.dev_refs = {}
        self.allowed = {"unknowns": max_unknowns}
        self.colors = {"unknowns": "red", "contigs": "green",
                       "distance": "purple", "assembly_size": "orange"}
        self.genomes = [Genome.Genome(genome, self.assembly_summary)
                        for genome in self.genome_paths]
        self.assess_tree()

    def __str__(self):
        self.message = [
            "Species: {}".format(self.name),
            "Maximum Unknown Bases:  {}".format(self.max_unknowns),
            "Acceptable Deviations,",
            "Contigs, {}".format(self.contigs),
            "Assembly Size, {}".format(self.assembly_size),
            "MASH: {}".format(self.mash)]
        return '\n'.join(self.message)

    def assess(f):
        # TODO: This can have a more general application if the pickling
        # functionality is implemented elsewhere
        @functools.wraps(f)
        def wrapper(self):
            try:
                assert self.stats is not None
                assert os.path.isfile(self.allowed_path)
                assert (sorted(self.genome_ids().tolist()) ==
                        sorted(self.stats.index.tolist()))
                self.complete = True
                with open(self.allowed_path, 'rb') as p:
                    self.allowed = pickle.load(p)
                self.log.info('Already complete')
            except AssertionError:
                self.complete = False
                f(self)
        return wrapper

    def assess_tree(self):
        try:
            assert self.tree is not None
            assert self.stats is not None
            leaf_names = [re.sub(".fasta", "", i) for i in
                          self.tree.get_leaf_names()]
            assert (sorted(leaf_names) ==
                    sorted(self.stats.index.tolist()) ==
                    sorted(self.genome_ids().tolist()))
            self.tree_complete = True
            self.log.info("Tree already complete")
        except AssertionError:
            self.tree_complete = False

    @property
    def genome_paths(self, ext="fasta"):
        # Why doesn't this work when importing at top of file?
        """Returns a generator for every file ending with `ext`

        :param ext: File extension of genomes in species directory
        :returns: Generator of Genome objects for all genomes in species dir
        :rtype: generator
        """
        return [os.path.join(self.path, genome) for genome in os.listdir(self.path)
                if genome.endswith(ext)]

    # @property
    # def genomes(self):
    #     """Returns a generator for every file ending with `ext`

    #     :param ext: File extension of genomes in species directory
    #     :returns: Generator of Genome objects for all genomes in species dir
    #     :rtype: generator
    #     """
    #     return (Genome.Genome(genome, self.assembly_summary) for genome in self.genome_paths)

    @property
    def total_genomes(self):
        return len(list(self.genomes))

    def sketches(self):
        return (i.msh for i in self.genomes)

    def genome_ids(self):
        ids = [i.name for i in self.genomes]
        return pd.Index(ids)

    # may be redundant. see genome_ids attrib
    @property
    def accession_ids(self):
        ids = [i.accession_id for i in self.genomes
               if i.accession_id is not None]
        return ids

    def mash_paste(self):
        if os.path.isfile(self.paste_file):
            os.remove(self.paste_file)
        sketches = os.path.join(self.qc_dir, "*msh")
        cmd = "mash paste {} {}".format(self.paste_file, sketches)
        Popen(cmd, shell="True", stderr=DEVNULL).wait()
        self.log.info("MASH paste completed")
        if not os.path.isfile(self.paste_file):
            self.log.error("MASH paste failed")
            self.paste_file = None

    def mash_dist(self):
        cmd = "mash dist -p {} -t '{}' '{}' > '{}'".format(
            self.ncpus, self.paste_file, self.paste_file, self.dmx_path)
        Popen(cmd, shell="True", stderr=DEVNULL).wait()
        self.log.info("MASH distance completed")
        self.dmx = pd.read_csv(self.dmx_path, index_col=0, sep="\t")
        # Make distance matrix more readable
        names = [os.path.splitext(i)[0].split('/')[-1] for i in self.dmx.index]
        self.dmx.index = names
        self.dmx.columns = names
        self.dmx.to_csv(self.dmx_path, sep="\t")
        self.log.info("dmx.csv created")

    def sketch_genomes(self):
        """Sketch all genomes"""
        with Pool(ncpus=self.ncpus) as pool:
            self.log.info("{} cpus in pool".format(pool.ncpus))
            pool.map(Genome.sketch_genome, self.genome_paths)
        self.log.info("All genomes sketched")

    def get_tree(self):
        # Use decorator instead of if statement
        if self.tree_complete is False:
            from ete3.coretype.tree import TreeError
            import numpy as np
            # import matplotlib as mpl
            # mpl.use('TkAgg')
            from skbio.tree import TreeNode
            from scipy.cluster.hierarchy import weighted
            ids = ['{}.fasta'.format(i) for i in self.dmx.index.tolist()]
            triu = np.triu(self.dmx.as_matrix())
            hclust = weighted(triu)
            t = TreeNode.from_linkage_matrix(hclust, ids)
            nw = t.__str__().replace("'", "")
            self.tree = Tree(nw)
            # midpoint root tree
            try:
                self.tree.set_outgroup(self.tree.get_midpoint_outgroup())
            except TreeError as e:
                self.log.exception()
            self.tree.write(outfile=self.nw_path)

    def get_stats(self):
        """Get stats for all genomes. Concat the results into a DataFrame"""
        dmx_mean = [self.dmx.mean()] * len(self.genome_paths)
        with Pool(ncpus=self.ncpus) as pool:
            results = pool.map(Genome.mp_stats, self.genome_paths, dmx_mean)
        self.stats = pd.concat(results)
        self.stats.to_csv(self.stats_path)
        self.log.info("Generated stats and wrote to disk")

    def MAD(self, df, col):
        """Get the median absolute deviation for col"""
        MAD = abs(df[col] - df[col].median()).mean()
        return MAD

    def MAD_ref(MAD, tolerance):
        """Get the reference value for median absolute deviation"""
        dev_ref = MAD * tolerance
        return dev_ref

    def bound(df, col, dev_ref):
        lower = df[col].median() - dev_ref
        upper = df[col].median() + dev_ref
        return lower, upper

    def filter_unknown_bases(self):
        """Filter out genomes with too many unknown bases."""
        self.failed["unknowns"] = self.stats.index[
            self.stats["unknowns"] > self.tolerance["unknowns"]]
        self.passed = self.stats.drop(self.failed["unknowns"])
        self.log.info("Analyzed unknowns")

    def check_passed_count(f):
        """
        Count the number of genomes in self.passed.
        Commence with filtering only if self.passed has more than five genomes.
        """
        @functools.wraps(f)
        def wrapper(self, *args):
            if len(self.passed) > 5:
                f(self, *args)
            else:
                self.allowed[args[0]] = ''
                self.failed[args[0]] = ''
                self.log.info("Stopped filtering after {}".format(f.__name__))
        return wrapper

    @check_passed_count
    def filter_contigs(self, criteria):
        """
        Only look at genomes with > 10 contigs to avoid throwing off the
        median absolute deviation.
        Median absolute deviation - Average absolute difference between
        number of contigs and the median for all genomes
        Extract genomes with < 10 contigs to add them back in later.
        Add genomes with < 10 contigs back in
        """
        eligible_contigs = self.passed.contigs[self.passed.contigs > 10]
        not_enough_contigs = self.passed.contigs[self.passed.contigs <= 10]
        # TODO Define separate function for this
        med_abs_dev = abs(eligible_contigs - eligible_contigs.median()).mean()
        self.med_abs_devs["contigs"] = med_abs_dev
        # Define separate function for this
        # The "deviation reference"
        dev_ref = med_abs_dev * self.contigs
        self.dev_refs["contigs"] = dev_ref
        self.allowed["contigs"] = eligible_contigs.median() + dev_ref
        self.failed["contigs"] = eligible_contigs[
            abs(eligible_contigs - eligible_contigs.median()) > dev_ref].index
        eligible_contigs = eligible_contigs[
            abs(eligible_contigs - eligible_contigs.median()) <= dev_ref]
        eligible_contigs = pd.concat([eligible_contigs, not_enough_contigs])
        eligible_contigs = eligible_contigs.index
        self.passed = self.passed.loc[eligible_contigs]
        self.log.info("Analyzed contigs")

    @check_passed_count
    def filter_MAD_range(self, criteria):
        """
        Filter based on median absolute deviation.
        Passing values fall within a lower and upper bound.
        """
        # Get the median absolute deviation
        med_abs_dev = abs(self.passed[criteria] -
                          self.passed[criteria].median()).mean()
        dev_ref = med_abs_dev * self.tolerance[criteria]
        lower = self.passed[criteria].median() - dev_ref
        upper = self.passed[criteria].median() + dev_ref
        allowed_range = (str(int(x)) for x in [lower, upper])
        allowed_range = '-'.join(allowed_range)
        self.allowed[criteria] = allowed_range
        self.failed[criteria] = self.passed[
            abs(self.passed[criteria] -
                self.passed[criteria].median()) > dev_ref].index
        self.passed = self.passed[
            abs(self.passed[criteria] -
                self.passed[criteria].median()) <= dev_ref]
        self.log.info("Filtered based on median absolute deviation range")

    @check_passed_count
    def filter_MAD_upper(self, criteria):
        """
        Filter based on median absolute deviation.
        Passing values fall under the upper bound.
        """
        # Get the median absolute deviation
        med_abs_dev = abs(self.passed[criteria] -
                          self.passed[criteria].median()).mean()
        dev_ref = med_abs_dev * self.tolerance[criteria]
        upper = self.passed[criteria].median() + dev_ref
        self.failed[criteria] = self.passed[
            self.passed[criteria] > upper].index
        self.passed = self.passed[
            self.passed[criteria] <= upper]
        upper = "{:.4f}".format(upper)
        self.allowed[criteria] = upper
        self.log.info("Filtered based on MAD upper bound")

    def base_node_style(self):
        from ete3 import NodeStyle, AttrFace
        nstyle = NodeStyle()
        nstyle["shape"] = "sphere"
        nstyle["size"] = 2
        nstyle["fgcolor"] = "black"
        for n in self.tree.traverse():
            n.set_style(nstyle)
            if re.match('.*fasta', n.name):
                nf = AttrFace('name', fsize=8)
                nf.margin_right = 150
                nf.margin_left = 3
                n.add_face(nf, column=0)
        self.log.info("Applied base node style")

    # Might be better in a layout function
    def style_and_render_tree(self, file_types=["svg"]):
        from ete3 import TreeStyle, TextFace, CircleFace
        ts = TreeStyle()
        title_face = TextFace(self.name.replace('_', ' '), fsize=20)
        title_face.margin_bottom = 10
        ts.title.add_face(title_face, column=0)
        ts.branch_vertical_margin = 10
        ts.show_leaf_name = False
        # Legend
        ts.legend.add_face(TextFace(""), column=1)
        for category in ["Allowed", "Tolerance", "Filtered", "Color"]:
            category = TextFace(category, fsize=8, bold=True)
            category.margin_bottom = 2
            category.margin_right = 40
            ts.legend.add_face(category, column=1)
        for i, criteria in enumerate(self.criteria, 2):
            title = criteria.replace("_", " ").title()
            title = TextFace(title, fsize=8, bold=True)
            title.margin_bottom = 2
            title.margin_right = 40
            cf = CircleFace(4, self.colors[criteria], style="sphere")
            cf.margin_bottom = 5
            filtered_count = len(list(
                filter(None, self.failed_report.criteria == criteria)))
            filtered = TextFace(filtered_count, fsize=8)
            filtered.margin_bottom = 5
            allowed = TextFace(self.allowed[criteria], fsize=8)
            allowed.margin_bottom = 5
            allowed.margin_right = 25
            tolerance = TextFace(self.tolerance[criteria], fsize=8)
            tolerance.margin_bottom = 5
            ts.legend.add_face(title, column=i)
            ts.legend.add_face(allowed, column=i)
            ts.legend.add_face(tolerance, column=i)
            ts.legend.add_face(filtered, column=i)
            ts.legend.add_face(cf, column=i)
        for f in file_types:
            out_tree = os.path.join(self.qc_results_dir, 'tree.{}'.format(f))
            self.tree.render(out_tree, tree_style=ts)
            self.log.info("tree.{} generated".format(f))

    def color_tree(self):
        from ete3 import NodeStyle
        self.base_node_style()
        for genome in self.failed_report.index:
            n = self.tree.get_leaves_by_name(genome+".fasta").pop()
            nstyle = NodeStyle()
            nstyle["fgcolor"] = self.colors[self.failed_report.loc[genome, 'criteria']]
            nstyle["size"] = 9
            n.set_style(nstyle)
        self.style_and_render_tree()

    def filter(self):
        self.filter_unknown_bases()
        self.filter_contigs("contigs")
        self.filter_MAD_range("assembly_size")
        self.filter_MAD_upper("distance")
        with open(self.allowed_path, 'wb') as p:
            pickle.dump(self.allowed, p)
            self.log.info("Pickled results of filtering")
        self.summary()
        self.write_failed_report()

    def write_failed_report(self):
        from itertools import chain
        if os.path.isfile(self.failed_path):
            os.remove(self.failed_path)
        ixs = chain.from_iterable([i for i in self.failed.values()])
        self.failed_report = pd.DataFrame(index=ixs, columns=["criteria"])
        for criteria in self.failed.keys():
            if type(self.failed[criteria]) == pd.Index:
                self.failed_report.loc[self.failed[criteria],
                                       'criteria'] = criteria
        self.failed_report.to_csv(self.failed_path)
        self.log.info("Wrote failed report")

    def summary(self):
        summary = [
            self.name,
            "Unknown Bases",
            "Allowed: {}".format(self.allowed["unknowns"]),
            "Tolerance: {}".format(self.tolerance["unknowns"]),
            "Filtered: {}".format(len(self.failed["unknowns"])),
            "\n",
            "Contigs",
            "Allowed: {}".format(self.allowed["contigs"]),
            "Tolerance: {}".format(self.tolerance["contigs"]),
            "Filtered: {}".format(len(self.failed["contigs"])),
            "\n",
            "Assembly Size",
            "Allowed: {}".format(self.allowed["assembly_size"]),
            "Tolerance: {}".format(self.tolerance["assembly_size"]),
            "Filtered: {}".format(len(self.failed["assembly_size"])),
            "\n",
            "MASH",
            "Allowed: {}".format(self.allowed["distance"]),
            "Tolerance: {}".format(self.tolerance["distance"]),
            "Filtered: {}".format(len(self.failed["distance"])),
            "\n"]
        summary = '\n'.join(summary)
        with open(os.path.join(self.summary_path), "w") as f:
            f.write(summary)
            self.log.info("Wrote QC summary")
        return summary

    def link_genomes(self):
        if not os.path.exists(self.passed_dir):
            os.mkdir(self.passed_dir)
        for genome in self.passed.index:
            fname = "{}.fasta".format(genome)
            src = os.path.join(self.path, fname)
            dst = os.path.join(self.passed_dir, fname)
            try:
                os.link(src, dst)
            except FileExistsError:
                pass
        self.log.info("Links created for genomes that passed QC")

    @assess
    def qc(self):
        if not os.path.isdir(self.qc_dir):
            os.mkdir(self.qc_dir)
        if not os.path.isdir(self.qc_results_dir):
            os.mkdir(self.qc_results_dir)
        self.sketch_genomes()
        self.mash_paste()
        self.mash_dist()
        self.get_stats()
        self.filter()
        self.link_genomes()
        self.get_tree()
        self.color_tree()
        self.log.info("qc command completed")

    def metadata(self):
        metadata = []
        for genome in self.genomes:
            if genome.accession_id in self.metadata_df.index:
                continue
            genome.get_metadata()
            metadata.append(genome.metadata)
        self.metadata_df = pd.concat([self.metadata_df,
                                      pd.DataFrame(metadata).set_index("accession")])
        self.metadata_df.to_csv(self.metadata_path)
        self.log.info("Completed metadata command")
