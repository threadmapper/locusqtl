#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103, W1514, C0301, C0413, R0902


"""

- Find significant markers
- boxplots

"""

# from dataclassy import dataclass

import collections
from dataclasses import dataclass, field
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests


@dataclass
class Marker:
    """ For storing point-wise marker attributes  """

    name: str = field(init=False)
    linkage_group: str = field(init=False)
    cm: float = field(init=False)
    genotype: List[str] = field(init=False)
    T: float = field(init=False)
    p_value: float = field(init=False)
    padj_bonferroni: float = field(init=False)
    accession_pheno: None


    def __str__(self):
        return f'{self.name} {self.linkage_group}{self.cm}{"".join(self.genotype)}{self.T} {self.p_value}'

    def box_plot(self):
        """ plotting adjusted P values against each marker """

        table = []
        for x in self.accession_pheno['J']:
            table.append(['J', x])
        for x in self.accession_pheno['C']:
            table.append(['C', x])

        df = pd.DataFrame(table, columns=['parent', 'phenotype'])

        plt.close('all')
        boxplot = sns.boxplot(x='parent', y='phenotype', data=df, order=["J", "C"])

        # adding data points
        sns.stripplot(x='parent', y='phenotype', data=df, color="grey")

        boxplot.axes.set_title(f"{self.name} Padj: {self.padj_bonferroni:.2e}", fontsize=16)
        boxplot.set_xlabel("Parent", fontsize=14)
        boxplot.set_ylabel("Phenotype", fontsize=14)

        plt.savefig(self.name + '.png')
        plt.savefig(self.name + '.pdf')
        plt.savefig(self.name + '.svg')
        plt.close('all')


@dataclass
class Linkage:
    """ Marker tested in each linkage group against phenotype """

    num_tests: int = 0
    markers: List[Marker] = field(default_factory=list)
    rils: List[float] = field(default_factory=list)
    phenotype: List[float] = field(default_factory=list)

    def __post_init__(self):
        """ load the genotypes by markers """
        if not self.markers:
            self.markers = self.get_genotype(geno_filer='CameorF7-pheno.csv')

    def get_genotype(self, geno_filer='CameorF7-pheno.csv'):
        """"
        parse
        return genotype strings
        """
        marker_list = []
        num = 0
        with open(geno_filer) as inp:
            for line in inp:
                marker, linkage_group, cm, *genotype = line.strip().split(',')
                assert len(genotype) == 61
                if marker == 'Pheno':
                    self.phenotype = [float(_) for _ in genotype]
                    continue
                if marker == 'marker':
                    self.rils = genotype
                    continue

                m = Marker(accession_pheno=None)
                m.name = marker
                m.linkage_group = linkage_group
                m.cm = cm
                m.genotype = genotype

                m.accession_pheno = collections.defaultdict(list)
                for gen, pheno in zip(m.genotype, self.phenotype):
                    m.accession_pheno[gen].append(pheno)

                m.T, m.p_value = stats.ttest_ind(m.accession_pheno['C'],
                                                 m.accession_pheno['J'])  # alternative='greater',

                marker_list.append(m)
                num += 1
                # if num > 4: break

        self.num_tests = len(marker_list)

        return marker_list

    def manhattan_plot(self, p_adjusted_vector):
        """"  plot all linkage groups  """

        df = pd.DataFrame({
            'marker': [m.name for m in self.markers],
            'LG': [m.linkage_group for m in self.markers],
            'cm': [float(m.cm) for m in self.markers],
            't_pvalues': [m.p_value for m in self.markers],
            'T': [m.T for m in self.markers],
            'padj': p_adjusted_vector,
        })

        df['-log10(padj)'] = -np.log10(df['padj'] + 1e-10)  # minimum p value 4.52657012417506E-09

        running_pos = 0
        cumulative_pos = []
        for _, group_df in df.groupby('LG'):
            cumulative_pos.append(group_df['cm'] + running_pos)
            running_pos += group_df['cm'].max()

        df['cumulative_pos'] = pd.concat(cumulative_pos)

        g = sns.relplot(
            data=df,
            x='cumulative_pos',
            y='-log10(padj)',
            aspect=4,
            hue='LG',
            palette='tab10',
            linewidth=0,
            s=6,
            legend='full'
        )

        for ax in g.axes.flat:
            ax.set_xlabel('LG')
            ax.set_xticks(df.groupby('LG')['cumulative_pos'].median())
            ax.set_xticklabels(df['LG'].unique())
            ax.axhline(-np.log10(0.05), linestyle='--', linewidth=0.3)

        plt.savefig('CameorF7-manhattan-BON.svg')
        plt.savefig('CameorF7-manhattan-BON.pdf')
        plt.savefig('CameorF7-manhattan-BON.png')
        plt.close('all')

        df.to_csv('CameorF7-manhattan-BON.csv', encoding='utf-8', index=False)

        print(df.tail())


if __name__ == '__main__':

    qtl = Linkage()

    t_pvalues = [mark.p_value for mark in qtl.markers]

    # control p-values
    p_adjusted_vect = multipletests(t_pvalues, method='bonferroni')[1].tolist()
    # p_adjusted_vect = multipletests(t_pvalues, method='fdr_bh')[1].tolist()

    print(f'number of tests: {len(p_adjusted_vect)}')
    qtl.manhattan_plot(p_adjusted_vect)

    for mark, padj in zip(qtl.markers, p_adjusted_vect):
        mark.padj_bonferroni = padj
        if mark.name in ['AX-183567168', 'AX-183567200', 'AX-183861105', 'AX-183578342']:
            mark.box_plot()

    print('DONE')
