import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn
from scipy.stats import fisher_exact


def make_heatmap(significant_var_dfs, row_orders, name):
    concated_heatmap = pd.DataFrame()
    annots_concat = pd.DataFrame()
    for i, variant_df in enumerate(significant_var_dfs):
        heatmap_data = pd.pivot_table(variant_df, values='Score',
                                      index=['Variant'],
                                      columns='Sample')
        annots = pd.pivot_table(variant_df, values='perc_reads',
                                index=['Variant'],
                                columns='Sample')
        heatmap_data = heatmap_data.reindex(row_orders[i])
        annots = annots.reindex(row_orders[i])

        concated_heatmap = pd.concat([concated_heatmap, heatmap_data])
        annots_concat = pd.concat([annots_concat, annots])

    sample_order = ['OC1', 'NC1', 'NC2', 'NC3', 'NC4', 'NC5', 'NC6', 'NC7', 'NC8', 'RC1', 'RC2', 'RC3',
                    'RC4', 'RC5', 'RC6', 'RC7', 'RC8', 'RC9', 'RC10', 'RC11', 'RC12', 'RC13', 'RC14', 'RC15', 'RC16']
    for sample in sample_order:
        if sample not in concated_heatmap:
            concated_heatmap[sample] = [np.nan]*len(concated_heatmap)
            annots_concat[sample] = [np.nan]*len(annots_concat)
    concated_heatmap = concated_heatmap[sample_order]
    annots_concat = annots_concat[sample_order]

    # export data to excel
    if 'NC' in name:
        path = './analysis/NC-OC/'
    elif 'OC' in name:
        path = './analysis/RC/'
    else:
        path = './analysis/'

    concated_heatmap.to_excel(path+name+".xlsx")

    # plot heatmap
    seaborn.set(context="paper", palette="Paired", style="darkgrid",
                font_scale=2, rc={'ytick.left': True, 'xtick.bottom': True})
    fig, ax = plt.subplots(figsize=(36, 16))

    ax.xaxis.tick_top()
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    cmap = "Reds"
    heat = seaborn.heatmap(concated_heatmap, vmin=15, vmax=45, cmap=cmap, ax=ax, cbar_kws={'pad': .15,
                                                                                           'label': 'CADD Deleteriousness Score'})  # , annot=annots_concat, fmt=".2f")  # annot=True, fmt=".0f", )
    plt.yticks(rotation=0)
    heat.set_xlabel('')
    plt.subplots_adjust(bottom=0.28)
    heat.get_figure().savefig(path + name + '-heatmap.svg', dpi=200)
    plt.clf()


def make_indiv_heatmap(significant_var_df, row_order, name):
    sample_order = ['OC1', 'NC1', 'NC2', 'NC3', 'NC4', 'NC5', 'NC6', 'NC7', 'NC8', 'RC1', 'RC2', 'RC3',
                    'RC4', 'RC5', 'RC6', 'RC7', 'RC8', 'RC9', 'RC10', 'RC11', 'RC12', 'RC13', 'RC14', 'RC15', 'RC16']
    heatmap_data = pd.pivot_table(significant_var_df, values='Score',
                                  index=['Variant'],
                                  columns='Sample')
    annots = pd.pivot_table(significant_var_df, values='perc_reads',
                            index=['Variant'],
                            columns='Sample')

    heatmap_data = heatmap_data.reindex(row_order)
    annots = annots.reindex(row_order)

    for sample in sample_order:
        if sample not in heatmap_data:
            heatmap_data[sample] = [np.nan]*len(heatmap_data)
            annots[sample] = [np.nan]*len(annots)
    heatmap_data = heatmap_data[sample_order]
    annots = annots[sample_order]
    if len(heatmap_data) > 1:
        heatmap_data.drop(heatmap_data.tail(1).index, inplace=True)
        annots.drop(annots.tail(1).index, inplace=True)

    # export data to excel
    path = './analysis/'

    # plot heatmap
    seaborn.set(context="paper", palette="Paired", style="darkgrid",
                font_scale=2, rc={'ytick.left': True, 'xtick.bottom': True})
    fig, ax = plt.subplots(figsize=(36, 24))

    ax.xaxis.tick_top()
    cmap = "Reds"
    heat = seaborn.heatmap(heatmap_data, square=True, vmin=15, vmax=45, cmap=cmap, ax=ax, cbar_kws={
                           'label': 'CADD Deleteriousness Score'}, annot=annots, fmt=".2f")  # square=True, annot=True, fmt=".0f", )

    plt.yticks(rotation=0)
    heat.set_xlabel('')
    heat.get_figure().savefig(path + name + '-heatmap.pdf', dpi=200)
    plt.clf()


def calc_p_value(variant_samples, sample_groups):
    num_NC_var = len(sample_groups['NC-OC'].intersection(variant_samples))
    num_NC_ref = len(sample_groups['NC-OC']) - num_NC_var
    perc_NC_var = num_NC_var / len(sample_groups['NC-OC'])

    num_RC_var = len(sample_groups['RC'].intersection(variant_samples))
    num_RC_ref = len(sample_groups['RC']) - num_RC_var

    if perc_NC_var > 0:
        p_val = 1
    else:
        odds, p_val = fisher_exact(
            [[num_NC_ref, num_NC_var], [num_RC_ref, num_RC_var]])

    return p_val


def create_significant_variant_df(variants, sample_groups, name, do_not_filter=False):
    df = pd.DataFrame()
    num_variants_evaluated = 0
    for idx, row in variants.iterrows():
        # creates list of "clean" sample names
        variant_samples = [x.split('.')[0] for x in row['Samples'].split(';')]
        alt_reads = [float(x) for x in row['Alternate_reads'].split(';')]
        total_reads = [float(x) for x in row['Total_reads'].split(';')]

        # necessary bc missing indel deleteriousness
        score = 15 if math.isnan(row['Phred']) else row['Phred']
        if score < 15:
            continue
        # Filter variants
        filtered_samples = []
        for i, sample in enumerate(variant_samples):
            if (total_reads[i] >= 4) and ((alt_reads[i]/total_reads[i]) >= 0.10):
                filtered_samples.append(sample)
            elif total_reads[i] < 4:
                filtered_samples.append(sample)
        if len(filtered_samples) > 0:
            num_variants_evaluated += 1

        p_val = calc_p_value(filtered_samples, sample_groups)
        if p_val <= .05 or do_not_filter:
            for i, sample in enumerate(filtered_samples):
                if 'inframe' in row['Sequence_Ontology']:
                    continue
                if type(row['Protein_Change']) == str:
                    if 'ins' in row['Protein_Change']:
                        var_id = row['Sequence_Ontology']
                    else:
                        var_id = row['Protein_Change']
                elif type(row['cDNA_change']) == str:
                    var_id = row['cDNA_change']
                elif type(row['Sequence_Ontology']) == str:
                    var_id = row['Sequence_Ontology']
                else:
                    var_id = row['Chrom'] + "-" + str(row["Position"])

                entry = pd.DataFrame({'Sample': [sample], 'Variant': [str(row['Gene']) + " - " + var_id], 'Present': [
                                     1], 'Score': [score], 'p-val': [p_val], 'perc_reads': [(alt_reads[i]/total_reads[i])]})
                df = pd.concat([df, entry])

    if len(df) == 0:
        entry = pd.DataFrame({'Sample': ['nan'], 'Variant': ['No Variants'], 'Present': [
                             0], 'Score': [0], 'p-val': [1], 'perc_reads': 0})
        df = pd.concat([df, entry])
        
    # Create spacer row
    entry = pd.DataFrame({'Sample': ['nan'], 'Variant': [''], 'Present': [
                         0], 'Score': [0], 'p-val': [1], 'perc_reads': 0})
    df = pd.concat([df, entry])

    print(name, "Num Variants:", num_variants_evaluated)

    return df


def main():
    samples = ['NC1', 'NC4', 'NC7', 'RC1', 'RC12', 'RC15', 'RC3', 'RC6', 'RC9',
               'NC2', 'NC5', 'NC8', 'RC10', 'RC13', 'RC16', 'RC4', 'RC7',
               'NC3', 'NC6', 'OC1', 'RC11', 'RC14', 'RC2', 'RC5', 'RC8']

    # Separate out NC/OC and RC analysis
    sample_groups = {'NC-OC': set(filter(lambda x: 'NC' in x or 'OC' in x, samples)),
                     'RC': set(filter(lambda x: 'RC' in x, samples))}

    # Load all gene sets
    cadd_variants = pd.read_csv(
        './annotated_variant_files/cadd_phred_over_15.tsv', sep="\t", header=[6])

    variable_genes_variants = pd.read_csv(
        './annotated_variant_files/top500variable_coding.tsv', sep="\t", header=[6])

    deg_variants = pd.read_csv(
        './annotated_variant_files/degFC2p0.01_coding.tsv', sep="\t", header=[6])

    clone_gene_variants = pd.read_csv(
        './annotated_variant_files/clonegenesFC1p0.05_coding.tsv', sep="\t", header=[6])

    genes_2017_coding_variants = pd.read_csv(
        './annotated_variant_files/2017_genes_coding.tsv', sep="\t", header=[6])

    resistance_genes_coding_variants = pd.read_csv(
        './annotated_variant_files/resistance_genes_coding.tsv', sep="\t", header=[6])

    epigenetic_genes_coding_variants = pd.read_csv(
        './annotated_variant_files/epigenetic_modifiers_coding.tsv', sep="\t", header=[6])

    # Select only variants that are significantly different between NC/OC and RC
    names = ['Genes from 2017 Paper variants', 'DEG FC2 p0.01 variants', 'Clone Genes FC1 p0.05 variants',
             'Top 500 Variable Genes variants', 'Epigenetic Modifiers', 'CADD >15 variants', 'Resistance+Pos Control Genes variants']
    significant_variant_dfs = []
    row_orders = []
    for i, variant_data in enumerate([genes_2017_coding_variants, deg_variants, clone_gene_variants, variable_genes_variants, epigenetic_genes_coding_variants, cadd_variants, resistance_genes_coding_variants]):
        if i == 6:  # Resistance genes
            significant_variant_dfs.append(create_significant_variant_df(
                variant_data, sample_groups, names[i], do_not_filter=True))
            row_orders.append(significant_variant_dfs[i]['Variant'].unique())
        else:
            significant_variant_dfs.append(create_significant_variant_df(
                variant_data, sample_groups, names[i]))
            row_orders.append(significant_variant_dfs[i].sort_values(
                by=['p-val'], ascending=True)['Variant'].unique())

    # Add in Vemurafenib resistance genes w/o p-val filtering

    # Make heatmaps with gene x variant x sample data
    make_heatmap(significant_variant_dfs[:-1],
                 row_orders[:-1], "genetic-resistance")

    make_indiv_heatmap(significant_variant_dfs[-1], row_orders[-1], names[-1])


main()
