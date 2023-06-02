import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn


def make_indiv_heatmap(significant_var_df, row_order, name):
    sample_order = ['OC1', 'NC1', 'NC2', 'NC3',
                    'NC4', 'NC5', 'NC6', 'NC7', 'NC8', ]
    heatmap_data = pd.pivot_table(significant_var_df, values='Score',  # Score
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

    # export data to excel
    path = './analysis/'

    # plot heatmap
    seaborn.set(context="paper", palette="Paired", style="darkgrid",
                font_scale=2, rc={'ytick.left': True, 'xtick.bottom': True})
    fig, ax = plt.subplots(figsize=(36, 24))

    ax.xaxis.tick_top()
    cmap = "Reds"
    heat = seaborn.heatmap(heatmap_data, vmin=15, vmax=45, cmap=cmap, ax=ax, cbar_kws={
                           'label': 'CADD Deleteriousness Score'}, annot=annots, fmt=".2f")  # square=True, )
    plt.yticks(rotation=0)
    heat.set_xlabel('')

    heat.get_figure().savefig(path + name + '-heatmap.pdf', dpi=200)
    plt.clf()


def calc_p_value(variant_samples, sample_groups):
    num_NC_var = len(sample_groups['NC-OC'].intersection(variant_samples))
    perc_NC_var = num_NC_var / len(sample_groups['NC-OC'])

    p_val = .01

    return p_val, perc_NC_var


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
            if 'RC' in sample:
                continue

            if (total_reads[i] >= 4) and ((alt_reads[i]/total_reads[i]) >= 0.10):
                filtered_samples.append(sample)
            elif total_reads[i] < 4:
                filtered_samples.append(sample)

        if len(filtered_samples) > 0:
            num_variants_evaluated += 1

        p_val, perc_nc = calc_p_value(filtered_samples, sample_groups)
        if p_val <= .05 or do_not_filter:
            for i, sample in enumerate(filtered_samples):
                if 'inframe' in row['Sequence_Ontology']:
                    continue
                if type(row['cDNA_change']) == str:
                    var_id = row['cDNA_change']
                elif type(row['Sequence_Ontology']) == str:
                    var_id = row['Sequence_Ontology']
                else:
                    var_id = row['Chrom'] + "-" + str(row["Position"])

                entry = pd.DataFrame({'Sample': [sample], 'Variant': [str(row['Gene']) + " - " + var_id], 'Present': [1], 'Score': [
                                     score], 'p-val': [p_val], 'perc_samp': [perc_nc], 'perc_reads': [(alt_reads[i]/total_reads[i])]})
                df = pd.concat([df, entry])

    if len(df) == 0:
        entry = pd.DataFrame({'Sample': ['nan'], 'Variant': [
                             'No Variants'], 'Present': [0], 'Score': [0], 'p-val': [1]})
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

    epigenetic_genes_coding_variants = pd.read_csv(
        './annotated_variant_files/epigenetic_modifiers_coding.tsv', sep="\t", header=[6])

    # Select only variants that are significantly different between NC/OC and RC
    name = 'Epigenetic Modifiers NC Only'
    significant_variant_df = create_significant_variant_df(
        epigenetic_genes_coding_variants, sample_groups, name)
    row_order = significant_variant_df.sort_values(
        by=['perc_samp', 'Variant'], ascending=[False, True])['Variant'].unique()
    print(row_order)

    make_indiv_heatmap(significant_variant_df, row_order, name)


main()
