import matplotlib.pyplot as plt
import pandas as pd
import seaborn
from scipy.stats import hypergeom
from matplotlib.colors import LogNorm


def make_heatmap(concated_data, samples, name):
    heatmap_data = pd.pivot_table(concated_data, values='p_val',
                                  index=['Sample 2'],
                                  columns='Sample 1')
    heatmap_data = heatmap_data[samples]

    heatmap_data = heatmap_data.reindex(samples)

    # export data to excel
    path = './analysis/'

    heatmap_data.to_excel(path+name+".xlsx")

    # plot heatmap
    seaborn.set(context="paper", palette="Paired", style="darkgrid",
                font_scale=2, rc={'ytick.left': True, 'xtick.bottom': True})
    fig, ax = plt.subplots(figsize=(36, 16))

    # ax.xaxis.tick_top()
    heat = seaborn.heatmap(heatmap_data, ax=ax, cbar_kws={
                           'label': 'Hypergeometric Test P-Value'}, annot=False, fmt=".0f", norm=LogNorm())
    plt.yticks(rotation=0)
    plt.subplots_adjust(bottom=0.28)
    heat.get_figure().savefig(path + name + '-heatmap.svg', dpi=200)
    plt.clf()

    return


def process_variants(variants, samples):
    variant_sets = {}
    all_variant_dict = {}

    for sample in samples:
        variant_sets[sample] = set()

    for index, row in variants.iterrows():
        # creates list of "clean" sample names
        variant_samples = [x.split('.')[0] for x in row['Samples'].split(';')]
        if len(variant_samples) == 25:
            continue
        var_id = "_".join([row['Chrom'], str(row['Position']),
                           row['Ref_Base'], row['Alt_Base']])

        for sample in variant_samples:
            variant_sets[sample].add(var_id)
            if var_id in all_variant_dict:
                all_variant_dict[var_id] += 1
            else:
                all_variant_dict[var_id] = 1

    return variant_sets, all_variant_dict


def generate_sample_overlap_stats(samples, variant_sets, all_variants):
    all_vars_list = list(all_variants.keys())

    percent_overlap_dict = {'Sample 1': [],
                            'Sample 2': [], 'Percent Overlap': [], 'p_val': []}
    for i in range(len(samples)):
        for j in range(i, len(samples)):
            percent_overlap_dict['Sample 1'].append(samples[i])
            percent_overlap_dict['Sample 2'].append(samples[j])
            overlap = 100 * len(variant_sets[samples[i]].intersection(variant_sets[samples[j]])) / len(
                variant_sets[samples[i]].union(variant_sets[samples[j]]))
            percent_overlap_dict['Percent Overlap'].append(overlap)

            p_val = hypergeom.sf(len(variant_sets[samples[i]].intersection(variant_sets[samples[j]]))-1, len(
                all_vars_list), len(variant_sets[samples[i]]), len(variant_sets[samples[j]]))

            percent_overlap_dict['p_val'].append(p_val)

    return pd.DataFrame(percent_overlap_dict)


def main():
    samples = ['OC1', 'NC1', 'NC2', 'NC3', 'NC4', 'NC5', 'NC6', 'NC7', 'NC8', 'RC1', 'RC2',
               'RC3', 'RC4', 'RC5', 'RC6', 'RC7', 'RC8', 'RC9', 'RC10', 'RC11', 'RC12', 'RC13',
               'RC14', 'RC15', 'RC16']

    # Load all gene sets
    variant_annot_path = './annotated_variant_files/cadd_phred_over_15.tsv'
    variants = pd.read_csv(variant_annot_path, sep="\t", header=[6])

    # Process variant information, sort into sets
    variant_sets, all_variant_dict = process_variants(variants, samples)

    # Find overlapping genes
    overlap_df = generate_sample_overlap_stats(
        samples, variant_sets, all_variant_dict)

    # Plot heatmap
    make_heatmap(overlap_df, samples, 'Overlap')


main()
