import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns
from matplotlib.lines import Line2D


def well_number_to_label(well_number):
    rows = 'ABCDEFGH'
    col = (well_number - 1) % 12 + 1
    row = rows[(well_number - 1) // 12]
    return f"{row}{col}"


def layout2png(layout_file, filename, n):
    global s_name
    df = layout_file
    df.columns = df.columns.astype(int)
    base_name = "plate_layout"
    output_file = f"{filename.rsplit('.', 1)[0]}.png"
    plate_rows = list("ABCDEFGH")
    plate_cols = list(range(1, 13))

    well_contents = df.values.flatten()
    gene_targets = sorted(set(
        item.split('_')[1] for item in well_contents
        if item not in ["OMIT", "NC"] and "_" in item
    ))

    colors = sns.color_palette("hsv", len(gene_targets))
    color_dict = dict(zip(gene_targets, colors))

    fig, axes = plt.subplots(8, 12, figsize=(16, 10))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    for row_idx, row in enumerate(plate_rows):
        for col_idx, col in enumerate(plate_cols):
            sample = df.at[row, col]
            well_label = f"{row}{col}"

            if sample == "OMIT":
                color = "white"
            elif sample == "NC":
                color = "gray"
            else:
                gene_target = sample.split('_')[1]
                s_name = sample.split('_')[0]
                color = color_dict.get(gene_target, "white")

            ax = axes[row_idx, col_idx]
            ax.set_facecolor(color)
            ax.text(0.5, 0.5, s_name, ha='center', va='center', fontsize=10, color='black', fontweight='bold')
            ax.set_title(well_label, fontsize=10, fontweight='bold')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect('equal')

    legend_patches = [plt.Rectangle((0, 0), 1, 1, facecolor=color, edgecolor='black') for color in colors]
    ax.legend(legend_patches, gene_targets, title="Gene Targets", loc='upper right', bbox_to_anchor=(1.2, 1))

    if n == 999:
        huts = f"{output_file}"
    else:
        huts = f"{base_name}_{n + 1}"
    plt.suptitle(huts, fontsize=14, fontweight='bold')

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def generate_sequential_layout_1(samples, gene_targets, num_replicates):
    plate_rows = list("ABCDEFGH")
    plate_cols = list(range(1, 13))

    wells_per_sample = len(gene_targets) * num_replicates
    max_wells = 96 - wells_per_sample  # Reserve 4 wells for NC
    wells_needed = len(samples) * wells_per_sample

    num_plates = (wells_needed // max_wells) + (1 if wells_needed % max_wells else 0)
    plate_designs = []

    sample_index = 0
    for plate in range(num_plates):
        current_plate = pd.DataFrame("OMIT", index=plate_rows, columns=plate_cols)
        wells_left = max_wells
        col_main = 1
        row = 0
        for i in range(0, 180, 1):
            if sample_index >= len(samples):
                col = col_main
                sample = "NC"
                for target in gene_targets:
                    for rep in range(1, num_replicates + 1):
                        current_plate.at[plate_rows[row], col] = f"{sample}_{target}_{rep}"
                        col += 1
                break
            wells_left -= wells_per_sample

            if wells_left >= 0:
                sample = samples[sample_index]
                sample_index += 1
            else:
                sample = "NC"
            col = col_main
            for target in gene_targets:
                for rep in range(1, num_replicates + 1):
                    current_plate.at[plate_rows[row], col] = f"{sample}_{target}_{rep}"
                    col += 1
            row += 1
            if row > 7:
                row = 0
                col_main += num_replicates * len(gene_targets)
                if col_main > 12:
                    break

        plate_designs.append(current_plate)

    return plate_designs


def plater(input_file, design, gene_targets, num_replicates):
    if input_file.endswith('.tsv'):
        layout = pd.read_csv(input_file, sep='\t', index_col=0, header=0)
        layout2png(layout, input_file, 999)

    elif input_file.endswith(".txt"):
        if not gene_targets:
            raise ValueError("Please specify gene targets")
        with open(input_file, 'r') as f:
            samples = [line.strip() for line in f.readlines()]
        if design == "sequential_1":
            plate_designs = generate_sequential_layout_1(samples, gene_targets, num_replicates)
        else:
            raise ValueError("Unkown design, use --help for options")
        for i, plate in enumerate(plate_designs):
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            filename = f"{base_name}_layout_{i + 1}.tsv"
            plate.to_csv(filename, sep='\t', index=True)
            print(f"Saved: {filename}")
            layout2png(plate, filename, i)

    else:
        raise ValueError("Provide either a .txt sample list or .tsv precomputed layout as input")


def process_ct(ct_file, layout_file):
    if not layout_file:
        raise ValueError("Provide a layout file through -i_lay for Ct analysis.")
    with open(layout_file, 'r') as file:
        lines = file.readlines()
    df_ct = pd.read_csv(ct_file, skiprows=2, header=0)

    plate_dict = {}
    well_number = 1
    well_id_dict = {}
    rows = list("ABCDEFGH")
    columns = list(range(1, 13))

    for i, line in enumerate(lines[1:]):
        row = line.strip().split('\t')[1:]  # Skip row label
        for j, sample_name in enumerate(row):
            well_id = f"{rows[i]}{columns[j]}"
            plate_dict[well_number] = sample_name
            well_id_dict[well_number] = well_id
            well_number += 1
    df_ct['Well'] = df_ct['Well'].astype(int)
    df_ct['SampleName'] = df_ct['Well'].map(plate_dict).fillna("Unknown")
    df_ct['WellID'] = df_ct['Well'].map(well_id_dict).fillna("Unknown")

    return df_ct


def analyze_ct(ct_dataframe, method, gene_housekeeping, file_name):
    ct_dataframe['Sample'] = ct_dataframe['SampleName'].apply(
        lambda x: x.split('_')[0] if '_' in x else x)
    ct_dataframe['Type'] = ct_dataframe['SampleName'].apply(
        lambda x: x.split('_')[1] if len(x.split('_')) > 1 else 'Unknown')

    result = ct_dataframe.groupby(['Sample', 'Type'])['Ct'].apply(list).reset_index()

    result_pivot = result.pivot(index='Sample', columns='Type', values='Ct').reset_index()

    def average_ct(values):
        tmp_list = []
        val_list = values.iloc[0]
        if isinstance(val_list, list):
            for value in val_list:
                if isinstance(value, (int, float)) and not pd.isna(value):
                    tmp_list.append(value)
            if len(tmp_list) == 0:
                return float('nan')
            return np.mean(tmp_list)
        else:
            return 'Invalid input type'

    primers = [col for col in result_pivot.columns if col != 'Sample']
    for primer in primers:
        result_pivot[f'{primer}_Avg'] = result_pivot[[primer]].apply(average_ct, axis=1)

    if not method:
        raise ValueError("Please specify method for CT analysis through --m")

    if method == "dct":
        if not gene_housekeeping:
            raise ValueError("Please specify housekeeping gene through --hg")

        if gene_housekeeping not in primers:
            raise ValueError(f"Housekeeping gene {gene_housekeeping} not found in primers: {primers}")

        for primer in primers:
            if primer != gene_housekeeping:
                result_pivot[f'dCt_{primer}'] = result_pivot.apply(
                    lambda x: x[f'{primer}_Avg'] - x[f'{gene_housekeeping}_Avg']
                    if isinstance(x[f'{primer}_Avg'], (int, float)) and isinstance(x[f'{gene_housekeeping}_Avg'],
                                                                                   (int, float))
                    else None,
                    axis=1
                )

                result_pivot[f'2^-dCt_{primer}'] = result_pivot.apply(
                    lambda x: 2 ** (-x[f'dCt_{primer}']) if x[f'dCt_{primer}'] is not None else None,
                    axis=1
                )

        df_sorted = result_pivot.sort_values(by='Sample')
        df_filtered = df_sorted[~df_sorted['Sample'].isin(['NC', 'PC'])]
        fig, ax = plt.subplots(figsize=(10, 6))
        df_filtered = df_filtered.set_index("Sample")
        df_filtered = df_filtered.fillna(0)
        df_filtered[[f'2^-dCt_{primer}' for primer in primers if primer != gene_housekeeping]].plot(kind='bar', ax=ax)
        ax.set_xticks(range(len(df_filtered)))
        ax.set_xticklabels(df_filtered.index, rotation=45, ha="right")
        plt.xlabel('Sample')
        plt.ylabel(f'2^−(Ct[Target]−Ct[{gene_housekeeping}])')
        plt.tight_layout()
        plt.show()
        df_sorted.to_csv(f"{file_name}_ct_analyzed.tsv", sep='\t', index=False)
        print(f"Saved: {file_name}_ct_analyzed.tsv")


def process_spec(spec_dataframe, cycles, df_ct, embed_ct):
    df_spec = pd.read_csv(spec_dataframe, header=0)

    filter_colors = {
        'Filter A': 'blue',
        'Filter B': 'red',
        'Filter C': 'green',
        'Filter D': 'purple'
    }

    filtered_df = df_spec[df_spec['Reading'] <= cycles]
    y_min = filtered_df[['Filter A', 'Filter B', 'Filter C', 'Filter D']].min().min()
    y_max = filtered_df[['Filter A', 'Filter B', 'Filter C', 'Filter D']].max().max()

    if df_ct is None and embed_ct:
        raise ValueError("No CT file provided while embed_ct flag is set to True (Default).")

    fig, axes = plt.subplots(8, 12, figsize=(16, 10))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    if df_ct is None or embed_ct == False:
        for well in range(1, 97):
            ax = axes[(well - 1) // 12, (well - 1) % 12]
            well_label = well_number_to_label(well)
            well_data = filtered_df[filtered_df['Well'] == well]

            for filter_name, color in filter_colors.items():
                ax.plot(well_data['Reading'], well_data[filter_name], color=color, linewidth=1.5)

            ax.set_title(well_label, fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim([1, cycles])
            ax.set_ylim([y_min, y_max])

    if df_ct is not None:
        if embed_ct:
            df_ct['Target'] = df_ct['SampleName'].str.split('_').str[1]
            unique_targets = df_ct['Target'].dropna().unique()
            cmap = sns.color_palette("husl", len(unique_targets))
            target_colors = {target: cmap[i] for i, target in enumerate(unique_targets)}
            legend_labels = list(target_colors.keys())
            legend_colors = list(target_colors.values())

            for well in range(1, 97):
                ax = axes[(well - 1) // 12, (well - 1) % 12]
                well_label = well_number_to_label(well)
                well_data = filtered_df[filtered_df['Well'] == well]
                ct_data = df_ct[df_ct['WellID'] == well_label]

                if not ct_data.empty:
                    new_label = ct_data['SampleName'].str.split('_').str[0].iloc[0]
                else:
                    new_label = ""

                bg_color = "black"
                if not ct_data.empty:
                    target = ct_data.iloc[0]['Target']
                    if target in target_colors:
                        bg_color = target_colors[target]

                for filter_name, color in filter_colors.items():
                    ax.plot(well_data['Reading'], well_data[filter_name], color=color, linewidth=1.5)

                ax.set_facecolor(bg_color)

                if not ct_data.empty:
                    ct_value = ct_data.iloc[0]['Ct']
                    if not pd.isna(ct_value):
                        ax.text(0.05, 0.85, f"{ct_value:.2f}", transform=ax.transAxes, fontsize=8,
                                color="white", bbox=dict(facecolor='black', alpha=0.5, edgecolor='none'))

                ax.set_title(new_label, fontsize=8, color='black')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xlim([1, cycles])
                ax.set_ylim([y_min, y_max])

            custom_legend = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10)
                             for color in legend_colors]
            plt.legend(custom_legend, legend_labels, title="Targets", loc='center left', bbox_to_anchor=(1, 0.5),
                       fontsize=8)

    filename = f"{os.path.splitext(spec_dataframe)[0]}.png"
    plt.savefig(filename, dpi=1200)
    plt.close()
    print(f"Saved: {filename}")


def process_dis(dis_dataframe):
    def extract_data(sect_df, index):
        df = pd.DataFrame(columns=["Well", "Data"])
        for zline in sect_df[index][1:]:
            list_item = zline.split(',')
            well_num = list_item[0]
            relevant_data = list_item[3:-1]
            relevant_data = [float(x) for x in relevant_data]
            new_row = pd.DataFrame({"Well": [well_num], "Data": [relevant_data]})
            df = pd.concat([df, new_row], ignore_index=True)
        return df

    with open(dis_dataframe, 'r') as file:
        lines = file.readlines()

    sections = []
    current_section = []
    headers_detected = False

    for line in lines:
        if line.startswith("Well") and headers_detected:
            sections.append(current_section)
            current_section = []
        current_section.append(line)
        headers_detected = True

    sections.append(current_section)  # Append the last section

    if len(sections) < 3:
        raise ValueError(
            "CSV file does not contain the expected sections for Temperature, Raw Data, and Derivative Data.")

    temperature_df = extract_data(sections, 0)

    derivative_df = extract_data(sections, 2)

    # Find global min and max for y-limits
    y_min = float('inf')
    y_max = float('-inf')

    # First, determine the global y limits
    for _, row in temperature_df.iterrows():
        well = int(row["Well"])
        try:
            well_deriv = derivative_df[derivative_df["Well"] == str(well)]["Data"].iloc[0]
            y_min = min(y_min, min(well_deriv))
            y_max = max(y_max, max(well_deriv))
        except ValueError:
            print(f"Skipping invalid well number: {well}")

    # Now plot each well with the same y limits
    fig, axes = plt.subplots(8, 12, figsize=(16, 10))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    for i in range(1, 97):
        well = i
        try:
            row_idx = (well - 1) // 12
            col_idx = (well - 1) % 12

            if 0 <= row_idx < 8 and 0 <= col_idx < 12:
                ax = axes[row_idx, col_idx]
                well_label = well_number_to_label(well)

                well_temp_data = temperature_df[temperature_df["Well"] == str(well)]["Data"]
                if not well_temp_data.empty:
                    well_temp = well_temp_data.iloc[0]
                else:
                    well_temp = None
                well_deriv_data = derivative_df[derivative_df["Well"] == str(well)]["Data"]
                if not well_deriv_data.empty:
                    well_deriv = well_deriv_data.iloc[0]
                else:
                    well_deriv = None

                if well_temp and well_deriv:
                    ax.plot(well_temp, well_deriv, color='red', linewidth=1.5, label='Derivative')

                ax.set_title(well_label, fontsize=8)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xticklabels([])
                ax.set_yticklabels([])

                # Set the same y limits for each subplot
                ax.set_ylim([y_min, y_max])

        except ValueError:
            print(f"Skipping invalid well number: {well}")

    filename = f"{os.path.splitext(dis_dataframe)[0]}_derivative_plot.png"
    plt.savefig(filename, dpi=1200)
    plt.show()
    print(f"Saved: {filename}")


def analyzer(input_file_ct, input_file_spec, input_file_dis, layout_file, method, housekeeping_gene, show_cycles,
             embed_ct):
    count = 0
    if not input_file_ct:
        print("No CT scores provided.")
        count += 1
        df_ct = None
    else:
        print("Processing Ct values...")
        df_ct = process_ct(input_file_ct, layout_file)
        output_base_file = os.path.splitext(input_file_ct)[0]
        analyze_ct(df_ct, method, housekeeping_gene, output_base_file)
    if not input_file_spec:
        print("No raw spectra provided.")
        count += 1
    else:
        print("Processing spectra...")
        process_spec(input_file_spec, show_cycles, df_ct, embed_ct)
    if not input_file_dis:
        print("No dissociation curve spectra provided.")
        count += 1
    else:
        print("Processing dissociation curves...")
        process_dis(input_file_dis)
    if count >= 3:
        raise ValueError("No input file is provided, no analysis will be performed.")


def main():
    parser = argparse.ArgumentParser(description="qPCR Plate Designer and Analyzer")
    subparsers = parser.add_subparsers(dest='command')

    plater_parser = subparsers.add_parser("plater", help="Design a qPCR plate.")
    plater_parser.add_argument("-i", "--input", required=True,
                               help="Sample list (.txt) or prepared layout file (.tsv).")
    plater_parser.add_argument("-d", "--design", required=False, default='sequential_1', type=str,
                               choices=["sequential_1"],
                               help="Plate design for automatic well assignment.")
    plater_parser.add_argument("-t", "--target_genes", required=False, nargs='+',
                               help="List containing 3-letter primer gene target IDs.")
    plater_parser.add_argument("-n", "--number_technical_replicates", required=False, default=2, type=int,
                               help="Number of technical replicates")

    analyzer_parser = subparsers.add_parser("analyzer", help="Analyze qPCR results")
    analyzer_parser.add_argument("-i_ct", "--input_ct", required=False, help="CSV file with ct results.")
    analyzer_parser.add_argument("-i_spec", "--input_spec", required=False, help="CSV file with spectra results.")
    analyzer_parser.add_argument("-i_dis", "--input_dis", required=False,
                                 help="CSV file with dissociation curve results.")
    analyzer_parser.add_argument("-i_lay", "--input_layout", required=False, help="TSV file with plate layout.")

    analyzer_parser.add_argument("-m", "--method", required=False, choices=["dct"], default="dct",
                                 help="Method for analysis")
    analyzer_parser.add_argument("-hg", "--housekeeping_gene", required=False, type=str,
                                 help="Housekeeping target ID used in layout.")
    analyzer_parser.add_argument("-c", "--show_cycles", required=False, type=int, default=40,
                                 help="Number of cycles to show in raw spectra traces.")
    analyzer_parser.add_argument("-e", "--embed_ct", required=False, type=lambda x: x.lower() == 'true', default=True,
                                 help="Whether or not to embed the ct data and layout information into spectra output. "
                                      "Set to 'True' or 'False'.")

    args = parser.parse_args()

    if args.command == "plater":
        plater(args.input, args.design, args.target_genes, args.number_technical_replicates)
    elif args.command == "analyzer":
        analyzer(args.input_ct, args.input_spec, args.input_dis, args.input_layout, args.method, args.housekeeping_gene,
                 args.show_cycles, args.embed_ct)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
