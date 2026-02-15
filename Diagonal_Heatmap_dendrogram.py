from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from numpy import *
import pandas as pd
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import dendrogram, linkage
import argparse


def read_data(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
        if ',' in first_line:
            delimiter = ','
        elif '\t' in first_line:
            delimiter = '\t'
        else:
            raise ValueError("Unable to determine the file format. Please specify the delimiter.")

    df = pd.read_csv(file_path, delimiter=delimiter)
    return df


def draw_rectangle(ax, center_x, center_y, value=None, cmap='viridis', df=None):
    n = 4
    t = arange(0, 360 + (360 / n), 360 / n)
    x = center_x + 10 * sin(radians(t))
    y = center_y + 10 * cos(radians(t))
    ax.text(center_x, center_y, f"{value:.2f}", ha='center', va='center', fontsize=8)
    ax.fill(x, y, edgecolor='black',
            facecolor=get_color(value, cmap, df))


def get_color(value, cmap, df):
    norm = Normalize(vmin=df['relation'].min(),
                     vmax=df['relation'].max())
    return plt.cm.viridis(norm(value))


def draw_pyramid_of_squares(n, wide_df, df):
    for i in range(n):
        for j in range(i + 1):
            if (i % 2 == 0 and j % 2 == 0) or (i % 2 != 0 and j % 2 != 0):
                x = i * 10
                y = j * 10
                value = wide_df.iloc[(n) - (i - j) // 2,
                                     (i + j) // 2]

                draw_rectangle(plt.gca(), x, y,
                               value, cmap='viridis', df=df)

                if j != 0:
                    value = wide_df.iloc[(n) - (i + j) // 2,
                                         (i - j) // 2]

                    draw_rectangle(plt.gca(), x, -y,
                                   value, cmap='viridis', df=df)


def create_wide_df(data):
    unique_vals = pd.concat([data['x1'], data['x2']]).unique()
    wide_df = pd.DataFrame(100,
                           columns=unique_vals,
                           index=unique_vals)

    for index, row in data.iterrows():
        wide_df.loc[row['x1'], row['x2']] = row['relation']
        wide_df.loc[row['x2'], row['x1']] = row['relation']

    return wide_df


def diagonal_heatmap(data,
                     include_dendrogram=False,
                     cmap='viridis',
                     show_cbar=False):

    n = len(pd.concat([data['x1'],
                       data['x2']]).unique())

    wide_df = create_wide_df(data)

    fig = plt.figure(figsize=(12, 10))

    left_heatmap = 0.05
    bottom_heatmap = 0.1
    width_heatmap = 0.5
    height_heatmap = 0.7

    left_dendrogram = 0.8
    bottom_dendrogram = 0.12
    width_dendrogram = 0.5
    height_dendrogram = 0.65

    ax_heatmap = fig.add_axes([left_heatmap,
                               bottom_heatmap,
                               width_heatmap,
                               height_heatmap])

    ax_heatmap.set_axis_off()

    corr_data = wide_df.values

    row_linkage = linkage(corr_data,
                          method='ward',
                          metric='euclidean')

    col_linkage = linkage(corr_data.T,
                          method='ward',
                          metric='euclidean')

    row_order = dendrogram(row_linkage,
                           labels=wide_df.index,
                           no_plot=True,
                           orientation='right')['leaves']

    col_order = dendrogram(col_linkage,
                           no_plot=True)['leaves']

    corr_data = corr_data[row_order][:, col_order]
    df_order = wide_df.iloc[row_order, col_order]

    draw_pyramid_of_squares(n - 1,
                            df_order,
                            data)

    if include_dendrogram:
        ax_dendrogram = fig.add_axes([left_dendrogram,
                                      bottom_dendrogram,
                                      width_dendrogram,
                                      height_dendrogram])

        dendrogram(row_linkage,
                   labels=wide_df.index,
                   orientation='right')
        plt.gca().collections[0].set_linewidth(0.8)

    if show_cbar:
        cbar_ax = fig.add_axes([0.1, 0.05,
                                0.4, 0.02])

        cbar = plt.colorbar(
            plt.cm.ScalarMappable(
                norm=Normalize(vmin=data['relation'].min(),
                               vmax=data['relation'].max()),
                cmap='viridis'),
            cax=cbar_ax,
            orientation='horizontal')

        cbar.set_label('Relation Value')

    plt.savefig("heatmap_and_dendrogram.png",
                dpi=300)

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Generate heatmap with dendrogram.')

    parser.add_argument('file_path',
                        type=str,
                        help='Path to the input data file.')

    args = parser.parse_args()

    df = read_data(args.file_path)

    print("Columns detected:", df.columns)

    # ✅ If using OrthoANI_results.tsv
    # Automatically rename if needed
    if 'Genome1' in df.columns:
        df = df.rename(columns={
            'Genome1': 'x1',
            'Genome2': 'x2',
            'Reciprocal_ANI': 'relation'
        })

    diagonal_heatmap(df,
                     include_dendrogram=True,
                     cmap='viridis',
                     show_cbar=True)


if __name__ == "__main__":
    main()
