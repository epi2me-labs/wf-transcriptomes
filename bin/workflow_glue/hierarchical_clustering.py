"""Hierarchical clustering heatmap plots."""

from bokeh.layouts import column as bokeh_column, row as bokeh_row
from bokeh.models import (
    ColorBar,
    ColumnDataSource,
    Div,
    FixedTicker,
    HoverTool,
    LinearColorMapper,
    Range1d,
    Spacer,
)
from bokeh.palettes import Category10, Category20, RdBu11, Turbo256
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from dominate.util import raw
from ezcharts.plots import BokehPlot
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform


HEATMAP_CELL_HEIGHT = 2
TOP_DENDROGRAM_HEIGHT = 42
TOP_DENDROGRAM_HEADROOM = 0.15
HEATMAP_WIDTH = 200


def clustering_info(data_dtype):
    """Get formatted info text that describe plots."""
    return (
        raw(
            "<b>(Left)</b> Hierarchical clustering heatmap of the top 150 most "
            f"variable {data_dtype} across all samples. Rows represent {data_dtype} "
            "(Z-scored and log2 transformed fold changes), "
            " columns represent samples. Dendrograms show "
            "clustering of both genes and samples. "

            "<b>(Middle)</b> Principal component analysis (PCA showing the first two "
            f"principal components) of sample {data_dtype} expression profiles. "
            " Each point represents a sample, coloured by condition. "
            " Samples that cluster together have similar overall "
            " expression profiles. "

            "<b>(Right)</b> Sample-to-sample Euclidean distance matrix calculated "
            "from log2-transformed fold change values. Lower values "
            " (darker blue) indicate more similar expression "
            " profiles between samples."
            )
    )


def _expression_matrix(data, id_column, samples):
    """Return an expression matrix from a CPM-style table."""
    labels = data[id_column].astype(str).to_numpy()

    metadata = samples.copy()
    sample_columns = metadata['sample'].astype(str)
    sample_columns = [column for column in sample_columns if column in data.columns]

    matrix = data[sample_columns].apply(pd.to_numeric, errors="coerce")
    matrix = matrix.replace([np.inf, -np.inf], np.nan)
    keep = matrix.notna().all(axis=1).to_numpy()
    return matrix.loc[keep].to_numpy(dtype=float), labels[keep], sample_columns


def _top_variable_rows(matrix, labels, top_n):
    """Filter an expression matrix to the top N rows by variance."""
    if len(matrix) <= top_n:
        return matrix, labels
    variances = np.var(matrix, axis=1)
    keep = np.argsort(variances)[-top_n:]
    return matrix[keep], labels[keep]


def _row_zscore(matrix):
    """Z-score each row for heatmap colour scaling."""
    centered = matrix - matrix.mean(axis=1, keepdims=True)
    scale = matrix.std(axis=1, keepdims=True)
    scale[scale == 0] = 1
    return centered / scale


def _cluster_order(matrix):
    """Cluster rows and return linkage data plus leaf order."""
    if len(matrix) < 2:
        return None, np.arange(len(matrix))
    linked = linkage(matrix, method="average", metric="correlation")
    dendro = dendrogram(linked, no_plot=True)
    return linked, np.array(dendro["leaves"])


def _scale_dendrogram_distances(distances):
    """Spread small dendrogram distances for clearer display."""
    return np.sqrt(np.asarray(distances))


def _dendrogram_limit(linked):
    """Return the plotted dendrogram distance limit."""
    if linked is None:
        return 1
    return float(_scale_dendrogram_distances(np.max(linked[:, 2])))


def _dendrogram_source(linked, orientation):
    """Build Bokeh multi-line coordinates from scipy dendrogram output."""
    if linked is None:
        return ColumnDataSource({"xs": [], "ys": []})
    dendro = dendrogram(linked, no_plot=True)
    distances = [
        _scale_dendrogram_distances(segment).tolist()
        for segment in dendro["dcoord"]
    ]
    if orientation == "top":
        xs = [[(x - 5) / 10 for x in segment] for segment in dendro["icoord"]]
        ys = distances
    else:
        xs = distances
        ys = [[(y - 5) / 10 for y in segment] for segment in dendro["icoord"]]
    return ColumnDataSource({"xs": xs, "ys": ys})


def _empty_plot(message):
    """Return a BokehPlot containing a message instead of a heatmap."""
    plot = BokehPlot()
    plot._fig = Div(text=message)
    return plot


def _sample_pca(matrix, sample_names):
    """Project samples onto the first two principal components."""
    sample_matrix = np.asarray(matrix, dtype=float).T
    centered = sample_matrix - sample_matrix.mean(axis=0, keepdims=True)
    if centered.shape[0] < 2 or centered.shape[1] == 0:
        return None

    u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    scores = u * singular_values
    if scores.shape[1] < 2:
        scores = np.column_stack([scores[:, 0], np.zeros(scores.shape[0])])

    total_variance = np.square(singular_values).sum()
    explained = np.square(singular_values[:2]) / total_variance
    return pd.DataFrame({
        "sample": sample_names,
        "pc1": scores[:, 0],
        "pc2": scores[:, 1],
        "pc1_label": f"PC1 ({explained[0] * 100:.1f}%)",
        "pc2_label": f"PC2 ({explained[1] * 100:.1f}%)",
    })


def _sample_distance_data(matrix, sample_names):
    """Return long-form sample-sample Euclidean distances."""
    sample_matrix = np.asarray(matrix, dtype=float).T
    distance_matrix = squareform(pdist(sample_matrix, metric="euclidean"))
    n_samples = len(sample_names)
    x_values = np.tile(np.arange(n_samples), n_samples)
    y_values = np.repeat(np.arange(n_samples), n_samples)
    return pd.DataFrame({
        "x": x_values,
        "y": y_values,
        "sample_x": [sample_names[index] for index in x_values],
        "sample_y": [sample_names[index] for index in y_values],
        "distance": distance_matrix.flatten(),
    })


def _category_palette(size):
    """Return a categorical palette sized for the number of classes."""
    if size <= 10:
        return Category10[10][:size]
    if size <= 20:
        return Category20[20][:size]
    steps = np.linspace(0, len(Turbo256) - 1, num=size, dtype=int)
    return [Turbo256[index] for index in steps]


def _add_colour(samples, condition_column):
    """Align sample metadata and attach a colour for each class."""
    classes = samples[condition_column].drop_duplicates().tolist()
    color_map = dict(zip(classes, _category_palette(len(classes))))
    samples["contrast_color"] = samples[condition_column].map(color_map)
    return samples


def _condition_strip_plot(sample_metadata, col_order, x_range):
    """Build a compact condition strip aligned to the heatmap columns."""
    smeta = sample_metadata.iloc[col_order]
    condition_height = int(TOP_DENDROGRAM_HEIGHT * 0.5)
    strip_source = ColumnDataSource({
        "x": np.arange(len(smeta)),
        "sample": smeta["sample"].tolist(),
        "condition": smeta["condition"].tolist(),
        "sample_color": smeta["contrast_color"].tolist(),
    })
    strip = figure(
        sizing_mode="stretch_width",
        height=condition_height,
        x_range=x_range,
        y_range=Range1d(0, 1),
        tools="",
        toolbar_location=None,
        min_border_left=40,
        min_border_right=0,
        min_border_top=0,
        min_border_bottom=0,
    )
    strip.rect(
        x="x",
        y=0.5,
        width=1,
        height=1,
        source=strip_source,
        fill_color="sample_color",
        line_color=None,
    )
    strip.add_tools(HoverTool(tooltips=[
        ("Sample", "@sample"),
        ("Condition", "@condition"),
    ]))
    strip.axis.visible = False
    strip.grid.visible = False
    strip.outline_line_color = None
    return strip


def _condition_legend_plot(color_map):
    """Build a single legend panel for condition colours."""
    conditions = list(color_map.keys())
    legend_source = ColumnDataSource({
        "x": np.arange(len(conditions)),
        "label": conditions,
        "color": [color_map[condition] for condition in conditions],
    })
    legend_plot = BokehPlot(
        width=max(260, 150 * len(conditions)),
        height=56,
        x_range=Range1d(-0.8, len(conditions) - 0.2),
        y_range=Range1d(0, 1),
        tools=""
    )
    legend = legend_plot._fig
    legend.title.text_font_size = "7pt"
    legend.title.align = "center"
    legend.rect(
        x="x",
        y=0.5,
        width=0.2,
        height=0.36,
        source=legend_source,
        fill_color="color",
        line_color=None,
    )
    legend.text(
        x=np.arange(len(conditions)) + 0.18,
        y=[0.5] * len(conditions),
        text=conditions,
        text_align="left",
        text_baseline="middle",
        text_font_size="10pt",
    )
    legend.axis.visible = False
    legend.grid.visible = False
    legend.outline_line_color = None

    return legend_plot


def _create_title_figure(title_text, height=40):
    """Create a title figure with centered text."""
    title_fig = figure(
        height=height,
        tools="",
        toolbar_location=None,
    )
    title_fig.outline_line_color = None
    title_fig.text(
        x=[5], y=[0.5], text=[title_text],
        text_align="center", text_baseline="middle", text_font_size="11pt")
    title_fig.axis.visible = False
    title_fig.grid.visible = False
    return title_fig


def _blue_white_palette(size=256):
    """Generate a blue to white palette with the specified number of colors."""
    palette = []
    for i in range(size):
        ratio = i / (size - 1) if size > 1 else 0.5
        # Interpolate between blue (#0000FF) and white (#FFFFFF)
        # by increasing red and green from 0 to 255 while keeping blue at 255
        r = int(255 * ratio)
        g = int(255 * ratio)
        b = 255
        palette.append(f'#{r:02x}{g:02x}{b:02x}')
    return palette


def distance_plot(matrix, col_order, sample_names, height):
    """Build the sample-distance heatmap figure."""
    distance_data = _sample_distance_data(matrix[:, col_order], sample_names)
    distance_source = ColumnDataSource(distance_data)
    plot = BokehPlot(
        frame_height=height,
        x_range=Range1d(-0.5, len(sample_names) - 0.5),
        y_range=Range1d(len(sample_names) - 0.5, -0.5),
        tools="",
        toolbar_location=None,
    )
    distance_fig = plot._fig
    palette = _blue_white_palette(256)[::-1]
    distance_mapper = LinearColorMapper(
        palette=palette,
        low=float(distance_data["distance"].min()),
        high=float(distance_data["distance"].max()),
    )
    distance_fig.rect(
        x="x",
        y="y",
        width=1,
        height=1,
        source=distance_source,
        fill_color={"field": "distance", "transform": distance_mapper},
        line_color=None,
    )
    distance_fig.xaxis.ticker = FixedTicker(ticks=list(range(len(sample_names))))
    distance_fig.xaxis.major_label_overrides = {
        index: sample for index, sample in enumerate(sample_names)
    }
    distance_fig.xaxis.major_label_orientation = 1.0
    distance_fig.yaxis.ticker = FixedTicker(ticks=list(range(len(sample_names))))
    distance_fig.yaxis.major_label_overrides = {
        index: sample for index, sample in enumerate(sample_names)
    }
    distance_fig.grid.visible = False
    distance_fig.add_tools(HoverTool(tooltips=[
        ("Sample 1", "@sample_x"),
        ("Sample 2", "@sample_y"),
        ("Distance", "@distance{0.000}"),
    ]))

    min_dist = float(distance_data["distance"].min())
    max_dist = float(distance_data["distance"].max())

    mapper = linear_cmap(field_name='', palette=palette, low=min_dist, high=max_dist)

    color_bar = ColorBar(
        color_mapper=mapper['transform'], width=250, height=15, location=(0, 0))
    color_bar.title = "log2 fold change Euclidian distance"
    color_bar.title_text_font_size = "9pt"

    distance_fig.add_layout(color_bar, 'below')

    title_fig = _create_title_figure("Sample distance")
    plot._fig = bokeh_column(
        title_fig,
        distance_fig,
        sizing_mode="stretch_width"
    )
    return plot


def heatmap_plot(
    z_matrix, labels, sample_names, row_linkage, col_linkage, sample_metadata, col_order
        ):
    """Build the clustered heatmap and row dendrogram layout."""
    n_rows, n_cols = z_matrix.shape
    x_values = np.tile(np.arange(n_cols), n_rows)
    y_values = np.repeat(np.arange(n_rows), n_cols)
    heat_source = ColumnDataSource({
        "x": x_values,
        "y": y_values,
        "sample": [sample_names[index] for index in x_values],
        "feature": [labels[index] for index in y_values],
        "value": z_matrix.flatten(),
    })
    mapper = LinearColorMapper(
        palette=list(reversed(RdBu11)),
        low=-2,
        high=2,
    )
    heatmap_height = max(160, HEATMAP_CELL_HEIGHT * n_rows)

    heatmap = figure(
        sizing_mode="stretch_width",
        frame_height=heatmap_height,
        x_range=Range1d(-0.5, n_cols - 0.5),
        y_range=Range1d(n_rows - 0.5, -0.5),
        tools="",
        toolbar_location=None,
    )
    heatmap.rect(
        x="x",
        y="y",
        width=1,
        height=1,
        source=heat_source,
        fill_color={"field": "value", "transform": mapper},
        line_color=None,
    )
    heatmap.add_tools(HoverTool(tooltips=[
        ("Feature", "@feature"),
        ("Sample", "@sample"),
        ("Row z-score", "@value{0.000}"),
    ]))
    heatmap.xaxis.ticker = FixedTicker(ticks=list(range(n_cols)))
    heatmap.xaxis.major_label_overrides = {
        index: sample for index, sample in enumerate(sample_names)
    }
    heatmap.xaxis.major_label_orientation = 1.0
    heatmap.yaxis.ticker = FixedTicker(ticks=list(range(n_rows)))
    heatmap.yaxis.visible = False
    heatmap.min_border_left = 40
    heatmap.min_border_right = 0
    heatmap.min_border_top = 0
    heatmap.grid.visible = False

    feature_dendro_fig = figure(
        width=84,
        frame_height=heatmap_height,
        x_range=Range1d(0, _dendrogram_limit(row_linkage)),
        y_range=heatmap.y_range,
        toolbar_location=None,
        tools="",
        min_border_left=0,
        min_border_right=0,
        min_border_top=0,
        min_border_bottom=0,
    )
    feature_dendro_fig.multi_line(
        xs="xs",
        ys="ys",
        source=_dendrogram_source(row_linkage, "right"),
        line_color="#333333",
        line_width=1,
    )
    feature_dendro_fig.axis.visible = False
    feature_dendro_fig.grid.visible = False

    # top sample dendro plot
    top_limit = _dendrogram_limit(col_linkage)
    sample_dendro_fig = figure(
        frame_height=TOP_DENDROGRAM_HEIGHT,
        x_range=heatmap.x_range,
        y_range=Range1d(0, top_limit * (1 + TOP_DENDROGRAM_HEADROOM)),
        toolbar_location=None,
        min_border_left=40,
        min_border_right=0,
        min_border_top=0,
        min_border_bottom=0,
    )
    sample_dendro_fig.multi_line(
        xs="xs",
        ys="ys",
        source=_dendrogram_source(col_linkage, "top"),
        line_color="#333333",
        line_width=1,
    )
    sample_dendro_fig.axis.visible = False
    sample_dendro_fig.grid.visible = False

    final_fig = BokehPlot()
    if sample_metadata is not None:
        strip_height = int(TOP_DENDROGRAM_HEIGHT * 0.5)
        strip = _condition_strip_plot(sample_metadata, col_order, heatmap.x_range)
        left_column = bokeh_column(
            sample_dendro_fig, strip, heatmap, sizing_mode="stretch_width", spacing=0)
        right_column = bokeh_column(
            Spacer(height=TOP_DENDROGRAM_HEIGHT + strip_height, width=84),
            feature_dendro_fig,
            spacing=0,
        )
        title_fig = _create_title_figure("Hierarchical clustering")
        final_fig._fig = bokeh_column(
            title_fig,
            bokeh_row(
                left_column, right_column, sizing_mode="stretch_width", spacing=0),
            sizing_mode="stretch_width"
        )
    else:
        final_fig._fig = bokeh_row(
            heatmap, feature_dendro_fig, sizing_mode="stretch_width", spacing=0)
    return final_fig, heatmap_height


def pca_plot(matrix, col_order, sample_names, sample_metadata, condition_column):
    """Build the sample PCA figure."""
    pca_data = _sample_pca(matrix[:, col_order], sample_names)
    if sample_metadata is not None:
        pca_data = pca_data.merge(sample_metadata, on="sample", how="left")
    else:
        pca_data["contrast_color"] = "#4C78A8"

    plot = BokehPlot(
        tools="",
        height=360
    )
    pca = plot._fig
    pca_source = ColumnDataSource(pca_data)
    tooltips = [
        ("Sample", "@sample"),
        ("PC1", "@pc1{0.000}"),
        ("PC2", "@pc2{0.000}"),
    ]
    if "condition" in pca_data.columns:
        tooltips.insert(1, ("Condition", "@condition"))
    pca.scatter(
        x="pc1",
        y="pc2",
        size=7,
        source=pca_source,
        color="contrast_color",
        line_color="contrast_color",
        fill_alpha=0.65,
    )
    pca.add_tools(HoverTool(tooltips=tooltips))
    x_min = float(pca_data["pc1"].min())
    x_max = float(pca_data["pc1"].max())
    y_min = float(pca_data["pc2"].min())
    y_max = float(pca_data["pc2"].max())
    x_pad = max((x_max - x_min) * 0.08, 0.1)
    y_pad = max((y_max - y_min) * 0.08, 0.1)
    pca.x_range = Range1d(x_min - x_pad, x_max + x_pad)
    pca.y_range = Range1d(y_min - y_pad, y_max + y_pad)
    pca.xaxis.axis_label = pca_data["pc1_label"].iat[0]
    pca.yaxis.axis_label = pca_data["pc2_label"].iat[0]
    pca.grid.grid_line_alpha = 0.3

    title_fig = _create_title_figure("Sample PCA")

    color_map = dict(zip(
        sample_metadata[condition_column].tolist(),
        sample_metadata['contrast_color'].tolist()),
    )
    legend = _condition_legend_plot(color_map)
    pca_and_legend = BokehPlot()
    pca_and_legend._fig = bokeh_column(
        title_fig, pca, legend._fig, sizing_mode="stretch_width", spacing=5)

    return pca_and_legend


def hierarchical(
    data,
    id_column,
    samples,
    condition_column,
    top_n=150,
):
    """Build a clustered expression heatmap with dendrograms and sample PCA."""
    samples.rename(columns={'alias': 'sample'}, inplace=True)
    log2_matrix, labels, sample_names = _expression_matrix(
        data,
        id_column=id_column,
        samples=samples,
    )

    if log2_matrix.shape[0] < 2 or log2_matrix.shape[1] < 2:
        return _empty_plot(
            "Hierarchical clustering requires at least two features and "
            "two numeric sample columns."
        )

    log2_matrix = np.log2(log2_matrix + 1)
    log2_matrix, labels = _top_variable_rows(log2_matrix, labels, top_n=top_n)
    log2_z_matrix = _row_zscore(log2_matrix)

    row_linkage, row_order = _cluster_order(log2_z_matrix)
    col_linkage, col_order = _cluster_order(log2_z_matrix.T)
    log2_z_matrix = log2_z_matrix[row_order][:, col_order]
    labels = labels[row_order]
    sample_names = [sample_names[index] for index in col_order]
    meta = _add_colour(
        samples,
        condition_column
    )

    heatmap_plt, hm_height = heatmap_plot(
        z_matrix=log2_z_matrix,
        labels=labels,
        sample_names=sample_names,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        sample_metadata=meta,
        col_order=col_order
    )
    pca_plt = pca_plot(
        matrix=log2_matrix,
        col_order=col_order,
        sample_names=sample_names,
        sample_metadata=meta,
        condition_column=condition_column,
    )
    distance_plt = distance_plot(
        matrix=log2_matrix,
        col_order=col_order,
        sample_names=sample_names,
        height=hm_height
    )

    return heatmap_plt, pca_plt, distance_plt
