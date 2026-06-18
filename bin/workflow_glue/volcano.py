"""Interactive volcano plot visualization for differential expression analysis."""

from bokeh.events import DocumentReady, Reset, Tap
from bokeh.layouts import column as bokeh_column, row as bokeh_row
from bokeh.models import (
    AutocompleteInput,
    Button,
    ColumnDataSource,
    CustomJS,
    CustomJSHover,
    DataTable,
    HoverTool,
    HTMLTemplateFormatter,
    LabelSet,
    NumberFormatter,
    Range1d,
    Slider,
    Span,
    TableColumn,
    TextInput,
    Toggle
)
from bokeh.plotting import figure
from ezcharts.plots import BokehPlot
import numpy as np
import pandas as pd

SIGNIFICANCE_CLASSES = {
    "significant_high_effect": {
        "color": "#d62728",
        "marker": "circle",
        "symbol": "&#9679;",
        "label": "Significant, high effect",
    },
    "significant_low_effect": {
        "color": "#1f77b4",
        "marker": "square",
        "symbol": "&#9632;",
        "label": "Significant, low effect",
    },
    "non_significant_high_effect": {
        "color": "#2ca02c",
        "marker": "triangle",
        "symbol": "&#9650;",
        "label": "Not significant, high effect",
    },
    "non_significant_low_effect": {
        "color": "#8a8a8a",
        "marker": "x",
        "symbol": "&#10005;",
        "label": "Not significant, low effect",
    },
}


def _class_attr(attr,):
    return {k: v[attr] for k, v in SIGNIFICANCE_CLASSES.items()}


def _class_attr_list(attr):
    return [v[attr] for k, v in SIGNIFICANCE_CLASSES.items()]


def _padded_range(start, end, padding_fraction=0.02, min_padding=0.1):
    """Return a lightly padded numeric range for plot axes."""
    span = end - start
    padding = max(min_padding, span * padding_fraction)
    return start - padding, end + padding


def _padded_log_range(start, end, padding_fraction=0.02):
    """Return a lightly padded positive range for log-scaled axes."""
    factor = 1 + padding_fraction
    return start / factor, end * factor


def _volcano_source_data(data, fold_threshold=1, p_threshold=0.05):
    """Prepare the minimal browser payload for the volcano/MA plot."""
    data = data.copy()
    input_columns = set(data.columns)
    data["log2FoldChange"] = pd.to_numeric(data["log2FoldChange"], errors="coerce")
    data["padj"] = pd.to_numeric(data["padj"], errors="coerce")

    data.replace([np.inf, -np.inf], np.nan, inplace=True)
    data.dropna(subset=["log2FoldChange", "padj"], inplace=True)

    is_transcript_plot = "featureID" in data.columns
    if is_transcript_plot:
        data.rename(columns={
            "featureID": "TXNAME",
            "groupID": "GENEID",
            "exonBaseMean": "mean_expression"
            },
            inplace=True)
        identifier_column = "TXNAME"

    else:
        identifier_column = "GENEID"
        data.rename(columns={
            "baseMean": "mean_expression"
        }, inplace=True)

    if identifier_column not in data.columns:
        raise ValueError(
            "Volcano plot data must contain 'TXNAME' for transcript plots or "
            "'GENEID' for gene plots."
        )

    if "mean_expression" not in data.columns:
        raise ValueError(
            "Volcano/MA plot data must contain 'baseMean' or 'exonBaseMean' column."
        )

    data["mean_expression"] = pd.to_numeric(
        data['mean_expression'], errors="coerce"
    )

    has_gene_name = "gene_name" in data.columns
    if not has_gene_name:
        data["gene_name"] = ""

    data["gene_group"] = data["GENEID"]
    data["selected_label"] = ""
    data["label_x_offset"] = 7
    data["label_align"] = "left"

    # Transformations and thresholds
    data["neg_log10_padj"] = -np.log10(data["padj"])
    y_threshold = -np.log10(p_threshold)
    outside_fold_threshold = data["log2FoldChange"].abs() >= fold_threshold
    above_p_threshold = data["neg_log10_padj"] >= y_threshold

    # Setup initial significance classes
    data["volcano_class"] = "non_significant_low_effect"
    data.loc[
        outside_fold_threshold & (data["neg_log10_padj"] < y_threshold),
        "volcano_class",
    ] = "non_significant_high_effect"
    data.loc[
        outside_fold_threshold & above_p_threshold,
        "volcano_class",
    ] = "significant_high_effect"
    data.loc[
        ~outside_fold_threshold & above_p_threshold,
        "volcano_class",
    ] = "significant_low_effect"
    data["color"] = data["volcano_class"].map(_class_attr("color"))
    data["marker"] = data["volcano_class"].map(_class_attr("marker"))

    # All the columns we need to send to the browser for interactivity and export
    source_columns = [
        identifier_column,
        "gene_group",
        "gene_name",
        "log2FoldChange",
        "mean_expression",
        "marker",
        "neg_log10_padj",
        "padj",
        "selected_label",
        "label_x_offset",
        "label_align",
        "color",
        "volcano_class",
    ]

    if is_transcript_plot:
        source_columns.insert(1, "GENEID")

    source_data = {}
    for column in source_columns:
        values = data[column]
        if column in (
            "log2FoldChange",
            "mean_expression",
            "padj",
            "neg_log10_padj",
        ):
            values = values.astype(np.float32)
        source_data[column] = values.reset_index(drop=True)

    export_columns = [identifier_column]
    if is_transcript_plot:
        export_columns.append("GENEID")
    if "gene_name" in input_columns:
        export_columns.append("gene_name")
    export_columns.extend([
        "log2FoldChange",
        "mean_expression",
        "padj",
        "volcano_class",
    ])

    return source_data, export_columns


def _tap_selection_callback_code(x_field, y_field):
    """Build the shared JS tap-selection callback for volcano-style plots."""
    return f"""
        const data = source.data;
        const tapX = cb_obj.sx;
        const tapY = cb_obj.sy;
        const maxPixelDistanceSquared = 20 * 20;
        function findPlotView(view, plotId) {{
            if (!view) {{
                return null;
            }}
            if (view.model && view.model.id === plotId && view.frame) {{
                return view;
            }}
            const childViews = view.child_views || view._child_views;
            if (!childViews) {{
                return null;
            }}
            const children = childViews.values ? childViews.values() : childViews;
            for (const child of children) {{
                const found = findPlotView(child, plotId);
                if (found) {{
                    return found;
                }}
            }}
            return null;
        }}
        function selectNearestPointByPixelDistance() {{
            let plotView = Bokeh.index[plot.id];
            if (!plotView || !plotView.frame) {{
                for (const rootView of Object.values(Bokeh.index)) {{
                    plotView = findPlotView(rootView, plot.id);
                    if (plotView) {{
                        break;
                    }}
                }}
            }}
            if (
                !plotView
                || !Number.isFinite(tapX)
                || !Number.isFinite(tapY)
            ) {{
                return;
            }}
            const xScale = plotView.frame.x_scale;
            const yScale = plotView.frame.y_scale;
            let bestIndex = -1;
            let bestDistance = Infinity;
            let i = 0;
            while (i !== data.{x_field}.length) {{
                const pointX = xScale.compute(data.{x_field}[i]);
                const pointY = yScale.compute(data.{y_field}[i]);
                if (!Number.isFinite(pointX) || !Number.isFinite(pointY)) {{
                    i += 1;
                    continue;
                }}
                const dx = pointX - tapX;
                const dy = pointY - tapY;
                const distance = dx * dx + dy * dy;
                if (Math.sign(bestDistance - distance) === 1) {{
                    bestDistance = distance;
                    bestIndex = i;
                }}
                i += 1;
            }}
            if (bestIndex !== -1 && bestDistance <= maxPixelDistanceSquared) {{
                if (gene_select_toggle.active) {{
                    const geneGroup = data.gene_group[bestIndex];
                    const selected = [];
                    let j = 0;
                    while (j !== data.gene_group.length) {{
                        if (data.gene_group[j] === geneGroup) {{
                            selected.push(j);
                        }}
                        j += 1;
                    }}
                    selection_state.data.indices[0] = selected;
                }} else {{
                    const selected = new Set(selection_state.data.indices[0]);
                    selected.add(bestIndex);
                    selection_state.data.indices[0] = Array.from(selected).sort(
                        function(a, b) {{
                        return a - b;
                    }});
                }}
                update_callback.execute(source);
            }}
        }}
        selectNearestPointByPixelDistance();
        """


def volcano(data, fold_threshold=1, p_threshold=0.05):
    """Build an interactive volcano plot with selection and filtering widgets."""
    input_size = len(data)
    source_data = data.dropna(subset=["log2FoldChange", "padj"])
    n_filtered = input_size - len(source_data)
    source_data, original_columns = _volcano_source_data(
        data, fold_threshold=fold_threshold, p_threshold=p_threshold
    )
    is_transcript_plot = "TXNAME" in source_data
    identifier_col = "TXNAME" if is_transcript_plot else "GENEID"
    has_gene_name = "gene_name" in data.columns
    counts = pd.Series(source_data["volcano_class"]).value_counts().to_dict()
    x_min = float(source_data["log2FoldChange"].min())
    x_max = float(source_data["log2FoldChange"].max())
    y_max = float(source_data["neg_log10_padj"].max())
    mean_min = float(source_data["mean_expression"].min())
    mean_max = float(source_data["mean_expression"].max())
    if x_min == x_max:
        x_min -= 1
        x_max += 1
    x_range_min, x_range_max = _padded_range(x_min, x_max)
    if y_max <= 0:
        y_max = 1
    _, y_axis_max = _padded_range(0, y_max, min_padding=0.2)
    if mean_min == mean_max:
        mean_min *= 0.9
        mean_max *= 1.1
    # Log scale requires strictly positive lower bound; baseMean can be 0
    _pos_means = source_data["mean_expression"][source_data["mean_expression"] > 0]
    mean_min_log = float(_pos_means.min()) if len(_pos_means) > 0 else 0.01
    ma_x_range_min, ma_x_range_max = _padded_log_range(mean_min_log, mean_max)
    y_threshold = -np.log10(p_threshold)
    responsive_slider_title_stylesheet = """
    @media screen and (max-width: 1100px) {
        .bk-slider-title {
            display: block;
            max-width: 100%;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }
        .bk-slider-title:hover {
            overflow: visible;
            text-overflow: clip;
            white-space: normal;
            position: relative;
            z-index: 1;
            background: white;
        }
    }
    """

    source = ColumnDataSource(source_data)
    selection_state = ColumnDataSource({"indices": [[]]})
    highlight_source = ColumnDataSource({
        "log2FoldChange": [],
        "neg_log10_padj": [],
        "mean_expression": [],
        "gene_name": [],
    })
    selected_source = ColumnDataSource({
        "source_index": [],
        "owner_id": [],
        "remove": [],
        "GENEID": [],
        "TXNAME": [],
        "gene_name": [],
        "log2FoldChange": [],
        "padj": [],
        "volcano_class": [],
    })
    legend_source = ColumnDataSource({
        "symbol": _class_attr_list("symbol"),
        "color": _class_attr_list("color"),
        "label": _class_attr_list("label"),
        "count": [counts.get(k, 0) for k in SIGNIFICANCE_CLASSES],
    })

    hover_content_formatter = CustomJSHover(
        args=dict(
            source=source,
            identifier_col=identifier_col,
            is_transcript_plot=is_transcript_plot,
            has_gene_name=has_gene_name,
        ),
        code="""
        const n = special_vars.indices.length;
        if (n > 3) {
            if (special_vars.index !== special_vars.indices[0]) { return ''; }
            return (
                '<em style="color:#888;">' +
                n + ' points \u2014 zoom in to resolve' +
                '</em>'
            );
        }
        if (!n) { return ''; }
        const data = source.data;

        function fmtPadj(v) {
            const num = Number(v);
            if (!isFinite(num)) { return String(v); }
            if (num < 0.0001) { return num.toExponential(4); }
            return num.toFixed(4);
        }

        const th = 'style="text-align:right;padding-right:4px;color:#666;"';
        function row(label, value, valStyle) {
            const td = valStyle
                ? '<td style="' + valStyle + '">' + value + '</td>'
                : '<td>' + value + '</td>';
            return '<tr><th ' + th + '>' + label + '</th>' + td + '</tr>';
        }
        function renderPoint(i) {
            const idLabel = is_transcript_plot ? 'TX' : 'Gene';
            const idValue = is_transcript_plot
                ? String(data.TXNAME ? data.TXNAME[i] || '' : '')
                : (has_gene_name
                    ? String(data.gene_name[i] || '')
                    : String(data.GENEID ? data.GENEID[i] || '' : ''));
            let h = '<table style="border-spacing:2px 1px">';
            h += row(idLabel, idValue, 'color:#6bb5d6;font-weight:bold;');
            if (is_transcript_plot) {
                const gene = has_gene_name
                    ? String(data.gene_name[i] || '')
                    : String(data.GENEID ? data.GENEID[i] || '' : '');
                h += row('Gene', gene);
            }
            h += row('log2FC', Number(data.log2FoldChange[i]).toFixed(3));
            h += row('Mean expr', Number(data.mean_expression[i]).toFixed(3));
            h += row('padj', fmtPadj(data.padj[i]));
            h += row('Class', String(data.volcano_class[i] || ''));
            h += '</table>';
            return h;
        }

        return renderPoint(special_vars.indices[0]);
        """,
    )
    hover = HoverTool(
        tooltips=f"@{identifier_col}{{custom}}",
        formatters={f"@{identifier_col}": hover_content_formatter},
    )

    vol_fig = figure(
        x_axis_label="log2 fold change",
        y_axis_label="-log10 adjusted p-value",
        height=480,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,reset,save",
        x_range=Range1d(
            x_range_min, x_range_max, bounds=(x_range_min, x_range_max)),
        y_range=Range1d(0, y_axis_max, bounds=(0, y_axis_max)),
    )
    vol_fig.toolbar.logo = None

    points_renderer = vol_fig.scatter(
        x="log2FoldChange",
        y="neg_log10_padj",
        marker="marker",
        source=source,
        color="color",
        size=7,
        fill_alpha=0.25,
        line_alpha=0.40,
        selection_color="#111111",
        selection_alpha=0.9,
        nonselection_alpha=0.12,
    )
    hover.renderers = [points_renderer]
    hover.callback = CustomJS(
        args=dict(selected_source=selected_source),
        code="""
        const hovered = cb_data.index.indices;
        if (!hovered.length || hovered.length > 3) {
            selected_source.selected.indices = [];
            selected_source.selected.change.emit();
            selected_source.change.emit();
            return;
        }
        const sourceIndex = hovered[0];
        const selectedRows = [];
        let i = 0;
        while (i !== selected_source.data.source_index.length) {
            if (selected_source.data.source_index[i] === sourceIndex) {
                selectedRows.push(i);
                break;
            }
            i += 1;
        }
        selected_source.selected.indices = selectedRows;
        selected_source.selected.change.emit();
        selected_source.change.emit();
        """,
    )
    vol_fig.add_tools(hover)
    vol_fig.circle(
        x="log2FoldChange",
        y="neg_log10_padj",
        source=highlight_source,
        size=15,
        fill_alpha=0,
        line_alpha=1,
        line_color="#111111",
        line_width=2,
    )
    labels = LabelSet(
        x="log2FoldChange",
        y="neg_log10_padj",
        text="selected_label",
        source=source,
        x_offset="label_x_offset",
        y_offset=7,
        text_align="label_align",
        text_font_size="10px",
        text_color="#111111",
    )
    vol_fig.add_layout(labels)

    ma_hover = HoverTool(
        tooltips=f"@{identifier_col}{{custom}}",
        formatters={f"@{identifier_col}": hover_content_formatter},
    )
    ma_fig = figure(
        x_axis_label="mean expression",
        y_axis_label="log2 fold change",
        height=480,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,reset,save",
        x_axis_type="log",
        x_range=Range1d(
            ma_x_range_min, ma_x_range_max,
            bounds=(ma_x_range_min, ma_x_range_max),
        ),
        y_range=Range1d(x_range_min, x_range_max, bounds=(x_range_min, x_range_max)),
    )
    ma_fig.toolbar.logo = None
    ma_fig.visible = False
    ma_points_renderer = ma_fig.scatter(
        x="mean_expression",
        y="log2FoldChange",
        marker="marker",
        source=source,
        color="color",
        size=7,
        fill_alpha=0.25,
        line_alpha=0.40,
        selection_color="#111111",
        selection_alpha=0.9,
        nonselection_alpha=0.12,
    )
    ma_hover.renderers = [ma_points_renderer]
    ma_hover.callback = hover.callback
    ma_fig.add_tools(ma_hover)
    ma_fig.circle(
        x="mean_expression",
        y="log2FoldChange",
        source=highlight_source,
        size=15,
        fill_alpha=0,
        line_alpha=1,
        line_color="#111111",
        line_width=2,
    )
    ma_labels = LabelSet(
        x="mean_expression",
        y="log2FoldChange",
        text="selected_label",
        source=source,
        x_offset=7,
        y_offset=7,
        text_font_size="10px",
        text_color="#111111",
    )
    ma_fig.add_layout(ma_labels)
    ma_zero_line = Span(
        location=0,
        dimension="width",
        line_dash="dashed",
        line_color="#444444",
        line_width=1,
    )
    ma_fig.add_layout(ma_zero_line)

    # Setup the selection table
    table_columns = [
        TableColumn(
            field="remove",
            title="",
            width=28,
            sortable=False,
            formatter=HTMLTemplateFormatter(
                template=(
                    '<button class="volcano-row-delete" '
                    'data-selected-source-id="<%= owner_id %>" '
                    'data-source-index="<%= source_index %>" '
                    'style="border:0;background:transparent;color:#666;'
                    'cursor:pointer;font-size:16px;line-height:1;'
                    'padding:0 4px;" title="Remove row">&times;</button>'
                ),
            ),
        ),
    ]
    if is_transcript_plot:
        table_columns.extend([
            TableColumn(field="TXNAME", title="TXNAME"),
            TableColumn(field="GENEID", title="GENEID"),
        ])
    else:
        table_columns.append(TableColumn(field="GENEID", title="GENEID"))
    if has_gene_name:
        table_columns.append(TableColumn(field="gene_name", title="gene_name"))
    table_columns.extend([
        TableColumn(
            field="log2FoldChange",
            title="log2FC",
            formatter=NumberFormatter(format="0.000"),
        ),
        TableColumn(
            field="padj",
            title="padj",
            formatter=HTMLTemplateFormatter(
                template="""
                <% if (value < 0.001 && value !== 0) { %>
                    <span><%= value.toExponential(2) %></span>
                <% } else { %>
                    <span><%= value.toFixed(4) %></span>
                <% } %>
                """
                ),
        ),
    ])
    selected_table = DataTable(
        source=selected_source,
        columns=table_columns,
        autosize_mode="force_fit",
        height=355,
        index_position=None,
        visible=False,
    )
    legend_table = DataTable(
        source=legend_source,
        columns=[
            TableColumn(
                field="label",
                title="Class",
                width=150,
                formatter=HTMLTemplateFormatter(
                    template=(
                        '<span style="color:<%= color %>;font-size:13px;'
                        'margin-right:4px;"><%= symbol %></span><%= value %>'
                    ),
                ),
            ),
            TableColumn(field="count", title="n", width=50),
        ],
        height=130,
        autosize_mode="force_fit",
        index_position=None,
        editable=False,
        sortable=False,
        styles={"padding-bottom": "4px"},
    )
    view_toggle = Toggle(
        label="Toggle highlights",
        button_type="default",
        active=True,
        width=120,
    )
    gene_select_toggle = Toggle(
        label="Select by gene",
        button_type="default",
        active=False,
        width=120,
        visible=is_transcript_plot,
    )
    plot_mode_toggle = Button(
        label="Show MA plot",
        button_type="default",
        width=120,
    )
    left_fold_line = Span(
        location=-fold_threshold,
        dimension="height",
        line_dash="dashed",
        line_color="#444444",
        line_width=1,
    )
    right_fold_line = Span(
        location=fold_threshold,
        dimension="height",
        line_dash="dashed",
        line_color="#444444",
        line_width=1,
    )
    p_line = Span(
        location=y_threshold,
        dimension="width",
        line_dash="dashed",
        line_color="#444444",
        line_width=1,
    )
    vol_fig.add_layout(left_fold_line)
    vol_fig.add_layout(right_fold_line)
    vol_fig.add_layout(p_line)

    fold_slider = Slider(
        title="Absolute log2 fold-change threshold",
        value=fold_threshold,
        start=0,
        end=max(5, np.ceil(np.abs(source_data["log2FoldChange"]).max())),
        step=0.1,
        show_value=False,
        sizing_mode="stretch_width",
        stylesheets=[responsive_slider_title_stylesheet],
    )
    fold_input = TextInput(
        title="",
        value=f"{fold_threshold:.1f}",
        width=80,
    )
    p_slider = Slider(
        title="Adjusted p-value threshold",
        value=y_threshold,
        start=0,
        end=max(5, np.ceil(y_max)),
        step=0.1,
        show_value=False,
        sizing_mode="stretch_width",
        stylesheets=[responsive_slider_title_stylesheet],
    )
    p_input = TextInput(
        title="",
        value=f"{p_threshold:.3f}" if p_threshold >= 0.001 else f"{p_threshold:.2e}",
        width=80,
    )
    y_axis_slider = Slider(
        title="y-axis maximum",
        value=y_axis_max,
        start=0.1,
        end=max(5, np.ceil(y_axis_max)),
        step=0.1,
        show_value=False,
        orientation="horizontal",
        sizing_mode="stretch_width",
        stylesheets=[responsive_slider_title_stylesheet],
    )
    # Setup callbacks for interactivity
    update_callback = CustomJS(
        args=dict(
            source=source,
            vol_fig=vol_fig,
            selection_state=selection_state,
            highlight_source=highlight_source,
            selected_source=selected_source,
            selected_table=selected_table,
            legend_source=legend_source,
            fold_slider=fold_slider,
            fold_input=fold_input,
            p_slider=p_slider,
            p_input=p_input,
            left_fold_line=left_fold_line,
            right_fold_line=right_fold_line,
            p_line=p_line,
            colors=_class_attr("color"),
            markers=_class_attr("marker"),
            labels=_class_attr("label"),
            view_toggle=view_toggle,
            identifier_col=identifier_col,
        ),
        code="""
        const data = source.data;
        const fold = fold_slider.value;
        const yThreshold = p_slider.value;
        const pThreshold = Math.pow(10, -yThreshold);
        fold_input.value = fold.toFixed(1);
        const pThresholdLabel = pThreshold >= 0.001
            ? pThreshold.toFixed(3)
            : pThreshold.toExponential(2);
        p_input.value = pThresholdLabel;
        const log2fc = data.log2FoldChange;
        const negLog10Padj = data.neg_log10_padj;
        const color = data.color;
        const marker = data.marker;
        const volcanoClass = data.volcano_class;
        const selectedLabel = data.selected_label;
        const labelXOffset = data.label_x_offset;
        const labelAlign = data.label_align;
        const counts = {
            significant_high_effect: 0,
            significant_low_effect: 0,
            non_significant_high_effect: 0,
            non_significant_low_effect: 0,
        };

        let i = 0;
        while (i !== log2fc.length) {
            const outsideFold = Math.sign(Math.abs(log2fc[i]) - fold) !== -1;
            const aboveP = Math.sign(negLog10Padj[i] - yThreshold) !== -1;
            let klass;
            if (aboveP) {
                if (outsideFold) {
                    klass = "significant_high_effect";
                } else {
                    klass = "significant_low_effect";
                }
            } else {
                if (outsideFold) {
                    klass = "non_significant_high_effect";
                } else {
                    klass = "non_significant_low_effect";
                }
            }
            volcanoClass[i] = klass;
            counts[klass] += 1;
            color[i] = colors[klass];
            marker[i] = markers[klass];
            i += 1;
        }
        const selected = Array.from(
            new Set(selection_state.data.indices[0])
        ).sort(function(a, b) {
            return a - b;
        });
        selection_state.data.indices[0] = selected;
        source.selected.indices = [];

        const selectedData = {
            source_index: [],
            owner_id: [],
            remove: [],
            GENEID: [],
            TXNAME: [],
            gene_name: [],
            log2FoldChange: [],
            padj: [],
            volcano_class: [],
        };
        i = 0;
        while (i !== selectedLabel.length) {
            selectedLabel[i] = "";
            labelXOffset[i] = 7;
            labelAlign[i] = "left";
            i += 1;
        }
        const xRangeStart = vol_fig.x_range.start;
        const xRangeEnd = vol_fig.x_range.end;
        const xSpan = xRangeEnd - xRangeStart;
        const labelFlipThreshold = xRangeEnd - (xSpan * 0.12);
        selected.forEach(function(index) {
            if (view_toggle.active) {
                selectedLabel[index] = identifier_col === "TXNAME"
                    ? data.TXNAME[index]
                    : (data.gene_name[index] || data[identifier_col][index]);
                if (data.log2FoldChange[index] >= labelFlipThreshold) {
                    labelXOffset[index] = -7;
                    labelAlign[index] = "right";
                }
            }
            selectedData.source_index.push(index);
            selectedData.owner_id.push(selected_source.id);
            selectedData.remove.push("x");
            selectedData.GENEID.push(data.GENEID ? data.GENEID[index] : "");
            selectedData.TXNAME.push(data.TXNAME ? data.TXNAME[index] : "");
            selectedData.gene_name.push(data.gene_name[index]);
            selectedData.log2FoldChange.push(data.log2FoldChange[index]);
            selectedData.padj.push(data.padj[index]);
            selectedData.volcano_class.push(data.volcano_class[index]);
        });
        selected_source.data = selectedData;
        selected_table.visible = selected.length > 0;
        selected_table.height = window.innerWidth <= 1100 ? 180 : 355;
        if (view_toggle.active) {
            const highlightData = {
                log2FoldChange: [],
                neg_log10_padj: [],
                mean_expression: [],
                gene_name: [],
            };
            selected.forEach(function(index) {
                highlightData.log2FoldChange.push(data.log2FoldChange[index]);
                highlightData.neg_log10_padj.push(data.neg_log10_padj[index]);
                highlightData.mean_expression.push(data.mean_expression[index]);
                highlightData.gene_name.push(data.gene_name[index]);
            });
            highlight_source.data = highlightData;
        }

        left_fold_line.location = -fold;
        right_fold_line.location = fold;
        p_line.location = yThreshold;
        legend_source.data.count = [
            counts.significant_high_effect,
            counts.significant_low_effect,
            counts.non_significant_high_effect,
            counts.non_significant_low_effect,
        ];
        legend_source.change.emit();
        source.change.emit();
        selection_state.change.emit();
        highlight_source.change.emit();
        selected_source.change.emit();
        """,
    )
    vol_fig.js_on_event(Tap, CustomJS(
        args=dict(
            plot=vol_fig,
            source=source,
            selection_state=selection_state,
            gene_select_toggle=gene_select_toggle,
            update_callback=update_callback,
        ),
        code=_tap_selection_callback_code("log2FoldChange", "neg_log10_padj"),
    ))
    ma_fig.js_on_event(Tap, CustomJS(
        args=dict(
            plot=ma_fig,
            source=source,
            selection_state=selection_state,
            gene_select_toggle=gene_select_toggle,
            update_callback=update_callback,
        ),
        code=_tap_selection_callback_code("mean_expression", "log2FoldChange"),
    ))
    fold_slider.js_on_change("value", update_callback)
    fold_input.js_on_change("value", CustomJS(
        args=dict(
            fold_input=fold_input,
            fold_slider=fold_slider,
            update_callback=update_callback,
            source=source,
        ),
        code="""
        const rawValue = fold_input.value.trim();
        if (!rawValue) {
            return;
        }
        const parsed = Number(rawValue);
        if (!Number.isFinite(parsed) || parsed < 0) {
            return;
        }
        fold_slider.value = parsed;
        update_callback.execute(source);
        """,
    ))
    p_slider.js_on_change("value", update_callback)
    p_input.js_on_change("value", CustomJS(
        args=dict(
            p_input=p_input,
            p_slider=p_slider,
            update_callback=update_callback,
            source=source,
        ),
        code="""
        const rawValue = p_input.value.trim();
        if (!rawValue) {
            return;
        }
        const parsed = Number(rawValue);
        if (!Number.isFinite(parsed) || parsed <= 0 || parsed > 1) {
            return;
        }
        p_slider.value = -Math.log10(parsed);
        update_callback.execute(source);
        """,
    ))
    y_axis_slider.js_on_change("value", CustomJS(
        args=dict(vol_fig=vol_fig, y_axis_slider=y_axis_slider),
        code="""
        vol_fig.y_range.end = y_axis_slider.value;
        """,
    ))
    vol_fig.js_on_event(Reset, CustomJS(
        args=dict(vol_fig=vol_fig, y_axis_slider=y_axis_slider),
        code="""
        y_axis_slider.value = vol_fig.y_range.end;
        """,
    ))
    view_toggle.js_on_change("active", CustomJS(
        args=dict(
            source=source,
            vol_fig=vol_fig,
            selection_state=selection_state,
            highlight_source=highlight_source,
            view_toggle=view_toggle,
            identifier_col=identifier_col,
        ),
        code="""
        const selectedLabel = source.data.selected_label;
        const labelXOffset = source.data.label_x_offset;
        const labelAlign = source.data.label_align;
        const xRangeStart = vol_fig.x_range.start;
        const xRangeEnd = vol_fig.x_range.end;
        const xSpan = xRangeEnd - xRangeStart;
        const labelFlipThreshold = xRangeEnd - (xSpan * 0.12);
        let i = 0;
        while (i !== selectedLabel.length) {
            selectedLabel[i] = "";
            labelXOffset[i] = 7;
            labelAlign[i] = "left";
            i += 1;
        }
        source.selected.indices = [];
        if (view_toggle.active) {
            const highlightData = {
                log2FoldChange: [],
                neg_log10_padj: [],
                mean_expression: [],
                gene_name: [],
            };
            selection_state.data.indices[0].forEach(function(index) {
                selectedLabel[index] = identifier_col === "TXNAME"
                    ? source.data.TXNAME[index]
                    : (source.data.gene_name[index] || \
                        source.data[identifier_col][index]);
                if (source.data.log2FoldChange[index] >= labelFlipThreshold) {
                    labelXOffset[index] = -7;
                    labelAlign[index] = "right";
                }
                highlightData.log2FoldChange.push(source.data.log2FoldChange[index]);
                highlightData.neg_log10_padj.push(source.data.neg_log10_padj[index]);
                highlightData.mean_expression.push(source.data.mean_expression[index]);
                highlightData.gene_name.push(source.data.gene_name[index]);
            });
            highlight_source.data = highlightData;
        } else {
            highlight_source.data = {
                log2FoldChange: [],
                neg_log10_padj: [],
                mean_expression: [],
                gene_name: [],
            };
        }
        source.change.emit();
        highlight_source.change.emit();
        """,
    ))
    vol_fig.x_range.js_on_change("start", update_callback)
    vol_fig.x_range.js_on_change("end", update_callback)
    gene_select_toggle.js_on_change("active", CustomJS(
        args=dict(
            source=source,
            selection_state=selection_state,
            gene_select_toggle=gene_select_toggle,
            update_callback=update_callback,
        ),
        code="""
        if (!gene_select_toggle.active) {
            return;
        }
        const selected = selection_state.data.indices[0];
        if (!selected.length) {
            return;
        }
        const geneGroup = source.data.gene_group[selected[0]];
        const groupIndices = [];
        let i = 0;
        while (i !== source.data.gene_group.length) {
            if (source.data.gene_group[i] === geneGroup) {
                groupIndices.push(i);
            }
            i += 1;
        }
        selection_state.data.indices[0] = groupIndices;
        update_callback.execute(source);
        """,
    ))
    plot_mode_toggle.js_on_click(CustomJS(
        args=dict(
            volcano_fig=vol_fig,
            y_axis_slider=y_axis_slider,
            ma_fig=ma_fig,
            plot_mode_toggle=plot_mode_toggle,
        ),
        code="""
        const showingVolcano = volcano_fig.visible;
        volcano_fig.visible = !showingVolcano;
        y_axis_slider.visible = !showingVolcano;
        ma_fig.visible = showingVolcano;
        plot_mode_toggle.label = showingVolcano
            ? "Show volcano plot"
            : "Show MA plot";
        """,
    ))
    selected_source.selected.js_on_change("indices", CustomJS(
        args=dict(
            source=source,
            selected_source=selected_source,
            highlight_source=highlight_source,
            view_toggle=view_toggle,
        ),
        code="""
        if (view_toggle.active) {
            return;
        }
        const rows = selected_source.selected.indices;
        if (rows.length) {
            const row = rows[0];
            const sourceIndex = selected_source.data.source_index[row];
            highlight_source.data = {
                log2FoldChange: [source.data.log2FoldChange[sourceIndex]],
                neg_log10_padj: [source.data.neg_log10_padj[sourceIndex]],
                mean_expression: [source.data.mean_expression[sourceIndex]],
                gene_name: [source.data.gene_name[sourceIndex]],
            };
        } else {
            highlight_source.data = {
                log2FoldChange: [],
                neg_log10_padj: [],
                mean_expression: [],
                gene_name: [],
            };
        }
        highlight_source.change.emit();
        """,
    ))
    selected_source.js_on_change("data", CustomJS(
        args=dict(
            source=source,
            selection_state=selection_state,
            selected_source=selected_source,
            highlight_source=highlight_source,
            gene_select_toggle=gene_select_toggle,
            update_callback=update_callback,
        ),
        code="""
        const listenerKey = "volcano_delete_listener_" + selected_source.id;
        if (!window[listenerKey]) {
            window[listenerKey] = true;
            document.addEventListener("click", function(event) {
                const path = event.composedPath ? event.composedPath() : [];
                const deleteButton = path.find(function(element) {
                    return element.classList
                        && element.classList.contains("volcano-row-delete");
                });
                if (!deleteButton) {
                    return;
                }
                if (deleteButton.dataset.selectedSourceId !== selected_source.id) {
                    return;
                }
                event.preventDefault();
                event.stopPropagation();

                const sourceIndex = Number(deleteButton.dataset.sourceIndex);
                if (gene_select_toggle.active) {
                    selection_state.data.indices[0] = [];
                } else {
                    const kept = [];
                    const selected = selection_state.data.indices[0];
                    let i = 0;
                    while (i !== selected.length) {
                        const index = selected[i];
                        if (index !== sourceIndex) {
                            kept.push(index);
                        }
                        i += 1;
                    }
                    selection_state.data.indices[0] = kept;
                }
                selected_source.selected.indices = [];
                highlight_source.data = {
                    log2FoldChange: [],
                    neg_log10_padj: [],
                    mean_expression: [],
                    gene_name: [],
                };
                update_callback.execute(source);
            });
        }
        """,
    ))

    search_input = AutocompleteInput(
        placeholder=(
            "Search by gene/transcript…" if is_transcript_plot else "Search by gene…"),
        completions=[],
        min_characters=1,
        height=32,
        width=225,
    )
    search_input.js_on_change("value_input", CustomJS(
        args=dict(
            source=source,
            search_input=search_input,
        ),
        code="""
        if (search_input.completions.length) {
            return;
        }
        const data = source.data;
        const completions = new Set();
        const fields = [data.gene_name, data.GENEID];
        if (data.TXNAME) {
            fields.push(data.TXNAME);
        }
        fields.forEach(function(field) {
            let i = 0;
            while (i !== field.length) {
                const value = field[i];
                if (value) {
                    completions.add(String(value));
                }
                i += 1;
            }
        });
        search_input.completions = Array.from(completions).sort();
        """,
    ))
    search_callback = CustomJS(
        args=dict(
            source=source,
            selection_state=selection_state,
            search_input=search_input,
            update_callback=update_callback,
        ),
        code="""
        const query = search_input.value.trim().toLowerCase();
        if (!query) {
            return;
        }
        const data = source.data;
        const currentSelected = new Set(selection_state.data.indices[0]);
        let i = 0;
        while (i !== data.gene_name.length) {
            const fields = [data.gene_name, data.GENEID];
            if (data.TXNAME) {
                fields.push(data.TXNAME);
            }
            let matched = false;
            let fi = 0;
            while (fi !== fields.length) {
                const val = fields[fi][i];
                if (val && val.toLowerCase() === query) {
                    matched = true;
                    break;
                }
                fi += 1;
            }
            if (matched) {
                currentSelected.add(i);
            }
            i += 1;
        }
        selection_state.data.indices[0] = Array.from(currentSelected).sort(
            function(a, b) { return a - b; }
        );
        search_input.value = "";
        update_callback.execute(source);
        """,
    )
    search_input.js_on_change("value", search_callback)

    select_significant_button = Button(
        label="Select significant",
        button_type="primary",
        width=150,
    )
    select_significant_button.js_on_click(CustomJS(
        args=dict(
            source=source,
            selection_state=selection_state,
            fold_slider=fold_slider,
            p_slider=p_slider,
            update_callback=update_callback,
        ),
        code="""
        const data = source.data;
        const fold = fold_slider.value;
        const yThreshold = p_slider.value;
        const passing = [];
        let i = 0;
        while (i !== data.log2FoldChange.length) {
            const outsideFold = Math.abs(data.log2FoldChange[i]) >= fold;
            const aboveP = data.neg_log10_padj[i] >= yThreshold;
            if (outsideFold && aboveP) {
                passing.push(i);
            }
            i += 1;
        }
        selection_state.data.indices[0] = passing;
        update_callback.execute(source);
        """,
    ))

    clear_button = Button(
        label="Clear selection",
        button_type="default",
        height=32,
        width=120,
    )
    clear_button.js_on_click(CustomJS(
        args=dict(
            source=source,
            selection_state=selection_state,
            update_callback=update_callback,
        ),
        code="""
        selection_state.data.indices[0] = [];
        update_callback.execute(source);
        """,
    ))

    export_button = Button(
        label="\u2B07 TSV",
        button_type="light",
        height=32,
        width=55
    )

    export_button.js_on_click(CustomJS(
        args=dict(
            source=source,
            selection_state=selection_state,
            original_columns=original_columns,
            filename=(
                "selected_transcripts.tsv"
                if is_transcript_plot
                else "selected_genes.tsv"
            ),
        ),
        code="""
        const indices = selection_state.data.indices[0];
        if (!indices.length) {
            return;
        }
        const cols = original_columns;
        const data = source.data;
        const rows = [cols.join('\\t')];
        indices.forEach(function(i) {
            rows.push(cols.map(function(col) {
                const v = data[col] !== undefined ? data[col][i] : '';
                return v !== null && v !== undefined ? String(v) : '';
            }).join('\\t'));
        });
        const blob = new Blob([rows.join('\\n') + '\\n'], {
            type: 'text/tab-separated-values',
        });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);
        """,
    ))

    volcano_ma_plot = BokehPlot()
    controls = bokeh_row(
        fold_slider,
        fold_input,
        p_slider,
        p_input,
        y_axis_slider,
        sizing_mode="stretch_width",
    )
    toggle_row = bokeh_row(
        plot_mode_toggle,
        view_toggle,
        gene_select_toggle,
        select_significant_button,
        sizing_mode="stretch_width",
        styles={
            "flex-wrap": "wrap",
            "row-gap": "6px",
        },
    )
    volcano_ma_plot._fig = bokeh_column(
        bokeh_row(controls, sizing_mode="stretch_width"),
        bokeh_row(
            bokeh_column(
                toggle_row,
                vol_fig,
                ma_fig,
                sizing_mode="stretch_width",
            ),
            sizing_mode="stretch_width",
        ),
        sizing_mode="stretch_width",
    )

    classes_table = BokehPlot()
    classes_table._fig = bokeh_column(
        legend_table,
        sizing_mode="stretch_width",
        styles={
            "padding-top": "4px",
        },
    )
    selected_plot = BokehPlot()
    selected_plot._fig = bokeh_column(
        bokeh_row(
            bokeh_row(
                search_input, clear_button, export_button)
        ),
        selected_table,
        sizing_mode="stretch_width",
        styles={
            "padding-top": "8px",
        })
    selected_plot._fig.js_on_event(DocumentReady, CustomJS(
        args=dict(
            selected_source=selected_source,
            selected_table=selected_table,
        ),
        code="""
        const listenerKey = "volcano_selected_table_resize_" + selected_source.id;
        function updateSelectedTableLayout() {
            selected_table.visible = selected_source.data.source_index.length > 0;
            selected_table.height = window.innerWidth <= 1100 ? 180 : 355;
        }
        updateSelectedTableLayout();
        if (!window[listenerKey]) {
            window[listenerKey] = true;
            window.addEventListener("resize", updateSelectedTableLayout);
        }
        """,
    ))
    message = None
    if n_filtered > 0:
        tx_file = (
            "results_dtu_transcript.tsv" if is_transcript_plot else "results_dge.tsv")
        message = (
            f"{n_filtered} features were omitted from the volcano plot because "
            "they do not have plottable log2 fold-change or adjusted p-values. "
            "<br>This is expected for low-information features filtered by the "
            "differential analysis method. Full results are available in "
            f"<i>{tx_file}</i>."
        )
    return volcano_ma_plot, classes_table, selected_plot, message
