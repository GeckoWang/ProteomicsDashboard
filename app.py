import itertools
import dash
from dash import dcc, html, Input, Output
import plotly.express as px
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# load data
df = pd.read_csv("data/TPD_genename6686_20230618.csv")
# replace index to actual name
mapping_h = {"N": "Normal", "M": "Multinodular goiter", "A": "Follicular adenoma",
             "C": "Follicular thyroid carcinoma", "P": "Papillary thyroid carcinoma"}
mapping_c = {"M": "Malignant", "B": "Benign"}
df["Histopathology_type"] = df["Histopathology_type"].replace(mapping_h)
df["Classification_type"] = df["Classification_type"].replace(mapping_c)
protein_name = df.columns[1:]


# median of the whole type
df_r2 = pd.read_csv("data/df_r2.csv")
bar_fig = px.bar(
    data_frame=df_r2,
    x="rank",
    y="value",
    color="Location",
    category_orders={"Location": ["Extracellular Space", "Cytoplasm", "Nucleus", "Plasma Membrane", "Other", "Not Available"
   ]},
    hover_name="index", hover_data=["Biomarker Application(s)"],
)


bar_fig.update_layout(
    bargap=0,
    # barmode="overlay",
    xaxis=dict(
        tickmode="linear",  # 线性刻度
        tick0=0,            # 从 0 开始
        dtick=1000),
    yaxis=dict(title_text="Log2(protein abundances)", range=[10, None]),
    plot_bgcolor="white",
)

# median of the select type
df_r3raw = pd.read_csv("data/df_r3.csv")


app = dash.Dash(__name__)


def make_break(num_breaks):
    br_list = [html.Br()] * num_breaks
    return br_list


app.layout = html.Div(children=[
    html.Div(children=[
        html.Div(children=
                 [html.H1("Thyroid website")]),
        html.Div(children=
                 [html.H1("Protein")]),
        html.Div(children=[
            html.Div(children=[
                html.H4("Please select your protein of interest"),
                dcc.Dropdown(
                    id="dropdown",
                    options=[
                        {"label":protein, "value":protein} for protein in protein_name
                    ],
                    placeholder="Please choose a protein",
                ),
                *make_break(2),
                html.H4("Label type"),
                dcc.RadioItems(
                    id="label_type",
                    options=[
                        {"label": "Histopathology", "value": "Histopathology_type"},
                        {"label": "Classification", "value": "Classification_type"},
                    ],
                ),
                *make_break(2),
                html.H4("Display jitter?"),
                dcc.RadioItems(
                    id="jitter",
                    options=[
                        {"label": "Yes", "value": "True"},
                        {"label": "No", "value": "False"},
                    ],
                ),
                *make_break(2),
                html.H4("Conduct t test for each group?"),
                dcc.RadioItems(
                    id="t_test",
                    options=[
                        {"label": "Raw t-test", "value": "t_test"},
                        {"label": "B-H adjusted t-test", "value": "bh_t_test"},
                        {"label": "No t-test", "value": "False"}
                    ],
                ),
                *make_break(2),
                html.H4("Conduct ANOVA for each group?"),
                dcc.RadioItems(
                    id="anova",
                    options=[
                        {"label": "Yes", "value": "True"},
                        {"label": "No", "value": "False"},
                    ],
                ),
            ], style={"width": "30%", "backgroundColor": "#f2f2f2", "padding": "20px"}),
            dcc.Graph(id="proteomics_graph", style={"width": "70%", "height": "auto", "max-height": "700px"})
        ], style={"display": "flex", "margin": "25px 0px"}),
    ]),
    html.Div(children=[
        html.Div(children=[
            html.H4("Select your protein of interest"),
            dcc.Dropdown(
                id="drop_rp",
                options=[
                    {"label":protein, "value":protein} for protein in protein_name
                ],
                placeholder="Please choose a protein",
            ),
            ], style={"width": "30%", "backgroundColor": "#f2f2f2", "padding": "20px"}),
        dcc.Graph(id="rank_graph2", figure=bar_fig, style={"width": "70%", "height": "400px"})
        ], style={"display": "flex", "margin": "25px 0px"}),
    html.Div(children=[
        html.Div(children=[
            html.H4("Protein expression in different tissue types"),
            dcc.Dropdown(
                id="drop_rt",
                options=[
                    {"label": "Normal", "value": "Normal"},
                    {"label": "Multinodular goiter", "value":"Multinodular goiter"},
                    {"label": "Follicular adenoma", "value": "Follicular adenoma"},
                    {"label": "Follicular thyroid carcinoma", "value": "Follicular thyroid carcinoma"},
                    {"label": "Papillary thyroid carcinoma", "value": "Papillary thyroid carcinoma"},
                    {"label": "Malignant", "value": "Malignant"},
                    {"label": "Benign", "value": "Benign"}
                    ],
                placeholder="Please choose a histopathology/classification type",
            ),
        ], style={"width": "30%", "backgroundColor": "#f2f2f2", "padding": "20px"}),
        dcc.Graph(id="rank_graph3", style={"width": "70%", "height": "400px"})
    ], style={"display": "flex", "margin": "25px 0px"}),

])


@app.callback(
    # Output("text_output", "children"),
    Output("proteomics_graph", "figure"),
    Input("dropdown","value"),
    Input("label_type", "value"),
    Input("t_test", "value"),
    Input("jitter", "value"),
    # Input("t_step", "value"),
    Input("anova", "value"),
)
def update_label_type(select_protein, select_label, select_t_test, select_jitter, select_anova):
    df_2 = df.copy(deep=True)

    if select_protein and select_label:
        select_data = df_2[[select_protein, select_label]]
        valid_data = select_data.dropna(subset=[select_protein]).groupby(select_label)[select_protein]
        valid1_groups = [valid1_group.tolist() for _, valid1_group in valid_data if len(valid1_group) > 1]
        sample_counts = valid_data.count()
        label_type_name = select_data[select_label].unique()
        # protein_data = select_data.groupby(select_label).sum().reset_index()
        # fig = px.box(protein_data, x=select_label, color=select_label, points="all")
        # protein_data = select_data.groupby(select_label, as_index=False).apply(lambda group: group)
        # protein_data = protein_data.reset_index(drop=True)
        annotations = []
        for i, (category, n_count) in enumerate(sample_counts.items()):
            n_label = label_type_name.tolist().index(category)
            if n_count > 0:
                annotations.append(
                    dict(
                        x=n_label, y=-0.02,
                        text=f"n={n_count}",
                        showarrow=False,
                        font=dict(size=10, color="black"),
                        yanchor="top", yref="paper"
                    )
                )
        fig = px.box(select_data, x=select_label, y=select_protein, color=select_label, points="all")
        fig.update_layout(xaxis=dict(tickangle=0, tickfont=dict(size=10), showline=True, linewidth=1, linecolor='black')
                          , yaxis=dict(showline=True, linewidth=1, linecolor="black",
                                       title_text="Log2(protein abundances)"),
                          plot_bgcolor="white",
                          annotations=annotations, margin=dict(b=80))

        if select_t_test in ["t_test", "bh_t_test"]:
            # t_test combinations

            label_type_combinations = list(itertools.combinations(label_type_name, 2))
            # adjust y_axis for t_test
            y_max = select_data[select_protein].max()
            # step increase of t test label
            y_offset = 1

            p_values = []
            pair_labels = []
            for group1, group2 in label_type_combinations:
                data1 = select_data[select_data[select_label] == group1][select_protein].dropna()
                data2 = select_data[select_data[select_label] == group2][select_protein].dropna()

                # skip cases with insufficient samples
                if data1.empty or data2.empty or len(data1) <= 1 or len(data2) <= 1:
                    continue

                t_stat, p_value = stats.ttest_ind(data1, data2)
                p_values.append(p_value)
                pair_labels.append((group1, group2))

            # adjust different type of t_test
            if select_t_test == "bh_t_test":
                # for B-H-adjusted t-test
                _, adj_p_values, _, _ = multipletests(p_values, method="fdr_bh")
            else:
                adj_p_values = p_values

            for i, ((group1, group2), p_value) in enumerate(zip(pair_labels, adj_p_values)):
                p_value_f = f"{p_value:.4e}"
                # significance signs, no longer use
                # significance = get_significance(p_value)

                # t_test place
                y_offset = y_offset
                y_pos = y_max + (i + 1) * y_offset
                # add t_test line
                fig.add_shape(
                    type="line",
                    x0=group1, x1=group2, y0=y_pos, y1=y_pos,
                    line=dict(color="black", width=1)
                )
                # add t_test significance
                fig.add_annotation(
                    x=(label_type_name.tolist().index(group1) + label_type_name.tolist().index(group2)) / 2,
                    y=y_pos + 0.6,
                    text=p_value_f,
                    showarrow=False,
                    font=dict(size=12)
                )

        if select_anova == "True":
            _, anova_value = stats.f_oneway(*valid1_groups)
            if len(valid1_groups) > 1:
                fig.add_annotation(
                    x=0,
                    y=1,
                    yref="paper",
                    text=f"ANOVA p-value: {anova_value:.4e}",
                    showarrow=False
                )

        if select_jitter == "True":
            fig.update_traces(marker=dict(size=2))

        return fig

    else:
        return {}


@app.callback(
    Output("rank_graph2", "figure"),
    Input("drop_rp", "value"),
)
def update_rank_annotation(select_rp):
    fig2 = bar_fig
    fig2.layout.annotations = []
    if select_rp:
        protein_info = df_r2[df_r2["index"] == select_rp]

        x_value = protein_info["rank"].values[0]
        y_value = protein_info["value"].values[0]
        annotation_text = f"{select_rp}<br>{y_value:.2f}"

        fig2.add_annotation(
            x=x_value, y=y_value,
            text=annotation_text,
            showarrow=True, arrowhead=2,
            font=dict(size=12, color="black"),
            bgcolor="lightgrey",
            # bordercolor="black"
        )

    return fig2


@app.callback(
    Output("rank_graph3", "figure"),
    Input("drop_rt", "value")
)
def update_rank_type(select_rt):
    df_r3 = df_r3raw.copy(deep=True)
    if select_rt:

        df_r3 = df_r3.dropna(subset=[select_rt])
        df_r3_sort = df_r3.sort_values(by=select_rt, ascending=False)

        # df_r3_sort["Protein"] = df_r3_sort.index
        df_r3_sort.reset_index(inplace=True)
        df_r3_sort.rename(columns={"Unnamed: 0": "Protein"}, inplace=True)
        df_r3_sort["rank"] = range(1, df_r3_sort[select_rt].count() + 1)

        fig3 = px.bar(
            data_frame=df_r3_sort,
            x="rank",
            y=select_rt,
            color="Location",
            category_orders={"Location": ["Extracellular Space", "Cytoplasm", "Nucleus", "Plasma Membrane", "Other", "Not Available"
                             ]},
            hover_data=["Protein", "Biomarker Application(s)"],
        )

        fig3.update_layout(
            bargap=0,
            # barmode="overlay",
            xaxis=dict(
                tickmode="linear",  # 线性刻度
                tick0=0,  # 从 0 开始
                dtick=1000),
            yaxis=dict(title_text="Log2(protein abundances)", range=[10, None]),
            plot_bgcolor="white",
        )

        return fig3
    else:
        return {}


if __name__ == "__main__":
    app.run_server(debug=True)

