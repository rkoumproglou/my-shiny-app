from shiny import App, ui, render, reactive
import plotly.graph_objs as go
import numpy as np
from scipy.stats import chisquare
import re
from collections import Counter

models = {
    "3:1": [3, 1],
    "1:2:1": [1, 2, 1],
    "9:3:4": [9, 3, 4],
    "12:3:1": [12, 3, 1],
    "13:3": [13, 3],
    "9:7": [9, 7],
    "15:1": [15, 1],
    "9:6:1": [9, 6, 1]
}

model_explanations = {
    "3:1": "This occurs when a single gene has one dominant and one recessive allele.",
    "1:2:1": "This ratio arises when there is incomplete dominance or codominance.",
    "9:3:4": "The recessive allele of one gene masks the expression of the other gene (recessive epistasis).",
    "12:3:1": "The dominant allele of one gene masks the expression of the other gene (dominant epistasis).",
    "13:3": "A dominant inhibitor or homozygous recessive at the other gene suppresses phenotype expression.",
    "9:7": "Both genes must carry at least one dominant allele for the phenotype to be expressed (duplicate recessive epistasis).",
    "15:1": "A dominant allele in either gene is sufficient to produce the phenotype (duplicate dominant epistasis).",
    "9:6:1": "Both genes contribute additively to the phenotype, producing enhanced, intermediate, or basic traits."
}

app_ui = ui.page_fluid(
    ui.panel_title("Chi-square Genetic Model Fit Tool"),
    ui.input_text_area("observed", "Paste your data (phenotypes or counts):", placeholder="e.g. red"),
    ui.output_ui("result"),
    ui.output_plot("plot")
)

def server(input, output, session):
    @reactive.calc
    def parsed_observed():
        input_str = input.observed().strip()
        if not input_str:
            return None, None

        # Normalize line breaks and tabs
        normalized_input = re.sub(r'[\t\r]+', '\n', input_str)
        lines = [line.strip() for line in normalized_input.split('\n') if line.strip()]

        # Determine if it's individual phenotypes or counts
        if all(re.match(r'^[a-zA-Z0-9 _-]+$', item) for item in lines):
            counter = Counter(lines)
            phenotype_labels = list(counter.keys())
            observed_counts = list(counter.values())
        else:
            normalized_input = re.sub(r'[\n\t]+', ',', input_str)
            normalized_input = re.sub(r'\s*,\s*', ',', normalized_input)
            try:
                observed_counts = list(map(int, normalized_input.split(',')))
                phenotype_labels = [f"Class {i+1}" for i in range(len(observed_counts))]
            except ValueError:
                return None, None

        return observed_counts, phenotype_labels

    @reactive.calc
    def best_model():
        observed, _ = parsed_observed()
        if not observed:
            return None

        total = sum(observed)
        min_chi2 = float('inf')
        best_fit = None
        best_p = None
        best_expected = []

        for model, ratio in models.items():
            if len(ratio) != len(observed):
                continue
            expected = np.array(ratio) / sum(ratio) * total
            chi2, p = chisquare(f_obs=observed, f_exp=expected)
            if chi2 < min_chi2:
                min_chi2 = chi2
                best_fit = model
                best_p = p
                best_expected = expected

        return best_fit, best_p, best_expected

    @output
    @render.ui
    def result():
        observed, labels = parsed_observed()
        if not observed:
            return ui.div("Please enter valid phenotype labels or counts.")

        best_fit, p, expected = best_model()
        p_cutoff = 0.05

        interpretation = (
            f"The observed ratio does not significantly deviate from the expected ratio based on the Chi-square test, "
            f"so it is considered consistent with the best-fit model (*{best_fit}*)."
            if p > p_cutoff else
            "The observed ratio significantly deviates from all expected ratios. No model fits the data well."
        )

        components = [
            ui.div(
                ui.h4("Best-Fit Model Report"),
                class_="p-3 mb-3",
                style="background-color: #d4edda; border-radius: 6px;"
            ),
            ui.card(
                ui.markdown(f"**Interpretation:** {interpretation}"),
                style="box-shadow: 0 4px 8px rgba(0,0,0,0.2); padding: 15px;"
            )
        ]

        if p > p_cutoff:
            components += [
                ui.card(
                    ui.markdown(f"**Best-fit model:** *{best_fit}*"),
                    ui.markdown(f"**Expected ratio:** {models[best_fit]}"),
                    ui.markdown(f"**Chi-square:** {round(chisquare(f_obs=observed, f_exp=expected)[0], 2)}"),
                    ui.markdown(f"**p-value:** {round(p, 4)}"),
                    style="box-shadow: 0 4px 8px rgba(0,0,0,0.2); padding: 15px;"
                ),
                ui.card(
                    ui.markdown(f"**Explanation for *{best_fit}* model:** {model_explanations[best_fit]}"),
                    style="box-shadow: 0 4px 8px rgba(0,0,0,0.2); padding: 15px;"
                )
            ]
        else:
            components += [
                ui.card(
                    ui.markdown(f"**Best-fit model:** *{best_fit}*"),
                    ui.markdown(f"**Expected ratio:** {models[best_fit]}"),
                    ui.markdown(f"**Chi-square:** {round(chisquare(f_obs=observed, f_exp=expected)[0], 2)}"),
                    ui.markdown(f"**p-value:** {round(p, 4)}"),
                    style="box-shadow: 0 4px 8px rgba(0,0,0,0.2); padding: 15px;"
                )
            ]

        return ui.div(*components)

    @output
    @render.plot
    def plot():
        observed, labels = parsed_observed()
        if not observed:
            return None

        best_fit, p, expected = best_model()
        p_cutoff = 0.05

        trace1 = go.Bar(
            x=labels,
            y=observed,
            name='Observed',
            marker_color='steelblue',
            offsetgroup=0
        )
        trace2 = go.Bar(
            x=labels,
            y=expected,
            name='Expected' + (f" ({best_fit})" if p > p_cutoff else ""),
            marker_color='orange',
            offsetgroup=1
        )

        layout = go.Layout(
            title='Observed vs Expected Counts',
            barmode='group',
            bargap=0.2,
            showlegend=True,
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
        )
        fig = go.Figure(data=[trace1, trace2], layout=layout)
        fig.update_layout(
            legend=dict(bgcolor="white", bordercolor="black", borderwidth=1),
        )

        return fig

app = App(app_ui, server)
