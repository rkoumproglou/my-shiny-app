from shiny import App, ui, render, reactive
import matplotlib.pyplot as plt
import matplotlib
import plotly.graph_objs as go
import pandas as pd
import numpy as np
from collections import Counter
import re
from scipy.stats import chisquare

matplotlib.use('Agg')

models = {
    "3:1": [0.75, 0.25],
    "1:2:1": [0.25, 0.5, 0.25],
    "9:3:4": [9/16, 3/16, 4/16],
    "12:3:1": [12/16, 3/16, 1/16],
    "13:3": [13/16, 3/16],
    "9:7": [9/16, 7/16],
    "15:1": [15/16, 1/16],
    "9:6:1": [9/16, 6/16, 1/16]
}

model_descriptions = {
    "3:1": "*3:1 Ratio: This occurs when a single gene has one dominant and one recessive allele. Result: 3 = dominant phenotype, 1 = recessive phenotype.",
    "1:2:1": "*1:2:1 Ratio: This arises when there is incomplete dominance or codominance. Result: 1 = homozygous dominant, 2 = heterozygous, 1 = homozygous recessive.",
    "9:3:4": "*Recessive Epistasis (9:3:4): The recessive allele of one gene masks the expression of the other gene.",
    "12:3:1": "*Dominant Epistasis (12:3:1): The dominant allele of one gene masks the expression of the other gene.",
    "13:3": "*Dominant and Recessive Epistasis (13:3): A dominant inhibitor or a recessive genotype suppresses the expression of the phenotype.",
    "9:7": "*Duplicate Recessive Epistasis (9:7): Both genes must carry at least one dominant allele for the phenotype to be expressed.",
    "15:1": "*Duplicate Dominant Epistasis (15:1): A dominant allele in either gene is sufficient to produce the phenotype.",
    "9:6:1": "*Polymeric Gene Interaction (9:6:1): Both genes contribute additively to the phenotype."
}

def chi_square_test(observed, expected):
    expected_counts = [sum(observed) * p for p in expected]
    stat, p = chisquare(f_obs=observed, f_exp=expected_counts)
    return stat, p

def get_best_model(observed):
    best_model = None
    best_p = -1
    best_stat = None
    for model, expected in models.items():
        if len(expected) == len(observed):
            stat, p = chi_square_test(observed, expected)
            if p > best_p:
                best_model = model
                best_p = p
                best_stat = stat
    return best_model, best_stat, best_p

app_ui = ui.page_fluid(
    ui.panel_title("Genetic Segregation Analyzer"),
    ui.layout_sidebar(
        ui.panel_sidebar(
            ui.input_text_area("observed", "Paste Observed Data (counts or phenotype names)",
                               placeholder="e.g. red\nred\nwhite\nred"),
            ui.input_action_button("analyze", "Analyze")
        ),
        ui.panel_main(
            ui.output_ui("results"),
            ui.output_plot("bar_plot")
        )
    )
)

def server(input, output, session):

    @reactive.event(input.analyze)
    def analyze_data():
        input_str = input.observed().strip()
        normalized_input = re.sub(r'[\t\r]+', '\n', input_str)
        lines = [line.strip() for line in normalized_input.split('\n') if line.strip()]

        if all(re.match(r'^[a-zA-Z0-9 _-]+$', item) for item in lines):
            counter = Counter(lines)
            phenotype_labels = list(counter.keys())
            observed_counts = list(counter.values())
        else:
            normalized_input = re.sub(r'[\n\t]+', ',', input_str)
            normalized_input = re.sub(r'\s*,\s*', ',', normalized_input)
            observed_counts = list(map(int, normalized_input.split(',')))
            phenotype_labels = [f"Class {i+1}" for i in range(len(observed_counts))]

        best_model, best_stat, best_p = get_best_model(observed_counts)

        result_items = []
        if best_p > 0.05:
            result_items.append(ui.h4("Best-Fit Model Report", style="background-color: lightgreen; padding: 0.5em;"))
            result_items.append(ui.card(
                ui.h5(f"The observed ratio does not significantly deviate from the expected ratio based on the Chi-square test (p = {best_p:.3f})."),
                ui.p(f"Best fit model: *{best_model}"),
                ui.p(model_descriptions[best_model]),
                ui.p(f"Chi-square statistic = {best_stat:.3f}, p = {best_p:.3f}"),
                style="box-shadow: 2px 2px 10px rgba(0,0,0,0.2); padding: 1em;"
            ))
        else:
            result_items.append(ui.h4("Best-Fit Model Report", style="background-color: lightgreen; padding: 0.5em;"))
            result_items.append(ui.card(
                ui.h5(f"None of the tested models fit the observed segregation pattern (best p = {best_p:.3f})."),
                ui.p(f"Best attempted model: *{best_model}"),
                ui.p(f"Chi-square statistic = {best_stat:.3f}, p = {best_p:.3f}"),
                style="box-shadow: 2px 2px 10px rgba(0,0,0,0.2); padding: 1em;"
            ))

        return result_items, observed_counts, best_model, best_p, phenotype_labels

    @output
    @render.ui
    def results():
        result, _, _, _, _ = analyze_data()
        return result

    @output
    @render.plot
    def bar_plot():
        _, observed_counts, best_model, best_p, phenotype_labels = analyze_data()
        expected = models[best_model]
        expected_counts = [sum(observed_counts) * p for p in expected]

        x = np.arange(len(observed_counts))
        width = 0.35

        fig, ax = plt.subplots()
        rects1 = ax.bar(x - width/2, observed_counts, width, label='Observed', color='skyblue')
        rects2 = ax.bar(x + width/2, expected_counts, width, label=f'Expected ({best_model}*)' if best_p > 0.05 else 'Expected', color='orange')

        ax.set_ylabel('Counts')
        ax.set_title('Observed vs Expected Segregation')
        ax.set_xticks(x)
        ax.set_xticklabels(phenotype_labels)
        ax.legend()
        fig.patch.set_facecolor('white')
        fig.subplots_adjust(bottom=0.2)
        fig.tight_layout()

        return fig

app = App(app_ui, server)
