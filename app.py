from shiny import App, reactive, render, ui
import numpy as np
from scipy.stats import chisquare
from collections import Counter
import re
import pandas as pd
from plotnine import ggplot, aes, geom_bar, labs, position_dodge, theme_minimal, scale_fill_manual, ggtitle

# Segregation models
models = {
    "3:1": [3, 1],
    "1:2:1": [1, 2, 1],
    "1:1:1:1": [1, 1, 1, 1],
    "9:3:3:1": [9, 3, 3, 1],
    "12:3:1": [12, 3, 1],
    "9:7": [9, 7],
    "15:1": [15, 1],
    "9:3:4": [9, 3, 4],
    "13:3": [13, 3],
    "9:6:1": [9, 6, 1]
}

# Model explanations
model_explanations = {
    "3:1": "Occurs when a single gene has one dominant and one recessive allele. F₂ ratio: 3 dominant : 1 recessive phenotype.",
    "1:2:1": "Arises from incomplete dominance or codominance. F₂ ratio: 1 homozygous dominant : 2 heterozygous : 1 homozygous recessive phenotype.",
    "1:1:1:1": "Typical of dihybrid test crosses with independent assortment of two genes.",
    "9:3:3:1": "Expected from a dihybrid F₂ with independent assortment of two genes each with dominant and recessive alleles.",
    "12:3:1": "Dominant epistasis: a dominant allele of one gene masks expression of another gene.",
    "9:7": "Duplicate recessive epistasis: both genes must carry a dominant allele for the trait to appear.",
    "15:1": "Duplicate dominant epistasis: a dominant allele at either gene is enough to produce the trait.",
    "9:3:4": "Recessive epistasis: a recessive allele at one gene masks expression of another gene.",
    "13:3": "Dominant and recessive epistasis: suppression occurs when an inhibitor is present or the other gene is homozygous recessive.",
    "9:6:1": "Polymeric gene interaction: both genes additively contribute to phenotype; both dominant alleles produce enhanced trait."
}

def test_segregation(observed, ratios):
    observed = np.array(observed)
    ratios = np.array(ratios)
    total = np.sum(observed)
    expected = ratios / np.sum(ratios) * total
    chi_stat, p_value = chisquare(f_obs=observed, f_exp=expected)
    return {
        'observed': observed.tolist(),
        'expected': expected.tolist(),
        'chi_square_stat': chi_stat,
        'p_value': p_value,
        'best_fit_ratio': ratios.tolist()
    }

def compare_models(observed_counts):
    results = []
    for name, ratio in models.items():
        if len(ratio) == len(observed_counts):
            result = test_segregation(observed_counts, ratio)
            result['model'] = name
            results.append(result)
    if not results:
        return None
    best_result = max(results, key=lambda x: x['p_value'])
    return best_result

# UI
app_ui = ui.page_fluid(
    ui.h2("Genetic Segregation Tester"),
    ui.input_text_area("phenotypes", "Enter phenotypes (comma or newline-separated)",
                       placeholder="e.g. red, red, orange\nor\nred\nred\norange", rows=5),
    ui.output_ui("result_ui"),
    ui.output_plot("result_plot")
)

# Server
def server(input, output, session):

    @reactive.Calc
    def observed_data():
        try:
            raw_input = input.phenotypes().strip()
            values = re.split(r"[\n,]+", raw_input)
            values = [v.strip().lower() for v in values if v.strip()]
            count_dict = Counter(values)
            labels = sorted(count_dict.keys())
            counts = [count_dict[label] for label in labels]
            return counts, labels
        except Exception:
            return None, None

    @reactive.Calc
    def model_result():
        counts, labels = observed_data()
        if not counts:
            return None
        result = compare_models(counts)
        if result:
            result["labels"] = labels
        return result

    @output
    @render.ui
    def result_ui():
        result = model_result()
        if result is None:
            return ui.p("Please enter valid phenotype data.")

        explanation = model_explanations.get(result['model'], "Explanation not available.")
        return ui.panel_well(
            ui.h4(f"Best Fitting Model: {result['model']}"),
            ui.p(f"Phenotype Labels: {result['labels']}"),
            ui.p(f"Observed Counts: {result['observed']}"),
            ui.p(f"Expected Counts: {np.round(result['expected'], 2).tolist()}"),
            ui.p(f"Chi-square Statistic: {result['chi_square_stat']:.4f}"),
            ui.p(f"P-value: {result['p_value']:.4f}"),
            ui.p(f"Tested Ratio: {result['best_fit_ratio']}"),
            ui.hr(),
            ui.p(ui.strong("Model Explanation:")),
            ui.p(explanation)
        )

    @output
    @render.plot
    def result_plot():
        result = model_result()
        if result is None:
            return None

        labels = result["labels"]
        df = pd.DataFrame({
            "Phenotype": labels * 2,
            "Count": result["observed"] + result["expected"],
            "Type": ["Observed"] * len(labels) + ["Expected"] * len(labels)
        })

        plot = (
            ggplot(df, aes(x="Phenotype", y="Count", fill="Type")) +
            geom_bar(stat="identity", position=position_dodge(width=0.7)) +
            scale_fill_manual(values={"Observed": "#FF9999", "Expected": "#9999FF"}) +
            labs(y="Count", fill="Legend") +
            ggtitle(f"Observed vs Expected Counts ({result['model']})") +
            theme_minimal()
        )
        return plot

# Run app
app = App(app_ui, server)
