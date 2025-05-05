from shiny import App, reactive, render, ui
import numpy as np
from scipy.stats import chisquare
from collections import Counter

# Define segregation models
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

# UI layout
app_ui = ui.page_fluid(
    ui.h2("Genetic Segregation Tester (Phenotype Input)"),
    ui.input_text("phenotypes", "Enter phenotypes (comma-separated)", placeholder="e.g. red, red, orange, red"),
    ui.output_ui("result_ui")
)

# Server logic
def server(input, output, session):

    @reactive.Calc
    def observed_counts():
        try:
            values = [x.strip().lower() for x in input.phenotypes().split(",")]
            count_dict = Counter(values)
            labels = sorted(count_dict.keys())  # consistent ordering
            counts = [count_dict[label] for label in labels]
            return counts, labels
        except Exception:
            return None, None

    @output
    @render.ui
    def result_ui():
        counts, labels = observed_counts()
        if not counts or not labels:
            return ui.p("Please enter valid comma-separated phenotype values (e.g., red, orange, red).")

        result = compare_models(counts)
        if result is None:
            return ui.p(f"No genetic model matches the number of phenotype categories: {len(counts)}.")

        return ui.panel_well(
            ui.h4(f"Best Fitting Model: {result['model']}"),
            ui.p(f"Phenotype Labels: {labels}"),
            ui.p(f"Observed Counts: {result['observed']}"),
            ui.p(f"Expected Counts: {np.round(result['expected'], 2).tolist()}"),
            ui.p(f"Chi-square Statistic: {result['chi_square_stat']:.4f}"),
            ui.p(f"P-value: {result['p_value']:.4f}"),
            ui.p(f"Tested Ratio: {result['best_fit_ratio']}")
        )

# Create app
app = App(app_ui, server)
