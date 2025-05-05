from shiny import App, reactive, render, ui
import numpy as np
from scipy.stats import chisquare

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
    
    result = {
        'observed': observed.tolist(),
        'expected': expected.tolist(),
        'chi_square_stat': chi_stat,
        'p_value': p_value,
        'best_fit_ratio': ratios.tolist()
    }
    return result

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
    ui.h2("Genetic Segregation Ratio Tester"),
    ui.input_text("counts", "Enter observed counts (comma-separated)", placeholder="e.g. 90, 30"),
    ui.output_ui("result_ui")
)

# Server logic
def server(input, output, session):

    @reactive.Calc
    def observed_counts():
        try:
            counts = [int(x.strip()) for x in input.counts().split(",")]
            return counts
        except Exception:
            return None

    @output
    @render.ui
    def result_ui():
        counts = observed_counts()
        if not counts:
            return ui.p("Please enter valid comma-separated integers.")
        
        result = compare_models(counts)
        if result is None:
            return ui.p("No model matches the length of the observed data.")
        
        return ui.panel_well(
            ui.h4(f"Best Fitting Model: {result['model']}"),
            ui.p(f"Observed: {result['observed']}"),
            ui.p(f"Expected: {np.round(result['expected'], 2).tolist()}"),
            ui.p(f"Chi-square Statistic: {result['chi_square_stat']:.4f}"),
            ui.p(f"P-value: {result['p_value']:.4f}"),
            ui.p(f"Tested Ratio: {result['best_fit_ratio']}")
        )

# Create app
app = App(app_ui, server)