# Corrected and complete Shiny for Python script with class count selection

from shiny import App, reactive, render, ui
import numpy as np
from scipy.stats import chisquare
import plotly.graph_objects as go
from collections import Counter

# Define possible Mendelian models
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

# Chi-square test function
def test_segregation(observed, ratios):
    observed = np.array(observed)
    ratios = np.array(ratios)
    total = np.sum(observed)
    expected = (ratios / np.sum(ratios)) * total
    if len(observed) != len(expected):
        raise ValueError("The length of observed and expected counts must match.")
    chi_stat, p_value = chisquare(f_obs=observed, f_exp=expected)
    return {
        'observed': observed.tolist(),
        'expected': expected.tolist(),
        'chi_square_stat': chi_stat,
        'p_value': p_value,
        'best_fit_ratio': ratios.tolist()
    }

# Compare all models matching number of classes
def compare_models(observed_counts):
    results = []
    num_classes = len(observed_counts)
    print("Testing models for", num_classes, "classes")  # Debug

    for name, ratio in models.items():
        if len(ratio) != num_classes:
            continue  # Skip models that don't match the observed class count
        try:
            result = test_segregation(observed_counts, ratio)
            result['model'] = name
            results.append(result)
        except Exception as e:
            print(f"Error testing model {name}: {e}")
            continue

    if not results:
        return None

    # Choose the model with the highest p-value
    best_result = max(results, key=lambda x: x['p_value'])
    return best_result

# UI layout
app_ui = ui.page_fluid(
    ui.h2("Genetic Segregation Ratio Tester"),
    ui.input_select("class_count", "Number of phenotypic classes", {
        "2": "2 classes",
        "3": "3 classes",
        "4": "4 classes"
    }),
    ui.input_text_area("counts", "Enter observed phenotypic values (one per line)", placeholder="Paste one observation per line"),
    ui.output_ui("result_ui"),
    ui.output_ui("plot_ui")
)

# Server logic
def server(input, output, session):

  @reactive.Calc
def observed_counts():
    try:
        raw_data = input.counts().split("\n")
        cleaned_data = [x.strip().lower() for x in raw_data if x.strip()]
        counts = Counter(cleaned_data)
        print("Observed phenotypes:", counts)  # Debug
        return counts
    except Exception as e:
        print("Error parsing input:", e)
        return None

    @output
    @render.ui
    def result_ui():
        counts = observed_counts()
        if not counts:
            return ui.p("Please enter valid data.")

        try:
            sorted_keys = sorted(counts.keys())
            sorted_counts = [counts[key] for key in sorted_keys]
            selected_class_count = int(input.class_count())
            if len(sorted_counts) != selected_class_count:
                return ui.p(f"Mismatch: You selected {selected_class_count} classes, but data contains {len(sorted_counts)} unique classes.")

            result = compare_models(sorted_counts, selected_class_count)
            if result is None:
                return ui.p("No suitable model found.")

            best_fit_msg = "No model explains the observed segregation" if result['p_value'] < 0.05 else f"The best fit model is {result['model']}."
            explanation = "The observed ratio does not significantly deviate from the expected ratio." if result['p_value'] > 0.05 else "The observed ratio significantly deviates from the expected ratio."

            return ui.panel_well(
                ui.h4("Best-Fit Model Report"),
                ui.p(best_fit_msg),
                ui.p(explanation),
                ui.p(f"Observed: {result['observed']}"),
                ui.p(f"Expected: {np.round(result['expected'], 2).tolist()}"),
                ui.p(f"Chi-square Statistic: {result['chi_square_stat']:.4f}"),
                ui.p(f"P-value: {result['p_value']:.4f}"),
                ui.p(f"Tested Ratio: {result['best_fit_ratio']}"),
                ui.p(f"* {result['model']} ratio: {models[result['model']]}") if result['p_value'] > 0.05 else None
            )
        except Exception as e:
            return ui.p(f"Error: {e}")

    @output
    @render.ui
    def plot_ui():
        counts = observed_counts()
        if not counts:
            return ui.p("No plot available, please enter valid data.")

        try:
            sorted_keys = sorted(counts.keys())
            sorted_counts = [counts[key] for key in sorted_keys]
            selected_class_count = int(input.class_count())
            if len(sorted_counts) != selected_class_count:
                return ui.p("Mismatch in class count and observed data.")

            result = compare_models(sorted_counts, selected_class_count)
            if result is None or result['p_value'] < 0.05:
                return ui.p("No plot available for this model as it does not fit well.")

            fig = go.Figure(data=[
                go.Bar(name="Observed", x=sorted_keys, y=result['observed'], marker=dict(color="blue")),
                go.Bar(name="Expected", x=sorted_keys, y=result['expected'], marker=dict(color="red"))
            ])
            fig.update_layout(
                barmode='group',
                title="Observed vs Expected Counts",
                xaxis_title="Phenotypic Class",
                yaxis_title="Count",
                legend_title="Legend"
            )
            return ui.HTML(fig.to_html(full_html=False))
        except Exception as e:
            return ui.p(f"Plot error: {e}")

app = App(app_ui, server)
