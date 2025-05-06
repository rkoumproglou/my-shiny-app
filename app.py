from shiny import App, reactive, render, ui
import numpy as np
import plotly.graph_objects as go
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

# UI layout with CSS for drop shadow
app_ui = ui.page_fluid(
    ui.h2("Genetic Segregation Ratio Tester"),
    ui.input_textarea("counts", "Paste phenotypic values (one per line)", placeholder="e.g.\nred\nred\nred\norange\norange"),
    ui.output_ui("result_ui"),
    ui.tags.style("""
        .graph-container {
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
        }
    """)
)

# Server logic
def server(input, output, session):

    @reactive.Calc
    def observed_counts():
        try:
            # Process the input as a list of phenotypic values
            phenotypes = [x.strip() for x in input.counts().splitlines() if x.strip()]
            # Count the occurrences of each phenotype
            unique_phenotypes = set(phenotypes)
            counts = [phenotypes.count(phenotype) for phenotype in unique_phenotypes]
            return counts, list(unique_phenotypes)
        except Exception:
            return None, None

    @output
    @render.ui
    def result_ui():
        counts, unique_phenotypes = observed_counts()
        if not counts:
            return ui.p("Please enter valid phenotypic values (one per line).")
        
        result = compare_models(counts)
        if result is None:
            return ui.p("No model matches the length of the observed data.")
        
        # Determine p-value interpretation
        if result['p_value'] < 0.05:
            interpretation = f"The observed ratio does not significantly deviate from the expected ratio based on the Chi-square test, so it is considered consistent with the best fit model {result['model']}."
        else:
            interpretation = "No model explains the observed segregation ratio based on the Chi-square test."
        
        # Construct the output
        output = ui.panel_well(
            ui.h4("Best-Fit Model Report", style="background-color: lightgreen; padding: 10px;"),
            ui.p(interpretation),
            ui.p(f"Best Fit Model: {result['model']}"),
            ui.p(f"Chi-square Statistic: {result['chi_square_stat']:.4f}"),
            ui.p(f"P-value: {result['p_value']:.4f}"),
            ui.p(f"Tested Ratio: {result['best_fit_ratio']}")
        )

        # Add graph if p-value is greater than 0.05
        if result['p_value'] >= 0.05:
            output += ui.card(
                ui.plotly({
                    "data": [
                        go.Bar(
                            x=unique_phenotypes,
                            y=result['observed'],
                            name='Observed',
                            marker=dict(color="rgb(204,204,255)")
                        ),
                        go.Bar(
                            x=unique_phenotypes,
                            y=result['expected'],
                            name='Expected',
                            marker=dict(color="rgb(255,204,204)")
                        ),
                    ],
                    "layout": {
                        "title": "Observed vs. Expected Distribution",
                        "showlegend": True,
                        "xaxis": {"title": "Phenotypes"},
                        "yaxis": {"title": "Counts"},
                        "bargap": 0.2  # To prevent the bars from touching
                    }
                }),
                style="box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);"
            )
        
        return output

# Create app
app = App(app_ui, server)
