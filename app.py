from shiny import App, reactive, render, ui
import numpy as np
from scipy.stats import chisquare
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
    ui.input_text_area("counts", "Enter observed counts (paste column from Excel, one per line)", placeholder="Paste one observation per line"),
    ui.output_ui("result_ui"),
    ui.output_ui("plot_ui")  # Change this to a generic output_ui for Plotly rendering
)

# Server logic
def server(input, output, session):
    
    @reactive.Calc
    def observed_counts():
        try:
            # Split the input by newlines and convert to integers
            counts = [int(x.strip()) for x in input.counts().split("\n") if x.strip()]
            return counts
        except Exception:
            return None

    @output
    @render.ui
    def result_ui():
        counts = observed_counts()
        if not counts:
            return ui.p("Please enter valid data.")
        
        result = compare_models(counts)
        if result is None:
            return ui.p("No model matches the length of the observed data.")
        
        # Best-fit model message
        best_fit_msg = "No model explains the observed segregation" if result['p_value'] < 0.05 else f"The best fit model is {result['model']}."

        # Explanation based on p-value
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
            ui.p(f"* {result['model']} ratio: {models[result['model']]}")
        )
    
    @output
    @render.ui
    def plot_ui():
        counts = observed_counts()
        if not counts:
            return ui.p("No plot available, please enter valid data.")
        
        result = compare_models(counts)
        if result is None or result['p_value'] < 0.05:
            return ui.p("No plot available for this model as it does not fit well.")
        
        # Create a bar chart
        fig = go.Figure(data=[
            go.Bar(name="Observed", x=list(range(len(counts))), y=result['observed'], marker=dict(color="blue")),
            go.Bar(name="Expected", x=list(range(len(result['expected']))), y=result['expected'], marker=dict(color="red"))
        ])
        
        fig.update_layout(
            barmode='group',  # Bars within each category won't touch each other
            title="Observed vs Expected Counts",
            xaxis_title="Phenotypic Class",
            yaxis_title="Count",
            legend_title="Model",
            legend=dict(
                itemsizing='constant',
                traceorder='normal',
                font=dict(size=12)
            )
        )

        # Embed Plotly graph into HTML div
        plot_html = fig.to_html(full_html=False)
        return ui.HTML(plot_html)  # Use ui.HTML to embed the plotly graph as HTML

# Create app
app = App(app_ui, server)
