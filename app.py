from shiny import App, ui, render, reactive
from plotly import graph_objs as go
from scipy.stats import chisquare
import math

# Define models and expected ratios
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

model_definitions = {
    "3:1": "This occurs when a single gene has one dominant and one recessive allele.",
    "1:2:1": "This ratio arises when there is incomplete dominance or codominance.",
    "9:3:4": "The recessive allele of one gene masks the expression of the other gene.",
    "12:3:1": "The dominant allele of one gene masks the expression of the other gene.",
    "13:3": "A dominant inhibitor or a homozygous recessive condition suppresses the phenotype.",
    "9:7": "Both genes must carry at least one dominant allele for the phenotype to be expressed.",
    "15:1": "A dominant allele in either gene is sufficient to produce the phenotype.",
    "9:6:1": "Both genes contribute additively to the phenotype."
}

def parse_input(text_input: str) -> list[int]:
    # Replace commas with newlines just in case both are used
    normalized = text_input.replace(",", "\n")
    # Split by newlines and remove empty lines/spaces
    lines = [line.strip() for line in normalized.strip().splitlines() if line.strip()]
    try:
        counts = [int(value) for value in lines]
        return counts
    except ValueError:
        return []

app_ui = ui.page_fluid(
    ui.panel_title("Genetic Segregation Model Fit"),
    ui.input_text_area("observed", "Paste Observed Counts (comma-separated or Excel column)", rows=5, placeholder="e.g. 100, 35, 15"),
    ui.output_ui("results")
)

def server(input, output, session):
    @reactive.Calc
    def best_model():
        observed = parse_input(input.observed())
        if not observed:
            return None

        total = sum(observed)
        best_p = -1
        best_fit = None
        best_expected = []
        stats_summary = ""

        for name, ratio in models.items():
            factor = total / sum(ratio)
            expected = [r * factor for r in ratio]
            if len(expected) != len(observed):
                continue
            stat, p = chisquare(f_obs=observed, f_exp=expected)
            if p > best_p:
                best_p = p
                best_fit = name
                best_expected = expected
                stats_summary = f"Chi-square: {stat:.2f}, p-value: {p:.4f}"

        return {
            "model": best_fit,
            "p": best_p,
            "expected": best_expected,
            "observed": observed,
            "stats": stats_summary
        }

    @output
    @render.ui
    def results():
        result = best_model()
        if not result:
            return ui.markdown("**Please enter valid observed counts.**")

        interpretation = ""
        if result["p"] > 0.05:
            interpretation = (
                f"<div style='background-color: #d4edda; padding: 10px; border-radius: 5px;'><strong>Best-Fit Model Report</strong></div>"
                f"<div class='card shadow p-3 mb-3'>"
                f"<h5><strong>*{result['model']}*</strong> is the best fit model.</h5>"
                f"<p>The observed ratio does not significantly deviate from the expected ratio based on the Chi-square test, so it is considered consistent with the model.</p>"
                f"<p>{model_definitions[result['model']]}</p>"
                f"<p><strong>Statistics:</strong> {result['stats']}</p>"
                f"</div>"
            )
        else:
            interpretation = (
                f"<div class='card shadow p-3 mb-3'>"
                f"<h5><strong>No model adequately explains the observed segregation ratio.</strong></h5>"
                f"<p>Although the best fit was <strong>*{result['model']}*</strong>, the observed values significantly deviate from the expected values (p = {result['p']:.4f}).</p>"
                f"<p><strong>Statistics:</strong> {result['stats']}</p>"
                f"</div>"
            )

        # Build the graph
        data = []
        categories = [f"Class {i+1}" for i in range(len(result["observed"]))]
        bar_width = 0.35
        x = list(range(len(categories)))

        data.append(go.Bar(
            x=[i - bar_width/2 for i in x],
            y=result["observed"],
            name="Observed",
            marker=dict(color="rgba(55, 83, 109, 0.7)"),
            width=bar_width
        ))

        data.append(go.Bar(
            x=[i + bar_width/2 for i in x],
            y=result["expected"],
            name="Expected" + (f" ({result['model']})" if result["p"] > 0.05 else ""),
            marker=dict(color="rgba(26, 118, 255, 0.7)"),
            width=bar_width
        ))

        layout = go.Layout(
            title="Observed vs Expected Counts",
            xaxis=dict(title="Phenotypic Classes", tickvals=x, ticktext=categories),
            yaxis=dict(title="Count"),
            barmode='overlay',
            bargap=0.2,
            showlegend=True,
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            margin=dict(l=40, r=40, t=40, b=40),
            legend=dict(bgcolor='rgba(255,255,255,0.5)', bordercolor='rgba(0,0,0,0.1)', borderwidth=1),
        )

        fig = go.Figure(data=data, layout=layout)
        fig.update_layout(
            shapes=[dict(type='rect', xref='paper', yref='paper',
                         x0=0, y0=0, x1=1, y1=1,
                         line=dict(width=0), fillcolor='rgba(255,255,255,0)', layer='below')]
        )

        return ui.div(
            ui.HTML(interpretation),
            ui.output_plot("plot")
        )

    @output
    @render.plot
    def plot():
        result = best_model()
        if not result:
            return

        categories = [f"Class {i+1}" for i in range(len(result["observed"]))]
        bar_width = 0.35
        x = list(range(len(categories)))

        fig = go.Figure()
        fig.add_bar(
            x=[i - bar_width/2 for i in x],
            y=result["observed"],
            name="Observed",
            marker_color="rgba(55, 83, 109, 0.7)",
            width=bar_width
        )
        fig.add_bar(
            x=[i + bar_width/2 for i in x],
            y=result["expected"],
            name="Expected" + (f" ({result['model']})" if result["p"] > 0.05 else ""),
            marker_color="rgba(26, 118, 255, 0.7)",
            width=bar_width
        )
        fig.update_layout(
            title="Observed vs Expected Counts",
            xaxis=dict(title="Phenotypic Classes", tickvals=x, ticktext=categories),
            yaxis=dict(title="Count"),
            barmode='overlay',
            bargap=0.2,
            showlegend=True,
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            margin=dict(l=40, r=40, t=40, b=40),
            legend=dict(bgcolor='rgba(255,255,255,0.5)', bordercolor='rgba(0,0,0,0.1)', borderwidth=1),
        )
        return fig

app = App(app_ui, server)
