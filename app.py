from shiny import App, ui, render, reactive
import pandas as pd
from collections import Counter
from scipy.stats import chisquare

# Define Mendelian models
MENDELIAN_MODELS = {
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

# Define UI
app_ui = ui.page_fluid(
    ui.h2("Mendelian Ratio Chi-square Tester"),
    ui.markdown("**Paste your phenotypic data (one value per line):**"),
    ui.input_text_area("phenodata", "", rows=15),
    ui.hr(),
    ui.output_ui("result_ui")
)

# Define server logic
def server(input, output, session):

    @reactive.calc
    def parsed_data():
        lines = input.phenodata().strip().splitlines()
        cleaned = [line.strip().lower() for line in lines]
        return Counter(cleaned)

    @reactive.calc
    def matched_models():
        count = len(parsed_data())
        return {
            name: ratio
            for name, ratio in MENDELIAN_MODELS.items()
            if len(ratio) == count
        }

    @reactive.calc
    def test_results():
        obs_counts = parsed_data()
        total = sum(obs_counts.values())
        sorted_categories = sorted(obs_counts.keys())

        results = []

        for name, ratio in matched_models().items():
            expected = sorted([total * r / sum(ratio) for r in ratio])
            observed = sorted([obs_counts.get(cat, 0) for cat in sorted_categories])
            
            if len(observed) != len(expected):
                continue
            chi2, p = chisquare(f_obs=observed, f_exp=expected)
            results.append((name, chi2, p, observed, expected))

        return sorted(results, key=lambda x: -x[2])  # sort by p-value desc

    @output
    @render.ui
    def result_ui():
        if not parsed_data():
            return ui.p("Paste data to begin.")

        if not matched_models():
            return ui.p("❌ No Mendelian model matches the number of phenotypic categories.")
        print(test_results())
        best = test_results()[0]
        name, chi2, p, obs, exp = best
        return ui.div(
            ui.h4("✔️ Best-Fitting Mendelian Model"),
            ui.p(f"Model: **{name}**"),
            ui.p(f"Observed counts: {obs}"),
            ui.p(f"Expected counts: {[round(e, 2) for e in exp]}"),
            ui.p(f"Chi-square statistic: {chi2:.4f}"),
            ui.p(f"P-value: {p:.4f}")
        )

# Create and run the app
app = App(app_ui, server)
