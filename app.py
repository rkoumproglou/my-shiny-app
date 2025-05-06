from shiny import App, reactive, render, ui
import numpy as np
import pandas as pd
from scipy.stats import chisquare
import plotnine as p9

# Define segregation models and their interpretations
models = {
    "3:1": {
        "ratios": [3, 1],
        "interpretation": "Single gene with dominant/recessive alleles (monohybrid cross)."
    },
    "1:2:1": {
        "ratios": [1, 2, 1],
        "interpretation": "Incomplete dominance or codominance."
    },
    "1:1:1:1": {
        "ratios": [1, 1, 1, 1],
        "interpretation": "Dihybrid test cross."
    },
    "9:3:3:1": {
        "ratios": [9, 3, 3, 1],
        "interpretation": "Two independently assorting genes."
    },
    "12:3:1": {
        "ratios": [12, 3, 1],
        "interpretation": "Dominant epistasis."
    },
    "9:7": {
        "ratios": [9, 7],
        "interpretation": "Duplicate recessive epistasis."
    },
    "15:1": {
        "ratios": [15, 1],
        "interpretation": "Duplicate dominant epistasis."
    },
    "9:3:4": {
        "ratios": [9, 3, 4],
        "interpretation": "Recessive epistasis."
    },
    "13:3": {
        "ratios": [13, 3],
        "interpretation": "Dominant and recessive epistasis."
    },
    "9:6:1": {
        "ratios": [9, 6, 1],
        "interpretation": "Polymeric gene interaction."
    }
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
        'p_value': p_value
    }
    return result

def compare_models(observed_counts):
    results = []
    for name, model in models.items():
        ratio = model['ratios']
        if len(ratio) == len(observed_counts):
            result = test_segregation(observed_counts, ratio)
            result['model'] = name
            result['interpretation'] = model['interpretation']
            result['ratio'] = ratio
            results.append(result)
    if not results:
        return None
    best_result = max(results, key=lambda x: x['p_value'])
    return best_result

app_ui = ui.page_fluid(
    ui.h2("Genetic Segregation Ratio Tester"),
    ui.input_text_area("phenotypes", "Enter phenotype data (comma-separated or pasted from Excel column)", placeholder="e.g. red, red, orange, red"),
    ui.output_ui("result_ui")
)

def server(input, output, session):

    @reactive.Calc
    def observed_counts():
        try:
            raw_input = input.phenotypes()
            cleaned = [x.strip() for x in raw_input.replace("\n", ",").split(",") if x.strip() != ""]
            values, counts = np.unique(cleaned, return_counts=True)
            return values.tolist(), counts.tolist()
        except Exception:
            return None, None

    @output
    @render.ui
    def result_ui():
        phenotypes, counts = observed_counts()
        if not counts:
            return ui.p("Please enter valid phenotype data.")

        result = compare_models(counts)
        if result is None:
            return ui.p("No model matches the number of observed phenotypic classes.")

        if result['p_value'] < 0.05:
            message = ui.div(
                ui.card(
                    ui.h4("Best-Fit Model Report", style="background-color: lightgreen; padding: 10px;"),
                    ui.p("None of the tested genetic models adequately explains the observed segregation ratio (p < 0.05)."),
                    style="box-shadow: 2px 2px 10px #888;"
                )
            )
            return message

        # Prepare plot data
        df = pd.DataFrame({
            "Phenotype": phenotypes * 2,
            "Count": result['observed'] + result['expected'],
            "Type": ["Observed"] * len(result['observed']) + ["Expected"] * len(result['expected'])
        })

        plot = (
            p9.ggplot(df, p9.aes(x='Phenotype', y='Count', fill='Type')) +
            p9.geom_bar(stat='identity', position=p9.position_dodge(width=0.8)) +
            p9.theme_minimal() +
            p9.ggtitle(f"Observed vs Expected Counts ({result['model']})") +
            p9.scale_fill_manual(values=["#FF9999", "#9999FF"]) if result['p_value'] >= 0.05 else
            p9.ggplot(df, p9.aes(x='Phenotype', y='Count', fill='Type')) +
            p9.geom_bar(stat='identity', position=p9.position_dodge(width=0.8)) +
            p9.theme_minimal() +
            p9.ggtitle("Observed vs Expected Counts") +
            p9.scale_fill_manual(values=["#FF9999", "#9999FF"])
        )

        content = ui.card(
            ui.h4("Best-Fit Model Report", style="background-color: lightgreen; padding: 10px;"),
            ui.h5(f"Best Fitting Model: {result['model']} - {result['interpretation']}"),
            ui.p("The observed ratio does not significantly deviate from the expected ratio based on the Chi-square test (p ≥ 0.05), so it is considered consistent with the best fit model."),
            ui.p(f"Observed counts: {result['observed']}"),
            ui.p(f"Expected counts: {np.round(result['expected'], 2).tolist()}"),
            ui.p(f"Chi-square Statistic: {result['chi_square_stat']:.4f}"),
            ui.p(f"P-value: {result['p_value']:.4f}"),
            ui.output_plot("plot"),
            style="box-shadow: 2px 2px 10px #888;"
        )
        return content

    @output
    @render.plot
    def plot():
        phenotypes, counts = observed_counts()
        result = compare_models(counts)
        if result is None or result['p_value'] < 0.05:
            return None

        df = pd.DataFrame({
            "Phenotype": phenotypes * 2,
            "Count": result['observed'] + result['expected'],
            "Type": ["Observed"] * len(result['observed']) + ["Expected"] * len(result['expected'])
        })

        p = (
            p9.ggplot(df, p9.aes(x='Phenotype', y='Count', fill='Type')) +
            p9.geom_bar(stat='identity', position=p9.position_dodge(width=0.8)) +
            p9.theme_minimal() +
            p9.ggtitle(f"Observed vs Expected Counts ({result['model']})") +
            p9.scale_fill_manual(values=["#FF9999", "#9999FF"])
        )
        return p

app = App(app_ui, server)
